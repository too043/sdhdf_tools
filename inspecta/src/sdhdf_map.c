//  Copyright (C) 2019, 2020, 2021, 2022 George Hobbs

/*
 *    This file is part of sdhdfProc. 
 * 
 *    sdhdfProc is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    sdhdfProc is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with sdhdfProc.  If not, see <http://www.gnu.org/licenses/>. 
 */

//
// Software to plot a map
//
// Usage:
// sdhdf_map <filenames>


//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"
#include "TKfit.h"
#include <cpgplot.h>


#define MAX_CHANS 262144
#define VNUM "v0.1"

double haversine(double centre_long,double centre_lat,double src_long,double src_lat);

void help()
{
  printf("sdhdf_map: routine to map a region of the sky\n\n");
  printf("-b <n>     Select band number <n>\n");
  printf("-h         This help\n");
  printf("-highlight <X> <Y> highlight point at X and Y coordinate\n");
  printf("-poly <n>  Remove n'th order polynomial\n");
  printf("\n\n");
  printf("filenames are given on the command line, e.g.:\n");
  printf("sdhdf_map -b 5 uwl*.hdf\n");
}

int main(int argc,char *argv[])
{
  int        i,j,k,l,ii,jj,kk;
  int iband=0;
  char       fname[MAX_FILES][64];
  FILE *fout;
  sdhdf_fileStruct *inFile;
  int        idump,ibeam;
  int av=0;
  int sump=0;
  int polPlot=1;
  float highlightX[1024],highlightY[1024];
  int nHighlight=0;
  float size;
  float val1,val2;
  long writepos;
  float mx,my;
  float minx,maxx,miny,maxy;
  char key;
  int allocateMemory=0;
  int inFiles;
  char grDev[128];
  float ra,dec;
  double sum[1024];  // Should be MAX_DUMPS
  double xfit[1024];
  int   nc = 0;

  int   nx = 64;
  int   ny = 32;


  // float x0 = 270.2;
  //  float x1 = 271.6;
  //  float y0 = -30.7;
  //  float y1 = -29.4;

  float x0 = 60.835;
  float x1 = 63.328;
  float y0 = -66.3;
  float y1 = -65.3;

  /*
  int nx = 32;
  int ny = 16;
  float x0 = 270.6;
  float x1 = 271.4;
  float y0 = -30.0;
  float y1 = -29.6;
  */

  float angle;
  float *arr;
  float tr[6];
  float sigma;
  float min,max;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float *wt;
  float weight;
  float f0,lambda;
  float telDiameter = 64;
  float fwhm;

  int flagDump;
  
  int nPos=0;
  float raVal[4096],decVal[4096];

  int removeFunction = 0;
  int nd;
  

  
  // Defaults
  idump = iband = ibeam = 0;
  strcpy(grDev,"/xs");
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);

  inFiles=0;
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-poly")==0)
	sscanf(argv[++i],"%d",&removeFunction);
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-b")==0)
	sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-highlight")==0)
	{
	  sscanf(argv[++i],"%f",&highlightX[nHighlight]);
	  sscanf(argv[++i],"%f",&highlightY[nHighlight++]);	  
	}
      else
	strcpy(fname[inFiles++],argv[i]);
    }

  arr = (float *)calloc(sizeof(float),nx*ny);
  wt  = (float *)calloc(sizeof(float),nx*ny);
  tr[0] = x0; tr[1] = (x1-x0)/(float)nx; tr[2] = 0;
  tr[3] = y0; tr[4] = 0; tr[5] = (y1-y0)/(float)ny; 

  min = 1e30;
  max = -1e30;
  
  for (i=0;i<inFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_loadMetaData(inFile);

      // Process each band
      //      for (j=0;j<inFile->beam[ibeam].nBand;j++)
      //      j=12;
      //      j=8;
      //      j = 14;
      j=iband;
      //      j=16;
      sdhdf_loadBandData(inFile,ibeam,j,1);
      f0     = inFile->beam[ibeam].bandHeader[j].fc;
      lambda = 3.0e8/(f0*1e6);
      fwhm   = 1.22*lambda*180/M_PI/telDiameter; 
      sigma  = fwhm/2.35;
      printf("Sigma = %g, fwhm = %g, lambda = %g, fc = %g\n",sigma,fwhm,lambda,f0);
      {
	printf("Processing band: %d file = %s\n",j,fname[i]);

	// Process each dump
	nd=0;
	for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l++) 
	  {
	    flagDump=0;

	    if (flagDump==0)
	      {
		xfit[nd] = l;
		sum[nd]  = 0;
		nc=0;
		for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		  {
		    // FIX ME: using [0] for frequency dump
		    f0 = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
		    // Note: not flagged (YET) on a per-dump basis
		    //		if (f0 > 1419 && f0 < 1422)
		    //		if (f0 > 1370 && f0 < 1410)
		    //		if (f0 > 2000 && f0 < 2100)
		    {
		      if (inFile->beam[ibeam].bandData[j].astro_data.dataWeights[k] != 0)
			{
			  val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			  val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			  sum[nd]+=val1+val2;			  
			  nc++;
			}
		    }
		  }

		if (nc == 0)
		  {
		    printf("ERROR: all data flagged\n");
		    exit(1);
		  }
		
		sum[nd]/=(float)nc;
		nd++;
	      }
	  }
      
	printf("Removing function\n");
	if (removeFunction > 0)
	  TKremovePoly_d(xfit,sum,inFile->beam[ibeam].bandHeader[j].ndump,removeFunction);

	for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l++)
	  {
	    flagDump=0;
	    if (strcmp(fname[i],"uwl_200720_024617.hdf.autoflag")==0 && (l>=105)){flagDump=1;}
	      //&& (l>=20)) 
	    
	    if (flagDump==0)
	      {
		ra  = inFile->beam[ibeam].bandData[j].astro_obsHeader[l].raDeg;
		dec = inFile->beam[ibeam].bandData[j].astro_obsHeader[l].decDeg;
		
		raVal[nPos]  = ra;
		decVal[nPos] = dec;
		nPos++;
		
		for (jj=0;jj<ny;jj++)
		  {
		    for (ii=0;ii<nx;ii++)
		      {
			// Should use a lookup table?
			angle = haversine(ra,dec,x0+(x1-x0)*ii/(float)nx,y0+(y1-y0)*jj/(float)ny);
			// Should the weight be squared or not?
			// weight = pow(exp(-(angle*angle)/2.0/sigma/sigma),2);
			weight = exp(-(angle*angle)/2.0/sigma/sigma);
			arr[jj*nx+ii] += sum[l]*weight;
			wt[jj*nx+ii]  += weight;
		      }
		  }
	      }
	  }
      }
      
      sdhdf_closeFile(inFile);     
    }
  fout = fopen("weightVals.dat","w");
  for (jj=0;jj<ny;jj++)
    {
      for (ii=0;ii<nx;ii++)
	{	  
	  fprintf(fout,"%d %d %g\n",ii,jj,wt[jj*nx+ii]*10);
	}
      fprintf(fout,"\n");
    }
  fclose(fout);
  for (ii=0;ii<nx*ny;ii++)
    {
      if (wt[ii] > 0)
	arr[ii]/=wt[ii];
	//			arr[ii]/=1.;
      else
	{
	  printf("WARNING: weighting = 0\n");
	}
      if (min > arr[ii]) min = arr[ii];
      if (max < arr[ii]) max = arr[ii];
    }
  printf("min = %g, max = %g\n",min,max);

  cpgbeg(0,grDev,1,1);
  cpgsch(1.4);
  cpgslw(2);
  cpgscf(2);
  cpgask(0);
  minx = x0;  miny = y0;
  maxx = x1;  maxy = y1;
  do {
    cpgenv(minx,maxx,miny,maxy,0,1);
    cpglab("Right ascension (deg)","Declination (deg)","");
    
    cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
    cpgimag(arr,nx,ny,1,nx,1,ny,min,max,tr);
    cpgsci(3);
    cpgpt(nPos,raVal,decVal,15);
    cpgsci(1);
    cpgsci(4);
    cpgpt(nHighlight,highlightX,highlightY,18);
    cpgsci(1);
    cpgcurs(&mx,&my,&key);
    if (key=='s')
      {
	printf("Enter min (%f) ",min);
	scanf("%f",&min);
	printf("Enter max (%f) ",max);
	scanf("%f",&max);
      }
    else if (key=='z')
      {
	float mx2,my2;
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	if (mx != mx2 && my != my2)
	  {
	    if (mx < mx2)
	      {minx = mx; maxx = mx2;}
	    else
	      {minx = mx2; maxx = mx;}
	    
	    if (my < my2)
	      {miny = my; maxy = my2;}
	    else
	      {miny = my2; maxy = my;}
	  }
      }
  } while (key!='q');
  cpgend();
  
  free(arr);
  free(wt);
}


double haversine(double centre_long,double centre_lat,double src_long,double src_lat)
{
  double dlon,dlat,a,c;
  double deg2rad = M_PI/180.0;

  centre_long*=deg2rad;
  centre_lat*=deg2rad;
  src_long*=deg2rad;
  src_lat*=deg2rad;
  
  /* Apply the Haversine formula */
  dlon = (src_long - centre_long);
  dlat = (src_lat  - centre_lat);
  a = pow(sin(dlat/2.0),2) + cos(centre_lat) *
    cos(src_lat)*pow(sin(dlon/2.0),2);
  if (a==1)
    c = M_PI;
  else
    c = 2.0 * atan2(sqrt(a),sqrt(1.0-a));
  return c/deg2rad;
}
