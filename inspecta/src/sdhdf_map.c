//  Copyright (C) 2019, 2020, 2021, 2022 George Hobbs

/*
 *    This file is part of INSPECTA. 
 * 
 *    INSPECTA is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    INSPECTA is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with INSPECTA.  If not, see <http://www.gnu.org/licenses/>. 
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
  printf("-xrange <x0> <x1> <dx> X-range for map\n");
  printf("-yrange <y0> <y1> <dy> Y-range for map\n");
  printf("-singleCh <chOn> <chOff> Obtain power from a single channel and baseline from chOff\n");
  printf("-b <n>     Select band number <n>\n");
  printf("-bl <f0> <f1> Remove baseline based on this frequency range\n");
  printf("-bl_min       Remove baseline based on minimum value\n");
  printf("-divideWeights Divide by the weightings\n");
  printf("-g <grdev> Set the graphics device\n");
  printf("-h         This help\n");
  printf("-highlight <X> <Y> <label> highlight point at X and Y coordinate\n");
  printf("-label     Add beam labels\n");
  printf("-log       Take log of z-values\n");
  printf("-overlay <filename> Overlay positions of sources given in the data file: format <source> <raj (deg)> <decj (deg)>\n");
  printf("-pol <val> Where val = 0 (default) for AA+BB, = 1 for AA and = 2 for BB\n");
  printf("-poly <n>  Remove n'th order polynomial from the time series produced over all spectral dumps\n");
  printf("-stat <n>  0 [default] = sum the spectrum, 1 = take the highest value\n");
  printf("-title <str> Title for plot\n");
  printf("\n\n");
  printf("filenames are given on the command line, e.g.:\n");
  printf("sdhdf_map -b 5 uwl*.hdf\n");
}

int main(int argc,char *argv[])
{
  int        i,j,k,l,ii,jj,kk;
  int iband=0;
  char       fname[MAX_FILES][64];
  char overlayName[1024]="unset";  
  FILE *fout;
  sdhdf_fileStruct *inFile;
  int        idump,ibeam;
  int divideWeights=0;
  int av=0;
  int stat=0;
  int sCh0=-1;
  int sCh1=-1;
  int sump=0;
  int polPlot=1;
  float highlightX[1024],highlightY[1024];
  char highlightText[1024][32];
  int nHighlight=0;
  float size;
  char grDev[128];
  float val1,val2;
  long writepos;
  float mx,my;
  float omx,omy;
  float minx,maxx,miny,maxy;
  char key;
  int allocateMemory=0;
  int inFiles;
  float ra,dec;
  double sum[1024];  // Should be MAX_DUMPS
  double xfit[1024];
  int   nc = 0;

  int   nx = 64;
  int   ny = 64;


  // float x0 = 270.2;
  //  float x1 = 271.6;
  //  float y0 = -30.7;
  //  float y1 = -29.4;

//  float x0 = 60.835;
//  float x1 = 63.328;
//  float y0 = -66.3;
//  float y1 = -65.3;

  float x0 = 230.5;
  float x1 = 234;
  float y0 = -57.5;
  float y1 = -55.5;

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
  int setMinMax=0;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float *wt;
  float weight;
  float f0,lambda;
  float telDiameter = 64;
  float fwhm;

  float bl0,bl1;
  float val_bl0,val_bl1;
  int bl_count=0;
  
  int flagDump;
  
  int nPos=0;
  float raVal[4096],decVal[4096];
  int pol=0;
  
  int removeFunctionDump = 0;
  int nd;
  int log=0;

  int label=0;
  char labelText[1024];
  char title[1024]="";
  float bl_min=-1;
  
  bl0=bl1=-1;
  
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
  min = 1e30;
  max = -1e30;

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-polyDump")==0)
	sscanf(argv[++i],"%d",&removeFunctionDump);
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-stat")==0)
	sscanf(argv[++i],"%d",&stat);
      else if (strcmp(argv[i],"-overlay")==0)
	strcpy(overlayName,argv[++i]);
      else if (strcmp(argv[i],"-singleCh")==0)
	{
	  sscanf(argv[++i],"%d",&sCh0);
	  sscanf(argv[++i],"%d",&sCh1);
	}
      else if (strcmp(argv[i],"-divideWeights")==0)
	divideWeights=1;
      else if (strcmp(argv[i],"-title")==0)
	strcpy(title,argv[++i]);
      else if (strcmp(argv[i],"-xrange")==0)
	{
	  sscanf(argv[++i],"%f",&x0);
	  sscanf(argv[++i],"%f",&x1);
	  sscanf(argv[++i],"%d",&nx);
	}
      else if (strcmp(argv[i],"-label")==0)
	label=1;
      else if (strcmp(argv[i],"-yrange")==0)
	{
	  sscanf(argv[++i],"%f",&y0);
	  sscanf(argv[++i],"%f",&y1);
	  sscanf(argv[++i],"%d",&ny);
	}
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-b")==0)
	sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-log")==0)
	log=1;
      else if (strcmp(argv[i],"-pol")==0)
	sscanf(argv[++i],"%d",&pol);
      else if (strcmp(argv[i],"-bl")==0)
	{
	  sscanf(argv[++i],"%f",&bl0);
	  sscanf(argv[++i],"%f",&bl1);
	}
      else if (strcmp(argv[i],"-bl_min")==0)
	bl_min = 1;
      else if (strcmp(argv[i],"-highlight")==0)
	{
	  sscanf(argv[++i],"%f",&highlightX[nHighlight]);
	  sscanf(argv[++i],"%f",&highlightY[nHighlight]);
	  strcpy(highlightText[nHighlight++],argv[++i]);
	}
      else if (strcmp(argv[i],"-zrange")==0)
	{
	  sscanf(argv[++i],"%f",&min);
	  sscanf(argv[++i],"%f",&max);
	  setMinMax=1;
	}      
      else
	strcpy(fname[inFiles++],argv[i]);
    }

  if (strcmp(overlayName,"unset")!=0)
    {
      FILE *fin;
      char stmp[1024];
      if (!(fin = fopen(overlayName,"r")))
	  printf("Unable to open overlay file: %s\n",overlayName);
      else
	{
	  while (!feof(fin))
	    {
	      if (fscanf(fin,"%s %f %f\n",highlightText[nHighlight],&highlightX[nHighlight],&highlightY[nHighlight])==3)
		nHighlight++;
	    }
	}
    }
  arr = (float *)calloc(sizeof(float),nx*ny);
  wt  = (float *)calloc(sizeof(float),nx*ny);
  tr[0] = x0; tr[1] = (x1-x0)/(float)nx; tr[2] = 0;
  tr[3] = y0; tr[4] = 0; tr[5] = (y1-y0)/(float)ny; 

  
  for (i=0;i<inFiles;i++)
  //  i=30;
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
		// Remove baseline if needed
		val_bl0=val_bl1=0;
		bl_count=0;
		if (sCh0>-1) // Use single channels
		  {
		    // No need to remove a baseline
		  }
		else if (bl_min == 1)
		  {
		    val_bl0 = val_bl1 = 1e30;
		    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		      {
			if (val_bl0 > inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k]) val_bl0 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			if (val_bl1 > inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k]) val_bl1 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			
		      }
		    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		      {
			inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k]-=(val_bl0);
			inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k]-=(val_bl1);
		      }		    		    
		  }		
		else if (bl0 > -1)
		  {
		    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		      {
			f0 = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
			if (f0 > bl0 && f0 < bl1)
			  {
			    val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			    val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];

			    val_bl0+=val1;
			    val_bl1+=val2;
			    bl_count++;
			  }
		      }
		    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		      {
			inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k]-=(val_bl0/bl_count);
			inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k]-=(val_bl1/bl_count);
		      }		    
		  }
		xfit[nd] = l;
		sum[nd]  = 0;
		nc=0;		
	      
		if (sCh0>-1) // Use single channels
		  {
		    nc = 1;
		    val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+sCh0] - inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+sCh1];
		    val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+sCh0] - inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+sCh1];
		    if (pol==0)
		      sum[nd]+=(val1+val2);
		    else if (pol==1)
		      sum[nd]+=(val1);
		    else if (pol==2)
		      sum[nd]+=(val2);	   
		  }
		else
		  {
		    if (stat==1)
		      sum[nd] = -1e30;
		    
		    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		      {
			f0 = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
			{
			  if (inFile->beam[ibeam].bandData[j].astro_data.dataWeights[k] != 0)
			    {
			      val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			      val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			      //			  printf("Loading %g %g\n",val1,val2);
			      if (stat==0)
				{
				  if (pol==0)
				    sum[nd]+=(val1+val2);
				  else if (pol==1)
				    sum[nd]+=(val1);
				  else if (pol==2)
				    sum[nd]+=(val2);
				}
			      else
				{
				  if (sum[nd] < (val1+val2)) // This has only been implmented for pol=0
				      sum[nd] = val1+val2;
				}
			      nc++;
			    }
			}
		      }
		  }
		if (nc == 0)
		  {
		    printf("ERROR: all data flagged\n");
		    exit(1);
		  }
		if (stat==0)
		  {
		    sum[nd]/=(float)nc;
		    if (log==1)
		      sum[nd] = log10(sum[nd]);
		  }
		nd++;
	      }
	  }
     

	if (removeFunctionDump > 0)
	  {
	    printf("Removing function from spectral dumps\n");
	    TKremovePoly_d(xfit,sum,inFile->beam[ibeam].bandHeader[j].ndump,removeFunctionDump);
	  }
	for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l++)
	  {
	    flagDump=0;
	    
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
			weight = pow(exp(-(angle*angle)/2.0/sigma/sigma),2);
			//			weight = exp(-(angle*angle)/2.0/sigma/sigma);
			
			//			printf("test: %d %d %g %g\n",ii,jj,angle,weight);
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
      if (divideWeights==1)
	{
	  
	  if (wt[ii] > 0)
	    {	  
	      arr[ii]/=wt[ii];
	      //			arr[ii]/=1.;
	    }
	  else
	    {
	      printf("WARNING: weighting = 0\n");
	    }
	}
      if (setMinMax == 0)
	{
	  if (min > arr[ii]) min = arr[ii];
	  if (max < arr[ii]) max = arr[ii];
	}
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
    cpglab("Right ascension (deg)","Declination (deg)",title);
    
    cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
    cpgimag(arr,nx,ny,1,nx,1,ny,min,max,tr);
    cpgsci(3);
    cpgpt(nPos,raVal,decVal,15);
    if (label==1)
      {
	cpgsch(1.0);
	for (i=0;i<nPos;i++)
	  {		
	    sprintf(labelText,"%d",i);
	    cpgtext(raVal[i],decVal[i],labelText);
	  }
	cpgsch(1.4);
      }
    cpgsci(1);
    cpgsci(4);
    cpgpt(nHighlight,highlightX,highlightY,18);
    cpgsch(0.8);
    for (i=0;i<nHighlight;i++)
      cpgtext(highlightX[i],highlightY[i],highlightText[i]);
    cpgsci(1);
    cpgsch(1.4);




    if (strcmp(grDev,"/xs")==0)
      {
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
		if (mx < mx2) {minx = mx; maxx = mx2;}
		else          {minx = mx2; maxx = mx;}	    
		if (my < my2) {miny = my; maxy = my2;}
		else          {miny = my2; maxy = my;}
	      }
	  }
	else if (key=='u')
	  {
	    minx = x0;  miny = y0;
	    maxx = x1;  maxy = y1;
	    
	  }
      }
    else
      key='q';
    
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
