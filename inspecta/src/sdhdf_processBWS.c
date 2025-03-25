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
// Software that processes the spectral line data from the Baade's window survey
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"
#include <cpgplot.h>
#include "TKfit.h"

#define VNUM "v0.1"
#define MAX_CHANS 262144
#define MAX_DUMPS 1024
#define MAX_TRANSMITTERS 1024

void loadBaselineRegions(char *fname, float *baselineX1,float *baselineX2,int *nBaseline);

int main(int argc,char *argv[])
{
  int        i,j,k,l,ii,kk;
  char       fname[MAX_FILES][64];
  sdhdf_fileStruct *inFile;
  int        idump,iband,ibeam;
  int av=0;
  int sump=0;
  int nx=1;
  int ny=1;
  int polPlot=1;
  float *highResX,*highResY1,*highResY2;
  float *zoomX,*zoomY1,*zoomY2;
  float *zoomY12;
  double *bX,*bY1,*bY2;
  int *bI;
  int nB;
  long writepos;
  int allocateMemory=0;
  float miny,maxy;
  float f0,f1;
  int inFiles;
  char text[1024];
  int sumPol_nPlot;
  char grDev[128];
  int nZoom;
  float xp1,xp2,yp1,yp2;
  float val1,val2;
  float baselineX1[1024],baselineX2[1024];
  int nBaseline;
  double params1[16],params2[16],v[16];
  int nfit;
  float fx[2],fy[2];
  float fc;
  double sx,sx2,sdev;
  int noBaseline=0,nc;
  
  loadBaselineRegions("cleanRegions.dat",baselineX1,baselineX2,&nBaseline);
  
  highResX = (float *)malloc(sizeof(float)*MAX_CHANS);
  highResY1 = (float *)calloc(sizeof(float),MAX_CHANS);
  highResY2 = (float *)calloc(sizeof(float),MAX_CHANS);
  zoomX = (float *)malloc(sizeof(float)*MAX_CHANS);
  zoomY1 = (float *)calloc(sizeof(float),MAX_CHANS);
  zoomY2 = (float *)calloc(sizeof(float),MAX_CHANS);
  zoomY12 = (float *)calloc(sizeof(float),MAX_CHANS);
  bI = (int *)malloc(sizeof(int)*MAX_CHANS);
  bX = (double *)malloc(sizeof(double)*MAX_CHANS);
  bY1 = (double *)calloc(sizeof(double),MAX_CHANS);
  bY2 = (double *)calloc(sizeof(double),MAX_CHANS);

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
      strcpy(fname[inFiles++],argv[i]);
    }
  
  for (i=0;i<inFiles;i++)
    {
      sprintf(grDev,"%s.bws.ps/cps",fname[i]);
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);

      sdhdf_loadMetaData(inFile);

      cpgbeg(0,grDev,1,1);
      cpgsch(0.8);
      cpgslw(2);
      cpgscf(2);
      cpgask(0);

      //      cpgsvp(0.15,0.95,0.15,0.95);
      //      cpgswin(0,1,0,1);
      //      cpgpage();
      
      for (j=0;j<inFile->beam[ibeam].nBand;j++)
	{
	  sdhdf_loadBandData(inFile,ibeam,j,1);
	}
      

      // Process each band
      for (j=0;j<inFile->beam[ibeam].nBand;j++)
      // j=0;
      //      j=5;
      //      j=17;
      {
	  printf("Processing band: %d\n",j);

	  // Process each dump
	  //	  for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l++) 
	    l=2;
	  {	      	      	      
	    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
	      {
		// FIX ME: using [0] for frequency dump
		highResX[k] = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
		val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
		val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
		
		highResY1[k] = val1;
		highResY2[k] = val2;
	      }

	    xp1 = 0.15;
	    yp1 = 0.15;
	    xp2 = xp1+0.2;
	    yp2 = yp1+0.15;
	    for (k=0;k<16;k++) // 16 = Number of zoom bands
	      {
		nZoom=0;
		f0 = 704+j*128 + k/16.*128;
		f1 = f0 + 128.0/16.0;
		fc = (f0+f1)/2.;
		nB=0;
		for (kk=0;kk<inFile->beam[ibeam].bandHeader[j].nchan;kk++)
		  {
		    if (highResX[kk] >= f0 && highResX[kk] < f1)
		      {
			zoomX[nZoom]  = highResX[kk];
			zoomY1[nZoom] = highResY1[kk];
			zoomY2[nZoom] = highResY2[kk];
			
			
			// Check for baseline region
			for (ii=0;ii<nBaseline;ii++)
			  {
			    if (zoomX[nZoom] >= baselineX1[ii] && zoomX[nZoom] < baselineX2[ii])
			      {
				bI[nB] = nZoom;
				bX[nB] =  (double)zoomX[nZoom]-fc; 
				bY1[nB] = (double)zoomY1[nZoom];
				bY2[nB] = (double)zoomY2[nZoom];
				nB++;
				break;
			      }
			  }
			nZoom++;				 
		      }	      
		  }
		if (nB==0) // No baseline region
		  {
		    noBaseline=1;
		    for (kk=0;kk<inFile->beam[ibeam].bandHeader[j].nchan;kk++)
		      {
			if (highResX[kk] >= f0 && highResX[kk] < f1)
			  {
			    bI[nB] = nZoom;
			    bX[nB] =  (double)zoomX[nZoom]-fc; 
			    bY1[nB] = (double)zoomY1[nZoom];
			    bY2[nB] = (double)zoomY2[nZoom];
			    nB++;			    
			  }
		      }
		  }
		else
		  noBaseline=0;
		// Remove the baseline

		if (j==0 && (k==0 || k==1 || k==3 || k==4 || k==7))
		  nfit = 2;
		else
		  nfit = 8;
		
		TKleastSquares_svd_noErr(bX,bY1,nB,params1,nfit,TKfitPoly);
		TKleastSquares_svd_noErr(bX,bY2,nB,params2,nfit,TKfitPoly);
		printf("Params1 are %g %g %g\n",params1[0],params1[1],params1[2]);
		printf("Params2 are %g %g %g\n",params2[0],params2[1],params2[2]);
		
		miny = 1e99; maxy = -1e99;
		for (ii=0;ii<nZoom;ii++)
		  {
		    TKfitPoly(zoomX[ii]-fc,v,nfit);
		    for (kk=0;kk<nfit;kk++)
		      {
			zoomY1[ii] -= v[kk]*params1[kk]; 
			zoomY2[ii] -= v[kk]*params2[kk];
		      }
		    zoomY12[ii] = zoomY1[ii] + zoomY2[ii];
		    if (miny > zoomY12[ii]) miny = zoomY12[ii];
		    if (maxy < zoomY12[ii]) maxy = zoomY12[ii];
		    /*
		      if (miny > zoomY1[ii]) miny = zoomY1[ii];
		      if (miny > zoomY2[ii]) miny = zoomY2[ii];
		      if (maxy < zoomY1[ii]) maxy = zoomY1[ii];
		      if (maxy < zoomY2[ii]) maxy = zoomY2[ii];
		    */
		    
		  }


		// Get variance in the baseline region		
		sx=sx2=0;
		nc = 0;
		for (kk=0;kk<nB;kk++)
		  {
		    sx  += zoomY12[bI[kk]];
		    sx2 += pow(zoomY12[bI[kk]],2);
		    nc++;
		  }	       		
		
		sdev = sqrt(1./(double)nc*sx2 - pow(1.0/(double)nc * sx,2));
		
		printf("From %g %g %g %g (nB = %d, nZoom = %d)\n",xp1,xp2,yp1,yp2,nB,nZoom);
		printf("min max = %g %g\n",miny,maxy);
		cpgsvp(xp1,xp2,yp1,yp2);
		//		cpgswin(f0,f1,miny,maxy);
		cpgswin(f0,f1,-5*sdev,10*sdev);
		cpgbox("ABCTSN",0,0,"ABC",0,0);
		if (noBaseline==1)
		  {
		    cpgsci(2); cpgline(nZoom,zoomX,zoomY12); cpgsci(1);
		  }
		else
		  {
		    cpgsci(4); cpgline(nZoom,zoomX,zoomY12); cpgsci(1);
		  }
		//		  cpgsci(4); cpgline(nZoom,zoomX,zoomY1); cpgsci(1); 
		//		  cpgsci(2); cpgline(nZoom,zoomX,zoomY2); cpgsci(1);
		
		fx[0] = f0; fx[1] = f1;
		fy[0] = fy[1] = 0;
		cpgsci(3); cpgsls(2); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
		xp1 += 0.2;
		if (xp1 > 0.8)
		  {
		    xp1 = 0.15;
		    yp1 += 0.2;
		  }
		xp2 = xp1 + 0.2;
		yp2 = yp1 + 0.15;
	      }	      
	    cpgpage();
	  }
      }
     

      cpgend();  
 
      sdhdf_closeFile(inFile);


    }

  
  free(highResX); free(highResY1); free(highResY2);
  free(zoomX); free(zoomY1); free(zoomY2);
  free(bX); free(bY1); free(bY2);
  free(bI);
  free(zoomY12);
}


void loadBaselineRegions(char *fname, float *baselineX1,float *baselineX2,int *nBaseline)
{
  FILE *fin;

  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open >%s<\n",fname);
      exit(1);
    }
  *nBaseline = 0;
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f",&baselineX1[*nBaseline],&baselineX2[*nBaseline])==2)
	(*nBaseline)++;
    }
  fclose(fin);
}
