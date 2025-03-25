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
// Software that calculates statistical properties of spectra
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
#include "T2toolkit.h"
#include "TKfit.h"

#define VNUM "v0.1"

void loadBaselineRegions(char *fname, float *baselineX1,float *baselineX2,int *nBaseline);

void help()
{
  printf("sdhdf_statistics\n\n");
  printf("-h            This help\n");
  printf("-sb <sb>      Select sub-band number sb\n");
  printf("-fr <f0> <f1> Define frequency range (can also use -freqRange)\n");
  printf("-nfit <val>   Set the polynomial order for least squares fitting\n");
  printf("-scaleRMS     Produces the rms as a function of N where N = number of dumps averaged\n");
  printf("\n\n");
  printf("Example\n");
  printf("sdhdf_statistics -sb 5 -fr 1400 1410 *.hdf\n");
  exit(1);
}


int main(int argc,char *argv[])
{
  int        i,j,k,l,ii,kk,ll;
  int nAv;
  char       fname[MAX_FILES][64];
  sdhdf_fileStruct *inFile;
  int        idump,iband,ibeam;
  int av=0;
  int sump=0;
  int nx=1;
  int ny=1;
  int polPlot=1;
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
  int nfit=0;
  float fx[2],fy[2];
  float fc;
  double sx,sx2,sx_2,sx2_2,sdev1,sdev2,mean1,mean2,min1,min2,max1,max2;
  float *avData1,*avData2;
  double freq;
  int np,t;
  int noBaseline=0;
  int scaleRMS=0;
  int nc;
  float freq0;
  float freq1;
  int setFreq=0;
  
  //  loadBaselineRegions("cleanRegions.dat",baselineX1,baselineX2,&nBaseline);
  
  // Defaults
  idump = iband = ibeam = 0;
  strcpy(grDev,"/xs");
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  //iband = 5;
  iband=0;
  //  iband = 5;
  
  sdhdf_initialiseFile(inFile);

  inFiles=0;

  // Have options for mean, median, peak
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-sb")==0)  // FIX ME -- WORK OUT BAND FROM FREQUENCY RANGE
	{
	  sscanf(argv[++i],"%d",&iband);
	  printf("Setting iband to %d\n",iband);
	}
	  else if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcasecmp(argv[i],"-freqRange")==0 || strcasecmp(argv[i],"-fr")==0)
	{
	  sscanf(argv[++i],"%f",&freq0);
	  sscanf(argv[++i],"%f",&freq1);
	  setFreq=1;
	}
      else if (strcmp(argv[i],"-nfit")==0) 
	sscanf(argv[++i],"%d",&nfit);
      else if (strcasecmp(argv[i],"-scaleRMS")==0)
	scaleRMS=1;
      else
	strcpy(fname[inFiles++],argv[i]);
    }
  if (setFreq==1)
    printf("\n Results for frequency range %g to %g MHz\n\n\n",freq0,freq1);
  printf("------------------------------------------------------------------------------------------------------\n");
  printf("Label File Source Dump Mean_pol1 Mean_pol2 Sdev_pol1 Sdev_pol2 Min_pol1 Max_pol1 Min_pol2 Max-Pol2 MJD\n");
  printf("------------------------------------------------------------------------------------------------------\n");
  
  for (i=0;i<inFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);

      sdhdf_loadMetaData(inFile);
      if (setFreq==0)
	{
	  freq0 = inFile->beam[ibeam].bandHeader[iband].f0;
	  freq1 = inFile->beam[ibeam].bandHeader[iband].f1;
	  //	  printf("Setting frequency to %g %g\n",freq0,freq1);
	}
      
      sdhdf_loadBandData(inFile,ibeam,iband,1);
      // Process each band
      //      for (j=0;j<inFile->beam[ibeam].nBand;j++)
      j=iband; // FIX ME -- UPDATE TO PROCESS EACH BAND
      {
	// Process each dump
	if (scaleRMS == 1)
	  {
	    avData1 = (float *)calloc(sizeof(float),inFile->beam[ibeam].bandHeader[j].nchan);
	    avData2 = (float *)calloc(sizeof(float),inFile->beam[ibeam].bandHeader[j].nchan);
	  }
	for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l++) 
	  {
	    sx=sx2=sx_2=sx2_2=0.0;
	    np=0;
	    nc=0;
	    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
	      {
		// FIX ME: using [0] for frequency dump
		freq = inFile->beam[ibeam].bandData[j].astro_data.freq[idump*inFile->beam[ibeam].bandHeader[j].nchan+k];
		if (freq > freq0 && freq < freq1)
		  {
		    val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[(l)*inFile->beam[ibeam].bandHeader[j].nchan+k];
		    val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[(l)*inFile->beam[ibeam].bandHeader[j].nchan+k];
		    if (scaleRMS==1)
		      {
			avData1[nc] += val1;
			avData2[nc] += val2;			
		      }
		    nc++;
		    if (np==0)
		      {
			min1=max1 = val1;
			min2=max2 = val2;
		      }
		    else
		      {
			if (min1 > val1) min1 = val1;
			if (max1 < val1) max1 = val1;
			if (min2 > val2) min2 = val2;
			if (max2 < val2) max2 = val2;
		      }
		    sx+=val1;
		    sx2+=val2;
		    sx_2 += pow(val1,2);
		    sx2_2 += pow(val2,2);
		    np++;
		  }
	      }
	    mean1=sx/(double)np;
	    mean2=sx2/(double)np;

	    sdev1 = sqrt(1.0/(double)np * sx_2 - pow(1.0/(double)np * sx,2));
	    sdev2 = sqrt(1.0/(double)np * sx2_2 - pow(1.0/(double)np * sx2,2));
	    
	    printf("[stats] %s %s %d %g %g %g %g %g %g %g %g %.5f\n",inFile->fname,inFile->beamHeader[0].source,l,mean1,mean2,sdev1,sdev2,min1,max1,min2,max2,inFile->beam[0].bandData[j].astro_obsHeader[l].mjd);	  
	    if (scaleRMS == 1)
	      {
		float calcVals1[nc],calcVals2[nc];
		float fx[nc];
		float rms1,rms2;
		printf("Scaling with nc = %d\n",nc);
		for (i=0;i<nc;i++)
		  {
		    fx[i] = i; // Probably should use the frequency instead here -- FIX ME
		    calcVals1[i] = avData1[i]/(double)(l+1);
		    calcVals2[i] = avData2[i]/(double)(l+1);
		  }
		if (nfit > 0)
		  {
		    TKremovePoly_f(fx,calcVals1,nc,nfit);
		    TKremovePoly_f(fx,calcVals2,nc,nfit);
		  }
		rms1 = TKfindRMS_f(calcVals1,nc);
		rms2 = TKfindRMS_f(calcVals2,nc);
		printf("[rmsScale] %d %g %g\n",l+1,rms1,rms2); // Note not checking all dumps are the same length -- FIX ME
	      }
	  }
	if (scaleRMS == 1)
	  {
	    free(avData1);
	    free(avData2);
	  }
	
      }
           
      sdhdf_closeFile(inFile);

    }

  
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
