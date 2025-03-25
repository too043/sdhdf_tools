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
#include <cpgplot.h>


#define MAX_CHANS 262144
#define MAX_SPECTRA 262144 // Should be NFILES * NDUMP

#define VNUM "v0.1"

double haversine(double centre_long,double centre_lat,double src_long,double src_lat);

int main(int argc,char *argv[])
{
  int        i,j,k,l,ii,jj,kk;
  char       fname[MAX_FILES][64];
  sdhdf_fileStruct *inFile;
  int        idump,ibeam;
  int av=0;
  int sump=0;
  int polPlot=1;
  float **highResX,**highResY1,**highResY2;
  float *avSpectraX,*avSpectraY1,*avSpectraY2;
  float *plotX,*plotX3,*plotY1,*plotY2,*plotY3;
  int nAv = 16;
  int avCount=0;
  float size;
  float val1,val2;
  long writepos;
  int allocateMemory=0;
  float miny,maxy,miny2,maxy2;
  int inFiles;
  char grDev[128];
  float ra[MAX_SPECTRA],dec[MAX_SPECTRA];
  int sdump[MAX_SPECTRA];
  char  **filename;
  float f0[MAX_SPECTRA];
  int nChan[MAX_SPECTRA];

  
  float sum;
  float ra0,dec0;
  float angle;
  float tr[6];
  float sigma;
  float min,max;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float weight;
  float freq,lambda;
  float telDiameter = 64;
  float fwhm;
  float sumWeight;
  int iband=0;
  int nSpec=0;
  int nPos=0;
  FILE *fout;
  char outName[1024]="weightedSpectrum.dat";
  int nOff=0;

  float fLow  = 0;
  float fHigh = 5000;
  int   nc    = 0;
  int nSum=0;
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
      if (strcmp(argv[i],"-ra")==0)
	sscanf(argv[++i],"%f",&ra0);
      else if (strcmp(argv[i],"-dec")==0)
	sscanf(argv[++i],"%f",&dec0);
      else if (strcmp(argv[i],"-band")==0)
	sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outName,argv[++i]);
      else
	strcpy(fname[inFiles++],argv[i]);
    }

  printf("Number of files = %d\n",inFiles);

  filename = (char **)malloc(sizeof(char *)*MAX_SPECTRA);
  for (i=0;i<MAX_SPECTRA;i++)
    filename[i] = (char *)malloc(sizeof(char)*128);


  
  avSpectraX  = (float *)calloc(sizeof(float),MAX_CHANS);
  avSpectraY1 = (float *)calloc(sizeof(float),MAX_CHANS);
  avSpectraY2 = (float *)calloc(sizeof(float),MAX_CHANS);
  
  highResX    = (float **)malloc(sizeof(float *)*MAX_SPECTRA);
  highResY1   = (float **)malloc(sizeof(float *)*MAX_SPECTRA);
  highResY2   = (float **)malloc(sizeof(float *)*MAX_SPECTRA);

  for (i=0;i<MAX_SPECTRA;i++)
    {
      highResX[i]  = (float *)malloc(sizeof(float)*MAX_CHANS);
      highResY1[i] = (float *)malloc(sizeof(float)*MAX_CHANS);
      highResY2[i] = (float *)malloc(sizeof(float)*MAX_CHANS);      
    }
  plotX  = (float *)malloc(sizeof(float)*MAX_CHANS);
  plotX3 = (float *)malloc(sizeof(float)*MAX_CHANS);
  plotY1 = (float *)malloc(sizeof(float)*MAX_CHANS);
  plotY2 = (float *)malloc(sizeof(float)*MAX_CHANS);
  plotY3 = (float *)malloc(sizeof(float)*MAX_CHANS);


  printf("Complete allocation\n");

  // Load in all the spectra
  //  float ra[MAX_SPECTRA],dec[MAX_SPECTRA];
  for (i=0;i<inFiles;i++)
    {
      printf("Loading data from file: %s\n",fname[i]);
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_loadMetaData(inFile);

      
      sdhdf_loadBandData(inFile,ibeam,iband,1);
      // Process each dump
      for (l=0;l<inFile->beam[ibeam].bandHeader[iband].ndump;l++) 
	{
	  strcpy(filename[nSpec],fname[i]);
	  sdump[nSpec] = l;
	  f0[nSpec] = inFile->beam[ibeam].bandHeader[iband].fc;

	  ra[nSpec] = inFile->beam[ibeam].bandData[iband].astro_obsHeader[l].raDeg;
	  dec[nSpec] = inFile->beam[ibeam].bandData[iband].astro_obsHeader[l].decDeg;
	  nc=0;
	  for (k=0;k<inFile->beam[ibeam].bandHeader[iband].nchan;k++)
	    {
	      // FIX ME: using [0] for frequency dump
	      if (inFile->beam[ibeam].bandData[iband].astro_data.freq[k] > fLow && inFile->beam[ibeam].bandData[iband].astro_data.freq[k] <= fHigh)
		{
		  highResX[nSpec][nc]  = inFile->beam[ibeam].bandData[iband].astro_data.freq[k];
		  highResY1[nSpec][nc] = inFile->beam[ibeam].bandData[iband].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[iband].nchan+k];
		  highResY2[nSpec][nc] = inFile->beam[ibeam].bandData[iband].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[iband].nchan+k];
		  nc++;
		}
	    }
	  nChan[nSpec] = nc;

	  nSpec++;
	  if (nSpec == MAX_SPECTRA)
	    {
	      printf("ERROR: need to increase MAX_SPECTRA in sdhdf_weightSpectra.c\n");
	      exit(1);
	    }
	  
	}
      printf("Closing file\n");
      sdhdf_closeFile(inFile);
      printf("Next file\n");
    }
  printf("Data loaded\n");
  
  // Now form the weighted spectra

  // Currently assuming all bands are the same for each spectra
  freq = f0[0];
  lambda = 3.0e8/(freq*1e6);
  fwhm = 1.22*lambda*180/M_PI/telDiameter; 
  sigma = fwhm/2.35;
  printf("Sigma = %g, fwhm = %g, lambda = %g, fc = %g\n",sigma,fwhm,lambda,freq);
  
  sumWeight=0;
  for (k=0;k<nSpec;k++)
    {	     
      angle = haversine(ra[k],dec[k],ra0,dec0);      
      weight = exp(-(angle*angle)/2.0/sigma/sigma);
      printf("%s dump = %d k = %d, angle = %g, weight = %g ra/dec = (%g,%g), ra0/dec0 = (%g,%g)\n",filename[k],sdump[k],k,angle,weight,ra[k],dec[k],ra0,dec0);
      if (weight > 0.01)
	{
	  //		  printf("Angle = %g, weight = %g\n",angle,weight);
	  sumWeight+=weight;
	  for (l=0;l<nChan[k];l++)
	    {
	      avSpectraX[l]   = highResX[k][l]; // Should check if frequencies change etc. because of Doppler correction * NEEDS FIXING
	      avSpectraY1[l] += (weight*highResY1[k][l]);
	      avSpectraY2[l] += (weight*highResY2[k][l]);
	      // Should record the weightings somehow
	    }
	}
    }
  printf("Producing weighted spectra\n");
  if (sumWeight == 0)
    {
      printf("WARNING: sumWeight = 0\n");
      exit(1);
    }

  if (!(fout = fopen(outName,"w")))
    {
      printf("Unable to open file >%s<\n",outName);
      exit(1);
    }
  for (l=0;l<nChan[0];l++) // CHECK 0 HERE
    {
      avSpectraY1[i] /= sumWeight;
      avSpectraY2[i] /= sumWeight;
      fprintf(fout,"%.6f %g %g\n",avSpectraX[l],avSpectraY1[l],avSpectraY2[l]);
    }
  fclose(fout);
  
  for (i=0;i<MAX_SPECTRA;i++)
    {
      free(highResX[i]);
      free(highResY1[i]);
      free(highResY2[i]);
    }
  free(highResX); free(highResY1); free(highResY2);
  free(avSpectraX); free(avSpectraY1); free(avSpectraY2);
  free(plotX); free(plotY1); free(plotY2); free(plotY3);
  free(plotX3);

  for (i=0;i<MAX_SPECTRA;i++)
    free(filename[i]);
  free(filename);

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
