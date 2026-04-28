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
// Software to measure, analyse and remove baselines from SDHDF data sets
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
int help()
{
  printf("-e <str>        Set the file extension\n");
  printf("-nfit <integer> Number of polynomial coefficients when fitting the baseline\n");
  printf("-excludeFreq <f0> <f1> When fitting the baseline ignore frequencies between f0 and f1\n");
  exit(1);
}


int main(int argc,char *argv[])
{
  int        i,j,k,l,ii,kk;
  char       fname[MAX_FILES][64];
  sdhdf_fileStruct *inFile,*outFile;
  int inFiles;
  int        idump,iband,ibeam;
  double *baselineX,*baselineY1,*baselineY2;
  float *outArr;
  float freq;
  int nBaseline;
  int nfit = 1;
  double params1[16],params2[16],v[16];
  float excludeFreq0,excludeFreq1;
  double fc;
  char oname[1024]="baseline.hdf";
  float val1,val2;
  int nchan,npol,ndump;
  char ext[1024]="baseline";  
  // Defaults
  idump = iband = ibeam = 0;
  excludeFreq0=excludeFreq1 = -1;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
  sdhdf_initialiseFile(inFile);
    if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for the output file\n");
      exit(1);
    }
  sdhdf_initialiseFile(outFile);
  
  inFiles=0;
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-e")==0)
	strcpy(ext,argv[++i]);    
      else if (strcmp(argv[i],"-nfit")==0)
	sscanf(argv[++i],"%d",&nfit);
      else if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-excludeFreq")==0)
	{
	  sscanf(argv[++i],"%f",&excludeFreq0);
	  sscanf(argv[++i],"%f",&excludeFreq1);
	}
      else
	strcpy(fname[inFiles++],argv[i]);
    }
  for (i=0;i<inFiles;i++)
    {
      printf("Processing file: %s\n",fname[i]);
      sdhdf_initialiseFile(inFile);
      sdhdf_initialiseFile(outFile);
      sdhdf_formOutputFilename(fname[i],ext,oname);
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_openFile(oname,outFile,3);

      sdhdf_loadMetaData(inFile);
     for (j=0;j<inFile->beam[ibeam].nBand;j++)
       sdhdf_loadBandData(inFile,ibeam,j,1);
           
     sdhdf_copyRemainder(inFile,outFile,0);
     // Now do the modelling
     // Process each band
      for (j=0;j<inFile->beam[ibeam].nBand;j++)
      {
	  printf("Processing band: %d\n",j);
	  
	  // Process each dump
	  for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l++) 
	    {
	      baselineX = (double *)malloc(sizeof(double)*inFile->beam[ibeam].bandHeader[j].nchan);
	      baselineY1 = (double *)malloc(sizeof(double)*inFile->beam[ibeam].bandHeader[j].nchan);
	      baselineY2 = (double *)malloc(sizeof(double)*inFile->beam[ibeam].bandHeader[j].nchan);
	      nBaseline=0;
	      fc = 0.5*(inFile->beam[ibeam].bandData[j].astro_data.freq[0]+inFile->beam[ibeam].bandData[j].astro_data.freq[inFile->beam[ibeam].bandHeader[j].nchan-1]);
	      for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		{
		  val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
		  val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];		
		  freq = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
		  if (excludeFreq0 < freq && excludeFreq1 >= freq)
		    {
		      // Do nothing
		    }
		  else
		    {
		      baselineX[nBaseline] = freq-fc;
		      baselineY1[nBaseline] = val1;
		      baselineY2[nBaseline] = val2;
		      nBaseline++;
		    }
		}
	      printf("Number of channels used to form the baseline: %d\n",nBaseline);
	      
	      TKleastSquares_svd_noErr(baselineX,baselineY1,nBaseline,params1,nfit,TKfitPoly);
	      TKleastSquares_svd_noErr(baselineX,baselineY2,nBaseline,params2,nfit,TKfitPoly);
	      printf("Baseline parameters for polarisation 1 and 2 are \n");
	      for (k=0;k<nfit;k++)
		printf(" [%d] %g %g\n",k, params1[k],params2[k]);
	      
	    // Remove the baseline
	    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
	      {
		freq = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
		TKfitPoly(freq-fc,v,nfit);
		for (kk=0;kk<nfit;kk++)
		  {		   
		    inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k]-=v[kk]*params1[kk];
		    inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k]-=v[kk]*params2[kk];		
		  }
		//		printf("output: %.6f %g %g\n",freq,inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k],inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k]);
	      }	    
	    	    
	    free(baselineX);
	    free(baselineY1);
	    free(baselineY2);
	    }
	  printf("Replacing %s\n",inFile->beamHeader[ibeam].label);
	  
	  // Get a single array
	  nchan = inFile->beam[ibeam].bandHeader[j].nchan;
	  npol = 4 ; // FIX ME
	  ndump = inFile->beam[ibeam].bandHeader[j].ndump;
	  
	  outArr = (float *)malloc(sizeof(float)*nchan*npol*ndump);
	  for (k=0;k<ndump;k++)
	    {
	      for (kk=0;kk<nchan;kk++)
		{
		  outArr[k*npol*nchan + 0*nchan + kk] = inFile->beam[ibeam].bandData[j].astro_data.pol1[kk];
		  outArr[k*npol*nchan + 1*nchan + kk] = inFile->beam[ibeam].bandData[j].astro_data.pol2[kk];
		  outArr[k*npol*nchan + 2*nchan + kk] = inFile->beam[ibeam].bandData[j].astro_data.pol3[kk];
		  outArr[k*npol*nchan + 3*nchan + kk] = inFile->beam[ibeam].bandData[j].astro_data.pol4[kk];
		}
	    }
	  sdhdf_replaceSpectrumData(outFile,inFile->beam[ibeam].bandHeader[j].label,inFile->beamHeader[ibeam].label,j,outArr,ndump,npol,nchan);
	  free(outArr);
      }
    



      sdhdf_closeFile(inFile);
      sdhdf_closeFile(outFile);
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
