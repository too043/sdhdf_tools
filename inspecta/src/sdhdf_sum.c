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
// Software to sum spectra together.
//
// Usage:
// sdhdf_sum -o <outfilename> ... list of input filenames ...
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"

#define VNUM "v0.1"

void help()
{
  printf("sdhdf_sum:        %s\n",VNUM);
  printf("sdhfProc version: %s\n",SOFTWARE_VER);
  printf("author:           George Hobbs\n");
  printf("\n");
  printf("Software to sum spectra together\n");
  printf("\n\nCommand line arguments:\n\n");
  printf("-divN             Divide by the number of summed observations\n");
  printf("-o <file>         Output filename\n");
  printf("-h                This help\n");
  

  printf("---------------------\n");
}

int main(int argc,char *argv[])
{
  int i,j,k,p,b;
  int timeThrough=1;
  sdhdf_fileStruct *inFile;
  sdhdf_fileStruct *file0;
  sdhdf_fileStruct *outFile;
  //  spectralDumpStruct spectrum;
  float *out_data;
  float in_pol1;
  float in_pol2;
  float in_pol3;
  float in_pol4;
  double tav=0;
  sdhdf_bandHeaderStruct *inBandParams;  
  long nchan,npol,ndump;
  int divN=0;

  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;

  
  char outFileName[MAX_STRLEN];
  int nFiles=0;
  char **fname;
  FILE *fout;
  double scaleFactor=1;

  printf("Starting sdhdf_sum\n");
  //  help();
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >onFile<\n");
      exit(1);
    }
  if (!(file0 = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >file0<\n");
      exit(1);
    }
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);
  sdhdf_initialiseFile(file0);
  sdhdf_initialiseFile(outFile);

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outFileName,argv[++i]);
      else if (strcmp(argv[i],"-divN")==0) // NEED TO UPDATE COMMAND LINE ARGUMENT LIST BELOW AS WELL -- SHOULD FIX
	divN=1;
      else if (strcmp(argv[i],"-h")==0) // Should actually free memory
	{help(); exit(1);}
      else
	nFiles++;
    }
  printf("Loaded command line arguments\n");
  fname = (char **)malloc(sizeof(char *)*nFiles);
  for (i=0;i<nFiles;i++)
    fname[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);

  nFiles=0;
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	i++;
      else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-divN")==0)
	{} // Do nothing
      else
	strcpy(fname[nFiles++],argv[i]);
    }
  printf("Have %d input files\n",nFiles);
  // Get information from the first file

  sdhdf_openFile(fname[0],file0,1);
  sdhdf_loadMetaData(file0);
  printf("Number of beams = %d\n",file0->nBeam);
  // Open the output file
  sdhdf_openFile(outFileName,outFile,3);
  printf("File opened for output >%s<\n",outFileName);
  
  for (b=0;b<file0->nBeam;b++)
    {
      inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*file0->beam[b].nBand);      
      sdhdf_copyBandHeaderStruct(file0->beam[b].bandHeader,inBandParams,file0->beam[b].nBand);
                
      for (i=0;i<file0->beam[b].nBand;i++)
	{
	  tav=0;
	  sdhdf_loadBandData(file0,b,i,1);
	  
	  printf("Processing band: %s\n",file0->beam[b].bandHeader[i].label);
	  nchan = file0->beam[b].bandHeader[i].nchan;
	  npol  = file0->beam[b].bandHeader[i].npol;
	  
	  // GH: FIX ME
	  //      sdhdf_loadFrequencyAttributes(file0,i);
	  //      sdhdf_loadDataAttributes(file0,i);
	  //      strcpy(outFile->frequency_attr.frame,file0->frequency_attr.frame);
	  //      strcpy(outFile->frequency_attr.unit,file0->frequency_attr.unit);
	  //      strcpy(outFile->data_attr.unit,file0->data_attr.unit);
	  	  
	  ndump = 1;
	  out_data = (float *)calloc(sizeof(float),nchan*npol*ndump);      
	  
	  // OPEN AND CLOSE ALL THE FILES LOTS OF TIMES -- CAN IMPROVE THE CODING HERE
	  for (j=0;j<nFiles;j++)
	    {
	      printf(" ... %s\n",fname[j]);
	      sdhdf_initialiseFile(inFile);	  
	      printf("Opening\n");
	      sdhdf_openFile(fname[j],inFile,1);
	      printf("Complete opening\n");
	      sdhdf_loadMetaData(inFile);
	      printf("Complete loading metadata\n");
	      
	      nchan = inFile->beam[b].bandHeader[i].nchan;
	      npol  = inFile->beam[b].bandHeader[i].npol;
	      ndump = 1;
	      printf("nchan = %d, npol = %d\n",(int)nchan,(int)npol);
	      
	      sdhdf_loadBandData(inFile,b,i,1);
	      tav += inFile->beam[b].bandHeader[i].dtime;
	      printf("Loaded the data\n");
	      scaleFactor=1;
	      
	      for (k=0;k<nchan;k++)
		{
		  // 0 here should be the spectral dump number
		  if (npol == 1)
		    {
		      in_pol1 = inFile->beam[b].bandData[i].astro_data.pol1[k];
		      if (divN==1)
			out_data[k+0*npol*nchan]         += scaleFactor*(in_pol1)/nFiles;
		      else
			out_data[k+0*npol*nchan]         += scaleFactor*(in_pol1);
		    }
		  else if (npol==2)
		    {
		      in_pol1 = inFile->beam[b].bandData[i].astro_data.pol1[k];
		      in_pol2 = inFile->beam[b].bandData[i].astro_data.pol2[k];
		      if (divN==1)
			{
			  out_data[k+0*npol*nchan]         += scaleFactor*(in_pol1)/nFiles;
			  out_data[k+0*npol*nchan+nchan]   += scaleFactor*(in_pol2)/nFiles;		      
			}
		      else
			{
			  out_data[k+0*npol*nchan]         += scaleFactor*(in_pol1);
			  out_data[k+0*npol*nchan+nchan]   += scaleFactor*(in_pol2);		      
			}
		    }
		  else
		    {
		      in_pol1 = inFile->beam[b].bandData[i].astro_data.pol1[k];
		      in_pol2 = inFile->beam[b].bandData[i].astro_data.pol2[k];
		      in_pol3 = inFile->beam[b].bandData[i].astro_data.pol3[k];
		      in_pol4 = inFile->beam[b].bandData[i].astro_data.pol4[k];
		      if (divN==1)
			{
			  out_data[k+0*npol*nchan]         += scaleFactor*(in_pol1)/nFiles;
			  out_data[k+0*npol*nchan+nchan]   += scaleFactor*(in_pol2)/nFiles;
			  out_data[k+0*npol*nchan+2*nchan] += scaleFactor*(in_pol3)/nFiles;
			  out_data[k+0*npol*nchan+3*nchan] += scaleFactor*(in_pol4)/nFiles;		      
			}
		      else
			{
			  out_data[k+0*npol*nchan]         += scaleFactor*(in_pol1);
			  out_data[k+0*npol*nchan+nchan]   += scaleFactor*(in_pol2);
			  out_data[k+0*npol*nchan+2*nchan] += scaleFactor*(in_pol3);
			  out_data[k+0*npol*nchan+3*nchan] += scaleFactor*(in_pol4);		      
			}
		    }
		}
	  
	      if (j==0) // Note that this is getting the attributes from the first file -- really should check if they are the same in all the files
		{
		  
		  // Need to do this separately for each band *** FIX ME **** CURRENTLY IT WILL ONLY TAKE THE LAST ONE
		  //		  sdhdf_copyAttributes(inFile->beam[b].bandData[i].astro_obsHeaderAttr,inFile->beam[b].bandData[i].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
		  //		  sdhdf_copyAttributes(inFile->beam[b].bandData[i].astro_obsHeaderAttr_freq,inFile->beam[b].bandData[i].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);	  
		}

	      printf("Releasing data\n");
	      sdhdf_releaseBandData(inFile,b,i,1);
	      printf("Closing file\n");
	      sdhdf_closeFile(inFile);
	      printf("File closed\n");
	    }
	  printf("Writing the output\n");
	  inBandParams[i].ndump=1;
	  inBandParams[i].dtime = tav;
	  // NOTE THAT THE ATTRIBUTES ARE INCORRECT HERE -- SEE COMMENT ABOVE ABOUT ONLY TAKING THE LAST ONE ***** FIX ME

	  sdhdf_writeSpectrumData(outFile,file0->beamHeader[b].label,file0->beam[b].bandHeader[i].label,b,i,out_data,file0->beam[b].bandData[i].astro_data.freq,file0->beam[b].bandData[i].astro_data.nFreqDumps,nchan,1,npol,ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	  // GH: FIX ME
	  //      sdhdf_writeFrequencyAttributes(outFile,file0->bandHeader[i].label);
	  //      sdhdf_writeDataAttributes(outFile,file0->bandHeader[i].label);

	  free(out_data);
	}
      sdhdf_writeBandHeader(outFile,inBandParams,file0->beamHeader[b].label,file0->beam[b].nBand,1,inFile->primary[0].hdr_defn_version);
      free(inBandParams);
    }
  sdhdf_writeHistory(outFile,file0->history,file0->nHistory);  
  sdhdf_copyRemainder(file0,outFile,0);
  
  sdhdf_closeFile(outFile);
  sdhdf_closeFile(file0);

  free(inFile);
  free(outFile);
  free(file0);

  for (i=0;i<nFiles;i++)
    free(fname[i]);
  free(fname);

  printf("Complete: have written file %s\n",outFileName);
}

