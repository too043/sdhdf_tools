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
// Software to provide a quick-look at the data within a SDHDF file.
//
// Usage:
// sdhdf_onoff -on <filename> -off <filename>
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"

#define VNUM "v2.0"
#define MAX_FILES_BATCH 1024

void help()
{
  printf("\nsdhdf_onoff       %s\n",VNUM);
  printf("INSPECTA version: %s\n",SOFTWARE_VER);
  printf("Author:           George Hobbs\n");
  printf("Software to produce (ON-OFF)/OFF spectra\n");

  printf("\n\nCommand line arguments:\n\n");
  printf("-h                This help\n");
  printf("-batch <filename> Read filenames from the file specified\n");
  printf("-on <filename>    SDHDF file corresponding to 'ON' source pointing\n");
  printf("-off <filename>   SDHDF file corresponding to 'OFF' source pointing\n");
  printf("-o <filename>     Output SDHDF file containing (ON-OFF)/OFF\n");
  printf("-subtract         Format ON-OFF, instead of (ON-OFF)/OFF\n");
  
  printf("\nExample:\n\n");
  printf("sdhdf_onoff -on on.hdf.T -off off.hdf.T -o diff.hdf\n");
  
  exit(1);
}

void determineOffSourceScale(float *freqOffScl,float *offSclA,float *offSclB,int n,float freq,float *sclA2,float *sclB2);

int main(int argc,char *argv[])
{
  int i,j,k,p,ii,b;
  char fnameOn[MAX_FILES_BATCH][MAX_STRLEN];
  char fnameOff[MAX_FILES_BATCH][MAX_STRLEN];
  char outFileName[MAX_FILES_BATCH][MAX_STRLEN];
  char defineOn[MAX_STRLEN];
  char defineOff[MAX_STRLEN];
  char defineOutFileName[MAX_STRLEN];
  char batchFileName[MAX_STRLEN];
  int batch=0;
  int nFiles=0;
  sdhdf_fileStruct *onFile;
  sdhdf_fileStruct *offFile;
  sdhdf_fileStruct *outFile;
  sdhdf_bandHeaderStruct *inBandParams;  
  //  spectralDumpStruct spectrumOn;
  //  spectralDumpStruct spectrumOff;
  float *out_data;
  float on_pol1,off_pol1;
  float on_pol2,off_pol2;
  float on_pol3,off_pol3;
  float on_pol4,off_pol4;
  float *freq;
  long nchan,npol,ndump;
  float sclA=1.0;
  float sclB=1.0;
  char scaleOff[128]="NULL";


  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;

  
  float freqOffScl[8192];
  float offSclA[8192];
  float offSclB[8192];
  int nOffScl=0;
  float sclA2=1.0,sclB2=1.0;

  int procType=1;
  
  FILE *fout;

  if (argc==1)
    help();

  printf("Starting\n");
  // help();

  if (!(onFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >onFile<\n");
      exit(1);
    }
  if (!(offFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >offFile<\n");
      exit(1);
    }
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }

  for (i=1;i<argc;i++)
    {       
      if (strcmp(argv[i],"-on")==0)
	strcpy(fnameOn[0],argv[++i]);	
      else if (strcmp(argv[i],"-h")==0)
	{
	  help();
	  exit(1);
	    
	}
      else if (strcmp(argv[i],"-off")==0)
	strcpy(fnameOff[0],argv[++i]);	
      else if (strcmp(argv[i],"-sclA")==0)
	sscanf(argv[++i],"%f",&sclA);
      else if (strcmp(argv[i],"-sclB")==0)
	sscanf(argv[++i],"%f",&sclB);
      else if (strcmp(argv[i],"-subtract")==0)
	procType=2;
      else if (strcmp(argv[i],"-1div")==0)
	procType=3;
      else if (strcmp(argv[i],"-div")==0)
	procType=4;
      else if (strcmp(argv[i],"-scaleOff")==0)
	strcpy(scaleOff,argv[++i]);
      else if (strcmp(argv[i],"-h")==0) // Should actually free memory
	exit(1);
      else if (strcmp(argv[i],"-batch")==0)
	{
	  batch=1;
	  strcpy(batchFileName,argv[++i]);
	}
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outFileName[0],argv[++i]);
    }

  if (strcmp(scaleOff,"NULL")!=0)
    {
      FILE *fin;
      fin = fopen(scaleOff,"r");
      nOffScl=0;
      while (!feof(fin))
	{
	  if (fscanf(fin,"%f %f %f",&freqOffScl[nOffScl],&offSclA[nOffScl],&offSclB[nOffScl])==3)
	    nOffScl++;
	}
      fclose(fin);
      printf("Loaded %d values for the off source system noise/temperature\n",nOffScl);
    }

  
  if (batch==0)
    nFiles=1;
  else
    {
      char line[4096];
      FILE *fin;
      nFiles=0;
      fin = fopen(batchFileName,"r");
      while (!feof(fin))
	{
	  if (fgets(line,4096,fin)!=NULL)
	    //	  if (fscanf(fin,"%s %s %s",fnameOn[nFiles],fnameOff[nFiles],outFileName[nFiles])==3)
	    {
	      sscanf(line,"%s %s %s",fnameOn[nFiles],fnameOff[nFiles],outFileName[nFiles]);
	      nFiles++;
	      if (nFiles == MAX_FILES_BATCH)
		{
		  printf("ERROR: Need to increase MAX_FILES_BATCH in sdhdf_onoff.c\n");
		  printf("Currently the maximum number of files is set to %d\n",MAX_FILES_BATCH);
		  exit(1);
		}
	    }
	}
      fclose(fin);
      printf("Loaded %d entires from %s\n",nFiles,batchFileName);
    }
  
  for (ii=0;ii<nFiles;ii++)
    {
      printf("Processing files: %s %s %s\n",fnameOn[ii],fnameOff[ii],outFileName[ii]);
      sdhdf_initialiseFile(onFile);
      sdhdf_initialiseFile(offFile);
      sdhdf_initialiseFile(outFile);
      sdhdf_openFile(outFileName[ii],outFile,3);         
      sdhdf_openFile(fnameOn[ii],onFile,1);
      sdhdf_openFile(fnameOff[ii],offFile,1);
      sdhdf_loadMetaData(onFile);
      sdhdf_loadMetaData(offFile);
      // Should do checks for: (NOT DOING YET)
      // 1) ASSUMING NDUMP = 1 --- NEED TO FIX THIS!!!
      // 2) number of channels same between files
      
      //  inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*onFile->nBands);      
      //  sdhdf_loadBandMetadataStruct(onFile,inBandParams);
      for (b=0;b<onFile->nBeam;b++)
	{
	  for (i=0;i<onFile->beam[b].nBand;i++)
	    {
	      printf("Processing band: %d\n",i);
	      nchan = onFile->beam[b].bandHeader[i].nchan;
	      npol  = onFile->beam[b].bandHeader[i].npol;
	      ndump = 1; //onFile->beam[b].bandHeader[i].ndump;
	      out_data = (float *)malloc(sizeof(float)*nchan*npol*ndump);
	      freq = (float *)malloc(sizeof(float)*nchan);      
	      sdhdf_loadBandData(onFile,b,i,1);
	      sdhdf_loadBandData(offFile,b,i,1);
	      /* FIX ME
		 sdhdf_loadFrequencyAttributes(onFile,i);
		 sdhdf_loadDataAttributes(onFile,i);
		 strcpy(outFile->frequency_attr.frame,onFile->frequency_attr.frame);
		 strcpy(outFile->frequency_attr.unit,onFile->frequency_attr.unit);
		 strcpy(outFile->data_attr.unit,onFile->data_attr.unit);
	      */
	      for (k=0;k<nchan;k++)		
		freq[k] = onFile->beam[b].bandData[i].astro_data.freq[k];  // FIX ME: DUMP FOR FREQ AXIS
	      for (k=0;k<nchan;k++)
		{
		  // 0 here should be the spectral dump number
		  if (npol==1)	   
		    {
		      on_pol1  = onFile->beam[b].bandData[i].astro_data.pol1[k];
		      off_pol1 = offFile->beam[b].bandData[i].astro_data.pol1[k];

		      if (procType==1)			
			out_data[k+0*npol*nchan]         = sclA*(on_pol1-off_pol1)/off_pol1;
		      else if (procType==2)
			out_data[k+0*npol*nchan]         = sclA*(on_pol1-off_pol1);
		      else if (procType==3)
			out_data[k+0*npol*nchan]         = sclA*on_pol1/(off_pol1-on_pol1);
		      else
			out_data[k+0*npol*nchan]         = sclA*on_pol1/(off_pol1);
			
		    }
		  else if (npol == 2)
		    {
		      on_pol1 = onFile->beam[b].bandData[i].astro_data.pol1[k];
		      on_pol2 = onFile->beam[b].bandData[i].astro_data.pol2[k];
		      off_pol1 = offFile->beam[b].bandData[i].astro_data.pol1[k];
		      off_pol2 = offFile->beam[b].bandData[i].astro_data.pol2[k];

		      if (procType==1)
			{
			  out_data[k+0*npol*nchan]         = sclA*(on_pol1-off_pol1)/off_pol1;
			  out_data[k+0*npol*nchan+nchan]   = sclB*(on_pol2-off_pol2)/off_pol2;
			}
		      else if (procType==2)
			{
			  out_data[k+0*npol*nchan]         = sclA*(on_pol1-off_pol1);
			  out_data[k+0*npol*nchan+nchan]   = sclB*(on_pol2-off_pol2);
			}
		      else if (procType==3)
			{
			  out_data[k+0*npol*nchan]         = sclA*on_pol1/(off_pol1-on_pol1);
			  out_data[k+0*npol*nchan+nchan]   = sclB*on_pol2/(off_pol2-on_pol2);
			}
		      else
			{
			  out_data[k+0*npol*nchan]         = sclA*on_pol1/(off_pol1);
			  out_data[k+0*npol*nchan+nchan]   = sclB*on_pol2/(off_pol2);
			}
			
		      
		    }
		  else
		    {
		      on_pol1 = onFile->beam[b].bandData[i].astro_data.pol1[k];
		      on_pol2 = onFile->beam[b].bandData[i].astro_data.pol2[k];
		      on_pol3 = onFile->beam[b].bandData[i].astro_data.pol3[k];
		      on_pol4 = onFile->beam[b].bandData[i].astro_data.pol4[k];

		      off_pol1 = offFile->beam[b].bandData[i].astro_data.pol1[k];
		      off_pol2 = offFile->beam[b].bandData[i].astro_data.pol2[k];
		      off_pol3 = offFile->beam[b].bandData[i].astro_data.pol3[k];
		      off_pol4 = offFile->beam[b].bandData[i].astro_data.pol4[k];
		      sclA2=sclB2 = 1.0;
		      if (nOffScl>0)
			determineOffSourceScale(freqOffScl,offSclA,offSclB,nOffScl,freq[k],&sclA2,&sclB2);
		      
		      if (procType==1)
			{
			  out_data[k+0*npol*nchan]         = sclA2*sclA*(on_pol1-off_pol1)/off_pol1;
			  out_data[k+0*npol*nchan+nchan]   = sclB2*sclB*(on_pol2-off_pol2)/off_pol2;
			  out_data[k+0*npol*nchan+2*nchan] = (on_pol3-off_pol3)/off_pol3;
			  out_data[k+0*npol*nchan+3*nchan] = (on_pol4-off_pol4)/off_pol4;
			}
		      else if (procType==2)
   			{
			  out_data[k+0*npol*nchan]         = sclA2*sclA*(on_pol1-off_pol1);
			  out_data[k+0*npol*nchan+nchan]   = sclB2*sclB*(on_pol2-off_pol2);
			  out_data[k+0*npol*nchan+2*nchan] = (on_pol3-off_pol3);
			  out_data[k+0*npol*nchan+3*nchan] = (on_pol4-off_pol4);
			}
		      else if (procType==3)
			{
			  out_data[k+0*npol*nchan]         = sclA2*sclA*on_pol1/(off_pol1-on_pol1);
			  out_data[k+0*npol*nchan+nchan]   = sclB2*sclB*on_pol2/(off_pol2-on_pol2);
			  out_data[k+0*npol*nchan+2*nchan] = on_pol3/(off_pol3-on_pol3);
			  out_data[k+0*npol*nchan+3*nchan] = on_pol4/(off_pol4-on_pol4);
			}
		      else
			{
			  out_data[k+0*npol*nchan]         = sclA2*sclA*on_pol1/(off_pol1);
			  out_data[k+0*npol*nchan+nchan]   = sclB2*sclB*on_pol2/(off_pol2);
			  out_data[k+0*npol*nchan+2*nchan] = on_pol3/(off_pol3);
			  out_data[k+0*npol*nchan+3*nchan] = on_pol4/(off_pol4);
			}
		    }
		      //	  printf("DIFF =  %g %g %g %g\n",freq[k],on_pol1,off_pol1,(on_pol1-off_pol1)/off_pol1);
	      
	      // 0 should be replaced by dump number
		}

	      // Note that this is getting the attributes from the ON file -- really should check if they are the same as in the OFF file
	      sdhdf_copyAttributes(onFile->beam[b].bandData[i].astro_obsHeaderAttr,onFile->beam[b].bandData[i].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
	      sdhdf_copyAttributes(onFile->beam[b].bandData[i].astro_obsHeaderAttr_freq,onFile->beam[b].bandData[i].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);	  

	      sdhdf_releaseBandData(onFile,b,i,1);
	      sdhdf_releaseBandData(offFile,b,i,1);

	      // FIX ME: Only sending 1 frequency channel through
	      sdhdf_writeSpectrumData(outFile,onFile->beamHeader[b].label,onFile->beam[b].bandHeader[i].label,b,i,out_data,freq,1,nchan,1,npol,ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	      free(out_data);
	      free(freq);
	    }
	}
      sdhdf_writeHistory(outFile,onFile->history,onFile->nHistory);  

      //
      // Need to copy time stamps etc. from the ON source data
      //
      sdhdf_copyRemainder(onFile,outFile,0);
      //  free(inBandParams);
      
      
      sdhdf_closeFile(onFile);
      sdhdf_closeFile(offFile);
      sdhdf_closeFile(outFile);
    }	
    
  free(onFile);
  free(offFile);
  free(outFile);



}

void determineOffSourceScale(float *freqOffScl,float *offSclA,float *offSclB,int n,float freq,float *sclA2,float *sclB2)
{
  int i;

  for (i=0;i<n-1;i++)
    {
      if (freq <= freqOffScl[i])
	{
	  *sclA2 = offSclA[i]; // Should do an interpolation
	  *sclB2 = offSclB[i];
	  break;
	}
    }
}
