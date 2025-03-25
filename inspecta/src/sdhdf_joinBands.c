//
// This joins a single band over different time chunks together into one file
// -- joining in time
//

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
// Software to join files together
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"

#define MAX_BANDS 26     // FIX THIS

//
// Note that the spectral dump number may be different for different subbands
// Should use the band extract first
//

void help()
{
  printf("sdhdf_join: routine to join sdhdf files together\n");
  
}

int main(int argc,char *argv[])
{
  int ii,i,j,k,b;
  sdhdf_fileStruct *inFile,*inFile0,*outFile;
  sdhdf_bandHeaderStruct *outBandParams,*outCalBandParams,*inBandParams,*inCalBandParams;
  sdhdf_obsParamsStruct  *outObsParams,*outCalObsParams,*inObsParams,*inCalObsParams;
  
  char fname[MAX_FILES][128];
  int nFiles=0;
  char oname[MAX_STRLEN]="join.hdf";
  int nBeam;
  int nDumpIn,nDumpOut,ndump,nFreqDumpOut;
  int nFreqDumps;
  int dumpPos,freqPos,dataPos;
  int nOutBands;
  float *out_data,*out_freq,*out_wts;
  unsigned char *out_flags;
  int haveWeights,haveFlags=0;
  
  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;
  char args[MAX_ARGLEN]="";
  
  sdhdf_storeArguments(args,MAX_ARGLEN,argc,argv);
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(oname,argv[++i]);
      else
	strcpy(fname[nFiles++],argv[i]);
    }

  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

    if (!(inFile0 = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(outFile);
  sdhdf_openFile(oname,outFile,3);
  
  
  sdhdf_openFile(fname[0],inFile0,1); 
  sdhdf_loadMetaData(inFile0);
  nBeam = inFile0->nBeam;
  sdhdf_allocateBeamMemory(outFile,nBeam);  

  
  for (b=0;b<nBeam;b++)
    {
      printf("Number of bands = %d\n",inFile0->beam[b].nBand);
      nOutBands = inFile0->beam[b].nBand;
      outBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nOutBands);
      sdhdf_copyBandHeaderStruct(inFile0->beam[b].bandHeader,outBandParams,inFile0->beam[b].nBand);  

      for (i=0;i<inFile0->beam[b].nBand;i++)
	{
	  sdhdf_loadDataFreqAttributes(inFile0,b,i,dataAttributes,&nDataAttributes,freqAttributes,&nFreqAttributes);
	  nDumpIn = inFile0->beam[b].bandHeader[i].ndump;

	  nDumpOut = 0;
	  dumpPos = freqPos = dataPos = 0;

	  for (ii=0;ii<nFiles;ii++)
	    {
	      sdhdf_openFile(fname[ii],inFile,1); 
	      sdhdf_loadMetaData(inFile); // Could do more checks here to check that we're properly combining everything
	      nDumpOut     += inFile->beam[b].bandHeader[i].ndump;
	      sdhdf_closeFile(inFile);
	    }
	  printf("beam %d, band %d, number of output spectral dumps = %d\n",b,i,nDumpOut);
	  //
	  outBandParams[i].ndump = nDumpOut; // FIX ME: SHOULD CHECK ALL THE DUMP TIMES ARE THE SAME
	  out_data = (float *)malloc(sizeof(float)*nDumpOut*inFile0->beam[b].bandHeader[i].nchan * inFile0->beam[b].bandHeader[i].npol);
	  out_freq = (float *)malloc(sizeof(float)*nDumpOut*inFile0->beam[b].bandHeader[i].nchan);
	  out_wts = (float *)malloc(sizeof(float)*nDumpOut*inFile0->beam[b].bandHeader[i].nchan);
	  out_flags = (unsigned char *)malloc(sizeof(unsigned char)*nDumpOut*inFile0->beam[b].bandHeader[i].nchan);
	  
	  outObsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*nDumpOut);
	  for (ii=0;ii<nFiles;ii++)
	    {
	      sdhdf_openFile(fname[ii],inFile,1); 
	      sdhdf_loadMetaData(inFile); // Could do more checks here to check that we're properly combining everything

	      sdhdf_loadBandData2Array(inFile,b,i,1,out_data+dataPos);
	      sdhdf_loadFrequency2Array(inFile,b,i,out_freq+freqPos,&nFreqDumps);
	      haveWeights = sdhdf_loadWeights2Array(inFile,b,i,out_wts+freqPos);
	      haveFlags = sdhdf_loadFlags2Array(inFile,b,i,out_flags+freqPos);


	      ndump = inFile->beam[b].bandHeader[i].ndump;

	      if (nFreqDumps != ndump)
		{
		  if (nFreqDumps == 1)
		    {
		      for (j=0;j<ndump;j++) // Make multiple copies
			memcpy(out_freq+freqPos+j*inFile->beam[b].bandHeader[i].nchan,out_freq+freqPos,sizeof(float)*inFile->beam[b].bandHeader[i].nchan);
		    }
		  else
		    {printf("ERROR: Not sure what to do\n"); exit(1);}
		}
	      dataPos += (inFile->beam[b].bandHeader[i].nchan * inFile->beam[b].bandHeader[i].npol * inFile->beam[b].bandHeader[i].ndump);
	      freqPos += (inFile->beam[b].bandHeader[i].nchan * inFile->beam[b].bandHeader[i].ndump);
	      //	      freqDumpPos += inFile->beam[b].bandData[i].astro_data.nFreqDumps;

	      for (k=0;k<ndump;k++)
		sdhdf_copySingleObsParams(inFile,b,i,k,outObsParams+k+dumpPos);
	      dumpPos += ndump;
	      sdhdf_closeFile(inFile);
	    }
	  sdhdf_writeSpectrumData(outFile,inFile0->beamHeader[b].label,inFile0->beam[b].bandHeader[i].label,b,i,out_data,
				  out_freq,nDumpOut,inFile0->beam[b].bandHeader[i].nchan,1,inFile0->beam[b].bandHeader[i].npol,
				  nDumpOut,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	  if (haveWeights==1)
	    sdhdf_writeDataWeights(outFile,b,i,out_wts,inFile0->beam[b].bandHeader[i].nchan,nDumpOut,inFile0->beamHeader[b].label,outBandParams[i].label);
	  if (haveFlags==1)
	    sdhdf_writeFlags(outFile,b,i,out_flags,inFile0->beam[b].bandHeader[i].nchan,nDumpOut,inFile->beamHeader[b].label,outBandParams[i].label);

	  sdhdf_writeObsParams(outFile,outBandParams[i].label,inFile0->beamHeader[b].label,i,outObsParams,nDumpOut,1,inFile->primary[0].hdr_defn_version);
	  free(outObsParams);	  
	}
      sdhdf_writeBandHeader(outFile,outBandParams,inFile0->beamHeader[b].label,nOutBands,1,inFile->primary[0].hdr_defn_version);
      free(outBandParams);
    }
  sdhdf_addHistory(inFile0->history,inFile0->nHistory,"sdhdf_join","INSPECTA software to combine files",args);
  printf("args = >%s<\n",args);
  inFile0->nHistory++;
  sdhdf_writeHistory(outFile,inFile0->history,inFile0->nHistory);
  sdhdf_copyRemainder(inFile0,outFile,0);
  
  sdhdf_closeFile(inFile0);

  sdhdf_closeFile(outFile);
}


