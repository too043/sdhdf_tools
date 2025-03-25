//
// Combines hdf files into a single file
// This will only change metadata - no data table will be updated
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


void help()
{
  printf("sdhdf_joinFiles: routine to join sdhdf files together\n");
  
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
  sdhdf_beamHeaderStruct *beamHeader;
  sdhdf_primaryHeaderStruct *primaryHeader;
  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;
  char args[MAX_ARGLEN]="";
  int bandParamsPos=0;
  char newBandName[1024];
  char groupName1[1024];
  char groupName2[1024];
  char tempStr[1024];
  
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
  sdhdf_allocateBeamMemory(outFile,nBeam);   // Assuming the same number of beams in each file

  beamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nBeam);  
  primaryHeader = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct));
  strcpy(beamHeader[0].label,inFile0->beamHeader[0].label);
  strcpy(beamHeader[0].source,inFile0->beamHeader[0].source);
  primaryHeader = inFile0->primary;
  
  for (b=0;b<nBeam;b++)
    {

      nOutBands = 0;

      for (ii=0;ii<nFiles;ii++)
	{
	  sdhdf_openFile(fname[ii],inFile,1); 
	  sdhdf_loadMetaData(inFile); // Could do more checks here to check that we're properly combining everything
	  nOutBands  += inFile->beam[b].nBand;
	  sdhdf_closeFile(inFile);
	}
      beamHeader[b].nBand = nOutBands;
      printf("Number of output bands = %d\n",nOutBands);
      
      outBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nOutBands);
      outCalBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nOutBands);

      bandParamsPos=0;
      for (ii=0;ii<nFiles;ii++)
	{
	  printf("Opening file %s\n",fname[ii]);
	  sdhdf_openFile(fname[ii],inFile,1);
	  sdhdf_loadMetaData(inFile);
	  sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,outBandParams+bandParamsPos,inFile->beam[b].nBand);

	  if (strcmp(inFile->primary[0].cal_mode,"ON")==0)
	    sdhdf_copyBandHeaderStruct(inFile->beam[b].calBandHeader,outCalBandParams+bandParamsPos,inFile->beam[b].nBand);
	  
	  // Replace band names
	  sprintf(newBandName,"band_%03d",ii); // NEEDS FIXING -- NEED TO USE ALL BANDS IN THAT FILE

	  strcpy(outBandParams[bandParamsPos].label,newBandName);
	  if (strcmp(inFile->primary[0].cal_mode,"ON")==0)
	    strcpy(outCalBandParams[bandParamsPos].label,newBandName);
	 
	  bandParamsPos += inFile->beam[b].nBand;
	  printf("Copying in file %d\n",ii);
	  sprintf(groupName1,"%s/%s/astronomy_data",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[0].label); // FIX FOR MORE THAN ONE BAND
	  sprintf(tempStr,"%s",inFile->beamHeader[b].label);
	  if (sdhdf_checkGroupExists(outFile,tempStr) == 1)
	    {
	      hid_t group_id;
	      herr_t status=0;
	      
	      group_id = H5Gcreate2(outFile->fileID,tempStr,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	      status   = H5Gclose(group_id);
	    }
	  sprintf(tempStr,"%s/%s",inFile->beamHeader[b].label,newBandName);
	  if (sdhdf_checkGroupExists(outFile,tempStr) == 1)
	    {
	      hid_t group_id;
	      herr_t status=0;
	      group_id = H5Gcreate2(outFile->fileID,tempStr,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	      status   = H5Gclose(group_id);
	    }

	  sprintf(groupName2,"%s/%s/astronomy_data",inFile->beamHeader[b].label,newBandName); // FIX FOR MORE THAN ONE BAND	  
	  printf("Copying from %s to %s\n",groupName1,groupName2);
	  sdhdf_copyEntireGroupDifferentLabels(groupName1,inFile,groupName2,outFile);

	  sprintf(groupName1,"%s/%s/metadata",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[0].label); // FIX FOR MORE THAN ONE BAND
	  sprintf(groupName2,"%s/%s/metadata",inFile->beamHeader[b].label,newBandName); // FIX FOR MORE THAN ONE BAND	  
	  printf("Copying from %s to %s\n",groupName1,groupName2);
	  sdhdf_copyEntireGroupDifferentLabels(groupName1,inFile,groupName2,outFile);

	  if (strcmp(inFile->primary[0].cal_mode,"ON")==0)
	    {
	      sprintf(groupName1,"%s/%s/calibrator_data",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[0].label); // FIX FOR MORE THAN ONE BAND
	      sprintf(groupName2,"%s/%s/calibrator_data",inFile->beamHeader[b].label,newBandName); // FIX FOR MORE THAN ONE BAND	  
	      printf("Copying from %s to %s\n",groupName1,groupName2);
	      sdhdf_copyEntireGroupDifferentLabels(groupName1,inFile,groupName2,outFile);
	    }

	  printf("Complete copy\n");
	  // Need to use: sdhdf_copyEntireGroupDifferentLabels
	  //	  sdhdf_copyRemainder(inFile,outFile,0);
	}
      sdhdf_writeBeamHeader(outFile,beamHeader,nBeam,inFile->primary[0].hdr_defn_version);
      for (i=0;i<nOutBands;i++)
	printf("band label = %s\n",outBandParams[i].label);
      sdhdf_writeBandHeader(outFile,outBandParams,inFile0->beamHeader[b].label,nOutBands,1,inFile->primary[0].hdr_defn_version);

      if (strcmp(inFile->primary[0].cal_mode,"ON")==0)
	sdhdf_writeBandHeader(outFile,outCalBandParams,inFile0->beamHeader[b].label,nOutBands,2,inFile->primary[0].hdr_defn_version);
      free(outBandParams);
      free(outCalBandParams);
    }
  sdhdf_addHistory(inFile0->history,inFile0->nHistory,"sdhdf_joinFiles","INSPECTA software to combine files",args);
  printf("args = >%s<\n",args);
  inFile0->nHistory++;
  printf("Writing history\n");
  sdhdf_writeHistory(outFile,inFile0->history,inFile0->nHistory);
  printf("Writing primary Header\n");
  sdhdf_writePrimaryHeader(outFile,primaryHeader);
  printf("Free beamHeader\n");
  //  sdhdf_copyRemainder(inFile0,outFile,0);
  free(beamHeader);
  printf("close file 1\n");
  //  free(primaryHeader);
  sdhdf_closeFile(inFile0);
  printf("close file 2\n");
  sdhdf_closeFile(outFile);
  printf("Ending the code\n");
}


