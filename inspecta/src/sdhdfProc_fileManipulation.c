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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sdhdfProc.h"
#include <ctype.h>

//
// Initialise the sdhdf_fileStruct structure
// This should be run as soon as the memory for the file structure has been allocated
//
void sdhdf_initialiseFile(sdhdf_fileStruct *inFile)
{
  inFile->fileOpen = 0;
  strcpy(inFile->fname,"unset");
  inFile->nPrimary = -1;
  inFile->primaryAllocatedMemory = 0;
  inFile->nSoftware = -1;
  inFile->softwareAllocatedMemory = 0;
  inFile->nHistory = -1;
  inFile->historyAllocatedMemory = 0;
  inFile->nBeam = -1;
  inFile->beamAllocatedMemory = 0;  
}

// Open an SDHDF file
// openType =
//            1 Open for read only
//            2 Open for read/write
//            3 Open for read/write as truncate
// Returns -1 if unsuccessfully opened
//
int sdhdf_openFile(char *fname,sdhdf_fileStruct *inFile,int openType)
{
  if (openType==1)
    inFile->fileID  = H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT);
  else if (openType==2)
    inFile->fileID  = H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
  else if (openType==3)
    inFile->fileID  = H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if (inFile->fileID != -1)
    {
      strcpy(inFile->fname,fname);
      inFile->fileOpen = 1;
    }
  return inFile->fileID;
}
// Routine to close a file that is already open
// it also frees memory used by the sdhdf_fileStruct
//
void sdhdf_closeFile(sdhdf_fileStruct *inFile)
{
  herr_t status;
  int i,j,k;
  if (inFile->fileOpen==1)
    {
      if (inFile->beamAllocatedMemory > 0)
	{
	  for (i=0;i<inFile->nBeam;i++)
	    {
	      if (inFile->beam[i].bandAllocatedMemory > 0)
		{
		  for (j=0;j<inFile->beam[i].nBand;j++)
		    {
		      if (inFile->beam[i].bandData[j].astro_data.freqAllocatedMemory > 0)
			{
			  free(inFile->beam[i].bandData[j].astro_data.freq);
			  inFile->beam[i].bandData[j].astro_data.freqAllocatedMemory = 0;
			}
		      if (inFile->beam[i].bandData[j].astro_data.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.flag); inFile->beam[i].bandData[j].astro_data.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].astro_data.dataWeightsAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.dataWeights); inFile->beam[i].bandData[j].astro_data.dataWeightsAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].astro_data.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.pol1); inFile->beam[i].bandData[j].astro_data.pol1AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].astro_data.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.pol2); inFile->beam[i].bandData[j].astro_data.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].astro_data.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.pol3); inFile->beam[i].bandData[j].astro_data.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].astro_data.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.pol4); inFile->beam[i].bandData[j].astro_data.pol4AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].cal_on_data.freqAllocatedMemory > 0)
			{
			  free(inFile->beam[i].bandData[j].cal_on_data.freq); inFile->beam[i].bandData[j].cal_on_data.freqAllocatedMemory = 0;
			}
		      if (inFile->beam[i].bandData[j].cal_on_data.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.flag); inFile->beam[i].bandData[j].cal_on_data.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.dataWeightsAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.dataWeights); inFile->beam[i].bandData[j].cal_on_data.dataWeightsAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.pol1); inFile->beam[i].bandData[j].cal_on_data.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.pol2); inFile->beam[i].bandData[j].cal_on_data.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.pol3); inFile->beam[i].bandData[j].cal_on_data.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.pol4); inFile->beam[i].bandData[j].cal_on_data.pol4AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.freqAllocatedMemory > 0)
			{
			  free(inFile->beam[i].bandData[j].cal_off_data.freq); inFile->beam[i].bandData[j].cal_off_data.freqAllocatedMemory = 0;
			}
		      if (inFile->beam[i].bandData[j].cal_off_data.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.flag); inFile->beam[i].bandData[j].cal_off_data.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.dataWeightsAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.dataWeights); inFile->beam[i].bandData[j].cal_off_data.dataWeightsAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.pol1); inFile->beam[i].bandData[j].cal_off_data.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.pol2); inFile->beam[i].bandData[j].cal_off_data.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.pol3); inFile->beam[i].bandData[j].cal_off_data.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.pol4); inFile->beam[i].bandData[j].cal_off_data.pol4AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.freqAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.freq); inFile->beam[i].bandData[j].cal_proc_tsys.freqAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.flag); inFile->beam[i].bandData[j].cal_proc_tsys.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.pol1); inFile->beam[i].bandData[j].cal_proc_tsys.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.pol2); inFile->beam[i].bandData[j].cal_proc_tsys.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.pol3); inFile->beam[i].bandData[j].cal_proc_tsys.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.pol4); inFile->beam[i].bandData[j].cal_proc_tsys.pol4AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.freqAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.freq); inFile->beam[i].bandData[j].cal_proc_diff_gain.freqAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.flag); inFile->beam[i].bandData[j].cal_proc_diff_gain.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.pol1); inFile->beam[i].bandData[j].cal_proc_diff_gain.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.pol2); inFile->beam[i].bandData[j].cal_proc_diff_gain.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.pol3); inFile->beam[i].bandData[j].cal_proc_diff_gain.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.pol4); inFile->beam[i].bandData[j].cal_proc_diff_gain.pol4AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.freqAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.freq); inFile->beam[i].bandData[j].cal_proc_diff_phase.freqAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.flag); inFile->beam[i].bandData[j].cal_proc_diff_phase.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.pol1); inFile->beam[i].bandData[j].cal_proc_diff_phase.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.pol2); inFile->beam[i].bandData[j].cal_proc_diff_phase.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.pol3); inFile->beam[i].bandData[j].cal_proc_diff_phase.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.pol4); inFile->beam[i].bandData[j].cal_proc_diff_phase.pol4AllocatedMemory = 0; }
				      
		      if (inFile->beam[i].bandData[j].astro_obsHeaderAllocatedMemory == 1)
			free(inFile->beam[i].bandData[j].astro_obsHeader);
		      if (inFile->beam[i].bandData[j].cal_obsHeaderAllocatedMemory == 1)
			free(inFile->beam[i].bandData[j].cal_obsHeader);

		    }
		  free(inFile->beam[i].bandData);
		  free(inFile->beam[i].bandHeader);
		  if (inFile->beam[i].calBandAllocatedMemory==1)
		    {
		      free(inFile->beam[i].calBandHeader);
		      inFile->beam[i].calBandAllocatedMemory = 0;
		    }
		  inFile->beam[i].bandAllocatedMemory = 0;
		  
		}
	    }
	  free(inFile->beam);
	  inFile->beamAllocatedMemory = 0;
	}

      if (inFile->historyAllocatedMemory > 0)
	{free(inFile->history); inFile->historyAllocatedMemory = 0;}

      if (inFile->softwareAllocatedMemory > 0)
	{free(inFile->software); inFile->softwareAllocatedMemory = 0;}
      if (inFile->primaryAllocatedMemory > 0)
	{free(inFile->primary); inFile->primaryAllocatedMemory = 0;}
      status = H5Fclose(inFile->fileID);
      inFile->fileOpen = 0;
    }
}


// Type = 1: astro_data
// Type = 2: cal_on
// Type = 3: cal_off
void sdhdf_loadBandData2Array(sdhdf_fileStruct *inFile,int beam,int band,int type,float *arr)
{
  int i,j,n;
  int nchan;
  int ndump;
  int npol;
  char dataName[MAX_STRLEN];
  hid_t dataset_id;
  herr_t status;
  float *allData;

  // astro_data
  if (type==1)
    {
      printf("Loading attributes with version: %d\n",inFile->primary[0].hdr_prior4);
      sprintf(dataName,"%s/%s/astronomy_data/data",inFile->beamHeader[beam].label,inFile->beam[beam].bandHeader[band].label);      
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes = sdhdf_getNattributes(inFile,dataName);
      for (j=0;j< inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes;j++)
	sdhdf_readAttributeFromNum(inFile,dataName,j,&(inFile->beam[beam].bandData[band].astro_obsHeaderAttr[j]));
      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,arr);  
      status = H5Dclose(dataset_id);
    }
  
}

// Type = 1: astro_data
// Type = 2: cal_on
// Type = 3: cal_off
void sdhdf_loadQuantisedBandData2Array(sdhdf_fileStruct *inFile,int beam,int band,int type,unsigned char *arr)
{
  int i,j,n;
  int nchan;
  int ndump;
  int npol;
  char dataName[MAX_STRLEN];
  hid_t dataset_id;
  herr_t status;
  float *allData;

  // astro_data
  if (type==1)
    {
      sprintf(dataName,"%s/%s/astronomy_data/data",inFile->beamHeader[beam].label,inFile->beam[beam].bandHeader[band].label);      
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes = sdhdf_getNattributes(inFile,dataName);
      for (j=0;j< inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes;j++)
	sdhdf_readAttributeFromNum(inFile,dataName,j,&(inFile->beam[beam].bandData[band].astro_obsHeaderAttr[j]));
      status = H5Dread(dataset_id,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,arr);  
      status = H5Dclose(dataset_id);
    }
  
}

// Note loads in all dumps (if present)
void sdhdf_loadFrequency2Array(sdhdf_fileStruct *inFile,int beam,int band,float *arr,int *nFreqDump)
{
  int i,j,n,ndims;
  char dataName[MAX_STRLEN];
  hid_t dataset_id,space;
  herr_t status;
  float *allData;
  hsize_t dims[2];
  
  sprintf(dataName,"%s/%s/astronomy_data/frequency",inFile->beamHeader[beam].label,inFile->beam[beam].bandHeader[band].label);      
  dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
  inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq = sdhdf_getNattributes(inFile,dataName);



  printf("Loading attributes with version: %d\n",inFile->primary[0].hdr_prior4);
  for (j=0;j<inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq;j++)
    sdhdf_readAttributeFromNum(inFile,dataName,j,&(inFile->beam[beam].bandData[band].astro_obsHeaderAttr_freq[j]));

  space      = H5Dget_space(dataset_id);
  ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
  status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,arr);  
  *nFreqDump = dims[0];
  status = H5Dclose(dataset_id);

}


int sdhdf_loadWeights2Array(sdhdf_fileStruct *inFile,int beam,int band,float *arr)
{
  int i,j,n,ndims;
  char dataName[MAX_STRLEN];
  hid_t dataset_id,space;
  herr_t status;
  float *allData;
  hsize_t dims[2];
  
  sprintf(dataName,"%s/%s/astronomy_data/weights",inFile->beamHeader[beam].label,inFile->beam[beam].bandHeader[band].label);      
  if (sdhdf_checkGroupExists(inFile,dataName)==0)
    {      
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq = sdhdf_getNattributes(inFile,dataName);
      
      space      = H5Dget_space(dataset_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,arr);  
      printf("Reading status = %d\n",status);
      status = H5Dclose(dataset_id);
      return 1;
    }
  return 0;
}

// FIX ME: SHOULD RETURN IF DON'T HAVE FLAGS
int sdhdf_loadFlags2Array(sdhdf_fileStruct *inFile,int beam,int band,unsigned char *arr)
{
  int i,j,n,ndims;
  char dataName[MAX_STRLEN];
  hid_t dataset_id,space;
  herr_t status;
  float *allData;
  hsize_t dims[2];

  if (sdhdf_checkGroupExists(inFile,dataName)==0)
    {      
      sprintf(dataName,"%s/%s/astronomy_data/flags",inFile->beamHeader[beam].label,inFile->beam[beam].bandHeader[band].label);      
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq = sdhdf_getNattributes(inFile,dataName);
      
      space      = H5Dget_space(dataset_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
      status = H5Dread(dataset_id,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,arr);  
      status = H5Dclose(dataset_id);
      return 1;
    }
  else
    return 0;
}


// Type = 1: astro_data
// Type = 2: cal_on
// Type = 3: cal_off
void sdhdf_loadBandData(sdhdf_fileStruct *inFile,int beam,int band,int type)
{
  int i,j,n;
  int nchan;
  int ndump;
  int npol;
  char dataName[MAX_STRLEN];
  hid_t dataset_id;
  herr_t status;
  float *allData;
  char beamLabel[MAX_STRLEN];
  int calName=1;
  
  //  printf("LOADING DATA type = %d\n",type);
  // astro_data
  if (type==1)
    {
      strcpy(beamLabel,inFile->beamHeader[beam].label);
      nchan = inFile->beam[beam].bandHeader[band].nchan;
      ndump = inFile->beam[beam].bandHeader[band].ndump;
      npol  = inFile->beam[beam].bandHeader[band].npol;
      sdhdf_allocateBandData(&(inFile->beam[beam].bandData[band].astro_data),nchan,ndump,npol);

      {
	hid_t space;
	hsize_t dims[2];
	int ndims;
	      
	//	printf("Loading frequency\n");
	sprintf(dataName,"%s/%s/astronomy_data/frequency",beamLabel,inFile->beam[beam].bandHeader[band].label);      
	dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);

	space      = H5Dget_space(dataset_id);
	ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
	status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].astro_data.freq);  
	if (ndims==1)
	  inFile->beam[beam].bandData[band].astro_data.nFreqDumps = 1;
	else
	  inFile->beam[beam].bandData[band].astro_data.nFreqDumps = dims[0];

	// Increase array if more actual dumps than frequency dumps
	if (inFile->beam[beam].bandData[band].astro_data.nFreqDumps == 1 && inFile->beam[beam].bandHeader[band].ndump > 1)
	  {
	    for (i=1;i<inFile->beam[beam].bandHeader[band].ndump;i++)
	      memcpy(inFile->beam[beam].bandData[band].astro_data.freq+i*nchan, inFile->beam[beam].bandData[band].astro_data.freq,sizeof(float)*nchan);
	    inFile->beam[beam].bandData[band].astro_data.nFreqDumps = inFile->beam[beam].bandHeader[band].ndump;
	  }
	//	printf("Done loading frequency, ndims = %d\n",ndims);
      }
      // Read attributes
            printf("Reading attributes\n");
      inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq = sdhdf_getNattributes(inFile,dataName);
            printf("Number of attributes = %d FREQUENCY SECTION\n",inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq);
      for (j=0;j<inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq;j++)
      	{
	  //	  	  printf("key GEORGE loading: j = %d\n",j);
	  sdhdf_readAttributeFromNum(inFile,dataName,j,&(inFile->beam[beam].bandData[band].astro_obsHeaderAttr_freq[j]));
	  //	  	  printf("Loaded %d >%s<\n",j,inFile->beam[beam].bandData[band].astro_obsHeaderAttr_freq[j].value);
      	}
      //                  printf("Done reading freq attributes\n");
      sprintf(dataName,"%s/%s/astronomy_data/data",beamLabel,inFile->beam[beam].bandHeader[band].label);      

      // Read attributes
      inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes = sdhdf_getNattributes(inFile,dataName);
      //            printf("Reading data attributes\n");
	    //            printf("We have %d\n",inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes);

      //    printf("Reading data attributes\n");
      for (j=0;j< inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes;j++)
      	sdhdf_readAttributeFromNum(inFile,dataName,j,&(inFile->beam[beam].bandData[band].astro_obsHeaderAttr[j]));
	  

      //      printf("Done reading data attributes\n");
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      //      printf("Allocating memory with %d %d %d\n",nchan,npol,ndump);
      allData = (float *)malloc(sizeof(float)*nchan*npol*ndump);      
      if (allData==NULL)
	{
          printf("ERROR: Failure in sdhdfProc_fileManipulation to allocate memory\n");
          printf("ERROR: Cannot recover. Sorry\n");
          exit(1);
	}
      //      printf("Now reading the data %s\n",dataName);
      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,allData);  
      //      printf("Done reading the data\n");
      sdhdf_extractPols(&(inFile->beam[beam].bandData[band].astro_data),allData,nchan,ndump,npol);

      free(allData);
      status = H5Dclose(dataset_id);
      // Flags
      //      printf("Loading flags\n");
      sprintf(dataName,"%s/%s/astronomy_data/flags",beamLabel,inFile->beam[beam].bandHeader[band].label);      
      if (sdhdf_checkGroupExists(inFile,dataName)==0)
	{
	  dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
	  status = H5Dread(dataset_id,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].astro_data.flag);        
	  status = H5Dclose(dataset_id);
	}
      else
	{
	  printf("No flags in the SDHDF file. Setting to 0\n");
	  for (i=0;i<nchan*ndump;i++)
	    inFile->beam[beam].bandData[band].astro_data.flag[i] = 0;
	}

      // Data weights
      //            printf("Loading weights\n");
      sprintf(dataName,"%s/%s/astronomy_data/weights",beamLabel,inFile->beam[beam].bandHeader[band].label);      
      if (sdhdf_checkGroupExists(inFile,dataName)==0)
	{
	  //	  printf("Loading data weights\n");
	  dataset_id = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
	  status     = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].astro_data.dataWeights);
	  if (inFile->beam[beam].bandData[band].astro_data.dataWeights[0] == -1)
	    {
	      int kk,ll;
	      printf("The first weight is -1. This seems incorrect and so setting all weights to +1\n");
	      for (ll=0;ll<ndump;ll++)
		{
		  for (kk=0;kk<inFile->beam[beam].bandHeader[band].nchan;kk++)
		    inFile->beam[beam].bandData[band].astro_data.dataWeights[ll*nchan+kk]=1;
		}
	    }
	  status     = H5Dclose(dataset_id);
	}
      else
	{
	  double chbw,tint;
	  printf("No weightings in the SDHDF file. Setting to the integration time x channel bandwidth\n");

	  tint = inFile->beam[beam].bandHeader[band].dtime;
	  chbw = 1e6*fabs(inFile->beam[beam].bandHeader[band].f0-inFile->beam[beam].bandHeader[band].f1)/(double)inFile->beam[beam].bandHeader[band].nchan;

	  for (i=0;i<nchan*ndump;i++)
	    inFile->beam[beam].bandData[band].astro_data.dataWeights[i] = tint*chbw; 
	}
    }
  // calibrator data
  else if (type==2 || type==3)
    {
      strcpy(beamLabel,inFile->beamHeader[beam].label);
      nchan = inFile->beam[beam].calBandHeader[band].nchan;
      ndump = inFile->beam[beam].calBandHeader[band].ndump;
      npol  = inFile->beam[beam].calBandHeader[band].npol;

      if (type==2)
	sdhdf_allocateBandData(&(inFile->beam[beam].bandData[band].cal_on_data),nchan,ndump,npol);
      else if (type==3)
	sdhdf_allocateBandData(&(inFile->beam[beam].bandData[band].cal_off_data),nchan,ndump,npol);
      
      sprintf(dataName,"%s/%s/calibrator_data/cal_frequency",beamLabel,inFile->beam[beam].bandHeader[band].label);      
      if (sdhdf_checkGroupExists(inFile,dataName)==0)
	{printf("Using cal_frequency to determine the noise source frequency information\n"); calName=1;}
      else
	{sprintf(dataName,"%s/%s/calibrator_data/frequency",beamLabel,inFile->beam[beam].bandHeader[band].label); calName=2;}
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      // FIX ME: freq[0]
      {
	hid_t space;
	hsize_t dims[2];
	int ndims;

	space      = H5Dget_space(dataset_id);
	ndims      = H5Sget_simple_extent_dims(space,dims,NULL);


	if (type==2)
	  status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].cal_on_data.freq);
	else if (type==3)
	  status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].cal_off_data.freq);

	printf("Number of dimentions = %d\n",ndims);
	
	if (ndims==1)
	  {
 	    if (type==2) inFile->beam[beam].bandData[band].cal_on_data.nFreqDumps = 1;
	    else if (type==3)inFile->beam[beam].bandData[band].cal_off_data.nFreqDumps = 1;
	  }	 
	else
	  {
 	    if (type==2) inFile->beam[beam].bandData[band].cal_on_data.nFreqDumps = dims[0];
	    else if (type==3)inFile->beam[beam].bandData[band].cal_off_data.nFreqDumps = dims[0];
	  }
	// FIX ME: SHOULD INCREASE NUMBER OF FREQDUMPS TO NDUMP
      }
      if (calName==1)
	{
	  if (type==2)
	    sprintf(dataName,"%s/%s/calibrator_data/cal_data_on",beamLabel,inFile->beam[beam].bandHeader[band].label);      
	  else if (type==3)
	    sprintf(dataName,"%s/%s/calibrator_data/cal_data_off",beamLabel,inFile->beam[beam].bandHeader[band].label);
	}
      else
	{
	  if (type==2)
	    sprintf(dataName,"%s/%s/calibrator_data/cal_on",beamLabel,inFile->beam[beam].bandHeader[band].label);      
	  else if (type==3)
	    sprintf(dataName,"%s/%s/calibrator_data/cal_off",beamLabel,inFile->beam[beam].bandHeader[band].label);

	}
	  
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      allData = (float *)malloc(sizeof(float)*nchan*npol*ndump);
      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,allData);  
      if (type==2)
	{
	  //	  printf("Extrating with %d %d %d\n",nchan,ndump,npol);
	  sdhdf_extractPols(&(inFile->beam[beam].bandData[band].cal_on_data),allData,nchan,ndump,npol);
	}
	  else if (type==3)
	sdhdf_extractPols(&(inFile->beam[beam].bandData[band].cal_off_data),allData,nchan,ndump,npol);
      free(allData);

      status = H5Dclose(dataset_id);
    }
}

void sdhdf_loadDataFreqAttributes(sdhdf_fileStruct *inFile,int beam,int band,sdhdf_attributes_struct *dataAttributes,int *nDataAttributes,
				  sdhdf_attributes_struct *freqAttributes,int *nFreqAttributes)
{
  int j;
  char dataName[MAX_STRLEN];
  char beamLabel[MAX_STRLEN];
  
  strcpy(beamLabel,inFile->beamHeader[beam].label);
  sprintf(dataName,"%s/%s/astronomy_data/frequency",beamLabel,inFile->beam[beam].bandHeader[band].label);      
  *nFreqAttributes = sdhdf_getNattributes(inFile,dataName);
  //      printf("Number of attributes = %d FREQUENCY SECTION\n",inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq);
  for (j=0;j<*nFreqAttributes;j++)
      sdhdf_readAttributeFromNum(inFile,dataName,j,&freqAttributes[j]);
  
  sprintf(dataName,"%s/%s/astronomy_data/data",beamLabel,inFile->beam[beam].bandHeader[band].label);      

  // Read attributes
  *nDataAttributes = sdhdf_getNattributes(inFile,dataName);
  for (j=0;j< *nDataAttributes; j++)
    sdhdf_readAttributeFromNum(inFile,dataName,j,&dataAttributes[j]);
  

}

void sdhdf_extractPols(sdhdf_spectralDumpsStruct *spec,float *in,int nchan,int ndump,int npol)
{
  unsigned long int i,j;
  unsigned long int ii;
  
  //  printf("Extracting polarisations\n");

  // SHOULD USE MEMCPY -- FIX ME
  // Note: nchan*ndump*npol can be greater than INT_MAX
  //
  for (i=0;i<ndump;i++)
    {
      //      printf("dump = %d\n",i);
      for (j=0;j<nchan;j++)
	{
	  //	  printf("Looking for %d\n",j+nchan*i*npol);
	  spec->pol1[(unsigned long int)j+(unsigned long int)nchan*i] = in[(unsigned long int)j+(unsigned long int)nchan*i*(unsigned long int)npol];
	  if (npol > 1)
	    {
	      
	      spec->pol2[j+(unsigned long int)nchan*i] = in[j+(unsigned long int)nchan*i*(unsigned long int)npol+(unsigned long int)nchan];
	      if (npol > 2)
		{
		  spec->pol3[j+(unsigned long int)nchan*i] = in[j+(unsigned long int)nchan*i*(unsigned long int)npol+2*(unsigned long int)nchan];
		  spec->pol4[j+(unsigned long int)nchan*i] = in[j+(unsigned long int)nchan*i*(unsigned long int)npol+3*(unsigned long int)nchan];
		}
	    }
	}
    }

}

void sdhdf_allocateBandData(sdhdf_spectralDumpsStruct *spec,int nchan,int ndump,int npol)
{
  int k;
  if (spec->freqAllocatedMemory == 0)
    {
      spec->freq = (float *)malloc(sizeof(float *)*ndump*nchan);
      spec->freqAllocatedMemory = 1;
      spec->nFreqDumps = ndump;
    }
  if (spec->flagAllocatedMemory == 0)
    {
      spec->flag = (unsigned char *)calloc(sizeof(unsigned char),nchan*ndump); // Initialise to 0
      spec->flagAllocatedMemory = 1;
    }
  if (spec->dataWeightsAllocatedMemory == 0)
    {
      int i;

      spec->dataWeights = (float *)malloc(sizeof(float)*nchan*ndump); // Initialise to -1
      for (i=0;i<nchan*ndump;i++)
	spec->dataWeights[i]=-1;
      spec->dataWeightsAllocatedMemory = 1;
    }
  if (spec->pol1AllocatedMemory == 0)
    {
      spec->pol1 = (float *)malloc(sizeof(float)*nchan*npol*ndump); 
      spec->pol1AllocatedMemory = 1;
    }
  if (npol > 1 && spec->pol2AllocatedMemory == 0)
    {
      spec->pol2 = (float *)malloc(sizeof(float)*nchan*npol*ndump); 
      spec->pol2AllocatedMemory = 1;
    }
  if (npol > 3 && spec->pol3AllocatedMemory == 0 && spec->pol4AllocatedMemory == 0)
    {
      spec->pol3 = (float *)malloc(sizeof(float)*nchan*ndump);
      spec->pol4 = (float *)malloc(sizeof(float)*nchan*ndump); 
      spec->pol3AllocatedMemory = 1;
      spec->pol4AllocatedMemory = 1;
    }
}
  

// Type = 1: astro_data
// FIX ME __ OTHER TYPES
void sdhdf_releaseBandData(sdhdf_fileStruct *inFile,int beam,int band,int type)
{
  int k;
  
  if (type==1)
    {
      if (inFile->beam[beam].bandData[band].astro_data.freqAllocatedMemory == 1)
	{
	  free(inFile->beam[beam].bandData[band].astro_data.freq);
	  inFile->beam[beam].bandData[band].astro_data.freqAllocatedMemory = 0;
	}
      if (inFile->beam[beam].bandData[band].astro_data.flagAllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.flag); inFile->beam[beam].bandData[band].astro_data.flagAllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.dataWeightsAllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.dataWeights); inFile->beam[beam].bandData[band].astro_data.dataWeightsAllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.pol1AllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.pol1); inFile->beam[beam].bandData[band].astro_data.pol1AllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.pol2AllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.pol2); inFile->beam[beam].bandData[band].astro_data.pol2AllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.pol3AllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.pol3); inFile->beam[beam].bandData[band].astro_data.pol3AllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.pol4AllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.pol4); inFile->beam[beam].bandData[band].astro_data.pol4AllocatedMemory = 0;}	      
    }
}

void sdhdf_add1arg(char *args,char *add)
{
  if (strlen(args) < MAX_STRLEN - strlen(add) - 1)
    {
      strcat(args,add);
      strcat(args," ");
    }
}


void sdhdf_add2arg(char *args,char *add1,char *add2)
{
  if (strlen(args) < MAX_STRLEN - strlen(add1) - 1)
    {
      strcat(args,add1);
      strcat(args," ");
    }
  if (strlen(args) < MAX_STRLEN - strlen(add2) - 1)
    {
      strcat(args,add2);
      strcat(args," ");
    }
}


// ** SHOULD COPY FLAG TABLE AS WELL **

// type = 0: update DATA and frequency, but copy the cal datasets
// type = 1: update data and frequency only
// type = 2: update cal_on and frequency only
// type = 3: update cal_off only
// type = 4: update flag, data and frequency only


// GEORGE: NEED TO ADD ATTRIBUTES HERE -- PASS THEM THROUGH
void sdhdf_writeSpectrumData(sdhdf_fileStruct *outFile,char *beamLabel,char *blabel, int ibeam,int iband,  float *out,float *freq,int nFreqDump,long nchan,long nbin,long npol,long nsub,int type,sdhdf_attributes_struct *dataAttributes,int nDataAttributes,sdhdf_attributes_struct *freqAttributes,int nFreqAttributes)
{
  int i;
  hid_t dset_id,datatype_id,group_id,ocpl_id;
  herr_t status=0;
  hid_t dataspace_id,stid;
  char groupName[MAX_STRLEN];
  char dsetName[MAX_STRLEN],label[MAX_STRLEN];
  hsize_t dims[5];

  
  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"%s",beamLabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status   = H5Gclose(group_id);
	}

      sprintf(groupName,"%s/%s",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}

      sprintf(groupName,"%s/%s/astronomy_data",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
    }
  else if (type==2 || type==3)
    {
      sprintf(groupName,"%s",beamLabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status   = H5Gclose(group_id);
	}

      sprintf(groupName,"%s/%s",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
      sprintf(groupName,"%s/%s/calibrator_data",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
    }
  dims[0] = nsub;
  dims[1] = npol;
  dims[2] = nchan;
  dims[3] = nbin;
  dataspace_id = H5Screate_simple(4,dims,NULL);

  /* FIX ME -- SHOULD COPY CAL INFO
  if (type==0)    
    {
      ocpl_id = H5Pcreate(H5P_OBJECT_COPY);
      sprintf(label,"%s/cal_data_on",groupName);
      status = H5Ocopy(inFile->fileID, label, outFile->fileID, label, ocpl_id, H5P_DEFAULT);
      sprintf(label,"%s/cal_data_off",groupName);
      status = H5Ocopy(inFile->fileID, label, outFile->fileID, label, ocpl_id, H5P_DEFAULT);
      sprintf(label,"%s/cal_frequency",groupName);
      status = H5Ocopy(inFile->fileID, label, outFile->fileID, label, ocpl_id, H5P_DEFAULT);

      H5Pclose(ocpl_id);
    }
  */

  // FIX ME
  /*
  if (type==4)
    {
      ocpl_id = H5Pcreate(H5P_OBJECT_COPY);
      sprintf(label,"%s/flag",groupName);
      status = H5Ocopy(inFile->fileID, label, outFile->fileID, label, ocpl_id, H5P_DEFAULT);
      H5Pclose(ocpl_id);
    }
  */

  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"%s/%s/astronomy_data",beamLabel,blabel);
      sprintf(dsetName,"%s/data",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
      for (i=0;i<nDataAttributes;i++)
	{
	  //	  printf("GEORGE: in fileManipulation: writing data attribute\n");
	  //	  printf("GEORGE: in fileManipulation: %d %d\n",i,nDataAttributes);
	  //	  printf("GEORGE: in fileManipulation: %s\n",dataAttributes[i].key);
	  sdhdf_writeAttribute(outFile,dsetName,&dataAttributes[i]); //dataAttributes[i].key,dataAttributes[i].value);
	  //	  printf("GEORGE: complete writing the attribute\n");
	}

    }
  else if (type==2)
    {
      sprintf(groupName,"%s/%s/calibrator_data",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
      sprintf(dsetName,"%s/cal_data_on",groupName);
      printf("CREATING %s\n",dsetName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
    }
  else if (type==3)
    {
      sprintf(groupName,"%s/%s/calibrator_data",beamLabel,blabel);
      sprintf(dsetName,"%s/cal_data_off",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
    }
  
  status = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,out);  
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);

  if (type==2 || type==3)
    {
      dims[0] = 1; // FIX ME IF NEEDED **
      dims[1] = nchan;
    }
  else
    {
      dims[0] = nFreqDump;
      dims[1] = nchan;
    }
      dataspace_id = H5Screate_simple(2,dims,NULL);

  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"%s/%s/astronomy_data",beamLabel,blabel);
      sprintf(dsetName,"%s/frequency",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);     
      status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,freq);  
    
      // Now write attributes
      //      printf("Number of attributes to write out = %d\n",outFile->beam[ibeam].bandData[iband].nAstro_obsHeaderAttributes_freq);
      for (i=0;i<nFreqAttributes;i++)
	sdhdf_writeAttribute(outFile,dsetName,&freqAttributes[i]); //freqAttributes[i].key,freqAttributes[i].value);
      
      status  = H5Dclose(dset_id);
    }
  else if (type==2)
    {
      sprintf(groupName,"%s/%s/calibrator_data",beamLabel,blabel);
      sprintf(dsetName,"%s/cal_frequency",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);     
      status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,freq);  
      status  = H5Dclose(dset_id);      
    }
  status = H5Sclose(dataspace_id);
  if (type == 0 || type == 1 || type == 4)
    status = H5Gclose(group_id);


}

// GEORGE: NEED TO ADD ATTRIBUTES HERE -- PASS THEM THROUGH
// Writes out unsigned characters
void sdhdf_writeQuantisedSpectrumData(sdhdf_fileStruct *outFile,char *beamLabel,char *blabel, int ibeam,int iband,  unsigned char *out,float *freq,int nFreqDump,long nchan,long nbin,long npol,long nsub,int type,sdhdf_attributes_struct *dataAttributes,int nDataAttributes,sdhdf_attributes_struct *freqAttributes,int nFreqAttributes)
{
  int i;
  hid_t dset_id,datatype_id,group_id,ocpl_id;
  herr_t status=0;
  hid_t dataspace_id,stid;
  char groupName[MAX_STRLEN];
  char dsetName[MAX_STRLEN],label[MAX_STRLEN];
  hsize_t dims[5];

  
  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"%s",beamLabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status   = H5Gclose(group_id);
	}

      sprintf(groupName,"%s/%s",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}

      sprintf(groupName,"%s/%s/astronomy_data",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
    }
  else if (type==2 || type==3)
    {
      sprintf(groupName,"%s",beamLabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status   = H5Gclose(group_id);
	}

      sprintf(groupName,"%s/%s",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
      sprintf(groupName,"%s/%s/calibrator_data",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
    }
  dims[0] = nsub;
  dims[1] = npol;
  dims[2] = nchan;
  dims[3] = nbin;
  dataspace_id = H5Screate_simple(4,dims,NULL);


  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"%s/%s/astronomy_data",beamLabel,blabel);
      sprintf(dsetName,"%s/data",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_UCHAR,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
      for (i=0;i<nDataAttributes;i++)
	{
	  sdhdf_writeAttribute(outFile,dsetName,&dataAttributes[i]); //dataAttributes[i].key,dataAttributes[i].value);
	}

    }
  else if (type==2)
    {
      sprintf(groupName,"%s/%s/calibrator_data",beamLabel,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
      sprintf(dsetName,"%s/cal_data_on",groupName);
      printf("CREATING %s\n",dsetName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_UCHAR,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
    }
  else if (type==3)
    {
      sprintf(groupName,"%s/%s/calibrator_data",beamLabel,blabel);
      sprintf(dsetName,"%s/cal_data_off",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_UCHAR,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
    }
  
  status = H5Dwrite(dset_id,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,out);  
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);

  if (type==2 || type==3)
    {
      dims[0] = 1; // FIX ME IF NEEDED **
      dims[1] = nchan;
    }
  else
    {
      dims[0] = nFreqDump;
      dims[1] = nchan;
    }
      dataspace_id = H5Screate_simple(2,dims,NULL);

  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"%s/%s/astronomy_data",beamLabel,blabel);
      sprintf(dsetName,"%s/frequency",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);     
      status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,freq);  
    
      // Now write attributes
      //      printf("Number of attributes to write out = %d\n",outFile->beam[ibeam].bandData[iband].nAstro_obsHeaderAttributes_freq);
      for (i=0;i<nFreqAttributes;i++)
	sdhdf_writeAttribute(outFile,dsetName,&freqAttributes[i]); //freqAttributes[i].key,freqAttributes[i].value);
      
      status  = H5Dclose(dset_id);
    }
  else if (type==2)
    {
      sprintf(groupName,"%s/%s/calibrator_data",beamLabel,blabel);
      sprintf(dsetName,"%s/cal_frequency",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);     
      status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,freq);  
      status  = H5Dclose(dset_id);      
    }
  status = H5Sclose(dataspace_id);
  if (type == 0 || type == 1 || type == 4)
    status = H5Gclose(group_id);


}


void sdhdf_replaceSpectrumData(sdhdf_fileStruct *outFile,char *blabel, int ibeam,int iband,  float *out,int nsub,int npol,int nchan)
{
  int i;
  hid_t dset_id,datatype_id,group_id,ocpl_id;
  herr_t status=0;
  hid_t dataspace_id,stid;
  char groupName[MAX_STRLEN];
  char dsetName[MAX_STRLEN],label[MAX_STRLEN];
  hsize_t dims[5];

  sprintf(groupName,"beam_%d/%s/astronomy_data",ibeam,blabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  dims[0] = nsub;
  dims[1] = 1;
  dims[2] = npol;
  dims[3] = nchan;
  dims[4] = 1;
  // FIX ME -- NOW 4 DIMENSIONS
  dataspace_id = H5Screate_simple(5,dims,NULL);
  
  sprintf(dsetName,"%s/data",groupName);
  dset_id = H5Dopen2(outFile->fileID,dsetName,H5P_DEFAULT);
  status = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,out);  
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);

  //  status = H5Gclose(group_id);


}

// type = 0: copy all
// type = 1: do not copy beam data
void sdhdf_copyRemainder(sdhdf_fileStruct *inFile,sdhdf_fileStruct *outFile,int type)
{
  char groupName[1024];
  hid_t group_id;
  herr_t status;
  int b,i;
  
  // File information
  sprintf(groupName,"/metadata");
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    sdhdf_copyEntireGroup(groupName,inFile,outFile);
  else
    {
      if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)       
	sprintf(groupName,"/metadata/beam_parameters");
      else
	sprintf(groupName,"/metadata/beam_params");
      if (sdhdf_checkGroupExists(outFile,groupName) == 1) sdhdf_copyEntireGroup(groupName,inFile,outFile);
      sprintf(groupName,"/metadata/history");
      if (sdhdf_checkGroupExists(outFile,groupName) == 1) sdhdf_copyEntireGroup(groupName,inFile,outFile);
      sprintf(groupName,"/metadata/primary_header");
      if (sdhdf_checkGroupExists(outFile,groupName) == 1) sdhdf_copyEntireGroup(groupName,inFile,outFile);
      sprintf(groupName,"/metadata/software_versions");
      if (sdhdf_checkGroupExists(outFile,groupName) == 1) sdhdf_copyEntireGroup(groupName,inFile,outFile);
    }
  
  // Config
  sprintf(groupName,"/config");
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    sdhdf_copyEntireGroup(groupName,inFile,outFile);


  // Beams
  if (type!=1)
    {
      for (b=0;b<inFile->nBeam;b++)
	{
	  sprintf(groupName,"/%s",inFile->beamHeader[b].label);
	  //	  printf("HERE: GROUP: %s %d\n",groupName,sdhdf_checkGroupExists(outFile,groupName));
	  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	    {
	      sdhdf_copyEntireGroup(groupName,inFile,outFile);
	    }
	  else
	    {
	      sprintf(groupName,"/%s/metadata",inFile->beamHeader[b].label);
	      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
		sdhdf_copyEntireGroup(groupName,inFile,outFile);
	      else
		{
		  if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)       
		    sprintf(groupName,"/%s/metadata/band_parameters",inFile->beamHeader[b].label);
		  else
		    sprintf(groupName,"/%s/metadata/band_params",inFile->beamHeader[b].label);
		  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
		    sdhdf_copyEntireGroup(groupName,inFile,outFile);
		  if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)       
		    sprintf(groupName,"/%s/metadata/cal_band_parameters",inFile->beamHeader[b].label);
		  else
		    sprintf(groupName,"/%s/metadata/cal_band_params",inFile->beamHeader[b].label);

		  if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
		    sdhdf_copyEntireGroup(groupName,inFile,outFile);
		}
	      
	      for (i=0;i<inFile->beam[b].nBand;i++)
		{
		  sprintf(groupName,"/%s/%s",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
		  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
		    sdhdf_copyEntireGroup(groupName,inFile,outFile);
		  else
		    {
		      // Metadata
		      sprintf(groupName,"/%s/%s/metadata",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
		      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
			sdhdf_copyEntireGroup(groupName,inFile,outFile);
		      else
			{
			  if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)       
			    sprintf(groupName,"/%s/%s/metadata/cal_observation_parameters",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
			  else
			    sprintf(groupName,"/%s/%s/metadata/cal_obs_params",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
			    
			  if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
			    sdhdf_copyEntireGroup(groupName,inFile,outFile);
			  if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)       
			    sprintf(groupName,"/%s/%s/metadata/observation_parameters",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
			  else
			    sprintf(groupName,"/%s/%s/metadata/obs_params",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
			  if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
			    sdhdf_copyEntireGroup(groupName,inFile,outFile);
			  
			}
		      // calibrator_data
		      sprintf(groupName,"/%s/%s/calibrator_data",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
		      if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
			sdhdf_copyEntireGroup(groupName,inFile,outFile);
		      
		      // astronomy_data
		      sprintf(groupName,"/%s/%s/astronomy_data",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
		      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
			sdhdf_copyEntireGroup(groupName,inFile,outFile);
		      
		      // Flagging information
		      sprintf(groupName,"%s/%s/astronomy_data/flags",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
		      if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
			sdhdf_copyEntireGroup(groupName,inFile,outFile);

		      // Data weights information
		      sprintf(groupName,"%s/%s/astronomy_data/weights",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
		      if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
			sdhdf_copyEntireGroup(groupName,inFile,outFile);
		    }
		}
	    }
	}
    }

}

void sdhdf_loadCalProc(sdhdf_fileStruct *inFile,int ibeam,int iband,char *cal_label,float *vals)
{
  int i;
  hid_t dataset_id;
  herr_t status;
  char dataName[1024];

  sprintf(dataName,"/%s/%s/calibrator_proc/%s",inFile->beamHeader[ibeam].label,inFile->beam[ibeam].bandHeader[iband].label,cal_label);
  dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
  status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,vals);  
  status = H5Dclose(dataset_id);
}

void sdhdf_writeCalProc(sdhdf_fileStruct *outFile,int ibeam,int iband,char *band_label,char *cal_label,float *vals,int nchan,int npol,int ndumps)
{
  int i;
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  char groupName[MAX_STRLEN];
  char dsetName[MAX_STRLEN],label[MAX_STRLEN];
  hsize_t dims[5];

  // Create group if needed
  sprintf(groupName,"/beam_%d/%s/calibrator_proc",ibeam,band_label);

  // Check if group exists
  status = H5Eset_auto1(NULL,NULL);
  status = H5Gget_objinfo(outFile->fileID,groupName,0,NULL);
  if (status == 0)  
    {
      //printf("cal_proc already exists\n");
    }
  else
    {
      status=0;
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    }
  dims[0] = ndumps;
  dims[1] = 1;
  dims[2] = npol;
  dims[3] = nchan;
  dims[4] = 1;
  dataspace_id = H5Screate_simple(5,dims,NULL);
  
  sprintf(dsetName,"beam_%d/%s/calibrator_proc/%s",ibeam,band_label,cal_label);
  dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
  status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,vals);  
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);
  status = H5Gclose(group_id);      

}


void sdhdf_copyEntireGroup(char *bandLabel,sdhdf_fileStruct *inFile,sdhdf_fileStruct *outFile)
{
  hid_t loc_id,ocpl_id,lcpl_id;
  herr_t status;

  loc_id = inFile->fileID;

  ocpl_id = H5Pcreate(H5P_OBJECT_COPY);
  lcpl_id = H5Pcreate(H5P_LINK_CREATE);


  status = H5Ocopy(loc_id, bandLabel, outFile->fileID, bandLabel, ocpl_id, lcpl_id);

  H5Pclose(ocpl_id);
  H5Pclose(lcpl_id);
}

void sdhdf_copyEntireGroupDifferentLabels(char *bandLabelIn,sdhdf_fileStruct *inFile,char *bandLabelOut,sdhdf_fileStruct *outFile)
{
  hid_t loc_id,ocpl_id,lcpl_id;
  herr_t status;

  loc_id = inFile->fileID;

  ocpl_id = H5Pcreate(H5P_OBJECT_COPY);
  lcpl_id = H5Pcreate(H5P_LINK_CREATE);


  status = H5Ocopy(loc_id, bandLabelIn, outFile->fileID, bandLabelOut, ocpl_id, lcpl_id);

  H5Pclose(ocpl_id);
  H5Pclose(lcpl_id);
}

void sdhdf_writeFlags(sdhdf_fileStruct *outFile,int ibeam,int iband,unsigned char *flag,int nchan,int ndump,char *beamLabel,char *bandLabel)
{
  char dSetName[MAX_STRLEN];
  char groupName[MAX_STRLEN];
  hid_t dataspace_id,stid,dset_id,datatype_id,group_id;
  hsize_t dims[2];
  herr_t status;
  
  dims[0] = ndump; 
  dims[1] = nchan;
  dataspace_id = H5Screate_simple(2,dims,NULL);

  sprintf(dSetName,"%s/%s/astronomy_data/flags",beamLabel,bandLabel);
  if (sdhdf_checkGroupExists(outFile,dSetName) == 1)
    dset_id = H5Dcreate2(outFile->fileID,dSetName,H5T_NATIVE_UCHAR,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
  else
    dset_id = H5Dopen2(outFile->fileID,dSetName,H5P_DEFAULT);
  
  status  = H5Dwrite(dset_id,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,flag);  
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);
}

void sdhdf_writeDataWeights(sdhdf_fileStruct *outFile,int ibeam,int iband,float *flag,int nchan,int ndump,char *beamLabel,char *bandLabel)
{
  char dSetName[MAX_STRLEN];
  char groupName[MAX_STRLEN];
  hid_t dataspace_id,stid,dset_id,datatype_id,group_id;
  hsize_t dims[2];
  herr_t status;


  printf("In here and writing weights for ndump and nchan = %d and %d\n",ndump,nchan);
  dims[0] = ndump; 
  dims[1] = nchan;
  dataspace_id = H5Screate_simple(2,dims,NULL);

  sprintf(dSetName,"%s/%s/astronomy_data/weights",beamLabel,bandLabel);
  printf("Writing to %s\n",dSetName);
  if (sdhdf_checkGroupExists(outFile,dSetName) == 1)
    dset_id = H5Dcreate2(outFile->fileID,dSetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
  else
    dset_id = H5Dopen2(outFile->fileID,dSetName,H5P_DEFAULT);
  if (dset_id < 0) // Returns negative number if not successful 
    {
      printf("WARNING: CANNOT WRITE WEIGHTS\n");
    }
  else
    {  
      status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,flag);  
    }
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);
}

void sdhdf_loadRestFrequencies(sdhdf_restfrequency_struct *restFreq,int *nRestFreq)
{
  char runtimeDir[1024];
  char fname[1024];
  FILE *fin;
  int i;
  char line[4096];
  char *tok;
  char trim[4096];
  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=========================================================\n");
      printf("Error: we require that the SDHDF_RUNTIME directory is set\n");
      printf("=========================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  sprintf(fname,"%s/spectralLines/spectralLines.dat",runtimeDir);
  
  *nRestFreq=0;
  if (!(fin = fopen(fname,"r")))
    printf("ERROR: unable to open file: %s .... trying to continue \n",fname);
  else
    {
      while (!feof(fin))
	{
	  if (fgets(line,4096,fin)!=NULL)
	    {
	      strcpy(trim,sdhdf_trim(line)); strcpy(line,trim);
	      if (line[0]=='#' || strlen(line)<2) // Comment line
		{}
	      else
		{
		  tok = strtok(line," ");
		  strcpy(trim,sdhdf_trim(tok));
		  sscanf(trim,"%lf",&(restFreq[*nRestFreq].f0));
		  tok = strtok(NULL," ");
		  strcpy(trim,sdhdf_trim(tok));
		  sscanf(trim,"%d",&(restFreq[*nRestFreq].flag));
		  tok = strtok(NULL,"\n");
		  strcpy(trim,sdhdf_trim(tok));
		  strcpy(restFreq[*nRestFreq].label,trim);
		  (*nRestFreq)++;
		}		  
	    }					 
	}
    }
}

char *sdhdf_ltrim(char *s)
{
    while(isspace(*s)) s++;
    return s;
}

char *sdhdf_rtrim(char *s)
{
    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}

char *sdhdf_trim(char *s)
{
    return sdhdf_rtrim(sdhdf_ltrim(s)); 
}
