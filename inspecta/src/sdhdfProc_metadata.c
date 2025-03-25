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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include "inspecta.h"
#include <locale.h>
#include <wchar.h>
#include <time.h>
#include <sys/utsname.h>

void sdhdf_loadConfig(sdhdf_fileStruct *inFile);

void sdhdf_storeArguments(char *args,int maxLen,int argc,char *argv[])
{
  int i;
  int len1,len2;
  for (i=1;i<argc;i++)
    {
      len1 = strlen(args);
      len2 = strlen(argv[i]);
      if (len1 + len2 > maxLen-5)
	{strcat(args,"..."); break;}
      else
	{
	  strcat(args,argv[i]);
	  strcat(args," ");
	}
    }
}


// Sort out the filename when adding on an extension
// If ending of filename = hdf then remove that, add the extension and then add .hdf at the end
// If ending of filename = hdf5 then remove that, add the extension and then add .hdf5 at the end
// Otherwise just append to the end
void sdhdf_formOutputFilename(char *inFile,char *extension,char *oname)
{
  if (strcmp(inFile+strlen(inFile)-4,".hdf")==0)
    {
      char tmp[1024];
      strcpy(tmp,inFile);
      tmp[strlen(tmp)-4]='\0';
      sprintf(oname,"%s.%s.hdf",tmp,extension);
    }
  else if (strcmp(inFile+strlen(inFile)-5,".hdf5")==0)
    {
      char tmp[1024];
      strcpy(tmp,inFile);
      tmp[strlen(tmp)-5]='\0';
      sprintf(oname,"%s.%s.hdf",tmp,extension);
    }
  else if (strcmp(inFile+strlen(inFile)-6,".sdhdf")==0)
    {
      char tmp[1024];
      strcpy(tmp,inFile);
      tmp[strlen(tmp)-6]='\0';
      sprintf(oname,"%s.%s.hdf",tmp,extension);
    }
  else
    sprintf(oname,"%s.%s",inFile,extension);
}


// Fix underscores for GIZA
void sdhdf_fixUnderscore(char *input,char *output)
{
  char *tok;
  char tstr[1024];
  int pos=0;
  strcpy(tstr,input);
  tok = strtok(tstr,"_");
  strcpy(output,"");
  while (tok!=NULL)
    {
      if (pos>0)
	strcat(output,"\\_");
      strcat(output,tok);	      
      tok = strtok(NULL,"_");
      pos=1;
    } 
}

//
// Load all available metadata from the SDHDF file
//
void sdhdf_loadMetaData(sdhdf_fileStruct *inFile)  // Include loading attributes
{
  sdhdf_loadPrimaryHeader(inFile);
  sdhdf_loadConfig(inFile);
  sdhdf_loadBeamHeader(inFile);
  sdhdf_loadBandHeader(inFile,1);
  sdhdf_loadObsHeader(inFile,1);
  if (inFile->beam[0].bandData[0].astro_obsHeader[0].dtime == -1) // Check if dump time not in the observation parameters
    {
      int i,j,b;
      for (b=0;b<inFile->nBeam;b++)
	{
	  for (i=0;i<inFile->beam[b].nBand;i++)
	    {
	      for (j=0;j<inFile->beam[b].bandHeader[i].ndump;j++)
		inFile->beam[b].bandData[i].astro_obsHeader[j].dtime = inFile->beam[b].bandHeader[i].dtime;
	    }
	}
    }

  if (strcmp(inFile->primary[0].cal_mode,"ON")==0)
    {
      sdhdf_loadBandHeader(inFile,2);
      sdhdf_loadObsHeader(inFile,2);
    }
  sdhdf_loadHistory(inFile);
  sdhdf_loadSoftware(inFile);  
}

void sdhdf_loadConfig(sdhdf_fileStruct *inFile)
{
  int i,j;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  sdhdf_backendConfigStruct *configVals;
  char groupName[128];
  
  sprintf(groupName,"config");
  if (sdhdf_checkGroupExists(inFile,groupName) == 1)
    {
      //      printf("Warning: No configuration table in SDHDF file\n");
    }
  else
    {  
      header_id  = H5Dopen2(inFile->fileID,"config/backend_config",H5P_DEFAULT);
      headerT    = H5Dget_type(header_id);
      space      = H5Dget_space(header_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
      
      configVals = (sdhdf_backendConfigStruct *)malloc(sizeof(sdhdf_backendConfigStruct)*dims[0]);
      
      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_backendConfigStruct));
      stid    = H5Tcopy(H5T_C_S1);
      status  = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,"BACKEND_PHASE",HOFFSET(sdhdf_backendConfigStruct,backend_phase),stid);
      H5Tinsert(val_tid,"CAL_FREQ",HOFFSET(sdhdf_backendConfigStruct,cal_freq),stid);
      H5Tinsert(val_tid,"CAL_EPOCH",HOFFSET(sdhdf_backendConfigStruct,cal_epoch),stid);
      H5Tinsert(val_tid,"CAL_DUTY_CYCLE",HOFFSET(sdhdf_backendConfigStruct,cal_duty_cycle),stid);
      H5Tinsert(val_tid,"CAL_PHASE",HOFFSET(sdhdf_backendConfigStruct,cal_phase),stid);
      status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,configVals);
      
      strcpy(inFile->cal_epoch,configVals[0].cal_epoch);
      sscanf(configVals[0].cal_freq,"%lf",&(inFile->cal_freq));
      sscanf(configVals[0].cal_phase,"%lf",&(inFile->cal_phase));
      sscanf(configVals[0].cal_duty_cycle,"%lf",&(inFile->cal_duty_cycle));
      
      //  printf("Loaded %s %s %s\n",configVals[0].backend_phase,configVals[0].cal_freq,configVals[0].cal_epoch);
      free(configVals);
      
      status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(header_id);
    }
}

//
// Load the primary header meta data in the metadata group
//
void sdhdf_loadPrimaryHeader(sdhdf_fileStruct *inFile)
{
  int i,j;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  int ver4plus=0;

  if (sdhdf_checkGroupExists(inFile,"metadata/primary_header") == 1)
    {
      printf("WARNING: We do not have metadata/primary_header in the data file %s\n",inFile->fname);
    }
  else
    {
      header_id  = H5Dopen2(inFile->fileID,"metadata/primary_header",H5P_DEFAULT);
      headerT    = H5Dget_type(header_id);
      space      = H5Dget_space(header_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
      
      if (inFile->primaryAllocatedMemory == 0)
	{
	  inFile->primary = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct)*dims[0]);
	  inFile->nPrimary = dims[0];
	  inFile->primaryAllocatedMemory = 1;
	  if (inFile->nPrimary != 1)
	    {
	      printf("ERROR: UNEXPECTED SIZE IN PRIMARY HEADER\n");
	      exit(1);
	    }
	}
      
      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_primaryHeaderStruct));
      stid    = H5Tcopy(H5T_C_S1);
      status  = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h

      // We need to first read the head version
      if (H5Aexists(header_id,"HEADER_DEFINITION") > 0)
	{
	  ver4plus=1;
	  inFile->primary[0].hdr_prior4 = 0;	  
	}
      else
	{
	  ver4plus=0;
	  inFile->primary[0].hdr_prior4 = 1;	  
	}
      //      status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->primary);

      
      H5Tinsert(val_tid,"DATE",HOFFSET(sdhdf_primaryHeaderStruct,date),stid);

      // Need to determine whether this is stored as HDR_DEFN or HDR_DEFINITION
      if (ver4plus==1)
	{
	  H5Tinsert(val_tid,"HEADER_DEFINITION",HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn),stid);
	  H5Tinsert(val_tid,"HEADER_DEFINITION_VERSION",HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn_version),stid);

	  H5Tinsert(val_tid,"FILE_FORMAT",HOFFSET(sdhdf_primaryHeaderStruct,file_format),stid);
	  H5Tinsert(val_tid,"FILE_FORMAT_VERSION",HOFFSET(sdhdf_primaryHeaderStruct,file_format_version),stid);
	  H5Tinsert(val_tid,"SCHEDULE_ID",HOFFSET(sdhdf_primaryHeaderStruct,sched_block_id),H5T_NATIVE_INT);
	  H5Tinsert(val_tid,"CALIBRATION_MODE",HOFFSET(sdhdf_primaryHeaderStruct,cal_mode),stid);
	  H5Tinsert(val_tid,"INSTRUMENT",HOFFSET(sdhdf_primaryHeaderStruct,instrument),stid);
	  H5Tinsert(val_tid,"OBSERVER",HOFFSET(sdhdf_primaryHeaderStruct,observer),stid);
	  H5Tinsert(val_tid,"PROJECT_ID",HOFFSET(sdhdf_primaryHeaderStruct,pid),stid);
	  H5Tinsert(val_tid,"RECEIVER",HOFFSET(sdhdf_primaryHeaderStruct,rcvr),stid);
	  
	  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
	  H5Tinsert(val_tid,"TELESCOPE",HOFFSET(sdhdf_primaryHeaderStruct,telescope),stid);
	  H5Tinsert(val_tid,"UTC_START",HOFFSET(sdhdf_primaryHeaderStruct,utc0),stid);
	  H5Tinsert(val_tid,"NUMBER_OF_BEAMS",HOFFSET(sdhdf_primaryHeaderStruct,nbeam),H5T_NATIVE_INT);
	  status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->primary);
	  inFile->nBeam    = inFile->primary[0].nbeam;

	}
      else
	{
	  H5Tinsert(val_tid,"HDR_DEFN",HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn),stid);
	  H5Tinsert(val_tid,"HDR_DEFN_VERSION",HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn_version),stid);
	  H5Tinsert(val_tid,"FILE_FORMAT",HOFFSET(sdhdf_primaryHeaderStruct,file_format),stid);
	  H5Tinsert(val_tid,"FILE_FORMAT_VERSION",HOFFSET(sdhdf_primaryHeaderStruct,file_format_version),stid);
	  H5Tinsert(val_tid,"SCHED_BLOCK_ID",HOFFSET(sdhdf_primaryHeaderStruct,sched_block_id),H5T_NATIVE_INT);
	  H5Tinsert(val_tid,"CAL_MODE",HOFFSET(sdhdf_primaryHeaderStruct,cal_mode),stid);
	  H5Tinsert(val_tid,"INSTRUMENT",HOFFSET(sdhdf_primaryHeaderStruct,instrument),stid);
	  H5Tinsert(val_tid,"OBSERVER",HOFFSET(sdhdf_primaryHeaderStruct,observer),stid);
	  H5Tinsert(val_tid,"PID",HOFFSET(sdhdf_primaryHeaderStruct,pid),stid);
	  H5Tinsert(val_tid,"RECEIVER",HOFFSET(sdhdf_primaryHeaderStruct,rcvr),stid);
	  
	  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
	  H5Tinsert(val_tid,"TELESCOPE",HOFFSET(sdhdf_primaryHeaderStruct,telescope),stid);
	  H5Tinsert(val_tid,"UTC_START",HOFFSET(sdhdf_primaryHeaderStruct,utc0),stid);
	  
	  H5Tinsert(val_tid,"N_BEAMS",HOFFSET(sdhdf_primaryHeaderStruct,nbeam),H5T_NATIVE_INT);
	  status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->primary);
	  inFile->nBeam    = inFile->primary[0].nbeam;

	}
      
      // Load attributes
      for (i=0;i<inFile->nPrimary;i++)
	{
	  inFile->nPrimaryAttributes = sdhdf_getNattributes(inFile,"metadata/primary_header");
	  for (j=0;j<inFile->nPrimaryAttributes;j++)
	    sdhdf_readAttributeFromNum(inFile,"metadata/primary_header",j,&(inFile->primaryAttr[j]));
	}
      status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(header_id);
    }
}


void sdhdf_allocateBeamMemory(sdhdf_fileStruct *inFile,int nbeam)
{
  int i;
  if (inFile->beamAllocatedMemory == 0)
    {
      inFile->beam = (sdhdf_beamStruct *)malloc(sizeof(sdhdf_beamStruct)*nbeam);
      inFile->beamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nbeam);
      inFile->nBeam = nbeam; 
      inFile->beamAllocatedMemory = 1;
      inFile->nBeamHeaderAttributes=0;
      for (i=0;i<nbeam;i++)
	{
	  inFile->beam[i].nBand=0;
	  inFile->beam[i].bandAllocatedMemory=0;
	  inFile->beam[i].nCalBand=0;
	  inFile->beam[i].calBandAllocatedMemory=0;

	  inFile->beam[i].nBandHeaderAttributes=0;
	  inFile->beam[i].nBandDataAttributes=0;
	}
    }
}

//
// Load the beam header meta data in the beam_xx group
//
void sdhdf_loadBeamHeader(sdhdf_fileStruct *inFile)
{
  int   i,j;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  int nbeam;
  char label[MAX_STRLEN];

  nbeam = inFile->nBeam;
 
  if (inFile->beamAllocatedMemory == 0)
    {
      if (nbeam < 0)
	{
	  printf("ERROR: Unexpected number of beams >%d<.\n",nbeam);
	  exit(1);
	}
      inFile->beam       = (sdhdf_beamStruct *)malloc(sizeof(sdhdf_beamStruct)*nbeam);
      inFile->beamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nbeam);
      inFile->nBeam      = nbeam; 
      inFile->beamAllocatedMemory = 1;
    }

  
  if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
    sprintf(label,"metadata/beam_parameters");
  else
    sprintf(label,"metadata/beam_params");
    
  header_id  = H5Dopen2(inFile->fileID,label,H5P_DEFAULT);
  headerT    = H5Dget_type(header_id);
  space      = H5Dget_space(header_id);
  ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

  if (dims[0] != nbeam)
    {
      printf("ERROR: wrong number of beams in beam_params\n");
      printf("HDR version = %s\n",inFile->primary[0].hdr_defn_version);
      printf("nbeam = %d\n",nbeam);
      printf("ndims = %d\n",ndims);
      exit(1);      
    }

  val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_beamHeaderStruct));
  stid    = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(stid,MAX_STRLEN); // Should set to value defined in sdhdf_v1.9.h
  
  H5Tinsert(val_tid,"LABEL",HOFFSET(sdhdf_beamHeaderStruct,label),stid);
  if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
    H5Tinsert(val_tid,"NUMBER_OF_BANDS",HOFFSET(sdhdf_beamHeaderStruct,nBand),H5T_NATIVE_INT);
  else
    H5Tinsert(val_tid,"N_BANDS",HOFFSET(sdhdf_beamHeaderStruct,nBand),H5T_NATIVE_INT);
  H5Tinsert(val_tid,"SOURCE",HOFFSET(sdhdf_beamHeaderStruct,source),stid);
  

  status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,(inFile->beamHeader));
  for (i=0;i<nbeam;i++)
    {  
      inFile->beam[i].bandAllocatedMemory = 0;
      inFile->beam[i].calBandAllocatedMemory = 0;
      inFile->beam[i].nBand = inFile->beamHeader[i].nBand;
    }

  status = H5Tclose(val_tid);
  status = H5Tclose(stid);
  status = H5Dclose(header_id);
  
  // Load attributes

  inFile->nBeamHeaderAttributes = sdhdf_getNattributes(inFile,label);
  for (j=0;j<inFile->nBeamHeaderAttributes;j++)
    sdhdf_readAttributeFromNum(inFile,label,j,&(inFile->beamHeaderAttr[j]));


}


// Set defaults for metadata
void sdhdf_setMetadataDefaults(sdhdf_primaryHeaderStruct *primaryHeader,sdhdf_beamHeaderStruct *beamHeader,
			       sdhdf_bandHeaderStruct *bandHeader,sdhdf_softwareVersionsStruct *softwareVersions,sdhdf_historyStruct *history,
			       int nbeam,int nband)
{
  int i;

  // History
  strcpy(history->date,"UNKNOWN");
  strcpy(history->proc_name,"UNKNOWN");
  strcpy(history->proc_descr,"UNKNOWN");
  strcpy(history->proc_args,"UNKNOWN");
  strcpy(history->proc_host,"UNKNOWN");
  
  // Software versions
  strcpy(softwareVersions->proc_name,"UNKNOWN");
  strcpy(softwareVersions->software,"UNKNONW");
  strcpy(softwareVersions->software_descr,"UNKNOWN");
  strcpy(softwareVersions->software_version,"UNKNOWN");
  
  // Primary header
  strcpy(primaryHeader->date,"UNKNOWN");
  strcpy(primaryHeader->hdr_defn,"UNKNOWN");
  strcpy(primaryHeader->hdr_defn_version,"UNKNOWN");
  strcpy(primaryHeader->file_format,"UNKNOWN");
  strcpy(primaryHeader->file_format_version,"UNKNOWN");
  primaryHeader->sched_block_id=-1;
  strcpy(primaryHeader->cal_mode,"UNKNOWN");
  strcpy(primaryHeader->instrument,"UNKNOWN");
  strcpy(primaryHeader->observer,"UNKNOWN");
  strcpy(primaryHeader->pid,"UNKNOWN");
  strcpy(primaryHeader->rcvr,"UNKNOWN");
  strcpy(primaryHeader->telescope,"UNKNOWN");
  strcpy(primaryHeader->utc0,"UNKNOWN");
  primaryHeader->nbeam=nbeam;

  // Should setup multiple beams correctly
  strcpy(beamHeader->label,"beam0");
  beamHeader->nBand=1;
  strcpy(beamHeader->source,"SOURCE NAME");

  strcpy(bandHeader->label,"band0");
  bandHeader->fc = 0; 
  bandHeader->f0 = 0; 
  bandHeader->f1 = 0; 
  bandHeader->nchan = 0;
  bandHeader->npol = 0;
  strcpy(bandHeader->pol_type,"UNKNOWN");
  bandHeader->dtime = 0; 
  bandHeader->ndump = 0;
 
}

void sdhdf_writeSoftwareVersions(sdhdf_fileStruct *outFile,sdhdf_softwareVersionsStruct *outParams)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];
			   
  dims[0] = 1;

  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_softwareVersionsStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h

  H5Tinsert(datatype_id,"PROC",HOFFSET(sdhdf_softwareVersionsStruct,proc_name),stid);
  H5Tinsert(datatype_id,"SOFTWARE",HOFFSET(sdhdf_softwareVersionsStruct,software),stid);
  H5Tinsert(datatype_id,"SOFTWARE_DESCR",HOFFSET(sdhdf_softwareVersionsStruct,software_descr),stid);
  H5Tinsert(datatype_id,"SOFTWARE_VERSION",HOFFSET(sdhdf_softwareVersionsStruct,software_version),stid);
  
  dataspace_id = H5Screate_simple(1,dims,NULL);

  sprintf(groupName,"metadata");
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(name,"metadata/software_versions");
  if (sdhdf_checkGroupExists(outFile,name) == 1)
      dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  else
      dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
  printf("Created software_versions\n");
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,outParams);
  status  = H5Dclose(dset_id);

  status  = H5Sclose(dataspace_id);  
  status  = H5Tclose(datatype_id);

  status  = H5Tclose(stid);
}


//
// Write the beam header meta data in the beam_xx group
//
void sdhdf_writeBeamHeader(sdhdf_fileStruct *outFile,sdhdf_beamHeaderStruct *beamHeader,int nBeams,char *ver)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];
			   
  dims[0] = nBeams;
  //  printf("outbands = %d\n",nBeams);
  
  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_beamHeaderStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,MAX_STRLEN); 
  printf("Setting size to %d with status %d\n",MAX_STRLEN,status);
  H5Tinsert(datatype_id,"LABEL",HOFFSET(sdhdf_beamHeaderStruct,label),stid);
  if (strcmp(ver,"4.0")==0)
    H5Tinsert(datatype_id,"NUMBER_OF_BANDS",HOFFSET(sdhdf_beamHeaderStruct,nBand),H5T_NATIVE_INT);
  else
    H5Tinsert(datatype_id,"N_BANDS",HOFFSET(sdhdf_beamHeaderStruct,nBand),H5T_NATIVE_INT);
  H5Tinsert(datatype_id,"SOURCE",HOFFSET(sdhdf_beamHeaderStruct,source),stid);
  printf("In beam write with %s\n",beamHeader[0].label);
  dataspace_id = H5Screate_simple(1,dims,NULL);

  sprintf(groupName,"metadata");
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  if (strcmp(ver,"4.0")==0)
    sprintf(name,"metadata/beam_parameters"); 
  else
    sprintf(name,"metadata/beam_params"); 

  if (sdhdf_checkGroupExists(outFile,name) == 1)
    dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  else
    {
      dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
      //      H5Dextend(dset_id,dims);
    }
  //
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,beamHeader);
  //  status  = H5Dwrite(dset_id,datatype_id,dataspace_id,dataspace_id,H5P_DEFAULT,beamHeader);
  status  = H5Dclose(dset_id);

  status  = H5Sclose(dataspace_id);  
  status  = H5Tclose(datatype_id);
  status  = H5Tclose(stid);  
}

//
// Write the beam header meta data in the beam_xx group
//
void sdhdf_writePrimaryHeader(sdhdf_fileStruct *outFile,sdhdf_primaryHeaderStruct *primaryHeader)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];
			   
  dims[0] = 1;

  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_primaryHeaderStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h



  H5Tinsert(datatype_id,"DATE",HOFFSET(sdhdf_primaryHeaderStruct,date),stid);
  H5Tinsert(datatype_id,"HDR_DEFN",HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn),stid);
  H5Tinsert(datatype_id,"HDR_DEFN_VERSION",HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn_version),stid);
  H5Tinsert(datatype_id,"FILE_FORMAT",HOFFSET(sdhdf_primaryHeaderStruct,file_format),stid);
  H5Tinsert(datatype_id,"FILE_FORMAT_VERSION",HOFFSET(sdhdf_primaryHeaderStruct,file_format_version),stid);
  H5Tinsert(datatype_id,"SCHED_BLOCK_ID",HOFFSET(sdhdf_primaryHeaderStruct,sched_block_id),H5T_NATIVE_INT);
  H5Tinsert(datatype_id,"CAL_MODE",HOFFSET(sdhdf_primaryHeaderStruct,cal_mode),stid);
  H5Tinsert(datatype_id,"INSTRUMENT",HOFFSET(sdhdf_primaryHeaderStruct,instrument),stid);
  H5Tinsert(datatype_id,"OBSERVER",HOFFSET(sdhdf_primaryHeaderStruct,observer),stid);
  H5Tinsert(datatype_id,"PID",HOFFSET(sdhdf_primaryHeaderStruct,pid),stid);
  H5Tinsert(datatype_id,"RECEIVER",HOFFSET(sdhdf_primaryHeaderStruct,rcvr),stid);

  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
  H5Tinsert(datatype_id,"TELESCOPE",HOFFSET(sdhdf_primaryHeaderStruct,telescope),stid);
  H5Tinsert(datatype_id,"UTC_START",HOFFSET(sdhdf_primaryHeaderStruct,utc0),stid);
  H5Tinsert(datatype_id,"N_BEAMS",HOFFSET(sdhdf_primaryHeaderStruct,nbeam),H5T_NATIVE_INT);

  
  dataspace_id = H5Screate_simple(1,dims,NULL);

  sprintf(groupName,"metadata");
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(name,"metadata/primary_header");
  if (sdhdf_checkGroupExists(outFile,name) == 1)
      dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  else
      dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
  printf("Created primary_header\n");
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,primaryHeader);
  status  = H5Dclose(dset_id);

  status  = H5Sclose(dataspace_id);  
  status  = H5Tclose(datatype_id);

  status  = H5Tclose(stid);

  
}

//
// Load the band header meta data in the beam_xx/beam_metadata group
// type = 1 for astronomy
// type = 2 for cal
//
void sdhdf_loadBandHeader(sdhdf_fileStruct *inFile,int type)
{
  int   i,j;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  int nbeam,nband;
  char label[MAX_STRLEN];
  char beamLabel[MAX_STRLEN];
  
  nbeam = inFile->nBeam;

  if (inFile->beamAllocatedMemory == 0 || nbeam < 1)
    {
      printf("ERROR: Trying to load bands, but no beams allocated\n");
      exit(1);
    }

  for (i=0;i<nbeam;i++)
    {  
      strcpy(beamLabel,inFile->beamHeader[i].label);
      nband = inFile->beam[i].nBand;
      if (nband < 1)
	{
	  printf("ERROR: No bands allocated\n");
	  exit(1);
 	}

      if (type==1)
	{
	  if (inFile->beam[i].bandAllocatedMemory == 0)
	    {
	      inFile->beam[i].bandData       = (sdhdf_bandStruct *)malloc(sizeof(sdhdf_bandStruct)*nband);
	      inFile->beam[i].bandHeader     = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nband);
	      inFile->beam[i].bandAllocatedMemory = 1;
	    }
	}
      else if (type==2)
	{
	  if (inFile->beam[i].calBandAllocatedMemory == 0)
	    {
	      inFile->beam[i].calBandHeader     = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nband);
	      inFile->beam[i].calBandAllocatedMemory = 1;
	    }
	}
      for (j=0;j<nband;j++)
	{
	  if (type==1)
	    {
	      sdhdf_initialise_bandHeader(&(inFile->beam[i].bandHeader[j]));	    
	      inFile->beam[i].bandData[j].astro_obsHeaderAllocatedMemory = 0;
	      inFile->beam[i].bandData[j].cal_obsHeaderAllocatedMemory   = 0;
	      inFile->beam[i].bandData[j].nAstro_obsHeader               = -1;
	      inFile->beam[i].bandData[j].nCal_obsHeader                 = -1;
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].astro_data));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_on_data));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_off_data));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_proc_tsys));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_proc_diff_gain));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_proc_diff_phase));
	    }
	  else if (type==2)
	    {
	      sdhdf_initialise_bandHeader(&(inFile->beam[i].calBandHeader[j]));	    
	    }
	}

      if (type==1)
	{
	  if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
	    sprintf(label,"%s/metadata/band_parameters",beamLabel);
	  else
	    sprintf(label,"%s/metadata/band_params",beamLabel);
	}
      else if (type==2)
	{
	  if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
	    sprintf(label,"%s/metadata/calibrator_band_parameters",beamLabel);
	  else
	    sprintf(label,"%s/metadata/cal_band_params",beamLabel);
	}
      header_id  = H5Dopen2(inFile->fileID,label,H5P_DEFAULT);
      headerT    = H5Dget_type(header_id);
      space      = H5Dget_space(header_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

      if (dims[0] != nband)
	{
	  printf("ERROR: incorrect number of bands.  In beam header have %d band. Dimension of data is %d\n",nband,(int)dims[0]);
	  printf("ndims = %d\n",ndims);
	  printf("Label = %s\n",label);
	  printf("File = %s\n",inFile->fname);
	  exit(1);
	}

      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_bandHeaderStruct));
      stid = H5Tcopy(H5T_C_S1);
      status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h

      if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
	{    
	  H5Tinsert(val_tid,"LABEL",HOFFSET(sdhdf_bandHeaderStruct,label),stid);
	  H5Tinsert(val_tid,"CENTRE_FREQUENCY",HOFFSET(sdhdf_bandHeaderStruct,fc),H5T_NATIVE_DOUBLE);
	  H5Tinsert(val_tid,"LOW_FREQUENCY",HOFFSET(sdhdf_bandHeaderStruct,f0),H5T_NATIVE_DOUBLE);
	  H5Tinsert(val_tid,"HIGH_FREQUENCY",HOFFSET(sdhdf_bandHeaderStruct,f1),H5T_NATIVE_DOUBLE);
	  H5Tinsert(val_tid,"NUMBER_OF_CHANNELS",HOFFSET(sdhdf_bandHeaderStruct,nchan),H5T_NATIVE_INT);
	  H5Tinsert(val_tid,"NUMBER_OF_POLARISATIONS",HOFFSET(sdhdf_bandHeaderStruct,npol),H5T_NATIVE_INT);
	  status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h
	  H5Tinsert(val_tid,"POLARISATION_TYPE",HOFFSET(sdhdf_bandHeaderStruct,pol_type),stid);
	  H5Tinsert(val_tid,"REQUESTED_INTEGRATION_TIME",HOFFSET(sdhdf_bandHeaderStruct,dtime),H5T_NATIVE_DOUBLE);
	  H5Tinsert(val_tid,"NUMBER_OF_INTEGRATIONS",HOFFSET(sdhdf_bandHeaderStruct,ndump),H5T_NATIVE_INT);
	}
      else
	{
	  H5Tinsert(val_tid,"LABEL",HOFFSET(sdhdf_bandHeaderStruct,label),stid);
	  H5Tinsert(val_tid,"CENTRE_FREQ",HOFFSET(sdhdf_bandHeaderStruct,fc),H5T_NATIVE_DOUBLE);
	  H5Tinsert(val_tid,"LOW_FREQ",HOFFSET(sdhdf_bandHeaderStruct,f0),H5T_NATIVE_DOUBLE);
	  H5Tinsert(val_tid,"HIGH_FREQ",HOFFSET(sdhdf_bandHeaderStruct,f1),H5T_NATIVE_DOUBLE);
	  H5Tinsert(val_tid,"N_CHANS",HOFFSET(sdhdf_bandHeaderStruct,nchan),H5T_NATIVE_INT);
	  H5Tinsert(val_tid,"N_POLS",HOFFSET(sdhdf_bandHeaderStruct,npol),H5T_NATIVE_INT);

	  status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h
	  H5Tinsert(val_tid,"POL_TYPE",HOFFSET(sdhdf_bandHeaderStruct,pol_type),stid);
	  H5Tinsert(val_tid,"DUMP_TIME",HOFFSET(sdhdf_bandHeaderStruct,dtime),H5T_NATIVE_DOUBLE);
	  H5Tinsert(val_tid,"N_DUMPS",HOFFSET(sdhdf_bandHeaderStruct,ndump),H5T_NATIVE_INT);

	}
      if (type==1)
	status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[i].bandHeader);
      else if (type==2)
	status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[i].calBandHeader);  

      
      status  = H5Tclose(val_tid);
      status  = H5Tclose(stid);
      status  = H5Dclose(header_id);
    }
}


void sdhdf_writeBandHeader(sdhdf_fileStruct *outFile,sdhdf_bandHeaderStruct *outBandParams,char *beamLabel,int outBands,int type,char *ver)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];

  // Need to allocate memory first
  //  outFile->beam[ibeam].nBand = outBands;  
  dims[0] = outBands;

  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_bandHeaderStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,12); 

  if (strcmp(ver,"4.0")==0)
    {    
      H5Tinsert(datatype_id,"LABEL",HOFFSET(sdhdf_bandHeaderStruct,label),stid);
      H5Tinsert(datatype_id,"CENTRE_FREQUENCY",HOFFSET(sdhdf_bandHeaderStruct,fc),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"LOW_FREQUENCY",HOFFSET(sdhdf_bandHeaderStruct,f0),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"HIGH_FREQUENCY",HOFFSET(sdhdf_bandHeaderStruct,f1),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"NUMBER_OF_CHANNELS",HOFFSET(sdhdf_bandHeaderStruct,nchan),H5T_NATIVE_INT);
      H5Tinsert(datatype_id,"NUMBER_OF_POLARISATIONS",HOFFSET(sdhdf_bandHeaderStruct,npol),H5T_NATIVE_INT);
      status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(datatype_id,"POLARISATION_TYPE",HOFFSET(sdhdf_bandHeaderStruct,pol_type),stid);
      H5Tinsert(datatype_id,"REQUESTED_INTEGRATION_TIME",HOFFSET(sdhdf_bandHeaderStruct,dtime),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"NUMBER_OF_INTEGRATIONS",HOFFSET(sdhdf_bandHeaderStruct,ndump),H5T_NATIVE_INT);
    }
  else
    {
      H5Tinsert(datatype_id,"LABEL",HOFFSET(sdhdf_bandHeaderStruct,label),stid);
      H5Tinsert(datatype_id,"CENTRE_FREQ",HOFFSET(sdhdf_bandHeaderStruct,fc),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"LOW_FREQ",HOFFSET(sdhdf_bandHeaderStruct,f0),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"HIGH_FREQ",HOFFSET(sdhdf_bandHeaderStruct,f1),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"N_CHANS",HOFFSET(sdhdf_bandHeaderStruct,nchan),H5T_NATIVE_INT);
      H5Tinsert(datatype_id,"N_POLS",HOFFSET(sdhdf_bandHeaderStruct,npol),H5T_NATIVE_INT);
      status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(datatype_id,"POL_TYPE",HOFFSET(sdhdf_bandHeaderStruct,pol_type),stid);
      H5Tinsert(datatype_id,"DUMP_TIME",HOFFSET(sdhdf_bandHeaderStruct,dtime),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"N_DUMPS",HOFFSET(sdhdf_bandHeaderStruct,ndump),H5T_NATIVE_INT);
    }
      dataspace_id = H5Screate_simple(1,dims,NULL);

  // Do we need to create the groups
  sprintf(groupName,"%s",beamLabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  
  sprintf(groupName,"%s/metadata",beamLabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  if (type==1)
    {
      if (strcmp(ver,"4.0")==0)
	sprintf(name,"%s/metadata/band_parameters",beamLabel);
      else
	sprintf(name,"%s/metadata/band_params",beamLabel);
      if (sdhdf_checkGroupExists(outFile,name) == 1)
	{
	  dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	}
      else
	dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
    }
  else
    {
      if (strcmp(ver,"4.0")==0)
	sprintf(name,"%s/metadata/calibrator_band_parameters",beamLabel);
      else
	sprintf(name,"%s/metadata/cal_band_params",beamLabel);
      if (sdhdf_checkGroupExists(outFile,name) == 1)
	dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      else
	dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
    }
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,outBandParams);
  status  = H5Dclose(dset_id);
  status  = H5Sclose(dataspace_id);  
  status  = H5Tclose(datatype_id);
  status  = H5Tclose(stid);
}

// 0 = already exists
int sdhdf_checkGroupExists(sdhdf_fileStruct *inFile,char *groupName)
{
  herr_t status;
  // Check if group exists
  status = H5Eset_auto1(NULL,NULL);
  status = H5Gget_objinfo(inFile->fileID,groupName,0,NULL);
  if (status==0)
    return 0;
  else
    return 1;
}

//
// Load the band observation meta data in the beam_xx/band_yy/band_metadata group
// type = 1: astronomy meta-data
// type = 2: cal meta-data
//
void sdhdf_loadObsHeader(sdhdf_fileStruct *inFile,int type)
{
  int   i,j,k;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  int nbeam,nband,ndump;
  char label[MAX_STRLEN];
  char beamLabel[MAX_STRLEN];
  
  nbeam = inFile->nBeam;

  if (inFile->beamAllocatedMemory == 0 || nbeam < 1)
    {
      printf("ERROR: Trying to load observation metadata, but no beams allocated\n");
      exit(1);
    }

  for (i=0;i<nbeam;i++)
    {  
      strcpy(beamLabel,inFile->beamHeader[i].label);
      nband = inFile->beam[i].nBand;

      if (inFile->beam[i].bandAllocatedMemory == 0 || nband < 1)
	{
	  printf("ERROR: Trying to load observation metadata, but no bands allocated\n");
	  exit(1);
	}

      for (j=0;j<nband;j++)
	{
	  if (type==1)
	    ndump = inFile->beam[i].bandHeader[j].ndump;
	  else if (type==2)
	    ndump = inFile->beam[i].calBandHeader[j].ndump;

	  if (ndump < 1)
	    {
	      printf("ERROR: Trying to load observation metadata, but ndump = 0\n");
	      exit(1);
	    }

	  if (type==1)
	    {
	      if (inFile->beam[i].bandData[j].astro_obsHeaderAllocatedMemory == 0)
		{
		  inFile->beam[i].bandData[j].astro_obsHeader = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
		  inFile->beam[i].bandData[j].astro_obsHeaderAllocatedMemory = 1;
		  inFile->beam[i].bandData[j].nAstro_obsHeader = ndump;
		}
	    }
	  else if (type==2)
	    {
	      if (inFile->beam[i].bandData[j].cal_obsHeaderAllocatedMemory == 0)
		{
		  inFile->beam[i].bandData[j].cal_obsHeader = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
		  inFile->beam[i].bandData[j].cal_obsHeaderAllocatedMemory = 1;
		  inFile->beam[i].bandData[j].nCal_obsHeader = ndump;
		}
	    }

	  for (k=0;k<ndump;k++)
	    {
	      if (type==1)
		sdhdf_initialise_obsHeader(&(inFile->beam[i].bandData[j].astro_obsHeader[k]));
	      else
		sdhdf_initialise_obsHeader(&(inFile->beam[i].bandData[j].cal_obsHeader[k]));
	    }
	 

	  if (type==1)
	    {
	      if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
		sprintf(label,"%s/%s/metadata/observation_parameters",beamLabel,inFile->beam[i].bandHeader[j].label);
	      else
		sprintf(label,"%s/%s/metadata/obs_params",beamLabel,inFile->beam[i].bandHeader[j].label);
	    }
	  else if (type==2)
	    {
	      if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
		sprintf(label,"%s/%s/metadata/calibrator_observation_parameters",beamLabel,inFile->beam[i].bandHeader[j].label);
	      else
		sprintf(label,"%s/%s/metadata/cal_obs_params",beamLabel,inFile->beam[i].bandHeader[j].label);
	    }
	  header_id  = H5Dopen2(inFile->fileID,label,H5P_DEFAULT);
	  headerT    = H5Dget_type(header_id);
	  space      = H5Dget_space(header_id);
	  ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
	  if (ndims == -1)
	    {
	      printf("ERROR: number of dimensions is -1 in %s\n",label);
	      exit(1);
	    }
	  if (dims[0] != ndump)
	    {
	      printf("ERROR: [%s] Missing data detected. In %s number of dumps %d (expected) !== %d (actual) (band = %d) [%s]. ndims = %d, dims[0] = %d, type = %d\n",inFile->fname,inFile->beam[i].bandHeader[j].label,ndump,(int)dims[0],j,label,ndims,(int)dims[0],(int)type);
	      printf("ATTEMPTING TO CONTINUE ....\n");
	    }
	  else
	    {
	      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_obsParamsStruct));
	      stid = H5Tcopy(H5T_C_S1);
	      status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h

	      if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
		{
		  H5Tinsert(val_tid,"ELAPSED_TIME",HOFFSET(sdhdf_obsParamsStruct,timeElapsed),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"INTEGRATION_TIME",HOFFSET(sdhdf_obsParamsStruct,dtime),H5T_NATIVE_DOUBLE);
		  //		  H5Tinsert(val_tid,"TIME_DB",HOFFSET(sdhdf_obsParamsStruct,timedb),stid);
		  H5Tinsert(val_tid,"MJD",HOFFSET(sdhdf_obsParamsStruct,mjd),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"UTC",HOFFSET(sdhdf_obsParamsStruct,utc),stid);
		  //		  H5Tinsert(val_tid,"UT_DATE",HOFFSET(sdhdf_obsParamsStruct,ut_date),stid);
		  H5Tinsert(val_tid,"LOCAL_TIME",HOFFSET(sdhdf_obsParamsStruct,local_time),stid);
		  H5Tinsert(val_tid,"RIGHT_ASCENSION",HOFFSET(sdhdf_obsParamsStruct,raStr),stid);
		  H5Tinsert(val_tid,"DECLINATION",HOFFSET(sdhdf_obsParamsStruct,decStr),stid);



		  //		  H5Tinsert(val_tid,"RA_DEG",HOFFSET(sdhdf_obsParamsStruct,raDeg),H5T_NATIVE_DOUBLE);
		  //		  H5Tinsert(val_tid,"DEC_DEG",HOFFSET(sdhdf_obsParamsStruct,decDeg),H5T_NATIVE_DOUBLE);
		  //		  H5Tinsert(val_tid,"RA_OFFSET",HOFFSET(sdhdf_obsParamsStruct,raOffset),H5T_NATIVE_DOUBLE);
		  //		  H5Tinsert(val_tid,"DEC_OFFSET",HOFFSET(sdhdf_obsParamsStruct,decOffset),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"GALACTIC_LONGITUDE",HOFFSET(sdhdf_obsParamsStruct,gl),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"GALACTIC_LATITUDE",HOFFSET(sdhdf_obsParamsStruct,gb),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"AZIMUTH_ANGLE",HOFFSET(sdhdf_obsParamsStruct,az),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"ZENITH_ANGLE",HOFFSET(sdhdf_obsParamsStruct,ze),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"ELEVATION_ANGLE",HOFFSET(sdhdf_obsParamsStruct,el),H5T_NATIVE_DOUBLE);
		  //		  H5Tinsert(val_tid,"AZ_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,az_drive_rate),H5T_NATIVE_DOUBLE);
		  //		  H5Tinsert(val_tid,"ZE_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,ze_drive_rate),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"HOUR_ANGLE",HOFFSET(sdhdf_obsParamsStruct,hourAngle),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"PARALLACTIC_ANGLE",HOFFSET(sdhdf_obsParamsStruct,paraAngle),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"WIND_DIRECTION",HOFFSET(sdhdf_obsParamsStruct,windDir),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"WIND_SPEED",HOFFSET(sdhdf_obsParamsStruct,windSpd),H5T_NATIVE_DOUBLE);
		}
	      else
		{
		  H5Tinsert(val_tid,"ELAPSED_TIME",HOFFSET(sdhdf_obsParamsStruct,timeElapsed),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"DUMP_TIME",HOFFSET(sdhdf_obsParamsStruct,dtime),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"TIME_DB",HOFFSET(sdhdf_obsParamsStruct,timedb),stid);
		  H5Tinsert(val_tid,"MJD",HOFFSET(sdhdf_obsParamsStruct,mjd),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"UTC",HOFFSET(sdhdf_obsParamsStruct,utc),stid);
		  H5Tinsert(val_tid,"UT_DATE",HOFFSET(sdhdf_obsParamsStruct,ut_date),stid);
		  H5Tinsert(val_tid,"LOCAL_TIME",HOFFSET(sdhdf_obsParamsStruct,local_time),stid);
		  H5Tinsert(val_tid,"RA_STR",HOFFSET(sdhdf_obsParamsStruct,raStr),stid);
		  H5Tinsert(val_tid,"DEC_STR",HOFFSET(sdhdf_obsParamsStruct,decStr),stid);
		  H5Tinsert(val_tid,"RA_DEG",HOFFSET(sdhdf_obsParamsStruct,raDeg),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"DEC_DEG",HOFFSET(sdhdf_obsParamsStruct,decDeg),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"RA_OFFSET",HOFFSET(sdhdf_obsParamsStruct,raOffset),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"DEC_OFFSET",HOFFSET(sdhdf_obsParamsStruct,decOffset),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"GL",HOFFSET(sdhdf_obsParamsStruct,gl),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"GB",HOFFSET(sdhdf_obsParamsStruct,gb),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"AZ",HOFFSET(sdhdf_obsParamsStruct,az),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"ZE",HOFFSET(sdhdf_obsParamsStruct,ze),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"EL",HOFFSET(sdhdf_obsParamsStruct,el),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"AZ_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,az_drive_rate),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"ZE_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,ze_drive_rate),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"HOUR_ANGLE",HOFFSET(sdhdf_obsParamsStruct,hourAngle),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"PARA_ANGLE",HOFFSET(sdhdf_obsParamsStruct,paraAngle),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"WIND_DIR",HOFFSET(sdhdf_obsParamsStruct,windDir),H5T_NATIVE_DOUBLE);
		  H5Tinsert(val_tid,"WIND_SPD",HOFFSET(sdhdf_obsParamsStruct,windSpd),H5T_NATIVE_DOUBLE);
		}
	      
	      if (type==1)
		status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[i].bandData[j].astro_obsHeader);
	      else if (type==2)
		status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[i].bandData[j].cal_obsHeader);

	      if (strcmp(inFile->primary[0].hdr_defn_version,"4.0")==0)
		{
		  int l;
		  // In this version we need to calculate RA_DEG and DEC_DEG
		  if (type==1)
		    {
		      for (l=0;l<inFile->beam[i].bandData[j].nAstro_obsHeader;l++)
			{
			  inFile->beam[i].bandData[j].astro_obsHeader[l].raDeg = hms_turn(inFile->beam[i].bandData[j].astro_obsHeader->raStr)*360;
			  inFile->beam[i].bandData[j].astro_obsHeader[l].decDeg = dms_turn(inFile->beam[i].bandData[j].astro_obsHeader->decStr)*360;
			}
		    }
		  else if (type==2)
		    {
		      for (l=0;l<inFile->beam[i].bandData[j].nCal_obsHeader;l++)
			{
			  inFile->beam[i].bandData[j].cal_obsHeader[l].raDeg = hms_turn(inFile->beam[i].bandData[j].cal_obsHeader->raStr)*360;
			  inFile->beam[i].bandData[j].cal_obsHeader[l].decDeg = dms_turn(inFile->beam[i].bandData[j].cal_obsHeader->decStr)*360;
			}

		    }
		  
		}

	      status  = H5Tclose(val_tid);
	      status  = H5Tclose(stid);
	    }
	  status  = H5Dclose(header_id);
	}
    }
}

void sdhdf_initialise_spectralDumps(sdhdf_spectralDumpsStruct *in)
{
  in->freqAllocatedMemory=0;
  in->flagAllocatedMemory=0;
  in->dataWeightsAllocatedMemory=0;
  in->pol1AllocatedMemory=0;
  in->pol2AllocatedMemory=0;
  in->pol3AllocatedMemory=0;
  in->pol4AllocatedMemory=0;
  in->nFreqAttributes=0;
  in->nDataAttributes=0;
}

void sdhdf_initialise_bandHeader(sdhdf_bandHeaderStruct *header)
{
  strcpy(header->label,"unset");
  header->fc = -1;
  header->f0 = -1;
  header->f1 = -1;
  header->nchan = -1;
  header->npol = -1;
  header->ndump = -1;
}

void sdhdf_initialise_obsHeader(sdhdf_obsParamsStruct *obs)
{
  obs->timeElapsed=-1;
  obs->dtime = -1;
  strcpy(obs->timedb,"unset");
  obs->mjd = -1;
  strcpy(obs->utc,"unset");
  strcpy(obs->ut_date,"unset");
  strcpy(obs->local_time,"unset");
  strcpy(obs->raStr,"unset");
  strcpy(obs->decStr,"unset");
  obs->raDeg = -1;
  obs->decDeg = -1;
  obs->raOffset = -1;
  obs->decOffset = -1;
  obs->gl = -1;
  obs->gb = -1;
  obs->az = -1;
  obs->ze = -1;
  obs->el = -1;
  obs->az_drive_rate = -1;
  obs->ze_drive_rate = -1;
  obs->hourAngle = -1;
  obs->paraAngle = -1;
  obs->windDir = -1;
  obs->windSpd = -1;
}

void sdhdf_writeHistory(sdhdf_fileStruct *outFile,sdhdf_historyStruct *outParams,int n)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char groupName[1024];
  
  dims[0] = n;
  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_historyStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,64); 
  dataspace_id = H5Screate_simple(1,dims,NULL);
  status = H5Tset_size(stid,20); H5Tinsert(datatype_id,"DATE",HOFFSET(sdhdf_historyStruct,date),stid);
  status = H5Tset_size(stid,64); H5Tinsert(datatype_id,"PROC",HOFFSET(sdhdf_historyStruct,proc_name),stid);
  status = H5Tset_size(stid,64); H5Tinsert(datatype_id,"PROC_DESCR",HOFFSET(sdhdf_historyStruct,proc_descr),stid);
  status = H5Tset_size(stid,64); H5Tinsert(datatype_id,"PROC_ARGS",HOFFSET(sdhdf_historyStruct,proc_args),stid);
  status = H5Tset_size(stid,64); H5Tinsert(datatype_id,"PROC_HOST",HOFFSET(sdhdf_historyStruct,proc_host),stid);

  sprintf(groupName,"/metadata");
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }

  if (sdhdf_checkGroupExists(outFile,"metadata/history") == 0)    
    {
      printf("THE HISTORY TABLE ALREADY EXISTS\n");
      dset_id  = H5Dopen2(outFile->fileID,"metadata/history",H5P_DEFAULT);
    }
  else
    dset_id = H5Dcreate2(outFile->fileID,"/metadata/history",datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,outParams);
  status  = H5Dclose(dset_id);
  status  = H5Sclose(dataspace_id);  
  status  = H5Tclose(datatype_id);
  status  = H5Tclose(stid);
}

//
// Load the history information
//
void sdhdf_loadHistory(sdhdf_fileStruct *inFile)
{
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];

  if (sdhdf_checkGroupExists(inFile,"metadata/history") == 1)
    {
      printf("WARNING: no history information\n");
      inFile->nHistory = 0;
      inFile->historyAllocatedMemory = 0;
    }
  else
    {
      header_id  = H5Dopen2(inFile->fileID,"metadata/history",H5P_DEFAULT);
      headerT    = H5Dget_type(header_id);
      space      = H5Dget_space(header_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
      
      if (inFile->historyAllocatedMemory == 0)
	{
	  inFile->history = (sdhdf_historyStruct *)malloc(sizeof(sdhdf_historyStruct)*MAX_HISTORY);
      if (dims[0] > MAX_HISTORY)
	{
	  printf("ERROR: need to increase MAX_HISTORY\n");
	  exit(1);
	}
      inFile->nHistory = dims[0]; 
      inFile->historyAllocatedMemory = 1;
	}
      
      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_historyStruct));
      stid = H5Tcopy(H5T_C_S1);
      
      status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,"DATE",HOFFSET(sdhdf_historyStruct,date),stid);
      status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,"PROC",HOFFSET(sdhdf_historyStruct,proc_name),stid);
      H5Tinsert(val_tid,"PROC_DESCR",HOFFSET(sdhdf_historyStruct,proc_descr),stid);
      H5Tinsert(val_tid,"PROC_ARGS",HOFFSET(sdhdf_historyStruct,proc_args),stid);
      H5Tinsert(val_tid,"PROC_HOST",HOFFSET(sdhdf_historyStruct,proc_host),stid);
      
      status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->history);
      
      status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(header_id);
    }
}

void sdhdf_addHistory(sdhdf_historyStruct *history,int n,char *procName,char *descr,char *args)
{
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  char timeDate[1024];
  char host[1024];
  
  struct utsname unameData;
  uname(&unameData);
  strcpy(host,unameData.nodename);
  
  sprintf(timeDate,"%d-%02d-%02d-%02d:%02d:%02d",tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
  strcpy(history[n].date,timeDate);
  strcpy(history[n].proc_name,procName);
  strcpy(history[n].proc_descr,descr);
  strcpy(history[n].proc_args,args);
  strcpy(history[n].proc_host,host);  
}



//
// Load the software information
//
void sdhdf_loadSoftware(sdhdf_fileStruct *inFile)
{
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];

  header_id  = H5Dopen2(inFile->fileID,"metadata/software_versions",H5P_DEFAULT);
  headerT    = H5Dget_type(header_id);
  space      = H5Dget_space(header_id);
  ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

  if (inFile->softwareAllocatedMemory == 0)
    {
      inFile->software = (sdhdf_softwareVersionsStruct *)malloc(sizeof(sdhdf_softwareVersionsStruct)*dims[0]);
      inFile->nSoftware = dims[0]; 
      inFile->softwareAllocatedMemory = 1;
    }

  val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_softwareVersionsStruct));
  stid = H5Tcopy(H5T_C_S1);

  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
  H5Tinsert(val_tid,"PROC",HOFFSET(sdhdf_softwareVersionsStruct,proc_name),stid);
  H5Tinsert(val_tid,"SOFTWARE",HOFFSET(sdhdf_softwareVersionsStruct,software),stid);
  H5Tinsert(val_tid,"SOFTWARE_DESCR",HOFFSET(sdhdf_softwareVersionsStruct,software_descr),stid);
  H5Tinsert(val_tid,"SOFTWARE_VERSION",HOFFSET(sdhdf_softwareVersionsStruct,software_version),stid);

  status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->software);
    
  status = H5Tclose(val_tid);
  status = H5Tclose(stid);
  status = H5Dclose(header_id);
}

int sdhdf_getNattributes(sdhdf_fileStruct *inFile,char *dataName)
{
  hid_t dataset_id;
  herr_t status;
  H5O_info_t object_info;
  //  printf("VERSION:   H5_VERS_RELEASE = %d %d\n",H5_VERS_MAJOR,H5_VERS_MINOR);
  dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
  #if H5_VERS_MINOR == 10
   status = H5Oget_info(dataset_id,&object_info);
  #else
   status = H5Oget_info(dataset_id,&object_info,H5O_INFO_NUM_ATTRS);
  #endif
  status = H5Dclose(dataset_id);
  return object_info.num_attrs;  
}

void sdhdf_readAttributeFromNum(sdhdf_fileStruct *inFile,char *dataName,int num,sdhdf_attributes_struct *attribute)
{
  hid_t dataset_id,attr_id,atype,atype_mem,aspace,space;
  H5T_class_t type_class;
  hid_t type,stid;
  herr_t status;
  int rank;
  hsize_t dims[5];
  int ndims;


  // Open the attribute based on its number
  attr_id      = H5Aopen_by_idx(inFile->fileID,dataName, H5_INDEX_NAME, H5_ITER_NATIVE,num,H5P_DEFAULT,H5P_DEFAULT);
  
  if (attr_id < 0)  // Unable to open attribute
    {
      strcpy(attribute->key,"uncertain");
      strcpy(attribute->value,"not defined");
      attribute->attributeType=0;
      return;
    }

  // Find out if this is a multi-dimensional array
  aspace       = H5Aget_space(attr_id);
  ndims        = H5Sget_simple_extent_dims(aspace,dims,NULL);
  //  printf("ndims = %d, aspace = %d, dims[0] = %d\n",ndims,aspace,dims[0]);
  if (ndims != 0) // 1 && ndims != 0)
    {
      hid_t dataspace_id,s1_tid;
      hsize_t storeSize;
      char *buf,*buffer;
      H5A_info_t attribute_info;
      sdhdf_attributes_struct attrStruct;
      char name[1024];
      //      printf("A\n");

      
      H5Aget_name(attr_id,MAX_STRLEN,name);
      //      printf("HERE with attribute->name = >%s<\n",name);

      // GEORGE: WORKING HERE ON THE ATTRIBUTES
      //      strcpy(attribute->value,"NOT SET");
      //      attribute->attributeType=0;
      //      printf("Cannot read attributes with multi dimensions >%s< >%s<\n",attribute->key,inFile->fname);

      buffer = (char *)malloc(sizeof(char)*MAX_STRLEN);
      
      buf = (char *)malloc(sizeof(char)*MAX_STRLEN);
      //      printf("ndims > 0.  ndims = %d (%s), dims[0] = %d\n",ndims,dataName,dims[0]);
      
      //      dataspace_id = H5Screate_simple(1,dims,NULL);
      //      printf("dataspace\n");
      //      storeSize = H5Aget_storage_size(attr_id);
      //      printf("Storage size - %d\n",storeSize);
      //      atype        = H5Aget_type(attr_id); // Retrieves a copy of the datatype for an attribute
      atype_mem = H5Tcopy(H5T_C_S1);
      status = H5Tset_size(atype_mem, MAX_STRLEN);
      s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(sdhdf_attributes_struct));
      
      //      H5Tinsert(s1_tid, "name", HOFFSET(sdhdf_attributes_struct, name), atype_mem);
      //      strcpy(attribute->name,attribute->key);
      H5Tinsert(s1_tid, "description", HOFFSET(sdhdf_attributes_struct, key), atype_mem);
      H5Tinsert(s1_tid, "unit", HOFFSET(sdhdf_attributes_struct, unit), atype_mem);
      H5Tinsert(s1_tid, "default", HOFFSET(sdhdf_attributes_struct, value), atype_mem);
      status = H5Aread(attr_id, s1_tid, attribute);
      strcpy(attribute->name,name);
      //      printf("Obtained >%s< >%s<\n",attribute->value,attribute->name);
      //      	strcpy(attribute->name, attr_name);
      //	printf("Obtained2 %s %s\n",attribute->value,attribute->name);
      //      printf("A status = %d\n",status);
      //      printf("... Attribute_encoding = %d:   H5T_CSET_UTF8 = %d  H5T_CSET_ASCII = %d\n",H5Tget_cset(atype_mem),H5T_CSET_UTF8,H5T_CSET_ASCII);
      //      strcpy(buffer,"");
      
      //      status = H5Aget_info(attr_id,&attribute_info);
      //      printf("Data size = %d\n",attribute_info.data_size);
      //      status = H5Tset_cset(atype_mem,H5T_CSET_ASCII);		  
		  //		  printf("Size = %d\n",attribute_info.data_size);
      //      status = H5Tset_size(atype_mem,attribute_info.data_sioze);		  
      //      printf("Setting size: status = %d\n",status);
      //      status = H5Aread(attr_id,atype_mem,buffer);
      //      printf("Buffer = %s, status = %d\n",buffer,status);
      // FIX ME: somehow it is not reading in the currect length
      //      buffer[attribute_info.data_size]='\0';
      
      //      strcpy(attribute->value,buffer);
      
      /*
      //      status    = H5Tset_cset(atype_mem,H5T_CSET_UTF8);
      status = H5Tset_cset(atype_mem,H5T_CSET_ASCII);		  
      //      status    = H5Tset_size(atype_mem,H5T_VARIABLE);
      status    = H5Tset_size(atype_mem,64);
      printf("B status = %d\n",status);
      //      status = H5Tset_size(atype_mem,storeSize);		  
      //      H5Tset_size (atype_mem, H5T_VARIABLE);
      status = H5Aread(attr_id,atype_mem,&buf);
      printf("LOADED >%s< %d\n",buf,status);
      */
      status = H5Sclose(dataspace_id);
      free(buf);
      free(buffer);
      // For now not read attributes like DIMENSION_LABEL and DIMENSION_LIST      
    }
  else // Here the number of dimensions is zero as storing a scalar
    {
      atype        = H5Aget_type(attr_id); // Retrieves a copy of the datatype for an attribute
      type_class   = H5Tget_class(atype);

      if (type_class == H5T_STRING)        // Is the attribute a string
	{
	  char string_out[MAX_STRLEN];
	  char *buffer;
	  buffer = (char *)malloc(sizeof(char)*MAX_STRLEN);
	    
	  H5Aget_name(attr_id,MAX_STRLEN,attribute->key);
	  //	  atype_mem = H5Tget_native_type(atype,H5T_DIR_ASCEND);
	  //	  printf("Attribute key = %s\n",attribute->key);
	  attribute->attributeType=0;
	  if (strcmp(attribute->key,"CLASS")==0)
	    {
	      //	      printf("Reading scalar CLASS\n");
	      strcpy(attribute->value,"FIX ME 2");
	    }
	  else 
	    {
	      atype_mem = H5Tcopy(H5T_C_S1);
	      if (H5Tget_cset(atype) == H5T_CSET_UTF8) // UTF-8 STRING
		{
		  status    = H5Tset_cset(atype_mem,H5T_CSET_UTF8);
		  status    = H5Tset_size(atype_mem,H5T_VARIABLE);

		  status = H5Aread(attr_id,atype_mem,&buffer);
		  strcpy(attribute->value,buffer);
		}
	      else // ASCII STRING
		{
		  H5A_info_t attribute_info;
		  //		  printf("Reading an ASCII string\n");
		  strcpy(buffer,"");

		  status = H5Aget_info(attr_id,&attribute_info);
		  status = H5Tset_cset(atype_mem,H5T_CSET_ASCII);		  
		  //		  printf("Size = %d\n",attribute_info.data_size);
		  status = H5Tset_size(atype_mem,attribute_info.data_size);		  
		  //		  printf("Setting size: status = %d\n",status);
		  status = H5Aread(attr_id,atype_mem,buffer);
		  //		  printf("Buffer = %s, status = %d\n",buffer,status);
		  // FIX ME: somehow it is not reading in the currect length
		  buffer[attribute_info.data_size]='\0';
		  
		  strcpy(attribute->value,buffer);


		}
	      
	      status = H5Tclose(atype_mem);	      
	      //	      printf("status at close = %d\n",status);
	    }
	  free(buffer);
	}
      else if (type_class == H5T_FLOAT)
	{
	  double fval;
	  attribute->attributeType=1;
	  H5Aget_name(attr_id,MAX_STRLEN,attribute->key);
	  status  = H5Aread(attr_id, atype, &fval);	  
	  attribute->fvalue = fval;
	}
      else if (type_class == H5T_INTEGER)
	{
	  printf("HAVE AN INTEGER TYPE ATTRIBUTE. ERROR: DON'T KNOW WHAT TO DO\n");
	  attribute->attributeType=2;
	  exit(1);
	}
      else
	{
	  printf("CANNOT READ ATTRIBUTE\n");
	  strcpy(attribute->key,"unset");
	  strcpy(attribute->value,"unset");
	  attribute->attributeType=1;
	}
      status = H5Tclose(atype);
    }
  //  printf("Closing off\n");
  status = H5Sclose(aspace);
  status = H5Aclose(attr_id);

}

void sdhdf_copyBandHeaderStruct(sdhdf_bandHeaderStruct *in,sdhdf_bandHeaderStruct *out,int n)
{
  int i;
  for (i=0;i<n;i++)
    {
      strcpy(out[i].label,in[i].label);
      out[i].fc = in[i].fc;
      out[i].f0 = in[i].f0;
      out[i].f1 = in[i].f1;
      out[i].nchan = in[i].nchan;
      out[i].npol = in[i].npol;
      strcpy(out[i].pol_type,in[i].pol_type);
      out[i].dtime = in[i].dtime;
      out[i].ndump = in[i].ndump;
    }
  
}

void sdhdf_writeObsParams(sdhdf_fileStruct *outFile,char *bandLabel,char *beamLabel,int iband,sdhdf_obsParamsStruct *obsParams,int ndump,int type,char *ver)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];

  // Need to allocate memory if I need the next line
  //  outFile->beam[ibeam].bandData[iband].nAstro_obsHeader = ndump;  
  dims[0] = ndump;
  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_obsParamsStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h

  if (strcmp(ver,"4.0")==0)
    {
      H5Tinsert(datatype_id,"ELAPSED_TIME",HOFFSET(sdhdf_obsParamsStruct,timeElapsed),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"INTEGRATION_TIME",HOFFSET(sdhdf_obsParamsStruct,dtime),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"MJD",HOFFSET(sdhdf_obsParamsStruct,mjd),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"UTC",HOFFSET(sdhdf_obsParamsStruct,utc),stid);
      H5Tinsert(datatype_id,"LOCAL_TIME",HOFFSET(sdhdf_obsParamsStruct,local_time),stid);
      H5Tinsert(datatype_id,"RIGHT_ASCENSION",HOFFSET(sdhdf_obsParamsStruct,raStr),stid);
      H5Tinsert(datatype_id,"DECLINATION",HOFFSET(sdhdf_obsParamsStruct,decStr),stid);
      H5Tinsert(datatype_id,"GALACTIC_LONGITUDE",HOFFSET(sdhdf_obsParamsStruct,gl),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"GALACTIC_LATITUDE",HOFFSET(sdhdf_obsParamsStruct,gb),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"AZIMUTH_ANGLE",HOFFSET(sdhdf_obsParamsStruct,az),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"ZENITH_ANGLE",HOFFSET(sdhdf_obsParamsStruct,ze),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"ELEVATION_ANGLE",HOFFSET(sdhdf_obsParamsStruct,el),H5T_NATIVE_DOUBLE);
      //      H5Tinsert(datatype_id,"AZ_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,az_drive_rate),H5T_NATIVE_DOUBLE);
      //      H5Tinsert(datatype_id,"ZE_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,ze_drive_rate),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"HOUR_ANGLE",HOFFSET(sdhdf_obsParamsStruct,hourAngle),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"PARALLACTIC_ANGLE",HOFFSET(sdhdf_obsParamsStruct,paraAngle),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"WIND_DIRECTION",HOFFSET(sdhdf_obsParamsStruct,windDir),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"WIND_SPEED",HOFFSET(sdhdf_obsParamsStruct,windSpd),H5T_NATIVE_DOUBLE);

    }
  else
    {
      H5Tinsert(datatype_id,"ELAPSED_TIME",HOFFSET(sdhdf_obsParamsStruct,timeElapsed),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"DUMP_TIME",HOFFSET(sdhdf_obsParamsStruct,dtime),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"TIME_DB",HOFFSET(sdhdf_obsParamsStruct,timedb),stid);
      H5Tinsert(datatype_id,"MJD",HOFFSET(sdhdf_obsParamsStruct,mjd),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"UTC",HOFFSET(sdhdf_obsParamsStruct,utc),stid);
      H5Tinsert(datatype_id,"UT_DATE",HOFFSET(sdhdf_obsParamsStruct,ut_date),stid);
      H5Tinsert(datatype_id,"LOCAL_TIME",HOFFSET(sdhdf_obsParamsStruct,local_time),stid);
      H5Tinsert(datatype_id,"RA_STR",HOFFSET(sdhdf_obsParamsStruct,raStr),stid);
      H5Tinsert(datatype_id,"DEC_STR",HOFFSET(sdhdf_obsParamsStruct,decStr),stid);
      H5Tinsert(datatype_id,"RA_DEG",HOFFSET(sdhdf_obsParamsStruct,raDeg),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"DEC_DEG",HOFFSET(sdhdf_obsParamsStruct,decDeg),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"RA_OFFSET",HOFFSET(sdhdf_obsParamsStruct,raOffset),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"DEC_OFFSET",HOFFSET(sdhdf_obsParamsStruct,decOffset),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"GL",HOFFSET(sdhdf_obsParamsStruct,gl),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"GB",HOFFSET(sdhdf_obsParamsStruct,gb),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"AZ",HOFFSET(sdhdf_obsParamsStruct,az),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"ZE",HOFFSET(sdhdf_obsParamsStruct,ze),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"EL",HOFFSET(sdhdf_obsParamsStruct,el),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"AZ_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,az_drive_rate),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"ZE_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,ze_drive_rate),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"PARA_ANGLE",HOFFSET(sdhdf_obsParamsStruct,paraAngle),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"WIND_DIR",HOFFSET(sdhdf_obsParamsStruct,windDir),H5T_NATIVE_DOUBLE);
      H5Tinsert(datatype_id,"WIND_SPD",HOFFSET(sdhdf_obsParamsStruct,windSpd),H5T_NATIVE_DOUBLE);
    }

  dataspace_id = H5Screate_simple(1,dims,NULL);
  // Do we need to create the groups

  sprintf(groupName,"%s",beamLabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(groupName,"%s/%s",beamLabel,bandLabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(groupName,"%s/%s/metadata",beamLabel,bandLabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  if (type==1)
    {
      if (strcmp(ver,"4.0")==0)
	sprintf(name,"%s/%s/metadata/observation_parameters",beamLabel,bandLabel);
      else
	sprintf(name,"%s/%s/metadata/obs_params",beamLabel,bandLabel);
    }
  else if (type==2)
    {
      if (strcmp(ver,"4.0")==0)
	sprintf(name,"%s/%s/metadata/calibration_observation_parameters",beamLabel,bandLabel);
      else
	sprintf(name,"%s/%s/metadata/cal_obs_params",beamLabel,bandLabel);
    }
      dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,obsParams);
  status  = H5Dclose(dset_id);

  status  = H5Sclose(dataspace_id);  
  status  = H5Tclose(datatype_id);

  status  = H5Tclose(stid);
}

void sdhdf_copySingleObsParams(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,sdhdf_obsParamsStruct *obsParam)
{
  obsParam->timeElapsed = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].timeElapsed;
  obsParam->dtime = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].dtime;
  strcpy(obsParam->timedb,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].timedb);
  obsParam->mjd = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].mjd;
  strcpy(obsParam->utc,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].utc);
  strcpy(obsParam->ut_date,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].ut_date);
  strcpy(obsParam->local_time,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].local_time);
  strcpy(obsParam->raStr,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].raStr);
  strcpy(obsParam->decStr,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].decStr);
  obsParam->raDeg = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].raDeg;
  obsParam->decDeg = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].decDeg;
  obsParam->raOffset = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].raOffset;
  obsParam->decOffset = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].decOffset;
  obsParam->gl = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].gl;
  obsParam->gb = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].gb;
  obsParam->az = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].az;
  obsParam->ze = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].ze;
  obsParam->el = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].el;
  obsParam->az_drive_rate = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].az_drive_rate;
  obsParam->ze_drive_rate = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].ze_drive_rate;
  obsParam->hourAngle = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].hourAngle;
  obsParam->paraAngle = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].paraAngle;
  obsParam->windDir = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].windDir;
  obsParam->windSpd = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].windSpd;
}

void sdhdf_copySingleObsParamsCal(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,sdhdf_obsParamsStruct *obsParam)
{
  obsParam->timeElapsed = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].timeElapsed;
  obsParam->dtime = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].dtime;
  strcpy(obsParam->timedb,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].timedb);
  obsParam->mjd = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].mjd;
  strcpy(obsParam->utc,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].utc);
  strcpy(obsParam->ut_date,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].ut_date);
  strcpy(obsParam->local_time,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].local_time);
  strcpy(obsParam->raStr,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].raStr);
  strcpy(obsParam->decStr,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].decStr);
  obsParam->raDeg = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].raDeg;
  obsParam->decDeg = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].decDeg;
  obsParam->raOffset = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].raOffset;
  obsParam->decOffset = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].decOffset;
  obsParam->gl = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].gl;
  obsParam->gb = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].gb;
  obsParam->az = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].az;
  obsParam->ze = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].ze;
  obsParam->el = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].el;
  obsParam->az_drive_rate = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].az_drive_rate;
  obsParam->ze_drive_rate = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].ze_drive_rate;
  obsParam->hourAngle = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].hourAngle;
  obsParam->paraAngle = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].paraAngle;
  obsParam->windDir = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].windDir;
  obsParam->windSpd = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].windSpd;
}

void sdhdf_copyAttributes(sdhdf_attributes_struct *in,int n_in,sdhdf_attributes_struct *out,int *n_out)
{
  int i;
  *n_out = n_in;
  for (i=0;i<n_in;i++)
    {
      //      printf("In sdhdf_copyAttributes: %s %s\n",in[i].key,in[i].value);
      strcpy(out[i].key,in[i].key);
      out[i].attributeType = in[i].attributeType;
      if (in[i].attributeType==0)
	strcpy(out[i].value,in[i].value);
      else if (in[i].attributeType==1)
	out[i].fvalue = in[i].fvalue;
      else if (in[i].attributeType==2)
	out[i].ivalue = in[i].ivalue;
    }
}


void sdhdf_writeAttribute(sdhdf_fileStruct *outFile,char *dataName,sdhdf_attributes_struct *attribute) //,char *attrName,char *result)
{
  hid_t memtype,dataset_id,attr,space,stid;
  herr_t status;
  hsize_t dims[1] = {1}; // Has only 1 string
  char *result;

  result = attribute->value;
  //  printf("In writeAttribute with dataName = >%s<, attrName = >%s< and result = >%s<\n",dataName,attribute->key,attribute->value);
  //  printf("AttributeType = %d\n",attribute->attributeType);
  // Python strings are written as variable length strings 


  dataset_id   = H5Dopen2(outFile->fileID,dataName,H5P_DEFAULT);

  // This produces a multi dimensional attribute - which is not the same as the python code
  /*
  stid    = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(stid,H5T_VARIABLE);     
  space   = H5Screate_simple(1,dims,NULL);
  attr    = H5Acreate(dataset_id,attrName,stid,space,H5P_DEFAULT,H5P_DEFAULT);
  printf("Trying to write >%s<\n",result);
  status = H5Awrite(attr,stid,&result);
  printf("The status after writing is: %d\n",status);
  H5Dclose(dataset_id);
  H5Aclose(attr);
  H5Tclose(memtype);
  */

  // This produces an attribute which is a string
  if (attribute->attributeType==0)
    {
      stid    = H5Tcopy(H5T_C_S1);
      //  status  = H5Tset_size(stid,H5T_VARIABLE);
      status  = H5Tset_size(stid,strlen(result));     
      space   = H5Screate(H5S_SCALAR);
      attr    = H5Acreate(dataset_id,attribute->key,stid,space,H5P_DEFAULT,H5P_DEFAULT);
      //      printf("writeAttribute: Trying to write >%s<\n",result);
      status  = H5Awrite(attr,stid,result);
      //      printf("writeAttribute: The status after writing is: %d\n",status);
      H5Dclose(dataset_id);
      H5Tclose(memtype);
      
      H5Aclose(attr);
    }
  else if (attribute->attributeType==1)
    {
      float val = 1.0;
      int ival= 1;
      
      //      stid    = H5Tcopy(H5T_C_S1);
      //  status  = H5Tset_size(stid,H5T_VARIABLE);
      //      status  = H5Tset_size(stid,strlen(result));     
      //      space   = H5Screate(H5S_SCALAR);

      space = H5Screate(H5S_SCALAR);

      //      attr    = H5Acreate(dataset_id,attribute->key,stid,space,H5P_DEFAULT,H5P_DEFAULT);
      attr    = H5Acreate(dataset_id,attribute->key,H5T_NATIVE_FLOAT,space,H5P_DEFAULT,H5P_DEFAULT);

      //      printf("writeAttribute: Trying to write >%s<\n",result);
      status  = H5Awrite(attr,H5T_NATIVE_FLOAT,&(attribute->fvalue));
      //      status  = H5Awrite(attr,stid,result);
      //      printf("writeAttribute: The status after writing is: %d\n",status);

      H5Dclose(dataset_id);
      H5Tclose(memtype);

      H5Aclose(attr);

      /*
      printf("GEORGE: WRITING A FLOATING VALUE INTO %s\n",attribute->key);
      space   = H5Screate(H5S_SCALAR);
      printf("GEORGE A %d\n",space);
      //      attr    = H5Acreate(dataset_id,attribute->key,H5T_NATIVE_FLOAT,space,H5P_DEFAULT,H5P_DEFAULT);
      //      attr    = H5Acreate(dataset_id,attribute->key,H5T_NATIVE_FLOAT,space,H5P_DEFAULT,H5P_DEFAULT);
      //      attr    = H5Acreate(dataset_id,attribute->key,H5T_FLOAT,space,H5P_DEFAULT,H5P_DEFAULT);
      attr    = H5Acreate(dataset_id,attribute->key,H5T_NATIVE_INT,space,H5P_DEFAULT,H5P_DEFAULT);
      printf("GEORGE B %f\n",attribute->fvalue);
      //      status  = H5Awrite(attr,H5T_NATIVE_FLOAT,&val);
      //      status  = H5Awrite(attr,H5T_FLOAT,&val);
      status  = H5Awrite(attr,H5T_NATIVE_INT,&ival);
      printf("GEORGE C >%d<\n",status);
      H5Aclose(attr);
      printf("COMPLETE WRITE WITH STATUS: %s\n",status);
      */
    }

  //  memtype = H5Tcopy(H5T_C_S1);
  //  status  = H5Tset_size(memtype,MAX_STRLEN);
  //  space   = H5Screate(H5S_SCALAR);
  
  //  status = H5Tset_size(memtype,H5T_VARIABLE);


  //  
  //  space = H5Screate_simple(1,dims,NULL);
  //  space = H5Screate(H5S_SCALAR);

  
  //  attr = H5Acreate(dataset_id,attrName,memtype,space,H5P_DEFAULT,H5P_DEFAULT);
  // attr = H5Acreate(dataset_id,attrName,memtype,space,H5P_DEFAULT,H5P_DEFAULT);

  //  status = H5Awrite(attr,space,result);
  

}

void sdhdf_loadPersistentRFI(sdhdf_rfi *rfi,int *nRFI,int maxRFI,char *tel)
{
  FILE *fin;
  char fname[MAX_STRLEN];
  char runtimeDir[MAX_STRLEN];
  char str[4096];
  char dirName[1024];

  *nRFI = 0;
  sdhdf_getTelescopeDirName(tel,dirName);
  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: require that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  
  sprintf(fname,"%s/observatory/%s/rfi/persistentRFI.dat",runtimeDir,dirName);
  if (!(fin = fopen(fname,"r")))
    {
      printf("ERROR: unable to open filename >%s<\n",fname);
    }
  else
    {
      while (!feof(fin))
	{
	  if (fgets(str,4096,fin)!=NULL)
	    {
	      if (str[0]!='#' && strlen(str) > 2) // Not a comment line
		{
		  sscanf(str,"%d %s %s %lf %lf %lf %lf",&(rfi[*nRFI].type),
			 rfi[*nRFI].observatory,rfi[*nRFI].receiver,&(rfi[*nRFI].f0),&(rfi[*nRFI].f1),
			 &(rfi[*nRFI].mjd0),&(rfi[*nRFI].mjd1));
		  (*nRFI)++;
		  // SHOULD CHCEK IF EXCEEDS MAXRFI -- FIX ME
		}
	    }
	}
    }
  fclose(fin);
}

void sdhdf_loadTransientRFI(sdhdf_transient_rfi *rfi,int *nRFI,int maxRFI,char *tel)
{
  FILE *fin;
  char fname[MAX_STRLEN];
  char runtimeDir[MAX_STRLEN];
  char str[4096];
  char dirName[1024];
  
  *nRFI = 0;


  sdhdf_getTelescopeDirName(tel,dirName);
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: require that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  
  sprintf(fname,"%s/observatory/%s/rfi/transientRFI.dat",runtimeDir,dirName);
  if (!(fin = fopen(fname,"r")))
    {
      printf("ERROR: unable to open filename >%s<\n",fname);
    }
  else
    {
      while (!feof(fin))
	{
	  if (fgets(str,4096,fin)!=NULL)
	    {
	      if (str[0]!='#' && strlen(str) > 2) // Not a comment line
		{
		  sscanf(str,"%d %s %s %lf %lf %lf %lf %lf %lf %lf",&(rfi[*nRFI].type),
			 rfi[*nRFI].observatory,rfi[*nRFI].receiver,&(rfi[*nRFI].f0),&(rfi[*nRFI].f1),
			 &(rfi[*nRFI].mjd0),&(rfi[*nRFI].mjd1),&rfi[*nRFI].threshold,&rfi[*nRFI].f2,&rfi[*nRFI].f3);
		  (*nRFI)++;
		  // SHOULD CHCEK IF EXCEEDS MAXRFI -- FIX ME
		}
	    }
	}
    }
  fclose(fin);
}

// Obtain the directory name corresponding to a telescope name
int sdhdf_getTelescopeDirName(char *tel,char *dir)
{
  char fname1[1024],fname2[1024];
  char loadLine[1024], observatoryDir[1024]="NULL";
  FILE *fin;
  char runtimeDir[1024];
  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("[sdhdf_getTelescopeDirName] ERROR: SDHDF_RUNTIME is not set\n");      
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  sprintf(fname1,"%s/observatory/observatories.list",runtimeDir);
  fin = fopen(fname1,"r");
  while (!feof(fin))
    {
      if (fgets(loadLine,1024,fin)!=NULL)
	{
	  if (loadLine[0]!='#')
	    {
	      if (strstr(loadLine,tel)!=NULL)
		{
		  sscanf(loadLine,"%s",observatoryDir);
		  break;
		}
	    }
	}	
    }
  fclose(fin);
  if (strcmp(observatoryDir,"NULL")==0)
    {
      printf("ERROR: in sdhdfProc_metadata.c - unable to find observatory %s in %s\n",tel,fname1);
      exit(1);
    }
  strcpy(dir,observatoryDir);

  return 0;
  
}
