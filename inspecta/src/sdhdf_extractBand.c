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
// Software to produce a new SDHDF file based on parts of an existing one
//
// Usage:
// sdhdf_extract -f <filename.hdf> -o <outputFile.hdf>
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"

#define MAX_BANDS 26     // FIX THIS

void help()
{
  printf("sdhdf_extractBand\n");
  printf("\n");
  printf("Command line arguments\n\n");
  printf("-e <ext>          file extension for output files\n");
  printf("-h                this help\n");
  printf("-zoom <f1> <f2>   produce zoom band between f1 and f2 MHz\n");
  printf("-b <bandLabel>    select sub-band with label (if the label is not found then it assumes a band number)\n");

  exit(1);
}

int main(int argc,char *argv[])
{
  char args[MAX_ARGLEN]="";
  int ii,i,j,k,kk,l,nchan,totNchan,b,nd;
  char fname[MAX_FILES][64];
  int nFiles=0;
  char ext[MAX_STRLEN]="extract";
  char oname[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*outFile;
  herr_t status;
  char selectBand[MAX_BANDS][MAX_STRLEN];
  float *outVals,*freqVals,*inData;
  sdhdf_obsParamsStruct  *outObsParams;
  int  nSelectBands=0;
  int  copyBand=0;
  float fSelect0[MAX_BANDS];
  float fSelect1[MAX_BANDS];
  int zoomBand=0;
  sdhdf_bandHeaderStruct *outBandParams,*outCalBandParams,*inBandParams,*inCalBandParams;
  int nBand=0;
  char zoomLabel[MAX_STRLEN];
  char groupName[MAX_STRLEN];
  char groupName1[MAX_STRLEN],groupName2[MAX_STRLEN];
  int cal=0;
  int selectBandID=-1;
  int npol;
  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;
  
  strcpy(oname,"sdhdf_extract_output.hdf");

  if (argc==1)
    help();
  sdhdf_storeArguments(args,MAX_ARGLEN,argc,argv);
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-h")==0)
	help();
    }
  printf("Allocating memory\n");
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }
  
  printf("Checking input\n");
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-e")==0)
	strcpy(ext,argv[++i]);
      else if (strcmp(argv[i],"-zoom")==0)
	{
	  sscanf(argv[++i],"%f",&fSelect0[zoomBand]);
	  sscanf(argv[++i],"%f",&fSelect1[zoomBand]);
	  zoomBand++;
	}
      else if (strcmp(argv[i],"-b")==0 || strcmp(argv[i],"-sb")==0 || strcmp(argv[i],"-band")==0)
	strcpy(selectBand[nSelectBands++],argv[++i]);
      else
	{
	  strcpy(fname[nFiles],argv[i]);
	  nFiles++;
	}
    }
  printf("Number of files = %d, number of zoom bands = %d\n",nFiles,zoomBand);
  for (ii=0;ii<nFiles;ii++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_initialiseFile(outFile);
      sdhdf_formOutputFilename(fname[ii],ext,oname);
      
      if (zoomBand>0)
	nSelectBands=zoomBand;
      
      sdhdf_openFile(fname[ii],inFile,1);
      sdhdf_openFile(oname,outFile,3);
      
      if (inFile->fileID!=-1) // Did we successfully open the file?
	{
	  sdhdf_loadMetaData(inFile);
	  if (strcmp(inFile->primary[0].cal_mode,"ON")==0)
	    cal=1;
	  else
	    cal=0;
	  
	  printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname[ii],inFile->primary[0].utc0,
		 inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[0].source,inFile->primary[0].telescope,
		 inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[0].nBand);
       
	  for (b=0;b<inFile->nBeam;b++)
	    {
	      printf("Setting memory %d %d\n",inFile->beam[b].nBand,nSelectBands);
	      inBandParams  = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);
	      outBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nSelectBands);

	      sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);

	      if (cal==1)
		{
		  inCalBandParams  = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);         
		  outCalBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nSelectBands);
		  sdhdf_copyBandHeaderStruct(inFile->beam[b].calBandHeader,inCalBandParams,inFile->beam[b].nBand);
		}
	      nBand=0;

	      if (zoomBand==0)
		{
		  copyBand=0;
		  printf("Number of bands = %d\n",inFile->beam[b].nBand);
		  for (j=0;j<nSelectBands;j++)
		    {
		      copyBand=0;
		      for (i=0;i<inFile->beam[b].nBand;i++)
			{
			  
			  if (strcmp(inFile->beam[b].bandHeader[i].label,selectBand[j])==0)
			    {
			      copyBand=1;
			      break;
			    }
			}
		      if (copyBand==0)
			{
			  int isel;
			  printf("WARNING: Unable to find band label >%s<\n",selectBand[j]);
			  if (sscanf(selectBand[j],"%d",&isel)==1)
			    {
			      printf("Using the band label as the band number instead\n");
			      strcpy(selectBand[j],inFile->beam[b].bandHeader[isel].label);
			    }
			  else
			    {
			      printf("ERROR: Unable to identify this band\n");
			      exit(1); // Should clean nicely! ** FIX ME
			    }
			}
		    }

		      

		  // Check if at least one band is to be copied
		  for (i=0;i<inFile->beam[b].nBand;i++)
		    {
		      copyBand=0;
		      for (j=0;j<nSelectBands;j++)
			{
			  if (strcmp(inFile->beam[b].bandHeader[i].label,selectBand[j])==0)
			    {
			      copyBand=1;
			      break;
			    }
			}
		      printf("Copying band: %d %d\n",i,copyBand);
		      if (copyBand==1)
			{
			  //			  sprintf(groupName,"beam_%d",b);
			  strcpy(groupName,inFile->beamHeader[b].label);
			  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
			    {
			      hid_t group_id;
			      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
			      status = H5Gclose(group_id);
			    }
			  sprintf(groupName,"%s/%s",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
			  sdhdf_copyEntireGroup(groupName,inFile,outFile);

			  strcpy(outBandParams[nBand].label,inBandParams[i].label);
			  outBandParams[nBand].fc = inBandParams[i].fc;
			  outBandParams[nBand].f0 = inBandParams[i].f0;
			  outBandParams[nBand].f1 = inBandParams[i].f1;
			  outBandParams[nBand].nchan = inBandParams[i].nchan;
			  outBandParams[nBand].dtime = inBandParams[i].dtime;
			  outBandParams[nBand].npol = inBandParams[i].npol;
			  strcpy(outBandParams[nBand].pol_type, inBandParams[i].pol_type);
			  outBandParams[nBand].ndump = inBandParams[i].ndump;
			  npol = inBandParams[i].npol;
			  
			  if (cal==1)
			    {
			      strcpy(outCalBandParams[nBand].label,inCalBandParams[i].label);
			      outCalBandParams[nBand].fc = inCalBandParams[i].fc;
			      outCalBandParams[nBand].f0 = inCalBandParams[i].f0;
			      outCalBandParams[nBand].f1 = inCalBandParams[i].f1;
			      outCalBandParams[nBand].nchan = inCalBandParams[i].nchan;
			      outCalBandParams[nBand].dtime = inCalBandParams[i].dtime;
			      outCalBandParams[nBand].npol = inCalBandParams[i].npol;
			      strcpy(outCalBandParams[nBand].pol_type, inCalBandParams[i].pol_type);
			      outCalBandParams[nBand].ndump = inCalBandParams[i].ndump;
			    }
			  nBand++;
			}
		    }
		}
	      else
		{
		  for (j=0;j<zoomBand;j++)
		    {

		      sprintf(zoomLabel,"band_zoom%03d",j);
		      printf("Processing zoom band: %d/%d (%s)\n",j,zoomBand-1,zoomLabel);
		      strcpy(outBandParams[j].label,zoomLabel);
		      outBandParams[j].fc = (fSelect0[j] + fSelect1[j])/2.;
		      outBandParams[j].f0 = fSelect0[j];
		      outBandParams[j].f1 = fSelect1[j];

		      // Find which band this corresponds to
		      //
		      // Count nchan (note that we may go across band edges)
		      totNchan=0;
		      selectBandID=-1;
		      for (k=0;k<inFile->beam[b].nBand;k++)
			{
			  sdhdf_loadBandData(inFile,b,k,1);
			  for (l=0;l<inFile->beam[b].bandHeader[k].nchan;l++)
			    {
			      if (inFile->beam[b].bandData[k].astro_data.freq[l] >= fSelect0[j] && // FIX ME :: FREQ FOR CORRECT DUMP
				  inFile->beam[b].bandData[k].astro_data.freq[l] <= fSelect1[j])
				{
				  if (selectBandID == -1) {selectBandID=k;}
				  totNchan++;
				}
			    }
			}
		      printf("Selecting band %d\n",selectBandID);
		      strcpy(outBandParams[j].pol_type,inBandParams[selectBandID].pol_type);
		      npol = inBandParams[selectBandID].npol;
		      outBandParams[j].npol = inBandParams[selectBandID].npol;
		      outBandParams[j].ndump = inBandParams[selectBandID].ndump; 
		      outObsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*outBandParams[j].ndump);

		      for (kk=0;kk<outBandParams[j].ndump;kk++)
			sdhdf_copySingleObsParams(inFile,b,selectBandID,kk,&outObsParams[kk]);
		    		 		      
		      if (cal==1)
			{
			  strcpy(outCalBandParams[j].label,zoomLabel);
			  // Copy the entire cal
			  outCalBandParams[j].fc = (fSelect0[j] + fSelect1[j])/2.;
			  outCalBandParams[j].f0 = fSelect0[j];
			  outCalBandParams[j].f1 = fSelect1[j];
			  strcpy(outCalBandParams[j].pol_type,inCalBandParams[selectBandID].pol_type);
			  outCalBandParams[j].nchan  = inCalBandParams[selectBandID].nchan;
			  outCalBandParams[j].npol   = inCalBandParams[selectBandID].npol;
			  outCalBandParams[j].ndump  = inCalBandParams[selectBandID].ndump;
			  outCalBandParams[j].dtime  = inCalBandParams[selectBandID].dtime;
			}
		    	      
		      outBandParams[j].nchan = totNchan;
		      outBandParams[j].dtime = inBandParams[selectBandID].dtime;
		      outVals  = (float *)malloc(sizeof(float)*outBandParams[j].nchan*4*outBandParams[j].ndump); // Fix 4 = npol, 1 = ndump
		      freqVals = (float *)malloc(sizeof(float)*outBandParams[j].nchan); 
		      printf("Total nchan = %d\n",totNchan);
		      
		      printf("Using original channel %d for number of spectral dumps, calibration solution etc.\n",selectBandID);
		      
		      for (nd=0;nd<inBandParams[selectBandID].ndump;nd++)
			{
			  nchan=0;
			  for (k=0;k<inFile->beam[b].nBand;k++)
			    {
			      for (l=0;l<inFile->beam[b].bandHeader[k].nchan;l++)
				{
				  // FIX ME: FREQ FOR CORRECT DUMP
				  if (inFile->beam[b].bandData[k].astro_data.freq[l] >= fSelect0[j] && inFile->beam[b].bandData[k].astro_data.freq[l] <= fSelect1[j])
				    {   
				      if (npol==4)
					{
					  outVals[nd*totNchan*4+nchan]            = inFile->beam[b].bandData[k].astro_data.pol1[l+nd*inFile->beam[b].bandHeader[k].nchan]; 
					  outVals[nd*totNchan*4+nchan+totNchan]   = inFile->beam[b].bandData[k].astro_data.pol2[l+nd*inFile->beam[b].bandHeader[k].nchan]; 
					  outVals[nd*totNchan*4+nchan+2*totNchan] = inFile->beam[b].bandData[k].astro_data.pol3[l+nd*inFile->beam[b].bandHeader[k].nchan];
					  outVals[nd*totNchan*4+nchan+3*totNchan] = inFile->beam[b].bandData[k].astro_data.pol4[l+nd*inFile->beam[b].bandHeader[k].nchan];
					}
				      else if (npol==2)
					{
					  outVals[nd*totNchan*2+nchan]            = inFile->beam[b].bandData[k].astro_data.pol1[l+nd*inFile->beam[b].bandHeader[k].nchan]; 
					  outVals[nd*totNchan*2+nchan+totNchan]   = inFile->beam[b].bandData[k].astro_data.pol2[l+nd*inFile->beam[b].bandHeader[k].nchan]; 
					}
				      else if (npol==1)
					outVals[nd*totNchan+nchan]            = inFile->beam[b].bandData[k].astro_data.pol1[l+nd*inFile->beam[b].bandHeader[k].nchan]; 
					
					  //				      printf("Seting freq\n");
				      // FIX ME: FREQ FOR CORRECT DUMP
				      if (nd==0) freqVals[nchan] = inFile->beam[b].bandData[k].astro_data.freq[l];
				      nchan++;
				    }
				}
			    }
			  //			  free(inData);
			}

		      sdhdf_copyAttributes(inFile->beam[b].bandData[selectBandID].astro_obsHeaderAttr,inFile->beam[b].bandData[selectBandID].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
		      sdhdf_copyAttributes(inFile->beam[b].bandData[selectBandID].astro_obsHeaderAttr_freq,inFile->beam[b].bandData[selectBandID].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);
		    
		      printf("Output nchan = %d (%d), npol = %d\n",nchan,totNchan,npol);
		      //		      sdhdf_writeSpectrumData(outFile,inFile,b,j,outVals,freqVals,nchan,4,1,0); // FIX 4,1,0

		      // FIX ME: Only sending 1 frequency channel through
		      sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,outBandParams[j].label,b,j,outVals,freqVals,1,totNchan,1,npol,outBandParams[j].ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
		      sdhdf_writeObsParams(outFile,outBandParams[j].label,inFile->beamHeader[b].label,j,outObsParams,outBandParams[j].ndump,1,inFile->primary[0].hdr_defn_version);
		    
		      sprintf(groupName1,"beam_%d/%s/metadata/cal_obs_params",b,inFile->beam[b].bandHeader[selectBandID].label);
		      sprintf(groupName2,"beam_%d/%s/metadata/cal_obs_params",b,outBandParams[j].label);
		      sdhdf_copyEntireGroupDifferentLabels(groupName1,inFile,groupName2,outFile);
		      //    sdhdf_writeObsParams(outFile,outBandParams[j].label,b,j,outCalObsParams,outCalBandParams[j].ndump,2); 

		      sprintf(groupName1,"%s/%s/calibrator_data",inFile->beamHeader[b].label,inFile->beam[b].bandHeader[selectBandID].label);
		      sprintf(groupName2,"%s/%s/calibrator_data",inFile->beamHeader[b].label,outBandParams[j].label);
		      printf("Copying %s to %s\n",groupName1,groupName2);
		      sdhdf_copyEntireGroupDifferentLabels(groupName1,inFile,groupName2,outFile);
		      
		      //		      sdhdf_writeNewBand(outFile,j,outVals,freqVals,&(outBandParams[j]));
		      printf("Completed writing out the spectrum for zoomband %d\n",j);
		   
		      free(outVals);
		      free(freqVals);
		      free(outObsParams);
		    }

		}
	      sdhdf_writeBandHeader(outFile,outBandParams,inFile->beamHeader[b].label,nSelectBands,1,inFile->primary[0].hdr_defn_version);

	      if (cal==1)
		sdhdf_writeBandHeader(outFile,outCalBandParams,inFile->beamHeader[b].label,nSelectBands,2,inFile->primary[0].hdr_defn_version);

	      
	      //	      sprintf(groupName,"beam_%d/metadata",b);
	      //	      sdhdf_copyEntireGroup(groupName,inFile,outFile);	      	      

	      inFile->beamHeader[b].nBand = nSelectBands;

	    }
       
	  // Copy other primary tables
	  sdhdf_copyEntireGroup("metadata",inFile,outFile);	      	      
	  sdhdf_writeBeamHeader(outFile,inFile->beamHeader,inFile->nBeam,inFile->primary[0].hdr_defn_version);
	  // Don't want to "copyRemainder" as have only selected specific bands
	  sdhdf_copyEntireGroup("config",inFile,outFile);
	  sdhdf_addHistory(inFile->history,inFile->nHistory,"sdhdf_extractBand","INSPECTA software to extractBands",args);
	  inFile->nHistory++;
	  printf("Writing history %d\n",inFile->nHistory);
	  sdhdf_writeHistory(outFile,inFile->history,inFile->nHistory);
	  printf("Done writing history\n");
	  free(outBandParams);
	  free(inBandParams);
	  if (cal==1)
	    {
	      free(outCalBandParams);
	      free(inCalBandParams);
	    }
	  sdhdf_closeFile(outFile);
	  sdhdf_closeFile(inFile);
	}
    }
  free(inFile);
  free(outFile);
}


