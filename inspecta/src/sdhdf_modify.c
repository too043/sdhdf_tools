
// Should have a normalise option
// Need to deal with wraps when averaging RAs/DECs/times etc.
// Should be able to do some beam or bin manipulation
// Should be able to manipulate the noise source
// History table
// lsr_regrid, bary_regrid, -regrid

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
// Software to produce a new SDHDF file based processing an existing one
//
// Usage:
// sdhdf_modify <filename.hdf> -o <outputFile.hdf>
//
//

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"
#include "TKfit.h"
#include <libgen.h>
#include "T2toolkit.h"

#define VNUM "v2.0"
#define MAX_COMMANDS 128

// Structure of commands
// type = 1 = time averaging, param1 = number of dumps to average, param2 = 1 = sum, = 2 = average
// type = 2 = frequency averaging, param2 = number of channels to average, param2 = 1 = sum, = 2 = average
// type = 3 = polarisation sum (AA + BB)  if param1 = 1, or average (AA+BB)/2 is param1 = 2
// type = 4 = convert to bary (param1 = 1) or lsr (param1 = 2) with (param2 = 2) or without (param2 = 1) regridding
typedef struct commandStruct {
  int type;
  float param1;
  float param2;
  int   iparam;
  int   iparam2;
} commandStruct;

// Structure allowing the data to be passed through multiple processing steps
typedef struct dataStruct {
  // Data ordering = dump,pol,channel
  float *data;
  float *data2; // For noise source
  double *freq; 
  int nFreqDump;
  float *wt;
  unsigned char *flag;
  int nchan;
  int npol;
  int ndump;
  int memoryAllocated;

  sdhdf_obsParamsStruct  *obsParams;

  int astroCal;
} dataStruct;


void processFile(char *fname,char *oname, commandStruct *commands, int nCommands,char *args,int astroCal,int singleFreqAxis,int verbose,int showVel,double mjdSet);
void processCommand(dataStruct *in,dataStruct *out,int iband, commandStruct *command,sdhdf_attributes_struct *freqAttributes,int nFreqAttributes,int verbose,int showVel,double mjdSet);
void timeAverage(dataStruct *in,dataStruct *out,int ndumpAv,int sum);
void frequencyAverage(dataStruct *in,dataStruct *out,int nfreqAv,int sum,int meanMedian);
void polarisationAverage(dataStruct *in,dataStruct *out,int sum);
void changeFrequencyAxis(dataStruct *in,dataStruct *out,int bary_lsr,int regrid,sdhdf_attributes_struct *freqAttributes,int nFreqAttributes,int verbose,int showVel,double mjdSet);
void scaleValues(dataStruct *in,dataStruct *out,int iband,  int useBand, int multDiv, float aaScale, float bbScale);
void allocateMemory(dataStruct *in);

//void timeAverage(float *in_data,int nchan,int npol,int ndump,int ntime,float *out_data,int *out_nchan,int *out_npol,int *out_ndump);

void help()
{
  printf("\nsdhdf_modify      %s\n",VNUM);
  printf("INSPECTA version: %s\n",SOFTWARE_VER);
  printf("Author:           George Hobbs\n");
  printf("Software to modify an SDHDF file\n");

  printf("\nCommand line arguments:\n\n");
  printf("-h                This help\n");
  printf("-bary_noregrid    Convert frequency axis to the barycentre\n");
  printf("-bary_regrid      Convert frequency axis to the barycentre and regrid back to the original frequency axis\n");
  printf("-cal              Process the noise source, not the astronomy data\n");
  printf("-e <ext>          Use output file extension <ext>\n");
  printf("-F or -Fsum       Completely sum frequency channels\n");
  printf("-Fav              Completely average frequency channels (weighted mean)\n");
  printf("-fav <val>        Average (weighted mean) val frequency channels together\n");
  printf("-favMedian <val>  Average (*un*weighted median) val frequency channels together\n");
  printf("-lsr_noregrid     Convert frequency axis to the LSR\n");
  printf("-lsr_regrid       Convert frequency axis to the LSR and regrid back to the original frequency axis\n");
  printf("-p1 or -p1_sum    Summation of P0 and P1\n");
  printf("-p1_av            (P0+P1)/2 \n");
  printf("-singleFreqAxis   Only write out a single frequency dimension (nchan) instead of (ndump x nchan)\n");
  printf("-T or -Tsum       Completely sum in time\n");
  printf("-Tav              Completely average (weighted mean) in time\n");
  printf("-tav <val>        Average (weighted mean) val time dumps together\n");
  printf("-tsum <val>       Sum val time dumps together\n");
  printf("-v                Verbose output\n");

  printf("\nExample:\n\n");
  printf("sdhdf_modify -T -e T *.hdf\n\n");

  exit(1);
}

int main(int argc,char *argv[])
{
  int i;
  int nCommands=0;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  char oname[MAX_STRLEN];
  char extension[MAX_STRLEN];
  commandStruct *commands;
  char args[MAX_ARGLEN]="";
  int verbose=0;
  int astroCal=0; // Modify the astronomy data
  int singleFreqAxis=0;
  int showVel=0;
  double mjdSet=-1;
  
  
  strcpy(extension,"modify");
  
  commands = (commandStruct *)malloc(sizeof(commandStruct)*MAX_COMMANDS);

  sdhdf_storeArguments(args,MAX_ARGLEN,argc,argv);
  if (argc==1)
    help();
	
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-Tsum")==0 || strcmp(argv[i],"-T")==0)        // Completely sum in time
	{commands[nCommands].type=1; commands[nCommands].param2 = 1; commands[nCommands++].param1 = 0;}
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcasecmp(argv[i],"-singleFreqAxis")==0)
	  singleFreqAxis=1;
      else if (strcmp(argv[i],"-v")==0)
	verbose=1;    
      else if (strcmp(argv[i],"-cal")==0)
	astroCal=1;
      else if (strcmp(argv[i],"-divide2pol")==0) 
	{commands[nCommands].type=5; commands[nCommands].iparam = 1; sscanf(argv[++i],"%d",&(commands[nCommands].iparam2)); sscanf(argv[++i],"%f",(&commands[nCommands].param1)); sscanf(argv[++i],"%f",(&commands[nCommands++].param2)); }
      else if (strcmp(argv[i],"-mult2pol")==0) 
	{commands[nCommands].type=5; commands[nCommands].iparam = 2; sscanf(argv[++i],"%d",&(commands[nCommands].iparam2)); sscanf(argv[++i],"%f",(&commands[nCommands].param1)); sscanf(argv[++i],"%f",(&commands[nCommands++].param2)); }
      else if (strcmp(argv[i],"-Tav")==0)    // Completely average in time
	{commands[nCommands].type=1; commands[nCommands].param2 = 2; commands[nCommands++].param1 = 0;}
      else if (strcmp(argv[i],"-tsum")==0)   // Sum <val> time dumps together
	{commands[nCommands].type=1; commands[nCommands].param2 = 1; sscanf(argv[++i],"%f",(&commands[nCommands++].param1)); }
      else if (strcmp(argv[i],"-Fsum")==0 || strcmp(argv[i],"-F")==0)   // Sum all> frequency channels together
	{commands[nCommands].type=2; commands[nCommands].iparam = 1; commands[nCommands].param2 = 1; commands[nCommands++].param1 = 0;}
      else if (strcmp(argv[i],"-Fav")==0)
	{commands[nCommands].type=2; commands[nCommands].iparam = 1; commands[nCommands].param2 = 2; commands[nCommands++].param1 = 0;}
      else if (strcmp(argv[i],"-fsum")==0)   // Sum <val> frequency channels together
	{commands[nCommands].type=2; commands[nCommands].iparam = 1; commands[nCommands].param2 = 1; sscanf(argv[++i],"%f",(&commands[nCommands++].param1));}
      else if (strcmp(argv[i],"-tav")==0)    // Average <val> time dumps together
	{commands[nCommands].type=1; commands[nCommands].param2 = 2; sscanf(argv[++i],"%f",(&commands[nCommands++].param1)); }
      else if (strcmp(argv[i],"-fav")==0)   // Average <val> frequency channels together
	{commands[nCommands].type=2; commands[nCommands].iparam = 1; commands[nCommands].param2 = 2; sscanf(argv[++i],"%f",(&commands[nCommands++].param1));}
      else if (strcmp(argv[i],"-favMedian")==0)   // Average <val> frequency channels together
	{commands[nCommands].type=2; commands[nCommands].iparam = 2; commands[nCommands].param2 = 2; sscanf(argv[++i],"%f",(&commands[nCommands++].param1));}
      else if (strcmp(argv[i],"-p1")==0 || strcmp(argv[i],"-p1_sum")==0)    // Form 1 polarisation (AA+BB)
	{commands[nCommands].type=3; commands[nCommands++].param1 = 1;}
      else if (strcmp(argv[i],"-p1_av")==0)    // Form 1 polarisation (AA+BB)/2
	{commands[nCommands].type=3; commands[nCommands++].param1 = 2;}
      else if (strcmp(argv[i],"-lsr_noregrid")==0)
	{commands[nCommands].type=4; commands[nCommands].param1 = 2; commands[nCommands++].param2 = 1;}
      else if (strcmp(argv[i],"-lsr2_noregrid")==0)
	{commands[nCommands].type=4; commands[nCommands].param1 = 3; commands[nCommands++].param2 = 1;}
      else if (strcmp(argv[i],"-bary_noregrid")==0)
	{commands[nCommands].type=4; commands[nCommands].param1 = 1; commands[nCommands++].param2 = 1;}
      else if (strcmp(argv[i],"-lsr_regrid")==0)
	{commands[nCommands].type=4; commands[nCommands].param1 = 2; commands[nCommands++].param2 = 2;}
      else if (strcmp(argv[i],"-lsr2_regrid")==0)
	{commands[nCommands].type=4; commands[nCommands].param1 = 3; commands[nCommands++].param2 = 2;}
      else if (strcmp(argv[i],"-bary_regrid")==0)
	{commands[nCommands].type=4; commands[nCommands].param1 = 1; commands[nCommands++].param2 = 2;}
      else if (strcasecmp(argv[i],"-showVel")==0)
	showVel=1;
      else if (strcasecmp(argv[i],"-mjdHack")==0)
	sscanf(argv[++i],"%lf",&mjdSet);
      else if (strcmp(argv[i],"-e")==0)
	{strcpy(extension,argv[++i]);}
      else
	strcpy(fname[nFiles++],argv[i]);
    }


  // List the commands
  printf("Processing: %d files with %d commands\n",nFiles,nCommands);

  for (i=0;i<nFiles;i++)
    {
      printf("Processing file: %s\n",fname[i]);
      sdhdf_formOutputFilename(fname[i],extension,oname);
      processFile(fname[i],oname,commands,nCommands,args,astroCal, singleFreqAxis,verbose,showVel,mjdSet);
    }


  free(commands);
}


void processFile(char *fname,char *oname, commandStruct *commands, int nCommands,char *args,int astroCal,int singleFreqAxis,int verbose,int showVel,double mjdSet)
{
  int ii,i,c,j,k;
  
  sdhdf_fileStruct *inFile,*outFile;
  sdhdf_bandHeaderStruct *inBandParams;
  sdhdf_bandHeaderStruct *outBandParams;


  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;

  
  int nBeam,b,ndump;
  dataStruct *dset;
  dataStruct *in,*out,*swp;  
  
  dset = (dataStruct *)malloc(sizeof(dataStruct)*2);
  dset[0].memoryAllocated = 0;
  dset[1].memoryAllocated = 0;
  in  = &dset[0];
  out = &dset[1]; 

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
     

  sdhdf_initialiseFile(inFile);
  sdhdf_initialiseFile(outFile);

  sdhdf_openFile(fname,inFile,1);
  sdhdf_openFile(oname,outFile,3);
  sdhdf_loadMetaData(inFile);

  nBeam = inFile->nBeam;
  sdhdf_allocateBeamMemory(outFile,nBeam);
  for (b=0; b<nBeam; b++)
    {
      printf("Processing beam %d/%d\n",b,nBeam-1);
      inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);
      outBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);      
      sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);
      sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,outBandParams,inFile->beam[b].nBand);

      // Each band is processed independently
      for (ii=0;ii<inFile->beam[b].nBand;ii++)
	{
	  printf(" ... Processing subband %d/%d\n",ii,inFile->beam[b].nBand-1);
	  if (astroCal==0)
	    {
	      in->nchan = inFile->beam[b].bandHeader[ii].nchan;
	      in->npol  = inFile->beam[b].bandHeader[ii].npol;
	      in->ndump = inFile->beam[b].bandHeader[ii].ndump;
	      sdhdf_loadBandData(inFile,b,ii,1);
	      in->nFreqDump = inFile->beam[b].bandData[ii].astro_data.nFreqDumps;
	      in->astroCal = 0;
	    }
	  else
	    {
	      in->nchan = inFile->beam[b].calBandHeader[ii].nchan;
	      in->npol  = inFile->beam[b].calBandHeader[ii].npol;
	      in->ndump = inFile->beam[b].calBandHeader[ii].ndump;
	      sdhdf_loadBandData(inFile,b,ii,2);
	      sdhdf_loadBandData(inFile,b,ii,3);
	      in->nFreqDump = in->nchan; // NFREQ_DUMPS ?
	      in->astroCal = 1;
	    }
	  // Load in the data
	  if (verbose==1) printf("Allocating memory\n");
	  allocateMemory(in);
	  if (verbose==1) printf("Complete\n");

	  sdhdf_copyAttributes(inFile->beam[b].bandData[ii].astro_obsHeaderAttr,inFile->beam[b].bandData[ii].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
	  sdhdf_copyAttributes(inFile->beam[b].bandData[ii].astro_obsHeaderAttr_freq,inFile->beam[b].bandData[ii].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);

	  if (astroCal==0)
	    {
	      if (verbose==1) printf("Copying memory\n");
	      for (i=0;i<in->nFreqDump;i++)
		{
		  memcpy(in->freq+i*in->nchan,inFile->beam[b].bandData[ii].astro_data.freq+i*in->nchan,sizeof(double)*in->nchan);
		}
	      if (verbose==1) printf("Complete memory\n");
	    }
	  else
	    memcpy(in->freq,inFile->beam[b].bandData[ii].cal_on_data.freq,sizeof(double)*in->nchan);
	  if (verbose==1) printf("Copying observation parameters\n");
	  for (i=0;i<in->ndump;i++)
	    {
	      if (astroCal==0)
		sdhdf_copySingleObsParams(inFile,b,ii,i,&in->obsParams[i]); 
	      else
		sdhdf_copySingleObsParamsCal(inFile,b,ii,i,&in->obsParams[i]); 
	    }
	  if (verbose==1) printf("Complete copying observation parameters\n");
	  //	  printf("POS A\n");
	  //	  printf("astroCal = %d\n",astroCal);
	  //	  printf("nchan =  %d\n",in->nchan);
	  //	  printf("ndump = %d\n",in->ndump);
	  if (in->npol==1)
	    {printf("ERROR: SORRY, for now only 2 pol data has been implemented not 1 pol\n"); exit(1);}
	    
	    
	  for (k=0;k<in->nchan;k++)
	    {
	      for (j=0;j<in->ndump;j++)
		{
		  // SHOULD CHECK NUMBER OF POLS -- FIX ME
		  // SHOULD HAVE A WAY TO FGET THE DATA DIRECTLY IN THE CORRECT FORMAT AS IT IS STORED IN THAT FORMAT ANYWAY
		  if (astroCal==0)
		    {
		      in->data[(long)((long)j*in->nchan*in->npol) + k] = inFile->beam[b].bandData[ii].astro_data.pol1[(long)k+(long)j*in->nchan];
		      in->data[(long)((long)j*in->nchan*in->npol) + (long)k + (long)in->nchan] =
			inFile->beam[b].bandData[ii].astro_data.pol2[(long)k+(long)j*in->nchan];
		      in->data[(long)((long)j*in->nchan*in->npol) + k + 2*in->nchan] = inFile->beam[b].bandData[ii].astro_data.pol3[k+(long)j*in->nchan];
		      in->data[(long)((long)j*in->nchan*in->npol) + k + 3*in->nchan] = inFile->beam[b].bandData[ii].astro_data.pol4[k+(long)j*in->nchan];
		      in->wt[(long)j*in->nchan + k] = inFile->beam[b].bandData[ii].astro_data.dataWeights[k+(long)j*in->nchan];
		      in->flag[(long)j*in->nchan + k] = inFile->beam[b].bandData[ii].astro_data.flag[k+(long)j*in->nchan];
		    }
		  else
		    {
		      in->data[j*in->nchan*in->npol + k] = inFile->beam[b].bandData[ii].cal_on_data.pol1[k+j*in->nchan];
		      in->data[j*in->nchan*in->npol + k + in->nchan] = inFile->beam[b].bandData[ii].cal_on_data.pol2[k+j*in->nchan];
		      in->data[j*in->nchan*in->npol + k + 2*in->nchan] = inFile->beam[b].bandData[ii].cal_on_data.pol3[k+j*in->nchan];
		      in->data[j*in->nchan*in->npol + k + 3*in->nchan] = inFile->beam[b].bandData[ii].cal_on_data.pol4[k+j*in->nchan];
		    
		      in->data2[j*in->nchan*in->npol + k] = inFile->beam[b].bandData[ii].cal_off_data.pol1[k+j*in->nchan];
		      in->data2[j*in->nchan*in->npol + k + in->nchan] = inFile->beam[b].bandData[ii].cal_off_data.pol2[k+j*in->nchan];
		      in->data2[j*in->nchan*in->npol + k + 2*in->nchan] = inFile->beam[b].bandData[ii].cal_off_data.pol3[k+j*in->nchan];
		      in->data2[j*in->nchan*in->npol + k + 3*in->nchan] = inFile->beam[b].bandData[ii].cal_off_data.pol4[k+j*in->nchan];
		      in->wt[j*in->nchan + k]   = 1;    // FIX ME HERE, currently do not FLAG or WEIGHT the noise source  and only 1 frequency dump
		      in->flag[j*in->nchan + k] = 0; 
		    }
		}
	    }
	  for (c=0;c<nCommands;c++)
	    {
	      processCommand(in,out,ii,&commands[c],freqAttributes,nFreqAttributes,verbose,showVel,mjdSet);
	      if (c!=nCommands-1)
		{ swp = in; in = out; out = swp;} 
	    }      
	  // Sort out the observation parameters

	  if (astroCal==0)
	    {
	      if (singleFreqAxis==0)
		sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[ii].label,b,ii,out->data,
					out->freq,out->nFreqDump,out->nchan,1,out->npol,out->ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	      else
		sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[ii].label,b,ii,out->data,
					out->freq,1,out->nchan,1,out->npol,out->ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
		
	      sdhdf_writeObsParams(outFile,inFile->beam[b].bandHeader[ii].label,inFile->beamHeader[b].label,ii,out->obsParams,out->ndump,1,inFile->primary[0].hdr_defn_version);
	    }
	  else
	    {
	      sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[ii].label,b,ii,out->data,
				      out->freq,out->nFreqDump,out->nchan,1,out->npol,out->ndump,2,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	      
	      sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[ii].label,b,ii,out->data2,
				      out->freq,out->nFreqDump,out->nchan,1,out->npol,out->ndump,3,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	      sdhdf_writeObsParams(outFile,inFile->beam[b].bandHeader[ii].label,inFile->beamHeader[b].label,ii,out->obsParams,out->ndump,2,inFile->primary[0].hdr_defn_version);
	    }
	  // Write out the obs_params file

	  sdhdf_writeDataWeights(outFile,b,ii,out->wt,out->nchan,out->ndump,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[ii].label);
	  sdhdf_writeFlags(outFile,b,ii,out->flag,out->nchan,out->ndump,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[ii].label);
	  	  
	  // Setup the band parameters
	  outBandParams[ii].nchan = out->nchan;
	  outBandParams[ii].npol  = out->npol;
	  // FIX ME: NOT PICKING DUMP HERE in FREQ[0]
	  outBandParams[ii].f0    = out->freq[0];
 	  outBandParams[ii].f1    = out->freq[out->nchan-1];
	  outBandParams[ii].ndump = out->ndump;
	  outBandParams[ii].dtime = out->obsParams[0].dtime;

	  if (astroCal==0)
	    sdhdf_releaseBandData(inFile,b,ii,1); 		      
	  else
	    {
	      sdhdf_releaseBandData(inFile,b,ii,2);
	      sdhdf_releaseBandData(inFile,b,ii,3); 		      
	    }
	}
      if (astroCal==0)
	sdhdf_writeBandHeader(outFile,outBandParams,inFile->beamHeader[b].label,inFile->beam[b].nBand,1,inFile->primary[0].hdr_defn_version);
      else
	sdhdf_writeBandHeader(outFile,outBandParams,inFile->beamHeader[b].label,inFile->beam[b].nBand,2,inFile->primary[0].hdr_defn_version);
      free(inBandParams);
      free(outBandParams);
    }
  // Clean up
  free(in->data);
  free(out->data);
  sdhdf_addHistory(inFile->history,inFile->nHistory,"sdhdf_modify","INSPECTA software to modify a file",args);
  inFile->nHistory++;
  sdhdf_writeHistory(outFile,inFile->history,inFile->nHistory);
  sdhdf_copyRemainder(inFile,outFile,0); 
  
  

  sdhdf_closeFile(inFile);
  sdhdf_closeFile(outFile);

  free(inFile); free(outFile);
  free(dset);
}

void processCommand(dataStruct *in,dataStruct *out,int iband,commandStruct *command,sdhdf_attributes_struct *freqAttributes,int nFreqAttributes,int verbose,int showVel,double mjdSet)
{
  if (verbose==1)
    printf("Processing command: %d\n",command->type);
  if (command->type==1)
    timeAverage(in,out,command->param1,command->param2);
  else if (command->type==2)
    frequencyAverage(in,out,command->param1,command->param2,command->iparam);
  else if (command->type==3)
    polarisationAverage(in,out,command->param1);
  else if (command->type==4)
    changeFrequencyAxis(in,out,command->param1,command->param2,freqAttributes,nFreqAttributes,verbose,showVel,mjdSet);
  else if (command->type==5)
    scaleValues(in,out,iband,command->iparam2,command->iparam,command->param1,command->param2);
}


void timeAverage(dataStruct *in,dataStruct *out,int ndumpAv,int sum)
{
  int i,j,k,kk,p;
  int include;
  int avDump,ndump;
  double av;
  float wt,wtVal;
  unsigned char flag;
  int nsum;
  double s_flag_wt;
  int flagVal;
  double dtime;
  double raDeg,decDeg,mjd,timeElapsed;
  double raOffset,decOffset,gl,gb,az,ze,el,az_drive_rate,ze_drive_rate,hourAngle,paraAngle,windDir,windSpd;
  char tempStr[1024];
  double turn;
  double turn_utc,turn_local;
  double turn_utc_av,turn_local_av;
  
  if (ndumpAv == 0)
    {ndump = 1; avDump = in->ndump;}
  else
    {ndump = (int)((double)in->ndump/(double)ndumpAv); avDump = ndumpAv;}


  out->nchan      = in->nchan;
  out->npol       =  in->npol;
  out->ndump      = ndump;
  out->nFreqDump  = ndump;
  printf("Averaging: %d dumps together giving %d dumps out of %d original dumps\n",avDump,ndump,in->ndump);

  allocateMemory(out);
  //  printf("A %d %d %d %d\n",out->nchan,out->npol,out->ndump,out->nchan*out->npol*out->ndump);
  //  printf("npol here = %d\n",out->npol);

  memcpy(out->freq,in->freq,sizeof(double)*out->nchan*out->nFreqDump);
  for (p=0;p<out->npol;p++)
    {
      for (j=0;j<out->nchan;j++)
	{
	  for (k=0;k<out->ndump;k++)
	    {
	      av=0;
	      wt=0;
	      flag = 1;
	      nsum=0;
	      s_flag_wt = 0;
	      if (p==0 && j==0)
		{
		  dtime = raDeg = decDeg = mjd = timeElapsed = 0;
		  raOffset = decOffset = gl = gb = az = ze = el = az_drive_rate = 0;
		  ze_drive_rate = hourAngle = paraAngle = windDir = windSpd = 0;
		  turn_utc_av = turn_local_av = 0;
		}

	      for (kk=0;kk<avDump;kk++)
		{
		  wtVal     = in->wt[((k*avDump)+kk)*in->nchan + j];
		  flagVal   = 1-in->flag[((k*avDump)+kk)*in->nchan + j];
		  av        += wtVal*flagVal*in->data[((k*avDump)+kk)*in->nchan*in->npol + p*in->nchan +j];		      
		  s_flag_wt += wtVal*flagVal;
		  if (p==0 && j==0)
		    {
		      dtime     += in->obsParams[(k*avDump)+kk].dtime;
		      raDeg     += in->obsParams[(k*avDump)+kk].raDeg;
		      decDeg    += in->obsParams[(k*avDump)+kk].decDeg;
		      mjd       += (in->obsParams[(k*avDump)+kk].mjd + in->obsParams[(k*avDump)+kk].fractional_mjd);
		      timeElapsed += in->obsParams[(k*avDump)+kk].timeElapsed; 
		      raOffset += in->obsParams[(k*avDump)+kk].raOffset;
		      decOffset += in->obsParams[(k*avDump)+kk].decOffset;
		      gl += in->obsParams[(k*avDump)+kk].gl;
		      gb += in->obsParams[(k*avDump)+kk].gb;
		      az += in->obsParams[(k*avDump)+kk].az;
		      ze += in->obsParams[(k*avDump)+kk].ze;
		      el += in->obsParams[(k*avDump)+kk].el;
		      turn_utc_av += hms_turn(in->obsParams[(k*avDump)+kk].utc);
		      turn_local_av += hms_turn(in->obsParams[(k*avDump)+kk].local_time);
		      az_drive_rate += in->obsParams[(k*avDump)+kk].az_drive_rate;
		      ze_drive_rate += in->obsParams[(k*avDump)+kk].ze_drive_rate;
		      hourAngle += in->obsParams[(k*avDump)+kk].hourAngle;
		      paraAngle += in->obsParams[(k*avDump)+kk].paraAngle;
		      windDir += in->obsParams[(k*avDump)+kk].windDir;
		      windSpd  += in->obsParams[(k*avDump)+kk].windSpd;
		    }
		  
		}

	      out->wt[k*out->nchan + j] =   s_flag_wt; 
	      if (s_flag_wt > 0)
		out->flag[k*out->nchan + j] = 0;
	      else
		out->flag[k*out->nchan + j] = 1;
	      
	      if (sum==1)
		out->data[k*out->nchan*out->npol + p*out->nchan + j] = av;
	      else
		{
		  if (s_flag_wt > 0)
		    out->data[k*out->nchan*out->npol + p*out->nchan + j] = av/s_flag_wt;
		  else
		    out->data[k*out->nchan*out->npol + p*out->nchan + j] = av;
		}

	      if (p==0 && j==0)
		{
		  // Determine the observation parameter output
		  strcpy(out->obsParams[k].timedb,"UNSET"); // FIX ME
		  strcpy(out->obsParams[k].ut_date,"UNSET"); // FIX ME
		  
		  // SHOULD FIX WRAPS HERE --- FIX ME
		  out->obsParams[k].raDeg = raDeg/avDump; // SHOULD BE WEIGHTED ** FIX ME
		  out->obsParams[k].decDeg = decDeg/avDump; // SHOULD BE WEIGHTED ** FIX ME
		  turn_hms(turn_utc_av/avDump,out->obsParams[k].utc);
		  turn_hms(turn_local_av/avDump,out->obsParams[k].local_time);
		  
		  
		  turn_hms(out->obsParams[k].raDeg/360.,out->obsParams[k].raStr);
		  turn_dms(out->obsParams[k].decDeg/360.,out->obsParams[k].decStr);
		  
		  
		  out->obsParams[k].timeElapsed = timeElapsed/avDump; // SHOULD BE WEIGHTED *** FIX ME
		  out->obsParams[k].dtime = dtime; // Note this is just a summation (not an average)
		  out->obsParams[k].mjd = mjd/avDump; // SHOULD BE WEIGHTED *** FIX ME
		  		  
		  out->obsParams[k].raOffset = raOffset/avDump;
		  out->obsParams[k].decOffset = decOffset/avDump;
		  out->obsParams[k].gl = gl/avDump;
		  out->obsParams[k].gb = gb/avDump;
		  out->obsParams[k].az = az/avDump;
		  out->obsParams[k].ze = ze/avDump;
		  out->obsParams[k].el = el/avDump;
		  out->obsParams[k].az_drive_rate = az_drive_rate/avDump;
		  out->obsParams[k].ze_drive_rate = ze_drive_rate/avDump;
		  out->obsParams[k].hourAngle = hourAngle/avDump;
		  out->obsParams[k].paraAngle = paraAngle/avDump;
		  out->obsParams[k].windDir = windDir/avDump;
		  out->obsParams[k].windSpd = windSpd/avDump;
		}
	    }
	}
    }
}

void frequencyAverage(dataStruct *in,dataStruct *out,int nfreqAv,int sum,int meanMedian)
{
  int i,j,k,kk,p;
  int include;
  int avFreq,nfreq;
  double av,avFreqVal,avFreqNoWt,av2;
  float wt,wtVal;
  unsigned char flag;
  int nsum;
  int nMedian;
  double s_flag_wt;
  int flagVal;
  float medVal[nfreqAv],medWt[nfreqAv];
  float medVal2[nfreqAv];
  if (nfreqAv == 0)
    {nfreq = 1; avFreq = in->nchan;}
  else
    {nfreq = (int)((double)in->nchan/(double)nfreqAv); avFreq = nfreqAv;}


  out->nchan = nfreq;
  out->npol =  in->npol;
  out->ndump = in->ndump;
  out->nFreqDump = in->nFreqDump;
  out->astroCal = in->astroCal;
  
  printf("Averaging: %d frequency channels together giving %d frequency channels out of %d original channels\n",avFreq,nfreq,in->nchan);
  allocateMemory(out);
  memcpy(out->obsParams,in->obsParams,sizeof(sdhdf_obsParamsStruct)*in->ndump); 

  if (meanMedian==1)
    printf("Calculating mean\n");
  else
    printf("Calculating median\n");
  
  for (p=0;p<out->npol;p++)
    {
      for (k=0;k<out->ndump;k++)
	{
	  for (j=0;j<out->nchan;j++)
	    {
	      s_flag_wt=0;
	      nsum=0;
	      av=av2=0;
	      avFreqVal=avFreqNoWt = 0;
	      
	      wt=0;
	      flag = 1;
	      nMedian=0;
	      for (kk=0;kk<avFreq;kk++)
		{
		  wtVal     = in->wt[((long)k*in->nchan) + ((long)j*avFreq)+kk];
		  flagVal   = 1-in->flag[((long)k*in->nchan) + ((long)j*avFreq)+kk];
		  if (meanMedian==1)
		    {
		      if (in->astroCal==0)
			{
			  avFreqVal  += wtVal*flagVal*in->freq[((long)j*avFreq)+kk+(long)k*in->nchan];
			  avFreqNoWt += in->freq[((long)j*avFreq)+kk+(long)k*in->nchan];
			  av += wtVal*flagVal*in->data[(long)k*in->nchan*in->npol + (long)p*in->nchan + ((long)j*avFreq)+kk];		      
			}
		      else
			{
			  avFreqVal  += wtVal*flagVal*in->freq[((long)j*avFreq)+kk];
			  avFreqNoWt += in->freq[((long)j*avFreq)+kk];
			  av  += wtVal*flagVal*in->data[(long)k*in->nchan*in->npol + p*in->nchan + ((long)j*avFreq)+kk];
			  av2 += wtVal*flagVal*in->data2[(long)k*in->nchan*in->npol + p*in->nchan + ((long)j*avFreq)+kk];		      
			  //			  printf("Using: %g %g\n",in->data[k*in->nchan*in->npol + p*in->nchan + (j*avFreq)+kk],in->data2[k*in->nchan*in->npol + p*in->nchan + (j*avFreq)+kk]);
			}

		      s_flag_wt+= wtVal*flagVal;
		    }
		  else
		    {
		      if (in->astroCal == 0)
			avFreqNoWt += in->freq[((long)j*avFreq)+kk+(long)k*in->nchan];
		      else
			avFreqNoWt += in->freq[((long)j*avFreq)+(long)kk];
		      if (flagVal == 1 && wtVal != 0)
			{
			  medVal[nMedian]  = in->data[(long)k*in->nchan*in->npol + (long)p*in->nchan + ((long)j*avFreq)+kk];
			  if (in->astroCal!=0)
			    medVal2[nMedian] = in->data2[(long)k*in->nchan*in->npol + (long)p*in->nchan + ((long)j*avFreq)+kk];
			  medWt[nMedian]  = in->wt[((long)k*in->nchan) + ((long)j*avFreq)+kk];
			  s_flag_wt+= wtVal;
			  nMedian++;
			}
		    }
		}
	      if (meanMedian==1) // Mean
		{
		  out->wt[(long)k*out->nchan + j] =   s_flag_wt; 
		  if (s_flag_wt > 0)
		    out->flag[(long)k*out->nchan + j] = 0;
		  else
		    out->flag[(long)k*out->nchan + j] = 1;
		  
		  if (sum==1)
		    {
		      if (in->astroCal==0)
			out->data[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = av;
		      else
			{
			  out->data[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = av;
			  out->data2[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = av2;
			}
		    }
		  else
		    {
		      if (in->astroCal==0)
			{
			  if (s_flag_wt > 0)
			    out->data[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = av/s_flag_wt;
			  else
			    out->data[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = av;
			}
		      else
			{

			  if (s_flag_wt > 0)
			    {
			      out->data[(long)k*out->nchan*out->npol + (long)p*out->nchan + j]  = av/s_flag_wt;
			      out->data2[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = av2/s_flag_wt;
			    }
			  else
			    {
			      out->data[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = av;
			      out->data2[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = av2;
			    }
			}
		    }
		      
		}
	      else
		{
		  // Currently doing an unweighted median ** FIX ME **
		  if (nMedian==0)
		    {
		      out->data[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = 0;
		      if (in->astroCal==1)
			out->data2[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = 0;
		      out->flag[(long)k*out->nchan + j] = 1;
		      out->wt[(long)k*out->nchan + j] = 0;
		    }
		  else
		    {
		      out->data[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = quick_select_float(medVal,nMedian);
		      if (in->astroCal==1)
			out->data2[(long)k*out->nchan*out->npol + (long)p*out->nchan + j] = quick_select_float(medVal2,nMedian);
		      out->flag[(long)k*out->nchan + j] = 0;
		      out->wt[(long)k*out->nchan + j]   = s_flag_wt/(float)nMedian; 
		    }
		}

	      // FIX ME **
	      // NOTE: Not doing weighted average as need to set for a single spectral dump
	      if (p==0)
		{
		  if (in->astroCal == 0)
		    out->freq[j+(long)k*out->nchan] = avFreqNoWt/avFreq; // FIX ME: NOT PICKING UP CORRECT BIT OF FREQUENCY ARRAY
		  else
		    out->freq[j+(long)k*out->nchan] = avFreqNoWt/avFreq;
		}
		  //	      else
		//		{
		  //		  printf("WARNING: channel %d (%g) has a weighting of zero\n",j,avFreqNoWt/avFreq);
		//		  out->freq[j] = avFreqNoWt/avFreq;
		//		}
	    }
	}
    }
}

void changeFrequencyAxis(dataStruct *in,dataStruct *out,int bary_lsr,int regrid,sdhdf_attributes_struct *freqAttributes,int nFreqAttributes,int verbose,int showVel,double mjdSet)
{
  int i,j,k,kk,p;
  int include;
  int avFreq,nfreq;
  double av,avFreqVal,avFreqNoWt;
  float wt,wtVal;
  unsigned char flag;
  int nsum;
  double s_flag_wt;
  int flagVal;
  double origFreq,newFreq;
  double iX1,iY1,iX2,iY2,m,c;
  double iX,iY;
  sdhdf_eopStruct *eop;
  int nEOP=0;
  double *mjdVals,*raDeg,*decDeg;
  double *vOverC;
  char ephemName[1024] = "DE436.1950.2050";
    
  // Load in Earth observation parameters
  eop = (sdhdf_eopStruct *)malloc(sizeof(sdhdf_eopStruct)*MAX_EOP_LINES);
  nEOP = sdhdf_loadEOP(eop);

  out->nchan = in->nchan;
  out->npol =  in->npol;
  out->ndump = in->ndump;
  out->nFreqDump = in->nFreqDump;
  mjdVals = (double *)malloc(sizeof(double)*out->ndump);
  raDeg = (double *)malloc(sizeof(double)*out->ndump);
  decDeg = (double *)malloc(sizeof(double)*out->ndump);
  vOverC = (double *)malloc(sizeof(double)*out->ndump);

  allocateMemory(out);

  // Should just memcopy all the data across
  memcpy(out->data,in->data,sizeof(float)*out->nchan*out->npol*out->ndump);
  memcpy(out->wt,in->wt,sizeof(float)*out->nchan*out->ndump);
  memcpy(out->flag,in->flag,sizeof(unsigned char)*out->nchan*out->ndump);
  if (out->nFreqDump != out->ndump)
    {
      printf("ERROR in sdhdf_modify: cannot form lsr/bary corrections if nFreqDump != ndump\n");
      exit(1);
    }
  for (i=0;i<out->ndump;i++)
    {
      if (mjdSet > 0)
	mjdVals[i] = mjdSet;
      else
	mjdVals[i] = (in->obsParams[i].mjd + in->obsParams[i].fractional_mjd);
      raDeg[i] = in->obsParams[i].raDeg;
      decDeg[i] = in->obsParams[i].decDeg;
  }
  memcpy(out->obsParams,in->obsParams,sizeof(sdhdf_obsParamsStruct)*in->ndump); 
  
  sdhdf_calcVoverC(mjdVals,raDeg,decDeg,out->ndump,vOverC,(char *)"Parkes",ephemName,eop,nEOP,bary_lsr);

  
  //  printf("WARNING: HAVE A TEST HACK IN PLACE\n");
  //  printf("RESETTING: %d %d\n",out->ndump,out->nchan);
  /*  for (i=0;i<out->ndump;i++)
    {
      for (j=0;j<out->nchan;j++)
      {
	//	in->freq[i*out->nchan+j] = (float)((long double)j*274.3518519e-6 -7.5445833334L/1000.0L);
	in->freq[i*out->nchan+j] = (double)((long double)j*274.34842250e-6L);
	//	in->freq[i*out->nchan+j] += 1662.000000L;
	in->freq[i*out->nchan+j] += 1660.000000L;
	  if (j < 5) printf("HAVE: %.6f\n",in->freq[i*out->nchan+j]);
	}
	} */
  //    printf("******************\n");
  //    printf("CHANNELS 0 and 1 = %.6f %.6f %.6f %6f\n",in->freq[0],in->freq[1],1662.000122 -7.5445833334/1000.0,1662.000122 + 274.3518519e-6 -7.5445833334/1000.0);
  for (i=0;i<out->ndump;i++)
    {
      if (verbose==1 || showVel==1)
	printf("velocity correction: mjd = %.6f idump = %d vOverC = %.8f ra/dec = %g/%g\n",mjdVals[i],i,vOverC[i],raDeg[i],decDeg[i]);
      if (regrid == 1) // (note 2 = regrid, 1 = don't regrid)
	{
	  for (j=0;j<out->nchan;j++)
	    {
	      newFreq = (double)in->freq[i*out->nchan+j]*(1.0L-vOverC[i]);
	      out->freq[i*out->nchan+j] = (double)newFreq;	     
	    }
	}
      else
	{
	  double newFreq,oldFreq;
	  int kk;
	  int deltaI;
	  double fracDeltaI;
	  double df;
	  double chanbw;

	  chanbw = in->freq[i*in->nchan+1] - in->freq[i*in->nchan];
	  //	  printf("CALCULATING CHANBW %.6f %.6f %.6f\n",in->freq[i*in->nchan+1],in->freq[i*in->nchan],chanbw);
	  //	  printf("===================================\n");
//	    printf("HAVE HACK IN PLACE - DO NOT USE\n");
//		 printf("===================================\n");
	  for (j=0;j<out->nchan;j++)
	    {
	      //	      in->freq[i*out->nchan+j] -= 10/1000.0; // FOR PAF ****
	   
	      newFreq = in->freq[i*out->nchan+j]*(1.0L-vOverC[i]);
	      oldFreq = in->freq[i*out->nchan+j];
	      out->freq[i*out->nchan+j] = oldFreq;
	      // Should first check if we're still in topocentric frequencies -- DO THIS ***  FIX ME
	      // FIX ME: this assumes equally sampled frequency channels
	      df = (double)(oldFreq-newFreq)/(double)chanbw; // Check if frequency channelisation changes - e.g., at subband boundaries ** FIX ME
	      if (df > 0)
		{
		  deltaI = (int)df;
		  fracDeltaI = df-deltaI;
		}
	      else
		{
		  // CHECK THIS -- FIX ME
		  //		  deltaI = -(int)(fabs(df)+0.5);
		  deltaI = -(int)(fabs(df)+1);
		  //		  fracDeltaI = deltaI - df;  // CHECK MINUS SIGN
		  fracDeltaI = df-deltaI;
		}
	      if (verbose==1)
		printf("Regridding df = %.6f %.6f %.6f %.6f chanbw = %g nchan1 = %d nchan2 = %d deltaI = %d fracDeltaI = %g\n",df,oldFreq,newFreq,oldFreq-newFreq,chanbw,in->nchan,out->nchan,deltaI,fracDeltaI);
	      if (j + deltaI >= 0 && j+deltaI < in->nchan-1)
		{
		  for (kk=0;kk<in->npol;kk++)
		    {
		      if (kk==4) // 5 pol
			out->data[i*out->nchan*out->npol + kk*out->nchan + j] = 0;
		      else
			{
			  //			  			  deltaI = 200;
			  //			  		  fracDeltaI = 0;
			  iX1 = 0; iX2 = 1;
			  //			  printf("Searching for %d %d\n",i*in->nchan*in->npol + kk*in->nchan + (j + deltaI),j);
			  iY1 = in->data[(long int)i*(long int)(in->nchan)*(long int)in->npol + (long int)kk*(long int)in->nchan + (long int)(j + deltaI)];
			  iY2 = in->data[(long int)i*(long int)(in->nchan)*(long int)in->npol + (long int)kk*(long int)in->nchan + (long int)(j + deltaI + 1)];
			  m   = (iY2-iY1)/(iX2-iX1);
			  c   = iY1;

			  iX  = fracDeltaI;
			  iY  = m*iX+c;

			  out->data[(long int)i*(long int)out->nchan*(long int)out->npol + (long int)kk*(long int)out->nchan + (long int)j] = iY;
			  //			  out->data[i*out->nchan*out->npol + kk*out->nchan + j] = iY1;
			  //
			}
		    }
		}
	    }
	}
    }
  // Fix up the attributes
  for (i=0;i<nFreqAttributes;i++)
    {
      if (strcmp(freqAttributes[i].key,"FRAME")==0)
	{
	  if (bary_lsr==1)
	    strcpy(freqAttributes[i].value,"barycentric");
	  else if (bary_lsr==2 || bary_lsr==3)
	    strcpy(freqAttributes[i].value,"LSR");				    
	}

    }
  free(eop);
  free(mjdVals);
  free(vOverC);
  free(raDeg);
  free(decDeg);
}


void allocateMemory(dataStruct *in)
{
  int k;
  if (in->memoryAllocated == 0)  // SHOULD CHECK NPOL HERE
    {
      if (!(in->data = (float *)malloc(sizeof(float)*in->nchan*in->ndump*in->npol)))
	{
	  printf("ERROR: Unable to allocate memory in sdhdf_modify\n");
	  exit(1);
	}
      if (in->astroCal==1)
	  in->data2 = (float *)malloc(sizeof(float)*in->nchan*in->ndump*in->npol);
      in->freq = (double *)malloc(sizeof(double)*in->nchan*in->nFreqDump);
      in->wt   = (float *)malloc(sizeof(float)*in->nchan*in->ndump);
      in->flag = (unsigned char *)malloc(sizeof(unsigned char)*in->nchan*in->ndump);
      in->obsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*in->ndump);      
      in->memoryAllocated = 1;
    }
  else
    {
      in->data = (float *)realloc(in->data,sizeof(float)*in->nchan*in->ndump*in->npol);
      if (in->astroCal==1)
	  in->data2 = (float *)realloc(in->data2,sizeof(float)*in->nchan*in->ndump*in->npol);

      in->freq = (double *)realloc(in->freq,sizeof(double)*in->nchan*in->nFreqDump);
      //      in->freq = (float *)realloc(in->freq,sizeof(float)*in->nchan);
      in->wt = (float *)realloc(in->wt,sizeof(float)*in->nchan*in->ndump);
      in->flag = (unsigned char *)realloc(in->flag,sizeof(unsigned char)*in->nchan*in->ndump);
      in->obsParams = (sdhdf_obsParamsStruct *)realloc(in->obsParams,sizeof(sdhdf_obsParamsStruct)*in->ndump);      
    }
}

void polarisationAverage(dataStruct *in,dataStruct *out,int sum)
{
  int i,j,k,kk,p;
  int include;
  int avDump,ndump;
  double av;
  float wt,wtVal;
  unsigned char flag;
  int nsum;
  double s_flag_wt;
  int flagVal;

  out->ndump = in->ndump;
  out->nchan = in->nchan;
  out->npol =  1;
  out->nFreqDump = in->nFreqDump;
  
  allocateMemory(out);
  
  memcpy(out->obsParams,in->obsParams,sizeof(sdhdf_obsParamsStruct)*in->ndump); 
  memcpy(out->freq,in->freq,sizeof(double)*out->nchan*out->nFreqDump);
  for (j=0;j<out->nchan;j++)
    {
      for (k=0;k<out->ndump;k++)
	{
	  for (p=0;p<out->npol;p++)
	    {	      
	      av=0;
	      wt=0;
	      flag = 1;
	      nsum=0;
	      s_flag_wt = 0;

	      av = in->data[k*in->nchan*in->npol +j] + in->data[k*in->nchan*in->npol + in->nchan +j];
	      if (sum==2) av/=2.0; // Average
	      
	      out->wt[k*out->nchan + j] =   in->wt[k*out->nchan+j]; 
	      out->flag[k*out->nchan + j] =   in->flag[k*out->nchan+j]; 
	      out->data[k*out->nchan + j] = av;
	      
	    }
	}
    }
}


void scaleValues(dataStruct *in,dataStruct *out,int iband, int useBand, int multDiv, float aaScale, float bbScale)
{
  int i,j,k,kk,p;
  double v;
  float wt,wtVal;
  unsigned char flag;
  int flagVal;

  out->ndump = in->ndump;
  out->nchan = in->nchan;
  out->npol =  in->npol;
  out->nFreqDump = in->nFreqDump;
  
  allocateMemory(out);
  
  memcpy(out->obsParams,in->obsParams,sizeof(sdhdf_obsParamsStruct)*in->ndump); 
  memcpy(out->freq,in->freq,sizeof(double)*out->nchan*out->nFreqDump);
  if (useBand != iband)
    memcpy(out->data,in->data,sizeof(float)*out->nchan*out->ndump*out->npol);
  else
    {
      printf("Scaling band %d by %g %g\n",iband,aaScale,bbScale);
      for (j=0;j<out->nchan;j++)
	{
	  for (k=0;k<out->ndump;k++)
	    {
	      for (p=0;p<out->npol;p++)
		{	      
		  v = in->data[k*in->nchan*in->npol + p*in->nchan + j];
		  if (multDiv==1)
		    {
		      if (p==0)      v/=aaScale;
		      else if (p==1) v/=bbScale;		  // No clear solution to the cross terms -- FIX ME
		    }
		  else
		    {
		      if (p==0)      v*=aaScale;
		      else if (p==1) v*=bbScale;		  // No clear solution to the cross terms -- FIX ME
		    }
		  if (p==0)
		    {
		      out->wt[k*out->nchan + j] =   in->wt[k*out->nchan+j]; 
		      out->flag[k*out->nchan + j] =   in->flag[k*out->nchan+j]; 
		    }
		  out->data[k*out->nchan + p*out->nchan + j] = v;
		  
		}
	    }
	}
    }
}
