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
  printf("sdhdf_extractDump\n\n");
  printf("Routine to extract spectral dumps from an sdhdf file\n\n");
  printf("-b <num>       Beam number\n");
  printf("-d <num>       Spectral dump number to extract (note this can be used multiple times)\n");
  printf("-e <ext>       Output file extension\n");
  printf("-h             This help\n");
  printf("\n\n");
  printf("Filenames are listed on the command line. For example\n\n");
  printf("sdhdf_extractDump -e extract -d 0 -d 5 -d 8 *.hdf\n");
}

int main(int argc,char *argv[])
{
  int ii,i,j,k,kk,jj,l,nchan,totNchan,b,nFreqDumps;
  char fname[MAX_FILES][64];
  int nFiles=0;
  char ext[MAX_STRLEN]="extract";
  char oname[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*outFile;
  herr_t status;
  int selectDump[MAX_BANDS];
  float *outVals,*freqVals,*inData;
  int  nSelectDumps=0;
  int  copySD=0;
  float fSelect0[MAX_BANDS];
  float fSelect1[MAX_BANDS];
  int zoomBand=0;
  sdhdf_obsParamsStruct *outObsParams;
  sdhdf_bandHeaderStruct *outBandParams,*outCalBandParams,*inBandParams,*inCalBandParams;
  int nBand=0;
  char zoomLabel[MAX_STRLEN];
  char groupName[MAX_STRLEN];
  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;
  int selectBeam=-1;
  
  long npol,ndump,out_ndump,out_ndump_num;
  float *out_data,*out_freq;
  
  strcpy(oname,"sdhdf_extract_output.hdf");
  
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
  
  
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-e")==0)
	strcpy(ext,argv[++i]);
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-d")==0)
	sscanf(argv[++i],"%d",&selectDump[nSelectDumps++]);
      else if (strcmp(argv[i],"-b")==0)
	sscanf(argv[++i],"%d",&selectBeam);
      else
	{
	  strcpy(fname[nFiles],argv[i]);
	  nFiles++;
	}
    }

  for (ii=0;ii<nFiles;ii++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_initialiseFile(outFile);

      sdhdf_formOutputFilename(fname[ii],ext,oname);
            
      sdhdf_openFile(fname[ii],inFile,1);
      sdhdf_openFile(oname,outFile,3);
      
      if (inFile->fileID!=-1) // Did we successfully open the file?
	{
	  sdhdf_loadMetaData(inFile);
	  printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname[ii],inFile->primary[0].utc0,
		 inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[0].source,inFile->primary[0].telescope,
		 inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[0].nBand);

	  sdhdf_allocateBeamMemory(outFile,inFile->nBeam);
	  printf("Allocated memory\n");
	  for (b=0;b<inFile->nBeam;b++)
	    {
	      if (selectBeam == -1 || b == selectBeam)
		{
		  printf("Processing beam %d\n",b);
		  inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);      
		  sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);
		  
		  for (i=0;i<inFile->beam[b].nBand;i++)
		    {
		      nchan = inFile->beam[b].bandHeader[i].nchan;
		      nFreqDumps = inFile->beam[b].bandData[i].astro_data.nFreqDumps;
		      npol  = inFile->beam[b].bandHeader[i].npol;
		      ndump = inFile->beam[b].bandHeader[i].ndump;
		      printf("Processing subband %d (number of spectral dumps = %d)\n",i,(int)ndump);
		      out_ndump = nSelectDumps;
		      outObsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*out_ndump);      
		      for (k=0;k<out_ndump;k++)
			sdhdf_copySingleObsParams(inFile,b,i,k,&outObsParams[k]);
		      
		      
		      out_data = (float *)malloc(sizeof(float)*nchan*npol*out_ndump);
		      // Should use nFreqDumps here -- FIX ME
		      out_freq = (float *)malloc(sizeof(float)*nchan);
		      
		      sdhdf_loadBandData(inFile,b,i,1);
		      for (j=0;j<nchan;j++)
			out_freq[j] = inFile->beam[b].bandData[i].astro_data.freq[j]; // FIX ME: DUMP FOR FREQ AXIS
		      out_ndump_num=0;
		      for (j=0;j<ndump;j++)
			{
			  copySD=0;
			  for (jj=0;jj<nSelectDumps;jj++)
			    {
			      if (selectDump[jj]==j)
				{
				  printf("Copying subband %d\n",j);
				  for (kk=0;kk<nchan;kk++)
				    {
				      if (npol==1)
					out_data[kk+out_ndump_num*nchan]         = inFile->beam[b].bandData[i].astro_data.pol1[kk+j*nchan];
				      else if (npol==2)
					{
					  out_data[kk+out_ndump_num*nchan*2]         = inFile->beam[b].bandData[i].astro_data.pol1[kk+j*nchan];
					  out_data[kk+nchan+out_ndump_num*nchan*2]   = inFile->beam[b].bandData[i].astro_data.pol2[kk+j*nchan];
					}
				      else if (npol==4)
					{
					  out_data[kk        +out_ndump_num*nchan*4] = inFile->beam[b].bandData[i].astro_data.pol1[kk+j*nchan];
					  out_data[kk+nchan  +out_ndump_num*nchan*4] = inFile->beam[b].bandData[i].astro_data.pol2[kk+j*nchan];
					  out_data[kk+2*nchan+out_ndump_num*nchan*4] = inFile->beam[b].bandData[i].astro_data.pol3[kk+j*nchan];
					  out_data[kk+3*nchan+out_ndump_num*nchan*4] = inFile->beam[b].bandData[i].astro_data.pol4[kk+j*nchan];
					}
				    }
				  out_ndump_num++;
				}
			    }
			}
		      //		  sdhdf_writeBandHeader(outFile,outBandParams,b,nSelectBands,1);
		      //		  sdhdf_writeBandHeader(outFile,outCalBandParams,b,nSelectBands,2);
		      printf("Writing spectral data\n");

		      // For now just picking one frequency dump for output! ** FIX ME
		      sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label,b,i,out_data,out_freq,1,nchan,1,npol,out_ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
		      printf("Writing obs params\n");
		      sdhdf_writeObsParams(outFile,inFile->beam[b].bandHeader[i].label,inFile->beamHeader[b].label,i,outObsParams,out_ndump,1,inFile->primary[0].hdr_defn_version);
		      printf("Releasing data\n");
		      sdhdf_releaseBandData(inFile,b,ii,1);		  

		      free(outObsParams);
		      inBandParams[i].ndump = out_ndump;
		    }
		  printf("Writing band header\n");
		  sdhdf_writeBandHeader(outFile,inBandParams,inFile->beamHeader[b].label,inFile->beam[b].nBand,1,inFile->primary[0].hdr_defn_version);
		  free(inBandParams);
		}
	    }
	  printf("Freeing\n");
	    

	      

	  printf("Writing other tables\n");
	  // Copy other primary tables
	  if (selectBeam > -1)
	    {
	      sdhdf_beamHeaderStruct *outBeamHeader;

	      outBeamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*1);
	      strcpy(outBeamHeader->label,inFile->beamHeader[selectBeam].label);
	      strcpy(outBeamHeader->source,inFile->beamHeader[selectBeam].source);
	      outBeamHeader->nBand = inFile->beamHeader[selectBeam].nBand;
	      inFile->primary[0].nbeam = 1;
	      sdhdf_writePrimaryHeader(outFile,inFile->primary);
	      sdhdf_writeBeamHeader(outFile,outBeamHeader,1,inFile->primary[0].hdr_defn_version);
	      free(outBeamHeader);
	    }
	  sdhdf_writeHistory(outFile,inFile->history,inFile->nHistory);
	  printf("Copy remainder\n");
	  //	  sdhdf_copyRemainder(inFile,outFile,1); // 1 = don't write extra beam information
	  sdhdf_copyRemainder(inFile,outFile,0); // 1 = don't write extra beam information
	  printf("Closing files\n");
	  sdhdf_closeFile(inFile);
	  sdhdf_closeFile(outFile);
	}
    }

  free(inFile);
  free(outFile);
}



