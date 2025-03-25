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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"

#define VNUM "v0.1"


int main(int argc,char *argv[])
{
  int i,j,k,p,ii,b;
  char fnameOn[MAX_FILES][MAX_STRLEN];
  char fnameOff[MAX_FILES][MAX_STRLEN];
  char defineOn[MAX_STRLEN];
  char defineOff[MAX_STRLEN];
  int nFiles=0;
  sdhdf_fileStruct *onFile;
  sdhdf_fileStruct *offFile;
  sdhdf_bandHeaderStruct *inBandParams;  
  //  spectralDumpStruct spectrumOn;
  //  spectralDumpStruct spectrumOff;
  float *out_data;
  float on_pol1,off_pol1;
  float on_pol2,off_pol2;
  float on_pol3,off_pol3;
  float on_pol4,off_pol4;
  float *freq;
  float *scalAA,*scalBB;
  float *av_calOn_srcOff_AA,*av_calOff_srcOff_AA;
  float *av_calOn_srcOff_BB,*av_calOff_srcOff_BB;

  float *av_calOn_srcOn_AA,*av_calOff_srcOn_AA;
  float *av_calOn_srcOn_BB,*av_calOff_srcOn_BB;
  long nchan,npol,ndump1,ndump2;
  float sclA=1.0;
  float sclB=1.0;

  double sModel;

  char outName[1024];
  FILE *fout;
  
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

  
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-on")==0)
	strcpy(fnameOn[0],argv[++i]);	
      else if (strcmp(argv[i],"-off")==0)
	strcpy(fnameOff[0],argv[++i]);	
      else if (strcmp(argv[i],"-sclA")==0)
	sscanf(argv[++i],"%f",&sclA);
      else if (strcmp(argv[i],"-sclB")==0)
	sscanf(argv[++i],"%f",&sclB);
      else if (strcmp(argv[i],"-h")==0) // Should actually free memory
	exit(1);
    }
  sprintf(outName,"scal.dat");
  fout = fopen(outName,"w");
  nFiles=1;
  
  for (ii=0;ii<nFiles;ii++)
    {
      sdhdf_initialiseFile(onFile);
      sdhdf_initialiseFile(offFile);

      sdhdf_openFile(fnameOn[ii],onFile,1);
      sdhdf_openFile(fnameOff[ii],offFile,1);
      sdhdf_loadMetaData(onFile);
      sdhdf_loadMetaData(offFile);

      for (b=0;b<onFile->nBeam;b++)
	{
	  for (i=0;i<onFile->beam[b].nBand;i++)
	    {	  
	      nchan = onFile->beam[b].calBandHeader[i].nchan;
	      npol  = onFile->beam[b].calBandHeader[i].npol;
	      // FIX ME: SHOULD CHECK IF NCHAN IS THE SAME FOR THE ON AND OFF FILE
	      ndump1  = onFile->beam[b].calBandHeader[i].ndump;
	      ndump2  = offFile->beam[b].calBandHeader[i].ndump;

	      freq = (float *)malloc(sizeof(float)*nchan);      
	      sdhdf_loadBandData(onFile,b,i,2);
	      sdhdf_loadBandData(onFile,b,i,3);
	      sdhdf_loadBandData(offFile,b,i,2);
	      sdhdf_loadBandData(offFile,b,i,3);


	      for (k=0;k<nchan;k++)
		freq[k] = onFile->beam[b].bandData[i].cal_on_data.freq[k];  // FIX ME -- FREQ AXIS FOR DUMP

	      av_calOn_srcOn_AA = (float *)calloc(sizeof(float),nchan);
	      av_calOn_srcOn_BB = (float *)calloc(sizeof(float),nchan);
	      av_calOff_srcOn_AA = (float *)calloc(sizeof(float),nchan);
	      av_calOff_srcOn_BB = (float *)calloc(sizeof(float),nchan);

	      av_calOn_srcOff_AA = (float *)calloc(sizeof(float),nchan);
	      av_calOn_srcOff_BB = (float *)calloc(sizeof(float),nchan);
	      av_calOff_srcOff_AA = (float *)calloc(sizeof(float),nchan);
	      av_calOff_srcOff_BB = (float *)calloc(sizeof(float),nchan);

	      scalAA = (float *)calloc(sizeof(float),nchan);
	      scalBB = (float *)calloc(sizeof(float),nchan);

	      for (j=0;j<ndump1;j++)
		{	      
		  for (k=0;k<nchan;k++)
		    {
		      av_calOn_srcOn_AA[k] += onFile->beam[b].bandData[i].cal_on_data.pol1[k+j*nchan]/ndump1;
		      av_calOn_srcOn_BB[k] += onFile->beam[b].bandData[i].cal_on_data.pol2[k+j*nchan]/ndump1;
		      av_calOff_srcOn_AA[k] += onFile->beam[b].bandData[i].cal_off_data.pol1[k+j*nchan]/ndump1;
		      av_calOff_srcOn_BB[k] += onFile->beam[b].bandData[i].cal_off_data.pol2[k+j*nchan]/ndump1;
		    }
		}

	      for (j=0;j<ndump2;j++)
		{	      
		  for (k=0;k<nchan;k++)
		    {
		      av_calOn_srcOff_AA[k] += offFile->beam[b].bandData[i].cal_on_data.pol1[k+j*nchan]/ndump2;
		      av_calOn_srcOff_BB[k] += offFile->beam[b].bandData[i].cal_on_data.pol2[k+j*nchan]/ndump2;
		      av_calOff_srcOff_AA[k] += offFile->beam[b].bandData[i].cal_off_data.pol1[k+j*nchan]/ndump2;
		      av_calOff_srcOff_BB[k] += offFile->beam[b].bandData[i].cal_off_data.pol2[k+j*nchan]/ndump2;
		    }
		}
	      
	      for (k=0;k<nchan;k++)
		{
		  // Set for 1934 -- FIX ME
		  // Note model in Jy and so output is Scal (**not Tcal**)
		  sModel = pow(10,-30.7667+26.4908*log10(freq[k]) - 7.0977*pow(log10(freq[k]),2) + 0.605334*pow(log10(freq[k]),3));
		  scalAA[k] = sModel*(av_calOn_srcOff_AA[k]-av_calOff_srcOff_AA[k])/(av_calOff_srcOn_AA[k]-av_calOff_srcOff_AA[k]);
		  scalBB[k] = sModel*(av_calOn_srcOff_BB[k]-av_calOff_srcOff_BB[k])/(av_calOff_srcOn_BB[k]-av_calOff_srcOff_BB[k]);
		  fprintf(fout,"%.6f %g %g\n",freq[k],scalAA[k],scalBB[k]);
		  printf("%.6f %g %g %g %g %g %g %g %g %g %g %g\n",
			 freq[k],
			 scalAA[k],
			 scalBB[k],
			 sModel,
			 av_calOn_srcOn_AA[k],
			 av_calOn_srcOn_BB[k],
			 av_calOff_srcOn_AA[k],
			 av_calOff_srcOn_BB[k],
			 av_calOn_srcOff_AA[k],
			 av_calOn_srcOff_BB[k],
			 av_calOff_srcOff_AA[k],
			 av_calOff_srcOff_BB[k]);
		}
	      sdhdf_releaseBandData(onFile,b,i,1);
	      sdhdf_releaseBandData(offFile,b,i,1);
	    
	      free(freq);
	      free(av_calOn_srcOn_AA);
	      free(av_calOff_srcOn_AA);
	      free(av_calOn_srcOn_BB);
	      free(av_calOff_srcOn_BB);
	      free(av_calOn_srcOff_AA);
	      free(av_calOff_srcOff_AA);
	      free(av_calOn_srcOff_BB);
	      free(av_calOff_srcOff_BB);
	      free(scalAA);
	      free(scalBB);

	    }
	}
      
      
      sdhdf_closeFile(onFile);
      sdhdf_closeFile(offFile);
    }	
  fclose(fout);    
  free(onFile);
  free(offFile);
}

