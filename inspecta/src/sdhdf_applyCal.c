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

// Code to apply a calibration solution to a data file
// for flux and polarisation calibration.  The output is a new, calibrated, file.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"


int main(int argc,char *argv[])
{
  int i,j,k,b,ii,i0,i1;
  sdhdf_fileStruct *inFile;
  sdhdf_fileStruct *outFile;
  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;

  float *cal_freq;
  float *tsys,*avTsys_AA,*avTsys_BB,fval;
  float *outData;
  float af0,af1,cbw;
  double sumAA,sumBB;
  double scaleAA,scaleBB;
  char fname[MAX_STRLEN];
  char oname[MAX_STRLEN];
  char extension[MAX_STRLEN];
  int nBeam;
  int nchanAst,npolAst,ndumpAst;
  
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
	strcpy(extension,argv[++i]);
      else
	{
	  strcpy(fname,argv[i]);
	  sdhdf_initialiseFile(inFile);
	  sdhdf_initialiseFile(outFile);

	  sdhdf_openFile(fname,inFile,1);
	  sprintf(oname,"%s.%s",fname,extension);
	  sdhdf_openFile(oname,outFile,3);

	  if (inFile->fileID!=-1) // Did we successfully open the file?
	    {
	      sdhdf_loadMetaData(inFile);

	      
	      printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname,inFile->primary[0].utc0,inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[0].source,inFile->primary[0].telescope,inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[0].nBand);
	      
	      nBeam = inFile->nBeam;
	      sdhdf_allocateBeamMemory(outFile,nBeam);
	      for (b=0; b<nBeam; b++)
		{
		  for (ii=0;ii<inFile->beam[b].nBand;ii++)
		    {
		      printf("Processing band %d, ndumpCal = %d\n",ii,inFile->beam[b].calBandHeader[ii].ndump);


		      // Calculate average Tsys values
		      tsys = (float *)calloc(sizeof(float),inFile->beam[b].calBandHeader[ii].ndump*inFile->beam[b].calBandHeader[ii].nchan*2); // * 2 = AA and BB
		      avTsys_AA = (float *)calloc(sizeof(float),inFile->beam[b].calBandHeader[ii].nchan);
		      avTsys_BB = (float *)calloc(sizeof(float),inFile->beam[b].calBandHeader[ii].nchan);
		      outData   = (float *)calloc(sizeof(float),inFile->beam[b].bandHeader[ii].nchan*inFile->beam[b].bandHeader[ii].ndump*inFile->beam[b].bandHeader[ii].npol);
		      sdhdf_loadBandData(inFile,b,ii,2);
		      sdhdf_loadBandData(inFile,b,ii,3);
		      sdhdf_loadCalProc(inFile,b,ii,"cal_proc_tsys",tsys);

		      for (k=0;k<inFile->beam[b].calBandHeader[ii].nchan;k++)
			{
			  for (j=0;j<inFile->beam[b].calBandHeader[ii].ndump;j++)
			    {
			      avTsys_AA[k] += tsys[k+j*2*inFile->beam[b].calBandHeader[ii].nchan];
			      avTsys_BB[k] += tsys[k+j*2*inFile->beam[b].calBandHeader[ii].nchan+inFile->beam[b].calBandHeader[ii].nchan];
			    }
			  avTsys_AA[k]/=inFile->beam[b].calBandHeader[ii].ndump;
			  avTsys_BB[k]/=inFile->beam[b].calBandHeader[ii].ndump;

			  //			  avTsys_AA[k] = 30.1954;
			  
			  //			  avTsys_AA[k]+=0.15*avTsys_AA[k];  // Used to scale Tsys
			  //			  avTsys_BB[k]+=0.15*avTsys_BB[k];
			}

		      // Now process the astronomy data
		      sdhdf_loadBandData(inFile,b,ii,1);
		      // FIX ME: using [0] for frequency dump
		      af0 = inFile->beam[b].bandData[ii].astro_data.freq[0];
		      af1 = inFile->beam[b].bandData[ii].astro_data.freq[inFile->beam[b].bandHeader[ii].nchan-1];

		      // FIXME: Assuming positive bandwidth
		      cbw = inFile->beam[b].bandData[ii].cal_on_data.freq[1] - inFile->beam[b].bandData[ii].cal_on_data.freq[0];
		      
		      for (k=0;k<inFile->beam[b].calBandHeader[ii].nchan;k++)
			{

			  // FIX ME: using [0] for frequency dump
			  fval = inFile->beam[b].bandData[ii].cal_on_data.freq[k];
			  //			  i0 = (int)(inFile->beam[b].bandHeader[ii].nchan*(fval-af0)/(af1-af0));
			  //			  i1 = (int)(inFile->beam[b].bandHeader[ii].nchan*(fval+cbw-af0)/(af1-af0));

			  //			  i0 = (int)(inFile->beam[b].bandHeader[ii].nchan*(fval-cbw/2.0-af0)/(af1-af0));
			  i0 = (int)(inFile->beam[b].bandHeader[ii].nchan*(fval-cbw-af0)/(af1-af0));
			  if (i0<0) i0=0;
			  //			  i1 = (int)(inFile->beam[b].bandHeader[ii].nchan*(fval+cbw/2.0-af0)/(af1-af0));
			  i1 = (int)(inFile->beam[b].bandHeader[ii].nchan*(fval-af0)/(af1-af0));
			  if (i1 > inFile->beam[b].bandHeader[ii].nchan) i1 = inFile->beam[b].bandHeader[ii].nchan;
			  sumAA=sumBB=0;
 			  for (j=i0;j<i1;j++)
			    {
			      sumAA+=inFile->beam[b].bandData[ii].astro_data.pol1[j];
			      sumBB+=inFile->beam[b].bandData[ii].astro_data.pol2[j];
			    }
			  sumAA/=(double)(i1-i0);
			  sumBB/=(double)(i1-i0);
			  scaleAA = sumAA/avTsys_AA[k];
			  scaleBB = sumBB/avTsys_BB[k];

			  //			  scaleAA=1;
			  //			  scaleBB=1;
						  
			  nchanAst = inFile->beam[b].bandHeader[ii].nchan;
			  npolAst  = inFile->beam[b].bandHeader[ii].npol;
			  ndumpAst  = inFile->beam[b].bandHeader[ii].ndump;

 			  for (j=i0;j<i1;j++)
			    {
			      outData[j+nchanAst*npolAst*0]              = inFile->beam[b].bandData[ii].astro_data.pol1[j]/scaleAA; // 0 -> DUMP (FIX ME)
			      outData[j+nchanAst*npolAst*0 + nchanAst*1] = inFile->beam[b].bandData[ii].astro_data.pol2[j]/scaleBB; // 0 -> DUMP (FIX ME)
			    }
			  // FIX ME: using [0] for frequency dump
			  printf("avcal: %.4f %g %g %g %g %g %g %d %d %g %g %g\n",inFile->beam[b].bandData[ii].cal_on_data.freq[k],avTsys_AA[k],avTsys_BB[k],sumAA,sumBB,scaleAA,scaleBB,i0,i1,af0,af1,cbw);

			}
		      
		      for (k=0;k<inFile->beam[b].bandHeader[ii].nchan;k++)
			{
			  // FIX ME: using [0] for frequency dump
			  printf("data: %.4f %g %g %g %g\n",inFile->beam[b].bandData[ii].astro_data.freq[k],
				 inFile->beam[b].bandData[ii].astro_data.pol1[k],
				 inFile->beam[b].bandData[ii].astro_data.pol2[k],
				 outData[k+inFile->beam[b].bandHeader[ii].nchan*0],outData[k+inFile->beam[b].bandHeader[ii].nchan*0+inFile->beam[b].bandHeader[ii].nchan*1]);
			}

		      sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[ii].label,b,ii,outData,inFile->beam[b].bandData[ii].astro_data.freq,inFile->beam[b].bandData[ii].astro_data.nFreqDumps,inFile->beam[b].bandHeader[ii].nchan,1,inFile->beam[b].bandHeader[ii].npol,inFile->beam[b].bandHeader[ii].ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
		    
		      free(tsys);
		      free(avTsys_AA);
		      free(avTsys_BB);
		      free(outData);
		    }
		}	     
	      sdhdf_copyRemainder(inFile,outFile,0);
	      sdhdf_closeFile(inFile);
	      sdhdf_closeFile(outFile);

	    }
	}
    }

  free(inFile);
  free(outFile);

}
