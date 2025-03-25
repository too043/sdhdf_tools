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
//
// Usage:
// sdhdf_calProc <filename.hdf> <filename2.hdf> ...
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"


int main(int argc,char *argv[])
{
  int ii,i,j,k,l,b;
  char fname[MAX_FILES][MAX_STRLEN];
  int nFiles=0;
  char oname[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*outFile;
  sdhdf_tcal_struct *tcalData;
  float *tsys,*gain,*phase;
  double tcalA,tcalB;
  double on1,off1,on2,off2;
  double on3,off3,on4,off4;
  long nTcal;
  double s1,s2,s3,s4;
  double d1,d2,d3,d4;
  double avTsys;
  
  tcalData = (sdhdf_tcal_struct *)malloc(sizeof(sdhdf_tcal_struct)*3328); // 3328 = number of channels in Tcal files  
  // HARDCODED
  printf("DO NOT USE THIS CODE\n");
  exit(1);
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }


  for (i=1;i<argc;i++)
    {
      strcpy(fname[nFiles++],argv[i]);
    }


  for (ii=0;ii<nFiles;ii++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_initialiseFile(outFile);
      
      sprintf(oname,"%s.%s",fname[ii],"calProc");
      printf("Processing %s\n",fname[ii]);
      sdhdf_openFile(fname[ii],inFile,1);
      sdhdf_openFile(oname,outFile,3);
      
      if (inFile->fileID!=-1) // Did we successfully open the file?
	{
	  sdhdf_loadMetaData(inFile);		  
	  sdhdf_copyRemainder(inFile,outFile,0);

	  for (b=0;b<inFile->nBeam;b++)
	    {
	      // Load and write the data
	      for (i=0;i<inFile->beam[b].nBand;i++)
		{
		  printf(" ... subband %d/%d\n",i,inFile->beam[b].nBand);
		  // 2 = Pol A and Pol B
		  tsys = (float *)malloc(sizeof(float)*inFile->beam[b].calBandHeader[i].nchan*2*inFile->beam[b].calBandHeader[i].ndump);
		  
		  // 1 pol for gain and phase
		  gain = (float *)malloc(sizeof(float)*inFile->beam[b].calBandHeader[i].nchan*inFile->beam[b].calBandHeader[i].ndump);
		  phase = (float *)malloc(sizeof(float)*inFile->beam[b].calBandHeader[i].nchan*inFile->beam[b].calBandHeader[i].ndump); 
		  

		  sdhdf_loadBandData(inFile,b,i,2);
		  sdhdf_loadBandData(inFile,b,i,3);
		  		  
		  // Now write the cal_proc group and table
		  for (j=0;j<inFile->beam[b].calBandHeader[i].ndump;j++)
		    {
		      for (k=0;k<inFile->beam[b].calBandHeader[i].nchan;k++)
			{

			  // FIX ME: using [0] for frequency dump
			  sdhdf_get_tcal(tcalData,nTcal,inFile->beam[b].bandData[i].cal_on_data.freq[k],&tcalA,&tcalB);

			  on1   = inFile->beam[b].bandData[i].cal_on_data.pol1[k+j*inFile->beam[b].calBandHeader[i].nchan];
			  off1  = inFile->beam[b].bandData[i].cal_off_data.pol1[k+j*inFile->beam[b].calBandHeader[i].nchan];
			  on2   = inFile->beam[b].bandData[i].cal_on_data.pol2[k+j*inFile->beam[b].calBandHeader[i].nchan];
			  off2  = inFile->beam[b].bandData[i].cal_off_data.pol2[k+j*inFile->beam[b].calBandHeader[i].nchan];
			  on3   = inFile->beam[b].bandData[i].cal_on_data.pol3[k+j*inFile->beam[b].calBandHeader[i].nchan];
			  off3  = inFile->beam[b].bandData[i].cal_off_data.pol3[k+j*inFile->beam[b].calBandHeader[i].nchan];
			  on4   = inFile->beam[b].bandData[i].cal_on_data.pol4[k+j*inFile->beam[b].calBandHeader[i].nchan];
			  off4  = inFile->beam[b].bandData[i].cal_off_data.pol4[k+j*inFile->beam[b].calBandHeader[i].nchan];

			  tsys[j*2*inFile->beam[b].calBandHeader[i].nchan + k] = tcalA*off1/(on1-off1);
			  tsys[j*2*inFile->beam[b].calBandHeader[i].nchan + inFile->beam[b].calBandHeader[i].nchan + k] = tcalB*off2/(on2-off2);
			  
			  d1 = on1 - off1;
			  d2 = on2 - off2;
			  d3 = on3 - off3;
			  d4 = on4 - off4;
			  
			  s1 = d1+d2;
			  s2 = d1-d2;
			  s3 = 2*d3;
			  s4 = 2*d4;
			  
			  gain[j*inFile->beam[b].calBandHeader[i].nchan + k] = 2*s2/s1;
			  phase[j*inFile->beam[b].calBandHeader[i].nchan + k] = atan2(s4,s3)*180.0/M_PI;
			
			  if (i==0 && j==0)
			    printf("result: %d %d %g %g %g %g %g %g\n",j,k,tsys[j*2*inFile->beam[b].calBandHeader[i].nchan + k], gain[j*inFile->beam[b].calBandHeader[i].nchan + k],phase[j*inFile->beam[b].calBandHeader[i].nchan + k],tcalA,on1,off1);
			}
		    }


		  // Calculate averages
		  for (k=0;k<inFile->beam[b].calBandHeader[i].nchan;k++)
		    {		      
		      avTsys=0;
		      for (j=0;j<inFile->beam[b].calBandHeader[i].ndump;j++)
			{
			  avTsys += tsys[j*2*inFile->beam[b].calBandHeader[i].nchan + k] ;
			}
		      //		      if (i==5)
		      //			printf("result: %.5f %g\n",inFile->beam[b].bandData[i].cal_on_data.freq[k],avTsys/inFile->beam[b].calBandHeader[i].ndump);
		    }
		  
		  sdhdf_writeCalProc(outFile,b,i,inFile->beam[b].bandHeader[i].label,"cal_proc_tsys",tsys,inFile->beam[b].calBandHeader[i].nchan,2,inFile->beam[b].calBandHeader[i].ndump);
		  sdhdf_writeCalProc(outFile,b,i,inFile->beam[b].bandHeader[i].label,"cal_proc_diff_gain",gain,inFile->beam[b].calBandHeader[i].nchan,1,inFile->beam[b].calBandHeader[i].ndump);
		  sdhdf_writeCalProc(outFile,b,i,inFile->beam[b].bandHeader[i].label,"cal_proc_diff_phase",phase,inFile->beam[b].calBandHeader[i].nchan,1,inFile->beam[b].calBandHeader[i].ndump);
		  
		  
		  free(tsys);
		  free(gain);
		  free(phase);
		}
	    }
	  
	  sdhdf_closeFile(inFile);
	  sdhdf_closeFile(outFile);
	}
    }
  free(inFile);
  free(outFile);

  free(tcalData);
}
