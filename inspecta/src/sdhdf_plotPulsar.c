//  Copyright (C) 2019, 2020, 2021, 2022, 2023 George Hobbs

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
// Software to plot a pulsar signal

//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"
#include <cpgplot.h>

#define VNUM "v0.1"
#define MAX_HLINE 16
#define MAX_VLINE 16

void help()
{
}


int main(int argc,char *argv[])
{
  int i,j,k;
  char       fname[MAX_STRLEN];
  sdhdf_fileStruct *inFile;
  int ibeam=0;
  int iband=0;
  float *allData;
  float *freqPhase;
  float *timePhase;
  float *plotProfile1;
  float px[4096];
  long nbin,nchan,ndump,npol;
  float tr[6];
  float minV,maxV;
  float minProf,maxProf;
  float fref;
  int nFreqDump;
  int i0;
  float dt;
  float f0,df;
  float mean,mean1,mean2;
  float dm=306.9;                // FIX ME
  float period=0.88776601200758; // FIX ME
  float binOffset;
  float *freqVals;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  //  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.2, 1.0, 1.0};
  //  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  //  float heat_b[] = {1.0, 0.8, 0.6, 0.4, 1.0};
  float heat_b[] = {1.0, 0.8, 0.6, 0.4, 1.0};
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);
  

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }
  sdhdf_initialiseFile(inFile);
  if (sdhdf_openFile(fname,inFile,1)==-1)
    {
      printf("Unable to open input file >%s<\n",fname);
      free(inFile);
      exit(1);
    }
  sdhdf_loadMetaData(inFile);
  printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname,inFile->primary[0].utc0, inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[ibeam].source,inFile->primary[0].telescope, inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[ibeam].nBand);

  cpgbeg(0,"?",1,1);
  cpgsch(1.2);
  cpgslw(2);
  cpgscf(2);

  // Plot colour plot of frequeny versus phase
  nbin = 1024; // FIX ME ******   NOT YET IN HEADER PARAMETERS
  nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
  ndump = inFile->beam[ibeam].bandHeader[iband].ndump;
  npol = inFile->beam[ibeam].bandHeader[iband].npol;
  printf("nbin = %d, nchan = %d\n",nbin,nchan);
  
  allData = (float *)malloc(sizeof(float)*nbin*nchan*npol*ndump);
  freqPhase = (float *)malloc(sizeof(float)*nbin*nchan);
  freqVals = (float *)malloc(sizeof(float)*nchan);
  timePhase = (float *)malloc(sizeof(float)*nbin*ndump);
  plotProfile1 = (float *)malloc(sizeof(float)*nchan);
  sdhdf_loadBandData2Array(inFile,ibeam,iband,1,allData);
  sdhdf_loadFrequency2Array(inFile,ibeam,iband,freqVals,&nFreqDump);
  f0 = freqVals[0]; df = freqVals[1]-freqVals[0];

  // Remove a mean from each spectrum
  for (i=0;i<ndump;i++)
    {
      for (j=0;j<nchan;j++)
	{
	  mean1=mean2=0.0;
	  for (k=0;k<nbin;k++)
	    {
	      mean1 += allData[i*nbin*nchan*npol + j*nbin + k];
	      mean2 += allData[i*nbin*nchan*npol + 1*nchan*nbin + j*nbin + k];
	    }
	  mean1/=(double)nbin;
	  mean2/=(double)nbin;
	  for (k=0;k<nbin;k++)
	    {
	      allData[i*nbin*nchan*npol + j*nbin + k] -= mean1;
	      allData[i*nbin*nchan*npol + 1*nbin*nchan + j*nbin + k] -= mean2;	  
	    }
	}
    }



  minV = 1e30; maxV = -1e30;
  for (k=0;k<nchan;k++)
    {
      for (i=0;i<nbin;i++)
	{
	  freqPhase[k*nbin+i]=0;
	  for (j=0;j<ndump;j++)
	    freqPhase[k*nbin+i] += (allData[j*nbin*nchan*npol + k*nbin + i] +
				    allData[j*nbin*nchan*npol + nchan*nbin + k*nbin + i]);
	}
    }

  
  for (i=0;i<nchan;i++)
    {
      for (j=0;j<nbin;j++)
	{
	  if (minV  > freqPhase[i*nbin+j]) minV = freqPhase[i*nbin+j];
	  if (maxV  < freqPhase[i*nbin+j]) maxV = freqPhase[i*nbin+j];
	}
    }
  printf("min/max = %g/%g\n",minV,maxV);
  tr[0] = 0; tr[1] = 1; tr[2] = 0;
  tr[3] = f0; tr[4] = 0; tr[5] = df;
  
  cpgsvp(0.1,0.9,0.10,0.45);
  cpgswin(0,nbin,f0,f0+nchan*df);
  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
  cpglab("Pulse phase bin","Frequency (MHz)","");
  cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
  cpgimag(freqPhase,nbin,nchan,1,nbin,1,nchan,minV,maxV,tr);

  // Time versus phase
  minV = 1e30; maxV = -1e30;  
  fref = 0.5*(freqVals[0]+freqVals[nchan-1]);
  for (k=0;k<ndump;k++)
    {
      for (i=0;i<nbin;i++)
	{
	  timePhase[k*nbin+i]=0;
	  for (j=0;j<nchan;j++)
	    {
	      dt = -4.15e-3*dm*(pow(fref/1000.0,-2)-pow(freqVals[j]/1000.0,-2));
	      binOffset = nbin*dt/period;
	      i0 = (int)(i + binOffset); // FIX ME 
	      //	      printf("dt = %g, binOffset = %g, i = %d i0 = %d\n",dt,binOffset,i,i0);
	      if (i0 < 0) i0+=nbin; // DO THIS PROPERLY FIX ME
	      if (i0 > nbin) i0-=nbin;

	      timePhase[k*nbin+i]+=(allData[k*nbin*nchan*npol + j*nbin + i0] +
				    allData[k*nbin*nchan*npol + nbin*nchan + j*nbin + i0]);
	    }
	}
    }
  for (i=0;i<nbin;i++)
    {
      px[i] = i;
      plotProfile1[i] = 0;
    }
  minProf = 1e30; maxProf = -1e30;
  for (i=0;i<ndump;i++)
    {
      mean=0;
      for (j=0;j<nbin;j++)
	mean+=timePhase[i*nbin+j];
      mean/=(double)nbin;
      for (j=0;j<nbin;j++)
	{
	  timePhase[i*nbin+j]-=mean;
	  plotProfile1[j] += timePhase[i*nbin+j];
	  if (minProf > plotProfile1[j]) minProf = plotProfile1[j];
	  if (maxProf < plotProfile1[j]) maxProf = plotProfile1[j];
	  if (minV > timePhase[i*nbin+j]) minV = timePhase[i*nbin+j];
	  if (maxV < timePhase[i*nbin+j]) maxV = timePhase[i*nbin+j];
	}
    }
  printf("time/phase min/max  %g/%g\n",minV,maxV);
  cpgsvp(0.1,0.9,0.45,0.8);
  cpgswin(0,nbin,0,ndump);
  cpgbox("ABCTS",0,0,"ABCTSN",0,0);
  cpglab("","Subintegration","");
  tr[0] = 0; tr[1] = 1; tr[2] = 0;
  tr[3] = 0; tr[4] = 0; tr[5] = 1;

  cpgimag(timePhase,nbin,ndump,1,nbin,1,ndump,minV,maxV,tr);

  // Plot average pulse profile with polarisation
  cpgsvp(0.1,0.9,0.8,0.98);
  cpgswin(0,nbin,minProf,maxProf);
  cpgbox("ABCTS",0,0,"ABCTS",0,0);
  cpgline(nbin,px,plotProfile1);

  cpgend();

  free(freqPhase);
  free(freqVals);
  free(timePhase);
  free(allData);




  sdhdf_closeFile(inFile);



  

  free(inFile);
}


