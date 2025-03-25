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


void oneBitFloat(int eight_bit_number, float *results, int *index);
void twoBitsFloat(int eight_bit_number, float *results, int *index);
void twoBitsFloat_unsigned(int eight_bit_number, float *results, int *index);
void fourBitsFloat(int eight_bit_number, float *results, int *index);
void eightBitsFloat(int eight_bit_number, float *results, int *index);
void pfits_bytesToFloats(int samplesperbyte,int n,unsigned char *cVals,float *out);


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
  float *plotData;
  unsigned char *loadData;
  float px[4096];
  long nsamples,nchan,npol,ndump,nsamplesPlot;
  int samplesperbyte;
  float tr[6];
  float minV,maxV;
  float minProf,maxProf;
  float fref;
  int nFreqDump;
  int i0;
  float dt;
  float f0,df;
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


  // Plot colour plot of frequeny versus phase
  nsamples = 4096*2051; // FIX ME ******   NOT YET IN HEADER PARAMETERS
  nsamplesPlot = 4096;
  samplesperbyte = 8;
  nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
  npol = inFile->beam[ibeam].bandHeader[iband].npol;
  ndump= inFile->beam[ibeam].bandHeader[iband].ndump;
  printf("nsamples = %d, nchan = %d ndump = %d\n",nsamples,nchan,ndump);
  loadData = (unsigned char *)malloc(sizeof(unsigned)*nsamples*nchan*npol/samplesperbyte);
  plotData = (float *)malloc(sizeof(float)*nsamplesPlot*nchan*npol);
  freqVals = (float *)malloc(sizeof(float)*nchan);
  printf("Loading freq\n");
  sdhdf_loadFrequency2Array(inFile,ibeam,iband,freqVals,&nFreqDump);
  f0 = freqVals[0]; df = freqVals[1]-freqVals[0];
  printf("Loading quantised\n");
  sdhdf_loadQuantisedBandData2Array(inFile,ibeam,iband,1,loadData);
  printf("Extracting bits\n");
  // Extract the bits from the bytes
  pfits_bytesToFloats(8,nchan*nsamplesPlot,loadData+(int)((long)421*(long)4096.0*nchan/samplesperbyte),plotData);
  //  for (i=0;i<nchan*nsamples;i++)
  //    printf("%f\n",plotData[i]);
  tr[0] = 412*4096*255e-6; tr[1] = 0; tr[2] = 255e-6;
  tr[3] = f0; tr[4] = df; tr[5] = 0;
  printf("Making the plot %d %g %g \n",nsamples,f0,f0+nchan*df);

  cpgbeg(0,"?",1,1);
  cpgsch(1.2);
  cpgslw(2);
  cpgscf(2);
  printf("Setting up cpgenv\n");
  cpgenv(412*4096*255e-6,412*4096*255e-6+255e-6*nsamplesPlot,f0,f0+nchan*df,0,1);
  cpglab("Time since start (seconds)","Frequency (MHz)","");
  cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
  cpgimag(plotData,nchan,nsamplesPlot,1,nchan,1,nsamplesPlot,0,1,tr);


  cpgend();

  free(loadData); free(plotData);
  free(freqVals);


  sdhdf_closeFile(inFile);



  

  free(inFile);
}


void pfits_bytesToFloats(int samplesperbyte,int n,unsigned char *cVals,float *out)
{
  int i,j;
  int pos=0;
  int index = 0;
  for (i = 0; i < n/samplesperbyte; i++)
    {
      void (*pointerFloat)(int, float *, int*);
      switch (samplesperbyte)
	{
	case 1:
	  pointerFloat = eightBitsFloat;
	  break;

	case 2:
	  pointerFloat = fourBitsFloat;
	  break;

	case 4:
	  //	  pointerFloat = twoBitsFloat;
	  pointerFloat = twoBitsFloat_unsigned;
	  break;

	case 8:
	  pointerFloat = oneBitFloat;
	  break;
	}
      pointerFloat(cVals[i], out, &index);
    }
}


void eightBitsFloat(int eight_bit_number, float *results, int *index)
{
  // 135 subtracted from the original number to make the range from -8 to something
  // 0.5 is then added to remove bias

  //results[(*index)++] = eight_bit_number - 135.5;
  results[(*index)++] = eight_bit_number; // -31.5; // -31.5;
}

void fourBitsFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant 4 bits with 0000 1111 will give the (signed) 4-bit number
  // shifting right 4 bits will produce the other (first) 4-bit numbers
  // 0.5 is added to each number to compensate for the bias, making the range -7.5 -> +7.5

  float tempResults[2];
  int i;
  for (i = 0; i < 2; i++) {
    int andedNumber = eight_bit_number & 15;
    tempResults[i] = andedNumber; // - 7.5;
    eight_bit_number = eight_bit_number >> 4;
  }

  for (i = 1; i >= 0; i--)
    results[(*index)++] = tempResults[i];
}

void twoBitsFloat_unsigned(int eight_bit_number, float *results, int *index)
{
  unsigned char uctmp = (unsigned char)eight_bit_number;
  int i;

  results[(*index)++] = ((uctmp >> 0x06) & 0x03);
  results[(*index)++] = ((uctmp >> 0x04) & 0x03);
  results[(*index)++] = ((uctmp >> 0x02) & 0x03);
  results[(*index)++] = ((uctmp & 0x03));
}

// Note DSPSR is using this technique, but then doesn't scale by ZERO_OFFS
void twoBitsFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant 2 bits with 0000 0011 will give the (signed 2-bit number
  // shifting right 2 bits will produce the next (previous) 2-bit number
  // the numbers are adjusted to compensate for bias and to give a rms of ~0.9
  //
  // 1 -> 2.0
  // 0 -> +0.5
  // -1 -> -0.5
  // -2 -> -2.0


  // NOTE: ACTUAL -2.5 below !!
  float tempResults[4];
  int i;

  for (i = 0; i < 4; i++) {
    int andedNumber = eight_bit_number & 3;

    switch (andedNumber) {
    case 0:
      //      tempResults[i] = 0; //-2.5;
      tempResults[i] = -2.5;
      break;
    case 1:
      //      tempResults[i] = 1; // -0.5;
      tempResults[i] = -0.5;      
      break;
    case 2:
      tempResults[i] = 0.5;
      //      tempResults[i] = 2; //0.5;
      break;
    case 3:
      tempResults[i] = 2.5;
      //      tempResults[i] = 3; // 2.5;
      break;
    }
    eight_bit_number = eight_bit_number >> 2;
  }

  for (i = 3; i >= 0; i--)
    {
      //      printf("Have %f\n",tempResults[i]);
      results[(*index)++] = tempResults[i];
    }
}

void oneBitFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant bit with 0000 0001 will give the (signed 1-bit number
  // shifting right 1 bit will produce the next (previous) 1-bit number
  //
  // 0.5 is added to each number to compensate for bias
  //

  float tempResults[8];
  int i;

  for (i = 0; i < 8; i++) {
    int andedNumber = eight_bit_number & 1;
    // tempResults[i] = andedNumber ? 0.5 : -0.5;
    tempResults[i] = andedNumber ? 0 : 1;
    eight_bit_number = eight_bit_number >> 1;
  }

  for (i = 7; i >= 0; i--)
    {
      results[(*index)++] = tempResults[i];
      //            printf("This pt: %d %g\n",*index,tempResults[i]);
    }
}
