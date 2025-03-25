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

void help()
{
}


int main(int argc,char *argv[])
{
  int i,j,k,p,ii,b;
  char fname[1024];
  char oname[1024] = "scaling.dat";
  char src[1024]="unset";
  int nFiles=0;
  sdhdf_fileStruct *inFile;
  sdhdf_bandHeaderStruct *inBandParams;  
  double freq;
  double pred;
  int nchan,npol;
  double sc1,sc2;
  FILE *fout;
  
  // help();

  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(oname,argv[++i]);
      else if (strcmp(argv[i],"-src")==0)
	strcpy(src,argv[++i]);
    }

  if (strcmp(src,"1934-638")==0)
    printf("Using primary flux calibrator: %s\n",src);
  else
    {
      printf("ERROR: Unknown primary flux calibrator. Use -src to set\n");
      exit(1);
    }
  sdhdf_initialiseFile(inFile);
  sdhdf_openFile(fname,inFile,1);
  sdhdf_loadMetaData(inFile);

  fout = fopen(oname,"w");
  
  for (b=0;b<inFile->nBeam;b++)
    {
      for (i=0;i<inFile->beam[b].nBand;i++)
	{	  
	  nchan = inFile->beam[b].bandHeader[i].nchan;
	  npol  = inFile->beam[b].bandHeader[i].npol;
	  sdhdf_loadBandData(inFile,b,i,1);

	  for (k=0;k<nchan;k++)
	    {
	      freq = inFile->beam[b].bandData[i].astro_data.freq[k]; // FIX ME -- FREQ AXIS FOR DUMP
	      pred = pow(10,-30.7667+26.4908*log10(freq) - 7.0977*pow(log10(freq),2) + 0.605334*pow(log10(freq),3));
	      sc1  =  inFile->beam[b].bandData[i].astro_data.pol1[k]/pred;
	      if (npol==2 || npol==4)
		sc2  =  inFile->beam[b].bandData[i].astro_data.pol2[k]/pred;

	      if (npol==1)
		fprintf(fout,"%.6f %g\n",freq,sc1);
	      else
		fprintf(fout,"%.6f %g %g\n",freq,sc1,sc2);
	    }
	  
	  sdhdf_releaseBandData(inFile,b,i,1);
	}
    }
  //  free(inBandParams);
  printf("Complete: scaling factors written to: %s\n",oname);
  
  sdhdf_closeFile(inFile);
    
  free(inFile);
  fclose(fout);


}

void determineOffSourceScale(float *freqOffScl,float *offSclA,float *offSclB,int n,float freq,float *sclA2,float *sclB2)
{
  int i;

  for (i=0;i<n-1;i++)
    {
      if (freq <= freqOffScl[i])
	{
	  *sclA2 = offSclA[i]; // Should do an interpolation
	  *sclB2 = offSclB[i];
	  break;
	}
    }
}
