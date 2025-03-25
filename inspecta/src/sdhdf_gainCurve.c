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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"
#include <cpgplot.h>

#define MAX_POINTS 4096
#define MAX_BANDS 26

float getMinMax2(float *in1,float *in2,int n,float *miny,float *maxy);
float getMinMax1(float *in1,int n,float *miny,float *maxy);

int main(int argc,char *argv[])
{
  int   i,j,k,b,l;
  float px1[MAX_BANDS][MAX_POINTS];
  float pyA[MAX_BANDS][MAX_POINTS];
  float pyB[MAX_BANDS][MAX_POINTS];
  int   npts[MAX_BANDS];
  float f0[MAX_BANDS],f1[MAX_BANDS];
  sdhdf_fileStruct *inFile;
  float *tsys,*avTsys_AA,*avTsys_BB,freq;
  int   nchan,ndump,npol,nband;
  char  fname[1024];
  int   ibeam=0;
  int   iband;
  float minx,maxx,miny,maxy;
  char title[1024];
  double mjd0=-1;
  int count;
  int xplot=2;
  
  // FIX ME: HARDCODED AND SET TO PARKES -- SHOULD USE FLAG VALUES
  f0[0] = 813;     f1[0] = 815;
  f0[1] = 938;     f1[1] = 942;
  f0[2] = 1007;     f1[2] = 1008;
  f0[3] = 1114;    f1[3] = 1119;
  f0[4] = 1305;    f1[4] = 1328;
  f0[5] = 1370;    f1[5] = 1410;
  f0[6] = 1510;    f1[6] = 1520;
  f0[7] = 1660;    f1[7] = 1680;
  f0[8] = 1775;    f1[8] = 1800;
  f0[9] = 1930;    f1[9] = 1970;
  f0[10] = 2010;   f1[10] = 2040;
  f0[11] = 2180;     f1[11] = 2210;
  f0[12] = 2270;     f1[12] = 2290;
  f0[13] = 2430;     f1[13] = 2445;
  f0[14] = 2520;     f1[14] = 2550;
  f0[15] = 2700;     f1[15] = 2720;
  f0[16] = 2780;     f1[16] = 2860;
  f0[17] = 2900;     f1[17] = 2990;
  f0[18] = 3030;     f1[18] = 3065;
  f0[19] = 3150;     f1[19] = 3250;
  f0[20] = 3290;     f1[20] = 3360;
  f0[21] = 3410;     f1[21] = 3440;
  f0[22] = 3580;     f1[22] = 3630;
  //f0[22] = 3605;     f1[22] = 3615;
  f0[23] = 3660;     f1[23] = 3750;
  f0[24] = 3850;     f1[24] = 3880;
  f0[25] = 3920;     f1[25] = 4000;
  
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  for (i=0;i<MAX_BANDS;i++)
    npts[i] = 0;
  
  
  for (i=1;i<argc;i++)
    {
      strcpy(fname,argv[i]);
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname,inFile,1);
      sdhdf_loadMetaData(inFile);

      nband = inFile->beam[ibeam].nBand;
      for (j=0;j<nband;j++)
	{
	  sdhdf_loadBandData(inFile,ibeam,j,2);
	  sdhdf_loadBandData(inFile,ibeam,j,3);

	  ndump = inFile->beam[ibeam].calBandHeader[j].ndump;
	  nchan = inFile->beam[ibeam].calBandHeader[j].nchan;
	  npol  = inFile->beam[ibeam].calBandHeader[j].npol;
	  tsys = (float *)malloc(sizeof(float)*ndump*nchan*2); // * 2 = AA and BB
	  avTsys_AA = (float *)calloc(sizeof(float),nchan);
	  avTsys_BB = (float *)calloc(sizeof(float),nchan);

	  sdhdf_loadCalProc(inFile,ibeam,j,"cal_proc_tsys",tsys);
	  for (l=0;l<nchan;l++)
	    {
	      avTsys_AA[l] = avTsys_BB[l] = 0;
	      for (k=0;k<ndump;k++)
		{
		  avTsys_AA[l] += tsys[l+k*2*nchan]/(double)ndump;
		  avTsys_BB[l] += tsys[l+k*2*nchan+nchan]/(double)ndump;
		  if (j==5 && k==ndump-1) printf("test: %g %g\n",tsys[l+k*2*nchan],avTsys_AA[l]);
		}
	    }
	  //	  px1[j][npts[j]] = i;

	  if (mjd0 < 0)
	    mjd0 = inFile->beam[ibeam].bandData[j].cal_obsHeader[0].mjd;

	  if (xplot==1)
	    px1[j][npts[j]] = (inFile->beam[ibeam].bandData[j].cal_obsHeader[0].mjd - mjd0)*86400.0/60.0/60.0; // FIX ME: Choose centre of dump 
	  else if (xplot==2)  // FIX ME: Choose centre of dump
	    px1[j][npts[j]] = inFile->beam[ibeam].bandData[j].cal_obsHeader[0].el;
	  else if (xplot==3)
	    px1[j][npts[j]] = inFile->beam[ibeam].bandData[j].cal_obsHeader[0].az;
	  else if (xplot==4)
	    px1[j][npts[j]] = inFile->beam[ibeam].bandData[j].cal_obsHeader[0].paraAngle;
	  else if (xplot==5)
	    px1[j][npts[j]] = inFile->beam[ibeam].bandData[j].cal_obsHeader[0].hourAngle;

	  pyA[j][npts[j]] = 0;
	  pyB[j][npts[j]] = 0;
	  count=0;
	  for (k=0;k<nchan;k++) // Should set properly
	    {
	      freq = inFile->beam[ibeam].bandData[j].cal_on_data.freq[k];  // FIX ME -- FREQ AXIS FOR DUMP
	      if (freq > f0[j] && freq <= f1[j])
		{
		  pyA[j][npts[j]]+=avTsys_AA[k];
		  pyB[j][npts[j]]+=avTsys_BB[k];
		  count++;
		}
	    }
	  pyA[j][npts[j]]/=(double)(count);
	  pyB[j][npts[j]]/=(double)(count);
	  if (j==5)
	    printf("tsys %f %f %f\n",px1[j][npts[j]],pyA[j][npts[j]],pyB[j][npts[j]]);
	  npts[j]++;
		  

	  
	  free(tsys);
	  free(avTsys_AA); free(avTsys_BB);
	}
      
      sdhdf_closeFile(inFile);
    }

  free(inFile);

  cpgbeg(0,"gainCurve.ps/cps",3,3);
  cpgsch(1.4);
  cpgslw(2);
  cpgscf(2);
  cpgask(0);
  iband=0;

  for (iband=0;iband<nband;iband++)
    {
      getMinMax2(pyA[iband],pyB[iband],npts[iband],&miny,&maxy);
      getMinMax1(px1[iband],npts[iband],&minx,&maxx);
      maxy = miny + 3;
      cpgenv(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,1);
      sprintf(title,"band: %d",iband);
      if (xplot==1)
	cpglab("Time (hours)","Total temperature (K)",title);
      else if (xplot==2)
	cpglab("Elevation (degrees)","Total temperature (K)",title);
      else if (xplot==3)
	cpglab("Azimuth (degrees)","Total temperature (K)",title);
      else if (xplot==4)
	cpglab("Parallactic angle (degrees)","Total temperature (K)",title);
      else if (xplot==5)
	cpglab("Hour angle (degrees)","Total temperature (K)",title);
      cpgline(npts[iband],px1[iband],pyA[iband]);
      cpgpt(npts[iband],px1[iband],pyA[iband],2);
      cpgsci(5);
      cpgline(npts[iband],px1[iband],pyB[iband]);
      cpgpt(npts[iband],px1[iband],pyB[iband],2);
      cpgsci(1);
    }

  
  cpgend();
  

}

float getMinMax2(float *in1,float *in2,int n,float *miny,float *maxy)
{
  int i;
  *miny = *maxy = in1[0];
  for (i=0;i<n;i++)
    {
      if (*miny > in1[i]) *miny = in1[i];
      if (*miny > in2[i]) *miny = in2[i];
      if (*maxy < in1[i]) *maxy = in1[i];
      if (*maxy < in2[i]) *maxy = in2[i];
    }
  return 0;
}

float getMinMax1(float *in1,int n,float *miny,float *maxy)
{
  int i;
  *miny = *maxy = in1[0];
  for (i=0;i<n;i++)
    {
      if (*miny > in1[i]) *miny = in1[i];
      if (*maxy < in1[i]) *maxy = in1[i];
    }
  return 0;
}
