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
// Software that carries out a quicklook at spectral line data
//
// Usage:
// sdhdf_plotSpectrum -f <filename>
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"
#include <cpgplot.h>

#define VNUM "v0.1"
#define MAX_CHANS 262144
#define MAX_DUMPS 1024
#define MAX_TRANSMITTERS 1024

void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY);
void showTransmitter(float freq,float bw,char *label,float miny,float maxy);
void plotSkyPositions(sdhdf_fileStruct *inFile,int ibeam);
void plotHI(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2);
void plotOH(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2);
void plotCH_SB0(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2);
void plotCH_SB19(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2);
void plotCH_SB20(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2);
void createPolSummedSpectrum(int band,int nPlot,float *px,float *py1,float *py2,float *cleanF1,
			     float *cleanF2,float *miny_sumPol,float *maxy_sumPol,float *sumPolX,float *sumPolY,
			     int *sumPolW,int *sumPol_nPlot);
void loadTransmitters(float *t1,float *t2,int *nTransmitters);
void setTransmitter(float freq,float bw,char *name,float *t1,float *t2);

void help()
{
  printf("sdhdf_quicklook: routine to have a quick look at a spectrum\n");  
  printf("-g                   Select graphics device\n");
  printf("-h                   This help");
  printf("-transmitters        Overplot transmitters\n\n");
  printf("\n\n");
  printf("Filenames are given on the command line\n");
  
}

int main(int argc,char *argv[])
{
  int        i,j,k,l,ii,kk;
  char       fname[MAX_FILES][64];
  sdhdf_fileStruct *inFile;
  int        idump,iband,ibeam;
  float f0=-1,f1=-1;
  int av=0;
  int sump=0;
  int nx=1;
  int ny=1;
  int polPlot=1;
  long nVals=0;
  float *px,*py1,*py2;
  float *highResX,*highResY1,*highResY2;
  float *sumPolX,*sumPolY;
  int *sumPolW;
  float *dump_x,*dump_y1,*dump_y2;
  long writepos;
  int allocateMemory=0;
  float miny,maxy,val1,val2,minx,maxx;
  float miny_sumPol,maxy_sumPol;
  float hmin,hmax;
  float miny1_dump,maxy1_dump;
  float miny2_dump,maxy2_dump;
  float ominy,omaxy,ominx,omaxx;
  float plotTransmitters=-1;
  int molecularLines=-1;
  int inFiles;
  int sb=-1;
  float mx,my;
  char key;
  char grDev[128]="";
  float xp,yp;
  int nav=512;
  int nPlot;
  int nSum;
  int nDumpAdd;
  float cleanF1[26],cleanF2[26];
  char text[1024];
  int sumPol_nPlot;
  
  px = (float *)malloc(sizeof(float)*MAX_CHANS);
  py1 = (float *)malloc(sizeof(float)*MAX_CHANS);
  py2 = (float *)malloc(sizeof(float)*MAX_CHANS);
  sumPolX = (float *)malloc(sizeof(float)*MAX_CHANS);
  sumPolY = (float *)malloc(sizeof(float)*MAX_CHANS);
  sumPolW = (int *)malloc(sizeof(int)*MAX_CHANS);


  
  highResX = (float *)malloc(sizeof(float)*MAX_CHANS);
  highResY1 = (float *)calloc(sizeof(float),MAX_CHANS);
  highResY2 = (float *)calloc(sizeof(float),MAX_CHANS);

  dump_x = (float *)malloc(sizeof(float)*MAX_DUMPS);
  dump_y1 = (float *)malloc(sizeof(float)*MAX_DUMPS);
  dump_y2 = (float *)malloc(sizeof(float)*MAX_DUMPS);

  cleanF1[0] = 740; cleanF2[0] = 750;
  cleanF1[1] = 857; cleanF2[1] = 863;
  cleanF1[2] = 985; cleanF2[2] = 991;
  cleanF1[3] = 1105; cleanF2[3] = 1115;
  cleanF1[4] = 1300; cleanF2[4] = 1308;
  cleanF1[5] = 1390; cleanF2[5] = 1403;
  cleanF1[6] = 1510; cleanF2[6] = 1520;
  cleanF1[7] = 1640; cleanF2[7] = 1660;
  cleanF1[8] = 1766; cleanF2[8] = 1784;
  cleanF1[9] = 1930; cleanF2[9] = 1960;
  cleanF1[10] = 2020; cleanF2[10] = 2040;
  cleanF1[11] = 2128; cleanF2[11] = 2136;
  cleanF1[12] = 2265; cleanF2[12] = 2285;
  cleanF1[13] = 2470; cleanF2[13] = 2477;
  cleanF1[14] = 2530; cleanF2[14] = 2550;
  cleanF1[15] = 2700; cleanF2[15] = 2730;
  cleanF1[16] = 2780; cleanF2[16] = 2810;
  cleanF1[17] = 2940; cleanF2[17] = 2960;
  cleanF1[18] = 3040; cleanF2[18] = 3060;
  cleanF1[19] = 1380; cleanF2[19] = 3220;
  cleanF1[20] = 3290; cleanF2[20] = 3320;
  cleanF1[21] = 3410; cleanF2[21] = 3435;
  cleanF1[22] = 3590; cleanF2[22] = 3600;
  cleanF1[23] = 3710; cleanF2[23] = 3740;
  cleanF1[24] = 3810; cleanF2[24] = 3830;
  cleanF1[25] = 3935; cleanF2[25] = 3970;


  
  printf("Starting\n");
  // Defaults
  idump = iband = ibeam = 0;
  strcpy(grDev,"/xs");
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);

  inFiles=0;
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-transmitters")==0)
	plotTransmitters=1;
      else if (strcmp(argv[i],"-h")==0)
	{help();exit(1);}
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else
	strcpy(fname[inFiles++],argv[i]);
    }
  
  for (i=0;i<inFiles;i++)
    {
      sprintf(grDev,"%s.quicklook.ps/cps",fname[i]);
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);

      sdhdf_loadMetaData(inFile);

      cpgbeg(0,grDev,1,1);
      cpgsch(1.0);
      cpgslw(2);
      cpgscf(2);
      cpgask(0);

      cpgsvp(0.15,0.95,0.15,0.95);
      cpgswin(0,1,0,1);
      cpgtext(0.1,0.97,fname[i]);
      sprintf(text,"Source: %s",inFile->beamHeader[0].source);
      cpgtext(0.1,0.93,text);
      sprintf(text,"Observer: %s",inFile->primary[0].observer);
      cpgtext(0.1,0.89,text);
      printf("Plotting sky\n");
      plotSkyPositions(inFile,ibeam);
      
      cpgpage();
      
      printf("Loading %s\n",fname[i]);
      
      
      nVals=0;
      for (j=0;j<inFile->beam[ibeam].nBand;j++)
	{
	  sdhdf_loadBandData(inFile,ibeam,j,1);
	  nVals+=(inFile->beam[ibeam].bandHeader[j].nchan);
	}
      

      for (j=0;j<inFile->beam[ibeam].nBand;j++)
	{
	  printf("Processing band: %d\n",j);
	  miny1_dump = 1e99; maxy1_dump = -1e99;
	  miny2_dump = 1e99; maxy2_dump = -1e99;
	  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
	    {
	      highResY1[k] = 0;
	      highResY2[k] = 0;
	    }
	  
	  for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l++) 
	    {	      
	      miny = 1e99; maxy = -1e99;
	      nPlot=0;

	      nDumpAdd=0;
	      dump_y1[l]=0.0;
	      dump_y2[l]=0.0;
	      
	      
	      for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		{
		  // FIX ME: using [0] for frequency dump
		  if (l==0) {highResX[k] = inFile->beam[ibeam].bandData[j].astro_data.freq[k];}
		  val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
		  val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
		       
		  highResY1[k] += val1;
		  highResY2[k] += val2;
		  
		  if (highResX[k] > cleanF1[j] && highResX[k] < cleanF2[j])
		    {
		      dump_y1[l] += val1;
		      dump_y2[l] += val2;
		      nDumpAdd++;
		    }
		}
	      dump_x[l] = l;		  
	      dump_y1[l]/=(float)nDumpAdd;
	      dump_y2[l]/=(float)nDumpAdd;
	      if (dump_y1[l] > maxy1_dump) maxy1_dump = dump_y1[l];
	      if (dump_y2[l] > maxy2_dump) maxy2_dump = dump_y2[l];
	      if (dump_y1[l] < miny1_dump) miny1_dump = dump_y1[l];
	      if (dump_y2[l] < miny2_dump) miny2_dump = dump_y2[l];		

	     	      
	
	    }
	      /*  
	      */
	  //	    }
    
	  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
	    {
	      highResY1[k]/=(float)inFile->beam[ibeam].bandHeader[j].ndump;
	      highResY2[k]/=(float)inFile->beam[ibeam].bandHeader[j].ndump;
	      //	      printf("vals : %g %g\n",highResY1[k],highResY2[k]);
	    }
	  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k+=nav)
	    {
	      py1[nPlot] = py2[nPlot] = 0;
	      px[nPlot] = 0;
	      nSum=0;
	      for (kk=k;kk<k+nav;kk++)
		{
		  if (kk < inFile->beam[ibeam].bandHeader[j].nchan)
		    {
		      val1 = highResY1[kk];
		      val2 = highResY2[kk];
		      highResX[kk] = highResX[kk];
		      highResY1[kk] += val1;
		      highResY2[kk] += val2;
		      py1[nPlot] += val1;
		      py2[nPlot] += val2;
		      // FIX ME: using [0] for frequency dump
		      px[nPlot] += inFile->beam[ibeam].bandData[j].astro_data.freq[kk];
		      nSum++;
		    }
		}
	      
	      py1[nPlot]/=(float)nSum;
	      py2[nPlot]/=(float)nSum;
	      px[nPlot]/=(float)nSum;
	      if (py1[nPlot] > maxy) maxy = py1[nPlot];
	      if (py2[nPlot] > maxy) maxy = py2[nPlot];
	      if (py1[nPlot] < miny) miny = py1[nPlot];
	      if (py2[nPlot] < miny) miny = py2[nPlot];		
	      nPlot++;
	      
	    }
	  
	
	  
	  // Plot the linear scaling of the spectrum
      /*
	cpgsvp(0.1,0.95,0.15,0.45);
	cpgswin(704+j*128,704+(j+1)*128,miny,maxy);
	cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
	cpglab("Frequency (MHz)","Spectrum (arb.)","");
	cpgline(nPlot,px,py1);
	cpgsci(2); cpgline(nPlot,px,py2);
	cpgsci(1);
      */
	  
      // Plot the logarithmic scaling
	  miny = 1e99; maxy = -1e99;
	  for (k=0;k<nPlot;k++)
	    {
	      py1[k] = log10(py1[k]);
	      py2[k] = log10(py2[k]);
	      if (py1[k] > maxy) maxy = py1[k];
	      if (py2[k] > maxy) maxy = py2[k];
	      if (py1[k] < miny) miny = py1[k];
	      if (py2[k] < miny) miny = py2[k];		
	    }
	  //	  cpgsvp(0.1,0.95,0.15,0.45);
	  cpgsvp(0.1,0.95,0.15,0.30);
	  cpgswin(704+j*128,704+(j+1)*128,miny,maxy);
	  cpglab("Frequency (MHz)","           Spectrum (arb.)","");
	  cpgbox("ABCTSN",0,0,"ABCTSNL",0,0);
	  cpgsci(1); cpgline(nPlot,px,py1);
	  cpgsci(2); cpgline(nPlot,px,py2); cpgsci(1);

	  // Plot pol-summed spectrum
	
	  createPolSummedSpectrum(j,nPlot,px,py1,py2,cleanF1,cleanF2,&miny_sumPol,&maxy_sumPol,sumPolX,sumPolY,sumPolW,&sumPol_nPlot);
	  cpgsvp(0.1,0.95,0.3,0.45);
	  cpgswin(704+j*128,704+(j+1)*128,miny_sumPol,maxy_sumPol);
	  cpglab("","","");
	  cpgbox("ABCTS",0,0,"ABCTSN",0,0);
	  cpgline(sumPol_nPlot,sumPolX,sumPolY);
	  for (k=0;k<sumPol_nPlot;k++)
	    {
	      if (sumPolW[k] == 0)
		{cpgsci(8); cpgpt(1,sumPolX+k,sumPolY+k,16); cpgsci(1);}
	    }
		 
	  // Make plot of the power as a function of dump in clean parts of the band
	  cpgsvp(0.1,0.45,0.55,0.80);
	  cpgswin(0,inFile->beam[ibeam].bandHeader[j].ndump,miny1_dump,maxy1_dump);
	  cpgsci(4);
	  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
	  cpgsch(0.8); cpglab("Spectral dump","Signal level",""); cpgsch(1.4);
	  cpgsci(1); cpgline(inFile->beam[ibeam].bandHeader[j].ndump,dump_x,dump_y1);
	  cpgswin(0,inFile->beam[ibeam].bandHeader[j].ndump,miny2_dump,maxy2_dump); // NOTE AXIS LABELS NOT CHANGING HERE
	  cpgsci(2); cpgline(inFile->beam[ibeam].bandHeader[j].ndump,dump_x,dump_y2); cpgsci(1);
	  cpgsci(1);

	  // Make a spectrum of interesting bit(s)
	  if (j==0)       plotCH_SB0(inFile,j,ibeam,highResX,highResY1,highResY2);
	  else if (j==5)  plotHI(inFile,j,ibeam,highResX,highResY1,highResY2);
	  else if (j==7)  plotOH(inFile,j,ibeam,highResX,highResY1,highResY2);
	  else if (j==19) plotCH_SB19(inFile,j,ibeam,highResX,highResY1,highResY2);
	  else if (j==20) plotCH_SB20(inFile,j,ibeam,highResX,highResY1,highResY2);
	  
	  // Provide some statistics
	  cpgsvp(0.1,0.45,0.83,0.98);
	  cpgswin(0,1,0,1);	 
	  cpgbox("ABC",0,0,"ABC",0,0);
	  cpgsch(0.7);
	  cpgtext(0.07,0.88,inFile->fname);
	  sprintf(text,"Baseline band: %.0f to %0.f MHz",cleanF1[j],cleanF2[j]);
	  cpgtext(0.07,0.76,inFile->beamHeader[0].source);
	  sprintf(text,"Baseline band: %.0f to %0.f MHz",cleanF1[j],cleanF2[j]);
	  cpgtext(0.07,0.64,text);
	  sprintf(text,"Spectral dumps: %d",inFile->beam[ibeam].bandHeader[j].ndump);
	  cpgtext(0.07,0.52,text);
	  cpgsch(1.0);
	    
	  cpgpage();
	}

      


      cpgend();  
 
      sdhdf_closeFile(inFile);


    }

  
  free(dump_x); free(dump_y1); free(dump_y2);
  free(highResX); free(highResY1); free(highResY2);
  free(px);
  free(py1);
  free(py2);
  free(inFile);
  free(sumPolX); free(sumPolY);
  free(sumPolW);
}

void plotHI(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2)
{
  float hmin,hmax;
  int i,k;
  float freq0,freq1;
  int k0=-1,k1=-1;
  
  freq0 = 1419;
  freq1 = 1422;
  
  cpgsvp(0.5,0.95,0.55,0.98);
  
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }
  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgsch(0.8);
  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
  cpglab("Frequency (MHz)","","");
  cpgsch(1.4);
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);
}

void plotOH(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2)
{
  float hmin,hmax;
  int i,k;
  float freq0,freq1;
  int k0=-1,k1=-1;

  cpgsch(0.8);
  
  freq0 = 1612;
  freq1 = 1613;
  k0=k1=-1;
  cpgsvp(0.48,0.7,0.55,0.75);
  
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpglab("Frequency (MHz)","","");
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);


  // Second transition
  k0=k1=-1;
  freq0 = 1664.8;
  freq1 = 1666;
  cpgsvp(0.73,0.95,0.55,0.75);
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpglab("Frequency (MHz)","","");
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);

  // Third transition
  k0=k1=-1;
  freq0 = 1667.15;
  freq1 = 1667.4;
  cpgsvp(0.48,0.7,0.8,0.99);
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);

  // Fourth transition
  k0=k1=-1;
  freq0 = 1720.2;
  freq1 = 1720.8;
  cpgsvp(0.73,0.95,0.8,0.99);
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);
  cpgsch(1);
  
}

void plotCH_SB0(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2)
{
  float hmin,hmax;
  int i,k;
  float freq0,freq1;
  int k0=-1,k1=-1;

  cpgsch(0.8);
  
  freq0 = 704;
  freq1 = 704.3;
  k0=k1=-1;
  cpgsvp(0.48,0.7,0.55,0.98);
  
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpglab("Frequency (MHz)","","");
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);


  // Second transition
  k0=k1=-1;
  freq0 = 722;
  freq1 = 722.8;
  cpgsvp(0.73,0.95,0.55,0.98);
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpglab("Frequency (MHz)","","");
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);

  
}

void plotCH_SB19(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2)
{
  float hmin,hmax;
  int i,k;
  float freq0,freq1;
  int k0=-1,k1=-1;

  cpgsch(0.8);
  
  freq0 = 3262;
  freq1 = 3264;
  k0=k1=-1;
  cpgsvp(0.48,0.7,0.55,0.98);
  
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpglab("Frequency (MHz)","","");
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);



}

void plotCH_SB20(sdhdf_fileStruct *inFile,int j, int ibeam, float *highResX,float *highResY1,float *highResY2)
{
  float hmin,hmax;
  int i,k;
  float freq0,freq1;
  int k0=-1,k1=-1;

  cpgsch(0.8);
  
  freq0 = 3334;
  freq1 = 3336;
  k0=k1=-1;
  cpgsvp(0.48,0.7,0.55,0.98);
  
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpglab("Frequency (MHz)","","");
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);


  // Second transition
  k0=k1=-1;
  freq0 = 3348;
  freq1 = 3350;
  cpgsvp(0.73,0.95,0.55,0.98);
  hmin = 1e99; hmax=-1e99;
  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
    {
      if (highResX[k] > freq0 && highResX[k] < freq1)
	{
	  if (k0 == -1) k0 = k;
	  if (hmin > highResY1[k]) hmin = highResY1[k];
	  if (hmin > highResY2[k]) hmin = highResY2[k];
	  if (hmax < highResY1[k]) hmax = highResY1[k];
	  if (hmax < highResY2[k]) hmax = highResY2[k];
	}
      else if (k0 != -1 && k1 == -1)
	k1 = k;
      
    }

  cpgswin(freq0,freq1,hmin,hmax);
  cpgsci(1);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpglab("Frequency (MHz)","","");
  cpgsci(1);
  
  cpgline(k1-k0,highResX+k0,highResY1+k0);
  cpgsci(2); cpgline(k1-k0,highResX+k0,highResY2+k0); cpgsci(1);

}



void showTransmitter(float freq,float bw,char *label,float miny,float maxy)
{
  float fx[2],fy[2];
  cpgptxt(freq,maxy-(maxy-miny)*0.05,0,0.5,label);
  cpgsls(4);
  fx[0] = fx[1] = freq-bw/2;
  fy[0] = miny;  fy[1] = maxy;
  cpgline(2,fx,fy);
  fx[0] = fx[1] = freq+bw/2;
  fy[0] = miny;  fy[1] = maxy;
  cpgline(2,fx,fy);
  cpgsls(1);
}

      

void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY)
{
  float fx[2],fy[2];

  fx[0] = fx[1] = freq;
  fy[0] = minY;
  fy[1] = maxY;
  cpgsci(3); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);

}

void plotSkyPositions(sdhdf_fileStruct *inFile,int ibeam)
{
  float *raDeg,*decDeg;
  int iband=0;
  int i;
  int ndump = inFile->beam[ibeam].bandHeader[iband].ndump;
  float minx,maxx,miny,maxy;
  float midRA,midDEC;
  float beamX[1024],beamY[1024];
  float angle;
  int nbeam;
  float x0,y0;
  float fwhm;
  
  minx = 1e99; maxx = -1e99;
  miny = 1e99; maxy = -1e99;
  
  printf("Setting ndump\n");
  raDeg  = (float *)malloc(sizeof(float)*ndump);
  decDeg = (float *)malloc(sizeof(float)*ndump);
  printf("IN HERE\n");
  for (i=0;i<ndump;i++)
    {
      raDeg[i]  = inFile->beam[ibeam].bandData[iband].astro_obsHeader[i].raDeg;
      decDeg[i] = inFile->beam[ibeam].bandData[iband].astro_obsHeader[i].decDeg;
      if (minx > raDeg[i])  minx = raDeg[i];
      if (maxx < raDeg[i])  maxx = raDeg[i];
      if (miny > decDeg[i]) miny = decDeg[i];
      if (maxy < decDeg[i]) maxy = decDeg[i];
      
      printf("Loading %g %g\n",raDeg[i],decDeg[i]);
    }
  midRA  = (minx+maxx)/2.0;
  midDEC = (miny+maxy)/2.0;

  for (i=0;i<ndump;i++)
    {
      raDeg[i]  = (raDeg[i]-midRA)*60;
      decDeg[i] = (decDeg[i]-midDEC)*60;
    }  

  minx = (minx-midRA)*60;
  maxx = (maxx-midRA)*60;
  miny = (miny-midDEC)*60;
  maxy = (maxy-midDEC)*60;\
  // Add a beam at the high end
  miny -= 8.0;
  maxy += 8.0;
  minx -= 8.0;
  maxx += 8.0;

  if (maxx > maxy) maxy = maxx;
  else maxx = maxy;
  if (minx < miny) miny = minx;
  else minx = miny;

  
  
  //  if (miny == maxy)
  //    {
  //      miny-=0.01;
  //      maxy+=0.01;
  //    }
  printf("min/max = %g %g %g %g\n",minx,maxx,miny,maxy);
  cpgsvp(0.2,0.8,0.2,0.8);
  cpgswin(minx,maxx,miny,maxy);
  cpgbox("BCTSN",0,0,"BCTSN",0,0);  
  cpglab("Offset in right ascension (min)","Offset in declination (min)","");
  cpgpt(ndump,raDeg,decDeg,17);

  // Plot a beam
  nbeam=1024;
  x0 = minx+(maxx-minx)*0.5; y0 = miny-(maxy-miny)*0.5; 
  fwhm = 0.11; // Degrees
  for (i=0;i<nbeam;i++)
    {
      angle = 360*(float)i/(float)nbeam;
      beamX[i] = x0 + fwhm/2.0*60.*cos(angle*M_PI/180.0);
      beamY[i] = x0 + fwhm/2.0*60.*sin(angle*M_PI/180.0);
    }
  cpgsci(2); cpgline(nbeam,beamX,beamY); cpgsci(1);

  x0 = minx+(maxx-minx)*0.5; y0 = miny-(maxy-miny)*0.5; 
  fwhm = 0.25; // Degrees
  for (i=0;i<nbeam;i++)
    {
      angle = 360*(float)i/(float)nbeam;
      beamX[i] = x0 + fwhm/2.0*60.*cos(angle*M_PI/180.0);
      beamY[i] = x0 + fwhm/2.0*60.*sin(angle*M_PI/180.0);
    }
  cpgsci(4); cpgline(nbeam,beamX,beamY); cpgsci(1);

  x0 = minx+(maxx-minx)*0.5; y0 = miny-(maxy-miny)*0.5; 
  fwhm = 0.43; // Degrees
  for (i=0;i<nbeam;i++)
    {
      angle = 360*(float)i/(float)nbeam;
      beamX[i] = x0 + fwhm/2.0*60.*cos(angle*M_PI/180.0);
      beamY[i] = x0 + fwhm/2.0*60.*sin(angle*M_PI/180.0);
    }
  cpgsci(5); cpgline(nbeam,beamX,beamY); cpgsci(1);

  free(raDeg);
  free(decDeg);
}

void createPolSummedSpectrum(int band,int nPlot,float *px,float *py1,float *py2,float *cleanF1,
			     float *cleanF2,float *miny_sumPol,float *maxy_sumPol,float *sumPolX,
			     float *sumPolY,int *sumPolW, int *sumPol_nPlot)
{
  int i,j;
  int np=0;
  float t1[MAX_TRANSMITTERS],t2[MAX_TRANSMITTERS];
  int nTransmitters=0;
  int transmitter=0;

  
  loadTransmitters(t1,t2,&nTransmitters);
  
  
  *miny_sumPol = 1e99;
  *maxy_sumPol = -1e99;
  //  printf("Loaded %d transmitters\n",nTransmitters);
  for (i=0;i<nPlot;i++)
    {
      transmitter=0;
      for (j=0;j<nTransmitters;j++)
	{
	  if (px[i] >= t1[j] && px[i] <= t2[j])
	    {transmitter=1; break;}
	}
      
      sumPolX[np] = px[i];
      sumPolY[np] = py1[i]+py2[i];

      if (transmitter == 0)
	{
	  sumPolW[np] = 1;
	  if (*miny_sumPol > sumPolY[np]) *miny_sumPol = sumPolY[np];
	  if (*maxy_sumPol < sumPolY[np]) *maxy_sumPol = sumPolY[np];	    
	}
      else
	{
	  sumPolW[np] = 0;
	}
      np++;
    }
  *sumPol_nPlot = np;
  
}


// Should use a more generic function for this .....
//
void loadTransmitters(float *t1,float *t2,int *nTransmitters)
{
  int nt=0;
  int i;
  setTransmitter(759,10,"Unknown",&t1[nt],&t2[nt]); nt++;
  
  setTransmitter(763,10,"Optus",&t1[nt],&t2[nt]); nt++;
  setTransmitter(778,20,"Telstra",&t1[nt],&t2[nt]); nt++;
  setTransmitter(847.8,0.4,"Parkes radio",&t1[nt],&t2[nt]); nt++;
  setTransmitter(849.5,0.2,"Police",&t1[nt],&t2[nt]); nt++;
  setTransmitter(872.5,5,"Vodafone",&t1[nt],&t2[nt]); nt++;
  setTransmitter(882.5,14.9,"Telstra",&t1[nt],&t2[nt]); nt++;
  setTransmitter(905,3,"Amateur radio/walkie talkie",&t1[nt],&t2[nt]); nt++;
  setTransmitter(930.5,0.2,"Police",&t1[nt],&t2[nt]); nt++;
  setTransmitter(947.6,8.4,"Optus",&t1[nt],&t2[nt]); nt++;
  setTransmitter(956.2,5,"Vodafone",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1018,0.5,"PKS airport",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1575.42,32.736,"Beidou",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1176.45,20.46,"Beidou",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1268.52,20.46,"Beidou",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1542,34,"Inmarsat",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1643.25,33.5,"Inmarsat",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1575,32,"Galileo",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1280,40,"Galileo",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1189,50,"Galileo",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1575.42,15.456,"GPS",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1227.6,11,"GPS",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1176.45,12.5,"GPS",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1626.25,0.5,"Iridium",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1542.5,35,"Thurya",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1597.21875,16.3125,"Glonass",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1245.78125,5.6875,"Glonass",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1815,20,"Telstra",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1835,20,"Telstra",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1857.5,25,"Optus",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2115,10,"Vodafone",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2167.5,4.8,"Vodafone",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2122.5,5,"Telstra",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2127.5,5,"Telstra",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2142.5,5,"Optus",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2147.5,5,"Optus",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2152.5,5,"Optus",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2160.0,9.9,"Telstra",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2312,19.2,"NBN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2331.2,19.2,"NBN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2350.4,19.2,"NBN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2369.6,19.2,"NBN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2386.7,15,"NBN",&t1[nt],&t2[nt]); nt++;
  
  setTransmitter(2650.0,40,"Telstra",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2680.0,20,"Optus",&t1[nt],&t2[nt]); nt++;
  
  setTransmitter(3455.0,19.9,"NBN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(3560.0,19.9,"NBN",&t1[nt],&t2[nt]); nt++;
  
  setTransmitter(921.5,13,"WiFi 802.11",&t1[nt],&t2[nt]); nt++;
  setTransmitter(2442.0,82,"WiFi",&t1[nt],&t2[nt]); nt++;
  
  setTransmitter(1090,50./1000.,"ADS-B",&t1[nt],&t2[nt]); nt++;

  setTransmitter(1732,8,"UNKNOWN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1980,22,"UNKNOWN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1024,2,"UNKNOWN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(1018,2,"UNKNOWN",&t1[nt],&t2[nt]); nt++;
  setTransmitter(3072,2,"UNKNOWN",&t1[nt],&t2[nt]); nt++; 

  for (i=0;i<26;i++)
    {setTransmitter(704+128*i,12,"Band edge",&t1[nt],&t2[nt]); nt++; }

  
  *nTransmitters = nt;
}

void setTransmitter(float freq,float bw,char *name,float *t1,float *t2)
{
  *t1 = freq-bw/2.;
  *t2 = freq+bw/2.;
}

