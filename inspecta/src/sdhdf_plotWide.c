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
// Software to plot a spectrum
//
// Usage:
// sdhdf_plotSpectrum -f <filename>
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>

#define VNUM "v0.1"
#define MAX_SHADE 128
#define MAX_FLAG 128
#define MAX_OVERLAY_COEFF 8


void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY);
void plotSpectrum(sdhdf_fileStruct *inFile,int iband,int idump,char *grDev,char *fname,float f0,float f1,int av,int sump,int nx,int ny,int polPlot);
void showTransmitter(float freq,float bw,char *label,float miny,float maxy);
void drawShades(int nShade,float *shadeF0,float *shadeF1,int *shadeCol,float miny,float maxy);
void drawBand(float f1,float f2,int nVals,float *px,float *pflag,float *py1,float *py2,
	      float *py3,float *py4,int plotPol,float *flagF0,float *flagF1,int nFlag,
	      int nShade,float *shadeF0,float *shadeF1,int *shadeCol,int log,int labelit,float miny,float maxy,int setMinMax,char *title,char *ylabel);

void help()
{
  printf("sdhdf_plotWide\n\n");
  printf("Command line arguments\n\n");
  printf("-4pol                 Plot 4 polarisations (defaults to 2)\n");
  printf("-ch <val>             Set the font size to val\n");
  printf("-dispAz               Display as a function of Azimuth\n");
  printf("-dispTime             Display a time label on a splitRF plot\n");
  printf("-flag <val1> <val2>   Flag data between val1 and val2 (MHz)\n");
  printf("-g <grDevice>         Set PGPLOT graphics device for plot\n");
  printf("-maxhold              Plot a 'max hold' spectrum\n");
  printf("-maxx <val>           Maximum frequency to plot (MHz)\n");
  printf("-maxy <val>           Maximumy-value for plot\n");
  printf("-minx <val>           Minimum frequency to plot (MHz)\n");
  printf("-miny <val>           Minimum y-value for plot\n");
  printf("-nolog                Do not take the logarithm of the values\n");
  printf("-oCoeff <val>         Overlay curve coefficient\n");
  printf("-sb <val>             Select specific sub-band\n");
  printf("-sd <val>             Select specifc spectral dump\n");
  printf("-shade <f0> <f1> <col> Shade region with colour 'col' between F0 and F1 in MHz\n");
  printf("-splitRF              Split the plot into the 3 RF bands for the Parkes UWL receiver\n");
  printf("-stokes               Form Stokes parameters\n");
  printf("-transmitters         Plot known RFI transmitters\n");
  printf("-waterfall            Plot as a waterfall plot\n");
  printf("-ylabel <str>         Label on the y-axis\n");
  
  printf("\n");
  printf("Example usage: sdhdf_plotWide -nolog onoff.hdf\n");

  printf("\n\n");
  printf("Interactive key presses\n\n");
  printf("left mouse click      Identify cursor position\n");
  printf("h                     This help\n");
  printf("m                     Toggle plotting the max hold spectrum\n");
  printf("M                     Toggle plotting rest frequencies for molecular lines\n");
  printf("q                     Quit\n");
  printf("t                     Plot toggling transmitters\n");
  printf("u                     Un-zoom back to original scale\n");
  printf("w                     Set new minimum and maximum values\n");
  printf("z                     Zoom into specific part of the spectrum\n");
  

}


int main(int argc,char *argv[])
{  
  int        i,j,k,l,ii;
  char       fname[MAX_FILES][64];
  char fnameFix[1024];
  sdhdf_fileStruct *inFile;
  int        idump,iband,ibeam;
  float f0=-1,f1=-1;
  float overlay[MAX_OVERLAY_COEFF];
  int nOverlay;
  int av=0;
  int sump=0;
  int nx=1;
  int ny=1;
  int stokes=0;
  int polPlot=1;
  long nVals=0;
  float *px,*py1,*py2,*pflag,*py1_max,*py2_max;
  float *py3,*py4;
  int plotPol=2;
  int *psum;
  long writepos;
  int allocateMemory=0;
  float miny,maxy,val1,val2,val3,val4,minx,maxx;
  float miny1,maxy1,miny2,maxy2,miny3,maxy3;
  int useAllMiny=0;
  float ominy,omaxy,ominx,omaxx;
  float setMinX=-1;
  float setMaxX=-1;
  float plotTransmitters=-1;
  int plotMaxHold=-1;
  int molecularLines=-1;
  int inFiles;
  int sb=-1;
  float mx,my;
  char key;
  char grDev[128]="/xs";
  char title[128]="";
  char ylabel[128]="Signal strength (arbitrary)";
  int splitRF=0;
  int waterfall=0;
  int wNchan=0;
  int wNdump=0;
  float wF0,wF1;
  float wMin,wMax;
  float wFileChange[MAX_FILES];
  int wMax_dump=2058; // How to set this properly?
  float owMin,owMax;
  float *arr;

  int nShade=0;
  int shadeCol[MAX_SHADE];
  float shadeF0[MAX_SHADE];
  float shadeF1[MAX_SHADE];
  int sd=-1;
  int s0,s1;
  int dispAz=0;
  int dispTime=0;
  float azVal,elVal;
  char utc[1024];
  int log=1;
  int setMinMax=0;
  
  int nFlag=0;
  float flagF0[MAX_FLAG],flagF1[MAX_FLAG];
  int haveFlagged=0;
  float charHeight=1.0;

  miny = maxy = -1.0;

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
      if (strcmp(argv[i],"-minx")==0)
	sscanf(argv[++i],"%f",&setMinX);
      else if (strcmp(argv[i],"-ch")==0)
	sscanf(argv[++i],"%f",&charHeight);
      else if (strcasecmp(argv[i],"-oCoeff")==0)
	{
	  sscanf(argv[++i],"%f",&overlay[nOverlay]);
	  nOverlay++;
	}
      else if (strcmp(argv[i],"-h")==0) {help(); exit(1);}
      else if (strcmp(argv[i],"-maxx")==0)
	sscanf(argv[++i],"%f",&setMaxX);
      else if (strcmp(argv[i],"-miny")==0)
	{sscanf(argv[++i],"%f",&miny); setMinMax=1;}
      else if (strcmp(argv[i],"-maxy")==0)
	{sscanf(argv[++i],"%f",&maxy); setMinMax=1;}
      else if (strcasecmp(argv[i],"-3bandrange")==0)
	{
	  sscanf(argv[++i],"%f",&miny1);
	  sscanf(argv[++i],"%f",&maxy1);
	  sscanf(argv[++i],"%f",&miny2);
	  sscanf(argv[++i],"%f",&maxy2);
	  sscanf(argv[++i],"%f",&miny3);
	  sscanf(argv[++i],"%f",&maxy3);


	  useAllMiny=1;
	}
      else if (strcmp(argv[i],"-stokes")==0)
	stokes=1;
      else if (strcmp(argv[i],"-transmitters")==0)
	plotTransmitters=1;
      else if (strcmp(argv[i],"-ylabel")==0)
	strcpy(ylabel,argv[++i]);
      else if (strcmp(argv[i],"-waterfall")==0)
	waterfall=1;
      else if (strcmp(argv[i],"-maxhold")==0)
	plotMaxHold=1;
      else if (strcmp(argv[i],"-4pol")==0)
	plotPol=4;
      else if (strcmp(argv[i],"-splitRF")==0)
	splitRF=1;
      else if (strcmp(argv[i],"-splitRF_rfi")==0)
	splitRF=2;
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-dispAz")==0)
	dispAz=1;
      else if (strcmp(argv[i],"-dispTime")==0)
	dispTime=1;
      else if (strcasecmp(argv[i],"-nolog")==0)
	log=0;
      else if (strcmp(argv[i],"-sd")==0)
	sscanf(argv[++i],"%d",&sd);
      else if (strcmp(argv[i],"-sb")==0)
	sscanf(argv[++i],"%d",&sb);
      else if (strcmp(argv[i],"-flag")==0)
      	{
      	  sscanf(argv[++i],"%f",&flagF0[nFlag]);
      	  sscanf(argv[++i],"%f",&flagF1[nFlag++]);
      	}
      else if (strcmp(argv[i],"-shade")==0)
	{
	  sscanf(argv[++i],"%f",&shadeF0[nShade]);
	  sscanf(argv[++i],"%f",&shadeF1[nShade]);
	  sscanf(argv[++i],"%d",&shadeCol[nShade]);
	  nShade++;
	}
      else
	strcpy(fname[inFiles++],argv[i]);
    }
  
  printf("Have inFiles=  %d\n",inFiles);
  for (i=0;i<inFiles;i++)
    {
      printf("Loading %s\n",fname[i]);
      
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);

      sdhdf_loadMetaData(inFile);
      
      if (allocateMemory ==0)
	{
	  nVals=0;
	}
      if (sb==-1)
	{
	  for (j=0;j<inFile->beam[ibeam].nBand;j++)
	    {
	      sdhdf_loadBandData(inFile,ibeam,j,1);
	      nVals+=(inFile->beam[ibeam].bandHeader[j].nchan);
	    }
	}
      else
	{
	  sdhdf_loadBandData(inFile,ibeam,sb,1);
	  nVals+=(inFile->beam[ibeam].bandHeader[sb].nchan);
	}
      if (allocateMemory==0)
	{
	  // NOTE: SHOULD BE ABLE TO LOAD MULTPLE FILES ...
	  allocateMemory=1;
	  printf("Allocating memory for %d data points\n",(int)nVals);
	  if (waterfall==1)
	    arr = (float *)malloc(sizeof(float)*inFile->beam[ibeam].bandHeader[sb].nchan*wMax_dump);
	  else
	    {
	      px      = (float *)malloc(sizeof(float)*nVals);
	      pflag   = (float *)malloc(sizeof(float)*nVals);
	      py1     = (float *)calloc(sizeof(float),nVals);
	      py2     = (float *)calloc(sizeof(float),nVals);
	      py3     = (float *)calloc(sizeof(float),nVals);
	      py4     = (float *)calloc(sizeof(float),nVals);
	      py1_max = (float *)calloc(sizeof(float),nVals);
	      py2_max = (float *)calloc(sizeof(float),nVals);
	      psum    = (int *)calloc(sizeof(int),nVals);

	      writepos=0;
	      if (sb==-1)
		{
		  for (j=0;j<inFile->beam[ibeam].nBand;j++)
		    {
		      for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
			{
			  // FIX ME: using [0] for frequency dump
			  px[writepos+k] = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
			  pflag[writepos+k] = inFile->beam[ibeam].bandData[j].astro_data.dataWeights[k];
			}
		      writepos+=inFile->beam[ibeam].bandHeader[j].nchan;
		    }
		}
	      else
		{
		  for (k=0;k<inFile->beam[ibeam].bandHeader[sb].nchan;k++)
		    {
		      // FIX ME: using [0] for frequency dump
		      px[writepos+k] = inFile->beam[ibeam].bandData[sb].astro_data.freq[k];
		      pflag[writepos+k] = inFile->beam[ibeam].bandData[sb].astro_data.dataWeights[k];
		    }
		  writepos+=inFile->beam[ibeam].bandHeader[sb].nchan;	      
		}
	    }
	}


      if (waterfall==1)
	{
	  if (sb==-1)
	    {
	      printf("For waterfall plot, must set -sb <sub-band>\n");
	      waterfall=0;
	    }
	  else
	    {
	      float val;
	      printf("Loading waterfall plot data\n");	      
	      printf("wNdump = %d\n",wNdump);
	      printf("New ndump = %d\n",inFile->beam[ibeam].bandHeader[sb].ndump);
	      printf("New nchan = %d\n",inFile->beam[ibeam].bandHeader[sb].nchan);
	      printf("Val 0 = %g\n",inFile->beam[ibeam].bandData[sb].astro_data.pol1[0*inFile->beam[ibeam].bandHeader[sb].nchan+0]);
	      for (k=0;k<inFile->beam[ibeam].bandHeader[sb].ndump;k++)
		{
		  for (j=0;j<inFile->beam[ibeam].bandHeader[sb].nchan;j++)
		    {
		      val = inFile->beam[ibeam].bandData[sb].astro_data.pol1[k*inFile->beam[ibeam].bandHeader[sb].nchan+j]
			+ inFile->beam[ibeam].bandData[sb].astro_data.pol2[k*inFile->beam[ibeam].bandHeader[sb].nchan+j];

		      arr[(k+wNdump)*inFile->beam[ibeam].bandHeader[sb].nchan+j] = val;
		      if (k==0 && j==0)
			  wMin = wMax = val;
		      else
			{
			  if (wMin > val) wMin = val;
			  if (wMax < val) wMax = val;
			}
		    }
		}
	      // FIX ME: using [0] for frequency dump
	      wF0 = inFile->beam[ibeam].bandData[sb].astro_data.freq[0];
	      wF1 = inFile->beam[ibeam].bandData[sb].astro_data.freq[inFile->beam[ibeam].bandHeader[sb].nchan-1];
	      
	      printf("Complete loading waterfall, min/max = %g/%g\n",wMin,wMax);
	      wNchan =inFile->beam[ibeam].bandHeader[sb].nchan;
	      wNdump+=inFile->beam[ibeam].bandHeader[sb].ndump;
	      wFileChange[i]=wNdump;

	      printf("Now wNdump = %d\n",wNdump);
	    }

	}
      else
	{
	  
	  printf("Total number of values to plot =  %d\n",(int)nVals);
	  writepos=0;
	  for (j=0;j<inFile->beam[ibeam].nBand;j++)
	    {
	      if (sb==-1 || j == sb)
		{
		  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		    {
		      if (sd<0)
			{
			  s0 = 0;
			  s1 = inFile->beam[ibeam].bandHeader[j].ndump;
			}
		      else
			{
			  s0 = sd;
			  s1 = sd+1;
			}
		      for (l=s0;l<s1;l++)
			{
			  val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			  val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			  val3 = inFile->beam[ibeam].bandData[j].astro_data.pol3[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			  val4 = inFile->beam[ibeam].bandData[j].astro_data.pol4[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
			  py1[writepos+k] += val1;
			  py2[writepos+k] += val2;
			  py3[writepos+k] += val3;
			  py4[writepos+k] += val4;
			  if (val1 > py1_max[writepos+k]) py1_max[writepos+k] = val1;
			  if (val2 > py2_max[writepos+k]) py2_max[writepos+k] = val2;
			  psum[writepos+k] ++;
			}		  
		    }
		  writepos+=inFile->beam[ibeam].bandHeader[j].nchan;
		}
	    }
	}
      
      if (sd >= 0)
	{
	  azVal = inFile->beam[0].bandData[0].astro_obsHeader[sd].az;
	  elVal = inFile->beam[0].bandData[0].astro_obsHeader[sd].el;      
	}
	  strcpy(utc,inFile->beam[0].bandData[0].astro_obsHeader[sd].utc);
      printf("Closing\n");
      sdhdf_closeFile(inFile);
      printf("Done close\n");
    }
  if (waterfall==0)
    {
      for (i=0;i<nVals;i++)
	{
	  py1[i]/=(float)psum[i];
	  py2[i]/=(float)psum[i];
	  py3[i]/=(float)psum[i];
	  py4[i]/=(float)psum[i];
	  if (stokes==1)
	    {
	      val1=py1[i]; val2=py2[i]; val3=py3[i]; val4=py4[i];
	      py1[i] = val1+val2;
	      py2[i] = val1-val2;
	      py3[i] = 2*val3;
	      py4[i] = 2*val4;
	    }

	  if (log==1)
	    {
	      if (py1[i] > 0) py1[i] = log10(py1[i]);
	      else py1[i] = -3; // FIX ME -- SHOULD SET SENSIBLY

	      if (py2[i] > 0) py2[i] = log10(py2[i]);
	      else py2[i] = -3; // FIX ME -- SHOULD SET SENSIBLY

	      if (py3[i] > 0) py3[i] = log10(py3[i]);
	      else py3[i] = -3; // FIX ME -- SHOULD SET SENSIBLY
	      
	      if (py3[i] > 0) py3[i] = log10(py3[i]);
	      else py3[i] = -3; // FIX ME -- SHOULD SET SENSIBLY

	      py1_max[i] = log10(py1_max[i]);
	      py2_max[i] = log10(py2_max[i]);
	    }
	}
    }
  
  if (miny == -1 && maxy == -1 && waterfall == 0)
    {
      int useThis;
      miny =  1e30;
      maxy = -1e30;
      for (i=0;i<nVals;i++)
	{
	  useThis=1;
	  if (setMinX != -1 && px[i] < setMinX)
	    useThis=0;
	  if (setMaxX != -1 && px[i] > setMaxX)
	    useThis=0;
	  
	  if (useThis==1)
	    {
	      if (py1[i] > maxy) maxy = py1[i];
	      if (py2[i] > maxy) maxy = py2[i];
	      if (py1_max[i] > maxy) maxy = py1_max[i];
	      if (py2_max[i] > maxy) maxy = py2_max[i];
	      if (py1[i] < miny) miny = py1[i];
	      if (py2[i] < miny) miny = py2[i];
	    }
	}
    }

  if (sb > -1 && setMinX == -1)
    {
      setMinX = 704+sb*128;
      setMaxX = 704+(sb+1)*128;
    }
  
  printf("GOT TO HERE\n");  
  if (setMinX >= 0)
    minx = setMinX;
  else
    minx = px[0];
  printf("and here nVals = %d\n",(int)nVals);
  if (setMaxX >= 0)
    maxx = setMaxX;
  else
    maxx = px[nVals-1];
  printf("and here 2\n");
  if (waterfall==1)
    {
      miny = 0;
      maxy = wNdump;	
    }

  owMin = wMin;
  owMax = wMax;
  printf("miny/maxy = %g/%g\n",miny,maxy);
  printf("minx/maxx = %g/%g\n",minx,maxx);
  
  ominy = miny;
  omaxy = maxy;
  ominx = minx;
  omaxx = maxx;

  printf("Starting plot >%s<\n",grDev);
  
  cpgbeg(0,grDev,1,1);  
  cpgsch(1.4);
  cpgslw(2);
  cpgscf(2);
  cpgask(0);

  printf("Making plot\n");
  
  do {
    if (waterfall==1)
      {
	float tr[6];
	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
	float tx[2],ty[2];
	int nx,ny;

	nx = wNchan;
	ny = wNdump;

	tr[0] = wF0; tr[1] = (wF1-wF0)/nx; tr[2] = 0;
	tr[3] = 0; tr[4] = 0; tr[5] = 1;
	  

	printf("Waterfall with %d and %d\n",nx,ny);
	
	cpgeras();
	cpgsci(1);
	cpgsch(1.0);
	cpgsvp(0.10,0.95,0.10,0.95);
	cpgswin(minx,maxx,miny,maxy);
	cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
	cpglab("Frequency (MHz)","Time","");
	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	cpgimag(arr,nx,ny,1,nx,1,ny,wMin,wMax,tr);
	cpgsls(4);
	for (i=0;i<inFiles;i++)
	  {
	    tx[0] = minx; tx[1] = maxx;
	    ty[0] = ty[1] = wFileChange[i];
	    cpgline(2,tx,ty);
	  }
	cpgsls(1);
      }
    else if (splitRF==0)
      {
	cpgeras();
	cpgsci(1);
	//	if (log==1)
	//	  cpgenv(minx,maxx,miny,maxy,0,20);
	//	else
	//	  cpgenv(minx,maxx,miny,maxy,0,1);
	//	cpglab("Observing frequency (MHz)","Signal strength","");
	//	cpgsci(2);
	cpgsch(charHeight);
	if (charHeight < 1.1)
	  cpgsvp(0.10,0.95,0.10,0.90);
	else
	  cpgsvp(0.10,0.95,0.15,0.90); 
	sdhdf_fixUnderscore(fname[0],fnameFix);
	if (sb >= 0 && sd >= 0)
	  sprintf(title,"%s (sub-band %d, spectral-dump %d)",fnameFix,sb,sd);
	else
	  sprintf(title,"%s (spectral-dump %d)",fnameFix,sd);
	drawBand(minx,maxx,nVals,px,pflag,py1,py2,py3,py4,plotPol,flagF0,flagF1,nFlag,nShade,shadeF0,shadeF1,shadeCol,log,3,miny,maxy,1,title,ylabel);
	if (nOverlay > 0) // Draw a curve overlay	  
	  {
	    int nv=128;
	    float xp[nv],yp[nv];
	    for (i=0;i<nv;i++)
	      {
		xp[i] = minx + i*(maxx-minx)/(double)nv;
		yp[i] = 0.0;
		for (j=0;j<nOverlay;j++)
		  yp[i] += overlay[j]*pow(log10(xp[i]/1000.0),j);
		yp[i] = pow(10,yp[i]);
		printf("Have %g %g\n",xp[i],yp[i]);


	      }
		cpgslw(3);
                cpgline(nv,xp,yp);
                cpgslw(1);
	  }
	//	cpgline(nVals,px,py1);
	if (plotMaxHold==1)
	  {
	    cpgsci(8);
	    cpgsls(4); cpgline(nVals,px,py1_max); cpgsls(1);
	  }
	cpgsci(4);
	//	cpgline(nVals,px,py2);
	if (plotMaxHold==1)
	  {
	    cpgsci(5);
	    cpgsls(4); cpgline(nVals,px,py2_max); cpgsls(1);
	  }
	cpgsci(1);
	if (dispAz==1)
	  {
	    char label[128];
	    cpgsch(1.0);
	    sprintf(label,"Azimuth: %.1f deg",azVal);
	    cpgtext(minx + (maxx-minx)*0.05,maxy-(maxy-miny)*0.1,label);
	    cpgsch(1.4);
	  }

      }
    else
      {
	float useMiny,useMaxy;
	int i0,i1,t=0;
	int flagit;
	int region;
	
        cpgsch(charHeight);
	cpgsvp(0.10,0.95,0.08,0.31);

	if (useAllMiny==1)
	  {
	    miny=miny1;
	    maxy=maxy1;
	    setMinMax=1;
	  }

	if (splitRF==1)
	  drawBand(704,1344,nVals,px,pflag,py1,py2,py3,py4,plotPol,flagF0,flagF1,nFlag,nShade,shadeF0,shadeF1,shadeCol,log,1,miny,maxy,setMinMax,title,ylabel);
	else
	  drawBand(704,960,nVals,px,pflag,py1,py2,py3,py4,plotPol,flagF0,flagF1,nFlag,nShade,shadeF0,shadeF1,shadeCol,log,1,miny,maxy,setMinMax,title,ylabel);
	  
	cpgsvp(0.10,0.95,0.36,0.61);
	if (useAllMiny==1)
	  {
	    miny=miny2;
	    maxy=maxy2;
	    setMinMax=1;
	    printf("Setting to %g %g\n",miny,maxy);
	  }

	if (splitRF==1)
	  drawBand(1344,2368,nVals,px,pflag,py1,py2,py3,py4,plotPol,flagF0,flagF1,nFlag,nShade,shadeF0,shadeF1,shadeCol,log,2,miny,maxy,setMinMax,title,ylabel);
	else
	  drawBand(960,1216,nVals,px,pflag,py1,py2,py3,py4,plotPol,flagF0,flagF1,nFlag,nShade,shadeF0,shadeF1,shadeCol,log,2,miny,maxy,setMinMax,title,ylabel);
	  
	cpgsvp(0.10,0.95,0.66,0.91);
	if (useAllMiny==1)
	  {
	    miny=miny3;
	    maxy=maxy3;
	    setMinMax=1;
	  }
	if (splitRF==1)
	  drawBand(2368,4096,nVals,px,pflag,py1,py2,py3,py4,plotPol,flagF0,flagF1,nFlag,nShade,shadeF0,shadeF1,shadeCol,log,0,miny,maxy,setMinMax,title,ylabel);
	else
	  drawBand(2368,2496,nVals,px,pflag,py1,py2,py3,py4,plotPol,flagF0,flagF1,nFlag,nShade,shadeF0,shadeF1,shadeCol,log,0,miny,maxy,setMinMax,title,ylabel);
	cpgsch(1.4);

	cpgsvp(0.10,0.95,0.91,0.99);
	cpgswin(0,1,0,1);
	//	cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
	
	//	if (dispAz==1)
	  {
	    char label[128];
	    cpgsch(1.0);
	    sprintf(label,"Azimuth: %.1f deg",azVal);
	    cpgtext(0.07,0.8,label);
	    sprintf(label,"Elevation: %.1f deg",elVal);
	    cpgtext(0.07,0.4,label);
	    cpgsch(1.4);
	  }
	  //	if (dispTime==1)
	  {
	    char label[128];
	    sdhdf_fixUnderscore(fname[0],fnameFix);
	    cpgsch(1.0);
	    sprintf(label,"UTC time: %s",utc);
	    cpgtext(0.3,0.8,label);
	    sprintf(label,"File: %s",fnameFix);
	    cpgtext(0.6,0.8,label);
	    sprintf(label,"Spectral dump: %d",sd);
	    cpgtext(0.6,0.4,label);


	    cpgsch(1.4);
	  }

      }
  
    if (molecularLines==1)
      {
	drawMolecularLine(1420.405752,"HI",minx,maxx,miny,maxy);
	drawMolecularLine(1612.2310,"OH",minx,maxx,miny,maxy);
	drawMolecularLine(1665.4018,"OH",minx,maxx,miny,maxy);
	drawMolecularLine(1667.3590,"OH",minx,maxx,miny,maxy);
	drawMolecularLine(1720.5300,"OH",minx,maxx,miny,maxy);
      }

    if (plotTransmitters==1)
      {
	cpgsch(0.8);
	showTransmitter(763,10,"Optus",miny,maxy);
	showTransmitter(778,20,"Telstra",miny,maxy);
	showTransmitter(847.8,0.4,"Parkes radio",miny,maxy);
	showTransmitter(849.5,0.2,"Police",miny,maxy);
	showTransmitter(872.5,5,"Vodafone",miny,maxy);
	showTransmitter(882.5,14.9,"Telstra",miny,maxy);
	showTransmitter(905,3,"Amateur radio/walkie talkie",miny,maxy);
	showTransmitter(930.5,0.2,"Police",miny,maxy);
	showTransmitter(947.6,8.4,"Optus",miny,maxy);
	showTransmitter(956.2,5,"Vodafone",miny,maxy);
	showTransmitter(1018,0.5,"PKS airport",miny,maxy);
	showTransmitter(1575.42,32.736,"Beidou",miny,maxy);
	showTransmitter(1176.45,20.46,"Beidou",miny,maxy);
	showTransmitter(1268.52,20.46,"Beidou",miny,maxy);
	showTransmitter(1542,34,"Inmarsat",miny,maxy);
	showTransmitter(1643.25,33.5,"Inmarsat",miny,maxy);
	showTransmitter(1575,32,"Galileo",miny,maxy);
	showTransmitter(1280,40,"Galileo",miny,maxy);
	showTransmitter(1189,50,"Galileo",miny,maxy);
	showTransmitter(1575.42,15.456,"GPS",miny,maxy);
	showTransmitter(1227.6,11,"GPS",miny,maxy);
	showTransmitter(1176.45,12.5,"GPS",miny,maxy);
	showTransmitter(1626.25,0.5,"Iridium",miny,maxy);
	showTransmitter(1542.5,35,"Thurya",miny,maxy);
	showTransmitter(1597.21875,16.3125,"Glonass",miny,maxy);
	showTransmitter(1245.78125,5.6875,"Glonass",miny,maxy);
	showTransmitter(1815,20,"Telstra",miny,maxy);
	showTransmitter(1835,20,"Telstra",miny,maxy);
	showTransmitter(1857.5,25,"Optus",miny,maxy);
	showTransmitter(2115,10,"Vodafone",miny,maxy);
	showTransmitter(2167.5,4.8,"Vodafone",miny,maxy);
	showTransmitter(2122.5,5,"Telstra",miny,maxy);
	showTransmitter(2127.5,5,"Telstra",miny,maxy);
	showTransmitter(2142.5,5,"Optus",miny,maxy);
	showTransmitter(2147.5,5,"Optus",miny,maxy);
	showTransmitter(2152.5,5,"Optus",miny,maxy);
	showTransmitter(2160.0,9.9,"Telstra",miny,maxy);
	showTransmitter(2312,19.2,"NBN",miny,maxy);
	showTransmitter(2331.2,19.2,"NBN",miny,maxy);
	showTransmitter(2350.4,19.2,"NBN",miny,maxy);
	showTransmitter(2369.6,19.2,"NBN",miny,maxy);
	showTransmitter(2386.7,15,"NBN",miny,maxy);

	showTransmitter(2650.0,40,"Telstra",miny,maxy);
	showTransmitter(2680.0,20,"Optus",miny,maxy);

	showTransmitter(3455.0,19.9,"NBN",miny,maxy);
	showTransmitter(3560.0,19.9,"NBN",miny,maxy);
	
	showTransmitter(921.5,13,"WiFi 802.11",miny,maxy);
	showTransmitter(2442.0,82,"WiFi",miny,maxy);

	showTransmitter(1090,50./1000.,"ADS-B",miny,maxy);

	cpgsch(1.4);
      }

    if (strcmp(grDev,"/xs")==0)
      {
	cpgcurs(&mx,&my,&key);
	
	if (key=='z')
	  {
	    float mx2,my2;
	    cpgband(2,0,mx,my,&mx2,&my2,&key);
	    printf("Zooming in on %g %g %g %g\n",mx,my,mx2,my2);
	    if (mx != mx2 && my != my2)
	      {
		if (mx < mx2)
		  {minx = mx; maxx = mx2;}
		else
		  {minx = mx2; maxx = mx;}
		
		if (my < my2)
		  {miny = my; maxy = my2;}
		else
		  {miny = my2; maxy = my;}
	      }
	    else
	      printf("Please press 'z' and then move somewhere and click left mouse button\n");
	  }
	else if (key=='h') help();
	else if (key=='A')
	  printf("Mouse press at (%.6f,%g)\n",mx,my);		 
	else if (key=='m')
	  plotMaxHold*=-1;
	else if (key=='w')
	  {
	    printf("Min/max = %g %g\n",wMin,wMax);
	    printf("Enter new values:\n");
	    scanf("%f %f",&wMin,&wMax);
	  }
	else if (key=='M')
	  molecularLines*=-1;
	else if (key=='t')
	  plotTransmitters*=-1;
	else if (key=='u')
	  {
	    miny = ominy;
	    maxy = omaxy;
	    minx = ominx;
	    maxx = omaxx;
	  }
      }
    else
      key='q';

  } while (key!='q');

  
  cpgend();  

  if (waterfall==1)
    free(arr);
  else
    {
      free(px);
      free(pflag);
      free(py1);
      free(py2);
      free(py3);
      free(py4);
      free(py1_max);
      free(py2_max);      
      free(psum);
    }
  free(inFile);
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


void drawShades(int nShade,float *shadeF0,float *shadeF1,int *shadeCol,float miny,float maxy)
{
  int i;
  cpgsfs(4);
  for (i=0;i<nShade;i++)
    {
      cpgsci(shadeCol[i]);
      cpgrect(shadeF0[i],shadeF1[i],miny,maxy);
      cpgsci(1);
    }
  cpgsfs(1);
}

void drawBand(float f1,float f2,int nVals,float *px,float *pflag,float *py1,float *py2,float *py3,float *py4,int plotPol,float *flagF0,float *flagF1,int nFlag,
	      int nShade,float *shadeF0,float *shadeF1,int *shadeCol,int log,int labelit,float miny,float maxy,int setMinMax,char *title,char *ylabel)
{
  int i,i0,i1,ii,flagit;
  float useMiny,useMaxy;
  int t,region;

  if (setMinMax==1)
    {
      useMiny = miny;
      useMaxy = maxy;
    }
  else
    {
      t=0;
      for (i=0;i<nVals;i++)
	{	    
	  if (px[i] >= f1 && px[i] < f2)
	    {
	      flagit=0;
	      if (pflag[i]==0)
		flagit=1;
	      else
		{
		  for (ii=0;ii<nFlag;ii++)
		    {
		      if (px[i] >= flagF0[ii] && px[i] < flagF1[ii])
			{
			  flagit=1;
			  break;
			}
		    }
		}
	      if (flagit==0)
		{
		  if (t==0)
		    {
		      useMiny = py1[i];
		      useMaxy = py1[i];
		      t=1;
		    }
		  else
		    {
		      if (useMiny > py1[i]) useMiny = py1[i];
		      if (useMiny > py2[i]) useMiny = py2[i];			
		      if (useMaxy < py1[i]) useMaxy = py1[i];
		      if (useMaxy < py2[i]) useMaxy = py2[i];			
		    }
		}
	    }
	}
    }
  printf("The minimum and maximum values are: %g/%g\n",useMiny,useMaxy);
  cpgswin(f1,f2,useMiny,useMaxy);

  if (log==1)
    cpgbox("ABCTSN",0,0,"ABCTSNL",0,0);
  else
    cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
  if (labelit==1)
    cpglab("Frequency (MHz)","",title);
  else if (labelit==2)
    cpglab("",ylabel,title);
  else if (labelit==3)
    cpglab("Frequency (MHz)",ylabel,title);
  drawShades(nShade,shadeF0,shadeF1,shadeCol,useMiny,useMaxy);

  i0=0;	
  region=0;
  printf("nVals = %d\n",nVals);
  for (i=0;i<nVals;i++)
    {
      flagit=0;
      if (pflag[i]==0)
	flagit=1;
      else	
	{
	  for (ii=0;ii<nFlag;ii++)
	    {
	      if (px[i] >= flagF0[ii] && px[i] < flagF1[ii])
		{
		  flagit=1;
		  break;
		}
	    }
	}
      //      printf("Here with %d %d\n",flagit,region);
      if (flagit==0 && region==2)
	{
	  region=0;
	  i0=i;
	}
      if (flagit==1 && region==0)
	{
	  i1 = i;
	  cpgsci(2);
	  cpgline(i1-i0,px+i0,py1+i0);
	  cpgsci(4);
	  cpgline(i1-i0,px+i0,py2+i0);
	  if (plotPol==4)
	    {
	      cpgsci(5);
	      cpgline(i1-i0,px+i0,py3+i0);
	      cpgsci(8);
	      cpgline(i1-i0,px+i0,py4+i0);
	    }
	  cpgsci(1);	  
	  region=2;
	  printf("Plotting %d\n",i1-i0);
	}
    }
  if (region==0)
    {
      i1 = i-1;
      cpgsci(2);
      cpgline(i1-i0,px+i0,py1+i0);
      cpgsci(4);
      cpgline(i1-i0,px+i0,py2+i0);
      if (plotPol==4)
	{
	  cpgsci(5);
	  cpgline(i1-i0,px+i0,py3+i0);
	  cpgsci(8);
	  cpgline(i1-i0,px+i0,py4+i0);
	}
      cpgsci(1);
      
    }
  

}
