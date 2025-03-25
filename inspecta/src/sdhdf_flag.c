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
// Software to flag data
//
// Usage:
// sdhdf_cal
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"
#include <cpgplot.h>
#include "TKfit.h"

#define MAX_RFI 512

void doPlot(sdhdf_fileStruct *inFile,int ibeam);
void saveFile(sdhdf_fileStruct *inFile,int ibeam);
void zapPersistent(sdhdf_fileStruct *inFile,int ibeam,int zapAll);
void autoZapTransmitters(sdhdf_fileStruct *inFile,int ibeam,int zapAll);
void autoZapDigitisers(sdhdf_fileStruct *inFile,int ibeam,int zapAll);
void autoZapAircraft(sdhdf_fileStruct *inFile,int ibeam,int zapAll);
void autoZapSatellites(sdhdf_fileStruct *inFile,int ibeam,int zapAll);
void autoZapWiFi(sdhdf_fileStruct *inFile,int ibeam,int zapAll);
void autoZapUnexplained(sdhdf_fileStruct *inFile,int ibeam,int zapAll);
void autoZapHandsets(sdhdf_fileStruct *inFile,int ibeam);
void flagChannels(sdhdf_fileStruct *inFile,int ibeam,float *f0,float *f1,int nT,int zapAll);

void help()
{
  printf("sdhdf_flag\n");
  printf("sdhfProc version:   %s\n",SOFTWARE_VER);
  printf("author:             George Hobbs\n");
  printf("\n");
  printf("Software to flag a spectrum in an interactive matter\n");
  printf("\n\nCommand line arguments:\n\n");

  printf("-f <string>          Input file name\n");
  printf("-h                   This help\n");
  printf("-setnband <number>   Set the maximum number of sub-bands to process\n");
  printf("-sb <number>         Set the band number for the initial plot (starting from 0)   [can also use -band instead of -sb]\n");

  printf("\n\n");
  printf("The plot is interactive and the following commands are available\n\n");
  
  printf("1        Select only polarisation channel 1\n");
  printf("2        Select only polarisation channel 2\n");
  printf("+        Move to next spectral dump\n");
  printf("-        Move to previous spectral dump\n");
  printf(">        Move to next sub-band\n");
  printf("<        Move to previous sub-band\n");
  printf("a        Toggle between flagging all the spectral dumps (default) or just the current one\n");
  printf("b        Flag 5 percent of the Parkes UWL band edges\n");
  printf("d        Toggle between showing and hiding (default) the flagged channels\n");
  printf("f (or Z) Define a region to flag specific channels\n");
  printf("h        This help\n");
  printf("H        Flag all channels where the value is above a specific value\n");
  printf("j        Jump to specific spectral dump\n");
  printf("k        Flag all the channels in the current zoom window for the current spectral dump and move on to the next\n");
  printf("K        Flag entire spectral dump and move on to the next\n");
  printf("p        Remove Parkes UWL persistent RFI\n");
  printf("q        Quit\n");
  printf("r        Remove flagging in specified region\n");
  printf("s        Save file\n");
  printf("t        Automatically identify and flag mobile handset RFI\n");
  printf("u        Unzoom - revert to default zoom level\n");
  printf("z        Zoom into a region using the mouse\n");
  printf("Z (or f) Define a region to flag specific channels\n");
}

int main(int argc,char *argv[])
{
  int i,j,k;
  char fname[MAX_STRLEN];
  sdhdf_fileStruct *inFile;
  int iband=0,ibeam=0,nchan,idump,nband;
  int ndump=0,totSize,chanPos;
  double tint,chbw;
  //  spectralDumpStruct spectrum;
  int npol=4;
  int setnband=-1;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);
    
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-f")==0)	       strcpy(fname,argv[++i]);	
      else if (strcmp(argv[i],"-setnband")==0) sscanf(argv[++i],"%d",&setnband);
      else if (strcmp(argv[i],"-sb")==0 || strcmp(argv[i],"-band")==0)       sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-h")==0) {help(); exit(1);}
    }
  
  sdhdf_openFile(fname,inFile,1);
  printf("File opened\n");
  sdhdf_loadMetaData(inFile);
  printf("Complete loading metadata\n");
  if (setnband > 0)
    inFile->beam[ibeam].nBand = setnband;
  
  
  totSize=0;
  nband = inFile->beam[ibeam].nBand;
   
  npol=4; // FIX HARDCODING OF NPOL
  for (i=0;i<nband;i++)
    {
      totSize+=(inFile->beam[ibeam].bandHeader[i].nchan*inFile->beam[ibeam].bandHeader[i].ndump*npol);
    }
  
  printf("Loading the data\n");
  for (i=0;i<nband;i++)
    {
      nchan = inFile->beam[ibeam].bandHeader[i].nchan;
      ndump = inFile->beam[ibeam].bandHeader[i].ndump;
      tint  = inFile->beam[ibeam].bandHeader[i].dtime;
      chbw  = 1e6*fabs(inFile->beam[ibeam].bandHeader[i].f0-inFile->beam[ibeam].bandHeader[i].f1)/(double)inFile->beam[ibeam].bandHeader[i].nchan;
      sdhdf_loadBandData(inFile,ibeam,i,1);
    }
  printf("Complete loading the data\n");
  printf("Preparing the plot\n");
  doPlot(inFile,ibeam);

  sdhdf_closeFile(inFile);

  free(inFile);
}

void doPlot(sdhdf_fileStruct *inFile,int ibeam)
{
  int plot=1,i,j,k;
  float *px,*py1,*py2,*pf1,*pf2;
  float pointX[2],pointY[2];
  float maxPosX;
  float mx,my,mx2,my2;
  float minx,maxx,miny,maxy;
  char key;
  char title[MAX_STRLEN];
  int iband=0;
  int idump=0;
  int nchan,ndump;
  long max_nchan;
  int npol=4;
  int recalc=1;
  int nFlagRegion=0;
  int n;
  int c0;
  int regionType=1;
  int zapAllDumps=1;
  int deleteFlagged=1;
  int selectPol=1;
  int firstThrough=1;
  int needSave=0;
  
  cpgbeg(0,"/xs",1,1);
  cpgslw(2);
  cpgscf(2);
  cpgsch(1.4);
  cpgask(0);
  printf("Doing plot\n");


  max_nchan=0;
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      if (max_nchan < inFile->beam[ibeam].bandHeader[i].nchan)
	max_nchan = inFile->beam[ibeam].bandHeader[i].nchan;
    }
  printf("Maximum channels = %d\n",(int)max_nchan);
  
  px  = (float *)malloc(sizeof(float)*max_nchan);
  py1 = (float *)malloc(sizeof(float)*max_nchan);
  py2 = (float *)malloc(sizeof(float)*max_nchan);
  pf1 = (float *)malloc(sizeof(float)*max_nchan);
  pf2 = (float *)malloc(sizeof(float)*max_nchan);
  
  do {    
    nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
    ndump = inFile->beam[ibeam].bandHeader[iband].ndump;

    if (recalc > 0)
      {
	if (zapAllDumps==1)
	  printf("Flagging will apply to all spectral dumps (press 'a' to toggle)\n");
	else
	  printf("Flagging will only apply to current spectral dump (press 'a' to toggle)\n");
	printf("Recalc ID = %d\n",recalc);
	regionType=1;
	nFlagRegion=0;
	c0 = 0;
	for (i=0;i<iband;i++)
	  c0+=(inFile->beam[ibeam].bandHeader[iband].nchan*inFile->beam[ibeam].bandHeader[iband].ndump*npol);

	if (recalc!=3)
	  {
	    if (recalc!=2) {minx = 1e99; maxx = -1e99;}
	    miny = 1e99; maxy = -1e99;
	  }
	for (i=0;i<nchan;i++)
	  {
	    // FIX ME: using [0] for frequency dump
	    px[i]  = inFile->beam[ibeam].bandData[iband].astro_data.freq[i];
	    py1[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol1[i+idump*nchan];
	    py2[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol2[i+idump*nchan];

	    /*
	    for (j=0;j<ndump;j++)
	      {
		py1[i] += inFile->beam[ibeam].bandData[iband].astro_data.pol1[i+j*nchan];
		py2[i] += inFile->beam[ibeam].bandData[iband].astro_data.pol2[i+j*nchan];
	      }
	    */
	    if (recalc!=2 && recalc!=3)
	      {
		if (px[i] > maxx) maxx = px[i];
		if (px[i] < minx) minx = px[i];
	      }
	  }
	for (i=0;i<nchan;i++)
	  {
	    if (recalc!=3)
	      {
		if (inFile->beam[ibeam].bandData[iband].astro_data.flag[i+nchan*idump] != 1)
		  {
		    if (px[i] > minx && px[i] < maxx)
		      {
			if (py1[i] > maxy) {maxy = py1[i]; maxPosX = px[i];}
			if (py2[i] > maxy) {maxy = py2[i]; maxPosX = px[i];}
			if (py1[i] < miny) {miny = py1[i];}
			if (py2[i] < miny) {miny = py2[i];}
		      }
		  }
	      }
	      
	    
	    if (regionType==1 && inFile->beam[ibeam].bandData[iband].astro_data.flag[i+nchan*idump] == 1) // No region set
	      {
		//		printf("Flagging region\n");
		regionType=2;
		pf1[nFlagRegion] = px[i];
	      }	     
	    else if (regionType==2 && inFile->beam[ibeam].bandData[iband].astro_data.flag[i+nchan*idump] == 1) 
	      {
		regionType=1;
		pf2[nFlagRegion++] = px[i];
	      }
	      
	  }
	if (regionType==2)
	    pf2[nFlagRegion++] = px[i-1];

	n = nchan;
	recalc=0;
      }
    printf("Number of flag regions = %d\n",nFlagRegion);
    cpgenv(minx,maxx,miny,maxy,0,1);
    sprintf(title,"Band: %s (spectral dump %d/%d)",inFile->beam[ibeam].bandHeader[iband].label,idump,ndump-1);
    cpglab("Frequency (MHz)","Signal strength (arbitrary)",title);

    // Update so we don't show zapped regions
    if (deleteFlagged==-1)
      {
	cpgline(n,px,py1);
	cpgsci(2); cpgline(n,px,py2); cpgsci(1);
      }
    else
      {
	int i0=0;
	int drawIt=-1;
	for (i=0;i<nchan;i++)
	  {
	    if (inFile->beam[ibeam].bandData[iband].astro_data.flag[i+nchan*idump] != 1 && drawIt==-1)
	      {
		drawIt=1;
		i0=i;
	      }
	    if (inFile->beam[ibeam].bandData[iband].astro_data.flag[i+nchan*idump] == 1 && drawIt==1)
	      {
		cpgsci(1); cpgline(i-1-i0,px+i0,py1+i0); cpgsci(1);
		cpgsci(2); cpgline(i-1-i0,px+i0,py2+i0); cpgsci(1);
		drawIt=-1;
		
	      }	      
	  }
	if (drawIt==1)
	  {
	    cpgsci(1); cpgline(i-1-i0,px+i0,py1+i0); cpgsci(1);
	    cpgsci(2); cpgline(i-1-i0,px+i0,py2+i0); cpgsci(1);
	  }
      }

    cpgsfs(3);
    cpgsci(5);
    for (i=0;i<nFlagRegion;i++)
      cpgrect(pf1[i],pf2[i],miny,maxy);
    cpgsci(1);
    cpgsfs(1);

    pointX[0] = maxPosX;
    pointY[0] = maxy;
    cpgsci(3);
    cpgpt(1,pointX,pointY,31);
    cpgsci(1);
    
    cpgcurs(&mx,&my,&key);
    if (key=='u')
      recalc=1;
    else if (key=='a')
      {
	zapAllDumps*=-1;
	recalc=1;
      }
    else if (key=='d')
      deleteFlagged*=-1;
    else if (key=='b') // Zap 5% of the band edges -- assuming Parkes UWL data
      {
	int ii,jj,zz;
	needSave=1;
	printf("Zapping band edges\n");
	for (ii=0;ii<inFile->beam[ibeam].nBand;ii++)
	  {
	    for (jj=0;jj<inFile->beam[ibeam].bandHeader[ii].nchan;jj++)
	      {
		for (zz=0;zz<=26;zz++)
		  {
		    // FIX ME: using [0] for frequency dump
		    if (inFile->beam[ibeam].bandData[ii].astro_data.freq[jj] >= 704-6.4+zz*128.0 &&
			inFile->beam[ibeam].bandData[ii].astro_data.freq[jj] < 704+6.4+zz*128.0)
		      {
			if (zapAllDumps==1)
			  {
			    for (k=0;k<ndump;k++)
			      inFile->beam[ibeam].bandData[ii].astro_data.flag[jj+k*inFile->beam[ibeam].bandHeader[ii].nchan] = 1;
			  }
			else
			  {
			    inFile->beam[ibeam].bandData[ii].astro_data.flag[jj+idump*inFile->beam[ibeam].bandHeader[ii].nchan] = 1;
			  }
		      }
		  }
	      }	
	  }
	recalc=1;
      }
    else if (key=='t') // Zap Parkes transient RFI
      {
	autoZapHandsets(inFile,ibeam);
	recalc=2;
      }
    else if (key=='p') // Zap for Parkes UWL data for persistent RFI
      {
	needSave=1;
	zapPersistent(inFile,ibeam,zapAllDumps);
	//	autoZapTransmitters(inFile,ibeam,zapAllDumps);
	//	autoZapDigitisers(inFile,ibeam,zapAllDumps);



	//	autoZapAircraft(inFile,ibeam,zapAllDumps);
	//	autoZapSatellites(inFile,ibeam,zapAllDumps);
	//	autoZapWiFi(inFile,ibeam,zapAllDumps);
	//	autoZapUnexplained(inFile,ibeam,zapAllDumps);
	//	autoZapHandsets(inFile,ibeam,zapAllDumps);
	recalc=1;
      }
    else if (key=='1')
      selectPol=1;
    else if (key=='2')
      selectPol=2;
    else if (key=='s')
      {
	saveFile(inFile,ibeam);
	needSave=0;
      }
    else if (key=='K') // Zap entire spectral dump and move to the next one
      {
	needSave=1;
	for (i=0;i<inFile->beam[ibeam].nBand;i++)
	  {
	    for (j=0;j<inFile->beam[ibeam].bandHeader[i].nchan;j++)
	      inFile->beam[ibeam].bandData[i].astro_data.flag[j+idump*inFile->beam[ibeam].bandHeader[i].nchan] = 1;      
	  }
	idump++;
	if (idump >= ndump)
	  {
	    idump = ndump-1;
	    recalc=1;
	  }
	else
	  recalc=2;
      }
    else if (key=='k') // Zap all channels in current spectral dump zoom region and move to the next one
      {
	for (i=0;i<inFile->beam[ibeam].nBand;i++)
	  {
	    for (j=0;j<inFile->beam[ibeam].bandHeader[i].nchan;j++)
	      { 
		// FIX ME: using [0] for frequency dump
		if (inFile->beam[ibeam].bandData[i].astro_data.freq[j] >= minx && inFile->beam[ibeam].bandData[i].astro_data.freq[j] <= maxx)
		  inFile->beam[ibeam].bandData[i].astro_data.flag[j+idump*inFile->beam[ibeam].bandHeader[i].nchan] = 1;      
	      }
	  }
	idump++;
	if (idump >= ndump)
	  {
	    idump = ndump-1;
	    recalc=1;
	  }
	else
	  recalc=2;
      }
    else if (key=='j') // Jump to dump
      {
	printf("Enter spectral dump number (starting from 0): ");
	scanf("%d",&idump);
	if (idump >= ndump)
	  idump = ndump-1;
	recalc=2;
      }
    else if (key=='+') 
      {
	idump++;
	if (idump >= ndump)
	  idump = ndump-1;
	recalc=2;
      }
    else if (key=='-') 
      {
	idump--;
	if (idump == -1)
	  idump=0;
	recalc=2;
      }
    else if (key=='>') // Note changed from previously +/-
      {
	iband++;
	if (iband >= inFile->beam[ibeam].nBand)
	  iband = inFile->beam[ibeam].nBand-1;
	recalc=1;
      }
    else if (key=='<')
      {
	iband--;
	recalc=1;
	if (iband < 0) iband=0;
      }
    else if (key=='z')
      {
	cpgband(2,0,mx,my,&mx2,&my2,&key);
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
    else if (key=='Z' || key=='f')
      {
	float lowX,highX;
	needSave=1;
	cpgband(4,0,mx,my,&mx2,&my2,&key);
	if (mx != mx2)
	  {
	    if (mx < mx2)
	      { lowX = mx; highX = mx2;}
	    else
	      { lowX = mx2; highX = mx;}
	    for (i=0;i<inFile->beam[ibeam].nBand;i++)
	      {
		for (j=0;j<inFile->beam[ibeam].bandHeader[i].nchan;j++)
		  {
		    // FIX ME: using [0] for frequency dump
		    if (inFile->beam[ibeam].bandData[i].astro_data.freq[j] >= lowX && inFile->beam[ibeam].bandData[i].astro_data.freq[j] <= highX)
		      {
			if (zapAllDumps==1)
			  {
			    for (k=0;k<ndump;k++)
			      inFile->beam[ibeam].bandData[i].astro_data.flag[j+k*inFile->beam[ibeam].bandHeader[i].nchan] = 1;
			  }
			else
			  inFile->beam[ibeam].bandData[i].astro_data.flag[j+idump*inFile->beam[ibeam].bandHeader[i].nchan] = 1;  
		      }
		  }
	      }
	    recalc=2;
	  }
      }
    else if (key=='r') // Un-flag data
      {
	float chbw;
	float lowX,highX;
	needSave=1;
	cpgband(4,0,mx,my,&mx2,&my2,&key);
	if (mx != mx2)
	  {
	    if (mx < mx2)
	      { lowX = mx; highX = mx2;}
	    else
	      { lowX = mx2; highX = mx;}
	    for (i=0;i<inFile->beam[ibeam].nBand;i++)
	      {
		for (j=0;j<inFile->beam[ibeam].bandHeader[i].nchan;j++)
		  { 
		    // FIX ME: using [0] for frequency dump
		    if (inFile->beam[ibeam].bandData[i].astro_data.freq[j] >= lowX && inFile->beam[ibeam].bandData[i].astro_data.freq[j] <= highX)
		      {
			chbw = 1e6*fabs(inFile->beam[ibeam].bandHeader[i].f0-inFile->beam[ibeam].bandHeader[i].f1)/(double)inFile->beam[ibeam].bandHeader[iband].nchan;
			if (zapAllDumps==1)
			  {
			    for (k=0;k<ndump;k++)
			      {				
				inFile->beam[ibeam].bandData[i].astro_data.flag[j+k*inFile->beam[ibeam].bandHeader[i].nchan] = 0;
			      }
			  }
			else
			  inFile->beam[ibeam].bandData[i].astro_data.flag[j+idump*inFile->beam[ibeam].bandHeader[i].nchan] = 0;
		      }
		  }
	      }
	    recalc=2;
	  }
      }
    else if (key=='h')
      help();
    else if (key=='H')
      {
	float lowX,highX;
	float *plotX,*plotY;
	int *plotI,*plotJ;
	int nPlot=0;
	int nFit = 2;
	float new_miny = miny;
	float new_maxy = maxy;
	needSave=1;
	printf("Polarisation selection = %d\n",selectPol);

	plotX = (float *)malloc(sizeof(float)*max_nchan);
	plotY = (float *)malloc(sizeof(float)*max_nchan);
	plotI = (int *)malloc(sizeof(int)*max_nchan);
	plotJ = (int *)malloc(sizeof(int)*max_nchan);

	for (i=0;i<inFile->beam[ibeam].nBand;i++)
	  {
	    for (j=0;j<inFile->beam[ibeam].bandHeader[i].nchan;j++)
	      {
		// FIX ME: using [0] for frequency dump
		if (inFile->beam[ibeam].bandData[i].astro_data.freq[j] >= minx && inFile->beam[ibeam].bandData[i].astro_data.freq[j] <= maxx)
		  {
		    // FIX ME: using [0] for frequency dump
		    plotX[nPlot] = inFile->beam[ibeam].bandData[i].astro_data.freq[j];
		    if (selectPol == 1) {plotY[nPlot] = inFile->beam[ibeam].bandData[i].astro_data.pol1[j+idump*nchan];}
		    else if (selectPol == 2) {plotY[nPlot] = inFile->beam[ibeam].bandData[i].astro_data.pol2[j+idump*nchan];}
		    plotI[nPlot] = i;
		    plotJ[nPlot] = j;
		    nPlot++;
		  }
	      }
	  }
	do {
	  cpgenv(minx,maxx,new_miny,new_maxy,0,1);
	  sprintf(title,"Band: %s (spectral dump %d/%d)",inFile->beam[ibeam].bandHeader[iband].label,idump,ndump-1);
	  cpglab("Frequency (MHz)","Signal strength (arbitrary)",title);
	  
	  {
	    int i0=0;
	    int drawIt=-1;
	    for (i=0;i<nPlot;i++)
	      {
		if (inFile->beam[ibeam].bandData[plotI[i]].astro_data.flag[plotJ[i]] != 1 && drawIt==-1)
		  {
		    drawIt=1;
		    i0=i;
		  }
		if (inFile->beam[ibeam].bandData[plotI[i]].astro_data.flag[plotJ[i]] == 1 && drawIt==1)
		  {
		    cpgsci(1); cpgline(i-1-i0,plotX+i0,plotY+i0); cpgsci(1);
		    drawIt=-1;
		    
		  }	      
	      }
	    if (drawIt==1)
	      {
		cpgsci(1); cpgline(i-1-i0,plotX+i0,plotY+i0); cpgsci(1);
	      }
	  }
	  
       
	  cpgband(5,0,mx,my,&mx2,&my2,&key);
	  if (key=='A') // Mouse click
	    {
	      for (i=0;i<nPlot;i++)
		{
		  if (plotY[i] > my2)
		    {
		      if (zapAllDumps==1)
			{
			  for (k=0;k<ndump;k++)
			    inFile->beam[ibeam].bandData[plotI[i]].astro_data.flag[plotJ[i]+k*inFile->beam[ibeam].bandHeader[plotI[i]].nchan] = 1;
			}
		      else
			inFile->beam[ibeam].bandData[plotI[i]].astro_data.flag[plotJ[i]+idump*inFile->beam[ibeam].bandHeader[plotI[i]].nchan] = 1;   
		    }
		}
	    }
	  else if (key=='1') nFit=1;
	  else if (key=='2') nFit=2;
	  else if (key=='3') nFit=3;
	  else if (key=='4') nFit=4;
	  else if (key=='5') nFit=5;
	  else if (key=='f')  // Fit a polynomial
	    {
	      printf("Nfit = %d\n",nFit);
	      TKremovePoly_f(plotX,plotY,nPlot,nFit);
	      new_miny = 1e30; new_maxy = -1e-30;
	      for (i=0;i<nPlot;i++)
		{
		  if (inFile->beam[ibeam].bandData[plotI[i]].astro_data.flag[plotJ[i]+idump*inFile->beam[ibeam].bandHeader[plotI[i]].nchan] != 1)
		    {
		      if (new_miny > plotY[i]) new_miny = plotY[i];
		      if (new_maxy < plotY[i]) new_maxy = plotY[i];
		    }
		}
	      printf("New min/max = %g %g\n",new_miny,new_maxy);
	    }
	} while (key!='q');
	key=' ';
	recalc=3;
	free(plotX); free(plotY); free(plotI); free(plotJ);
      }
  } while (key != 'q');

  if (needSave==1)
    {
      char yesno[128];
      printf("You have unsaved changes\n");
      printf("save (y/n) ");
      scanf("%s",yesno);
      if (strcmp(yesno,"y")==0)
	{
	  printf("Saving ...\n");
	  saveFile(inFile,ibeam);
	  printf("Finished saving\n");
	}
    }
  free(px);
  free(py1);
  free(py2);
  free(pf1);
  free(pf2);
}

void saveFile(sdhdf_fileStruct *inFile,int ibeam)
{
  char oname[MAX_STRLEN];
  char flagName[MAX_STRLEN];
  sdhdf_fileStruct *outFile;
  int i,j;
  hsize_t dims[1];
  hid_t dset_id,dataspace_id;
  herr_t status;
  int *outFlags;

  // Should check if already .flag extension
  //
  sdhdf_formOutputFilename(inFile->fname,"flag",oname);
  printf("Saving to %s, please wait ... \n",oname);

  
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(outFile);
  sdhdf_openFile(oname,outFile,3);

  sdhdf_copyRemainder(inFile,outFile,0);
  // Now add the flag table
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      sdhdf_writeFlags(outFile,ibeam,i,inFile->beam[ibeam].bandData[i].astro_data.flag,inFile->beam[ibeam].bandHeader[i].nchan,inFile->beam[ibeam].bandHeader[i].ndump,inFile->beamHeader[ibeam].label,inFile->beam[ibeam].bandHeader[i].label);
    }


  sdhdf_closeFile(outFile);

  free(outFile);
  printf("Save completed\n");
}

void autoZapTransmitters(sdhdf_fileStruct *inFile,int ibeam,int zapAll)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;
  int zapAllTransmitters = 0;
  //
  // Fixed mobile transmission towers
  // We see persistent emission from these at all telescope pointing angles
  //
  f0[nT] = 758; f1[nT++] = 768;       //  Optus  
  f0[nT] = 768; f1[nT++] = 788;       //  Telstra
  f0[nT] = 869.95; f1[nT++] = 875.05; //  Vodafone
  f0[nT] = 875.05; f1[nT++] = 889.95; //  Telstra
  f0[nT] = 915; f1[nT++] = 928;       //  Agricultural 920 MHz RFI
  f0[nT] = 943.4; f1[nT++] = 951.8;   //  Optus
  f0[nT] = 953.7; f1[nT++] = 958.7;   //  Vodafone
  f0[nT] = 1081.7; f1[nT++] = 1086.7; // Aliased signal
  f0[nT] = 1720; f1[nT++] = 1736.2;   // Aliased signal
  f0[nT] = 1805; f1[nT++] = 1825;     //  Telstra
  f0[nT] = 1845; f1[nT++] = 1865;     //  Optus
  f0[nT] = 1973.7; f1[nT++] = 1991.0; // Aliased signal
  f0[nT] = 2110; f1[nT++] = 2120;     //  Vodafone
  f0[nT] = 2140; f1[nT++] = 2145;     //  Optus
  f0[nT] = 2145; f1[nT++] = 2150;     //  Optus
  f0[nT] = 2164.9; f1[nT++] = 2170.1; //  Vodafone
  f0[nT] = 2302.05; f1[nT++] = 2321.95; //  NBN
  f0[nT] = 2322.05; f1[nT++] = 2341.95; //  NBN
  f0[nT] = 2342.05; f1[nT++] = 2361.95; //  NBN
  f0[nT] = 2362.05; f1[nT++] = 2381.95; //  NBN
  f0[nT] = 2487.00; f1[nT++] = 2496.00; //  NBN alias
  f0[nT] = 2670; f1[nT++] = 2690;     //  Optus
  f0[nT] = 3445.05; f1[nT++] = 3464.95; //  NBN
  f0[nT] = 3550.05; f1[nT++] = 3569.95; //  NBN

  if (zapAllTransmitters==1)
    {
      f0[nT] = 804.4; f1[nT++] = 804.6;               //  NSW Police Force
      f0[nT] = 2152.5-2.5; f1[nT++] = 2152.5+2.5;     //  Telstra - not always on
      f0[nT] = 847.8-0.2; f1[nT++] = 847.8+0.2;     //  Radio broadcast Parkes
      f0[nT] = 849.5-0.1; f1[nT++] = 849.5+0.1;               //  NSW Police Force
      f0[nT] = 848.6-0.230/2.; f1[nT++] = 848.6+0.230/2.; //  Radio broadcast Mount Coonambro
      
      // Note see Licence number 1927906/1 in the ACMA database
      f0[nT] = 2127.5-2.5; f1[nT++] = 2127.5+2.5; // Parkes: "Station open to official correspondence exclusively"

      f0[nT] = 3575; f1[nT++] = 3640;               //  Telstra from Orange or Dubbo
    }  
  printf("Number of fixed transmitters being removed = %d\n",nT);
  flagChannels(inFile,ibeam,f0,f1,nT,zapAll);
}

void autoZapDigitisers(sdhdf_fileStruct *inFile,int ibeam,int zapAllSub)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;
  int zapAll=1;
  
  //
  // Digitiser-related signals that we always see
  //
  f0[nT] = 1023; f1[nT++] = 1025;
  f0[nT] = 1919.9; f1[nT++] = 1920.1;
  f0[nT] = 3071.9; f1[nT++] = 3072.05;       
  
  printf("Number of digitiser-related signals being removed = %d\n",nT);

  flagChannels(inFile,ibeam,f0,f1,nT,zapAllSub);
}


void autoZapUnexplained(sdhdf_fileStruct *inFile,int ibeam,int zapAllSub)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;
  int zapAll=1;
  
  // We do not know what is causing this
  f0[nT] = 824.95; f1[nT++] = 825.05;

  f0[nT] = 1225; f1[nT++] = 1230; // This is potentially interesting - should look more into this one
  f0[nT] = 1399.9; f1[nT++] = 1400.15;
  f0[nT] = 1498.0; f1[nT++] = 1499.0;
  f0[nT] = 1499.8; f1[nT++] = 1500.2;
  f0[nT] = 1880; f1[nT++] = 1904;
  f0[nT] = 2032.0; f1[nT++] = 2033.0;
  f0[nT] = 2063.0; f1[nT++] = 2064.0;
  f0[nT] = 2077.0; f1[nT++] = 2078.0;
  f0[nT] = 2079.0; f1[nT++] = 2080.0;
  f0[nT] = 2093.0; f1[nT++] = 2094.0;
  f0[nT] = 2160.0; f1[nT++] = 2161;
  f0[nT] = 2191.0; f1[nT++] = 2192;
  f0[nT] = 2205.0; f1[nT++] = 2206.0;
  f0[nT] = 2207.0; f1[nT++] = 2208.0;
  f0[nT] = 2221.0; f1[nT++] = 2222.0;
  f0[nT] = 2226.3; f1[nT++] = 2226.7;
  printf("Number of unexplained signals being removed = %d\n",nT);

  flagChannels(inFile,ibeam,f0,f1,nT,zapAllSub);
}

void autoZapWiFi(sdhdf_fileStruct *inFile,int ibeam,int zapAllSub)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;
  int zapAll=1;
  
  //
  // Digitiser-related signals that we always see
  //
  if (zapAll==1)
    {
      f0[nT] = 2401; f1[nT++] = 2483; // Entire band
    }
  printf("Number of WiFi signals being removed = %d\n",nT);
  flagChannels(inFile,ibeam,f0,f1,nT,zapAllSub);
}

void autoZapSatellites(sdhdf_fileStruct *inFile,int ibeam,int zapAllSub)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;
  int zapAll=1;


    //
  // Satellites that we always see
  //
  f0[nT] = 1618; f1[nT++] = 1626.5; // Iridium
  if (zapAll==1)
    {
      f0[nT] = 1164; f1[nT++] = 1189;
      f0[nT] = 1189; f1[nT++] = 1214;
      f0[nT] = 1240; f1[nT++] = 1260;
      f0[nT] = 1260; f1[nT++] = 1300;
      f0[nT] = 1525; f1[nT++] = 1646.5;  // Inmarsat - this is too wide
    }
  
  printf("Number of satellite signals being removed = %d\n",nT);

  flagChannels(inFile,ibeam,f0,f1,nT,zapAllSub);
}

void autoZapHandsets(sdhdf_fileStruct *inFile,int ibeam)
{
  int maxTransmitters = 512;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  float cleanF0 = 746; // SHOULD READ FROM AN AUXILARLY FILE -- MUST BE WITHIN 1 BAND
  float cleanF1 = 748;
  float f0v,f1v;
  int nT=0;
  int i,j,k;
  int nzap=0;
  int zapAll=1;
  int bandNum=-1;
  int n =0;
  int n1=0;
  int nchan;
  
  float s1,s2;
  float s1_2,s2_2;
  float t1,t2,t1_2,t2_2;
  float *pwr_clean,pwr_test;
  float mean1,mean2,sdev1,sdev2;
  float meanT1,meanT2,sdevT1,sdevT2;
  float sigma1,sigma2,psigma1,psigma2;
  float peak1,peak2;

  pwr_clean = (float *)malloc(sizeof(float)*inFile->beam[ibeam].bandHeader[bandNum].ndump);

  
  nT=0;
  for (i=0;i<150;i++)
    {
      f0[nT] = 703+0.2*i; f1[nT] = 703+0.2*(i+1); 
      nT++;
    }  
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      printf("Checking %g %g\n",inFile->beam[ibeam].bandHeader[i].f0,inFile->beam[ibeam].bandHeader[i].f1);
      if (inFile->beam[ibeam].bandHeader[i].f0  < cleanF0 && inFile->beam[ibeam].bandHeader[i].f1 > cleanF1)
	{
	  bandNum = i;
	  break;
	}
    }

  //      nchan = 
  if (bandNum == -1)
    printf("WARNING: Cannot find the clean band. No flagging applied\n");
  else
    {
      printf("In band %d\n",bandNum);
      nchan = inFile->beam[ibeam].bandHeader[bandNum].nchan;
      n=0;
      s1 = s2 = s1_2 = s2_2 = 0;
      for (j=0;j<inFile->beam[ibeam].bandHeader[bandNum].ndump;j++)
	{
	  n1=0;
	  pwr_clean[j]=0;
	  for (i=0;i<nchan;i++)
	    {
	      // FIX ME: using [0] for frequency dump
	      if (inFile->beam[ibeam].bandData[bandNum].astro_data.freq[i] > cleanF0 && inFile->beam[ibeam].bandData[bandNum].astro_data.freq[i] < cleanF1)
		{
		  s1 += inFile->beam[ibeam].bandData[bandNum].astro_data.pol1[i+j*nchan];
		  s1_2 += pow(inFile->beam[ibeam].bandData[bandNum].astro_data.pol1[i+j*nchan],2);

		  s2 += inFile->beam[ibeam].bandData[bandNum].astro_data.pol2[i+j*nchan];
		  s2_2 += pow(inFile->beam[ibeam].bandData[bandNum].astro_data.pol2[i+j*nchan],2);
		  pwr_clean[j] += inFile->beam[ibeam].bandData[bandNum].astro_data.pol1[i+j*nchan]+inFile->beam[ibeam].bandData[bandNum].astro_data.pol2[i+j*nchan];
		  n++;
		  n1++;
		}
	    }
	  pwr_clean[j] /= (float)n1;
	}


      mean1 = s1/(float)n;
      mean2 = s2/(float)n;

      sdev1 = sqrt(1/(float)n*s1_2 - pow(1.0/(float)n*s1,2));
      sdev2 = sqrt(1/(float)n*s2_2 - pow(1.0/(float)n*s2,2));
      
      printf("In clean band. Mean values are %g and %g, sdev = %g and %g, n = %d, band = %d\n",mean1,mean2,sdev1,sdev2,n,bandNum);      
    }
  
  // Obtain signal in different bands
  nchan = inFile->beam[ibeam].bandHeader[bandNum].nchan;
  for (j=0;j<inFile->beam[ibeam].bandHeader[bandNum].ndump;j++)
  //  for (j=0;j<5;j++)
    {
      printf("Processing spectral dump %d/%d\n",j,inFile->beam[ibeam].bandHeader[bandNum].ndump);
      for (k=0;k<nT;k++)
	{
	  t1=t2=t1_2=t2_2=0;
	  peak1=peak2=0.0;
	  n=0;
	  for (i=0;i<nchan;i++)
	    {
	      // Could significantly speed this up by simply recording the start and end channels
		      // FIX ME: using [0] for frequency dump
	      if (inFile->beam[ibeam].bandData[bandNum].astro_data.freq[i] > f0[k] && inFile->beam[ibeam].bandData[bandNum].astro_data.freq[i] < f1[k])
		{
		  if (inFile->beam[ibeam].bandData[bandNum].astro_data.pol1[i+j*nchan] > peak1)
		    peak1 = inFile->beam[ibeam].bandData[bandNum].astro_data.pol1[i+j*nchan];

		  t1 += inFile->beam[ibeam].bandData[bandNum].astro_data.pol1[i+j*nchan];
		  t1_2 += pow(inFile->beam[ibeam].bandData[bandNum].astro_data.pol1[i+j*nchan],2);

		  if (inFile->beam[ibeam].bandData[bandNum].astro_data.pol2[i+j*nchan] > peak2)
		    peak2 = inFile->beam[ibeam].bandData[bandNum].astro_data.pol2[i+j*nchan];

		  t2 += inFile->beam[ibeam].bandData[bandNum].astro_data.pol2[i+j*nchan];
		  t2_2 += pow(inFile->beam[ibeam].bandData[bandNum].astro_data.pol2[i+j*nchan],2);
		  n++;

		}
	    }
	  pwr_test = (t1+t2)/(float)n;
	  meanT1 = t1/(float)n;
	  meanT2 = t2/(float)n;
	  
	  sdevT1 = sqrt(1/(float)n*t1_2 - pow(1.0/(float)n*t1,2));
	  sdevT2 = sqrt(1/(float)n*t2_2 - pow(1.0/(float)n*t2,2));

	  sigma1 = (meanT1-mean1)/sdev1;
	  sigma2 = (meanT2-mean2)/sdev2;
	  psigma1 = (peak1-mean1)/sdev1;
	  psigma2 = (peak2-mean2)/sdev2;
	  //	  if (psigma1 > 3 || psigma2 > 3)  // NOTE HARDCODING THRESHOLD HERE AND ONLY USING THE PEAK
	  printf("Checking %d %d %g %g [%g %g]\n",j,k,pwr_clean[j],pwr_test,f0[k],f1[k]);
	  if (pwr_test > pwr_clean[j]*1.5)
	    {
	      //	      printf("REMOVE Check %d %d %g %g %g %g sigma = %g %g, peak sigma = %g %g\n",j,k,meanT1,meanT2,sdevT1,sdevT2,sigma1,sigma2,psigma1,psigma2);
	      for (i=0;i<nchan;i++)
		{
		  // FIX ME: using [0] for frequency dump
		  if (inFile->beam[ibeam].bandData[bandNum].astro_data.freq[i] > f0[k] && inFile->beam[ibeam].bandData[bandNum].astro_data.freq[i] < f1[k])
		    inFile->beam[ibeam].bandData[bandNum].astro_data.flag[i+nchan*j] = 1;
		}  
	    }
	}
    }
  

  
  /*
  
  if (zapAll==1)
    {
      f0[nT] = 703;    f1[nT++] = 713;    // 4G Optus
      f0[nT] = 704.5;  f1[nT++] = 708;    // Alias
      f0[nT] = 713;    f1[nT++] = 733;
      f0[nT] = 825.1;  f1[nT++] = 829.9; // Vodafone 4G
      f0[nT] = 830.05; f1[nT++] = 844.95; // Telstra 3G
      f0[nT] = 847.6;  f1[nT++] = 848.0;
      f0[nT] = 898.4;  f1[nT++] = 906.4;  // Optus 3G
      f0[nT] = 906.8;  f1[nT++] = 915.0;  // Vodafone 3G
      f0[nT] = 953;    f1[nT++] = 960.1;  // Alias
      f0[nT] = 1710;   f1[nT++] = 1725;  // 4G Telstra
      f0[nT] = 1745;   f1[nT++] = 1755;  // 4G Optus
      f0[nT] = 2550;   f1[nT++] = 2570;  // 4G Optus
      
    }
  printf("Number of Handset signals being removed = %d\n",nT);

  flagChannels(inFile,ibeam,f0,f1,nT,zapAllSub);
  */

  free(pwr_clean);
}

void autoZapAircraft(sdhdf_fileStruct *inFile,int ibeam,int zapAllSub)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;
  int zapAll=1;

    //
  // Aircraft signals that we always see
  //
  f0[nT] = 1017; f1[nT++] = 1019;       // Parkes ground response
  f0[nT] = 1029; f1[nT++] = 1031;       // Parkes ground response

  if (zapAll==1)
    {
      f0[nT] = 1026.8; f1[nT++] = 1027.2;       // Unknown DME signal
      f0[nT] = 1027.8; f1[nT++] = 1028.2;       // Unknown DME signal
      f0[nT] = 1032.8; f1[nT++] = 1033.2;       // Unknown DME signal
      f0[nT] = 1040.8; f1[nT++] = 1041.2;       // Strong, unknown DME signal
      f0[nT] = 1061.8; f1[nT++] = 1062.2;       // Sydney DME
      f0[nT] = 1067.8; f1[nT++] = 1068.2;       // Richmond DME
      f0[nT] = 1071.8; f1[nT++] = 1072.2;       // Wagga Wagga DME
      f0[nT] = 1079.2; f1[nT++] = 1080.2;       // Unknown DME
      f0[nT] = 1080.8; f1[nT++] = 1081.2;       // Parkes DME
      f0[nT] = 1081.8; f1[nT++] = 1082.2;       // Sydney RW07 DME
      f0[nT] = 1103.8; f1[nT++] = 1104.2;       // Unknown DME
      f0[nT] = 1102.8; f1[nT++] = 1103.2;       // Unknown DME
      f0[nT] = 1120.8; f1[nT++] = 1121.2;       // Wagga Wagga DME
      f0[nT] = 1134.8; f1[nT++] = 1135.2;       // Nowra DME
      f0[nT] = 1137.6; f1[nT++] = 1138.4;       // Canberra DME
      f0[nT] = 1149.8; f1[nT++] = 1150.2;       // CHECK IF THIS IS A DME
      f0[nT] = 1150.8; f1[nT++] = 1151.2;       // CHECK IF THIS IS A DME
      
    }
  printf("Number of aircraft signals being removed = %d\n",nT);


  flagChannels(inFile,ibeam,f0,f1,nT,zapAllSub);
}


void flagChannels(sdhdf_fileStruct *inFile,int ibeam,float *f0,float *f1,int nT,int zapAll)
{
  int i,j,k,t;
  int nchan;
  float f;
  int zap=0;
  
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
  	  nchan = inFile->beam[ibeam].bandHeader[i].nchan; 
	  for (j=0;j<nchan;j++)
	    {
	      // FIX ME: using [0] for frequency dump
	      f = inFile->beam[ibeam].bandData[i].astro_data.freq[j];
	      zap=0;
	      for (k=0;k<nT;k++)
		{
		  if (f >= f0[k] && f <= f1[k])
		    {
		      zap=1;
		      break;
		    }
		}
	      if (zap==1)
		{
		  if (zapAll==0)
		    inFile->beam[ibeam].bandData[i].astro_data.flag[j]=1;
		  else
		    {
		      for (k=0;k<inFile->beam[ibeam].bandHeader[i].ndump;k++)
			inFile->beam[ibeam].bandData[i].astro_data.flag[j+k*inFile->beam[ibeam].bandHeader[i].nchan] = 1;
		    }
		  
	    }
	}
    }

  printf("Save completed\n");
}


void zapPersistent(sdhdf_fileStruct *inFile,int ibeam,int zapAll)
{
  sdhdf_rfi rfi[MAX_RFI];
  int nRFI;
  float f0[MAX_RFI],f1[MAX_RFI];
  int i;
  int n=0;
  float freq;
  double mjd;
  
  if (strcmp(inFile->primary[0].telescope,"Parkes")!=0)
    printf("WARNING: ONLY IMPLEMENTED PARKES OBSERVATORY\n"); // FIX ME	      
  sdhdf_loadPersistentRFI(rfi,&nRFI,MAX_RFI,"parkes");

  for (i=0;i<nRFI;i++)
    {
      if (rfi[i].type == 1) 
	{
	  mjd = inFile->beam[ibeam].bandData[0].astro_obsHeader[0].mjd; // FIX 0 FOR BAND AND DUMP HERE
	  if (mjd >= rfi[i].mjd0 && mjd <= rfi[i].mjd1)
	    {
	      f0[n] = rfi[i].f0;
	      f1[n] = rfi[i].f1;
	      n++;
	    }
	}
      else
	{
	  printf("WARNING: unknown RFI type in sdhdf_flag\n");
	}
    }
  
  flagChannels(inFile,ibeam,f0,f1,n,zapAll);  
}
