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
#include <cpgplot.h>
#include "TKfit.h"

#define MAX_SUBBANDS 1024
#define MAX_DUMPS 4096

void saveFile(sdhdf_fileStruct *inFile);
void drawVline(float t0, float t1,char *str);
int loadBands(char *bandFile,float *sbF0,float *sbF1,int *col,int *style,char **bandHeader);
float calcT(char *local,char *local0,int min);
void drawLabels(char *labelFile,char *t0,float miny,float maxy,int min);
void sort2d(int n,float *a,float *b);


void help()
{
  printf("sdhdf_rfiTime: routine to plot RFI as a function of time\n");
  printf("\n");
  printf("Command line arguments\n\n");
  printf("-b <filename>           Specify bands from file given by filename\n");
  printf("-div1                   Divide all values by first entry\n");
  printf("-gap <val>              Define gap size\n");
  printf("-h                      This help\n");
  printf("-ignoreAz <val1> <val2> Ignore azimuth ranges between val1 and val2\n");
  printf("-l <filename>           Load labels from file\n");
  printf("-label <x> <y> <col> <str> Place a label at x,y with colour col and the label is str\n");
  printf("-log                    Plot on logarithmic scale\n");
  printf("-miny <val>             Define minimum y-value for the plot\n");
  printf("-maxy <val>             Define maximum y-value for the plot\n");
  printf("-noEl                   Do not plot elevation\n");
  printf("-plotAz                 Plot RFI as a function of azimuth\n");
  printf("-showFiles              Write file information\n");
  printf("-sm <val>               Smoothing integer value\n");
  printf("-t0 <string>            Define local start time\n");
  printf("\n\n");
  printf("Filenames are given on the command line\n");
}

int main(int argc,char *argv[])
{
  int i,j,k,ii,jj,kk,b;
  sdhdf_fileStruct *inFile;
  int nchan,idump,nband;
  int totSize,chanPos;
  int ndump;
  //  spectralDumpStruct spectrum;
  int npol=4;
  int ibeam=0;
  int iband=0;
  int ichan;
  int setnband=-1;
  float tdump;
  int   np;
  int sdump=0;
  int sdump0=0;
  char local0[MAX_STRLEN]="00:00:00";
  
  double freq;
  float *timeVal;
  int min=1;
  float *elVal,*azVal;
  float **signalVal;

  FILE *fin;
  FILE *fout;
  FILE *fout2;
  char outName[1024];
  char header[1024];  
  char bandFile[1024] = "NULL";
  int  nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];

  int nsb=0;
  float *sbF0,*sbF1;
  int col[MAX_SUBBANDS];
  int style[MAX_SUBBANDS];
  char **bandHeader;
  int log=0;
  char labelFile[MAX_STRLEN] = "NULL";
  char xlabel[MAX_STRLEN];
  float px[MAX_DUMPS],py[MAX_DUMPS],ps[MAX_DUMPS];
  int nplot;
  float miny = -1,maxy = -1;
  float gap = 1800; // Gap length in seconds
  int i0=0;
  int sm=1;
  int nc;
  int div1=0;
  int sub0=0;
  
  int nLabel=0;
  char labelStr[8][MAX_STRLEN];
  float labelX[8];
  float labelY[8];
  int   labelCol[8];
  float fx[2],fy[2];

  float minAz = -5;
  float maxAz = 365;
  
  int plotType=1;
  int showFiles=0;
  int plotEl=1;
  
  int nIgnoreAz=0;
  float ignoreAz1[1024];
  float ignoreAz2[1024];
  int ignoreIt;
  
  timeVal = (float *)malloc(sizeof(float)*MAX_DUMPS);
  elVal   = (float *)malloc(sizeof(float)*MAX_DUMPS);
  azVal   = (float *)malloc(sizeof(float)*MAX_DUMPS);
  signalVal = (float **)malloc(sizeof(float *)*MAX_SUBBANDS);
  for (i=0;i<MAX_SUBBANDS;i++)
    signalVal[i] = (float *)malloc(sizeof(float)*MAX_DUMPS);

  
  // Allocate memory to record the required bands
  sbF0       = (float *)malloc(sizeof(float)*MAX_SUBBANDS);
  sbF1       = (float *)malloc(sizeof(float)*MAX_SUBBANDS);
  bandHeader = (char **)malloc(sizeof(char *)*MAX_SUBBANDS);
  for (i=0;i<MAX_SUBBANDS;i++)
    bandHeader[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);

  for (i=1;i<argc;i++)
    {
      printf("Checking %s\n",argv[i]);
      if (strcmp(argv[i],"-b")==0)
	strcpy(bandFile,argv[++i]);
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-log")==0)
	log=1;
      else if (strcmp(argv[i],"-t0")==0)
	strcpy(local0,argv[++i]);
      else if (strcmp(argv[i],"-plotAz")==0)
	plotType=2;
      else if (strcmp(argv[i],"-l")==0)
	strcpy(labelFile,argv[++i]);
      else if (strcasecmp(argv[i],"-noEl")==0)
	plotEl=0;
      else if (strcmp(argv[i],"-gap")==0)
	sscanf(argv[++i],"%f",&gap);
      else if (strcmp(argv[i],"-sm")==0)
	sscanf(argv[++i],"%d",&sm);
      else if (strcmp(argv[i],"-miny")==0)
	sscanf(argv[++i],"%f",&miny);
      else if (strcmp(argv[i],"-maxy")==0)
	sscanf(argv[++i],"%f",&maxy);
      else if (strcmp(argv[i],"-label")==0)
	{
	  sscanf(argv[++i],"%f",&labelX[nLabel]);
	  sscanf(argv[++i],"%f",&labelY[nLabel]);
	  sscanf(argv[++i],"%d",&labelCol[nLabel]);
	  strcpy(labelStr[nLabel++],argv[++i]);
	}
      else if (strcmp(argv[i],"-div1")==0)
	div1=1;
      else if (strcmp(argv[i],"-sub0")==0)
	sub0=1;
      else if (strcmp(argv[i],"-showFiles")==0)	
	showFiles=1;
      else if (strcmp(argv[i],"-ignoreAz")==0)
	{
	  sscanf(argv[++i],"%f",&ignoreAz1[nIgnoreAz]);
	  sscanf(argv[++i],"%f",&ignoreAz2[nIgnoreAz++]);
	}
      else
	strcpy(fname[nFiles++],argv[i]);
    }

  if (min == 1)
    {
      gap/=60.0;
    }
  printf("Will process bandfile >%s< and %d input files\n",bandFile,nFiles);
  nsb = loadBands(bandFile,sbF0,sbF1,col,style,bandHeader);
  printf("Loading %d sub-bands\n",nsb);
    
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdump=sdump0=0;
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_loadMetaData(inFile);

      // Find out which UWL bands we will need to load and process
      for (j=0;j<nsb;j++)
	{
	  printf("Processing SB %d (%s)\n",j,inFile->fname);
	  for (k=0;k<inFile->beam[ibeam].nBand;k++)
	    {
	      if (sbF0[j] < inFile->beam[ibeam].bandHeader[k].f1 &&
		  sbF1[j] > inFile->beam[ibeam].bandHeader[k].f0)
		{	      
		  printf("  Will use UWL band: %d\n",k);
		  
		  // *** Check that this doesn't reload if it doesn't need to
		  sdhdf_loadBandData(inFile,ibeam,k,1);
		  printf("      Loaded data\n");
		}
	    }
	  
	  // What to do if number of dumps is not the same in different bands?
	  ndump = inFile->beam[ibeam].bandHeader[0].ndump;
	  sdump = sdump0;
	  printf("    Number of spectral dumps = %d (sdump = %d)\n",ndump,sdump);
	  for (jj=0;jj<ndump;jj++)
	    {
	      if (j==0)
		{
		  //		  timeVal[sdump] = calcT(inFile->beam[ibeam].bandData[0].astro_obsHeader[jj].aest,aest0,min);
		  timeVal[sdump] = calcT(inFile->beam[ibeam].bandData[0].astro_obsHeader[jj].utc,local0,min);
		  azVal[sdump] = inFile->beam[ibeam].bandData[0].astro_obsHeader[jj].az;
		  elVal[sdump] = inFile->beam[ibeam].bandData[0].astro_obsHeader[jj].el;
		}
	      
	      ignoreIt=0;
	      for (kk=0;kk<nIgnoreAz;kk++)
		{
		  if ( inFile->beam[ibeam].bandData[0].astro_obsHeader[jj].az >= ignoreAz1[kk] && inFile->beam[ibeam].bandData[0].astro_obsHeader[jj].az < ignoreAz2[kk])
		    {
		      ignoreIt=1;
		      break;
		    }
		}
	      if (ignoreIt==0)
		{
		  signalVal[j][sdump] = 0.0;
		  np=0;
		  
		  for (k=0;k<inFile->beam[ibeam].nBand;k++)
		    {
		      if (sbF0[j] < inFile->beam[ibeam].bandHeader[k].f1 &&
			  sbF1[j] > inFile->beam[ibeam].bandHeader[k].f0)
			{
			  nchan = inFile->beam[ibeam].bandHeader[k].nchan;
			  for (ii=0;ii<nchan;ii++)
			    {
			      // FIX ME: using [0] for frequency dump
			      freq = inFile->beam[ibeam].bandData[k].astro_data.freq[ii];
			      if (freq >= sbF0[j] && freq <= sbF1[j])
				{
				  signalVal[j][sdump] += (inFile->beam[ibeam].bandData[k].astro_data.pol1[ii+jj*nchan]+inFile->beam[ibeam].bandData[k].astro_data.pol2[ii+jj*nchan]);
				  //			      printf("Complete set\n");
				  np++;
				}
			    }
			}
		    }
		  signalVal[j][sdump]/=(double)np;
		  sdump++;
		  if (sdump > MAX_DUMPS)
		    {
		      printf("ERROR: must increase MAX_DUMPS in sdhdf_rfiTime.c\n");
		      exit(1);
		    }

		}
	    }
	}
      sdump0+=(sdump-sdump0);

      sdhdf_closeFile(inFile);
    }

  sdump = sdump0;
  printf("Have loaded %d values\n",sdump);
  fout = fopen("output.dat","w");
  for (i=0;i<sdump;i++)
    {
      fprintf(fout,"%g %g %g ",timeVal[i],azVal[i],elVal[i]);
      for (j=0;j<nsb;j++)
	fprintf(fout,"%g ",signalVal[j][i]);
      fprintf(fout,"\n");
    }
  fclose(fout);

  if (sub0==1)
    {
      double v0[nsb];
      float minVal,maxVal;
      for (j=0;j<nsb;j++)
	v0[j] = signalVal[j][0];
      
      for (j=0;j<nsb;j++)
	{
	  minVal = 1e99;
	  maxVal = -1e99;

	  for (i=0;i<sdump;i++)
	    {
	      if (minVal > signalVal[j][i])
		minVal = signalVal[j][i];
	      if (maxVal < signalVal[j][i])
		maxVal = signalVal[j][i];
	    }
	  printf("%d min/max = %g %g\n",j,minVal,maxVal);
	  for (i=0;i<sdump;i++)
	    signalVal[j][i]= (signalVal[j][i]-minVal)/(maxVal-minVal);

	  /*
	    if (j==0)
	    {
	    for (i=0;i<sdump;i++)
		printf("Have %d %g\n",i,signalVal[j][i]);
	    }
	  */
	}
    }

    if (div1==1)
    {
      double v0[nsb];
      for (j=0;j<nsb;j++)
	v0[j] = signalVal[j][0];
      
      for (i=0;i<sdump;i++)
	{
	  for (j=0;j<nsb;j++)
	    {
	      if (log==1)
		signalVal[j][i]/=v0[j];
	      else
		signalVal[j][i]-=v0[j];
	    }
	}
    }

  
  // Convert to log
  if (log==1)
    {
      for (i=0;i<sdump;i++)
	{
	  for (j=0;j<nsb;j++)
	    signalVal[j][i]=log10(signalVal[j][i]);
	}
    }

  // Smooth, plot, label, log and highlight regions
  printf("PlotType= %d\n",plotType);
  if (plotType==1)
    {
      printf("sdump = %d\n",sdump);
      cpgbeg(0,"plot.ps/cps",1,1);
      cpgslw(2);
      cpgscf(2);
      cpgsch(1.0);
      // Plot Az/El
      cpgsch(1.0);
      cpgsci(2);
      cpgsvp(0.10,0.95,0.53,0.71);
      cpgswin(timeVal[0],timeVal[sdump-1],0,360);

      
      cpgbox("ABC",0,0,"ABCTSN",0,0);
      cpglab("","Az (deg)","");

      for (i=50;i<360;i+=50)
	{
	  fx[0] = timeVal[0]; fx[1] = timeVal[sdump-1];
	  fy[0] = fy[1] = i;
	  cpgsci(1); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(2);
	}
      
      cpgpt(sdump,timeVal,azVal,22);
      
      cpgsci(4);
      cpgsvp(0.10,0.95,0.715,0.895);
      cpgswin(timeVal[0],timeVal[sdump-1],15,90);
      cpgbox("ABC",0,0,"ABCTSN",0,0);
      cpglab("","El (deg)","");

      for (i=10;i<90;i+=10)
	{
	  fx[0] = timeVal[0]; fx[1] = timeVal[sdump-1];
	  fy[0] = fy[1] = i;
	  cpgsci(1); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(4);
	}
      
      
      cpgpt(sdump,timeVal,elVal,22);
      cpgsch(1.4);
      cpgsci(1);
      
      cpgsvp(0.10,0.95,0.12,0.5);
      
      if (miny == maxy)
	{
	  miny = 1e99; maxy = -1e99;
	  for (j=0;j<nsb;j++)
	    {
	      nplot=0;
	      i0=0;
	      for (i=0;i<sdump;i++)
		{
		  if (signalVal[j][i] > maxy) maxy = signalVal[j][i];
		  if (signalVal[j][i] < miny) miny = signalVal[j][i];
		}
	    }
	}
      
      
      if (log==1)
	{
	  //      cpgenv(timeVal[0],timeVal[sdump-1],miny,maxy,0,20);
	  cpgswin(timeVal[0],timeVal[sdump-1],miny,maxy);
	  cpgbox("ABCTSN",0,0,"ABCTSNL",0,0);
	}
      else
	{
	  //      cpgenv(timeVal[0],timeVal[sdump-1],miny,maxy,0,1);
	  cpgswin(timeVal[0],timeVal[sdump-1],miny,maxy);
	  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
	}
      if (min==1)
	sprintf(xlabel,"Minutes since UTC %s",local0);
	//	sprintf(xlabel,"Minutes since AEST %s",aest0);
      else
	sprintf(xlabel,"Seconds since UTC %s",local0);
	//	sprintf(xlabel,"Seconds since AEST %s",aest0);
      cpglab(xlabel,"Relative signal strength","");
      
      for (j=0;j<nsb;j++)
	{
	  nplot=0;
	  i0=0;
	  for (i=0;i<sdump;i++)
	    {
	      if (i+1 < sdump && timeVal[i+1]-timeVal[i] < gap)
		{
		  px[nplot] = timeVal[i];
		  py[nplot] = signalVal[j][i];
		  nplot++;
		}
	      else
		{
		  printf("Splitting here\n");
		  if (i+1 >= sdump)
		    {
		      px[nplot] = timeVal[i];
		      py[nplot] = signalVal[j][i];
		      nplot++;		  
		    }
		  
		  
		  // Apply some smoothing
		  for (ii=0;ii<nplot;ii++)
		    {
		      ps[ii] = 0;
		      nc=0;
		      for (jj=ii;jj<ii+sm;jj++)
			{
			  if (jj < nplot)
			    {
			      ps[ii] += py[jj];
			      nc++;
			    }
			}
		      ps[ii] /= (float)nc;
		    }
		  
		  cpgslw(4);
		  cpgsci(col[j]); cpgsls(style[j]); cpgline(nplot,px,ps);
		  cpgslw(2);
		  cpgsci(1); cpgsls(1);
		  nplot=0;
		}
	    }


	  cpgsci(1); cpgsls(1);
	}
      if (strcmp(labelFile,"NULL")!=0)
	{
	  drawLabels(labelFile,local0,miny,maxy,min);
	}
      cpgend();
      exit(1);

    }
  else if (plotType==2)
    {
      float lineX[1024],lineY[1024];
      int nl=0;
      int nc;
      float mean;
      
      cpgbeg(0,"plot_az.ps/cps",1,1);
      cpgslw(2);
      cpgscf(2);
      cpgsch(1.0);

      if (plotEl==1)
	{
	  cpgsci(4);
	  cpgsvp(0.10,0.95,0.715,0.895);
	  cpgswin(minAz,maxAz,15,90);
	  cpgbox("BC",0,0,"BCTSN",0,0);
	  cpglab("","El (deg)","");
	  for (i=10;i<90;i+=10)
	    {
	      fx[0] = minAz; fx[1] = maxAz;
	      fy[0] = fy[1] = i;
	      cpgsci(1); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(4);
	    }
	}
      
      cpgpt(sdump,azVal,elVal,22);
      cpgsch(0.8);
      cpgsci(1);

      if (plotEl==1)
	cpgsvp(0.10,0.95,0.6,0.7);
      else
	cpgsvp(0.10,0.95,0.8,0.9);
      cpgswin(minAz,maxAz,0,1);
      cpgptxt(0,0.5,0,0.5,"N");
      cpgptxt(360,0.5,0,0.5,"N");
      cpgptxt(180,0.5,0,0.5,"S");
      cpgptxt(90,0.5,0,0.5,"E");
      cpgptxt(270,0.5,0,0.5,"W");


      cpgsch(1.4);
      if (plotEl==1)
	cpgsvp(0.10,0.95,0.12,0.6);
      else
	cpgsvp(0.10,0.95,0.12,0.8);
      if (miny == maxy)
	{
	  miny = 1e99; maxy = -1e99;
	  for (j=0;j<nsb;j++)
	    {
	      nplot=0;
	      i0=0;
	      for (i=0;i<sdump;i++)
		{
		  if (signalVal[j][i] > maxy) maxy = signalVal[j][i];
		  if (signalVal[j][i] < miny) miny = signalVal[j][i];
		}
	    }
	  miny = miny-(maxy-miny)*0.1;
	}
      
      printf("min/max in plot is %g/%g\n",miny,maxy);
      if (log==1)
	{
	  //      cpgenv(timeVal[0],timeVal[sdump-1],miny,maxy,0,20);
	  cpgswin(minAz,maxAz,miny,maxy);
	  cpgbox("BCTSN",0,0,"BCTSNL",0,0);
	}
      else
	{
	  //      cpgenv(timeVal[0],timeVal[sdump-1],miny,maxy,0,1);
	  cpgswin(minAz,maxAz,miny,maxy);
	  cpgbox("BCTSN",0,0,"BCTSN",0,0);
	}


      sprintf(xlabel,"Azimuthal angle (deg)");
      cpglab(xlabel,"Signal strength","");
      
      fx[0] = fx[1] = 90; fy[0] = miny; fy[1]=maxy;
      cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
      fx[0] = fx[1] = 180; fy[0] = miny; fy[1]=maxy;
      cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
      fx[0] = fx[1] = 270; fy[0] = miny; fy[1]=maxy;
      cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
      fx[0] = fx[1] = 360; fy[0] = miny; fy[1]=maxy;
      cpgsls(4); cpgline(2,fx,fy); cpgsci(1);

	  
      for (j=0;j<nsb;j++)
	{
	  
	  mean=0;
	  nl=nc=0;
	  lineX[nl] = 0;
	  lineY[nl] = 0;

	  for (i=0;i<sdump;i++)
	    {
	      if (i < sdump-1 && (i==0 || fabs(azVal[i]-azVal[i-1]) < gap))
		{
		  mean+=signalVal[j][i];
		  nc++;
		}
	      else
		{
		  if (nc > 0)
		    {
		      mean/=(float)nc;
		      lineX[nl] = azVal[i-1];
		      lineY[nl] = mean;
		      mean=0;
		      nc=0;
		      nl++;
		    }
		}
	    }
	  cpgsci(col[j]);
	  cpgsls(1);

	  //	  cpgpt(sdump,azVal,signalVal[j],21);

	  // Put in to get the high resolution plot working ** check this
	  nl = sdump;
	  for (i=0;i<nl;i++)
	    {
	      lineX[i] = azVal[i];
	      lineY[i] = signalVal[j][i];
	    }
	  sort2d(nl,lineX,lineY);
	  cpgline(nl,lineX,lineY);

	  // HERE **
	  //	  cpgpt(nl,lineX,lineY,4);

	  //	  for (i=0;i<nl;i++)
	  //	    printf("PLOTTING: %g %g\n",lineX[i],lineY[i]);
	  cpgsci(1);

	}
    }

  // Put headers and legend

  cpgsvp(0.10,0.95,0.895,0.995);
  cpgswin(0,1,0,1);
  
  // Command line labels:
  for (i=0;i<nLabel;i++)
    {
      cpgsci(labelCol[i]);
      fx[0] = labelX[i];
      fx[1] = labelX[i]+0.05;
      fy[0] = fy[1] = labelY[i];
      printf("Plotting at %g %g %g %g\n",fx[0],fy[0],fx[1],fy[1]);
      cpgline(2,fx,fy);
      cpgsch(0.8);
      cpgptxt(fx[1]+0.01,fy[0],0,0.0,labelStr[i]);
      cpgsch(1.4);
      cpgsci(1);
    }
  
  // Write file information
  if (showFiles==1)
  {
    float xp=0.35;
    float yp=0.94;
    cpgsch(0.6);
    for (i=0;i<nFiles;i++)
      {
	cpgptxt(xp,yp,0,0.0,fname[i]);
	yp -= 0.15;
	if (yp < 0.15) {yp = 0.94; xp += 0.20;}
      }
    cpgsch(1.4);
  }


  
  free(sbF0);
  free(sbF1);
  for (i=0;i<MAX_SUBBANDS;i++)
    free(bandHeader[i]);
  free(bandHeader);
  for (i=0;i<MAX_SUBBANDS;i++)
    free(signalVal[i]);
  free(signalVal);
  free(timeVal);
  free(azVal);
  free(elVal);
}

void drawLabels(char *labelFile,char *t0,float miny,float maxy,int min)
{
  FILE *fin;
  char line[1024];
  char tstart[1024],tend[1024];
  float tval1,tval2,yfrac;
  int col;
  int style;
  float fx[2],fy[2];
  float xch,ych;
  
  fin = fopen(labelFile,"r");
  while (!feof(fin))
    {
      if (fgets(line,1024,fin)!=NULL)
	{
	  line[strlen(line)-1]='\0';
	  printf("%s\n",line);
	  sscanf(line,"%s %s %d %d %f",tstart,tend,&col,&style,&yfrac);
	  tval1 =  calcT(tstart,t0,min);
	  tval2 =  calcT(tend,t0,min);

	  cpgsci(col);
	  fx[0] = tval1; fx[1] = tval2;
	  if (style==0)
	    {
	      fy[0] = miny; fy[1] = maxy;
	      cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
	    }
	  else
	    {
	      fy[0] = fy[1] = miny + (maxy-miny)*yfrac;
	      cpgline(2,fx,fy);
	      printf("Label = %s\n",line+27);
	      cpgsch(0.9);

	      cpgqcs(4,&xch,&ych);
	      cpgptxt((fx[0]+fx[1])/2.,fy[0]+ych/2,0,0.5,line+27);
	      cpgsch(0.5);
	      cpgsah(1,45,0.3);
	      cpgarro(fx[0],fy[0],fx[1],fy[1]);
	      cpgarro(fx[1],fy[1],fx[0],fy[0]);
	      cpgsch(1.4);
	      cpgsls(4);
	      fy[0] = miny; fy[1] = maxy;
	      fx[0] = fx[1] = tval1;
	      cpgslw(3);
	      cpgline(2,fx,fy);
	      fx[0] = fx[1] = tval2;
	      cpgline(2,fx,fy);
	      cpgsls(1);
	      cpgslw(2);
	      
	    }


	  cpgsci(1);
	}
      
    }
  fclose(fin);
}

void drawVline(float t0, float t1,char *str)
{
  float fx[2],fy[2];

  fx[0] = fx[1] = t1-t0;
  fy[0] = 6; fy[1] = 11;
  cpgsls(4);
  cpgsci(5); cpgline(2,fx,fy); cpgsci(1);
  cpgptxt(fx[0],6.5,0,0.5,str);


  //  cpgptxt(fx[0],8.1,0,0.5,str);
  cpgsls(1);
}

int loadBands(char *bandFile,float *sbF0,float *sbF1,int *col,int *style,char **bandHeader)
{
  int nsb = 0;
  FILE *fin;

  fin = fopen(bandFile,"r");
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f %d %d %s",&sbF0[nsb],&sbF1[nsb],&col[nsb],&style[nsb],bandHeader[nsb])==5)
	nsb++;
    }
  fclose(fin);

  
  return nsb;
}

float calcT(char *local,char *local0,int min)
{
  float hr0,min0,sec0,t0;
  float hr1,min1,sec1,t1;
  sscanf(local0,"%f:%f:%f",&hr0,&min0,&sec0);
  sscanf(local,"%f:%f:%f",&hr1,&min1,&sec1);

  t0 = hr0*60*60+min0*60+sec0;
  t1 = hr1*60*60+min1*60+sec1;

  if (min==1)
    return (t1-t0)/60.0;
  else
    return t1-t0;
}

void sort2d(int n,float *a,float *b)
{
  int i;
  int swap=0;
  float t1,t2;
  
  do {
    swap=0;
    for (i=0;i<n-1;i++)
      {
	if (a[i] > a[i+1])
	  {
	    t1 = a[i];
	    t2 = b[i];
	    a[i] = a[i+1]; a[i+1] = t1;
	    b[i] = b[i+1]; b[i+1] = t2;
	    swap=1;
	  }
      }
  } while (swap==1);
}
