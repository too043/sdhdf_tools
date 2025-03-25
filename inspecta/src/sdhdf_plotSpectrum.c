//  Copyright (C) 2019, 2020, 2021, 2022 George Hobbs

/*
 *    This file is part of INSPECTA. 
 * 
 *    INSPECTA is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    INSPECTA is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with INSPECTA.  If not, see <http://www.gnu.org/licenses/>. 
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
#include "inspecta.h"
#include "hdf5.h"
#include <cpgplot.h>

void plotSpectrum(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,double fref,char *yUnit,char *freqFrame,char *freqUnit,
		  char *grDev);
void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY,int drawLabel);
void drawRecombLine(float minX,float maxX,float minY,float maxY);

#define VNUM "v2.0"
#define MAX_REST_FREQUENCIES 1024

void help()
{
  printf("\nsdhdf_plotSpectrum: %s\n",VNUM);
  printf("sdhfProc version:   %s\n",SOFTWARE_VER);
  printf("Author:             George Hobbs\n");
  printf("Software to plot a spectrum from an SDHDF file interactively\n");

  printf("\nCommand line arguments:\n\n");
  printf("-h                  This help\n");
  printf("-f <filename>       SDHDF file corresponding to observation\n");
  printf("-beam <num>         Set beam number\n");
  printf("-g <string>         PGPLOT graphical device\n");
  printf("-sb <num>           Set the sub-band number (starting from 0)\n");
  printf("-sd <num>           Set the spectral dump number (starting from 0)\n");
  printf("-yuint <string>     The text that will be written on the y-axis (by default this comes from the internal SDHDF attributes)\n");

  printf("\nExample:\n\n");
  printf("sdhdf_plotSpectrum -f file.hdf -sb 0\n");

  printf("\n\nThe plot is interactive and the following key/mouse strokes can be used\n\n");
  printf("Left mouse click or A  Display on the screen the current cursor position\n");
  printf(">                   Move to the next sub-band\n");
  printf("<                   Move to the previous sub-band\n");
  printf("+                   Move to the next spectral dump\n");
  printf("-                   Move to the previous spectral dump\n");
  printf("b                   Plot next beam\n");
  printf("B                   Plot previous beam\n");
  printf("f                   Toggle how the flagging is dealt with\n");
  printf("h                   This help\n");
  printf("l                   Toggle displaying the spectrum on a logarithmic scale\n");
  printf("L                   List to the screen the data (frequency, pol1) for the current plot\n");
  printf("m                   Toggle overlaying the rest frequency of the molecular lines\n");
  printf("o                   Write the data corresponding to the current zoom region to a local file\n");
  printf("p                   Toggle plotting 2 polarisations and 1\n");
  printf("q                   Quit. Exit the program\n");
  printf("r                   Toggle overlaying the rest frequency of recombination lines\n");
  printf("s                   Define a region using the mouse and then print the statistical properties of the data in that region\n");
  printf("u                   Unzoom - back to the original plot scales\n");
  printf("w                   Toggle dividing the spectra by the weights in the file\n");
  printf("x                   Toggle the x-axis\n");
  printf("z                   Start a zoom region from the current mouse position, move the mouse and click to set the zoom region\n\n");
}

int main(int argc,char *argv[])
{
  int i,j;
  int nAttr;
  char fname[MAX_STRLEN]="unset";
  sdhdf_fileStruct *inFile;
  double fref=-1;
  char yUnit[MAX_STRLEN] = "not set";
  char freqFrame[MAX_STRLEN] = "[unknown]";
  char freqUnit[MAX_STRLEN] = "unknown";
  char att_yUnit[MAX_STRLEN] = "arbitrary";
  int idump,iband;
  char dataName[MAX_STRLEN];
  int ibeam=0;
  char grDev[128]="/xs";
  int verbose=0;
  
  //  help();
  
  // Defaults
  idump = iband = 0;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

	if (argc==1) {
		help();
		  exit(1);
	  }

  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-f")==0)   	     strcpy(fname,argv[++i]);	
      else if (strcmp(argv[i],"-yunit")==0)  strcpy(yUnit,argv[++i]);	
      else if (strcmp(argv[i],"-fref")==0)   sscanf(argv[++i],"%lf",&fref);
      else if (strcmp(argv[i],"-h")==0)      {help(); exit(1);}
      else if (strcmp(argv[i],"-g")==0)      strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-beam")==0)   sscanf(argv[++i],"%d",&ibeam);
      else if (strcmp(argv[i],"-sd")==0)     sscanf(argv[++i],"%d",&idump);
      else if (strcmp(argv[i],"-v")==0)      verbose=1;
      else if (strcmp(argv[i],"-sb")==0) // FIX ME - ENABLE BAND DESCRIPTORS
	sscanf(argv[++i],"%d",&iband);
    }

  sdhdf_initialiseFile(inFile);
    
  if (strcmp(fname,"unset")==0)
    {
      printf("Please define an input file name using the -f option\n");
      free(inFile);
      exit(1);
    }
  
  if (sdhdf_openFile(fname,inFile,1)==-1)
    {
      printf("Unable to open input file >%s<\n",fname);
      free(inFile);
      exit(1);
    }
  if (verbose==1) printf("Loading meta data\n");
  sdhdf_loadMetaData(inFile);
  if (verbose==1) printf("Completed meta data\n");
  printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname,inFile->primary[0].utc0,
	 inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[ibeam].source,inFile->primary[0].telescope,
	 inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[ibeam].nBand);

  // FIX THIS
  /*
    for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-sb")==0)
	iband = sdhdf_getBandID(inFile,argv[++i]);
    }
  */
  
  // FIX THESE
  //  sdhdf_loadFrequencyAttributes(inFile,iband);
  //  sdhdf_loadDataAttributes(inFile,iband);

  //  strcpy(freqFrame,inFile->frequency_attr.frame);
  //  strcpy(freqUnit,inFile->frequency_attr.unit);
  //  strcpy(att_yUnit,inFile->data_attr.unit);
  
  
	//    strcpy(yUnit,att_yUnit);
  
  plotSpectrum(inFile,ibeam,iband,idump,fref,yUnit,freqFrame,freqUnit,grDev);

  
  sdhdf_closeFile(inFile);
  free(inFile);
}

void plotSpectrum(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,double fref,char *yUnit,char *freqFrame,char *freqUnit,
		  char *grDev)
{
  char key;
  float mx,my;
  int i,j;
  double wts;
  int flagVal;
  float *aa,*bb,*ab,*abs;
  float minx,maxx,miny,maxy,minz,maxz;
  float ominx,omaxx,ominy,omaxy;
  int   divWeights=-1;
  int t=0;
  int setLog=-1;
  char fnameFix[1024];
  char blabelFix[1024];
  char sourceFix[1024];
  char title[1024];
  char ylabel[1024];
  char xlabel[1024];
  char legend[1024];
  int underscoreFix=1;
  int plot=1;
  int molecularLines=0;
  int recombLines=-1;
  int nchan,maxNchan;
  float *pol1,*pol2,*pol3,*pol4;
  float *freq;
  int xplot=1;
  int npol=4;
  int reload=1;
  int flagIt=2;
  int nFreq;
  int allDataFlagged=1;
  
  sdhdf_restfrequency_struct *restFrequencies;
  restFrequencies = (sdhdf_restfrequency_struct *)malloc(sizeof(sdhdf_restfrequency_struct)*MAX_REST_FREQUENCIES);
  sdhdf_loadRestFrequencies(restFrequencies,&nFreq);
  


  cpgbeg(0,grDev,1,1);
  cpgask(0);
  cpgscf(2);
  cpgsch(1.4);
  cpgslw(2);

  maxNchan=0;
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      if (inFile->beam[ibeam].bandHeader[i].nchan > maxNchan)
	maxNchan = inFile->beam[ibeam].bandHeader[i].nchan;
    }

  if (iband > inFile->beam[ibeam].nBand)
    {
      printf("WARNING: Requested band (%d) > number of bands (%d)\n",iband,inFile->beam[ibeam].nBand);
      iband = inFile->beam[ibeam].nBand-1;
      printf("Setting band number to %d\n",iband);
    }
  
  pol1 = (float *)malloc(sizeof(float)*maxNchan);
  pol2 = (float *)malloc(sizeof(float)*maxNchan);
  pol3 = (float *)malloc(sizeof(float)*maxNchan);
  pol4 = (float *)malloc(sizeof(float)*maxNchan);
  freq = (float *)malloc(sizeof(float)*maxNchan);

  npol  = inFile->beam[ibeam].bandHeader[iband].npol;
  
  do
    {
      if (plot==1)
	{
	  nchan = inFile->beam[ibeam].bandHeader[iband].nchan;

	  if (reload==1) // Should reload if band or beam changes
	    {
	      // If already loaded then should release data **
	      // SHOULD ONLY LOAD IF NOT LOADED YET
	      if (inFile->beam[ibeam].bandData[iband].astro_data.pol1AllocatedMemory == 0)
		{
		  sdhdf_loadBandData(inFile,ibeam,iband,1);
		  if (strcmp(yUnit,"not set")==0)
		    {
		      int kk;
		      for (kk=0;kk<inFile->beam[ibeam].bandData[iband].nAstro_obsHeaderAttributes;kk++)
			{
			  // Change .key to .name
			  if (strcmp(inFile->beam[ibeam].bandData[iband].astro_obsHeaderAttr[kk].name,"UNIT")==0)
			    {strcpy(yUnit,inFile->beam[ibeam].bandData[iband].astro_obsHeaderAttr[kk].value); break;}
			}
		    }
		  //  char freqFrame[MAX_STRLEN] = "[unknown]";
		  //  char freqUnit[MAX_STRLEN] = "unknown";
		  //		  printf("freqUnit = %s <<, freqFrame = %s\n",freqUnit,freqFrame);

		  if (strcmp(freqFrame,"[unknown]")==0)
		    {
		      int kk;
		      for (kk=0;kk<inFile->beam[ibeam].bandData[iband].nAstro_obsHeaderAttributes_freq;kk++)
			{
			  //			  printf("Checking attribute %s\n",inFile->beam[ibeam].bandData[iband].astro_obsHeaderAttr_freq[kk].key);

			  // GEORGE: THIS .name was .key **
			  if (strcmp(inFile->beam[ibeam].bandData[iband].astro_obsHeaderAttr_freq[kk].name,"FRAME")==0)
			    {strcpy(freqFrame,inFile->beam[ibeam].bandData[iband].astro_obsHeaderAttr_freq[kk].value);}
			  else if (strcmp(inFile->beam[ibeam].bandData[iband].astro_obsHeaderAttr_freq[kk].name,"UNIT")==0)		    
			    {
			      //			      printf("Value = %s\n",inFile->beam[ibeam].bandData[iband].astro_obsHeaderAttr_freq[kk].value);
			      strcpy(freqUnit,inFile->beam[ibeam].bandData[iband].astro_obsHeaderAttr_freq[kk].value);
			    }
			}
		    }
		  
		}
	      reload=0;
	    }
	  allDataFlagged=1;
	  for (i=0;i<nchan;i++)
	    {
	      wts     = inFile->beam[ibeam].bandData[iband].astro_data.dataWeights[idump*nchan+i]; 
	      flagVal = inFile->beam[ibeam].bandData[iband].astro_data.flag[idump*nchan+i]; 

	      if (xplot==1)
		{
		  if (fref < 0)
		    freq[i] = inFile->beam[ibeam].bandData[iband].astro_data.freq[idump*nchan+i];
		  else
		    freq[i] = (1.0-inFile->beam[ibeam].bandData[iband].astro_data.freq[idump*nchan+i]/(fref))*SPEED_LIGHT/1000.; // km/s
		}
	      else
		freq[i] = i;
	      if (wts > 0 && flagVal == 0)
		{
		  //		  printf("Not flagged in channel %d\n",i);
		  allDataFlagged=0;
		  if (divWeights==1)
		    pol1[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol1[i+idump*nchan]/wts;		 
		  else
		    pol1[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol1[i+idump*nchan];
		  //printf("Loaded %f %f\n",freq[i],pol1[i]);
		  if (divWeights==1)
		    {
		      if (npol > 1) pol2[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol2[i+idump*nchan]/wts;
		    }
		  else
		    {
		      if (npol > 1) pol2[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol2[i+idump*nchan];
		    }
		  
		  if (setLog == 1)
		    {
		      pol1[i]=log10(pol1[i]);
		      if (npol > 1) pol2[i]=log10(pol2[i]);
		    }
		}
	      else
		{
		  pol1[i] = 0;
		  if (npol > 1) pol2[i] = 0;
		}
	    }

	
	  if (t==-1 || t==0 || t==1)
	    {
	      if (t==0 || t==1)
		{
		  int setMiny=0;
		  
		  for (i=0;i<nchan;i++)
		    {
		      //		      printf("IN HERE %g\n",pol1[i]);

		      wts     = inFile->beam[ibeam].bandData[iband].astro_data.dataWeights[idump*nchan+i];
		      flagVal = inFile->beam[ibeam].bandData[iband].astro_data.flag[idump*nchan+i];
		      if ((wts != 0 && flagVal == 0) || flagIt==0)
			{
			  if (setMiny==0)
			    {
			      miny = maxy = pol1[i];
			      setMiny=1;
			    }
			  else
			    {
			      if (pol1[i] > maxy) maxy = pol1[i];
			      if (pol1[i] < miny) miny = pol1[i];
			      
			      if (npol > 1)
				{
				  if (pol2[i] > maxy) maxy = pol2[i];
				  if (pol2[i] < miny) miny = pol2[i];
				}
			    }
			}
		    }
		  ominy = miny; omaxy = maxy;
		}
	   
	      if (t==0)
		{
		  if (xplot==1)
		    {
		      minx = freq[0]; maxx = freq[nchan-1];
		    }
		  else
		    {
		      minx = 0;
		      maxx = nchan;
		    }
		  ominx = minx; omaxx = maxx;
		}
	      else if (t==-1)
		{
		  float sminx,smaxx;
		  float f0,f1;
		  sminx = minx;
		  smaxx = maxx;
		  
		  if (xplot==1) // Convert channel to frequency
		    {
		      minx = freq[(int)sminx];
		      maxx = freq[(int)smaxx];
		      printf("Changing range from (%g,%g) to (%g,%g)\n",sminx,smaxx,minx,maxx);
		    }
		  else // Convert frequency to channel
		    {
		      // FIX ME: using [0] for frequency dump
		      f0 = inFile->beam[ibeam].bandData[iband].astro_data.freq[0];
		      f1 = inFile->beam[ibeam].bandData[iband].astro_data.freq[1];
		      minx = (int)((sminx-f0)/(f1-f0)+0.5);
		      maxx = (int)((smaxx-f0)/(f1-f0)+0.5);
		      printf("Changing range from (%g,%g) to (%g,%g)\n",sminx,smaxx,minx,maxx);
		    }
		}

	      t=2;
	      
	    }
	}
      
      if (minx == maxx || miny == maxy)
	{
	  cpgenv(0,1,0,1,0,0);
	  cpgtext(0.3,0.45,"No data to plot");
	}
      else
	{
	  if (setLog==-1)
	    cpgenv(minx,maxx,miny,maxy,0,0);
	  else
	    cpgenv(minx,maxx,miny,maxy,0,20);
	  sprintf(ylabel,"Signal strength (%s)",yUnit);
	}
      if (underscoreFix==1) // SHOULD HAVE A FUNCTION TO DO THIS
	{
	  sdhdf_fixUnderscore(inFile->fname,fnameFix);
	  sdhdf_fixUnderscore(inFile->beam[ibeam].bandHeader[iband].label,blabelFix);
	  sdhdf_fixUnderscore(inFile->beamHeader[ibeam].source,sourceFix);
	}
      else
	{
	  strcpy(fnameFix,inFile->fname);
	  strcpy(blabelFix,inFile->beam[ibeam].bandHeader[iband].label);
	  strcpy(sourceFix,inFile->beamHeader[ibeam].source);
	}
      if (inFile->beam[ibeam].bandHeader[iband].ndump==1)       
	sprintf(title,"%s, %s",fnameFix,blabelFix);
      else
	sprintf(title,"%s, %s, spectral dump %d",fnameFix,blabelFix,idump);
      if (ibeam > 0)
	{
	  char tStr[128];
	  sprintf(tStr,", beam %d",ibeam);
	  strcat(title,tStr);
	}
      if (fref < 0)
	{
	  sprintf(xlabel,"%s frequency (%s)",freqFrame,freqUnit);
	}
      else
	{
	  sprintf(xlabel,"Velocity (km s\\u-1\\d) [f\\dref\\u = %.2f MHz]",fref);
	}
      cpglab(xlabel,ylabel,title);
	  
      cpgsch(0.9);
      sprintf(legend,"Source: %s",sourceFix);
      cpgtext(minx+0.05*(maxx-minx),maxy-0.05*(maxy-miny),legend);
      sprintf(legend,"Tobs: %.2f sec",inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].dtime);
      cpgtext(minx+0.05*(maxx-minx),maxy-0.1*(maxy-miny),legend);
      sprintf(legend,"RA:  %s",inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].raStr);
      cpgtext(minx+0.05*(maxx-minx),maxy-0.15*(maxy-miny),legend);
      sprintf(legend,"Dec: %s",inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].decStr);
      cpgtext(minx+0.05*(maxx-minx),maxy-0.2*(maxy-miny),legend);
      cpgsch(1.4);
      
      if (molecularLines==1 || molecularLines==2)
	{
	  for (i=0;i<nFreq;i++)	    
	    drawMolecularLine(restFrequencies[i].f0,restFrequencies[i].label,minx,maxx,miny,maxy,molecularLines);
	}
      if (recombLines==1)
	{
	  drawRecombLine(minx,maxx,miny,maxy);
	}
      if (flagIt==2)
	{
	  int i0=0;
	  int drawIt=-1;
	  for (i=0;i<nchan;i++)
	    {
	      if ((inFile->beam[ibeam].bandData[iband].astro_data.dataWeights[i+idump*nchan] != 0 &&
		   inFile->beam[ibeam].bandData[iband].astro_data.flag[i+idump*nchan] == 0)
		   && drawIt==-1)
		{
		  drawIt=1;
		  i0=i;
		}
	      if ((inFile->beam[ibeam].bandData[iband].astro_data.dataWeights[i+idump*nchan] == 0  ||
		   inFile->beam[ibeam].bandData[iband].astro_data.flag[i+idump*nchan] == 1)
		  && drawIt==1)
		{
		  cpgsci(1); cpgline(i-1-i0,freq+i0,pol1+i0); cpgsci(1);
		  if (npol > 1)
		    {
		      cpgsci(2); cpgline(i-1-i0,freq+i0,pol2+i0); cpgsci(1);
		    }
		  drawIt=-1;

		}	      
	    }
	  if (drawIt==1)
	    {
	      cpgsci(1); cpgline(i-1-i0,freq+i0,pol1+i0); cpgsci(1);
	      if (npol > 1)
		{
		  cpgsci(2); cpgline(i-1-i0,freq+i0,pol2+i0); cpgsci(1);
		}
	    }
	}
      else
	{
	  cpgsci(1); cpgline(nchan,freq,pol1); cpgsci(1);
	  if (npol > 1)
	    {
	      cpgsci(2); cpgline(nchan,freq,pol2); cpgsci(1);
	    }
	}

      if (strstr(grDev,"/ps")!=NULL || strstr(grDev,"/cps")!=NULL || strstr(grDev,"/png")!=NULL)
	key='q';
      else
	{
	  cpgcurs(&mx,&my,&key);
	  if (key=='z')
	    {
	      float mx2,my2;
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
	  else if (key=='b')
	    {
	      if (ibeam < inFile->nBeam-1)
		{
		  ibeam++; reload=1;
		}
	    }
	  else if (key=='B')
	    {
	      if (ibeam > 0)
		{
		  ibeam--; reload=1;
		}
	    }
	  else if (key=='w')
	    {
	      divWeights*=-1;
	      if (divWeights==1)
		printf("Dividing spectra by weights\n");
	      else
		printf("Not dividing spectra by weights\n");
	      t=0;
	    }
	  else if (key=='h')
	    help();
	  else if (key=='L')
	    {
	      printf("Number of channels = %d\n",nchan);
	      for (i=0;i<nchan;i++)
		printf("%.6f %g\n",freq[i],pol1[i]);
	    }
	  else if (key=='s')
	    {
	      float mx2,my2;
	      float t;
	      int kk;
	      float sx[nchan],sy1[nchan],sy2[nchan];
	      double ssy,ssy_2,ssy2,ssy2_2,sdev1,sdev2;
	      double min1,max1;
	      double min2,max2;
	      double freqPeak1,freqPeak2;
	      int setMinMax=1;
	      int ns=0;
	      cpgband(4,0,mx,my,&mx2,&my2,&key);
	      if (mx != mx2 && my != my2)
		{
		  if (mx > mx2)
		    {t = mx; mx = mx2; mx2 = t;}

		  printf("Determining statistics between %.6f and %.6f\n",mx,mx2);
		  ssy = ssy_2 = ssy2 = ssy2_2 = 0;
		  for (kk=0;kk<nchan;kk++)
		    {
		      if (freq[kk] > mx && freq[kk] <= mx2)
			{
			  if (setMinMax==1)
			    {
			      min1 = max1 = pol1[kk];
			      min2 = max2 = pol2[kk];
			      freqPeak1 = freqPeak2 = freq[kk];
			      setMinMax=0;
			    }
			  if (min1 > pol1[kk]) min1 = pol1[kk];
			  if (min2 > pol2[kk]) min2 = pol2[kk];

			  if (max1 < pol1[kk]) {max1 = pol1[kk]; freqPeak1 = freq[kk];}
			  if (max2 < pol2[kk]) {max2 = pol2[kk]; freqPeak2 = freq[kk];}
			  sy1[ns]=pol1[kk];
			  sy2[ns]=pol2[kk];
			  ssy += pol1[kk];
			  ssy2 += pol2[kk];
			  ssy_2 += pow(pol1[kk],2);
			  ssy2_2 += pow(pol2[kk],2);
			  ns++;
			}
		    }
		  sdev1 = sqrt(1./(double)ns*ssy_2   - pow(1.0/(double)ns * ssy,2));
		  sdev2 = sqrt(1./(double)ns*ssy2_2 - pow(1.0/(double)ns * ssy2,2));
		  printf("# Pol   Npts   Mean    RMS   Min          Max   FreqMax\n");
		  printf("Pol 1: %d %g %g %g %g %.6f\n",ns,ssy/(double)ns,sdev1,min1,max1,freqPeak1);
		  printf("Pol 2: %d %g %g %g %g %.6f\n",ns,ssy2/(double)ns,sdev2,min2,max2,freqPeak2);

		}
	    }
	  else if (key=='f')
	    {
	      flagIt++;
	      if (flagIt==3) flagIt=0;
	      t=0;
	    }
	  else if (key=='A')
	    printf("Mouse cursor = (%g,%g)\n",mx,my);
	  else if (key=='1') {plot=1; t=0;} // DOES THIS DO ANYTHING?
	  else if (key=='2') {plot=2; t=1;} // DOES THIS DO ANYTHING?
	  else if (key=='l')
	    {setLog*=-1; t=0;}
	  else if (key=='p')
	    {
	      if (npol==1) npol=4;
	      else npol=1;
	    }
	  else if (key=='>')
	    {	  
	      iband++;
	      if (iband >= inFile->beam[ibeam].nBand)
		iband = inFile->beam[ibeam].nBand-1;
	      reload=1;
	      t=0;
	      
	    }
	  else if (key=='<')
	    {
	      iband--;
	      if (iband < 0) iband = 0;
	      reload=1;
	      t=0;
	    }
	  else if (key=='+')
	    {
	      if (idump < inFile->beam[ibeam].bandHeader[iband].ndump-1)
		idump++;
	    }
	  else if (key=='x')
	    {
	      xplot*=-1;
	      t=-1;
	    }
	  else if (key=='m')
	    {
	      molecularLines++;
	      if (molecularLines==3) molecularLines=0;
	    }
	    else if (key=='r')
	    recombLines*=-1;
	  else if (key=='-')
	    {
	      if (idump > 0)
		idump--;
	    }
	  else if (key=='o')
	    {
	      char fname[1024];
	      FILE *fout;
	      printf("Enter output filename ");
	      scanf("%s",fname);
	      fout = fopen(fname,"w");
	      for (i=0;i<nchan;i++)
		{
		  if (freq[i] > minx && freq[i] <= maxx)
		    {
		      if (npol > 1)
			fprintf(fout,"%.6f %g %g\n",freq[i],pol1[i],pol2[i]);
		      else
			fprintf(fout,"%.6f %g\n",freq[i],pol1[i]);
		    }
		}
	  fclose(fout);
	    }
	  else if (key=='u')
	    {
	      minx = ominx;
	      maxx = omaxx;
	      miny = ominy;
	      maxy = omaxy;
	    }
	}
  } while (key!='q');
  
  cpgend();

  free(pol1);
  free(pol2);
  free(pol3);
  free(pol4);
  free(restFrequencies);
  
  /*
  for (i=0;i<inFile->bandHeader[iband].nchan;i++)
    {
      printf("%d %d %d %.5f %g %g %g %g %g %g %g %g\n",iband,idump,i,inFile->bandHeader[iband].topoFreq[i],
	     spectrum.pol1[i].val,spectrum.pol2[i].val,spectrum.pol3[i].val,spectrum.pol4[i].val,
	     spectrum.pol1[i].weight,spectrum.pol2[i].weight,spectrum.pol3[i].weight,spectrum.pol4[i].weight);
    }
  */  

}

void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY,int drawLabel)
{
  float fx[2],fy[2];

  fx[0] = fx[1] = freq;
  fy[0] = minY;
  fy[1] = maxY;
  cpgsci(3); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);
  if (drawLabel==2)
    {
      cpgsch(0.8);
      cpgtext(fx[0],fy[1]-0.1*(fy[1]-fy[0]),label);
      cpgsch(1.4);
    }  
}

//
// See https://www.cv.nrao.edu/course/astr534/Recombination.html
//
void drawRecombLine(float minX,float maxX,float minY,float maxY)
{
  float fx[2],fy[2];
  double freq;
  double RM_c;
  double c = 299792458.0;
  int n;
  int dn;
  double Rinf_c = 3.28984e15; // Rydberg frequency (Hz)
  double me = 9.10938356e-31; // Kg
  double M = 1836.1 * me; // Hydrogen nucleus is a single proton of mass mp ~ 1836.1m_e
  char text[1024];
  char greek[6];
  
  //  n = 187;
  //  dn = 1;

  for (dn=1;dn<5;dn++)
    {
      if (dn==1)
	sprintf(greek,"\\ga");
      else if (dn==2)
	sprintf(greek,"\\gb");
      else if (dn==3)
	sprintf(greek,"\\gg");
      else if (dn==4)
	sprintf(greek,"\\gd");
      for (n=1;n<400;n++)
	{
	  RM_c = Rinf_c *pow(1+me/M,-1);  
	  freq = RM_c*(1./pow(n,2) - 1./pow(n+dn,2))/1e6; // MHz
	  if (freq > minX && freq < maxX)
	    {
	      fx[0] = fx[1] = freq;
	      fy[0] = minY;
	      fy[1] = maxY;
	      cpgsci(5);
	      cpgsls(2);
	      cpgline(2,fx,fy);
	      sprintf(text,"H %d%s",n,greek);
	      cpgsch(0.8);
	      cpgtext(fx[0],minY + (maxY-minY)*0.9,text);
	      cpgsch(1.4);
	      cpgsci(1);
	      cpgsls(1);
	      printf("%d %d Freq = %g\n",dn,n,freq);
	    }
	}
    }
  
}
