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
// Software to process and view the calibrator
//
// Usage:
// sdhdf_cal
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>


void doPlot(sdhdf_fileStruct *inFile,int beam,int totChan,int nScal,float *scalFreq,float *scalAA,float *scalBB,int av,int SorT);
void obtainScal(float freq,float *scalFreq,float *scalAA,float *scalBB,int nScal,float *retAA,float *retBB,int interp);

void help()
{
  printf("sdhdf_plotNoiseCal\n");
  printf("Provides ways to visualise the calibration signal/noise source\n");
  printf("\n\n");
  printf("-autoSsys           Automatically determine Ssys and print out to screen\n");
  printf("-autoSysGain        Automatically determine system gain and print out to screen\n");
  printf("-autoDiffPhase      Automatically determine differential phase and print out to screen\n");
  printf("-av                 Average the noise source information in time\n");
  printf("-f <filename>       Input filename\n");
  printf("-h                  This help\n");
  printf("-scal <filename>    Filename of Scal (Jy) information\n");
  printf("-scale <float>      Value to scale the calibrator values by\n");
  printf("\n\n");
  printf("The output is an interactive plot\n\n");
  printf("1                   Plot cal ON for polarisations 1 and 2\n");
  printf("2                   Plot cal ON and cal OFF for polarisation 1\n");
  printf("3                   Plot OFF/(ON-OFF) for both polarisations\n");
  printf("4                   Plot Ssys for both polaristaions\n");
  printf("5                   Plot differential gain");
  printf("6                   Plot differential phase");
  printf("7                   Plot ON,OFF for polarisations 1 and 2");
  printf("+                   Move to next time dump\n");
  printf("-                   Move to previous time dump\n");
  printf("a                   Toggle time averaging");
  printf("h                   This help\n");
  printf("q                   Quit\n");
  printf("u                   Unzoom\n");
  printf("w                   Write spectra to disk as a text file\n");
  printf("z                   Use mouse to zoom into a specific region\n");
}


int main(int argc,char *argv[])
{
  int i,j,k;
  char fname[MAX_STRLEN];
  sdhdf_fileStruct *inFile;
  int iband=0,nchan,idump;
  int ndump=0;
  int beam=0;
  int totChan=0;
  float *scalFreq,*scalAA,*scalBB;
  float retAA,retBB;
  char scalFname[1024]="NULL";
  char tcalFname[1024]="NULL";
  int SorT=0;
  int  nScal=-1;
  float scaleCal=1;
  int autoSsys=0;
  int autoSysGain=0;
  int autoDiffPhase=0;
  int av=0;
  int verbose=0;
  sdhdf_fluxCalibration *fluxCal;

  printf("Allocating memory\n");
  fluxCal = (sdhdf_fluxCalibration *)malloc(sizeof(sdhdf_fluxCalibration)*3328); //  Should change to MAX_FLUXCAL OR SIMILAR ** FIX ME
  
  scalFreq = (float *)malloc(sizeof(float)*3328);
  scalAA   = (float *)malloc(sizeof(float)*3328);
  scalBB   = (float *)malloc(sizeof(float)*3328);

  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
  printf("Initialising the data file\n");
  sdhdf_initialiseFile(inFile);
  printf("Reading the command line argument\n");
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);	
      else if (strcmp(argv[i],"-v")==0)
	verbose=1;
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-scal")==0)
	strcpy(scalFname,argv[++i]);
      else if (strcmp(argv[i],"-tcal")==0)
	strcpy(tcalFname,argv[++i]);
      else if (strcmp(argv[i],"-scale")==0)
	sscanf(argv[++i],"%f",&scaleCal); 
      else if (strcmp(argv[i],"-autoSsys")==0)
	autoSsys=1;
      else if (strcasecmp(argv[i],"-autoDiffPhase")==0)
	autoDiffPhase=1;
     else if (strcmp(argv[i],"-autoSysGain")==0)
	autoSysGain=1;
      else if (strcmp(argv[i],"-av")==0)
	av=1;
    }
  printf("Opening the data file\n");
  sdhdf_openFile(fname,inFile,1);
  printf("File opened\n");
  sdhdf_loadMetaData(inFile);
  printf("Loaded metadata\n");

  if (strcmp(scalFname,"NULL")!=0)
    {
      FILE *fin;
      char obsDir[1024];
      int nFluxCalChan,kk;
      sdhdf_getTelescopeDirName(inFile->primary[0].telescope,obsDir);
      sdhdf_loadFluxCal(fluxCal,&nFluxCalChan,obsDir,inFile->primary[0].rcvr,scalFname,0.0); 
      nScal=nFluxCalChan;
      for (kk=0;kk<nFluxCalChan;kk++)
	{
	  scalFreq[kk] = fluxCal[kk].freq;
	  scalAA[kk] = fluxCal[kk].scalAA * scaleCal;
	  scalBB[kk] = fluxCal[kk].scalBB * scaleCal;
	}
      SorT = 1;
      printf("Loaded %s, Number of SCAL values is %d\n",scalFname,nScal);
    }
  else if (strcmp(tcalFname,"NULL")!=0)
    {
      FILE *fin;
      char obsDir[1024];
      char fname[1024];
      int kk;
      
      sdhdf_getTelescopeDirName(inFile->primary[0].telescope,obsDir);
      //      sprintf(fname,"%s/observatory/%s/rfi/persistentRFI.dat",runtimeDir,dirName);
      
      sprintf(fname,"%s/observatory/%s/calibration/UWL/%s",getenv("SDHDF_RUNTIME"),obsDir,tcalFname);
      printf("Opening %s\n",fname);
      kk=0;
      fin = fopen(fname,"r");
      while (!feof(fin))
	{
	  if (fscanf(fin,"%f %f %f",&scalFreq[kk],&scalAA[kk],&scalBB[kk])==3)
	    kk++;
	}
      fclose(fin);
      nScal = kk;
      SorT=2;
      printf("Loaded %s, Number of TCAL values is %d\n",tcalFname,nScal);
    }
    

  ndump = inFile->beam[beam].calBandHeader[0].ndump;
  printf("ndumps = %d\n",ndump);
  for (i=0;i<inFile->beam[beam].nBand;i++)
    {
      printf("nchan = %d\n",inFile->beam[beam].calBandHeader[i].nchan);

      nchan = inFile->beam[beam].calBandHeader[i].nchan;
      totChan+=nchan;
      sdhdf_loadBandData(inFile,beam,i,2);
      sdhdf_loadBandData(inFile,beam,i,3);
    }
  if (autoSysGain==1 && av==0)
    {
      double freq,s1,s2;
      int np=0;
      float onP1,offP1,onP2,offP2;
      float onP3,offP3,onP4,offP4;
      
      for (k=0;k<inFile->beam[beam].calBandHeader[0].ndump;k++)
	{	  
	  np=0;
	  for (i=0;i<inFile->beam[beam].nBand;i++)
	    {
	      nchan = inFile->beam[beam].calBandHeader[i].nchan;
	      for (j=0;j<nchan;j++)
		{
		  freq = inFile->beam[beam].bandData[i].cal_on_data.freq[j];  // FIX ME -- FREQ AXIS FOR DUMP
		  onP1  = inFile->beam[beam].bandData[i].cal_on_data.pol1[j+k*nchan];
		  offP1 = inFile->beam[beam].bandData[i].cal_off_data.pol1[j+k*nchan];
		  onP2  = inFile->beam[beam].bandData[i].cal_on_data.pol2[j+k*nchan];
		  offP2 = inFile->beam[beam].bandData[i].cal_off_data.pol2[j+k*nchan];
		  onP3  = inFile->beam[beam].bandData[i].cal_on_data.pol3[j+k*nchan];
		  offP3 = inFile->beam[beam].bandData[i].cal_off_data.pol3[j+k*nchan];
		  onP4  = inFile->beam[beam].bandData[i].cal_on_data.pol4[j+k*nchan];
		  offP4 = inFile->beam[beam].bandData[i].cal_off_data.pol4[j+k*nchan];

		  //
		  // Should find the closest SCAL measurement to the frequency
		  //
		  obtainScal(freq,scalFreq,scalAA,scalBB,nScal,&retAA,&retBB,1);
		  s1    = (onP1-offP1)/retAA;
		  s2    = (onP2-offP2)/retBB;
		  printf("SysGain %.6f %g %g\n",freq,s1,s2);
		  np++;
		}
	    }
	  printf("SysGain\n");
	}
    }
  else if (autoSysGain==1 && av==1)
    {
      double freq,s1,s2;
      int np=0;
      float onP1,offP1,onP2,offP2;
      float onP3,offP3,onP4,offP4;
      
      np=0;
      for (i=0;i<inFile->beam[beam].nBand;i++)
	{
	  nchan = inFile->beam[beam].calBandHeader[i].nchan;
	  for (j=0;j<nchan;j++)
	    {
	      freq = inFile->beam[beam].bandData[i].cal_on_data.freq[j];  // FIX ME -- FREQ AXIS FOR DUMP
	      onP1=offP1=onP2=offP2=onP3=offP3=onP4=offP4=0;
	      for (k=0;k<inFile->beam[beam].calBandHeader[0].ndump;k++)
		{	  		  
		  onP1  += inFile->beam[beam].bandData[i].cal_on_data.pol1[j+k*nchan];
		  offP1 += inFile->beam[beam].bandData[i].cal_off_data.pol1[j+k*nchan];
		  onP2  += inFile->beam[beam].bandData[i].cal_on_data.pol2[j+k*nchan];
		  offP2 += inFile->beam[beam].bandData[i].cal_off_data.pol2[j+k*nchan];
		  onP3  += inFile->beam[beam].bandData[i].cal_on_data.pol3[j+k*nchan];
		  offP3 += inFile->beam[beam].bandData[i].cal_off_data.pol3[j+k*nchan];
		  onP4  += inFile->beam[beam].bandData[i].cal_on_data.pol4[j+k*nchan];
		  offP4 += inFile->beam[beam].bandData[i].cal_off_data.pol4[j+k*nchan];
		}
	      obtainScal(freq,scalFreq,scalAA,scalBB,nScal,&retAA,&retBB,1);
	      s1 = (onP1-offP1)/retAA;
	      s2 = (onP2-offP2)/retBB;

	      printf("AvSysGain %.6f %g %g\n",freq,s1,s2);
	      np++;
	    }
	}
    }
 else  if (autoSsys==1 && av==0)
    {
      double freq,s1,s2;
      int np=0;
      float onP1,offP1,onP2,offP2;
      float onP3,offP3,onP4,offP4;
      
      for (k=0;k<inFile->beam[beam].calBandHeader[0].ndump;k++)
	{	  
	  np=0;
	  for (i=0;i<inFile->beam[beam].nBand;i++)
	    {
	      nchan = inFile->beam[beam].calBandHeader[i].nchan;
	      for (j=0;j<nchan;j++)
		{
		  freq = inFile->beam[beam].bandData[i].cal_on_data.freq[j];  // FIX ME -- FREQ AXIS FOR DUMP
		  onP1  = inFile->beam[beam].bandData[i].cal_on_data.pol1[j+k*nchan];
		  offP1 = inFile->beam[beam].bandData[i].cal_off_data.pol1[j+k*nchan];
		  onP2  = inFile->beam[beam].bandData[i].cal_on_data.pol2[j+k*nchan];
		  offP2 = inFile->beam[beam].bandData[i].cal_off_data.pol2[j+k*nchan];
		  onP3  = inFile->beam[beam].bandData[i].cal_on_data.pol3[j+k*nchan];
		  offP3 = inFile->beam[beam].bandData[i].cal_off_data.pol3[j+k*nchan];
		  onP4  = inFile->beam[beam].bandData[i].cal_on_data.pol4[j+k*nchan];
		  offP4 = inFile->beam[beam].bandData[i].cal_off_data.pol4[j+k*nchan];
		  obtainScal(freq,scalFreq,scalAA,scalBB,nScal,&retAA,&retBB,1);
		  s1 = retAA*offP1/(onP1-offP1);
		  s2 = retBB*offP2/(onP2-offP2);
		  printf("Ssys %s %d %d %d %.6f %g %g %g %g\n",inFile->fname,k,i,j,freq,s1,s2,inFile->beam[beam].bandData[i].cal_obsHeader[k].el,inFile->beam[beam].bandData[i].cal_obsHeader[k].paraAngle);
		  np++;
		}
	    }
	  printf("Ssys\n");
	}
    }
  else if (autoSsys==1 && av==1)
    {
      double freq,s1,s2;
      int np=0;
      float onP1,offP1,onP2,offP2;
      float onP3,offP3,onP4,offP4;
      
      np=0;
      for (i=0;i<inFile->beam[beam].nBand;i++)
	{
	  nchan = inFile->beam[beam].calBandHeader[i].nchan;
	  for (j=0;j<nchan;j++)
	    {
	      freq = inFile->beam[beam].bandData[i].cal_on_data.freq[j];  // FIX ME -- FREQ AXIS FOR DUMP
	      onP1=offP1=onP2=offP2=onP3=offP3=onP4=offP4=0;
	      for (k=0;k<inFile->beam[beam].calBandHeader[0].ndump;k++)
		{	  		  
		  onP1  += inFile->beam[beam].bandData[i].cal_on_data.pol1[j+k*nchan];
		  offP1 += inFile->beam[beam].bandData[i].cal_off_data.pol1[j+k*nchan];
		  onP2  += inFile->beam[beam].bandData[i].cal_on_data.pol2[j+k*nchan];
		  offP2 += inFile->beam[beam].bandData[i].cal_off_data.pol2[j+k*nchan];
		  onP3  += inFile->beam[beam].bandData[i].cal_on_data.pol3[j+k*nchan];
		  offP3 += inFile->beam[beam].bandData[i].cal_off_data.pol3[j+k*nchan];
		  onP4  += inFile->beam[beam].bandData[i].cal_on_data.pol4[j+k*nchan];
		  offP4 += inFile->beam[beam].bandData[i].cal_off_data.pol4[j+k*nchan];
		}
	      obtainScal(freq,scalFreq,scalAA,scalBB,nScal,&retAA,&retBB,1);
	      s1 = retAA*offP1/(onP1-offP1);
	      s2 = retBB*offP2/(onP2-offP2);
	      printf("AvSsys %.6f %g %g\n",freq,s1,s2);
	      np++;
	    }
	}
    }
  else if (autoDiffPhase==1)
    {
      double freq,s1,s2;
      int np=0;
      float onP1,offP1,onP2,offP2;
      float onP3,offP3,onP4,offP4;
      
      np=0;
      for (i=0;i<inFile->beam[beam].nBand;i++)
	{
	  nchan = inFile->beam[beam].calBandHeader[i].nchan;
	  for (j=0;j<nchan;j++)
	    {
	      freq = inFile->beam[beam].bandData[i].cal_on_data.freq[j];  // FIX ME -- FREQ AXIS FOR DUMP
	      onP1=offP1=onP2=offP2=onP3=offP3=onP4=offP4=0;
	      for (k=0;k<inFile->beam[beam].calBandHeader[0].ndump;k++)
		{	  		  
		  onP1  += inFile->beam[beam].bandData[i].cal_on_data.pol1[j+k*nchan];
		  offP1 += inFile->beam[beam].bandData[i].cal_off_data.pol1[j+k*nchan];
		  onP2  += inFile->beam[beam].bandData[i].cal_on_data.pol2[j+k*nchan];
		  offP2 += inFile->beam[beam].bandData[i].cal_off_data.pol2[j+k*nchan];
		  onP3  += inFile->beam[beam].bandData[i].cal_on_data.pol3[j+k*nchan];
		  offP3 += inFile->beam[beam].bandData[i].cal_off_data.pol3[j+k*nchan];
		  onP4  += inFile->beam[beam].bandData[i].cal_on_data.pol4[j+k*nchan];
		  offP4 += inFile->beam[beam].bandData[i].cal_off_data.pol4[j+k*nchan];
		}
	      float ucal = 2*(onP3-offP3);
	      float vcal = 2*(onP4-offP4);
	      
	      printf("AvDiffPhase %.6f %g\n",freq,atan2(vcal,ucal)*180/M_PI);
	      np++;
	    }
	}
    }

  else
    doPlot(inFile,beam,totChan,nScal,scalFreq,scalAA,scalBB,av,SorT);

    
  sdhdf_closeFile(inFile);
  free(inFile);
  free(fluxCal);
  free(scalFreq); free(scalAA); free(scalBB);
}

void doPlot(sdhdf_fileStruct *inFile,int beam,int totChan,int nScal,float *scalFreq,float *scalAA,float *scalBB,int av,int SorT)
{
  float fx[totChan];
  float fy1[totChan];
  float fy2[totChan];
  float fy3[totChan];
  float fy4[totChan];
  int np=0;
  int i,j,k;
  float mx,my;
  char key;
  int idump=0;
  int nchan;
  int plot=1;
  int nLine=2;
  float onP1,offP1,onP2,offP2;
  float onP3,offP3,onP4,offP4;
  float miny,maxy;
  float uminx,umaxx;
  float uminy,umaxy;
  float set_miny=-1,set_maxy=-1;
  float set_minx=-1,set_maxx=-1;
  int log=0;
  char title[1024];
  //  int av=-1;
  float retAA,retBB;

  printf("nScal here is =%d\n",nScal);
  
  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  cpgslw(2);
  cpgscf(2);
  cpgsch(1.4);
  do {    
    np=0;
    for (i=0;i<inFile->beam[beam].nBand;i++)
      {
	nchan = inFile->beam[beam].calBandHeader[i].nchan;
	for (j=0;j<nchan;j++)
	  {
	    fx[np] = inFile->beam[beam].bandData[i].cal_on_data.freq[j];  // FIX ME -- FREQ AXIS FOR DUMP
	    if (av==0)
	      {
		onP1  = inFile->beam[beam].bandData[i].cal_on_data.pol1[j+idump*nchan];
		offP1 = inFile->beam[beam].bandData[i].cal_off_data.pol1[j+idump*nchan];
		onP2  = inFile->beam[beam].bandData[i].cal_on_data.pol2[j+idump*nchan];
		offP2 = inFile->beam[beam].bandData[i].cal_off_data.pol2[j+idump*nchan];
		onP3  = inFile->beam[beam].bandData[i].cal_on_data.pol3[j+idump*nchan];
		offP3 = inFile->beam[beam].bandData[i].cal_off_data.pol3[j+idump*nchan];
		onP4  = inFile->beam[beam].bandData[i].cal_on_data.pol4[j+idump*nchan];
		offP4 = inFile->beam[beam].bandData[i].cal_off_data.pol4[j+idump*nchan];
	      }
	    else
	      {
		onP1=offP1=onP2=offP2=onP3=offP3=onP4=offP4 = 0.0;
		for (k=0;k<inFile->beam[beam].calBandHeader[0].ndump;k++)
		  {
		    onP1  += inFile->beam[beam].bandData[i].cal_on_data.pol1[j+k*nchan];
		    offP1 += inFile->beam[beam].bandData[i].cal_off_data.pol1[j+k*nchan];
		    onP2  += inFile->beam[beam].bandData[i].cal_on_data.pol2[j+k*nchan];
		    offP2 += inFile->beam[beam].bandData[i].cal_off_data.pol2[j+k*nchan];
		    onP3  += inFile->beam[beam].bandData[i].cal_on_data.pol3[j+k*nchan];
		    offP3 += inFile->beam[beam].bandData[i].cal_off_data.pol3[j+k*nchan];
		    onP4  += inFile->beam[beam].bandData[i].cal_on_data.pol4[j+k*nchan];
		    offP4 += inFile->beam[beam].bandData[i].cal_off_data.pol4[j+k*nchan];
		  }

		/*
		  onP1  /= (double)inFile->beam[beam].calBandHeader[0].ndump;
		  offP1 /= (double)inFile->beam[beam].calBandHeader[0].ndump;
		  onP2  /= (double)inFile->beam[beam].calBandHeader[0].ndump;
		  offP2 /= (double)inFile->beam[beam].calBandHeader[0].ndump;
		  onP3  /= (double)inFile->beam[beam].calBandHeader[0].ndump;
		  offP3 /= (double)inFile->beam[beam].calBandHeader[0].ndump;
		  onP4  /= (double)inFile->beam[beam].calBandHeader[0].ndump;
		  offP4 /= (double)inFile->beam[beam].calBandHeader[0].ndump;
		*/						
	      }
	    //	    printf("cal: %.6f %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g \n",
	    //		   inFile->beam[beam].bandData[i].cal_on_data.freq[j],onP1,offP1,onP2,offP2,onP3,offP3,onP4,offP4);
	    if (plot==1)
	      {
		fy1[np] = onP1;
		fy2[np] = onP2;
		
		fy1[np] = log10(fy1[np]);
		fy2[np] = log10(fy2[np]);
		nLine=2;
		log=1;
	      }
	    else if (plot==2)
	      {
		fy1[np] = onP1;
		fy2[np] = offP1;
		fy1[np] = log10(fy1[np]);
		fy2[np] = log10(fy2[np]);
		nLine=2;
		log=1;
	      }
	    else if (plot==3)
	      {
		fy1[np] = offP1/(onP1-offP1);
		fy2[np] = offP2/(onP2-offP2);
		printf("Have %g %g\n",fy1[np],fy2[np]);
		nLine=2;
		log=0;
	      }
	    else if (plot==4)
	      {
		obtainScal(fx[np],scalFreq,scalAA,scalBB,nScal,&retAA,&retBB,1);
		fy1[np] = retAA*offP1/(onP1-offP1);
		fy2[np] = retBB*offP2/(onP2-offP2);
		//		printf("Have %.6f %g %g %g %g %g %g %g %g\n",fx[np],fy1[np],fy2[np],retAA,retBB,onP1,offP1,onP2,offP2);
		nLine=2;
		log=0;
	      }
	    else if (plot==5)
	      {
		float ical = (onP1-offP1) + (onP2-offP2);
		float qcal = (onP1-offP1) - (onP2-offP2);

		fy1[np] = 2*qcal/ical;
		printf("Have %g %g\n",fy1[np],fy2[np]);
		nLine=1;
		log=0;
	      }
	    else if (plot==6)
	      {
		float ucal = 2*(onP3-offP3);
		float vcal = 2*(onP4-offP4);
		//		fy1[np] = atan2(vcal,ucal)*180/M_PI;

		// Note PSRCHIVE has a - sign here
		fy1[np] = atan2(vcal,ucal)*180/M_PI;
		printf("Have %g %g\n",fy1[np],fy2[np]);
		nLine=1;
		log=0;
	      }
	    else if (plot==7)
	      {
		fy1[np] = onP1;
		fy2[np] = offP1;
		fy3[np] = onP2;
		fy4[np] = offP2;
		nLine=4;
		log=0;
	      }
	    np++;
	  }
      }
    for (i=0;i<np;i++)
      {
	if (i==0)
	  {
	    miny = maxy = fy1[i];
	  }
	if (miny > fy1[i]) miny = fy1[i];
	if (maxy < fy1[i]) maxy = fy1[i];
	if (nLine == 2)
	  {
	    if (miny > fy2[i]) miny = fy2[i];
	    if (maxy < fy2[i]) maxy = fy2[i];
	  }
	if (nLine == 4)
	  {
	    if (miny > fy2[i]) miny = fy2[i];
	    if (maxy < fy2[i]) maxy = fy2[i];
	    if (miny > fy3[i]) miny = fy3[i];
	    if (maxy < fy3[i]) maxy = fy3[i];
	    if (miny > fy4[i]) miny = fy4[i];
	    if (maxy < fy4[i]) maxy = fy4[i];
	  }
      }
    if (set_miny==-1 && set_maxy==-1)
      {
	uminy = miny;
	umaxy = maxy;
	uminx = 704;
	umaxx = 4032;
      }
    else
      {
	uminy = set_miny;
	umaxy = set_maxy;
	uminx = set_minx;
	umaxx = set_maxx;
      }
    if (log==1)
      cpgenv(uminx,umaxx,uminy,umaxy,0,20);
    else
      cpgenv(uminx,umaxx,uminy,umaxy,0,1);

    if (av==0)       sprintf(title,"spectral dump %d",idump);
    else              sprintf(title,"averaged spectral dump");

    if (plot==1)      cpglab("Frequency (MHz)","Cal On for pol 1 and 2",title);
    else if (plot==2) cpglab("Frequency (MHz)","Cal On and Cal off for pol 1",title);
    else if (plot==3) cpglab("Frequency (MHz)","OFF/(ON-OFF)",title);
    else if (plot==4)
      {
	if (SorT == 1) cpglab("Frequency (MHz)","S_sys (Jy)",title);
	else if (SorT == 2) cpglab("Frequency (MHz)","T_sys (K)",title);
	else cpglab("Frequency (MHz)","Sys (unknown)",title);
      }
    else if (plot==5) cpglab("Frequency (MHz)","Diff. gain",title);
    else if (plot==6) cpglab("Frequency (MHz)","Diff. phase",title);
    else              cpglab("Frequency (MHz)","",title);

    cpgsci(1); cpgline(np,fx,fy1);
    if (nLine == 2)
      {
	cpgsci(2); cpgline(np,fx,fy2);
      }
    else if (nLine==4)
      {
	cpgsci(2); cpgline(np,fx,fy2);
	cpgsci(3); cpgline(np,fx,fy3);
	cpgsci(4); cpgline(np,fx,fy4);
      }
    
    cpgsci(1);


    cpgcurs(&mx,&my,&key);
    if (key=='1') {plot=1; set_miny = set_maxy = -1;}
    else if (key=='2') {plot=2;set_miny = set_maxy = -1;}
    else if (key=='3') {plot=3;set_miny = set_maxy = -1;}
    else if (key=='4') {plot=4;set_miny = set_maxy = -1;}
    else if (key=='5') {plot=5;set_miny = set_maxy = -1;}
    else if (key=='6') {plot=6;set_miny = set_maxy = -1;}
    else if (key=='7') {plot=7;set_miny = set_maxy = -1;}
    else if (key=='a')
      {
	if (av==0) av=1;
	else av=0;
      }
	else if (key=='+')
      {
	if (idump < inFile->beam[beam].calBandHeader[0].ndump-1)
	  idump++;

      }
    else if (key=='-')
      {
	if (idump > 0)
	  idump--;
      }
    else if (key=='w')
      {
	FILE *fout;
	char outFile[128];
	printf("Loaded: %s\n",inFile->fname);
	printf("Enter output filename ");
	scanf("%s",outFile);
	fout = fopen(outFile,"w");
	for (i=0;i<np;i++)
	  {
	    if (nLine == 2)	
	      fprintf(fout,"%.6f %g %g\n",fx[i],fy1[i],fy2[i]);
	    else if (nLine == 4)	
	      fprintf(fout,"%.6f %g %g %g %g\n",fx[i],fy1[i],fy2[i],fy3[i],fy4[i]);
	    else
	      fprintf(fout,"%.6f %g\n",fx[i],fy1[i]);
	  }
	fclose(fout);
      }
    else if (key=='u')
      {
	set_miny = set_maxy = -1;
      }
    else if (key=='z')
      {
	float mx2,my2;
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	if (mx != mx2 && my != my2)
	  {
	    if (mx < mx2)
	      {set_minx = mx; set_maxx = mx2;}
	    else
	      {set_minx = mx2; set_maxx = mx;}
	    
	    if (my < my2)
	      {set_miny = my; set_maxy = my2;}
	    else
	      {set_miny = my2; set_maxy = my;}
	  }
	else
	  printf("Please press 'z' and then move somewhere and click left mouse button\n");
	
      }
  } while (key!='q');

  
}

//
// Obtain an SCAL determination from an array for a given frequency
// We assume that the frequency channels are increasing
//
void obtainScal(float freq,float *scalFreq,float *scalAA,float *scalBB,int nScal,float *retAA,float *retBB,int interp)
{
  float f0,f1;
  double df;
  int i,i0;

  f0 = scalFreq[0];
  f1 = scalFreq[nScal-1];

  i0=0;
  //  printf("nScal i = 0 to %d\n",nScal);
  for (i=0;i<nScal;i++)
    {
      //      printf("Have, nScal i = %d %g and %g\n",i,freq,scalFreq[i]);
      if (freq < scalFreq[i])
	{i0=i; break;}
    }

  if (interp==1)
    {
      *retAA = scalAA[i0];
      *retBB = scalBB[i0];
    }
  else
    {
      printf("ERROR: in obtainScal, interp style not implemented\n");
      *retAA = 0;
      *retBB = 0;
    }
}
  
 
