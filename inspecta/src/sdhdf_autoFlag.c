// From file *** MUST IMPLEMENT
// Shouldn't copy weights/flags if not requested
// Band edge RFI
// Fix issue about loading all band data
//

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
// Software to flag data automatically
//
// Usage:
// sdhdf_autoFlag
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>
#include "TKfit.h"

#define MAX_RFI 512

void flagDataPoint(sdhdf_fileStruct *inFile,int b,int i,int j,int s);
void flagPersistentRFI(sdhdf_fileStruct *inFile,sdhdf_rfi *rfi,int nRFI,int b,int i);
void flagTransientRFI(sdhdf_fileStruct *inFile,sdhdf_transient_rfi *transient_rfi,int nTransientRFI,int b,int i);
void flagAutoDump(sdhdf_fileStruct *inFile,float autoFlagDump_sigma,int b,int i);
void flagNarrowRFI(sdhdf_fileStruct *inFile,int b,int i);

void help()
{
  printf("sdhdf_autoFlag:   code to carry out automatic flagging of a data file\n");	 
  exit(1);
}


int main(int argc,char *argv[])
{
  sdhdf_fileStruct *inFile,*fromFile;
  int i,j,ii,b,k,s;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  char extension[MAX_STRLEN]="autoFlag";

  int nTransientRFI=0;
  int nRFI=0;
  int setWt = 0;
  int nManualRFI = 0;
  
  char oname[MAX_STRLEN];
  sdhdf_fileStruct *outFile;
  int nband;
  
  int fromFlag = 0; // = 1 => copy flag table from another file
  int bandEdgeFlag = 0;
  int persistentFlag = 0;
  int transientFlag = 0;
  int narrowFlag = 0;

  int resetWeights=0;
  
  double freq;
  int flagIt;
  int haveFlag;
  int copyType;
  
  float bandEdge=0;
  sdhdf_rfi rfi[MAX_RFI];
  sdhdf_rfi manual_rfi[MAX_RFI];
  sdhdf_transient_rfi transient_rfi[MAX_RFI];
  char args[MAX_ARGLEN]="";
  int autoFlagDump=0;
  float autoFlagDump_sigma=3;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }
  sdhdf_storeArguments(args,MAX_ARGLEN,argc,argv);
  for (i=1;i<argc;i++)
    {      
      //      if (strcmp(argv[i],"-from")==0)	{strcpy(fromFileName,argv[++i]); flagType=1; fromFlag=1;}
      if (strcmp(argv[i],"-e")==0) strcpy(extension,argv[++i]);
      else if (strcmp(argv[i],"-edge")==0) {sscanf(argv[++i],"%f",&bandEdge); bandEdgeFlag=1;}
      else if (strcmp(argv[i],"-persistent")==0) {persistentFlag=1;}
      else if (strcmp(argv[i],"-transient")==0) {transientFlag=1;}
      else if (strcasecmp(argv[i],"-resetWeights")==0) {resetWeights=1;}
      else if (strcmp(argv[i],"-narrow")==0) {narrowFlag=1;}
      else if (strcasecmp(argv[i],"-autoFlagDump")==0) {autoFlagDump=1; sscanf(argv[++i],"%f",&autoFlagDump_sigma);}
      else if (strcmp(argv[i],"-h")==0) help();
      else if (strcasecmp(argv[i],"-flag")==0)
	{
	  manual_rfi[nManualRFI].type=1;
	  sscanf(argv[++i],"%lf",&(manual_rfi[nManualRFI].f0));
	  sscanf(argv[++i],"%lf",&(manual_rfi[nManualRFI].f1));
	  nManualRFI++;
	}
      else
	strcpy(fname[nFiles++],argv[i]);
    }
  
  for (ii=0;ii<nFiles;ii++)
    {
      printf("Processing: %s\n",fname[ii]);
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[ii],inFile,1);
      sdhdf_loadMetaData(inFile);
      sdhdf_formOutputFilename(inFile->fname,extension,oname);
      printf(" .... output filename: %s\n",oname);

      // Copy the input file to the output file
      sdhdf_initialiseFile(outFile);
      sdhdf_openFile(oname,outFile,3);

      sdhdf_addHistory(inFile->history,inFile->nHistory,"sdhdf_autoFlag","INSPECTA software to add flags automatically",args);
      inFile->nHistory++;
      sdhdf_writeHistory(outFile,inFile->history,inFile->nHistory);
      
      sdhdf_copyRemainder(inFile,outFile,0);      
      sdhdf_loadMetaData(outFile);

      // Load the transient and persistent RFI
      if (transientFlag==1)
	sdhdf_loadTransientRFI(transient_rfi,&nTransientRFI,MAX_RFI,inFile->primary[0].telescope);
      if (persistentFlag==1)
	sdhdf_loadPersistentRFI(rfi,&nRFI,MAX_RFI,inFile->primary[0].telescope);

      printf("Have loaded %d persistent RFI signals\n",nRFI);
      printf("            %d transient RFI signals\n",nTransientRFI);
      printf("            %d manually entered RFI signals\n",nManualRFI);
      // Now apply the flag tables
      for (b=0;b<inFile->nBeam;b++)
	{
	  haveFlag=0;
	  nband = inFile->beam[b].nBand;
	  for (i=0;i<nband;i++)
	    {
	      // A waste as don't need always need to load in the data: FIX ME
	      // NOT APPLYING TO ALL SPECTRAL DUMPS **** FIX ME
	      sdhdf_loadBandData(inFile,b,i,1);

	      // Apply the persistent RFI
	      if (persistentFlag == 1)
		flagPersistentRFI(inFile,rfi,nRFI,b,i);
	      if (transientFlag == 1)
		flagTransientRFI(inFile,transient_rfi,nTransientRFI,b,i);
	      if (nManualRFI > 0)
		flagPersistentRFI(inFile,manual_rfi,nManualRFI,b,i);
	      if (narrowFlag == 1)
		flagNarrowRFI(inFile,b,i);
	      if (autoFlagDump == 1)
		flagAutoDump(inFile,autoFlagDump_sigma,b,i);
	      if (resetWeights==1)
		{
		  int k;
		  for (k=0;k<inFile->beam[b].bandHeader[i].nchan;k++)
		    inFile->beam[b].bandData[i].astro_data.dataWeights[k] = 1;
		}
	      sdhdf_writeDataWeights(outFile,b,i,inFile->beam[b].bandData[i].astro_data.dataWeights,inFile->beam[b].bandHeader[i].nchan,inFile->beam[b].bandHeader[i].ndump,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
	      sdhdf_writeFlags(outFile,b,i,inFile->beam[b].bandData[i].astro_data.flag,inFile->beam[b].bandHeader[i].nchan,inFile->beam[b].bandHeader[i].ndump,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[i].label);
	      sdhdf_releaseBandData(inFile,b,i,1);
	      printf("Complete band %d\n",i);
	      

	    }
	  printf("Complete beam %d\n",b);
	}
      printf("Closing input file\n");
      sdhdf_closeFile(inFile);
      printf("Closing output file\n");
      sdhdf_closeFile(outFile);	  
      
    }

  free(inFile);
  free(outFile);

  // FREE THE FROM FILE -- FIX ME!
  
}

void flagDataPoint(sdhdf_fileStruct *inFile,int b,int i,int j,int s)
{
  inFile->beam[b].bandData[i].astro_data.flag[j+s*inFile->beam[b].bandHeader[i].nchan]=1;
}

void flagAutoDump(sdhdf_fileStruct *inFile,float autoFlagDump_sigma,int b,int i)
{
  int s,j,k;
  float x,x2;
  float v,v2;
  float val;
  float mean,sdev;
  float meanDump,sdevDump;
  float meanVals[inFile->beam[b].bandHeader[i].ndump];
  int nc=0;
  
  // Calculate the mean and variance
  v=v2=0;
  for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
    {
      x=x2=0;
      nc=0;
      for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
	{
	  if (inFile->beam[b].bandData[i].astro_data.flag[j+s*inFile->beam[b].bandHeader[i].nchan] == 0)
	    {
	      val = (inFile->beam[b].bandData[i].astro_data.pol1[j+s*inFile->beam[b].bandHeader[i].nchan] 
		     + inFile->beam[b].bandData[i].astro_data.pol2[j+s*inFile->beam[b].bandHeader[i].nchan]);
	      
	      x += val;
	      x2 += pow(val,2);
	      nc++;
	    }
	}
      if (nc > 0)
	{
	  mean = x/(float)nc;
	  meanVals[s] = mean;
	  sdev = sqrt(x2/(float)nc - pow(mean,2));
	}
      v += mean;
      v2 += pow(mean,2);
    }
  meanDump = v/(float)inFile->beam[b].bandHeader[i].ndump;
  sdevDump = sqrt(v2/(float)inFile->beam[b].bandHeader[i].ndump - pow(meanDump,2));
  printf("Mean/sdev = %g/%g\n",meanDump,sdevDump);

  // Now determine what (if anything to flag)
  for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
    { 
      if (meanVals[s] > meanDump + autoFlagDump_sigma*sdevDump)
	{
	  printf("Flagging spectral dump %d in band %d because of noise level (deviates by %f sigma)\n",s,i,(meanVals[s]-meanDump)/sdevDump);
	  for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
	    inFile->beam[b].bandData[i].astro_data.flag[j+s*inFile->beam[b].bandHeader[i].nchan] = 1;
	}
    }
}


void flagPersistentRFI(sdhdf_fileStruct *inFile,sdhdf_rfi *rfi,int nRFI,int b,int i)
{
  int k,s,j;
  double freq;
  int haveFlag;
  
  for (k=0;k<nRFI;k++)
    {
      // Is this RFI within the band?
      if ((rfi[k].f0 > inFile->beam[b].bandHeader[i].f0 && rfi[k].f0 < inFile->beam[b].bandHeader[i].f1)
	  || (rfi[k].f1 > inFile->beam[b].bandHeader[i].f0 && rfi[k].f1 < inFile->beam[b].bandHeader[i].f1))
	{
	  // Now zap from each sub-band
	  printf("Processing persistent RFI: %s (%g to %g)\n",rfi[k].description,rfi[k].f0,rfi[k].f1);
	  haveFlag=0;
	  for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
	    {
	      for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
		{
		  freq = inFile->beam[b].bandData[i].astro_data.freq[s*inFile->beam[b].bandHeader[i].nchan+j];
		  if (freq > rfi[k].f0 && freq < rfi[k].f1)
		    {
		      flagDataPoint(inFile,b,i,j,s);
		      haveFlag=1;
		    }
		}				      
	    }
	  if (haveFlag==1)
	    printf(" ... have flagged this RFI\n");
	  
	}
    }
}

void flagTransientRFI(sdhdf_fileStruct *inFile,sdhdf_transient_rfi *transient_rfi,int nTransientRFI,int b,int i)
{
  int j,k,s;
  double v1,v2,v3;
  int n1,n2;
  int flagIt=0;
  double freq;
  double sumPol=0;
  // ii, b, i (band)
  for (k=0;k<nTransientRFI;k++)
    {
      // Is this RFI within the band?
      if ((transient_rfi[k].f0 > inFile->beam[b].bandHeader[i].f0 && transient_rfi[k].f0 < inFile->beam[b].bandHeader[i].f1)
	  || (transient_rfi[k].f1 > inFile->beam[b].bandHeader[i].f0 && transient_rfi[k].f1 < inFile->beam[b].bandHeader[i].f1))
	{
	  // Now zap from each sub-band
	  printf("Checking transient RFI: %s comparing %g -> %g with %g -> %g, type = %d\n",transient_rfi[k].description,
		 transient_rfi[k].f0,transient_rfi[k].f1,transient_rfi[k].f2,transient_rfi[k].f3,transient_rfi[k].type);
	  for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
	    {
	      v1 = v2 = v3 = 0;
	      
	      n1=n2=0;
	      flagIt=0;
	      for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
		{
		  freq = inFile->beam[b].bandData[i].astro_data.freq[s*inFile->beam[b].bandHeader[i].nchan+j];
		  if (freq > transient_rfi[k].f0 && freq < transient_rfi[k].f1)
		    {
		      sumPol = (inFile->beam[b].bandData[i].astro_data.pol1[j+s*inFile->beam[b].bandHeader[i].nchan] 
				+ inFile->beam[b].bandData[i].astro_data.pol2[j+s*inFile->beam[b].bandHeader[i].nchan]);
		      if (transient_rfi[k].type == 1 && v1 < sumPol) v1 = sumPol;
		      else if (transient_rfi[k].type == 2) v1 += sumPol;
		      n1 ++;
		    }
		  if (freq > transient_rfi[k].f2 && freq < transient_rfi[k].f3)
		    {
		      sumPol = (inFile->beam[b].bandData[i].astro_data.pol1[j+s*inFile->beam[b].bandHeader[i].nchan] 
				+ inFile->beam[b].bandData[i].astro_data.pol2[j+s*inFile->beam[b].bandHeader[i].nchan]);
		      v2 +=  sumPol;
		      n2 ++;
		    }
		}
	      if (transient_rfi[k].type == 2) v1/=(double)n1;
	      v2/=(double)n2; // SHOULD CHECK n2 > 0 --- FIX ME
	      
	      
	      if (v1 > transient_rfi[k].threshold*v2)
		{
		  flagIt=1;
		  printf("Flagging dump %d, comparing %g and %g (%d %d)\n",s,v1,v2,n1,n2);
		}
	      else
		printf("**Not** flagging dump %d, comparing %g and %g (%d %d)\n",s,v1,v2,n1,n2);
	      if (flagIt==1) // Do the flagging if necessary
		{
		
		  for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
		    {
		      // FIX ME: using [0] for frequency dump
		      freq = inFile->beam[b].bandData[i].astro_data.freq[s*inFile->beam[b].bandHeader[i].nchan+j];
		      if (freq > transient_rfi[k].f0 && freq < transient_rfi[k].f1)
			{
			  flagDataPoint(inFile,b,i,j,s);
			}
		    }					  
		}
	    }
	}
    }
}

void flagNarrowRFI(sdhdf_fileStruct *inFile,int b,int i)
{
  int j,k,s;
  int j0,j1,j2,j3;
  double v1,v2,v3;
  double mean,sdev;
  int n1,n2;
  int flagIt=0;
  double freq;
  double sumPol=0;
  int blSamples  = 4;
  int rfiSamples = 1;
  // ii, b, i (band)



  printf("Flagging narrow spectral features\n");
  

  for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
    {
      
      n1=n2=0;
      flagIt=0;
      for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
	{
	  sumPol = (inFile->beam[b].bandData[i].astro_data.pol1[j+s*inFile->beam[b].bandHeader[i].nchan] 
		    + inFile->beam[b].bandData[i].astro_data.pol2[j+s*inFile->beam[b].bandHeader[i].nchan]);
	  v1 = sumPol;
	  j0 = j-blSamples;
	  if (j0 < 0) j0 = 0;
	  j1 = j-1;
	  if (j1 < 0) j1 = 0;
	  j2 = j+rfiSamples;
	  j3 = j2+blSamples;
	  if (j2 > inFile->beam[b].bandHeader[i].nchan) j2 = inFile->beam[b].bandHeader[i].nchan;
	  if (j3 > inFile->beam[b].bandHeader[i].nchan) j3 = inFile->beam[b].bandHeader[i].nchan;
	  
	  // Need to account for flagged channels ** FIX ME
	  n1 = n2 = 0;
	  v2 = v3 = 0;

	  for (k=j0;k<j1;k++)
	    {
	      sumPol = (inFile->beam[b].bandData[i].astro_data.pol1[k+s*inFile->beam[b].bandHeader[i].nchan] 
			+ inFile->beam[b].bandData[i].astro_data.pol2[k+s*inFile->beam[b].bandHeader[i].nchan]);
	      v2 += sumPol;
	      v3 += pow(sumPol,2);
	      n1++;
	    }
	  for (k=j2;k<j3;k++)
	    {
	      sumPol = (inFile->beam[b].bandData[i].astro_data.pol1[k+s*inFile->beam[b].bandHeader[i].nchan] 
			+ inFile->beam[b].bandData[i].astro_data.pol2[k+s*inFile->beam[b].bandHeader[i].nchan]);
	      v2 += sumPol;
	      v3 += pow(sumPol,2);
	      n1++;
	    }
	  mean = v2/(double)n1;
	  sdev = sqrt(1./(double)n1*v3 - pow(1.0/(double)n1 * v2,2));
	  //	  printf("CHECKING NARROW: %g %g %g %d (V2 = %g, %g)\n",v1,mean,sdev,n1,v2,v3);
		
	  if (v1 > mean + 5*sdev) // FIX SDEV
	    flagDataPoint(inFile,b,i,j,s);	    
	}
    }
}
