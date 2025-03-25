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
// Software to output text files relating to the data sets in the file
//
// Usage:
// sdhdf_quickdump <filename.hdf>  <filename.hdf> ...
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"

#define VNUM "v2.0"

void help()
{
  printf("\nsdhdf_quickdump       %s\n",VNUM);
  printf("INSPECTA version:     %s\n",SOFTWARE_VER);
  printf("Author:               George Hobbs\n");
  printf("Software to output text files of the data in SDHDF files\n");

  printf("\nCommand line arguments:\n\n");
  printf("-h                     This help\n");
  printf("-bandRange <b0> <b1>   Output data in bands between b0 and b1\n");
  printf("-cal                   Output calibration data\n");
  printf("-cal_tav               Average calibration data in time\n");
  printf("-cal32_ch <ch>         Set the channel to 'ch' when outputting data for the 32-bin cal\n");
  printf("-cal32_tav             Average the 32-binned calibration signal in time before output\n");
  printf("-calOnOff              Calculate ON-OFF for the noise source parameters\n");
  printf("-dumpRange <s0> <s1>   Output data in specificed sub-integration range\n");
  printf("-e <ext>               Output to data file with defined extension\n");
  printf("-fref <fref> Use 'fref' as a reference frequency and output velocities as well as frequencies\n");
  printf("-freqRange <f0> <f1>   Output data between f0 and f1 (in MHz)\n");
  printf("-mjd                   Print dump time MJD in output\n");
  printf("-noflagged             Do not print out lines that have been flagged\n");
  printf("-stokes                Convert coherency products to Stokes (for calibrated data) or pseudo-Stokes (uncalibrated data)\n");
  printf("-tsys                  output system temperature measurements\n");

  printf("\nExample:\n\n");
  printf("sdhdf_quickdump file.hdf\n\n");

  exit(1);
}

void readPhaseResolvedCal(sdhdf_fileStruct *inFile,int band,int cal32_chN,int cal32_tav);

int main(int argc,char *argv[])
{
  int i,j,k,beam,band;
  int mjd=0;
  int display=0;
  int notFlagged=0;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  float freq;
  int   flag;
  float weight;
  float pol1,pol2,pol3,pol4;
  sdhdf_fileStruct *inFile;
  int npol=0;
  int stokes=0;
  int outFile=0;
  float *tsys;
  char extFileName[MAX_STRLEN];
  char outFileName[MAX_STRLEN];
  FILE *fout;
  float freq0,freq1;
  int   setFreqRange=0;
  int sd0,sd1;
  int   setDumpRange=0;
  int nchan;
  int dataType=1;
  int calOnOff = 0;
  int cal32_chN=-1;
  int cal32_tav=-1;
  int band0=-1,band1=-1;
  double fref=-1;
  int cal_tav=-1;
  
  if (argc==1)
    help();
	
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-mjd")==0)
	mjd=1;
      else if (strcasecmp(argv[i],"-notFlagged")==0)
	notFlagged=1;
      else if (strcmp(argv[i],"-stokes")==0)
	stokes=1;
      else if (strcmp(argv[i],"-tsys")==0)
	dataType=2;
      else if (strcmp(argv[i],"-fref")==0)
	sscanf(argv[++i],"%lf",&fref);
      else if (strcmp(argv[i],"-cal")==0)
	dataType=3;
      else if (strcmp(argv[i],"-cal_tav")==0)
	cal_tav=1;
      else if (strcmp(argv[i],"-calOnOff")==0)
	calOnOff=1;
      else if (strcmp(argv[i],"-cal32")==0)
	dataType=4;
      else if (strcmp(argv[i],"-cal32_ch")==0)
	sscanf(argv[++i],"%d",&cal32_chN);
      else if (strcmp(argv[i],"-cal32_tav")==0)
	cal32_tav=1;
      else if (strcmp(argv[i],"-freqRange")==0)
	{
	  sscanf(argv[++i],"%f",&freq0);
	  sscanf(argv[++i],"%f",&freq1);
	  setFreqRange=1;
	}
      else if (strcmp(argv[i],"-dumpRange")==0)
	{
	  sscanf(argv[++i],"%d",&sd0);
	  sscanf(argv[++i],"%d",&sd1);
	  setDumpRange=1;
	}
      else if (strcmp(argv[i],"-bandRange")==0)
	{
	  sscanf(argv[++i],"%d",&band0);
	  sscanf(argv[++i],"%d",&band1);
	  setDumpRange=1;
	}
      else if (strcmp(argv[i],"-e")==0)
	{
	  strcpy(extFileName,argv[++i]);
	  outFile=1;
	}
      else
	strcpy(fname[nFiles++],argv[i]);
    }

  //  printf("Starting quickdump\n");
  
  inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));
  
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      if (sdhdf_openFile(fname[i],inFile,1)==-1)
	printf("Warning: unable to open file >%s<. Skipping\n",fname[i]);
      else
	{
	  sdhdf_loadMetaData(inFile);
	  //	  printf("Have loaded metadata\n");
	  if (outFile==1)
	    {
	      sprintf(outFileName,"%s.%s",inFile->fname,extFileName);
	      if (!(fout = fopen(outFileName,"w")))
		{
		  printf("ERROR: Unable to open file >%s< for writing\n",outFileName);
		  exit(1);
		}
	    }
	  
	  for (beam=0;beam<inFile->nBeam;beam++)
	    {
	      if (band0 < 0)
		{
		  band0 = 0;
		  band1 = inFile->beam[beam].nBand;
		}
	      for (band=band0;band<band1;band++)
		{
		  npol = inFile->beam[beam].bandHeader[band].npol;
		  nchan = inFile->beam[beam].bandHeader[band].nchan;
		  //		  printf("Loading the data\n");
		  if (dataType!=3 && dataType!=4)
		    sdhdf_loadBandData(inFile,beam,band,1);
		  //		  printf("Loading calibration data\n");
		  if (dataType==2 || dataType==3)
		    {
		      sdhdf_loadBandData(inFile,beam,band,2);
		      sdhdf_loadBandData(inFile,beam,band,3);
		    }
		  //		  printf("Complete loading the data\n");

		  if (dataType==1) // Astronomy data
		    {
		      if (setDumpRange==0)
			{
			  sd0 = 0;
			  sd1 = inFile->beam[beam].bandHeader[band].ndump;
			}
		      
		      for (k=sd0;k<sd1;k++)
			{
			  //		      printf("Here with %d\n",k);
			  for (j=0;j<nchan;j++)
			    {
			      display=1;
			      //			  printf("Loading freq\n");
			      freq = inFile->beam[beam].bandData[band].astro_data.freq[j];// FIX ME: DUMP FOR FREQ AXIS
			      //			      printf("Loading pol1 %g\n",freq);
			      pol1 = inFile->beam[beam].bandData[band].astro_data.pol1[j+k*nchan];

			      if (setFreqRange==1 && (freq < freq0 || freq > freq1))
				display=0;
			      
			      if (npol > 1)
				{
				  pol2 = inFile->beam[beam].bandData[band].astro_data.pol2[j+k*nchan];
				  if (npol > 2)
				    {
				      pol3 = inFile->beam[beam].bandData[band].astro_data.pol3[j+k*nchan];
				      pol4 = inFile->beam[beam].bandData[band].astro_data.pol4[j+k*nchan];
				    }
				  else
				    pol3 = pol4 = 0;
				}
			      else
				pol2 = pol3 = pol4 = 0;

			      if (stokes==1)
				{
				  double sI,sQ,sU,sV;
				  sI = pol1 + pol2;
				  sQ = pol1 - pol2;
				  sU = 2*pol3;
				  sV = 2*pol4;

				  pol1 = sI;
				  pol2 = sQ;
				  pol3 = sU;
				  pol4 = sV;
				}
			      
			      weight = inFile->beam[beam].bandData[band].astro_data.dataWeights[k*nchan+j]*(1-inFile->beam[beam].bandData[band].astro_data.flag[k*nchan+j]);
			      if (notFlagged == 1 && weight == 0)
				display=0;
			      
			      if (display==1)
				{
				  if (outFile==1)
				    {
				      if (fref < 0)
					fprintf(fout,"%s %d %d %d %d %.6f %g %g %g %g %g\n",inFile->fname,beam,band,k,j,freq,pol1,pol2,pol3,pol4,weight);
				      else
					fprintf(fout,"%s %d %d %d %d %.6f %g %g %g %g %g %.6g\n",inFile->fname,beam,band,k,j,freq,pol1,pol2,pol3,pol4,weight,(1.0-freq/fref)*SPEED_LIGHT/1000.); // km/s);
				    }
				  else
				    {
				      if (fref < 0)
					printf("%s %d %d %d %d %.6f %g %g %g %g %g %.6f\n",inFile->fname,beam,band,k,j,freq,pol1,pol2,pol3,pol4,weight,
					       inFile->beam[beam].bandData[0].astro_obsHeader[k].mjd);
				      else
					printf("%s %d %d %d %d %.6f %g %g %g %g %g %.6f\n",inFile->fname,beam,band,k,j,freq,pol1,pol2,pol3,pol4,weight,(1.0-freq/fref)*SPEED_LIGHT/1000.); // km/s);
				    }
				}
			    }
			}
		    }
		  else if (dataType==2)
		    {
		      tsys = (float *)malloc(sizeof(float)*inFile->beam[beam].calBandHeader[band].ndump*inFile->beam[beam].calBandHeader[band].nchan*2); // * 2 = AA and BB
		      sdhdf_loadCalProc(inFile,beam,band,"cal_proc_tsys",tsys);
		      for (k=0;k<inFile->beam[beam].calBandHeader[band].ndump;k++)
			{
			  //		      printf("Here with %d\n",k);
			  for (j=0;j<inFile->beam[beam].calBandHeader[band].nchan;j++)
			    {
			      // FIX ME: DUMP FOR FREQ AXIS
			      if (outFile==1)
				fprintf(fout,"%s %d %d %d %d %.6f %g %g\n",inFile->fname,beam,band,k,j,inFile->beam[beam].bandData[band].cal_on_data.freq[j],
					tsys[k*2*inFile->beam[beam].calBandHeader[band].nchan + j],
					tsys[k*2*inFile->beam[beam].calBandHeader[band].nchan+inFile->beam[beam].calBandHeader[band].nchan + j]);
			      else
				printf("%s %d %d %d %d %.6f %g %g\n",inFile->fname,beam,band,k,j,inFile->beam[beam].bandData[band].cal_on_data.freq[j],
					tsys[k*2*inFile->beam[beam].calBandHeader[band].nchan + j],
					tsys[k*2*inFile->beam[beam].calBandHeader[band].nchan+inFile->beam[beam].calBandHeader[band].nchan + j]);

			    }
			  if (outFile==1)
			    fprintf(fout,"\n");
			}
		      free(tsys);
		    }
		  else if (dataType==3)
		    {
		      int nchanCal =  inFile->beam[beam].calBandHeader[band].nchan;
		      if (setDumpRange==0)
			{
			  sd0 = 0;
			  sd1 = inFile->beam[beam].calBandHeader[band].ndump;
			}
		      if (cal_tav==1)
			{
			  double cal_on_pol1,cal_off_pol1;
			  double cal_on_pol2,cal_off_pol2;
			  double cal_on_pol3,cal_off_pol3;
			  double cal_on_pol4,cal_off_pol4;
			  double mjdAv;
			  for (k=0;k<nchanCal;k++)
			    {
			      display=1;
			      //			  printf("Loading freq\n");
			      // FIX ME: DUMP FOR FREQ AXIS
			      freq = inFile->beam[beam].bandData[band].cal_on_data.freq[k];
			      if (setFreqRange==1 && (freq < freq0 || freq > freq1))
				display=0;
			      if (display==1)
				{
				  int nd;
				  cal_on_pol1=cal_off_pol1=0;
				  cal_on_pol2=cal_off_pol2=0;
				  cal_on_pol3=cal_off_pol3=0;
				  cal_on_pol4=cal_off_pol4=0;
				  mjdAv = 0;
				  nd=0;
				  for (j=sd0;j<sd1;j++)
				    {				      
				      cal_on_pol1+=inFile->beam[beam].bandData[band].cal_on_data.pol1[k+j*nchanCal];
				      cal_off_pol1+=inFile->beam[beam].bandData[band].cal_off_data.pol1[k+j*nchanCal];

				      cal_on_pol2+=inFile->beam[beam].bandData[band].cal_on_data.pol2[k+j*nchanCal];
				      cal_off_pol2+=inFile->beam[beam].bandData[band].cal_off_data.pol2[k+j*nchanCal];

				      cal_on_pol3+=inFile->beam[beam].bandData[band].cal_on_data.pol3[k+j*nchanCal];
				      cal_off_pol3+=inFile->beam[beam].bandData[band].cal_off_data.pol3[k+j*nchanCal];

				      cal_on_pol4+=inFile->beam[beam].bandData[band].cal_on_data.pol4[k+j*nchanCal];
				      cal_off_pol4+=inFile->beam[beam].bandData[band].cal_off_data.pol4[k+j*nchanCal];
				      mjdAv += inFile->beam[beam].bandData[band].cal_obsHeader[j].mjd;
				      nd++;
				    }
				  cal_on_pol1/=(double)nd;  cal_off_pol1/=(double)nd;
				  cal_on_pol2/=(double)nd;  cal_off_pol2/=(double)nd;
				  cal_on_pol3/=(double)nd;  cal_off_pol3/=(double)nd;
				  cal_on_pol4/=(double)nd;  cal_off_pol4/=(double)nd;
				  mjdAv /= (double)nd;
				  /*
				  printf("%s %d %d %d %d %.6f %g %g %g %g %g %g %g %g %.6f %.6f %.6f %.6f %s\n",inFile->fname,beam,band,k,0,freq,
					 cal_on_pol1,cal_off_pol1,cal_on_pol2,cal_off_pol2,cal_on_pol3,cal_off_pol3,cal_on_pol4,cal_off_pol4,
					 inFile->beam[beam].bandData[band].cal_obsHeader[j].mjd,inFile->beam[beam].bandData[band].cal_obsHeader[j].az,inFile->beam[beam].bandData[band].cal_obsHeader[j].el,inFile->beam[beam].bandData[band].cal_obsHeader[j].paraAngle,inFile->beamHeader[beam].source);
				  */
				  if (calOnOff==0)
				    printf("%s %d %d %d %d %g %g %g %g %g %g %g %g %g %.6f\n",inFile->fname,beam,band,k,0,freq,cal_on_pol1,cal_off_pol1
					   ,cal_on_pol2,cal_off_pol2,cal_on_pol3,cal_off_pol3,cal_on_pol4,cal_off_pol4,mjdAv);
				  else
				    printf("%s %d %d %d %d %g %g %g %g %g %.5f\n",inFile->fname,beam,band,k,0,freq,cal_on_pol1-cal_off_pol1
					   ,cal_on_pol2-cal_off_pol2,cal_on_pol3-cal_off_pol3,cal_on_pol4-cal_off_pol4,mjdAv);
				}
			    }

			}
		      else
			{
			  for (j=sd0;j<sd1;j++)
			    {
			      //		      printf("Here with %d\n",k);
			      for (k=0;k<nchanCal;k++)
				{
				  display=1;
				  //			  printf("Loading freq\n");
				  // FIX ME: DUMP FOR FREQ AXIS
				  freq = inFile->beam[beam].bandData[band].cal_on_data.freq[k];
				  if (setFreqRange==1 && (freq < freq0 || freq > freq1))
				    display=0;
				  if (display==1)
				    printf("%s %d %d %d %d %.6f %g %g %g %g %g %g %g %g %.6f %.6f %.6f %.6f %s\n",inFile->fname,beam,band,k,j,freq,inFile->beam[beam].bandData[band].cal_on_data.pol1[k+j*nchanCal],inFile->beam[beam].bandData[band].cal_off_data.pol1[k+j*nchanCal],inFile->beam[beam].bandData[band].cal_on_data.pol2[k+j*nchanCal],inFile->beam[beam].bandData[band].cal_off_data.pol2[k+j*nchanCal],inFile->beam[beam].bandData[band].cal_on_data.pol3[k+j*nchanCal],inFile->beam[beam].bandData[band].cal_off_data.pol3[k+j*nchanCal],inFile->beam[beam].bandData[band].cal_on_data.pol4[k+j*nchanCal],inFile->beam[beam].bandData[band].cal_off_data.pol4[k+j*nchanCal],inFile->beam[beam].bandData[band].cal_obsHeader[j].mjd,inFile->beam[beam].bandData[band].cal_obsHeader[j].az,inFile->beam[beam].bandData[band].cal_obsHeader[j].el,inFile->beam[beam].bandData[band].cal_obsHeader[j].paraAngle,inFile->beamHeader[beam].source);
				  
				}
			    }
			}
		    }
		  else if (dataType==4)
		    {
		      readPhaseResolvedCal(inFile,band,cal32_chN,cal32_tav);
		    }
		  
		  
		  sdhdf_releaseBandData(inFile,beam,band,1);
		  
		}
	    }
	  if (outFile==1)
	    fclose(fout);
	  sdhdf_closeFile(inFile);
	}
    }

  free(inFile);
}


void readPhaseResolvedCal(sdhdf_fileStruct *inFile,int band,int cal32_chN,int cal32_tav)
{
  int i,j,k,b;
  char dataName[1024];
  int beam=0; // FIX ME
  hid_t dataset_id;
  herr_t status;
  int ndump = inFile->beam[beam].calBandHeader[band].ndump;
  int npol = 4; // FIX ME
  int nchan = inFile->beam[beam].calBandHeader[band].nchan;
  int nbin = 32; // FIX ME
  float *data;
  int k0,k1;

  
  if (cal32_chN < 0)
    {
      k0 = 0;
      k1 = nchan;
    }
  else
    {
      k0 = cal32_chN;
      k1 = k0+1;
    }
  sprintf(dataName,"beam_%d/%s/calibrator_data/cal32_data",beam,inFile->beam[beam].bandHeader[band].label);
  dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
  data = (float *)malloc(sizeof(float)*nchan*npol*ndump*nbin);
  status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);  
  if (cal32_tav == 1)
    {
      float bin1[nbin],bin2[nbin],bin3[nbin],bin4[nbin];
      for (k=k0;k<k1;k++)
	{
	  for (b=0;b<nbin;b++)
	    bin1[b] = bin2[b] = bin3[b] = bin4[b] = 0;
	  for (i=0;i<ndump;i++)
	    {
	      for (b=0;b<nbin;b++)
		{
		  bin1[b] += data[i*npol*nchan*nbin + 0*nchan*nbin +k*nbin + b];
		  bin2[b] += data[i*npol*nchan*nbin + 1*nchan*nbin +k*nbin + b];
		  bin3[b] += data[i*npol*nchan*nbin + 2*nchan*nbin +k*nbin + b];
		  bin4[b] += data[i*npol*nchan*nbin + 3*nchan*nbin +k*nbin + b];
		}
	    }
	  for (b=0;b<nbin;b++)
	    {
	      bin1[b]/=(double)ndump;
	      bin2[b]/=(double)ndump;
	      bin3[b]/=(double)ndump;
	      bin4[b]/=(double)ndump;
	      printf("%d %d %d %g %g %g %g GEORGE\n",0,k,b,
		     bin1[b],bin2[b],bin3[b],bin4[b]);
	    }
	    }
      printf("\n");
      
    }
  else
    {
      for (i=0;i<ndump;i++)
	{
	  for (k=k0;k<k1;k++)
	    {
	      for (b=0;b<nbin;b++)
		printf("%d %d %d %g %g %g %g GEORGE\n",i,k,b,
		       data[i*npol*nchan*nbin + 0*nchan*nbin +k*nbin + b],
		       data[i*npol*nchan*nbin + 1*nchan*nbin +k*nbin + b],
		       data[i*npol*nchan*nbin + 2*nchan*nbin +k*nbin + b],
		       data[i*npol*nchan*nbin + 3*nchan*nbin +k*nbin + b]);
	    }
	  printf("\n");
	}
  }

  status = H5Dclose(dataset_id);
  free(data);
  //  exit(1);
}
