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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "inspecta.h"
#include "hdf5.h"

#define VNUM "v2.0"

//
// Code to convert from filterbank and PSRFITS format to SDHDF
//


// Defined from SIGPROC header.h
// define the signedness for the 8-bit data type
#define OSIGN 1
#define SIGNED OSIGN < 0



typedef struct header {
  //
  // Filterbank header information

  char fil_isign;
  int fil_telescope_id;
  int fil_machine_id;
  int fil_data_type;
  char fil_rawdatafile[80];
  char fil_source_name[80];
  int fil_barycentric;
  int fil_pulsarcentric;
  double fil_az_start;
  double fil_za_start;
  double fil_src_raj;
  double fil_src_decj;
  double fil_tstart;
  double fil_tsamp;
  int fil_nbits;
  int fil_nsamples;
  double fil_fch1;
  double fil_foff;
  char fil_freqStart;
  double fil_fchannel;
  char fil_freqEnd;
  int fil_nchans;
  int fil_nifs;
  double fil_refdm;
  double fil_period;
  int fil_ibeam;
  int fil_nbeams;
  int fil_nbins;
  int fil_npuls;

  // 
  char outputName[1024];
  int  nsblk;
} header;

void help()
{
  printf("\nsdhdf_convertFromPsr  %s\n",VNUM);
	printf("INSPECTA version:     %s\n",SOFTWARE_VER);
  printf("Author:               George Hobbs\n");
  printf("Software to convert a filterbank or PSRFITS file to SDHDF format\n");

  printf("\nCommand line arguments:\n\n");
	printf("-h                    This help\n");
	printf("-f <filename>         Filename to convert\n");
  printf("-cal                  \n");
  printf("-calOn					      \n");
	printf("-calOff					      \n");
  printf("-pulseOff		  	      \n");
  printf("-bin                  \n");
  printf("-type                 \n");
  printf("-nsblk                \n");
	printf("-o <filename>         Output filename\n");

	printf("\nExample:\n\n");
  printf("sdhdf_convertFromPsr -f file.hdf\n\n");

  exit(1);
}

void readFilHeader(char *fname,FILE *inputfile,header *head);
void get_string(FILE *fin,int *nbytes,char string[]);
int strings_equal (char *string1, char *string2);
void bytesToFloats(int samplesperbyte,int n,unsigned char *cVals,float *out);
void eightBitsFloat(int eight_bit_number, float *results, int *index);
void fourBitsFloat(int eight_bit_number, float *results, int *index);
void twoBitsFloat(int eight_bit_number, float *results, int *index);
void twoBitsFloat_unsigned(int eight_bit_number, float *results, int *index);
void oneBitFloat(int eight_bit_number, float *results, int *index);

int main(int argc,char *argv[])
{
  int i,j,k,b;

  int fileType=1; // 1 = fil, 2 = PSRFITS
  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;

  
  sdhdf_fileStruct *outFile;
  sdhdf_softwareVersionsStruct *softwareVersions;
  sdhdf_historyStruct *history;
  sdhdf_primaryHeaderStruct *primaryHeader;
  sdhdf_beamHeaderStruct *beamHeader;
  sdhdf_bandHeaderStruct *bandHeader;
  sdhdf_obsParamsStruct  *obsParams;
  sdhdf_obsParamsStruct  *calObsParams;
  sdhdf_bandHeaderStruct *calBandHeader;
  
  
  // Input FITS file
  int status=0;
  fitsfile *fptr;
  int maxdim = 5;
  int naxis;
  long naxes[maxdim];
  int type;
  long nrows;
  long maxsize;
  int    nchan;
  int    nchan_cal;
  int    npol;
  int    npol_cal;
  int    ndump;
  int    ndump_cal;
  long   long_ndump;
  double freq0;
  double chan_bw;
  float  *floatVals;
  float  *floatCalValsOn;
  float  *floatCalValsOff;
  float  *infloatVals;
  float  *freqVals;
  double *dfreqVals;
  float  *freqCalVals;
  double *dfreqCalVals;
  float *datOffs;
  float *datScl;
  double *doubleVals;
  short int *shortVals;
  short int n_sval=0;
  int    binVal = 0; 
  int    nVals;
  float  n_fval=0.;
  double n_dval=0.;
  long   n_lval=0;
  int    colnum;
  int    colnum_datoffs;
  int    colnum_datscl;

  int    initflag=0;
  double mean=0;

  char fname[1024];
  char outname[1024]="convert.hdf";

  int cal=0;
  char calFile[1024];
  int  calOn1,calOn2;
  int  calOff1,calOff2;

  double tdump,mjd;
  char srcName[1024];
  int nsblk = 2048;
  int pulseOff1,pulseOff2;
  double *baseline1,*baseline2,*baseline3,*baseline4;
  pulseOff1=pulseOff2=-1;
  
  primaryHeader    = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct));
  beamHeader       = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct));
  bandHeader       = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct));
  softwareVersions = (sdhdf_softwareVersionsStruct *)malloc(sizeof(sdhdf_softwareVersionsStruct));
  history          = (sdhdf_historyStruct *)malloc(sizeof(sdhdf_historyStruct));
  
  sdhdf_setMetadataDefaults(primaryHeader,beamHeader,bandHeader,softwareVersions,history,1,1);

  if (argc==1)
    help();
    
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-cal")==0)
	{
	  cal=1;
	  strcpy(calFile,argv[++i]);
	}
      else if (strcmp(argv[i],"-calOn")==0)
	{
	  sscanf(argv[++i],"%d",&calOn1);
	  sscanf(argv[++i],"%d",&calOn2);
	}
      else if (strcmp(argv[i],"-calOff")==0)
	{
	  sscanf(argv[++i],"%d",&calOff1);
	  sscanf(argv[++i],"%d",&calOff2);
	}
      else if (strcmp(argv[i],"-pulseOff")==0)
	{
	  sscanf(argv[++i],"%d",&pulseOff1);
	  sscanf(argv[++i],"%d",&pulseOff2);
	}
      else if (strcmp(argv[i],"-bin")==0)
	sscanf(argv[++i],"%d",&binVal);
      else if (strcmp(argv[i],"-type")==0)
	sscanf(argv[++i],"%d",&fileType);
      else if (strcmp(argv[i],"-nsblk")==0)
	sscanf(argv[++i],"%d",&nsblk);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outname,argv[++i]);
    }
  
  

  if (fileType==1) // Filterbank
    {    
      FILE *fin;
      header head;
      unsigned char *loadit; 
      float *fvals;
      int nread;
      int nbits;

      
      printf("Opening file >%s<\n",fname);
      if (!(fin = fopen(fname,"rb")))
	{
	  printf("Unable to open file %s\n",fname);
	  exit(1);
	}
      readFilHeader(fname,fin,&head);
      nbits = head.fil_nbits;
      printf("nbits = %d\n",nbits);
      nchan = head.fil_nchans;
      npol  = 1; // SHOULD CHECK
      loadit = (unsigned char *)malloc(sizeof(unsigned char)*nchan*nbits*nsblk/8);
      fvals = (float *)malloc(sizeof(float)*(long)nchan*(long)nsblk*(long)npol);
      ndump = 0; // nsamples isn't properly in the headers
      //      printf("Allocated space for %ld %d %d %d\n",(long)(nchan*npol*max_dumps),nchan,npol,max_dumps);
      if (!(floatVals   = (float *)calloc(sizeof(float),nchan*npol)))
	{
	  printf("Unable to allocate enough memory\n");
	  exit(1);
	}
      // Must now unpack the bits and store as floating point values
      printf("Number of dumps = %d %d\n",ndump,head.fil_nsamples);
      do
	{
	  if (nbits <= 8)
	    {
	      nread = fread(loadit,1,nbits*nchan*nsblk/8,fin); 
	      bytesToFloats(8/nbits,nchan*nsblk*npol,loadit,fvals);
	    }
	  else
	    nread = fread(fvals,sizeof(float),nchan,fin); 
	  if (nread > 0)
	    {
	      printf("%d Have read %d values, sizeof(float) = %d\n",(int)ndump,(int)nread,sizeof(float));
	      
	      printf("Complete conversion %d %d %d\n",nchan,npol,nsblk);
	      for (i=0;i<nchan;i++)
		{
		  for (j=0;j<npol;j++)
		    {
		      for (k=0;k<nsblk;k++)
			{
			  //			  printf("a %d %d %d %g %d\n",i,j,k,fvals[j*nchan+k*nchan*npol+i],ndump*nchan*npol+j*nchan+i);
			  floatVals[j*nchan + i] += fvals[j*nchan+k*nchan*npol+i];
			  //			  printf("b\n");
			}
		    }
		}
	      // Now should write out this dump
	      // ****** GOT TO HERE **********  CAN'T LOAD THE LOT IN ONE GO
	      
	      //	      printf("Got to here\n");
	      ndump++;
	    }
	  //	  if (ndump == 10)
	  //	    break;
	}
      while (nread > 0);
      printf("Loaded %d time intervals\n",ndump);
      printf("Complete conversion\n");

      freqVals    = (float *)malloc(sizeof(float)*nchan);
      dfreqVals    = (double *)malloc(sizeof(double)*nchan);
      for (i=0;i<nchan;i++)
	freqVals[i] = head.fil_fch1+i*head.fil_foff;
      tdump = head.fil_tsamp*nsblk;
      mjd   = head.fil_tstart;
      strcpy(srcName,head.fil_source_name);
      printf("mjd = %f\n",mjd);
      fclose(fin);
      free(loadit);
      free(fvals);
    }
  else if (fileType==2)
    {
      int nbin;
      int nbin_cal;
      int nOn,nOff;
      if (cal==1)
	{
	  calBandHeader  = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct));
	  calObsParams   = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);

	  strcpy(primaryHeader->cal_mode,"ON");
	  
	  
	  // First process cal
	  printf("Opening file >%s<\n",calFile);
	  fits_open_file(&fptr,calFile,READONLY,&status);
	  fits_report_error(stderr,status);      

	  printf("Reading PSRFITS file\n");
	  fits_movnam_hdu(fptr,BINARY_TBL,"SUBINT",1,&status);
	  fits_get_num_rows(fptr,&long_ndump,&status);
	  ndump_cal = (int)long_ndump;
	  
	  doubleVals = (double *)malloc(sizeof(double)*ndump_cal);
	  
	  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
	  fits_read_tdim(fptr,colnum,maxdim,&naxis,naxes,&status);
	  fits_report_error(stderr,status);
	  
	  nbin_cal  = naxes[0];
	  nchan_cal = naxes[1];
	  npol_cal  = naxes[2];
	  printf("nbin = %d, nchan = %d, npol = %d, status = %d\n",nbin_cal,nchan_cal,npol_cal,status);
	  
	  floatCalValsOn  = (float *)malloc(sizeof(float)*nchan_cal*npol_cal*ndump_cal);
	  floatCalValsOff = (float *)malloc(sizeof(float)*nchan_cal*npol_cal*ndump_cal);
	  freqCalVals    = (float *)malloc(sizeof(float)*nchan_cal);
	  dfreqCalVals    = (double *)malloc(sizeof(double)*nchan_cal);
	  shortVals   = (short int *)malloc(sizeof(short int)*nchan_cal*npol_cal*nbin_cal);
	  
	  fits_get_colnum(fptr,CASEINSEN,"DAT_OFFS",&colnum_datoffs,&status);
	  fits_get_colnum(fptr,CASEINSEN,"DAT_SCL",&colnum_datscl,&status);
	  datOffs = (float *)malloc(sizeof(float)*nchan_cal*npol_cal);
	  datScl = (float *)malloc(sizeof(float)*nchan_cal*npol_cal);	

	  
	  fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&colnum,&status);
	  fits_read_col(fptr,TFLOAT,colnum,1,1,nchan_cal,&n_fval,freqCalVals,&initflag,&status);


	  strcpy(calBandHeader->label,"band0");
	  calBandHeader->fc = (freqCalVals[0]+freqCalVals[nchan_cal-1])/2.0; 
	  calBandHeader->f0 = freqCalVals[0];
	  calBandHeader->f1 = freqCalVals[nchan_cal-1]; 
	  calBandHeader->nchan = nchan_cal;
	  calBandHeader->npol = npol_cal;
	  strcpy(calBandHeader->pol_type,"AABBCRCI");
	  calBandHeader->dtime = 100; // FIX
	  calBandHeader->ndump = ndump_cal;
	  
	  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
	  for (i=0;i<ndump_cal;i++)
	    {
	      fits_read_col(fptr,TFLOAT,colnum_datoffs,i+1,1,nchan_cal*npol_cal,&n_fval,datOffs,&initflag,&status);
	      fits_read_col(fptr,TFLOAT,colnum_datscl,i+1,1,nchan_cal*npol_cal,&n_fval,datScl,&initflag,&status);


	      calObsParams[i].timeElapsed = i; // FIX
	      strcpy(calObsParams[i].timedb,"UNKNOWN"); // FIX
	      calObsParams[i].mjd = 56000; // FIX
	      calObsParams[i].fractional_mjd = 0; // FIX
	      strcpy(calObsParams[i].utc,"UNKNOWN"); // FIX
	      strcpy(calObsParams[i].ut_date,"UNKNOWN"); // FIX
	      strcpy(calObsParams[i].local_time,"UNKNOWN"); // FIX  
	      strcpy(calObsParams[i].raStr,"UNKNOWN");
	      strcpy(calObsParams[i].decStr,"UNKNOWN");
	      calObsParams[i].raDeg = 0;
	      calObsParams[i].decDeg = 0;
	      calObsParams[i].raOffset = 0;
	      calObsParams[i].decOffset = 0;
	      calObsParams[i].gl = 0;
	      calObsParams[i].gb = 0;
	      calObsParams[i].az = 0;
	      calObsParams[i].ze = 0;
	      calObsParams[i].el = 0;
	      calObsParams[i].az_drive_rate = 0;
	      calObsParams[i].ze_drive_rate = 0;
	      calObsParams[i].hourAngle = 0;
	      calObsParams[i].paraAngle = 0;
	      calObsParams[i].windDir = 0;
	      calObsParams[i].windSpd = 0;      


	      
	      printf("Loading subintegration %d/%d\n",i,ndump_cal-1);
	      fits_read_col(fptr,TSHORT,colnum,i+1,1,nchan_cal*npol_cal*nbin_cal,&n_sval,shortVals,&initflag,&status);
	      // Select bin and scale by DAT_OFFS and DAT_SCL
	      for (j=0;j<npol_cal;j++)
		{
		  for (k=0;k<nchan_cal;k++)
		    {
		      floatCalValsOn[i*nchan_cal*npol_cal + j*nchan_cal + k] = 0;
		      floatCalValsOff[i*nchan_cal*npol_cal + j*nchan_cal + k] = 0;
		      nOn = nOff = 0;
		      for (b=0;b<nbin_cal;b++)
			{
			  if (b >= calOn1 && b <= calOn2)
			    {			    
			      floatCalValsOn[i*nchan_cal*npol_cal + j*nchan_cal + k] += (shortVals[k*nbin_cal+j*nchan_cal*nbin_cal + b])*datScl[j*nchan_cal+k]+datOffs[j*nchan_cal+k]; // WHAT ABOUT zeroOff?? -- CHECK
			      nOn++;
			    }
			  if (b>= calOff1 && b <= calOff2)
			    {
			      floatCalValsOff[i*nchan_cal*npol_cal + j*nchan_cal + k] += (shortVals[k*nbin_cal+j*nchan_cal*nbin_cal + b])*datScl[j*nchan_cal+k]+datOffs[j*nchan_cal+k]; // WHAT ABOUT zeroOff?? -- CHECK
			      nOff++;
			    }
			}
		      //		      printf("%d %d nOn = %d nOff = %d %d %d %d %d %g %g\n",j,k,nOn,nOff,calOn1,calOn2,calOff1,calOff2, floatCalValsOn[i*nchan_cal*npol_cal + j*nchan_cal + k],floatCalValsOff[i*nchan_cal*npol_cal + j*nchan_cal + k]);
		      floatCalValsOn[i*nchan_cal*npol_cal + j*nchan_cal + k] /= (double)nOn;
		      floatCalValsOff[i*nchan_cal*npol_cal + j*nchan_cal + k] /= (double)nOff;
		    }
		}
	    }
	
	  printf("Closing cal file\n");
	  fits_close_file(fptr,&status);
	  free(shortVals); free(datOffs); free(datScl);
	}
      // Now process the astronomy file
      
      printf("Opening file >%s<\n",fname);
      fits_open_file(&fptr,fname,READONLY,&status);
      fits_report_error(stderr,status);
      
      printf("Reading PSRFITS file\n");
      fits_movnam_hdu(fptr,BINARY_TBL,"SUBINT",1,&status);
      fits_get_num_rows(fptr,&long_ndump,&status);
      ndump = (int)long_ndump;

      doubleVals = (double *)malloc(sizeof(double)*ndump);
     
      fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
      fits_read_tdim(fptr,colnum,maxdim,&naxis,naxes,&status);
      fits_report_error(stderr,status);

      nbin  = naxes[0];
      nchan = naxes[1];
      npol  = naxes[2];
      printf("nbin = %d, nchan = %d, npol = %d, status = %d\n",nbin,nchan,npol,status);
      
      infloatVals = (float *)malloc(sizeof(float)*nchan*npol); // DON'T NEED THIS ONE .. BUT FREE AT THE END .. FIX ME
      floatVals   = (float *)malloc(sizeof(float)*nchan*npol*ndump*nbin);
      freqVals    = (float *)malloc(sizeof(float)*nchan);
      dfreqVals    = (double *)malloc(sizeof(double)*nchan);
      shortVals   = (short int *)malloc(sizeof(short int)*nchan*npol*nbin);
      baseline1   = (double *)malloc(sizeof(double)*nchan);
      baseline2   = (double *)malloc(sizeof(double)*nchan);
      baseline3   = (double *)malloc(sizeof(double)*nchan);
      baseline4   = (double *)malloc(sizeof(double)*nchan);
      
      
      fits_get_colnum(fptr,CASEINSEN,"DAT_OFFS",&colnum_datoffs,&status);
      fits_get_colnum(fptr,CASEINSEN,"DAT_SCL",&colnum_datscl,&status);
      datOffs = (float *)malloc(sizeof(float)*nchan*npol);
      datScl = (float *)malloc(sizeof(float)*nchan*npol);	
      
      fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&colnum,&status);
      fits_read_col(fptr,TFLOAT,colnum,1,1,nchan,&n_fval,freqVals,&initflag,&status);

      fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);

      // Do we want a baseline
      if (pulseOff1 >= 0 && pulseOff2 >= 0)
	{
	  int kk;
	  for (k=0;k<nchan;k++)
	    baseline1[k]=baseline2[k]=baseline3[k]=baseline4[k]=0;

	  for (i=0;i<ndump;i++)
	    {
	      fits_read_col(fptr,TFLOAT,colnum_datoffs,i+1,1,nchan*npol,&n_fval,datOffs,&initflag,&status);
	      fits_read_col(fptr,TFLOAT,colnum_datscl,i+1,1,nchan*npol,&n_fval,datScl,&initflag,&status);
	      
	      printf("Loading subintegration %d/%d\n",i,ndump-1);
	      fits_read_col(fptr,TSHORT,colnum,i+1,1,nchan*npol*nbin,&n_sval,shortVals,&initflag,&status);
	      // Select bin and scale by DAT_OFFS and DAT_SCL
	      for (k=0;k<nchan;k++)
		{
		  for (kk=pulseOff1;kk<pulseOff2;kk++)
		    {
		      j=0; baseline1[k] += (shortVals[k*nbin+j*nchan*nbin + kk])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		      j=1; baseline2[k] += (shortVals[k*nbin+j*nchan*nbin + kk])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		      j=2; baseline3[k] += (shortVals[k*nbin+j*nchan*nbin + kk])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		      j=3; baseline4[k] += (shortVals[k*nbin+j*nchan*nbin + kk])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		    }
		}
	    }

	  for (k=0;k<nchan;k++)
	    {
	      baseline1[k] /= (double)(ndump*(pulseOff2-pulseOff1));
	      baseline2[k] /= (double)(ndump*(pulseOff2-pulseOff1));
	      baseline3[k] /= (double)(ndump*(pulseOff2-pulseOff1));
	      baseline4[k] /= (double)(ndump*(pulseOff2-pulseOff1));
	      printf("baseline %d %g %g %g %g\n",k,baseline1[k],baseline2[k],baseline3[k],baseline4[k]);
	    }
	}
      else
	{
	  for (k=0;k<nchan;k++)
	    baseline1[k]=baseline2[k]=baseline3[k]=baseline4[k]=0;
	}

      for (i=0;i<ndump;i++)
	{
	  fits_read_col(fptr,TFLOAT,colnum_datoffs,i+1,1,nchan*npol,&n_fval,datOffs,&initflag,&status);
	  fits_read_col(fptr,TFLOAT,colnum_datscl,i+1,1,nchan*npol,&n_fval,datScl,&initflag,&status);

	  printf("Loading subintegration %d/%d\n",i,ndump-1);
	  fits_read_col(fptr,TSHORT,colnum,i+1,1,nchan*npol*nbin,&n_sval,shortVals,&initflag,&status);
	  // Select bin and scale by DAT_OFFS and DAT_SCL
	  for (j=0;j<npol;j++)
	    {
	      for (k=0;k<nchan;k++)
		{
		  floatVals[i*nchan*npol + j*nchan + k] = (shortVals[k*nbin+j*nchan*nbin + binVal])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		  if (j==0)      floatVals[i*nchan*npol + j*nchan + k] -= baseline1[k];
		  else if (j==1) floatVals[i*nchan*npol + j*nchan + k] -= baseline2[k];
		  else if (j==2) floatVals[i*nchan*npol + j*nchan + k] -= baseline3[k];
		  else if (j==3) floatVals[i*nchan*npol + j*nchan + k] -= baseline4[k];
		}
	    }
	}
      tdump = 1; // SET PROPERLY
      mjd = 56000; // SET PROPERLY 
      strcpy(srcName,"UNKNOWN"); // SET PROPERLY
      
      free(shortVals);
      printf("Closing file\n");
      fits_close_file(fptr,&status);
    }


  // Create the SDHDF file
  obsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
  printf("Here with ndump = %d\n",ndump);
  for (i=0;i<ndump;i++)
    {
      obsParams[i].timeElapsed = i*tdump; 
      strcpy(obsParams[i].timedb,"UNKNOWN"); // FIX
      obsParams[i].mjd = mjd;
      obsParams[i].fractional_mjd = 0.0; // FIX

      strcpy(obsParams[i].utc,"UNKNOWN"); // FIX
      strcpy(obsParams[i].ut_date,"UNKNOWN"); // FIX
      strcpy(obsParams[i].local_time,"UNKNOWN"); // FIX 
      strcpy(obsParams[i].raStr,"UNKNOWN");
      strcpy(obsParams[i].decStr,"UNKNOWN");
      obsParams[i].raDeg = 0;
      obsParams[i].decDeg = 0;
      obsParams[i].raOffset = 0;
      obsParams[i].decOffset = 0;
      obsParams[i].gl = 0;
      obsParams[i].gb = 0;
      obsParams[i].az = 0;
      obsParams[i].ze = 0;
      obsParams[i].el = 0;
      obsParams[i].az_drive_rate = 0;
      obsParams[i].ze_drive_rate = 0;
      obsParams[i].hourAngle = 0;
      obsParams[i].paraAngle = 0;
      obsParams[i].windDir = 0;
      obsParams[i].windSpd = 0;      
    }

  // Open the output HDF5 file
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }
  sdhdf_initialiseFile(outFile);

  // Set up the primary header information
      
  // Set up the beam information
  strcpy(beamHeader->source,srcName);
  strcpy(beamHeader->label,"beam_0");
  // Set up the band information
  strcpy(bandHeader->label,"band0");
  bandHeader->fc = (freqVals[0]+freqVals[nchan-1])/2.0; 
  bandHeader->f0 = freqVals[0];
  bandHeader->f1 = freqVals[nchan-1]; 
  bandHeader->nchan = nchan;
  bandHeader->npol = npol;
  if (npol==1)
    strcpy(bandHeader->pol_type,"AA+BB");
  else   
    strcpy(bandHeader->pol_type,"AABBCRCI");
  bandHeader->dtime = tdump;
  bandHeader->ndump = ndump;
  
  printf("Writing output file\n");
  sdhdf_openFile(outname,outFile,3);
  sdhdf_writePrimaryHeader(outFile,primaryHeader);
  sdhdf_writeBeamHeader(outFile,beamHeader,1,(char *)"4.0");
  sdhdf_writeBandHeader(outFile,bandHeader,beamHeader[0].label,1,1,(char *)"4.0");
  sdhdf_writeObsParams(outFile,bandHeader[0].label,beamHeader[0].label,0,obsParams,ndump,1,(char *)"4.0");
  // FIX ME: Only sending 1 frequency channel through
  for (i=0;i<nchan;i++)
    dfreqVals[i] = (double)freqVals[i];
  sdhdf_writeSpectrumData(outFile,beamHeader->label,bandHeader->label,0,0,floatVals,dfreqVals,1,nchan,1,npol,ndump,1,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);

  if (cal==1)
    {
      printf("Pos 4\n");
      sdhdf_writeBandHeader(outFile,calBandHeader,beamHeader[0].label,1,2,(char *)"4.0");
      sdhdf_writeObsParams(outFile,bandHeader[0].label,beamHeader[0].label,0,obsParams,ndump_cal,2,(char *)"4.0");
      // GEORGE HERE:
      // ... need to setup and then write the cal metadata
      // 

      // FIX ME: Only sending 1 frequency channel through
      for (i=0;i<nchan_cal;i++)
	dfreqCalVals[i] = (double)freqCalVals[i];
      sdhdf_writeSpectrumData(outFile,beamHeader->label,bandHeader->label,0,0,floatCalValsOn,dfreqCalVals,1,nchan_cal,1,npol_cal,ndump_cal,2,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
      sdhdf_writeSpectrumData(outFile,beamHeader->label,bandHeader->label,0,0,floatCalValsOff,dfreqCalVals,1,nchan_cal,1,npol_cal,ndump_cal,3,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	   }
  sdhdf_writeSoftwareVersions(outFile,softwareVersions);
  sdhdf_writeHistory(outFile,history,1);
  sdhdf_closeFile(outFile);
  free(outFile);
  if (cal==1)
    {
      free(floatCalValsOn);
      free(floatCalValsOff);
      free(freqCalVals);
      free(dfreqCalVals);
      free(calBandHeader);
      free(calObsParams);	    
    }
if (type==2)
  {
    free(baseline1); free(baseline2); free(baseline3); free(baseline4);
  }
  free(primaryHeader);
  free(beamHeader);
  free(bandHeader);
  free(obsParams);
  free(softwareVersions);
  free(history);
  free(floatVals);
  free(infloatVals);
  free(freqVals);
  free(dfreqVals);
  free(doubleVals);
} 


//
// Read the header information from the filterbank file
//

void readFilHeader(char *fname,FILE *inputfile,header *head)
{
  char string[80], message[80];
  int itmp,nbytes,totalbytes,expecting_rawdatafile=0,expecting_source_name=0; 
  int expecting_frequency_table=0,channel_index;
  /* added frequency table for use with non-contiguous data */
  double frequency_table[4096]; /* note limited number of channels */
  long int npuls; /* added for binary pulse profile format */
  
  head->fil_isign=0;


  /* try to read in the first line of the header */
  get_string(inputfile,&nbytes,string);
  if (!strings_equal(string,"HEADER_START")) {
	/* the data file is not in standard format, rewind and return */
	rewind(inputfile);
	exit(1);
	return;
  }
  /* store total number of bytes read so far */
  totalbytes=nbytes;

  /* loop over and read remaining header lines until HEADER_END reached */
  while (1) {
    get_string(inputfile,&nbytes,string);
    printf("STRING %s\n",string);
    if (strings_equal(string,"HEADER_END")) break;
    totalbytes+=nbytes;
    if (strings_equal(string,"rawdatafile")) {
      expecting_rawdatafile=1;
    } else if (strings_equal(string,"source_name")) {
      expecting_source_name=1;
    } else if (strings_equal(string,"FREQUENCY_START")) {
      expecting_frequency_table=1;
      channel_index=0;
    } else if (strings_equal(string,"FREQUENCY_END")) {
      expecting_frequency_table=0;
    } else if (strings_equal(string,"az_start")) {
      fread(&(head->fil_az_start),sizeof(head->fil_az_start),1,inputfile);
      totalbytes+=sizeof(head->fil_az_start);
    } else if (strings_equal(string,"za_start")) {
      fread(&(head->fil_za_start),sizeof(head->fil_za_start),1,inputfile);
      totalbytes+=sizeof(head->fil_za_start);
    } else if (strings_equal(string,"src_raj")) {
      fread(&(head->fil_src_raj),sizeof(head->fil_src_raj),1,inputfile);
      totalbytes+=sizeof(head->fil_src_raj);
    } else if (strings_equal(string,"src_dej")) {
      fread(&(head->fil_src_decj),sizeof(head->fil_src_decj),1,inputfile);
      totalbytes+=sizeof(head->fil_src_decj);
    } else if (strings_equal(string,"tstart")) {
      fread(&(head->fil_tstart),sizeof(head->fil_tstart),1,inputfile);
      totalbytes+=sizeof(head->fil_tstart);
    } else if (strings_equal(string,"tsamp")) {
      fread(&(head->fil_tsamp),sizeof(head->fil_tsamp),1,inputfile);
      totalbytes+=sizeof(head->fil_tsamp);
    } else if (strings_equal(string,"period")) {
      fread(&(head->fil_period),sizeof(head->fil_period),1,inputfile);
      totalbytes+=sizeof(head->fil_period);
    } else if (strings_equal(string,"fch1")) {
      fread(&(head->fil_fch1),sizeof(head->fil_fch1),1,inputfile);
      totalbytes+=sizeof(head->fil_fch1);
    } else if (strings_equal(string,"fchannel")) {
      fread(&frequency_table[channel_index++],sizeof(double),1,inputfile);
      totalbytes+=sizeof(double);
      head->fil_fch1=head->fil_foff=0.0; /* set to 0.0 to signify that a table is in use */
    } else if (strings_equal(string,"foff")) {
      fread(&(head->fil_foff),sizeof(head->fil_foff),1,inputfile);
      totalbytes+=sizeof(head->fil_foff);
    } else if (strings_equal(string,"nchans")) {
      fread(&(head->fil_nchans),sizeof(head->fil_nchans),1,inputfile);
      totalbytes+=sizeof(head->fil_nchans);
    } else if (strings_equal(string,"telescope_id")) {
      fread(&(head->fil_telescope_id),sizeof(head->fil_telescope_id),1,inputfile);
      totalbytes+=sizeof(head->fil_telescope_id);
    } else if (strings_equal(string,"machine_id")) {
      fread(&(head->fil_machine_id),sizeof(head->fil_machine_id),1,inputfile);
      totalbytes+=sizeof(head->fil_machine_id);
    } else if (strings_equal(string,"data_type")) {
      fread(&(head->fil_data_type),sizeof(head->fil_data_type),1,inputfile);
      totalbytes+=sizeof(head->fil_data_type);
    } else if (strings_equal(string,"ibeam")) {
    fread(&(head->fil_ibeam),sizeof(head->fil_ibeam),1,inputfile);
      totalbytes+=sizeof(head->fil_ibeam);
    } else if (strings_equal(string,"nbeams")) {
    fread(&(head->fil_nbeams),sizeof(head->fil_nbeams),1,inputfile);
      totalbytes+=sizeof(head->fil_nbeams);
    } else if (strings_equal(string,"nbits")) {
    fread(&(head->fil_nbits),sizeof(head->fil_nbits),1,inputfile);
      totalbytes+=sizeof(head->fil_nbits);
    } else if (strings_equal(string,"barycentric")) {
    fread(&(head->fil_barycentric),sizeof(head->fil_barycentric),1,inputfile);
      totalbytes+=sizeof(head->fil_barycentric);
    } else if (strings_equal(string,"pulsarcentric")) {
    fread(&(head->fil_pulsarcentric),sizeof(head->fil_pulsarcentric),1,inputfile);
      totalbytes+=sizeof(head->fil_pulsarcentric);
    } else if (strings_equal(string,"nbins")) {
    fread(&(head->fil_nbins),sizeof(head->fil_nbins),1,inputfile);
      totalbytes+=sizeof(head->fil_nbins);
    } else if (strings_equal(string,"nsamples")) {
      printf("IN HERE\n");
      fread(&(head->fil_nsamples),sizeof(head->fil_nsamples),1,inputfile);
      totalbytes+=sizeof(itmp);
    } else if (strings_equal(string,"nifs")) {
    fread(&(head->fil_nifs),sizeof(head->fil_nifs),1,inputfile);
      totalbytes+=sizeof(head->fil_nifs);
    } else if (strings_equal(string,"npuls")) {
    fread(&(head->fil_npuls),sizeof(head->fil_npuls),1,inputfile);
      totalbytes+=sizeof(head->fil_npuls);
    } else if (strings_equal(string,"refdm")) {
    fread(&(head->fil_refdm),sizeof(head->fil_refdm),1,inputfile);
      totalbytes+=sizeof(head->fil_refdm);
    } else if (strings_equal(string,"signed")) {
    fread(&(head->fil_isign),sizeof(head->fil_isign),1,inputfile);
      totalbytes+=sizeof(head->fil_isign);
    } else if (expecting_rawdatafile) {
      strcpy(head->fil_rawdatafile,string);
      expecting_rawdatafile=0;
    } else if (expecting_source_name) {
      strcpy(head->fil_source_name,string);
      expecting_source_name=0;
    } else {
      sprintf(message,"read_header - unknown parameter: %s\n",string);
      fprintf(stderr,"ERROR: %s\n",message);
      exit(1);
    } 
    if (totalbytes != ftell(inputfile)){
	    fprintf(stderr,"ERROR: Header bytes does not equal file position\n");
	    fprintf(stderr,"String was: '%s'\n",string);
	    fprintf(stderr,"       header: %d file: %d\n",totalbytes,(int)ftell(inputfile));
	    exit(1);
    }

    if (head->fil_isign < 0 && OSIGN > 0){
	    fprintf(stderr,"WARNING! You are reading unsigned numbers with a signed version of sigproc\n");
    }
    if (head->fil_isign > 0 && OSIGN < 0){
	    fprintf(stderr,"WARNING! You are reading signed numbers with a unsigned version of sigproc\n");
    }


  } 

  /* add on last header string */
  totalbytes+=nbytes;

  if (totalbytes != ftell(inputfile)){
	  fprintf(stderr,"ERROR: Header bytes does not equal file position\n");
	  fprintf(stderr,"       header: %d file: %d\n",totalbytes,(int)ftell(inputfile));
	  exit(1);
  }
  printf("total bytes = %d file position = %d\n",totalbytes,(int)ftell(inputfile));

}

// This routine is from read_header.c in sigproc
void get_string(FILE *fin,int *nbytes,char string[])
{
  int nchar;
  strcpy(string,"ERROR");
  fread(&nchar, sizeof(int),1,fin);
  printf("Have read in %d chars, sizeof(int) is %d\n",nchar,(int)sizeof(int));
  *nbytes = sizeof(int);

  if (feof(fin)) exit(0);
  if (nchar > 80 || nchar<1) return;
  fread(string,nchar,1,fin);
  string[nchar]='\0';
  *nbytes+=nchar;
}

// From sigproc
int strings_equal (char *string1, char *string2) /* includefile */
{
  if (!strcmp(string1,string2)) {
    return 1;
  } else {
    return 0;
  }
}

void bytesToFloats(int samplesperbyte,int n,unsigned char *cVals,float *out)
{
  int i,j;
  int pos=0;
  int index = 0;

  for (i = 0; i < n/samplesperbyte; i++)
    {
      void (*pointerFloat)(int, float *, int*);
      switch (samplesperbyte)
	{
	case 1:
	  pointerFloat = eightBitsFloat;
	  break;

	case 2:
	  pointerFloat = fourBitsFloat;
	  break;

	case 4:
	  //	  pointerFloat = twoBitsFloat;
	  pointerFloat = twoBitsFloat_unsigned;
	  break;

	case 8:
	  pointerFloat = oneBitFloat;
	  break;
	}
      pointerFloat(cVals[i], out, &index);
    }
}


void eightBitsFloat(int eight_bit_number, float *results, int *index)
{
  // 135 subtracted from the original number to make the range from -8 to something
  // 0.5 is then added to remove bias

  //results[(*index)++] = eight_bit_number - 135.5;
  results[(*index)++] = eight_bit_number; // -31.5; // -31.5;
}

void fourBitsFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant 4 bits with 0000 1111 will give the (signed) 4-bit number
  // shifting right 4 bits will produce the other (first) 4-bit numbers
  // 0.5 is added to each number to compensate for the bias, making the range -7.5 -> +7.5

  float tempResults[2];
  int i;
  for (i = 0; i < 2; i++) {
    int andedNumber = eight_bit_number & 15;
    tempResults[i] = andedNumber; // - 7.5;
    eight_bit_number = eight_bit_number >> 4;
  }

  for (i = 1; i >= 0; i--)
    results[(*index)++] = tempResults[i];
}

void twoBitsFloat_unsigned(int eight_bit_number, float *results, int *index)
{
  unsigned char uctmp = (unsigned char)eight_bit_number;
  int i;

  results[(*index)++] = ((uctmp >> 0x06) & 0x03);
  results[(*index)++] = ((uctmp >> 0x04) & 0x03);
  results[(*index)++] = ((uctmp >> 0x02) & 0x03);
  results[(*index)++] = ((uctmp & 0x03));
}

// Note DSPSR is using this technique, but then doesn't scale by ZERO_OFFS
void twoBitsFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant 2 bits with 0000 0011 will give the (signed 2-bit number
  // shifting right 2 bits will produce the next (previous) 2-bit number
  // the numbers are adjusted to compensate for bias and to give a rms of ~0.9
  //
  // 1 -> 2.0
  // 0 -> +0.5
  // -1 -> -0.5
  // -2 -> -2.0


  // NOTE: ACTUAL -2.5 below !!
  float tempResults[4];
  int i;

  for (i = 0; i < 4; i++) {
    int andedNumber = eight_bit_number & 3;

    switch (andedNumber) {
    case 0:
      //      tempResults[i] = 0; //-2.5;
      tempResults[i] = -2.5;
      break;
    case 1:
      //      tempResults[i] = 1; // -0.5;
      tempResults[i] = -0.5;      
      break;
    case 2:
      tempResults[i] = 0.5;
      //      tempResults[i] = 2; //0.5;
      break;
    case 3:
      tempResults[i] = 2.5;
      //      tempResults[i] = 3; // 2.5;
      break;
    }
    eight_bit_number = eight_bit_number >> 2;
  }

  for (i = 3; i >= 0; i--)
    {
      //      printf("Have %f\n",tempResults[i]);
      results[(*index)++] = tempResults[i];
    }
}

void oneBitFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant bit with 0000 0001 will give the (signed 1-bit number
  // shifting right 1 bit will produce the next (previous) 1-bit number
  //
  // 0.5 is added to each number to compensate for bias
  //

  float tempResults[8];
  int i;

  for (i = 0; i < 8; i++) {
    int andedNumber = eight_bit_number & 1;
    // tempResults[i] = andedNumber ? 0.5 : -0.5;
    tempResults[i] = andedNumber ? 0 : 1;
    eight_bit_number = eight_bit_number >> 1;
  }

  for (i = 7; i >= 0; i--)
    {
      results[(*index)++] = tempResults[i];
      //            printf("This pt: %d %g\n",*index,tempResults[i]);
    }
}

