//  Copyright (C) 2019, 2020, 2021, 2022, 2023 George Hobbs

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

void help()
{
  printf("\nsdhdf_convertPulsarSearchMode  %s\n",VNUM);
	printf("INSPECTA version:              %s\n",SOFTWARE_VER);
  printf("Author:                        George Hobbs\n");
  printf("Software to convert pulsar search-mode data to SDHDF format\n");

  printf("\nCommand line arguments:\n\n");
	printf("-h                This help\n");
	printf("-f <filename>            Filename to convert\n");
	printf("-o <filename>         Output filename\n");

	printf("\nExample:\n\n");
  printf("sdhdf_convertPulsarSearchMode -f file.sf\n\n");

  exit(1);
}

int main(int argc,char *argv[])
{
  int i,j,k,b;

  int fitsType=1; // 1 = SDFITS, 2 = PSRFITS
  
  sdhdf_fileStruct *outFile;
  sdhdf_softwareVersionsStruct *softwareVersions;
  sdhdf_historyStruct *history;
  sdhdf_primaryHeaderStruct *primaryHeader;
  sdhdf_beamHeaderStruct *beamHeader;
  sdhdf_bandHeaderStruct *bandHeader;
  sdhdf_obsParamsStruct  *obsParams;
  sdhdf_obsParamsStruct  *calObsParams;
  sdhdf_bandHeaderStruct *calBandHeader;

  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;
    
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
  unsigned char *storeVals;
  float  *freqVals;
  float *datOffs;
  float *datScl;

  unsigned char nval='0';
  float n_fval=0;
  int    colnum;
  int    colnum_datoffs;
  int    colnum_datscl;

  int    initflag=0;
  double mean=0;

  char fname[1024];
  char outname[1024]="convert.hdf";



  int nBeam;
  int nsblk;
  int samplesperbyte;
  
  nBeam = 1; // FIX ME
  
  primaryHeader    = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct));
  beamHeader       = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nBeam);
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
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outname,argv[++i]);
    }


  // Assuming PSRFITS search mode
  
  printf("Opening file >%s<\n",fname);
  fits_open_file(&fptr,fname,READONLY,&status);
  fits_report_error(stderr,status);
  
  printf("Reading PSRFITS file\n");
  fits_movnam_hdu(fptr,BINARY_TBL,"SUBINT",1,&status);
  fits_get_num_rows(fptr,&long_ndump,&status);
  ndump = (int)long_ndump; // Number of subints in the original file
    
  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
  fits_read_tdim(fptr,colnum,maxdim,&naxis,naxes,&status);
  fits_report_error(stderr,status);
  
  nchan = naxes[0];
  npol  = naxes[1];
  nsblk  = naxes[2];
  printf("nsblk = %d, nchan = %d, npol = %d, status = %d\n",nsblk,nchan,npol,status);


    samplesperbyte = 8; // FIX ME ****
  if (!(storeVals   = (unsigned char *)malloc(sizeof(unsigned char)*nchan*npol*ndump*nsblk/samplesperbyte)))
    {
      printf("ERROR: Unable to allocate sufficient memory\n");
      exit(1);
    }
  printf("data column num = %d\n",colnum);
  fits_get_colnum(fptr,CASEINSEN,"DAT_OFFS",&colnum_datoffs,&status);
  printf("Status = %d, colnum_datoffs = %d, nchan = %d, npol = %d\n",status,colnum_datoffs,nchan,npol);
  fits_get_colnum(fptr,CASEINSEN,"DAT_SCL",&colnum_datscl,&status);
  printf("Status = %d, colnum_datscl = %d, nchan = %d, npol = %d\n",status,colnum_datscl,nchan,npol);
  if (!(datOffs = (float *)malloc(sizeof(float)*nchan*npol)))
    {
      printf("ERROR: Unable to allocate sufficient memory\n");
      exit(1);
    }
  if (!(datScl = (float *)malloc(sizeof(float)*nchan*npol)))
    {
      printf("ERROR: Unable to allocate sufficient memory\n");
      exit(1);      
    }
  freqVals = (float *)malloc(sizeof(float)*nchan);
  printf("Got here\n");
  printf("colnum_datoffs = %d\n",colnum_datoffs);
  fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&colnum,&status);
  printf("colnum for DAT_FREQ = %d\n",colnum);
  fits_read_col(fptr,TFLOAT,colnum,1,1,nchan,&n_fval,freqVals,&initflag,&status);
  printf("Status = %d\n",status);
  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);  

  
  for (i=0;i<ndump;i++)
    {
      fits_read_col(fptr,TFLOAT,colnum_datoffs,i+1,1,nchan*npol,&n_fval,datOffs,&initflag,&status);
      fits_read_col(fptr,TFLOAT,colnum_datscl,i+1,1,nchan*npol,&n_fval,datScl,&initflag,&status);
      
      printf("Loading subintegration %d/%d\n",i,ndump-1);
      fits_read_col_byt(fptr,colnum,i+1,1,nchan*npol*nsblk/samplesperbyte,nval,storeVals+i*nchan*npol*nsblk/samplesperbyte,&initflag,&status);
    }

  printf("Closing file, status = %d\n",status);
  printf("ndump = %d\n",ndump);
  fits_close_file(fptr,&status);
  obsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
  for (i=0;i<ndump;i++)
    {
      obsParams[i].timeElapsed = i; // FIX
      strcpy(obsParams[i].timedb,"UNKNOWN"); // FIX
      obsParams[i].mjd = 56000; // FIX
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
  printf("nbeam = %d\n",nBeam);
  primaryHeader[0].nbeam = nBeam;
  
  // Set up the beam information

  // Set up the band information
  strcpy(bandHeader->label,"band0");
  bandHeader->fc = (freqVals[0]+freqVals[nchan-1])/2.0; 
  bandHeader->f0 = freqVals[0];
  bandHeader->f1 = freqVals[nchan-1]; 
  bandHeader->nchan = nchan;
  bandHeader->npol = npol;
  //  strcpy(bandHeader->pol_type,"AABBCRCI");
  strcpy(bandHeader->pol_type,"AABB");
  bandHeader->dtime = 100; // FIX
  bandHeader->ndump = 1;

  printf("Writing output file %s\n",outname);
  sdhdf_openFile(outname,outFile,3);

  sdhdf_writeBeamHeader(outFile,beamHeader,nBeam,(char *)"4.0");
  sdhdf_writePrimaryHeader(outFile,primaryHeader);
  for (b=0;b<nBeam;b++)
    {
      printf("Write band %s %s\n",beamHeader[b].label,bandHeader[b].label);
      sdhdf_writeBandHeader(outFile,bandHeader,beamHeader[b].label,1,1,(char *)"4.0");

      printf("Write obs\n");
      sdhdf_writeObsParams(outFile,bandHeader[0].label,beamHeader[b].label,0,obsParams,ndump,1,(char *)"4.0");
      // FIX ME: Only sending 1 frequency channel through
      printf("Writing data\n");
      sdhdf_writeQuantisedSpectrumData(outFile,beamHeader[b].label,bandHeader->label,b,0,storeVals,freqVals,1,nchan,nsblk*ndump/samplesperbyte,npol,1,1,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
      printf("Finished writing data\n");
    }
  sdhdf_writeSoftwareVersions(outFile,softwareVersions);
  sdhdf_writeHistory(outFile,history,1);
  
  sdhdf_closeFile(outFile);
  free(outFile);
  free(storeVals);

  if (nBeam > 0)
    free(beamHeader);
  free(primaryHeader);
  free(bandHeader);
  free(obsParams);
  free(softwareVersions);
  free(history);
  free(freqVals);
} 
