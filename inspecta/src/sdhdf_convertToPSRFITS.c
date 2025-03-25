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
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "fitsio.h"

#define VERSION "v0.5"

void createPrimaryHeader(fitsfile *cal,fitsfile *psr,  sdhdf_fileStruct *inFile);
void createSubintHeader(fitsfile *cal,fitsfile *psr,  sdhdf_fileStruct *inFile);
void createDataForCal(fitsfile *cal,sdhdf_fileStruct *inFile,int fixFreq);
void createDataForAstro(fitsfile *astro,sdhdf_fileStruct *inFile,int fixFreq);

void help()
{
  printf("sdhdf_convertToPSRFITS %s (INSPECTA %s)\n",VERSION,SOFTWARE_VER);
  printf("Authors: G. Hobbs\n");

  exit(1);
}



int main(int argc,char *argv[])
{
  int i,j,k;
  int fixFreq=0;
  sdhdf_fileStruct *inFile;
  char outName[1024];
  char outFname[1024];
  char fname[1024];
  fitsfile *fp_psr,*fp_cal;
  int status=0;
  
  if (argc==1) help();
    // Allocate memory for the input ifle
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outName,argv[++i]);
      else if (strcasecmp(argv[i],"-fixfreq")==0)
	fixFreq=1;
      else if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }

  // Open the input SDHDF file
  sdhdf_initialiseFile(inFile);
  sdhdf_openFile(fname,inFile,1);
  sdhdf_loadMetaData(inFile);
  printf("%s %s %s %s %s %s\n",inFile->fname,inFile->primary[0].pid,inFile->beamHeader[0].source,inFile->beam[0].bandData[0].astro_obsHeader[0].raStr,inFile->beam[0].bandData[0].astro_obsHeader[0].decStr,inFile->primary[0].utc0);
  
  // Create the PSRFITS output file

  sprintf(outFname,"!%s.cf(psrheader.fits)",outName);
  fits_create_file(&fp_cal,outFname,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  sprintf(outFname,"!%s.rf(psrheader.fits)",outName);
  fits_create_file(&fp_psr,outFname,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
      
  createPrimaryHeader(fp_cal,fp_psr,inFile);
  createSubintHeader(fp_cal,fp_psr,inFile);
  createDataForCal(fp_cal,inFile,fixFreq);
  createDataForAstro(fp_psr,inFile,fixFreq);
  sdhdf_closeFile(inFile);
  
  fits_close_file(fp_cal,&status);
  fits_close_file(fp_psr,&status);
  free(inFile);
  
  
}

void createDataForAstro(fitsfile *astro,sdhdf_fileStruct *inFile,int fixFreq)
{
  int i,j,k;
  int i0,j0,k0,l0,m0;
  int16_t *dataOut;
  float *dataOutF;
  int nchan;
  int npol=4;
  int nsub;
  int nbin=1; // Only 1 bin for continuum/spectral line data
  int status=0;
  int colnum_data;
  long naxes[4];
  int naxis=3;
  char tdim[16];
  double *raSub,*decSub;
  int colnum_datFreq, colnum_datWts,colnum_datScl,colnum_datOffs,colnum_raSub,colnum_decSub;
  int colnum_posAng,colnum_paraAng;
  int colnum_period; 
  int colnum_tsubint;
  char dataName[1024];
  hid_t dataset_id;
  float *datScl,*datOffs;
  float *data;
  float *floatArr;
  int ndump,nd;
  int nchanBand;
  long int n=0;
  int nband = inFile->beam[0].nBand;
  double sumFreq,chanbw,f0,f1;
  float cfreq,bw;
  float tfloat;
  
  double minShort,maxShort;
  double minVal,maxVal,tv;
  double period = 1; // FIX PERIOD
  minShort = INT16_MIN;
  maxShort = INT16_MAX;

  nchan=0;
  ndump = 0;
  for (j=0;j<nband;j++)
    {
      nchan += inFile->beam[0].bandHeader[j].nchan;
      ndump += inFile->beam[0].bandHeader[j].ndump;
    }
  nchanBand = inFile->beam[0].bandHeader[0].nchan; // ASSUME CONSTANT
  
  printf("astro nchan = %d, ndump = %d, nchanBand = %d\n",nchan,ndump,nchanBand);  
  dataOut = (int16_t *)malloc(sizeof(int16_t)*nchan*nbin*npol);
  dataOutF = (float *)malloc(sizeof(float)*nchan*nbin*npol*ndump);
  datScl = (float *)malloc(sizeof(float)*nchan*npol);
  datOffs = (float *)malloc(sizeof(float)*nchan*npol);
  fits_movnam_hdu(astro,BINARY_TBL,(char *)"SUBINT",0,&status);

  fits_get_colnum(astro, CASEINSEN, "RA_SUB", &colnum_raSub, &status);
  fits_get_colnum(astro, CASEINSEN, "DEC_SUB", &colnum_decSub, &status);
  fits_get_colnum(astro, CASEINSEN, "DATA", &colnum_data, &status);
  fits_get_colnum(astro, CASEINSEN, "TSUBINT", &colnum_tsubint, &status);
  fits_get_colnum(astro, CASEINSEN, "POS_ANG", &colnum_posAng, &status);
  fits_get_colnum(astro, CASEINSEN, "PAR_ANG", &colnum_paraAng, &status);
  fits_get_colnum(astro, CASEINSEN, "PERIOD", &colnum_period, &status);

  if (status) {fits_report_error(stderr, status); exit(1);}

  fits_modify_vector_len(astro,colnum_data,(long)((long)nchan*(long)npol*(long)nbin),&status);

  fits_get_colnum(astro, CASEINSEN, "DAT_FREQ", &colnum_datFreq, &status);
  fits_modify_vector_len(astro,colnum_datFreq,nchan,&status);
  fits_get_colnum(astro, CASEINSEN, "DAT_WTS", &colnum_datWts, &status);
  fits_modify_vector_len(astro,colnum_datWts,nchan,&status);
  fits_get_colnum(astro, CASEINSEN, "DAT_SCL", &colnum_datScl, &status);
  fits_modify_vector_len(astro,colnum_datScl,nchan*npol,&status);
  fits_get_colnum(astro, CASEINSEN, "DAT_OFFS", &colnum_datOffs, &status);
  fits_modify_vector_len(astro,colnum_datOffs,nchan*npol,&status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  
  naxes[0] = nbin;
  naxes[1] = nchan;
  naxes[2] = npol;
  
  sprintf(tdim,"TDIM%d",colnum_data);
  fits_delete_key(astro, tdim, &status);
  fits_write_tdim(astro, colnum_data, naxis, naxes, &status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  printf("Complete re-sizeing file\n");
  floatArr = (float *)malloc(sizeof(float)*nchan*npol);
  for (j0 = 0; j0 < nband; j0++)
    {
      sprintf(dataName,"beam_%d/%s/astronomy_data/data",0,inFile->beam[0].bandHeader[j0].label);
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      nd = inFile->beam[0].bandHeader[0].ndump;
      if (!(data = (float *)malloc(sizeof(float)*nchan*npol*nd*nbin)))
	{
	  printf("UNABLE to allocate enough memory in convertToPSRFITS\n");
	  exit(1);
	}

      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);  
      
      for (l0 = 0; l0 <  inFile->beam[0].bandHeader[j0].ndump; l0++)
	{
	  for (k0 = 0; k0 < inFile->beam[0].bandHeader[j0].nchan; k0++)
	    {
	      for (i0 = 0;i0 < npol;i0++)
		{
		  for (m0=0;m0<nbin;m0++)
		    {
		      dataOutF[l0*nchan*npol*nbin + i0*nchan*nbin + (j0*nchanBand+k0)*nbin + m0] =
			data[l0*nchanBand*npol*nbin + i0*nbin*nchanBand + k0*nbin + m0]; // Could be sped up
		    }
		}
	    }
	}
      //      printf("%g %d\n",data[i],dataOut[i]);
      status = H5Dclose(dataset_id);
      free(data);
    }
  printf("Writing astro data to PSRFITS\n");

  f0 = 1e6; f1 = -1e6;
  sumFreq = 0;
  raSub  = (double *)malloc(sizeof(double)*inFile->beam[0].bandHeader[0].ndump);
  decSub = (double *)malloc(sizeof(double)*inFile->beam[0].bandHeader[0].ndump);
  
  for (i = 0; i <  inFile->beam[0].bandHeader[0].ndump; i++) // Assume all bands have same number of dumps
    {
      raSub[i] = inFile->beam[0].bandData[0].astro_obsHeader[i].raDeg;
      decSub[i] = inFile->beam[0].bandData[0].astro_obsHeader[i].raDeg;
      printf("ra/dec = %g %g\n",raSub[i],decSub[i]);
      // Convert to SHORT and determine DAT_SCL and DAT_OFFS
      for (j0=0;j0<npol;j0++)
	{
	  for (i0=0;i0<nchan;i0++)
	    {
	      for (k0=0;k0<nbin;k0++)
		{
		  tv = dataOutF[i*nchan*nbin*npol + j0*nchan*nbin + i0*nbin + k0];
		  if (k0==0)
		    minVal = maxVal = tv;
		  else
		    {
		      if (minVal > tv) minVal = tv;
		      if (maxVal < tv) maxVal = tv;
		    }
		}
	      //	      printf("(%d,%d) = %g %g\n",j0,i0,minVal,maxVal);
	      datScl[j0*nchan+i0]  = 1; // (maxVal-minVal)/(maxShort-minShort); // WHAT SHOULD THIS BE WHEN ONLY 1 BIN
	      datOffs[j0*nchan+i0] = dataOutF[i*nchan*nbin*npol + j0*nchan*nbin + i0*nbin]; // maxVal-maxShort*datScl[j0*nchan+i0]; // NOTE BIN 0
	    }
	}
      for (j0=0;j0<npol;j0++)
	{
	  for (i0=0;i0<nchan;i0++)
	    {
	      for (k0=0;k0<nbin;k0++)
		{
		  dataOut[j0*nchan*nbin + i0*nbin + k0] = 1.0; //(int16_t)((dataOutF[i*nchan*nbin*npol + j0*nchan*nbin + i0*nbin+k0] - datOffs[j0*nchan+i0])/datScl[j0*nchan+i0]);
		  //		  printf("dataOut = %d %g %g %g\n", dataOut[j0*nchan*nbin + i0*nbin + k0],dataOutF[i*nchan*nbin*npol + j0*nchan*nbin + i0*nbin+k0],datOffs[j0*nchan+i0],datScl[j0*nchan+i0]);
		  //		  printf("From %g\n",dataOutF[i*nchan*nbin*npol + j0*nchan*nbin + i0*nbin+k0]);
		}
	    }
	}

      // dataOut = 32036 1.01508e+06 0 1
 //dataOut = -20391 503898 0 1
 //dataOut = -1840 522448 0 1
 //dataOut = -22997 501292 0 1
 //dataOut = -932 523356 0 1
 //dataOut = -25781 498507 0 1
 //dataOut = -5406 518882 0 1
 //dataOut = -785 523503 0 1
 //dataOut = -24491 499797 0 1
 //dataOut = -3861 520428 0 1
 //dataOut = -27973 496315 0 1
 //dataOut = -9010 515279 0 1

      
      printf("Writing new row: %d\n",i);
      //      if (i>0)
      //	fits_insert_rows(cal,i,1,&status);
      fits_write_col(astro,TSHORT,colnum_data,i+1,1,nchan*nbin*npol,dataOut,&status);  
      
      printf("Complete writing data: sttus = %d\n",status);
      if (status) {fits_report_error(stderr, status); exit(1);}
      fits_write_col(astro,TFLOAT,colnum_datScl,i+1,1,nchan*npol,datScl,&status);  
      fits_write_col(astro,TFLOAT,colnum_datOffs,i+1,1,nchan*npol,datOffs,&status);  
    
      for (j=0;j<nchan;j++) floatArr[j] = 1; // WTS
      fits_write_col(astro,TFLOAT,colnum_datWts,i+1,1,nchan,floatArr,&status);  
      printf("Status here = %d\n",status);
      
      // Frequency axis
      printf("Frequency axis\n");
      n=0;

      for (k=0;k<inFile->beam[0].nBand;k++)
	{
	  sdhdf_loadBandData(inFile,0,k,1);
	  
	  //	  printf("k = %d %d %d\n",k,n,inFile->beam[0].bandHeader[k].nchan);
	  for (j=0;j<inFile->beam[0].bandHeader[k].nchan;j++)
	    {
	      if (fixFreq==0)
		floatArr[n] = inFile->beam[0].bandData[k].astro_data.freq[j];  // FIX ME -- FREQ AXIS FOR DUMP 
	      else
		{
		  double cbw,fs,fe;




		  fs  = inFile->beam[0].bandData[k].astro_data.freq[0];  // FIX ME -- FREQ AXIS FOR DUMP
		  fe  = inFile->beam[0].bandData[k].astro_data.freq[inFile->beam[0].bandHeader[k].nchan-1];  // FIX ME -- FREQ AXIS FOR DUMP
		  cbw = (fe-fs)/inFile->beam[0].bandHeader[k].nchan;


		  floatArr[n] = (fs + 0.5*cbw) + j*cbw;
		}

	      sumFreq += floatArr[n];
	      if (f0 > floatArr[n]) f0 = floatArr[n];
	      if (f1 < floatArr[n]) f1 = floatArr[n];
	      n++;
	    }
	  if (k==0)
	    chanbw = floatArr[1] - floatArr[0]; // Assuming constant
	  sdhdf_releaseBandData(inFile,0,k,2);
	  sdhdf_releaseBandData(inFile,0,k,3); 		      

	}
      printf("nchan and n = %d %d\n",(int)nchan,(int)n);
      if (i==0)
	{
	  cfreq = sumFreq/n;
	  bw = nchan*chanbw;
	}
      fits_write_col(astro,TFLOAT,colnum_datFreq,i+1,1,nchan,floatArr,&status);  
      fits_write_col(astro,TDOUBLE,colnum_raSub,i+1,1,1,&raSub[i],&status);
      fits_write_col(astro,TDOUBLE,colnum_decSub,i+1,1,1,&raSub[i],&status);
      tfloat = inFile->beam[0].bandData[0].astro_obsHeader[i].paraAngle;
      fits_write_col(astro,TFLOAT,colnum_paraAng,i+1,1,1,&tfloat,&status);
      fits_write_col(astro,TDOUBLE,colnum_period,i+1,1,1,&period,&status);
      fits_write_col(astro,TDOUBLE,colnum_tsubint,i+1,1,1,&inFile->beam[0].bandHeader[0].dtime,&status);  
      printf("status at bottom = %d, cfreq = %g, bw = %g, chbw = %g\n",status,cfreq,bw,chanbw);
      
    }

  free(dataOut);
  free(dataOutF);
  free(datScl);
  free(datOffs);
  free(floatArr);
  free(raSub);
  free(decSub);
  fits_update_key(astro,TFLOAT,"CHAN_BW",&chanbw,NULL,&status);
  // Update the primary header
  fits_movabs_hdu(astro,1,NULL,&status);
  fits_update_key(astro,TFLOAT,"OBSFREQ",&cfreq,NULL,&status);
  fits_update_key(astro,TFLOAT,"OBSBW",&bw,NULL,&status);
  fits_update_key(astro,TINT,"OBSNCHAN",&nchan,NULL,&status);   
}


void createDataForCal(fitsfile *cal,sdhdf_fileStruct *inFile,int fixFreq)
{
  int i,j,k;
  int i0,j0,k0,l0,m0;
  int16_t *dataOut;
  float *dataOutF;
  int nchan;
  int npol=4;
  int nsub;
  int nbin=32;
  int status=0;
  int colnum_data;
  long naxes[4];
  int naxis=3;
  char tdim[16];
  double *raSub,*decSub;
  int colnum_datFreq, colnum_datWts,colnum_datScl,colnum_datOffs,colnum_raSub,colnum_decSub;
  int colnum_tsubint;
  char dataName[1024];
  hid_t dataset_id;
  float *datScl,*datOffs;
  float *data;
  float *floatArr;
  int ndump,nd;
  int nchanBand;
  long int n=0;
  int nband = inFile->beam[0].nBand;
  float sumFreq,chanbw,f0,f1;
  float cfreq,bw;

  double minShort,maxShort;
  double minVal,maxVal,tv;
  
  minShort = INT16_MIN;
  maxShort = INT16_MAX;

  printf("SHORT range = %g %g\n",minShort,maxShort);
  nchan=0;
  ndump = 0;
  for (j=0;j<nband;j++)
    {
      nchan += inFile->beam[0].calBandHeader[j].nchan;
      ndump += inFile->beam[0].calBandHeader[j].ndump;
    }
  nchanBand = inFile->beam[0].calBandHeader[0].nchan; // ASSUME CONSTANT
  
  printf("cal nchan = %d, ndump = %d, nchanBand = %d\n",nchan,ndump,nchanBand);  
  dataOut = (int16_t *)malloc(sizeof(int16_t)*nchan*nbin*npol);
  dataOutF = (float *)malloc(sizeof(float)*nchan*nbin*npol*ndump);
  datScl = (float *)malloc(sizeof(float)*nchan*npol);
  datOffs = (float *)malloc(sizeof(float)*nchan*npol);
  fits_movnam_hdu(cal,BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_get_colnum(cal, CASEINSEN, "RA_SUB", &colnum_raSub, &status);
  fits_get_colnum(cal, CASEINSEN, "DEC_SUB", &colnum_decSub, &status);
  fits_get_colnum(cal, CASEINSEN, "DATA", &colnum_data, &status);
  fits_get_colnum(cal, CASEINSEN, "TSUBINT", &colnum_tsubint, &status);
  if (status) {fits_report_error(stderr, status); exit(1);}

  fits_modify_vector_len(cal,colnum_data,(long)((long)nchan*(long)npol*(long)nbin),&status);

  fits_get_colnum(cal, CASEINSEN, "DAT_FREQ", &colnum_datFreq, &status);
  fits_modify_vector_len(cal,colnum_datFreq,nchan,&status);
  fits_get_colnum(cal, CASEINSEN, "DAT_WTS", &colnum_datWts, &status);
  fits_modify_vector_len(cal,colnum_datWts,nchan,&status);
  fits_get_colnum(cal, CASEINSEN, "DAT_SCL", &colnum_datScl, &status);
  fits_modify_vector_len(cal,colnum_datScl,nchan*npol,&status);
  fits_get_colnum(cal, CASEINSEN, "DAT_OFFS", &colnum_datOffs, &status);
  fits_modify_vector_len(cal,colnum_datOffs,nchan*npol,&status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  
  naxes[0] = nbin;
  naxes[1] = nchan;
  naxes[2] = npol;
  
  sprintf(tdim,"TDIM%d",colnum_data);
  fits_delete_key(cal, tdim, &status);
  fits_write_tdim(cal, colnum_data, naxis, naxes, &status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  printf("Complete re-sizeing file\n");
  floatArr = (float *)malloc(sizeof(float)*nchan*npol);
  for (j0 = 0; j0 < nband; j0++)
    {
      sprintf(dataName,"beam_%d/%s/calibrator_data/cal32_data",0,inFile->beam[0].bandHeader[j0].label);
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      nd = inFile->beam[0].calBandHeader[0].ndump;
      data = (float *)malloc(sizeof(float)*nchan*npol*nd*nbin);

      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);  
      
      for (l0 = 0; l0 <  inFile->beam[0].calBandHeader[j0].ndump; l0++)
	{
	  for (k0 = 0; k0 < inFile->beam[0].calBandHeader[j0].nchan; k0++)
	    {
	      for (i0 = 0;i0 < npol;i0++)
		{
		  for (m0=0;m0<nbin;m0++)
		    dataOutF[l0*nchan*npol*nbin + i0*nchan*nbin + (j0*nchanBand+k0)*nbin + m0] =
		      data[l0*nchanBand*npol*nbin + i0*nbin*nchanBand + k0*nbin + m0]; // Could be sped up
		}
	    }
	}
      //      printf("%g %d\n",data[i],dataOut[i]);
      status = H5Dclose(dataset_id);
      free(data);
    }
  printf("Writing cal data to PSRFITS\n");

  f0 = 1e6; f1 = -1e6;
  sumFreq = 0;
  raSub  = (double *)malloc(sizeof(double)*inFile->beam[0].calBandHeader[0].ndump);
  decSub = (double *)malloc(sizeof(double)*inFile->beam[0].calBandHeader[0].ndump);
  
  for (i = 0; i <  inFile->beam[0].calBandHeader[0].ndump; i++) // Assume all bands have same number of dumps
    {
      raSub[i] = inFile->beam[0].bandData[0].cal_obsHeader[i].raDeg;
      decSub[i] = inFile->beam[0].bandData[0].cal_obsHeader[i].raDeg;
      printf("ra/dec = %g %g\n",raSub[i],decSub[i]);
      // Convert to SHORT and determine DAT_SCL and DAT_OFFS
      for (j0=0;j0<npol;j0++)
	{
	  for (i0=0;i0<nchan;i0++)
	    {
	      for (k0=0;k0<nbin;k0++)
		{
		  tv = dataOutF[i*nchan*nbin*npol + j0*nchan*nbin + i0*nbin + k0];
		  if (k0==0)
		    minVal = maxVal = tv;
		  else
		    {
		      if (minVal > tv) minVal = tv;
		      if (maxVal < tv) maxVal = tv;
		    }
		}
	      //	      printf("(%d,%d) = %g %g\n",j0,j0,minVal,maxVal);
	      datScl[j0*nchan+i0]  = (maxVal-minVal)/(maxShort-minShort);
	      datOffs[j0*nchan+i0] = maxVal-maxShort*datScl[j0*nchan+i0];
	    }
	}
      for (j0=0;j0<npol;j0++)
	{
	  for (i0=0;i0<nchan;i0++)
	    {
	      for (k0=0;k0<nbin;k0++)
		dataOut[j0*nchan*nbin + i0*nbin + k0] = (int16_t)((dataOutF[i*nchan*nbin*npol + j0*nchan*nbin + i0*nbin+k0] - datOffs[j0*nchan+i0])/datScl[j0*nchan+i0]);
	    }
	}
      
      printf("Writing new row: %d\n",i);
      //      if (i>0)
      //	fits_insert_rows(cal,i,1,&status);
      fits_write_col(cal,TSHORT,colnum_data,i+1,1,nchan*nbin*npol,dataOut,&status);  
      
      printf("Complete writing data: sttus = %d\n",status);
      if (status) {fits_report_error(stderr, status); exit(1);}
      fits_write_col(cal,TFLOAT,colnum_datScl,i+1,1,nchan*npol,datScl,&status);  
      fits_write_col(cal,TFLOAT,colnum_datOffs,i+1,1,nchan*npol,datOffs,&status);  
    
      for (j=0;j<nchan;j++) floatArr[j] = 1; // WTS
      fits_write_col(cal,TFLOAT,colnum_datWts,i+1,1,nchan,floatArr,&status);  
      printf("Status here = %d\n",status);
      
      // Frequency axis
      printf("Frequency axis\n");
      n=0;

      for (k=0;k<inFile->beam[0].nBand;k++)
	{
	  sdhdf_loadBandData(inFile,0,k,2);
	  sdhdf_loadBandData(inFile,0,k,3);
	  
	  //	  printf("k = %d %d %d\n",k,n,inFile->beam[0].calBandHeader[k].nchan);
	  for (j=0;j<inFile->beam[0].calBandHeader[k].nchan;j++)
	    {
	      if (fixFreq==0)
		floatArr[n] = inFile->beam[0].bandData[k].cal_on_data.freq[j];  // FIX ME -- FREQ AXIS FOR DUMP
	      else
		{
		  double cbw,fs,fe;
		  fs  = inFile->beam[0].bandData[k].cal_on_data.freq[0];  // FIX ME -- FREQ AXIS FOR DUMP
		  fe  = inFile->beam[0].bandData[k].cal_on_data.freq[inFile->beam[0].calBandHeader[k].nchan-1];  // FIX ME -- FREQ AXIS FOR DUMP
		  cbw = (fe-fs)/inFile->beam[0].calBandHeader[k].nchan;
		  floatArr[n] = (fs + 0.5*cbw) + j*cbw;
		}
	      sumFreq += floatArr[n];
	      if (f0 > floatArr[n]) f0 = floatArr[n];
	      if (f1 < floatArr[n]) f1 = floatArr[n];
	      n++;
	    }
	  if (k==0)
	    chanbw = floatArr[1] - floatArr[0]; // Assuming constant
	  sdhdf_releaseBandData(inFile,0,k,2);
	  sdhdf_releaseBandData(inFile,0,k,3); 		      

	}
      printf("nchan and n = %d %d\n",nchan,n);
      if (i==0)
	{
	  cfreq = sumFreq/n;
	  bw = nchan*chanbw;
	}
      fits_write_col(cal,TFLOAT,colnum_datFreq,i+1,1,nchan,floatArr,&status);  
      fits_write_col(cal,TDOUBLE,colnum_raSub,i+1,1,1,&raSub[i],&status);
      fits_write_col(cal,TDOUBLE,colnum_decSub,i+1,1,1,&raSub[i],&status);
      fits_write_col(cal,TDOUBLE,colnum_tsubint,i+1,1,1,&inFile->beam[0].calBandHeader[0].dtime,&status);  
      printf("status at bottom = %d, cfreq = %g, bw = %g, chbw = %g\n",status,cfreq,bw,chanbw);
      
    }

  free(dataOut);
  free(dataOutF);
  free(datScl);
  free(datOffs);
  free(floatArr);
  free(raSub);
  free(decSub);
  fits_update_key(cal,TFLOAT,"CHAN_BW",&chanbw,NULL,&status);
  // Update the primary header
  fits_movabs_hdu(cal,1,NULL,&status);
  fits_update_key(cal,TFLOAT,"OBSFREQ",&cfreq,NULL,&status);
  fits_update_key(cal,TFLOAT,"OBSBW",&bw,NULL,&status);
  fits_update_key(cal,TINT,"OBSNCHAN",&nchan,NULL,&status);   
}

void createSubintHeader(fitsfile *cal,fitsfile *psr,  sdhdf_fileStruct *inFile)
{
  int i,j;
  int status=0;
  char name[128];
  char str[128];
  int iint;
  float ifloat;
  int nchan;
  
  strcpy(name,"SUBINT");
  fits_movnam_hdu(cal,BINARY_TBL,name,0,&status);fits_report_error(stderr,status);
  fits_movnam_hdu(psr,BINARY_TBL,name,0,&status);fits_report_error(stderr,status);
  iint = 4;
  fits_update_key(cal,TINT,"NPOL",&iint,NULL,&status);  fits_update_key(psr,TINT,"NPOL",&iint,NULL,&status);
  iint = 1;
  fits_update_key(cal,TINT,"NSBLK",&iint,NULL,&status);  fits_update_key(psr,TINT,"NSBLK",&iint,NULL,&status);
  fits_update_key(cal,TINT,"NBITS",&iint,NULL,&status);  fits_update_key(psr,TINT,"NBITS",&iint,NULL,&status);
  ifloat = 1; // FIX ME
  printf("TBIN IS WRONG!\n");
  fits_update_key(cal,TFLOAT,"TBIN",&ifloat,NULL,&status);
  fits_update_key(psr,TFLOAT,"TBIN",&ifloat,NULL,&status);
  ifloat = 0;
  fits_update_key(cal,TFLOAT,"DM",&ifloat,NULL,&status);  fits_update_key(psr,TFLOAT,"DM",&ifloat,NULL,&status);
  fits_update_key(cal,TFLOAT,"RM",&ifloat,NULL,&status);  fits_update_key(psr,TFLOAT,"RM",&ifloat,NULL,&status);

  fits_update_key(cal, TSTRING, (char *)"POL_TYPE", (char *)"AABBCRCI", NULL, &status);  fits_update_key(psr, TSTRING, (char *)"POL_TYPE", (char *)"AABBCRCI", NULL, &status);
  iint  =32; // For noise source
  fits_update_key(cal,TINT,"NBIN",&iint,NULL,&status);
  iint = 1; // For continuum source
  fits_update_key(psr,TINT,"NBIN",&iint,NULL,&status);


  
  // NCHAN for noise source
  nchan=0;
  for (j=0;j<inFile->beam[0].nBand;j++)
    nchan += inFile->beam[0].calBandHeader[j].nchan;
  
  iint = nchan;
  printf("cal nchan = %d\n",iint);  
  fits_update_key(cal,TINT,"NCHAN",&iint,NULL,&status);

  // NCHAN for astro
  nchan=0;
  for (j=0;j<inFile->beam[0].nBand;j++)
    nchan += inFile->beam[0].bandHeader[j].nchan;
  
  iint = nchan;
  printf("src nchan = %d\n",iint);  
  fits_update_key(psr,TINT,"NCHAN",&iint,NULL,&status);
}


void createPrimaryHeader(fitsfile *cal,fitsfile *psr,  sdhdf_fileStruct *inFile)
{
  int status=0;
  char str[128];
  int iint;
  float ifloat;
  int newHDUtype;
  
  // Delete unwanted tables
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"BANDPASS" , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"COHDDISP" , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"PSRPARAM" , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"POLYCO"   , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"T2PREDICT", 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"FLUX_CAL" , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"CAL_POLN" , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"FEEDPAR"  , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"DIG_STAT"  , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(cal, BINARY_TBL, (char *)"DIG_CNTS"  , 0, &status);   fits_delete_hdu(cal, &newHDUtype, &status);   fits_report_error(stdout,status);

  fits_movnam_hdu(psr, BINARY_TBL, (char *)"BANDPASS" , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"COHDDISP" , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"PSRPARAM" , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"POLYCO"   , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"T2PREDICT", 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"FLUX_CAL" , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"CAL_POLN" , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"FEEDPAR"  , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"DIG_STAT"  , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(psr, BINARY_TBL, (char *)"DIG_CNTS"  , 0, &status);   fits_delete_hdu(psr, &newHDUtype, &status);   fits_report_error(stdout,status);
  
  
  fits_movabs_hdu(cal,1,NULL,&status);
  fits_movabs_hdu(psr,1,NULL,&status);
  fits_write_date(cal,&status);
  fits_write_date(psr,&status);
  fits_update_key(cal, TSTRING, (char *)"OBSERVER", inFile->primary[0].observer, NULL, &status);  fits_update_key(psr, TSTRING, (char *)"OBSERVER", inFile->primary[0].observer, NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"PROJID", inFile->primary[0].pid, NULL, &status);    fits_update_key(psr, TSTRING, (char *)"PROJID", inFile->primary[0].pid, NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"TELESCOP", inFile->primary[0].telescope, NULL, &status);    fits_update_key(psr, TSTRING, (char *)"TELESCOP", inFile->primary[0].telescope, NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"BACKEND", (char *)"MEDUSA", NULL, &status);    fits_update_key(psr, TSTRING, (char *)"BACKEND", (char *)"MEDUSA", NULL, &status); // FIX MEDUSA HARDCODE
  fits_update_key(cal, TSTRING, (char *)"FRONTEND", inFile->primary[0].rcvr, NULL, &status);    fits_update_key(psr, TSTRING, (char *)"FRONTEND", inFile->primary[0].rcvr, NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"OBS_MODE", (char *)"CAL", NULL, &status);    fits_update_key(psr, TSTRING, (char *)"OBS_MODE", (char *)"PSR", NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"FD_POLN", (char *)"LIN", NULL, &status);    fits_update_key(psr, TSTRING, (char *)"FD_POLN", (char *)"LIN", NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"SRC_NAME", inFile->beamHeader[0].source, NULL, &status);    fits_update_key(psr, TSTRING, (char *)"SRC_NAME", inFile->beamHeader[0].source, NULL, &status);
  printf("RA = %s\n", inFile->beam[0].bandData[0].astro_obsHeader[0].raStr);
  fits_update_key(cal, TSTRING, (char *)"RA", inFile->beam[0].bandData[0].astro_obsHeader[0].raStr, NULL, &status);   fits_update_key(psr, TSTRING, (char *)"RA", inFile->beam[0].bandData[0].astro_obsHeader[0].raStr, NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"STT_CRD1", inFile->beam[0].bandData[0].astro_obsHeader[0].raStr, NULL, &status);   fits_update_key(psr, TSTRING, (char *)"STT_CRD1", inFile->beam[0].bandData[0].astro_obsHeader[0].raStr, NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"DEC", inFile->beam[0].bandData[0].astro_obsHeader[0].decStr, NULL, &status); fits_update_key(psr, TSTRING, (char *)"DEC", inFile->beam[0].bandData[0].astro_obsHeader[0].decStr, NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"STT_CRD2", inFile->beam[0].bandData[0].astro_obsHeader[0].decStr, NULL, &status); fits_update_key(psr, TSTRING, (char *)"STT_CRD2", inFile->beam[0].bandData[0].astro_obsHeader[0].decStr, NULL, &status);
  fits_update_key(cal, TSTRING, (char *)"COORD_MD", (char *)"J2000", NULL, &status); fits_update_key(psr, TSTRING, (char *)"COORD_MD", (char *)"J2000", NULL, &status);

  iint = (int)(inFile->beam[0].bandData[0].astro_obsHeader[0].mjd);
  fits_update_key(cal,TINT,"STT_IMJD",&iint,NULL,&status);  fits_update_key(psr,TINT,"STT_IMJD",&iint,NULL,&status);
  ifloat = (inFile->beam[0].bandData[0].astro_obsHeader[0].mjd - iint)*86400.0;
  fits_update_key(cal,TFLOAT,"STT_SMJD",&ifloat,NULL,&status);  fits_update_key(psr,TFLOAT,"STT_SMJD",&ifloat,NULL,&status);
  ifloat = 0.0;
  fits_update_key(cal,TFLOAT,"STT_OFFS",&ifloat,NULL,&status);  fits_update_key(psr,TFLOAT,"STT_OFFS",&ifloat,NULL,&status);
  ifloat = 2000.0;
  fits_update_key(cal,TFLOAT,"EQUINOX",&ifloat,NULL,&status);  fits_update_key(psr,TFLOAT,"EQUINOX",&ifloat,NULL,&status);

  // CAL specific
  fits_update_key(cal, TSTRING, (char *)"CAL_MODE", (char *)"SYNC", NULL, &status);
  ifloat = 100; // FIX ME
  fits_update_key(cal,TFLOAT,"CAL_FREQ",&ifloat,NULL,&status); 
  ifloat = 0.5; // FIX ME
  fits_update_key(cal,TFLOAT,"CAL_DCYC",&ifloat,NULL,&status); 
  ifloat = 0.0; // FIX ME
  fits_update_key(cal,TFLOAT,"CAL_PHS",&ifloat,NULL,&status); 

}
