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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "inspecta.h"
#include "hdf5.h"

#define VNUM "v0.1"

double dms_turn(char *line);
double hms_turn(char *line);
int    turn_hms(double turn, char *hms);
int    turn_dms(double turn, char *dms);
double turn_deg(double turn);
void convertGalactic(double raj,double decj,double *gl,double *gb);
void slaCaldj ( int iy, int im, int id, double *djm, int *j );
void slaCldj ( int iy, int im, int id, double *djm, int *j );

void help()
{
  printf("sdhdf_convertFromSDFITS\n\n");
  printf("Routine to convert SDFITS files to SDHDF\n");
  printf("\n\n");
  printf("-f <filename>        FITS filename for conversion\n");
  printf("-h                   This help\n");
  printf("-o <filename>        Output filename\n");

}

int main(int argc,char *argv[])
{
  int i,j,k,b;

  int fitsType=1; // 1 = SDFITS  
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

  char projectID[1024];
  char telescope[1024];
  char raStr[1024];
  char decStr[1024];
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
  float *spectrumData;
  float  *infloatVals;
  float  *freqVals;
  double *raVals,*decVals;
  double *timeVals;
  float *azimuth,*elevation;
  float *parAngle;
  float *exposure;
  float  *freqCalVals;
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

  int pulseOff1,pulseOff2;
  double *baseline1,*baseline2,*baseline3,*baseline4;

  int nbeam;
  int nband;
  int *beamNum;
  int *cycleNum;
  int *ifNum;
  char **object;
  char **dte;
  char n_str[16]="NULL";
  //
  int haveBeam=0;
  int n_ival;

  char **ctype;
  char readStr[1024];

  int chVal0;
  double val0;
  double dval;
  int colnum1,colnum2,colnum3;
  int nscan;

  int hr,min;
  double sec;
  int iy,im,id;
  double mjd0;
  int ret;
  
  // Read the command line arguments
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outname,argv[++i]);
    }


  // Read basic parameters from the SDFITS file
  printf("Loaded %s\n",fname);
  fits_open_file(&fptr,fname,READONLY,&status);
  fits_report_error(stderr,status);
  fits_movnam_hdu(fptr,BINARY_TBL,"SINGLE DISH",1,&status);

  // Read basic header information
  fits_read_key(fptr,TSTRING,"TELESCOP",&telescope,NULL,&status);
  fits_read_key(fptr,TSTRING,"PROJID",&projectID,NULL,&status);
  printf("Telescope = %s\n",telescope);
  fits_get_num_rows(fptr,&nrows,&status);
  printf("nrows = %d\n",nrows);

  beamNum = (int *)malloc(sizeof(int)*nrows);
  cycleNum = (int *)malloc(sizeof(int)*nrows);
  ifNum = (int *)malloc(sizeof(int)*nrows);
  object = (char **)malloc(sizeof(char *)*nrows);
  dte = (char **)malloc(sizeof(char *)*nrows);
  for (i=0;i<nrows;i++)
    {
      object[i] = (char *)malloc(sizeof(char)*16);
      dte[i]    = (char *)malloc(sizeof(char)*10);
    }
  // We now need to determine how those rows are divided
  fits_get_colnum(fptr,CASEINSEN,"CYCLE",&colnum,&status);
  fits_read_col(fptr,TINT,colnum,1,1,nrows,&n_ival,cycleNum,&initflag,&status);  

  fits_get_colnum(fptr,CASEINSEN,"BEAM",&colnum,&status);
  fits_read_col(fptr,TINT,colnum,1,1,nrows,&n_ival,beamNum,&initflag,&status);  

  fits_get_colnum(fptr,CASEINSEN,"IF",&colnum,&status);
  fits_read_col(fptr,TINT,colnum,1,1,nrows,&n_ival,ifNum,&initflag,&status);  

  fits_get_colnum(fptr,CASEINSEN,"OBJECT",&colnum,&status);
  fits_read_col(fptr,TSTRING,colnum,1,1,nrows,&n_str,object,&initflag,&status);  

  fits_get_colnum(fptr,CASEINSEN,"DATE-OBS",&colnum,&status);
  fits_read_col(fptr,TSTRING,colnum,1,1,nrows,&n_str,dte,&initflag,&status);  

  // Note: currently assuming a single object in the SDFITS file

  
  
  printf("SDFITS summary\n\n");
  printf("Scan  Cycle Beam IF  Object\n");
  nbeam=0;
  ndump=0;
  nband=0;
  for (i=0;i<nrows;i++)
    {
      printf("%-5.5d %-5.5d %-4.4d %-3.3d %-16.16s\n",i,cycleNum[i],beamNum[i],ifNum[i],object[i]);
      if (beamNum[i] > nbeam) nbeam=beamNum[i];  // This assumes that all beams are present in this file ** FIX ME
      if (cycleNum[i] > ndump) ndump=cycleNum[i];
      if (ifNum[i] > nband) nband =ifNum[i]; 
    }
  printf("==============================\n");
  printf("Number of spectral dumps/cycles = %d\n",ndump);
  printf("Number of beams = %d\n",nbeam);
  printf("Number of bands/IFs = %d\n",nband);
  printf("==============================\n");

  // Obtain the dimensions of the DATA array
  // Use the CTYPE to determine what the dimensions are
  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
  fits_read_tdim(fptr,colnum,maxdim,&naxis,naxes,&status);  
  fits_report_error(stderr,status);

  printf("Data array dimensions: %d\n",naxis);
  ctype = (char **)malloc(sizeof(char *)*naxis);
  for (i=0;i<naxis;i++)
    {
      ctype[i] = (char *)malloc(sizeof(char)*16);
      sprintf(readStr,"CTYPE%d",i+1);
      fits_read_key(fptr,TSTRING,readStr,ctype[i],NULL,&status);
    }

  nchan = -1;
  npol  = -1;
  for (i=0;i<naxis;i++)
    {
      if (strcmp(ctype[i],"FREQ")==0)
	nchan = naxes[i];
      else if (strcmp(ctype[i],"STOKES")==0)
	npol = naxes[i];
      printf(".... %d %d %s\n",i,naxes[i],ctype[i]);
    }
  if (nchan==-1)
    {
      printf("ERROR: Cannot read NCHAN\n");
      exit(1);
    }
  else if (npol==-1)
    {
      printf("ERROR: Cannot read NPOL\n");
      exit(1);
    }
  printf("nchan = %d\n",nchan);
  printf("npol = %d\n",npol);
  printf("==============================\n");

  //
  // Note can also read in further dimensions (e.g., N_RA, N_DEC) -- FIX ME
  //

  // Load in all the data
  floatVals   = (float *)malloc(sizeof(float)*nchan*npol*ndump*nbeam);  
  freqVals    = (float *)malloc(sizeof(float)*nchan*ndump*nbeam);
  timeVals    = (double *)malloc(sizeof(double)*ndump*nbeam); // Should include *nRA
  raVals    = (double *)malloc(sizeof(double)*ndump*nbeam); // Should include *nRA
  decVals    = (double *)malloc(sizeof(double)*ndump*nbeam); // Should include *nDEC
  azimuth    = (float *)malloc(sizeof(float)*ndump*nbeam); 
  elevation    = (float *)malloc(sizeof(float)*ndump*nbeam);
  parAngle    = (float *)malloc(sizeof(float)*ndump*nbeam);
  exposure    = (float *)malloc(sizeof(float)*ndump*nbeam); 


  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);

  // Assuming only NCHANxNPOL for each scan here -- FIX ME
  for (i=0;i<nrows;i++)
    fits_read_col(fptr,TFLOAT,colnum,i+1,1,nchan*npol,&n_fval,floatVals+i*nchan*npol,&initflag,&status);

  // Read frequency channels
  // Assuming 1 = frequency column -- FIX ME
  // Assuming frequency values in Hz -- FIX ME
  fits_get_colnum(fptr,CASEINSEN,"CRPIX1",&colnum1,&status);
  fits_get_colnum(fptr,CASEINSEN,"CRVAL1",&colnum2,&status);
  fits_get_colnum(fptr,CASEINSEN,"CDELT1",&colnum3,&status);
  for (i=0;i<nrows;i++)
    {      
      fits_read_col(fptr,TINT,colnum1,i+1,1,1,&n_ival,&chVal0,&initflag,&status);
      fits_read_col(fptr,TDOUBLE,colnum2,i+1,1,1,&n_dval,&val0,&initflag,&status);
      fits_read_col(fptr,TDOUBLE,colnum3,i+1,1,1,&n_dval,&dval,&initflag,&status);
      for (j=0;j<nchan;j++)
	freqVals[i*nchan+j] = (val0 + ((j+1)-chVal0)*dval)/1.0e6; // Convert to MHz
    }

  
  // Assume RA is simply CRVAL3 -- FIX ME
  fits_get_colnum(fptr,CASEINSEN,"CRVAL3",&colnum,&status);
  fits_read_col(fptr,TDOUBLE,colnum,1,1,nrows,&n_dval,raVals,&initflag,&status);
  
  // Assume DEC is simply CRVAL4 -- FIX ME
  fits_get_colnum(fptr,CASEINSEN,"CRVAL4",&colnum,&status);
  fits_read_col(fptr,TDOUBLE,colnum,1,1,nrows,&n_dval,decVals,&initflag,&status);

  fits_get_colnum(fptr,CASEINSEN,"TIME",&colnum,&status);
  fits_read_col(fptr,TDOUBLE,colnum,1,1,nrows,&n_dval,timeVals,&initflag,&status);

 
  fits_get_colnum(fptr,CASEINSEN,"AZIMUTH",&colnum,&status);
  fits_read_col(fptr,TFLOAT,colnum,1,1,nrows,&n_fval,azimuth,&initflag,&status);

  fits_get_colnum(fptr,CASEINSEN,"ELEVATIO",&colnum,&status);
  fits_read_col(fptr,TFLOAT,colnum,1,1,nrows,&n_fval,elevation,&initflag,&status);

  fits_get_colnum(fptr,CASEINSEN,"PARANGLE",&colnum,&status);
  fits_read_col(fptr,TFLOAT,colnum,1,1,nrows,&n_fval,parAngle,&initflag,&status);

  fits_get_colnum(fptr,CASEINSEN,"EXPOSURE",&colnum,&status);
  fits_read_col(fptr,TFLOAT,colnum,1,1,nrows,&n_fval,exposure,&initflag,&status);

  
  fits_close_file(fptr,&status);
  
  // -----------------------------------------------------------------------
  // Now produce the SDHDF file
  // -----------------------------------------------------------------------
  primaryHeader    = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct));
  beamHeader       = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nbeam);
  bandHeader       = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct));
  softwareVersions = (sdhdf_softwareVersionsStruct *)malloc(sizeof(sdhdf_softwareVersionsStruct));
  history          = (sdhdf_historyStruct *)malloc(sizeof(sdhdf_historyStruct));

  sdhdf_setMetadataDefaults(primaryHeader,beamHeader,bandHeader,softwareVersions,history,1,1);
  primaryHeader[0].nbeam = nbeam;
  hr = (int)(timeVals[0]/60./60.);
  min = (int)((timeVals[0]-hr*60*60.)/60.);
  sec = timeVals[0] - hr*60*60. - min*60;
  sprintf(primaryHeader[0].utc0,"%s-%02d:%02d:%02d",dte[0],hr,min,(int)sec); // Default to first dump/beam. 
  strcpy(primaryHeader[0].cal_mode,"OFF"); // FIX ME
  if (strcmp(telescope,"ATPKSMB")==0)
    {
      strcpy(primaryHeader[0].telescope,"Parkes");
      strcpy(primaryHeader[0].rcvr,"MB");
    }
  else    
    strcpy(primaryHeader[0].telescope,telescope);
  strcpy(primaryHeader[0].pid,projectID);
  
  for (i=0;i<nbeam;i++)
    {
      sprintf(beamHeader[i].label,"beam_%d",i);
      beamHeader[i].nBand = 1;
      strcpy(beamHeader[i].source,object[i]);	  
    }
  obsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);

  // Set up the band information
  strcpy(bandHeader->label,"band0");
  bandHeader->fc = (freqVals[0]+freqVals[nchan-1])/2.0;  // UPDATE
  bandHeader->f0 = freqVals[0]; // UPDATE
  bandHeader->f1 = freqVals[nchan-1];  // UPDATE
  bandHeader->nchan = nchan;
  bandHeader->npol = npol;
  strcpy(bandHeader->pol_type,"AABB"); // FIX ME
  bandHeader->dtime = exposure[0]; // FIX
  bandHeader->ndump = ndump;

  // Open the output HDF5 file
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }
  sdhdf_initialiseFile(outFile);

  sdhdf_openFile(outname,outFile,3);
  sdhdf_writePrimaryHeader(outFile,primaryHeader);
  sdhdf_writeBeamHeader(outFile,beamHeader,nbeam);
  spectrumData = (float *)malloc(sizeof(float)*nchan*npol*ndump);
  for (b=0;b<nbeam;b++)
    {
      printf("Creating observation parameters for beam %d\n",b);
      for (i=0;i<ndump;i++)
	{
	  obsParams[i].timeElapsed = (i+0.5)*exposure[0]; // FIX in case not all exposures are the same length
	  strcpy(obsParams[i].timedb,"UNKNOWN"); // FIX
	  sscanf(dte[i*nbeam+b],"%d-%d-%d",&iy,&im,&id);
	  slaCaldj (iy, im, id, &mjd0, &ret);
	  obsParams[i].mjd = mjd0 + timeVals[i*nbeam+b]/86400.; 
	  hr = (int)(timeVals[i*nbeam+b]/60./60.);
	  min = (int)((timeVals[i*nbeam+b]-hr*60*60.)/60.);
	  sec = timeVals[i*nbeam+b] - hr*60*60. - min*60;
	  sprintf(obsParams[i].utc,"%02d:%02d:%02d",hr,min,(int)sec);
	  strcpy(obsParams[i].ut_date,"UNKNOWN"); // FIX
	  strcpy(obsParams[i].local_time,"UNKNOWN"); // FIX 
	  
	  turn_hms(raVals[i*nbeam+b]/360.,raStr); // FIX: This assumes the specific structure of beams in the RA list
	  turn_dms(decVals[i*nbeam+b]/360.,decStr);	  
	  strcpy(obsParams[i].raStr,raStr);
	  strcpy(obsParams[i].decStr,decStr);
	  
	  obsParams[i].raDeg = raVals[i*nbeam+b]; // FIX: This assumes the specific structure of beams in the RA list
	  obsParams[i].decDeg = decVals[i*nbeam+b];
	  obsParams[i].raOffset = 0;
	  obsParams[i].decOffset = 0;
	  convertGalactic(obsParams[i].raDeg,obsParams[i].decDeg,&(obsParams[i].gl),&(obsParams[i].gb));
	  obsParams[i].az = azimuth[i*nbeam+b];
	  obsParams[i].ze = 90-elevation[i*nbeam+b];
	  obsParams[i].el = elevation[i*nbeam+b];
	  obsParams[i].az_drive_rate = 0;
	  obsParams[i].ze_drive_rate = 0;
	  obsParams[i].hourAngle = 0;
	  obsParams[i].paraAngle = parAngle[i*nbeam+b];
	  obsParams[i].windDir = 0;
	  obsParams[i].windSpd = 0;      
	}
      printf("Complete\n");

    
      sdhdf_writeBandHeader(outFile,bandHeader,beamHeader[b].label,1,1);
      sdhdf_writeObsParams(outFile,bandHeader[0].label,beamHeader[b].label,0,obsParams,ndump,1);

      nscan=0;
      for (k=0;k<nrows;k++)
	{
	  if (beamNum[k]==b+1)
	    {
	      for (i=0;i<nchan*npol;i++)
		spectrumData[nscan*nchan*npol+i] = floatVals[k*nchan*npol+i];	      
	      nscan++;
	    }
	}      
      sdhdf_writeSpectrumData(outFile,beamHeader[b].label,bandHeader->label,b,0,spectrumData,freqVals,nchan,npol,ndump,1,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
    }
  sdhdf_writeSoftwareVersions(outFile,softwareVersions);
  sdhdf_writeHistory(outFile,history,1);
  
  sdhdf_closeFile(outFile);
  free(outFile);

  

  free(spectrumData);
  free(floatVals);
  free(freqVals);
  free(timeVals);
  free(raVals);
  free(decVals);
  free(azimuth); free(elevation);
  free(parAngle);
  free(beamNum);
  free(cycleNum);
  free(ifNum);
  for (i=0;i<nrows;i++)
    {
      free(object[i]);
      free(dte[i]);
    }
  free(object);
  free(dte);
  for (i=0;i<naxis;i++)
    free(ctype[i]);
  free(ctype);
  exit(1);



  

   
  
  //      for (i=0;i<ndump;i++)
  // FIX ME
  ndump/=13;  // (as we have 13 beams)
  for (i=0;i<ndump;i++)
    {
      printf("Loading spectrum dump %d/%d\n",i,ndump-1);
      //	  fits_read_col(fptr,TFLOAT,colnum,i+1,1,nchan*npol,&n_fval,infloatVals,&initflag,&status);
      
      // FIX ME
      for (b=0;b<nbeam;b++)
	{
	  fits_read_col(fptr,TFLOAT,colnum,i*nbeam+b+1,1,nchan*npol,&n_fval,infloatVals,&initflag,&status);
	  // Rearrange ordering
	  for (j=0;j<npol;j++)
	    {
	      for (k=0;k<nchan;k++)
		floatVals[b*nchan*npol*ndump + i*nchan*npol + j*nchan + k]
		  = infloatVals[j*nchan + k];
	      
	      //floatVals[i*nchan*npol + j*nchan + k] = infloatVals[k*npol+j];
	    }
	}
    }

  printf("Closing file\n");


} 

double turn_deg(double turn){
 
  /* Converts double turn to string "sddd.ddd" */
  return turn*360.0;
}


int turn_dms(double turn, char *dms){
  
  /* Converts double turn to string "sddd:mm:ss.sss" */
  
  int dd, mm, isec;
  double trn, sec;
  char sign;
  
  sign=' ';
  if (turn < 0.){
    sign = '-';
    trn = -turn;
  }
  else{
    sign = '+';
    trn = turn;
  }
  dd = trn*360.;
  mm = (trn*360.-dd)*60.;
  sec = ((trn*360.-dd)*60.-mm)*60.;
  isec = (sec*1000. +0.5)/1000;
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        dd=dd+1;
      }
    }
  sprintf(dms,"%c%02d:%02d:%010.7f",sign,dd,mm,sec);
 
}


int turn_hms(double turn, char *hms){
 
  /* Converts double turn to string " hh:mm:ss.ssss" */
  
  int hh, mm, isec;
  double sec;

  hh = turn*24.;
  mm = (turn*24.-hh)*60.;
  sec = ((turn*24.-hh)*60.-mm)*60.;
  isec = (sec*10000. +0.5)/10000;
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        hh=hh+1;
        if(hh==24){
          hh=0;
        }
      }
    }

  sprintf(hms," %02d:%02d:%010.7f",hh,mm,sec);
 
}


double hms_turn(char *line){

  /* Converts string " hh:mm:ss.ss" or " hh mm ss.ss" to double turn */
  
  int i;int turn_hms(double turn, char *hms);
  double hr, min, sec, turn=0;
  char hold[MAX_STRLEN];

  strcpy(hold,line);

  /* Get rid of ":" */
  for(i=0; *(line+i) != '\0'; i++)if(*(line+i) == ':')*(line+i) = ' ';

  i = sscanf(line,"%lf %lf %lf", &hr, &min, &sec);
  if(i > 0){
    turn = hr/24.;
    if(i > 1)turn += min/1440.;
    if(i > 2)turn += sec/86400.;
  }
  if(i == 0 || i > 3)turn = 1.0;


  strcpy(line,hold);

  return turn;
}

double dms_turn(char *line){

  /* Converts string "-dd:mm:ss.ss" or " -dd mm ss.ss" to double turn */
  
  int i;
  char *ic, ln[40];
  double deg, min, sec, sign, turn=0;

  /* Copy line to internal string */
  strcpy(ln,line);

  /* Get rid of ":" */
  for(i=0; *(ln+i) != '\0'; i++)if(*(ln+i) == ':')*(ln+i) = ' ';

  /* Get sign */
  if((ic = strchr(ln,'-')) == NULL)
     sign = 1.;
  else {
     *ic = ' ';
     sign = -1.;
  }

  /* Get value */
  i = sscanf(ln,"%lf %lf %lf", &deg, &min, &sec);
  if(i > 0){
    turn = deg/360.;
    if(i > 1)turn += min/21600.;
    if(i > 2)turn += sec/1296000.;
    if(turn >= 1.0)turn = turn - 1.0;
    turn *= sign;
  }
  if(i == 0 || i > 3)turn =1.0;

  return turn;
}

void convertGalactic(double raj,double decj,double *gl,double *gb)
{
  double sinb,y,x,at;
  double rx,ry,rz,rx2,ry2,rz2;
  double deg2rad = M_PI/180.0;
  double gpoleRAJ = 192.85*deg2rad;
  double gpoleDECJ = 27.116*deg2rad;
  double rot[4][4];

  /* Note: Galactic coordinates are defined from B1950 system - e.g. must transform from J2000.0                      
                                                  
     equatorial coordinates to IAU 1958 Galactic coords */

  /* Convert to rectangular coordinates */
  rx = cos(raj)*cos(decj);
  ry = sin(raj)*cos(decj);
  rz = sin(decj);

  /* Now rotate the coordinate axes to correct for the effects of precession */
  /* These values contain the conversion between J2000 and B1950 and from B1950 to Galactic */
  rot[0][0] = -0.054875539726;
  rot[0][1] = -0.873437108010;
  rot[0][2] = -0.483834985808;
  rot[1][0] =  0.494109453312;
  rot[1][1] = -0.444829589425;
  rot[1][2] =  0.746982251810;
  rot[2][0] = -0.867666135858;
  rot[2][1] = -0.198076386122;
  rot[2][2] =  0.455983795705;

  rx2 = rot[0][0]*rx + rot[0][1]*ry + rot[0][2]*rz;
  ry2 = rot[1][0]*rx + rot[1][1]*ry + rot[1][2]*rz;
  rz2 = rot[2][0]*rx + rot[2][1]*ry + rot[2][2]*rz;

  /* Convert the rectangular coordinates back to spherical coordinates */
  *gb = asin(rz2);
  *gl = atan2(ry2,rx2);
  if (*gl < 0) (*gl)+=2.0*M_PI;
}

void slaCaldj ( int iy, int im, int id, double *djm, int *j ) /*includefile*/
/*
**  - - - - - - - - -
**   s l a C a l d j
**  - - - - - - - - -
**
**  Gregorian calendar to Modified Julian Date.
**
**  (Includes century default feature:  use slaCldj for years
**   before 100AD.)
**
**  Given:
**     iy,im,id   int      year, month, day in Gregorian calendar
**
**  Returned:
**     *djm       double   Modified Julian Date (JD-2400000.5) for 0 hrs
**     *j         int      status:
**                           0 = ok
**                           1 = bad year   (MJD not computed)
**                           2 = bad month  (MJD not computed)
**                           3 = bad day    (MJD computed)
**
**  Acceptable years are 00-49, interpreted as 2000-2049,
**                       50-99,     "       "  1950-1999,
**                       100 upwards, interpreted literally.
**
**  Called:  slaCldj
**
**  Last revision:   21 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int ny;

/* Default century if appropriate */
   if ( ( iy >= 0 ) && ( iy <= 49 ) )
      ny = iy + 2000;
   else if ( ( iy >= 50 ) && ( iy <= 99 ) )
      ny = iy + 1900;
   else
      ny = iy;

/* Modified Julian Date */
   slaCldj ( ny, im, id, djm, j );
}

void slaCldj ( int iy, int im, int id, double *djm, int *j ) /*includefile*/
/*
**  - - - - - - - -
**   s l a C l d j
**  - - - - - - - -
**
**  Gregorian calendar to Modified Julian Date.
**
**  Given:
**     iy,im,id     int    year, month, day in Gregorian calendar
**
**  Returned:
**     *djm         double Modified Julian Date (JD-2400000.5) for 0 hrs
**     *j           int    status:
**                           0 = OK
**                           1 = bad year   (MJD not computed)
**                           2 = bad month  (MJD not computed)
**                           3 = bad day    (MJD computed)
**
**  The year must be -4699 (i.e. 4700BC) or later.
**
**  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
**
**  Last revision:   29 August 1994
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   long iyL, imL;

/* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };



/* Validate year */
   if ( iy < -4699 ) { *j = 1; return; }

/* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { *j = 2; return; }

/* Allow for leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
             ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
             29 : 28;

/* Validate day */
   *j = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;

/* Lengthen year and month numbers to avoid overflow */
   iyL = (long) iy;
   imL = (long) im;

/* Perform the conversion */
   *djm = (double)
        ( ( 1461L * ( iyL - ( 12L - imL ) / 10L + 4712L ) ) / 4L
        + ( 306L * ( ( imL + 9L ) % 12L ) + 5L ) / 10L
        - ( 3L * ( ( iyL - ( 12L - imL ) / 10L + 4900L ) / 100L ) ) / 4L
        + (long) id - 2399904L );
}
