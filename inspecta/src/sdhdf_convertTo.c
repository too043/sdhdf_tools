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
#include "inspecta.h"
#include "hdf5.h"
#include "fitsio.h"



#define VNUM "v0.1"

#define MAX_SRC_NAMES 16

void help()
{
  printf("sdhdf_convertTo: %s\n",VNUM);
  printf("sdhfProc version:   %s\n",SOFTWARE_VER);
  printf("author:             George Hobbs\n");
  printf("\n");
  printf("\n\nCommand line arguments:\n\n");
  printf("-band <integer>     Sub-band selection (starting from 0 for the first subband in the file)\n");
  printf("-e <string>         Output file extension (defaults to sdfits)  (note that this is not used if the -o option is set)\n");
  printf("-h                  This help\n");
  printf("-o <string>         Output file name (this will automatically join input data sets together as one output file)\n");
  printf("-noXpol             Do not output cross polarisation terms\n");
  printf("-srcName <from> <to> Changes source names from <from> to <to>\n");
  printf("<filename>          SDHDF file corresponding to observations\n");
  printf("\n");
  printf("\nExample:\n\n");
  printf("sdhdf_convertTo <filename.hdf>\n");
  printf("---------------------\n");

  exit(1);
}

int main(int argc,char *argv[])
{
  int i,j,ii,k,kk;
  char fname[MAX_FILES][MAX_STRLEN];
  char oname[MAX_STRLEN];
  int  join=0;
  char requestName[MAX_STRLEN]="NULL";
  char ext[MAX_STRLEN] = "sdfits";

  int oneScanPerFile = 0;
  
  // SDHDF file
  sdhdf_fileStruct *inFile;
  int nFiles=0;

  // Parameters for fits file
  fitsfile *fptr;
  int       fitsStatus=0;
  float     fval;
  float    *fdata,*xdata;
  float     hr,min,sec;
  int       colnum;
  int       rowNum=0,ival;
  int       scan=0;
  int       scan0=0;
  long      naxes[4],naxes_xpol[4];
  float     intTime;  // Total integration time
  char    **strArray;
  char      timeStr[1024];

  int       haveXpol=0;
  int       useXpol=1;
  int       ibeam = 0;
  int       pnum = 0;

  int       selectBand=-1,i0,i1;
  long      nchan0;
  
  char      runtimeDir[1024];

  // Modification of source names as needed
  int       nSrcChange=0; 
  char      changeFrom[MAX_SRC_NAMES][64];
  char      changeTo[MAX_SRC_NAMES][64];
  char      srcName[64];

  
  if (argc==1) help();
  
  // Allocate memory for the input ifle
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  // Allocate memory for storing arrays of strings
  strArray = (char **)malloc(sizeof(char *));
  strArray[0] = (char *)malloc(sizeof(char)*1024);

  printf("WARNING: Various parameters not being set:\n"); //Should set FOCUSTAN correctly\n");
  printf("FOCUSROT, FOCUSTAN, SCANRATE, TSYS, CALFCTR, TCAL, TCALTIME\n");
  printf("TAMBIENT, PRESSURE and HUMIDTY being set to 0 as they currently are not recorded in the SDHDF file\n");
  for (i=1;i<argc;i++)
    {
      if (strcasecmp(argv[i],"-noXpol")==0)  // Ignore cross polarisation terms
	useXpol=0;
      else if (strcasecmp(argv[i],"-1scanPerFile")==0)   // 1 scan per file
	oneScanPerFile=1;
      else if (strcmp(argv[i],"-srcChange")==0)
	{
	  strcpy(changeFrom[nSrcChange],argv[++i]);
	  strcpy(changeTo[nSrcChange],argv[++i]);
	  nSrcChange++;	  
	}
      else if (strcmp(argv[i],"-e")==0)      // Output file extension
	strcpy(ext,argv[++i]);
      else if (strcmp(argv[i],"-band")==0)   // Selection of band
	sscanf(argv[++i],"%d",&selectBand);
      else if (strcmp(argv[i],"-o")==0)
	{strcpy(requestName,argv[++i]); join=1;}
      else if (strcmp(argv[i],"-h")==0)
	help();
      else
	strcpy(fname[nFiles++],argv[i]);
    }

  printf("Number of source name changes: %d\n",nSrcChange);

  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdf_convertTo requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  printf("Using runtime directory >%s<\n",runtimeDir);
  
  for (i=0;i<nFiles;i++)
    {
      printf("Processing: %s\n",fname[i]);
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_loadMetaData(inFile);
      if (strcmp(requestName,"NULL")==0)
	sprintf(oname,"!%s.%s(%s/fileConversion/sdHeader.fits)",fname[i],ext,runtimeDir);
      else
	sprintf(oname,"!%s(%s/fileConversion/sdHeader.fits)",requestName,runtimeDir);

      if (i==0 || join==0)
	{
	  fits_create_file(&fptr,oname,&fitsStatus);
	  fits_report_error(stderr, fitsStatus); 
      
	  // Primary header information in the SDFITS file
	  fits_movabs_hdu(fptr,1,NULL,&fitsStatus);
	  fits_write_date(fptr,&fitsStatus);
	  
	  // Header information in SINGLE_DISH table
	  fits_movnam_hdu(fptr, BINARY_TBL, (char *)"SINGLE DISH", 0, &fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"TELESCOP",inFile->primary[0].telescope,NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"INSTRUME",inFile->primary[0].instrument,NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"DATE-OBS",inFile->primary[0].utc0,NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"PROJID",inFile->primary[0].pid,NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"FRONTEND",inFile->primary[0].rcvr,NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"HDRVER",inFile->primary[0].hdr_defn_version,NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"OBSERVER",inFile->primary[0].observer,NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"EQUINOX",(char *)"2000",NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"RADESYS",(char *)"FK5",NULL,&fitsStatus);
	  fits_update_key(fptr,TSTRING,(char *)"SYSOBS",(char *)"TOPOCENT",NULL,&fitsStatus); // This assumes that we're working from the raw SDHDF file
      
	  // ASAP only allows [REST, LSRK, LSRD, BARY, GEO, TOPO, GALACTO, LGROUP, CMB, Undefined]
	  // Whereas the documentation allows (https://casa.nrao.edu/aips2_docs/notes/236/node14.html): -LSR,-HEL,-OBS,LSRK,-GEO,REST,-GAL
	  fits_update_key(fptr,TSTRING,(char *)"SPECSYS",(char *)"TOPO",NULL,&fitsStatus);
       
	  fits_update_key(fptr,TSTRING,(char *)"CTYPE2",(char *)"STOKES",NULL,&fitsStatus); // SHOULDN'T HARD CODE STOKES?? HERE
	  fval = 1.0; fits_update_key(fptr,TFLOAT,(char *)"CRPIX2",&fval,NULL,&fitsStatus);
	  fval = -5.0; fits_update_key(fptr,TFLOAT,(char *)"CRVAL2",&fval,NULL,&fitsStatus); // See https://casa.nrao.edu/aips2_docs/notes/236/node14.html
	  fval = -1.0; fits_update_key(fptr,TFLOAT,(char *)"CDELT2",&fval,NULL,&fitsStatus);
	  printf("Incorrectly setting CDELT3 and CDELT4\n");
	  fval = -1.0; fits_update_key(fptr,TFLOAT,(char *)"CDELT3",&fval,NULL,&fitsStatus);
	  fval = 1.0; fits_update_key(fptr,TFLOAT,(char *)"CDELT4",&fval,NULL,&fitsStatus);
	  
	  // These are taken from an old SDFITS file
	  //
	  fval = -4.554232087E+06; fits_update_key(fptr,TFLOAT,(char *)"OBSGEO-X",&fval,NULL,&fitsStatus);
	  fval =  2.816759046E+06; fits_update_key(fptr,TFLOAT,(char *)"OBSGEO-Y",&fval,NULL,&fitsStatus);
	  fval = -3.454035950E+06; fits_update_key(fptr,TFLOAT,(char *)"OBSGEO-Z",&fval,NULL,&fitsStatus);
	}

           
       
       //
       // ******************* Process each sub-band
       if (selectBand < 0)
	 {
	   i0 = 0;
	   i1 = inFile->beam[ibeam].nBand;
	 }
       else
	 {
	   printf("Selecting band %d\n",selectBand);
	   i0 = selectBand;
	   i1 = i0+1;
	 }
       
       for (ii=i0;ii<i1;ii++)
	 {
	   scan=0; // Set the scan number back to the start
      
	   printf("Processing sub-band %d (%s)\n",ii,inFile->beam[ibeam].bandHeader[ii].label);
	   printf("scan number = %d\n",scan+scan0);
	   if (ii==i0)
	     {
	       naxes[0] = inFile->beam[ibeam].bandHeader[i0].nchan; // **** FIX ME --- WE SHOULD ONLY SET THIS ONCE
	       if (useXpol == 1) // For some reason SDFITS seems to have DATA as (nchan,npol) and XDATA as (npol,nchan)!
		 {
		   naxes_xpol[0] = 2;
		   naxes_xpol[1] = inFile->beam[ibeam].bandHeader[i0].nchan;
		 }
	       if (inFile->beam[ibeam].bandHeader[i0].npol > 1)
		 {
		   naxes[1] = 2; // Xpol data in a different table
		   pnum = 2;
		   if (inFile->beam[ibeam].bandHeader[i0].npol == 4) // Do we have cross polar terms?
		     haveXpol = 1;
		 }
	       else
		 {
		   naxes[1] = 1;
		   pnum = 1;
		 }
	       naxes[2] = 1;
	       naxes[3] = 1;
	     }
	   if (strlen(inFile->beam[ibeam].bandHeader[ii].label)>0)
	     {
	       if (ii==i0)
		 nchan0=inFile->beam[ibeam].bandHeader[ii].nchan;
	       else
		 {
		   if (inFile->beam[ibeam].bandHeader[ii].nchan != nchan0)
		     {
		       printf("ERROR: SDFITS cannot accept different numbers of channels in different subbands.\n");
		       printf("sdhdf_convertTo therefore cannot convert this file with %d channels for subband 0 and %d channels for subband %d\n",(int)nchan0,(int)inFile->beam[ibeam].bandHeader[ii].nchan,(int)ii);
		       printf("Please re-run using the -band option to select a specific sub-band\n");
		       exit(1);
		     }
		 }
	       printf("Checking subintegrations\n");
	       printf("haveXpol = %d, useXpol = %d, pnum = %d\n",haveXpol,useXpol,pnum);
	       // ******************* Now process each subintegration
	       fdata = (float *)malloc(sizeof(float)*inFile->beam[ibeam].bandHeader[ii].nchan*pnum);
	       if (haveXpol==1 && useXpol==1)
		 xdata = (float *)malloc(sizeof(float)*inFile->beam[ibeam].bandHeader[ii].nchan*2);

	       if (ii==i0) 
		 {
		   char tdimName[128];

		   
		   fits_get_colnum(fptr,CASEINSEN,(char *)"DATA",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fits_modify_vector_len (fptr, colnum, inFile->beam[ibeam].bandHeader[ii].nchan*pnum, &fitsStatus); fits_report_error(stderr,fitsStatus);
		   sprintf(tdimName,"TDIM%d",colnum);
		   fits_delete_key(fptr, (char *)tdimName, &fitsStatus);
		   if (fitsStatus) { fits_report_error(stderr,fitsStatus); exit(1);}
		   fits_write_tdim(fptr,colnum,4,naxes,&fitsStatus);fits_report_error(stderr,fitsStatus);
		   if (fitsStatus) { fits_report_error(stderr,fitsStatus); exit(1);}
		   if (haveXpol == 1 && useXpol == 1)
		     {
		       fits_get_colnum(fptr,CASEINSEN,(char *)"XPOLDATA",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		       sprintf(tdimName,"TDIM%d",colnum);
		       fits_modify_vector_len (fptr, colnum, inFile->beam[ibeam].bandHeader[ii].nchan*pnum, &fitsStatus); fits_report_error(stderr,fitsStatus);
		       fits_delete_key(fptr, (char *)tdimName, &fitsStatus);
		       if (fitsStatus) { fits_report_error(stderr,fitsStatus); exit(1);}
		       fits_write_tdim(fptr,colnum,2,naxes_xpol,&fitsStatus);fits_report_error(stderr,fitsStatus);
		       if (fitsStatus) { fits_report_error(stderr,fitsStatus); exit(1);}
		     }
		 }


	       
	       sdhdf_loadBandData(inFile,ibeam,ii,1);
	      
	       intTime = inFile->beam[ibeam].bandHeader[ii].dtime * inFile->beam[ibeam].bandHeader[ii].ndump; /// GEORGE CHECK WHAT THIS SHOULD BE
	       
	       for (j=0;j<inFile->beam[ibeam].bandHeader[ii].ndump;j++)
		 {
		   for (k=0;k<inFile->beam[ibeam].bandHeader[ii].npol;k++)
		     {
		       for (kk=0;kk<inFile->beam[ibeam].bandHeader[ii].nchan;kk++)
			 {
			   if (k==0)
			     fdata[k*inFile->beam[ibeam].bandHeader[ii].nchan+kk] = inFile->beam[ibeam].bandData[ii].astro_data.pol1[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan];
			   else if (k==1)
			     fdata[k*inFile->beam[ibeam].bandHeader[ii].nchan+kk] = inFile->beam[ibeam].bandData[ii].astro_data.pol2[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan];
			   
			   // These are the cross polar terms
			   if (haveXpol==1 && useXpol==1)
			     {
			       // NOTE: -- The packing of the cross polarization terms seems to be (2,nchan), not (nchan, 2) and so this seems to be different from the data column above
			       //
			       if (k==2)
				 xdata[(kk*2)+(k-2)] = inFile->beam[ibeam].bandData[ii].astro_data.pol3[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan];
			       else if (k==3)
				 xdata[(kk*2)+(k-2)] = inFile->beam[ibeam].bandData[ii].astro_data.pol4[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan];
			     }
			 }
		     }

		   fits_insert_rows(fptr,rowNum,1,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"SCAN",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   ival = scan0+scan+1;  fits_write_col(fptr,TINT,colnum,rowNum+1,1,1,&ival,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CYCLE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   ival = j;  fits_write_col(fptr,TINT,colnum,rowNum+1,1,1,&ival,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"IF",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   ival = (ii-i0)+1;  fits_write_col(fptr,TINT,colnum,rowNum+1,1,1,&ival,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"BEAM",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   ival = 1;  fits_write_col(fptr,TINT,colnum,rowNum+1,1,1,&ival,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"DATE-OBS",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   strcpy(strArray[0],inFile->primary[ibeam].utc0);
		   fits_write_col(fptr,TSTRING,colnum,rowNum+1,1,1,strArray,&fitsStatus);
		
		   // TIME
		   strcpy(timeStr,inFile->primary[ibeam].utc0+11); // CHECK THAT THIS IS CORRECT *** GEORGE
		   sscanf(timeStr,"%f:%f:%f",&hr,&min,&sec);
		   fval = (hr*60*60+min*60.+sec);  
		   fits_get_colnum(fptr,CASEINSEN,(char *)"TIME",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"EXPOSURE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = intTime; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"OBJECT",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);

		   // Consider change source name
		   strcpy(srcName,inFile->beamHeader[ibeam].source);
		   for (kk=0;kk<nSrcChange;kk++)
		     {
		       if (strcmp(changeFrom[kk],srcName)==0)
			 {
			   printf("Changing source name: %s to %s\n",srcName,changeTo[kk]);
			   strcpy(srcName,changeTo[kk]);
			 }
		     }
		   strcpy(strArray[0],srcName);
		   fits_write_col(fptr,TSTRING,colnum,rowNum+1,1,1,strArray,&fitsStatus);
		   
		   // OBJ-RA
		   fits_get_colnum(fptr,CASEINSEN,(char *)"OBJ-RA",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].raDeg; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);		   
		
		   // OBJ-DEC
		   fits_get_colnum(fptr,CASEINSEN,(char *)"OBJ-DEC",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].decDeg; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // RESTFRQ
		   fits_get_colnum(fptr,CASEINSEN,(char *)"RESTFRQ",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 704e6; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   
		   // OBSMODE
		   strcpy(strArray[0],"POINTING");
		   fits_get_colnum(fptr,CASEINSEN,(char *)"OBSMODE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fits_write_col(fptr,TSTRING,colnum,rowNum+1,1,1,strArray,&fitsStatus);



		   // FREQRES
		   fits_get_colnum(fptr,CASEINSEN,(char *)"FREQRES",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 1e6*fabs(inFile->beam[ibeam].bandHeader[ii].f1-inFile->beam[ibeam].bandHeader[ii].f0)/inFile->beam[ibeam].bandHeader[ii].nchan; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // BANDWID
		   fits_get_colnum(fptr,CASEINSEN,(char *)"BANDWID",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 1e6*fabs(inFile->beam[ibeam].bandHeader[ii].f1-inFile->beam[ibeam].bandHeader[ii].f0); fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // CRPIX1 - must check if this is NCHAN - "array location of the reference pixel along axis 1"
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CRPIX1",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   //		   fval = (float)inFile->beam[ibeam].bandHeader[ii].nchan; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // CRVAL1
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CRVAL1",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 1e6*inFile->beam[ibeam].bandHeader[ii].f0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // CDELT1
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CDELT1",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 1e6*fabs(inFile->beam[ibeam].bandHeader[ii].f1-inFile->beam[ibeam].bandHeader[ii].f0)/inFile->beam[ibeam].bandHeader[ii].nchan; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // CRVAL3
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CRVAL3",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].raDeg; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // CRVAL4
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CRVAL4",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].decDeg; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // SCANRATE
		   // TSYS
		   // CALFCTR
		   // TCAL
		   // TCALTIME
		   // AZIMUTH
		   fits_get_colnum(fptr,CASEINSEN,(char *)"AZIMUTH",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].az; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   // ELEVATIO
		   fits_get_colnum(fptr,CASEINSEN,(char *)"ELEVATIO",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].el; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // PARANGLE
		   
		   fits_get_colnum(fptr,CASEINSEN,(char *)"PARANGLE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].paraAngle; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // FOCUSTAN
		   fits_get_colnum(fptr,CASEINSEN,(char *)"FOCUSTAN",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // FOCUSROT
		   fits_get_colnum(fptr,CASEINSEN,(char *)"FOCUSROT",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // TAMBIENT
		   fits_get_colnum(fptr,CASEINSEN,(char *)"TAMBIENT",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // PRESUURE
		   fits_get_colnum(fptr,CASEINSEN,(char *)"PRESSURE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // HUMIDITY
		   fits_get_colnum(fptr,CASEINSEN,(char *)"HUMIDITY",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   
		   // WINDSPEE
		   fits_get_colnum(fptr,CASEINSEN,(char *)"WINDSPEE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].windSpd; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // WINDDIRE
		   fits_get_colnum(fptr,CASEINSEN,(char *)"WINDDIRE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].windDir; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   
		   // Data		   		   
		   fits_get_colnum(fptr,CASEINSEN,(char *)"DATA",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,inFile->beam[ibeam].bandHeader[ii].nchan*pnum,fdata,&fitsStatus);
		 
		    
		   if (haveXpol == 1 && useXpol==1)
		     {
		       fits_get_colnum(fptr,CASEINSEN,(char *)"XPOLDATA",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		       fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,inFile->beam[ibeam].bandHeader[ii].nchan*2,xdata,&fitsStatus);		       
		     }
		   rowNum++;
		   if (oneScanPerFile==0)
		     scan++;  
		 }
	       free(fdata);
	       if (haveXpol == 1 && useXpol==1)
		 free(xdata);
	     }

	 }
       if (oneScanPerFile==1) scan=1;
       //       printf("Summing %d to scan0 %d\n",scan,scan0);
       scan0+=scan;
    	       
       if (i==nFiles-1 || join==0)
	 {
	   // Cross polarisation terms
	   if (useXpol==0) // Don't have or don't want Xpol informatino
	     {
	       fits_get_colnum(fptr,CASEINSEN,(char *)"XPOLDATA",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
	       fits_delete_col(fptr,colnum,&fitsStatus);
	       fits_report_error(stdout,fitsStatus);
	       fits_get_colnum(fptr,CASEINSEN,(char *)"XCALFCTR",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
	       fits_delete_col(fptr,colnum,&fitsStatus);
	       fits_report_error(stdout,fitsStatus);
	     }
	   
	   
	   
	   
	   
	   fits_report_error(stderr, fitsStatus);  /* print out any error messages */
	   fits_close_file(fptr,&fitsStatus);
	 }
       printf("Complete\n");
       
       sdhdf_closeFile(inFile);
    }
  free(inFile);
  free(strArray[0]);
  free(strArray);

}
