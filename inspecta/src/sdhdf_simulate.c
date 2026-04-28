//  Copyright (C) 2021, 2022 George Hobbs

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

/* To include

- Background sky from Haslam or equivalent
- HI from GASS or equivalent
- Able to include bright continuum sources (Haslam update or NVSS/SUMSS?)
*/
/*
# DON'T FORGET:
#
# BANDPASS - SLOPES
# RIPPLES
# GAIN VARIATIONS
# TONES AND COMBS
# PARALLACTIFY??


- Include bandpass shape/slope across band
- Correctly model gain in different pols and allow gain variations
- Include ripples
- Include RFI
      - impulsive (aircraft/satellites)
- Have ability to convert to counts (or Kelvin)
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "T2toolkit.h"
#include "hdf5.h"

#define MAX_INCLUDE_SPECTRA 128

typedef struct dumpParameterStruct {
  float tdump; // Time of spectral dump
  float raj0;   // Pointing position in degrees at centre of spectral dump (RA)
  float decj0;   // Pointing position in degrees at centre of spectral dump (DEC)
  
  double timeFromStart;
  double mjd;

  double gl;                  // Derived Galactic coordinate of beam position
  double gb;

} dumpParameterStruct;

typedef struct includeSpectrumStruct {
  int beam;
  int band;
  int dump;
  char fname[1024];
} includeSpectrumStruct;


typedef struct bandParameterStruct {
  char  label[MAX_STRLEN];    // Band label
  float f0;
  float f1;
  int   nchan;
} bandParameterStruct;

typedef struct beamParameterStruct {
  char  label[MAX_STRLEN];   // Beam label
  float tsys_aa;             // System temperature for this beam for AA polarisation (K)
  float tsys_bb;             // System temperature for this beam for BB polarisation (K)
  float delta_ra;            // Offset of beam from nominal pointing direction in RA (deg)
  float delta_dec;           // Offset of beam from nominal pointing direction in DEC (deg)
  char  src[MAX_STRLEN];     // Source name
} beamParameterStruct;

typedef struct parameterStruct {
  int nbeams;
  int nbands;
  int ndumps;
  int npols;
  char fname[MAX_STRLEN];
  double mjd0;
  char pid[MAX_STRLEN];
  char observer[MAX_STRLEN];
  char telescope[MAX_STRLEN];
  char rcvr[MAX_STRLEN];
  beamParameterStruct *beam;
  bandParameterStruct *band;
  dumpParameterStruct *dump;

  includeSpectrumStruct *includeSpectrum;
  int nIncludeSpectrum;
  
  int persistentRFI;
  int galacticHI;
} parameterStruct;


void convertGalactic(double raj,double decj,double *gl,double *gb);
double haversine(double centre_long,double centre_lat,double src_long,double src_lat);
float get_tsky(float *tsky_gl,float *tsky_gb,float *tsky_kelvin,int nTsky,float beamGl_deg,float beamGb_deg);
double beamScaling(double angle,float diameter,float freq);


void help()
{
  printf("sdhdf_simulate: Simulation code for sdhdf files\n");
  printf("-h           This help\n");
  printf("-p <file>    Input parameter file\n");

}

int main(int argc,char *argv[])
{
  int i,j,k,ii,jj;
  double beamRA_deg,beamDec_deg,beamGl_deg,beamGb_deg;
  long iseed = TKsetSeed();
  char fname[MAX_STRLEN]="unset";
  char paramFile[MAX_STRLEN]="unset";
  char labelStr[MAX_STRLEN];
  sdhdf_fileStruct          *outFile;
  sdhdf_beamHeaderStruct    *beamHeader;
  sdhdf_primaryHeaderStruct *primaryHeader;
  sdhdf_bandHeaderStruct    *bandHeader;
  sdhdf_obsParamsStruct     *obsParams;
  double dt;
  int nbeam=0;
  int ibeam=0;
  int ndump=0;
  int idump,iband;
  int nband;
  int nchan,nsub,npol;
  int nDataAttributes=0;
  int nFreqAttributes=0;
  sdhdf_attributes_struct *freqAttributes;
  sdhdf_attributes_struct *dataAttributes;
  double *freq;
  float *data;
  float ssys_aa = 22;
  float ssys_bb = 23;
  float *tsky_gl;
  float *tsky_gb;
  float *tsky_kelvin;
  int   nTsky =  3147775;
  float chbw;
  float req_chbw = 0.5e3; // in Hz
  float band_bw  = 192; // in MHz
  float rcvr_f0  = 700.0; // MHz
  float rcvr_f1  = 1852.0; // MHz. How will the band be split ??
  float signal;
  FILE *fin;
  char line[MAX_STRLEN];
  char word1[MAX_STRLEN];
  char word2[MAX_STRLEN];
  parameterStruct *params;
  double d_raj,d_decj;
  float tsky_model,tsky_model0;
  float ssky_model;
  float noiseScaleAA,noiseScaleBB;
  float ssrc_model;
  float *spectrumAdd;
  int addUserSpectrum=0;
  double cmb = 2.725;
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-p")==0)
	strcpy(paramFile,argv[++i]);
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
    }

  if (strcmp(paramFile,"unset")==0)
    {
      printf("ERROR: must set a parameter filename using the -p option\n");
      exit(1);
    }
  
  if (!(fin = fopen(paramFile,"r")))
    {
      printf("ERROR: unable to open parameter file >%s<\n",paramFile);
      exit(1);
    }
  params = (parameterStruct *)malloc(sizeof(parameterStruct));
  params->includeSpectrum = (includeSpectrumStruct *)malloc(sizeof(includeSpectrumStruct)*MAX_INCLUDE_SPECTRA);
  params->nIncludeSpectrum=0;
   
  // initialise
  params->persistentRFI=0;
  params->galacticHI=0;
  
  while (!feof(fin))
    {
      if (fgets(line,MAX_STRLEN,fin)!=NULL)
	{
	  if (line[0]!='#' && strlen(line) > 1) // Not a comment line
	    {
	      strcpy(word2,"NULL");
	      if (sscanf(line,"%s %s",word1,word2)==2)
		{
		  if (strcmp(word1,"nbeams:")==0)
		    {
		      sscanf(word2,"%d",&(params->nbeams));
		      // Allocate memory for this number of beams
		      params->beam = (beamParameterStruct *)malloc(sizeof(beamParameterStruct)*params->nbeams);
		      ibeam=0;
		    }
		  else if (strcmp(word1,"persistent_rfi:")==0 && strcmp(word2,"on")==0)
		    params->persistentRFI=1;
		  else if (strcmp(word1,"galactic_hi:")==0 && strcmp(word2,"gaussian")==0)
		    params->galacticHI=1;
		  else if (strcmp(word1,"inputSpectrum:")==0)
		    {
		      int nv = params->nIncludeSpectrum;
		      sscanf(line,"%s %d %d %d %s",word1,&(params->includeSpectrum[nv].beam),
			     &(params->includeSpectrum[nv].band),
			     &(params->includeSpectrum[nv].dump),
			     params->includeSpectrum[nv].fname);
		      (params->nIncludeSpectrum)++;
		    }
		  else if (strcmp(word1,"beam:")==0)
		    {
		      if (ibeam == params->nbeams)
			printf("WARNING: TOO MANY BEAMS SPECIFIED. IGNORING >%s<\n",line);
		      else
			{
			  sscanf(line,"%s %s %s %f %f %f %f",word1,params->beam[ibeam].label,params->beam[ibeam].src,
				 &(params->beam[ibeam].tsys_aa),
				 &(params->beam[ibeam].tsys_bb),&(params->beam[ibeam].delta_ra),&(params->beam[ibeam].delta_dec));			  
			  ibeam++;
			}
		    }
		  else if (strcmp(word1,"nbands:")==0)
		    {
		      sscanf(word2,"%d",&(params->nbands));
		      params->band = (bandParameterStruct *)malloc(sizeof(bandParameterStruct)*params->nbands);
		      iband=0;
		    }
		  else if (strcmp(word1,"band:")==0)
		    {
		      if (iband == params->nbands)
			{
			  printf("WARNING: TOO MANY BANDS SPECIFIED. IGNORING >%s<\n",line);
			  printf("word1 = %s, word2 = %s\n",word1,word2);
			  printf("iband = %d, nbands = %d\n",iband,params->nbands);
			}
			  else
			{
			  sscanf(line,"%s %s %f %f %d",word1,params->band[iband].label,
				 &(params->band[iband].f0), &(params->band[iband].f1),&(params->band[iband].nchan));
			  iband++;
			}
		    }
		  else if (strcmp(word1,"ndumps:")==0)
		    {
		      sscanf(word2,"%d",&(params->ndumps));
		      params->dump = (dumpParameterStruct *)malloc(sizeof(dumpParameterStruct)*params->ndumps);
		      idump=0;
		    }
		  else if (strcmp(word1,"dump:")==0)
		    {
		      if (idump == params->ndumps)
			  printf("WARNING: TOO MANY DUMPS SPECIFIED. IGNORING >%s<\n",line);
		      else
			{
			  sscanf(line,"%s %f %f %f",word1,
				 &(params->dump[idump].tdump), &(params->dump[idump].raj0),&(params->dump[idump].decj0));
			  idump++;
			}
		    }
		  else if (strcmp(word1,"npols:")==0)
		    sscanf(word2,"%d",&(params->npols));
		  else if (strcmp(word1,"start_mjd:")==0)
		    sscanf(word2,"%lf",&(params->mjd0));
		  else if (strcmp(word1,"output:")==0)
		    strcpy(params->fname,word2);
		  else if (strcmp(word1,"pid:")==0)
		    strcpy(params->pid,word2);
		  else if (strcmp(word1,"observer:")==0)
		    strcpy(params->observer,line+10);
		  else if (strcmp(word1,"receiver:")==0)
		    strcpy(params->rcvr,word2);
		  else if (strcmp(word1,"telescope:")==0)
		    strcpy(params->telescope,word2);
		}
	    }
	}
    }
  fclose(fin);

  // Load Tsky if required
  tsky_gl = (float *)malloc(sizeof(float)*nTsky);
  tsky_gb = (float *)malloc(sizeof(float)*nTsky);
  tsky_kelvin = (float *)malloc(sizeof(float)*nTsky);

  nTsky = 0;
  if (!(fin = fopen("haslam408_ds_Remazeilles2014.allSky.dat","r")))
    {
      printf("Unable to open haslam408_ds_Remazeilles2014.allSky.dat for the sky temperature model\n");      
    }
  else
    {
      while (!feof(fin))
	{
	  if (fscanf(fin,"%f %f %f",&tsky_gl[nTsky],&tsky_gb[nTsky],&tsky_kelvin[nTsky])==3)
	    nTsky++;
	}
      fclose(fin);
      printf("Loaded Tsky model\n");

    }

  // Calculate times and positions
  dt=0;
  for (i=0;i<params->ndumps;i++)
    {
      dt+=params->dump[i].tdump/2.;
      params->dump[i].timeFromStart = dt;
      params->dump[i].mjd = params->mjd0 + dt/86400.0;
      d_raj = params->dump[i].raj0;
      d_decj = params->dump[i].decj0;
      convertGalactic(d_raj*M_PI/180.,d_decj*M_PI/180.0,&(params->dump[i].gl),&(params->dump[i].gb));
      dt+=params->dump[i].tdump/2.;
    }

  
  
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
  sdhdf_initialiseFile(outFile);
  if (sdhdf_openFile(params->fname,outFile,3)==-1)
    {
      printf("Unable to open output file >%s<\n",fname);
      free(outFile);
      exit(1);
    }

  nbeam = params->nbeams;
  nband = params->nbands; 
  printf("Simulating %d bands\n",nband);
  
  beamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nbeam);
  primaryHeader = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct));

  for (i=0;i<nbeam;i++)
    {
      strcpy(beamHeader[i].label,params->beam[i].label);
      beamHeader[i].nBand = nband;
      strcpy(beamHeader[i].source,params->beam[i].src); 
    }
  sdhdf_writeBeamHeader(outFile,beamHeader,nbeam,(char *)"4.0");

  // Write the primary header
  strcpy(primaryHeader[0].date,"unknown");
  strcpy(primaryHeader[0].hdr_defn,"5");
  strcpy(primaryHeader[0].hdr_defn_version,"1.9");
  strcpy(primaryHeader[0].file_format,"unknown");
  strcpy(primaryHeader[0].file_format_version,"unknown");
  primaryHeader[0].sched_block_id = 1;
  strcpy(primaryHeader[0].cal_mode,"OFF");
  strcpy(primaryHeader[0].instrument,"simulate");
  strcpy(primaryHeader[0].observer,params->observer);
  strcpy(primaryHeader[0].pid,params->pid);
  strcpy(primaryHeader[0].rcvr,params->rcvr);
  strcpy(primaryHeader[0].telescope,params->telescope);
  strcpy(primaryHeader[0].utc0,"unknown");
  primaryHeader[0].nbeam = nbeam;
  sdhdf_writePrimaryHeader(outFile,primaryHeader);


  // Simulate data for each beam
  for (i=0;i<nbeam;i++)
    {
      printf("Processing beam %d out of %d\n",i+1,nbeam);
      bandHeader = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nband); // FIX ME -- may have different numbers of bands 
      for (j=0;j<nband;j++)
	{
	  printf("Processing band %d out of %d\n",j+1,nband);
	  ndump = params->ndumps;
	  npol  = params->npols;
	  
	  bandHeader[j].f0 = params->band[j].f0;
	  bandHeader[j].f1 = params->band[j].f1;
	  bandHeader[j].fc = (bandHeader[j].f0+bandHeader[j].f1)/2.0;

	  nchan = params->band[j].nchan;
	  spectrumAdd = (float *)malloc(sizeof(float)*nchan);
	  addUserSpectrum=0;
	  obsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
	  freq = (double *)malloc(sizeof(double)*nchan);
	  data = (float *)malloc(sizeof(float)*nchan*ndump*npol);
	  strcpy(bandHeader[j].label,params->band[j].label);

	  bandHeader[j].nchan = nchan;
	  bandHeader[j].npol = npol;
	  if (npol==4)
	    strcpy(bandHeader[j].pol_type,"AABBCRCI");
	  else if (npol==2)
	    strcpy(bandHeader[j].pol_type,"AA,BB");
	  else if (npol==1)
	    strcpy(bandHeader[j].pol_type,"AA+BB");
	  bandHeader[j].dtime = params->dump[0].tdump; 
	  bandHeader[j].ndump = ndump;	 
	  chbw = (bandHeader[j].f1 - bandHeader[j].f0)/(float)nchan;
	  // Setup the observation parameters for each band
	  for (k=0;k<nchan;k++)
	    freq[k] = bandHeader[j].f0+k*(bandHeader[j].f1-bandHeader[j].f0)/(float)nchan + chbw/2.0;
	  
	  for (k=0;k<ndump;k++)
	    {
	      printf("Processing spectral dump %d/%d\n",k+1,ndump);
	      beamRA_deg = params->dump[k].raj0; // SHOULD INCLUDE BEAM OFFSET ... FIX ME
	      beamDec_deg = params->dump[k].decj0; // SHOULD INCLUDE BEAM OFFSET ... FIX ME
	      beamGl_deg = params->dump[k].gl;     // SHOULD INCLUDE BEAM OFFSET ... FIX ME
	      beamGl_deg = params->dump[k].gb;    // SHOULD INCLUDE BEAM OFFSET ... FIX ME
	      
	      obsParams[k].timeElapsed = params->dump[k].timeFromStart;
	      strcpy(obsParams[k].timedb,"unknown");
	      obsParams[k].mjd = params->dump[k].mjd;
	      obsParams[k].fractional_mjd = 0;
	      strcpy(obsParams[k].utc,"unknown");
	      strcpy(obsParams[k].ut_date,"unknown");
	      strcpy(obsParams[k].local_time,"unknown");
	      strcpy(obsParams[k].raStr,"unknown");
	      strcpy(obsParams[k].decStr,"unknown");
	      obsParams[k].raOffset = 0;
	      obsParams[k].decOffset = 0;

	      obsParams[k].raDeg  = beamRA_deg;
	      obsParams[k].decDeg = beamDec_deg;
	      obsParams[k].gl     = beamGl_deg;
	      obsParams[k].gb     = beamGb_deg;
	      obsParams[k].az = 0;
	      obsParams[k].el = 0;
	      obsParams[k].az_drive_rate = 0;
	      obsParams[k].ze_drive_rate = 0;
	      obsParams[k].hourAngle = 0;
	      obsParams[k].paraAngle = 0;
	      obsParams[k].windDir = 0;
	      obsParams[k].windSpd = 0;		      

	      // Add in defined spectra
	      for (jj=0;jj<params->nIncludeSpectrum;jj++)
		{
		  if (params->includeSpectrum[jj].beam == i && params->includeSpectrum[jj].band == j && params->includeSpectrum[jj].dump == k)
		    {
		      fin = fopen(params->includeSpectrum[jj].fname,"r");
		      for (ii=0;ii<nchan;ii++)
			fscanf(fin,"%g",&spectrumAdd[ii]);
		      fclose(fin);
		      addUserSpectrum=1;
		    }
		}

	    
	      for (ii=0;ii<nchan;ii++)
		{
		  signal = 0;
		  if (params->galacticHI==1)
		    signal = 4*exp(-pow(freq[ii]-1421.0,2)/2./0.2/0.2);
		  if (params->persistentRFI==1)
		    {
		      if (freq[ii] > 754 && freq[ii] < 768)
			signal += 1e6;
		      else if (freq[ii] > 768 && freq[ii] < 788)
			signal += 8e5;
		      else if (freq[ii] > 869.95 && freq[ii] < 875.05)
			signal += 4e5;
		      else if (freq[ii] > 875.05 && freq[ii] < 889.95)
			signal += 3e5;
		      else if (freq[ii] > 943.4 && freq[ii] < 951.8)
			signal += 3e5;
		      else if (freq[ii] > 953.7 && freq[ii] < 960)
			signal += 3e5;
		      else if (freq[ii] > 1017 && freq[ii] < 1019)
			signal += 3e5;
		      else if (freq[ii] > 1023 && freq[ii] < 1025)
			signal += 9e5;
		      else if (freq[ii] > 1029 && freq[ii] < 1031)
			signal += 8e5;
		      else if (freq[ii] > 1805 && freq[ii] < 1865)
			signal += 9e5;
		    }
		  if (addUserSpectrum==1)
		    signal += spectrumAdd[ii]; 
		  // Build up ssys_aa and ssys_bb
		  tsky_model = 0;
		  if (ii==0)
		    {
		      tsky_model0 = get_tsky(tsky_gl,tsky_gb,tsky_kelvin,nTsky,beamGl_deg,beamGb_deg);
		      //		      printf("Value = %g\n",tsky_model0);
		    }
		  // Add in source
		  {
		    double angle;
		    //  62.0849118333333
		    // -65.7525223888889
		    angle = haversine(62.0849118333333,-65.7525223888889,beamRA_deg,beamDec_deg);
		    ssrc_model = 10*beamScaling(angle,64,freq[ii]);
		    if (ii==0)
		      printf("%d %g %g\n",k,angle,ssrc_model);
		  }
		  tsky_model = (tsky_model0 - cmb) * pow(freq[ii]/408.0,-2.6) + cmb;
		  ssky_model = tsky_model ; // FIX ME -- WITH GAINS ETC.
		
		  //		  printf("chbw = %g, dtime = %g ssky_model = %g\n",chbw,bandHeader[j].dtime,ssky_model);
		  noiseScaleAA = params->beam[i].tsys_aa/sqrt(chbw*1e6*bandHeader[j].dtime);
		  noiseScaleBB = params->beam[i].tsys_bb/sqrt(chbw*1e6*bandHeader[j].dtime);
		  //		  printf("scaleAA, scaleBB = %g %g\n",noiseScaleAA,noiseScaleBB);

		    
		  data[k*nchan*npol + ii]           = noiseScaleAA*TKgaussDev(&iseed) + ssys_aa + ssky_model + ssrc_model + signal;
		  data[k*nchan*npol + nchan + ii]   = noiseScaleBB*TKgaussDev(&iseed) + ssys_bb + ssky_model + ssrc_model + signal;
		  data[k*nchan*npol + 2*nchan + ii] = sqrt(noiseScaleAA*noiseScaleBB)*TKgaussDev(&iseed);
		  data[k*nchan*npol + 3*nchan + ii] = sqrt(noiseScaleAA*noiseScaleBB)*TKgaussDev(&iseed);
	    }
	}
	  sdhdf_writeObsParams(outFile,bandHeader[j].label,beamHeader[i].label,j,obsParams,ndump,1,(char *)"4.0");			       
      free(obsParams);
      free(spectrumAdd);
	  // SHOULD SET UP ATTRIBUTES
	  // FIX ME: Sending only one frequency channel through
      sdhdf_writeSpectrumData(outFile,beamHeader[i].label,bandHeader[j].label,i,j,data,freq,1,nchan,1,npol,ndump,1,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	  free(freq);
	  free(data);

	}
      sdhdf_writeBandHeader(outFile,bandHeader,beamHeader[i].label,nband,1,(char *)"4.0");
    
      free(bandHeader);

}


sdhdf_closeFile(outFile);
  free(outFile);
  free(primaryHeader);
  free(beamHeader);
  free(params->beam);
  free(params->band);
  free(params->dump);
  free(params->includeSpectrum);
  free(params);
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


// Should do a direct determination of the place in the file to look for -- FIX ME --- the following is slow
float get_tsky(float *tsky_gl,float *tsky_gb,float *tsky_kelvin,int nTsky,float beamGl_deg,float beamGb_deg)
{
  int i;
  double dist;

  for (i=0;i<nTsky;i++)
    {
      dist = haversine(tsky_gl[i],tsky_gb[i],beamGl_deg,beamGb_deg);
      if (dist < 0.1)
	{
	  return tsky_kelvin[i];
	  break;
	}
    }
  return -1;
}

double haversine(double centre_long,double centre_lat,double src_long,double src_lat)
{
  double dlon,dlat,a,c;
  double deg2rad = M_PI/180.0;

  centre_long*=deg2rad;
  centre_lat*=deg2rad;
  src_long*=deg2rad;
  src_lat*=deg2rad;
  
  /* Apply the Haversine formula */
  dlon = (src_long - centre_long);
  dlat = (src_lat  - centre_lat);
  a = pow(sin(dlat/2.0),2) + cos(centre_lat) *
    cos(src_lat)*pow(sin(dlon/2.0),2);
  if (a==1)
    c = M_PI;
  else
    c = 2.0 * atan2(sqrt(a),sqrt(1.0-a));
  return c/deg2rad;
}

// Calculates the reduction in gain caused by the source being offset from the beam pointing direction
// The input angle is in degrees, the diameter is the telescope diameter (m) and freq is the observing frequency (MHz).  The output value is the scaling factor
//

double beamScaling(double angle,float diameter,float freq)
{
  double radangle = angle*M_PI/180.0;
  double tt;
  double lambda = 3.0e8/(freq*1.0e6);
  double ang;
  ang = 1.22*lambda/diameter;
  
  tt = radangle*M_PI/ang;
  if (tt == 0)
    return 1;
  else
    return pow(sin(tt)/tt,2);
}
