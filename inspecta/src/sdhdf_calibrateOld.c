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


// To fix:
//
// Should check if the user actually wants polarisation and flux calibration, or just 1
// Need a way to read the calibration information from an earlier file - or inject an earlier cal into a different file
// Need a way to specify which fluxcal/pcm file to use (and how the user can know what options exist)
// Currrently hardcode alpha = -45 degrees in sdhdfProc_calibration.c -- FIX THIS!
// Currently always time-averaging the noise source signal
// FIX ME: IMPORTANT: ISSUE IN THAT SHOULDN'T INTERPOLATE USING BAD ENTIRES (AS THEY ARE ZERO) (important for noise source stokes)
// Must work with older-style, non-normalised files.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"

#define VERSION "v0.5"
#define MAX_POL_CAL_CHAN 4096    // FIX ME -- SHOULD SET DYNAMICALLY

void help()
{
  printf("sdhdf_calibrate %s (INSPECTA %s)\n",VERSION,SOFTWARE_VER);
  printf("Authors: G. Hobbs\n");
  printf("Purpose: to calibrate the data sets\n");
  printf("\n");
  printf("Command line arguments:\n\n");
  printf("-e <ext>    Output file extension\n");
  printf("-h          This help\n");
  printf("\n\n");
  printf("Example: sdhdf_calibrate -e cal uwl_*.hdf.T.f1024\n");
  exit(1);
}


int main(int argc,char *argv[])
{
  int b,i,j,k,ii,kk,beam,band;
  int nFiles=0;
  char fluxCalFile[1024] = "uwl_220705_132107.fluxcal";
  char fname[MAX_FILES][MAX_STRLEN];
  char runtimeDir[MAX_STRLEN];
  float pol1,pol2,pol3,pol4;
  int tcal = 0;  // Scal if tcal = 0.
  
  float *out_data,*out_freq;

  sdhdf_fileStruct *inFile,*outFile;
  sdhdf_calibration *polCal;
  sdhdf_fluxCalibration *fluxCal;
  int nPolCalChan=0;
  int nFluxCalChan=0;
  FILE *fout;

  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;

  sdhdf_bandHeaderStruct *inBandParams;
  sdhdf_obsParamsStruct  *outObsParams;
  float *dataWts;
  
  double stokesCalMeasured[4];
  double actualNoiseStokes[4];
  double coherencyCalMeasuredOn[4];
  double coherencyCalMeasuredOff[4];
  double complex rho[2][2];
  double complex Jast[2][2];
  double complex Rdag[2][2];
  double complex R_feed_pa[2][2];
  double complex R_feed_pa_dag[2][2];
  double complex finalJ[2][2];
  double aa,bb,rab,iab;
  double final_aa,final_bb,final_rab,final_iab;
  double fluxScale;
  double paraAng;
  double alpha = -45.0*M_PI/180.0; // FIX ME
  double freq;
  int nc=0;
  
  int npol,nchan,ndump,nchanCal,ndumpCal;
  int ichan;

  int normCal=0;    // Normalise the noise source counts
  int normAstro=0;  // Normalise the astronomy source counts
  double scaleFactor;
  double tdumpAstro=-1;
  int setTdumpAstro=0;
  double tdumpCal=-1;
  int nchanAstro=-1;
  int setNchanAstro=0;
  int nchanFreq;
  
  char extension[1024];
  char oname[1024];
  char args[MAX_STRLEN]="";
  float averageCalValuesF0=-1;
  float averageCalValuesF1=-1;
  int averageCal=0;
  int noFluxCal=0;
  char pcmFile[1024];
  int pcmClear=0;
  
  // Default PCM file // FIX ME - THIS IS FOR PARKES
  strcpy(pcmFile,"uwl_210621_145853.extzoom.cal.hdf");

  
  // Setup output defaults
  strcpy(oname,"sdhdf_calibrate_output.hdf");
  strcpy(extension,"calibrate");
  
  // Load in the file names
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-e")==0)
	{sdhdf_add2arg(args,argv[i],argv[i+1]); strcpy(extension,argv[++i]);}
      else if (strcmp(argv[i],"-tcal")==0)
	tcal=1;
      else if (strcmp(argv[i],"-averageCal")==0)
	{
	  averageCal=1;
	  sscanf(argv[++i],"%f",&averageCalValuesF0);
	  sscanf(argv[++i],"%f",&averageCalValuesF1);
	}
      else if (strcasecmp(argv[i],"-averageCalByBand")==0)
	averageCal=2;
      else if (strcmp(argv[i],"-norm")==0)
	{
	  normCal=1;
	  normAstro=1;
	}
      else if (strcasecmp(argv[i],"-noFluxCal")==0)
	noFluxCal=1;
      else if (strcmp(argv[i],"-tdumpAstro")==0)
	{sscanf(argv[++i],"%lf",&tdumpAstro); setTdumpAstro=1;}
      else if (strcmp(argv[i],"-nchanAstro")==0)
	{sscanf(argv[++i],"%d",&nchanAstro); setNchanAstro=1;}
      else if (strcmp(argv[i],"-pcmClear")==0)
	pcmClear=1;
      else if (strcmp(argv[i],"-fluxcal")==0)
	strcpy(fluxCalFile,argv[++i]);
      else if (strcasecmp(argv[i],"-pcmfile")==0)
	strcpy(pcmFile,argv[++i]);
      else if (strcmp(argv[i],"-h")==0)
	help();
      else
	strcpy(fname[nFiles++],argv[i]);
    }
  
  inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));
  polCal = (sdhdf_calibration *)malloc(sizeof(sdhdf_calibration)*MAX_POL_CAL_CHAN);
  fluxCal = (sdhdf_fluxCalibration *)malloc(sizeof(sdhdf_fluxCalibration)*MAX_POL_CAL_CHAN); //  Should change to MAX_FLUXCAL OR SIMILAR ** FIX ME
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }


  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdf_calibrate requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      if (sdhdf_openFile(fname[i],inFile,1)==-1)
	printf("Warning: unable to open file >%s<. Skipping\n",fname[i]);
      else
	{
	  sdhdf_formOutputFilename(fname[i],extension,oname);
	  printf("Opening file >%s<\n",oname);
	  sdhdf_openFile(oname,outFile,3);

	  sdhdf_loadMetaData(inFile);
	  sdhdf_allocateBeamMemory(outFile,inFile->nBeam);
	  
	  // Now process the astronomy data
	  for (b=0;b<inFile->nBeam;b++)
	    {
	      inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);      
	      sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);
	      
	      for (j=0;j<inFile->beam[b].nBand;j++)
		{
		  printf("Processing subband %d\n",j);
		  sdhdf_loadBandData(inFile,b,j,1);
		  nchan  = inFile->beam[b].bandHeader[j].nchan;
		  npol   = inFile->beam[b].bandHeader[j].npol;
		  ndump  = inFile->beam[b].bandHeader[j].ndump;
		  nchanFreq = inFile->beam[b].bandData[j].astro_data.nFreqDumps;

		  sdhdf_copyAttributes(inFile->beam[b].bandData[j].astro_obsHeaderAttr,inFile->beam[b].bandData[j].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
		  sdhdf_copyAttributes(inFile->beam[b].bandData[j].astro_obsHeaderAttr_freq,inFile->beam[b].bandData[j].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);
		  
		  
		  // Copy observation parameters
		  outObsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
		  for (kk=0;kk<ndump;kk++)
		    sdhdf_copySingleObsParams(inFile,b,j,kk,&outObsParams[kk]);
		  
		  out_freq  = (float *)malloc(sizeof(float)*nchan*nchanFreq);
		  out_data  = (float *)calloc(sizeof(float),nchan*npol*ndump);
		  dataWts   = (float *)calloc(sizeof(float),nchan*ndump);

		  if (setTdumpAstro==0)
		    tdumpAstro = inFile->beam[b].bandHeader[j].dtime;
		  if (setNchanAstro==0)
		    nchanAstro = nchan;

		  // Load in the calibration data for this band
		  	  // FIX ME: Should check if the user wants polarisation calibration	

		  if (averageCal==2)
		    {
		      averageCalValuesF0 = inFile->beam[b].bandHeader[j].f0;
		      averageCalValuesF1 = inFile->beam[b].bandHeader[j].f1;
		    }
		  
		  // Load the information within the PCM file	  
		  // This is assuming same number of channels in the PCM cal as in the noise source *** FIX ME
		  //
		  sdhdf_loadPCM(polCal,&nPolCalChan,"parkes","UWL",pcmFile,averageCal,averageCalValuesF0,averageCalValuesF1,pcmClear); // REMOVE HARDCODE
		  sdhdf_formPCM_response(polCal,nPolCalChan);
		  
		  
		  // Load the Scal information
		  nFluxCalChan=0;
		  if (tcal==0)
		    sdhdf_loadFluxCal(fluxCal,&nFluxCalChan,"parkes","UWL",fluxCalFile,0.0); // REMOVE HARDCODE
		  else
		    sdhdf_loadTcal(fluxCal,&nFluxCalChan,"parkes","UWL","tcal_noflag.dat");  // Note this is a text file -- FIX ME -- should make consistent with fluxcal
		  
		  
		  //	  for (j=0;j<nPolCalChan;j++)
		  //	    printf("%d %g %g %g %g %g %g %g %g %g %g %g\n",
		  //		   j,polCal[j].freq,polCal[j].noiseSource_QoverI,polCal[j].noiseSource_UoverI,polCal[j].noiseSource_VoverI,
		  //		   polCal[j].constant_gain,polCal[j].constant_diff_gain,polCal[j].constant_diff_phase,
		  //		   polCal[j].constant_b1,polCal[j].constant_b2,polCal[j].constant_r1,polCal[j].constant_r2);
		  
		  
		  // Obtain NOISE ON - NOISE OFF
		  // Should be able to do this from a different file ** FIX ME	  
		  
		  // Request CAL_ON-CAL_OFF for specfic frequency covered by the PCM file
		  beam=0; // FIX ME
		  sdhdf_set_stokes_noise_measured(inFile,beam,polCal,nPolCalChan,normCal,averageCal,averageCalValuesF0,averageCalValuesF1);
		  sdhdf_calculate_gain_diffgain_diffphase(polCal,nPolCalChan);
		  sdhdf_calculate_timedependent_response(polCal,nPolCalChan);   
		  
		  /*
		    for (j=0;j<nPolCalChan;j++)
		    {
		    printf("%g Have (%g,%g,%g,%g) and (%g,%g,%g,%g) gain = %g diffgain = %g diffphase = %g\n",polCal[j].freq,polCal[j].stokes_noise_measured[0],polCal[j].stokes_noise_measured[1],
		    polCal[j].stokes_noise_measured[2],polCal[j].stokes_noise_measured[3],
		    polCal[j].stokes_noise_actual[0],polCal[j].stokes_noise_actual[1],polCal[j].stokes_noise_actual[2],polCal[j].stokes_noise_actual[3],
		     polCal[j].gain,polCal[j].diff_gain,polCal[j].diff_phase);
		     
	    }
	  */

		  for (kk=0;kk<nPolCalChan;kk++)
		    printf("td_response = %g %g %g %g %g %g %g %g %g %g\n",polCal[kk].freq,polCal[kk].td_response[0][0],polCal[kk].td_response[0][1],polCal[kk].td_response[0][2],polCal[kk].td_response[1][0],polCal[kk].td_response[1][1],polCal[kk].td_response[1][2],polCal[kk].td_response[2][0],polCal[kk].td_response[2][1],polCal[kk].td_response[2][2]);

		  
		  for (k=0;k<inFile->beam[b].bandHeader[j].ndump;k++)
		    {
		      printf("Processing spectral dump %d\n",k);
		      //		      printf("Processing: %d %d %d\n",b,j,k);
		      for (ii=0;ii<nchan;ii++)
			{
			  freq = inFile->beam[b].bandData[j].astro_data.freq[k*nchan+ii];  // FIX ME -- FREQ AXIS FOR DUMP
			
			  out_freq[k*nchan+ii] = freq;  // NOTE: Assuming that the frequency is constant in all bands - not true if Doppler corrected -- FIX ME
			  dataWts[ii+k*nchan]  = inFile->beam[b].bandData[j].astro_data.dataWeights[k*nchan+ii];

			  // FIX ME: SHOULD ACCOUNT FOR WEIGHTING

			  if (normAstro==1)
			    {
			      //			      printf("Normalising using tdumpAstro = %g and nchanAstro = %d\n",tdumpAstro,nchanAstro);
			      //			      scaleFactor = 1.0/((double)nchanAstro/(double)tdumpAstro); 
  
			    }
			  else
			    scaleFactor=1;
			  scaleFactor=1;

			  aa  = scaleFactor*inFile->beam[b].bandData[j].astro_data.pol1[ii+k*nchan];
			  bb  = scaleFactor*inFile->beam[b].bandData[j].astro_data.pol2[ii+k*nchan];
			  rab = scaleFactor*inFile->beam[b].bandData[j].astro_data.pol3[ii+k*nchan];
			  iab = scaleFactor*inFile->beam[b].bandData[j].astro_data.pol4[ii+k*nchan];

	 		  sdhdf_complex_matrix_2x2(rho,aa,rab-I*iab,rab+I*iab,bb);

			  // Should obtain the relevant channel in the PCM structure
			  ichan = (int)((freq-polCal[0].freq)/(polCal[1].freq-polCal[0].freq)+0.5); // FIX ME ******
			  //			  ichan = 29; // FIX ME HARDCODE ******
			  // NOT SURE THIS IS TD_RESPNSE -- SHOULD BE R_PCM --- CHECK -- FIX ME
			  //			  sdhdf_copy_complex_matrix_2x2(Jast,polCal[ichan].td_response); // NEEDS FIXING WITH ICHAN

			  
			  // printf("%d td_response = ",ichan); sdhdf_display_complex_matrix_2x2(polCal[ichan].td_response);

			  sdhdf_copy_complex_matrix_2x2(Jast,polCal[ichan].td_response); 
			  sdhdf_multiply_complex_matrix_2x2(Jast,rho);

			  //  printf("%d rho = ",ii);  sdhdf_display_complex_matrix_2x2(rho);


			  // **** At this point we should have stabilised by the noise source ****
			  
			  // Form conjugate
			  sdhdf_complex_matrix_2x2_dagger(polCal[ichan].td_response,Rdag);
			  sdhdf_multiply_complex_matrix_2x2(Jast,Rdag);
			  // printf("%d Jast = ",ii);  sdhdf_display_complex_matrix_2x2(Jast);

			  // Obtain parallactic angle in radians
			  paraAng = inFile->beam[b].bandData[j].astro_obsHeader[k].paraAngle * M_PI/180.0;

 			  sdhdf_complex_matrix_2x2(R_feed_pa,
						   I*sin(alpha+paraAng),
						   -I*cos(alpha+paraAng),
						   -I*cos(alpha+paraAng),
						   -I*sin(alpha+paraAng));
			  // printf("%d Feed PA = ",ii); sdhdf_display_complex_matrix_2x2(R_feed_pa);
			  sdhdf_complex_matrix_2x2_dagger(R_feed_pa,R_feed_pa_dag);
			  // printf("%d Feed PA_dag = ",ii);  sdhdf_display_complex_matrix_2x2(R_feed_pa_dag);


   			  sdhdf_copy_complex_matrix_2x2(finalJ,R_feed_pa);
			  sdhdf_multiply_complex_matrix_2x2(finalJ,Jast);
			  // printf("%d Rfeed x Jast = ",ii);  sdhdf_display_complex_matrix_2x2(finalJ);

			  sdhdf_multiply_complex_matrix_2x2(finalJ,R_feed_pa_dag);
			  //			  printf("%d finalJ = ",ii); sdhdf_display_complex_matrix_2x2(finalJ);


			
			  // Multiply by Scal
			  // FIX ME -- NEED TO INTERPLATE -- USE DIFFERENT ICHAN (IN CASE NCHAN IS DIFFERENT WITH FLUXCAL
			  // printf("%d finalJ ",ii);  sdhdf_display_complex_matrix_2x2(finalJ);
			  if ((averageCal==1 || averageCal==2) && ii == 0)
			    {
			      int c,nv=0;
			      fluxScale=0;
			      for (c=0;c<nFluxCalChan;c++)
				{
				  if (fluxCal[c].freq >= averageCalValuesF0 && fluxCal[c].freq <= averageCalValuesF1)
				    {
				      fluxScale += (fluxCal[c].scalAA + fluxCal[c].scalBB);
				      nv++;
				    }
				}
			      fluxScale/=(double)nv;
			      printf("Using averaged flux scale of %g\n",fluxScale);
			    }
			  else if (averageCal==0)
			    fluxScale = fluxCal[ichan].scalAA + fluxCal[ichan].scalBB;

			  if (noFluxCal == 1)
			    fluxScale=1;
			  
			  //			  printf("flux scale = %g\n",fluxScale);
			  
			  if ((normAstro == 1 || normCal==1) &&
			      (((averageCal==1 || averageCal==2) && ii == 0) || averageCal == 0))// Currently normalising both
			    {
			      int nbinCal = 32; // WARNING HARDCODED
			      nchanCal = inFile->beam[b].calBandHeader[j].nchan;
			      tdumpCal = inFile->beam[b].calBandHeader[j].dtime;
			      //			      printf("Scaling factor = %g\n",((double)nchanAstro/(double)nchanCal * (double)1.0/(double)nbinCal * tdumpCal /tdumpAstro));
			      //			      printf("pre-fluxScale = %g\n",fluxScale);
			      fluxScale *= ((double)nchanAstro/(double)nchanCal * (double)1.0/(double)nbinCal * tdumpCal /tdumpAstro); 
			    }
			  //			  printf("Here fluxScale = %g\n",fluxScale);
			  // NORMALISATION
			  //			  fluxScale *= ((double)262144.0/(double)128. * (double)1.0/(double)32.0 * 5 /0.983);
			  //			  fluxScale *= ((double)32768.0/(double)128. * (double)1.0/(double)32.0 * 5 /4.997); 
			  //			  /262144./128.*5./0.983/32.


			  //			  fluxScale = 1;
			  final_aa  = fluxScale*creal(finalJ[0][0]);
			  final_bb  = fluxScale*creal(finalJ[1][1]);
			  final_rab = fluxScale*creal(finalJ[1][0]);
  			  final_iab = fluxScale*cimag(finalJ[0][1]);
			  //			  printf("Here with fluxScale = %g\n",fluxScale);

			  out_data[ii+nchan*k*npol]         = final_aa;
			  out_data[ii+nchan*k*npol+nchan]   = final_bb;
			  out_data[ii+nchan*k*npol+2*nchan] = final_rab;
			  out_data[ii+nchan*k*npol+3*nchan] = final_iab;

			  //			  printf("Output [%d] [%d] %.6f %g %g %g %g %g %g %g %g %d %g %g %g\n",j,k,inFile->beam[b].bandData[j].astro_data.freq[ii], aa,bb,rab,iab,final_aa,final_bb,final_rab,final_iab,ichan,fluxCal[ichan].scalAA,fluxCal[ichan].scalBB,fluxScale);
			  //			  exit(1);
			}
		    }
		
		  // Fix up the attributes
		  for (ii=0;ii<nDataAttributes;ii++)
		    {
		      if (strcmp(dataAttributes[ii].key,"UNIT")==0)
			{
			  if (tcal==0)
			    strcpy(dataAttributes[ii].value,"Jy");			  
			  else
			    strcpy(dataAttributes[ii].value,"K");
			}
		    }
		  
		  sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[j].label,b,j, out_data,out_freq,ndump,nchan,1,npol,ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
		  sdhdf_writeObsParams(outFile,inFile->beam[b].bandHeader[j].label,inFile->beamHeader[b].label,j,outObsParams,ndump,1,inFile->primary[0].hdr_defn_version);
		  sdhdf_writeDataWeights(outFile,b,j,dataWts,nchan,ndump,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[j].label);
		  sdhdf_releaseBandData(inFile,b,j,1); 		      
		  free(out_freq); free(out_data); free(dataWts);
		  free(outObsParams);

		}
	      sdhdf_writeBandHeader(outFile,inBandParams,inFile->beamHeader[b].label,inFile->beam[b].nBand,1,inFile->primary[0].hdr_defn_version);
	      free(inBandParams);

	    }
	  sdhdf_writeHistory(outFile,inFile->history,inFile->nHistory);
	  sdhdf_copyRemainder(inFile,outFile,0);

	  sdhdf_closeFile(inFile);
	  sdhdf_closeFile(outFile);
	}
    }

  free(inFile);
  free(outFile);
  free(fluxCal);
  free(polCal);
}



