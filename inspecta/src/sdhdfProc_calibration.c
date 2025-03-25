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
#include <math.h>
#include <string.h>
#include "inspecta.h"
#include "fitsio.h"
#include <complex.h> 

void sdhdf_set_pcm_response(sdhdf_calibration *polCal,double complex J1[2][2],double complex J2[2][2],double complex J3[2][2],double complex J4[2][2],
			    double complex J5[2][2],double complex J6[2][2],double complex J7[2][2]);


void sdhdf_calculate_timedependent_response(sdhdf_calibration *polCal,int nPolCalChan)
{
  int i;
  double complex gain[2][2];
  double complex diffgain[2][2];
  double complex diffphase[2][2];  
  double complex R1[2][2];
  double complex R1_inv[2][2];

  
  for (i=0;i<nPolCalChan;i++)
    {
      sdhdf_complex_matrix_2x2(gain,polCal[i].gain,0,0,polCal[i].gain);
      sdhdf_complex_matrix_2x2(diffgain,cosh(polCal[i].diff_gain) + sinh(polCal[i].diff_gain),0,0,cosh(polCal[i].diff_gain) - sinh(polCal[i].diff_gain));
      sdhdf_complex_matrix_2x2(diffphase,cos(polCal[i].diff_phase) + I*sin(polCal[i].diff_phase),0,0,cos(polCal[i].diff_phase) - I*sin(polCal[i].diff_phase));
      
      sdhdf_copy_complex_matrix_2x2(R1,gain);
      sdhdf_multiply_complex_matrix_2x2(R1,diffgain);
      sdhdf_multiply_complex_matrix_2x2(R1,diffphase);

      // Now we multiply this by the system response
      sdhdf_multiply_complex_matrix_2x2(R1,polCal[i].response_pcm);

      // Now we need the inverse of this
      if (sdhdf_complex_matrix_2x2_inverse(R1,R1_inv)==1) 
	polCal[i].bad=1;
      else
	{
	  // We are able to form an inverse
	  sdhdf_copy_complex_matrix_2x2(polCal[i].td_response,R1_inv);
	  //	  printf("%d set td_response = ",i);
	  //	  sdhdf_display_complex_matrix_2x2(polCal[i].td_response);

	}
    }
}


// Note: this averages the cal signal through the entire observation
//
void sdhdf_noiseSourceOnOff(sdhdf_fileStruct *inFile,int ibeam,int iband,float *freq,float *aa,float *bb,float *re,float *im)
{
  int i,j;
  int nchanCal = inFile->beam[ibeam].calBandHeader[iband].nchan;
  int nsubCal  = inFile->beam[ibeam].calBandHeader[iband].ndump;
  float aa_on,aa_off,bb_on,bb_off,re_on,re_off,im_on,im_off;
  
  for (i=0;i<nchanCal;i++)
    {
      // Should check NCHANFREQ ** FIX ME ** ASSUME 1
      freq[i] = inFile->beam[ibeam].bandData[iband].cal_on_data.freq[i];
      //      printf("Noise source frequency = %g\n",freq[i]);
      aa_on = aa_off = bb_on = bb_off = re_on = re_off = im_on = im_off = 0.;
      for (j=0;j<nsubCal;j++)
	{
	  aa_on += inFile->beam[ibeam].bandData[iband].cal_on_data.pol1[j*nchanCal+i];
	  bb_on += inFile->beam[ibeam].bandData[iband].cal_on_data.pol2[j*nchanCal+i];
	  re_on += inFile->beam[ibeam].bandData[iband].cal_on_data.pol3[j*nchanCal+i];
	  im_on += inFile->beam[ibeam].bandData[iband].cal_on_data.pol4[j*nchanCal+i];

	  aa_off += inFile->beam[ibeam].bandData[iband].cal_off_data.pol1[j*nchanCal+i];
	  bb_off += inFile->beam[ibeam].bandData[iband].cal_off_data.pol2[j*nchanCal+i];
	  re_off += inFile->beam[ibeam].bandData[iband].cal_off_data.pol3[j*nchanCal+i];
	  im_off += inFile->beam[ibeam].bandData[iband].cal_off_data.pol4[j*nchanCal+i];
	}
      aa_on /= (float)nsubCal;
      bb_on /= (float)nsubCal;
      re_on /= (float)nsubCal;    
      im_on /= (float)nsubCal;
      aa_off /= (float)nsubCal;
      bb_off /= (float)nsubCal;
      re_off /= (float)nsubCal;    
      im_off /= (float)nsubCal;

      aa[i] = aa_on-aa_off;
      bb[i] = bb_on-bb_off;
      re[i] = re_on-re_off;
      im[i] = im_on-im_off;
    }
  
}


void sdhdf_calculate_gain_diffgain_diffphase(sdhdf_calibration *polCal,int nPolCalChan)
{
  int i;
  double In,Qn,Un,Vn;
  double Im,Qm,Um,Vm;
  double G2;
  for (i=0;i<nPolCalChan;i++)
    {
      In = polCal[i].stokes_noise_actual[0];
      Qn = polCal[i].stokes_noise_actual[1];
      Un = polCal[i].stokes_noise_actual[2];
      Vn = polCal[i].stokes_noise_actual[3];

      Im = polCal[i].stokes_noise_measured[0];
      Qm = polCal[i].stokes_noise_measured[1];
      Um = polCal[i].stokes_noise_measured[2];
      Vm = polCal[i].stokes_noise_measured[3];

      polCal[i].diff_gain = 0.5*atanh((In*Qm-Im*Qn)/(In*Im - Qn*Qm));
      polCal[i].diff_phase = 0.5*atan2(Um*Vn - Vm*Un,Un*Um + Vn*Vm);

      G2 = sqrt((pow(Im,2) - pow(Qm,2))/(pow(In,2)-pow(Qn,2)));
      polCal[i].gain = sqrt(G2);
      printf("CALIBRATION: gain diffPhase: %g %g %g %g %d\n",polCal[i].freq,polCal[i].gain,polCal[i].diff_gain,polCal[i].diff_phase,nPolCalChan);
    }
    
}
void sdhdf_set_stokes_noise_measured(sdhdf_fileStruct *inFile,int ibeam,sdhdf_calibration *polCal,int nPolCalChan,int normalise,int av,float av1freq,float av2freq)
{
  int i,j;
  int ichan,nchanCal,ndumpCal;
  int iband=-1;
  double freq,f0,f1;
  double chbw;
  double aa,bb,rab,iab;
  double aa_on,bb_on,rab_on,iab_on;
  double aa_off,bb_off,rab_off,iab_off;
  double scaleFactor;
  double tdump;

  for (i=0;i<nPolCalChan;i++)
    {
      // Should set better or using interpolation -- FIX ME
      freq = polCal[i].freq;
      //      printf("At this point: checking freq %g\n",freq);
      iband=-1;
      for (j=0;j<inFile->beam[ibeam].nBand;j++)   // CAN SPEED THIS UP -- DON'T NEED TO CHECK SO OFTEN
	{   
	  f0   = inFile->beam[ibeam].calBandHeader[j].f0; // GEORGE: updated on 8th Sept 2021 to use the cal band header here
	  f1   = inFile->beam[ibeam].calBandHeader[j].f1;
	  //	  printf("At this point: checking j %d %g %g\n",j,f0,f1);
	  if (freq >= f0 && freq <= f1)
	    {
	      iband=j;
	      break;
	    }
	}
      if (iband==-1)
	polCal[i].bad=1;
      else
	{
	  sdhdf_loadBandData(inFile,ibeam,iband,2);
	  sdhdf_loadBandData(inFile,ibeam,iband,3);
	  f0   = inFile->beam[ibeam].bandData[iband].cal_on_data.freq[0];   // FIX ME FOR DUMP FREQ
	  chbw = inFile->beam[ibeam].bandData[iband].cal_on_data.freq[1]-f0;
	  ichan = (int)(((freq - f0)/(double)chbw)+0.5);  // FIX ME -- DOUBLE CHECK THIS AND INTERPOLATE IF NEEDED
	  aa_on = bb_on = rab_on = iab_on = 0.0;
	  aa_off = bb_off = rab_off = iab_off = 0.0;
	  nchanCal = inFile->beam[ibeam].calBandHeader[iband].nchan;
	  ndumpCal = inFile->beam[ibeam].calBandHeader[iband].ndump;

	  if (normalise==1)
	    {
	      // NOT DOING NORMLIASTION HERE NOW *** CHECK LOTS OF TIMES
	      //	      int calBin = 32; // FIX ME -- DON'T HARDCODE

	      //	      calBin = 32;
	      //	      calBin = 2;
	      // bin - 32 => 10^-6
	      // bin = 1 => 10^-4
	      //	      tdump = inFile->beam[ibeam].calBandHeader[iband].dtime;
	      // printf("tdump = %g\n",tdump);
	      //	      fluxScale *= ((double)32768.0/(double)128. * (double)1.0/(double)32.0 * 5 /4.997); 
	      //	      scaleFactor = 1.0/((double)nchanCal/(double)tdump*(double)calBin); // FIX ME -- SHOULD CHECK IF TDUMP IS CORRECT
	      //	      scaleFactor = 1.0/((double)nchanCal/(double)calBin); // FIX ME -- SHOULD CHECK IF TDUMP IS CORRECT
	      //	      scaleFactor = 1.0/((double)calBin*(double)nchanCal); // FIX ME -- SHOULD CHECK IF TDUMP IS CORRECT
	      //	      scaleFactor = 1.0/((double)calBin*(double)nchanCal); // FIX ME -- SHOULD CHECK IF TDUMP IS CORRECT
	      scaleFactor=1;
	    }
	  else
	    scaleFactor=1;

	  for (j=0;j<ndumpCal;j++)
	    {
	      // FIX ME -- SHOULD ACCOUNT FOR ANY WEIGHTING
	      aa_on += inFile->beam[ibeam].bandData[iband].cal_on_data.pol1[j*nchanCal+ichan];
	      aa_off += inFile->beam[ibeam].bandData[iband].cal_off_data.pol1[j*nchanCal+ichan];
	      bb_on += inFile->beam[ibeam].bandData[iband].cal_on_data.pol2[j*nchanCal+ichan];
	      bb_off += inFile->beam[ibeam].bandData[iband].cal_off_data.pol2[j*nchanCal+ichan];
	      rab_on += inFile->beam[ibeam].bandData[iband].cal_on_data.pol3[j*nchanCal+ichan];
	      rab_off += inFile->beam[ibeam].bandData[iband].cal_off_data.pol3[j*nchanCal+ichan];
	      iab_on += inFile->beam[ibeam].bandData[iband].cal_on_data.pol4[j*nchanCal+ichan];
	      iab_off += inFile->beam[ibeam].bandData[iband].cal_off_data.pol4[j*nchanCal+ichan];
	    }
	  sdhdf_releaseBandData(inFile,ibeam,iband,2);
	  sdhdf_releaseBandData(inFile,ibeam,iband,3);

	  aa_on /= (double)ndumpCal;	  aa_off /= (double)ndumpCal;
	  bb_on /= (double)ndumpCal;	  bb_off /= (double)ndumpCal;
	  rab_on /= (double)ndumpCal;	  rab_off /= (double)ndumpCal;
	  iab_on /= (double)ndumpCal;	  iab_off /= (double)ndumpCal;

	  // FIX ME -- should check for OFF > ON issues
	  aa = (aa_on-aa_off)*scaleFactor;
	  bb = (bb_on-bb_off)*scaleFactor;
	  rab = (rab_on-rab_off)*scaleFactor;
	  iab = (iab_on-iab_off)*scaleFactor;

	  // HARDCODE FIX
	  /*
	  aa = 23.3073;
	  bb = 19.66;
	  rab = -16.3943;
	  iab = -4.90985;
	  */
	  
	  polCal[i].stokes_noise_measured[0] = aa+bb;
	  polCal[i].stokes_noise_measured[1] = aa-bb;
	  polCal[i].stokes_noise_measured[2] = 2*rab;
	  polCal[i].stokes_noise_measured[3] = 2*iab;
	  printf("Calibration solution (Stokes noise measured): %g %g %g %g\n",polCal[i].stokes_noise_measured[0],polCal[i].stokes_noise_measured[1],polCal[i].stokes_noise_measured[2],polCal[i].stokes_noise_measured[3]);
	  //	  printf("%d At this point setting %g , iband = %d\n",i,polCal[i].stokes_noise_measured[0],iband);
	}
      //      printf("%d At this point\n",i);
    }
  
  
  // If requested should change to an averaged value
  if (av==1 || av==2)
    {
      double av1=0,av2=0,av3=0,av4=0;
      int n=0;
      for (i=0;i<nPolCalChan;i++)
	{
	  if (polCal[i].freq >= av1freq && polCal[i].freq <= av2freq)
	    {
	      av1 += polCal[i].stokes_noise_measured[0];
	      av2 += polCal[i].stokes_noise_measured[1];
	      av3 += polCal[i].stokes_noise_measured[2];
	      av4 += polCal[i].stokes_noise_measured[3];
	      n++;
	    }
	}
      av1 /= (double)n;
      av2 /= (double)n;
      av3 /= (double)n;
      av4 /= (double)n;
      printf("Setting input measurement of the noise source to (%g,%g,%g,%g)\n",av1,av2,av3,av4);
      for (i=0;i<nPolCalChan;i++)
	{
	  //	  printf("Changing %g to %g\n",polCal[i].stokes_noise_measured[0],av1);
	  polCal[i].stokes_noise_measured[0] = av1;
	  polCal[i].stokes_noise_measured[1] = av2;
	  polCal[i].stokes_noise_measured[2] = av3;
	  polCal[i].stokes_noise_measured[3] = av4;
	}
    }
}


// Find the best estimate of the noise source Stokes parameters for the specified frequency
// FIX ME: Currently only doing a linear interpolation
// FIX ME: IMPORTANT: ISSUE IN THAT SHOULDN'T INTERPOLATE USING BAD ENTIRES (AS THEY ARE ZERO)
void sdhdf_get_pcmcal_stokes(double freq,sdhdf_calibration *polCal,int nPolCalChan,double *actualNoiseStokes)
{
  int i;
  double m,c,y1,y2,x1,x2;
  int set=0;
  double aa,bb,cc,dd;

  aa=bb=cc=dd=0.0; 
  
  for (i=0;i<nPolCalChan-1;i++)
    {
      if (freq <= polCal[i].freq)
	{
	  if (i!=0)
	    {
	      // Stokes I
	      //
	      x1 = polCal[i-1].freq; x2 = polCal[i].freq;
	      y1 = polCal[i-1].noiseSource_I_postResponse; y2 = polCal[i].noiseSource_I_postResponse;
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      aa = m*freq+c; 

	      // Stokes Q
	      //
	      y1 = polCal[i-1].noiseSource_Q_postResponse; y2 = polCal[i].noiseSource_Q_postResponse;
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      bb = m*freq+c; 

	      // Stokes U
	      //
	      y1 = polCal[i-1].noiseSource_U_postResponse; y2 = polCal[i].noiseSource_U_postResponse;
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      cc = m*freq+c; 

	      // Stokes V
	      //
	      y1 = polCal[i-1].noiseSource_V_postResponse; y2 = polCal[i].noiseSource_V_postResponse;
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      dd = m*freq+c; 

	    }
	  else
	    {
	      // Stokes I
	      //
	      x1 = polCal[i+1].freq; x2 = polCal[i].freq;
	      y1 = polCal[i+1].noiseSource_I_postResponse; y2 = polCal[i].noiseSource_I_postResponse;
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      aa = m*freq+c; 

	      // Stokes Q
	      //
	      y1 = polCal[i+1].noiseSource_Q_postResponse; y2 = polCal[i].noiseSource_Q_postResponse;
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      bb = m*freq+c; 

	      // Stokes U
	      //
	      y1 = polCal[i+1].noiseSource_U_postResponse; y2 = polCal[i].noiseSource_U_postResponse;
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      cc = m*freq+c; 

	      // Stokes V
	      //
	      y1 = polCal[i+1].noiseSource_V_postResponse; y2 = polCal[i].noiseSource_V_postResponse;
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      dd = m*freq+c; 


	    }
	  set=1;
	      
	  break;
	}
    }
  if (set==0)
    {
      // Stokes I
      //
      x1 = polCal[i-2].freq; x2 = polCal[i-1].freq;
      y1 = polCal[i-2].noiseSource_I_postResponse; y2 = polCal[i-1].noiseSource_I_postResponse;
      m = (y2-y1)/(x2-x1); c = y1-m*x1;
      aa = m*freq+c; 
      
      // Stokes Q
      //
      y1 = polCal[i-2].noiseSource_Q_postResponse; y2 = polCal[i-1].noiseSource_Q_postResponse;
      m = (y2-y1)/(x2-x1); c = y1-m*x1;
      bb = m*freq+c; 
      
      // Stokes U
      //
      y1 = polCal[i-2].noiseSource_U_postResponse; y2 = polCal[i-1].noiseSource_U_postResponse;
      m = (y2-y1)/(x2-x1); c = y1-m*x1;
      cc = m*freq+c; 
      
      // Stokes V
      //
      y1 = polCal[i-2].noiseSource_V_postResponse; y2 = polCal[i-1].noiseSource_V_postResponse;
      m = (y2-y1)/(x2-x1); c = y1-m*x1;
      dd = m*freq+c; 


      
    }
  actualNoiseStokes[0] = aa;
  actualNoiseStokes[1] = bb;
  actualNoiseStokes[2] = cc;
  actualNoiseStokes[3] = dd;
  
}

//
// Routine to form the system response from a PCM parameterisation
//
void sdhdf_formPCM_response(sdhdf_calibration *polCal,int nPolCalChan)
{
  int i;

  // The response is formed from 7 Jones matrices. The elements of the Jones matrices are complex
  double complex J1[2][2];
  double complex J2[2][2];
  double complex J3[2][2];
  double complex J4[2][2];
  double complex J5[2][2];
  double complex J6[2][2];
  double complex J7[2][2];

  // Note this is "feed" in PSRCHIVE, but probably more relates to the injection of the noise source into the LNAs (for the Parkes UWL)
  double complex feed[2][2];
  double alpha = -45.0*M_PI/180.0;  // FIX ME -- DON'T HARDCODE

  // Hermition of response
  double complex Rdag[2][2];
  
  // Coherency products for the noise source
  double complex rho[2][2];
  double complex R_rho_Rdag[2][2];
  
  sdhdf_complex_matrix_2x2(feed,0-I*sin(alpha),0+I*cos(alpha),0+I*cos(alpha),0+I*sin(alpha));
  //  sdhdf_complex_matrix_2x2(feed,1,0,0,1);
  //  printf("Feed = \n");
  //  sdhdf_display_complex_matrix_2x2(feed);
    
  for (i=0;i<nPolCalChan;i++)
    {
      sdhdf_complex_matrix_2x2(J1,polCal[i].constant_gain,0,0,polCal[i].constant_gain);      
      sdhdf_complex_matrix_2x2(J2,(cosh(polCal[i].constant_diff_gain)+sinh(polCal[i].constant_diff_gain)),0,0,(cosh(polCal[i].constant_diff_gain)-sinh(polCal[i].constant_diff_gain)));
      //      sdhdf_display_complex_matrix_2x2(J2);

      sdhdf_complex_matrix_2x2(J3,cos(polCal[i].constant_diff_phase)+I*sin(polCal[i].constant_diff_phase),0,0,cos(polCal[i].constant_diff_phase)-I*sin(polCal[i].constant_diff_phase));
      sdhdf_complex_matrix_2x2(J4,cosh(polCal[i].constant_b1),sinh(polCal[i].constant_b1),sinh(polCal[i].constant_b1),cosh(polCal[i].constant_b1));
      sdhdf_complex_matrix_2x2(J5,cosh(polCal[i].constant_b2),-I*sinh(polCal[i].constant_b2),I*sinh(polCal[i].constant_b2),cosh(polCal[i].constant_b2));
      sdhdf_complex_matrix_2x2(J6,cos(polCal[i].constant_r1),I*sin(polCal[i].constant_r1),I*sin(polCal[i].constant_r1),cos(polCal[i].constant_r1));
      sdhdf_complex_matrix_2x2(J7,cos(polCal[i].constant_r2),sin(polCal[i].constant_r2),-sin(polCal[i].constant_r2),cos(polCal[i].constant_r2));
      sdhdf_set_pcm_response(&polCal[i],J1,J2,J3,J4,J5,J6,J7);

      // Add in feed response
      sdhdf_copy_complex_matrix_2x2(polCal[i].response_pcm_feed,polCal[i].response_pcm);
      sdhdf_multiply_complex_matrix_2x2(polCal[i].response_pcm_feed,feed); // FIX ME -- NOT RESPONSE_PCM


      sdhdf_complex_matrix_2x2_dagger(polCal[i].response_pcm_feed,Rdag);

      //      printf("%d Response_pcm = ",i);
      //      sdhdf_display_complex_matrix_2x2(polCal[i].response_pcm_feed);

      
      // Now we have the response of the system in terms of the noise source.
      // We now can apply that to the "actual" Stokes parameters of the noise source
      sdhdf_complex_matrix_2x2(rho,
			       0.5*(1.0+polCal[i].noiseSource_QoverI),
			       0.5*(polCal[i].noiseSource_UoverI - I*polCal[i].noiseSource_VoverI),
			       0.5*(polCal[i].noiseSource_UoverI + I*polCal[i].noiseSource_VoverI),
			       0.5*(1.0-polCal[i].noiseSource_QoverI));
      //      printf("%d rho = ",i);
      //      sdhdf_display_complex_matrix_2x2(rho);

      sdhdf_copy_complex_matrix_2x2(R_rho_Rdag,polCal[i].response_pcm_feed);
      sdhdf_multiply_complex_matrix_2x2(R_rho_Rdag,rho);
      //      printf("%d J rho = ",i);
      // This doesn't seem identical to the PSRCHIVE version -- CHECK?   FIX ME
      //      sdhdf_display_complex_matrix_2x2(R_rho_Rdag);

      sdhdf_multiply_complex_matrix_2x2(R_rho_Rdag,Rdag); // R * rho * R^dagger
      //      printf("%d Response_rho_Rdag = ",i);
      //      sdhdf_display_complex_matrix_2x2(R_rho_Rdag);

      polCal[i].stokes_noise_actual[0] = creal(R_rho_Rdag[0][0])+creal(R_rho_Rdag[1][1]);
      polCal[i].stokes_noise_actual[1] = creal(R_rho_Rdag[0][0])-creal(R_rho_Rdag[1][1]);
      polCal[i].stokes_noise_actual[2] = 2*creal(R_rho_Rdag[1][0]);
      polCal[i].stokes_noise_actual[3] = 2*cimag(R_rho_Rdag[0][1]);
      //      printf("Noise source postResponse = (%g,%g,%g,%g)\n",polCal[i].stokes_noise_actual[0],
      //	     polCal[i].stokes_noise_actual[1],polCal[i].stokes_noise_actual[2],polCal[i].stokes_noise_actual[3]);
      
    }
}

void sdhdf_set_pcm_response(sdhdf_calibration *polCal,double complex J1[2][2],double complex J2[2][2],double complex J3[2][2],double complex J4[2][2],
			    double complex J5[2][2],double complex J6[2][2],double complex J7[2][2])
{
  sdhdf_copy_complex_matrix_2x2(polCal->response_pcm,J1);
  sdhdf_multiply_complex_matrix_2x2(polCal->response_pcm,J2);  
  sdhdf_multiply_complex_matrix_2x2(polCal->response_pcm,J3);
  sdhdf_multiply_complex_matrix_2x2(polCal->response_pcm,J4);
  sdhdf_multiply_complex_matrix_2x2(polCal->response_pcm,J5);
  sdhdf_multiply_complex_matrix_2x2(polCal->response_pcm,J6);
  sdhdf_multiply_complex_matrix_2x2(polCal->response_pcm,J7);
}


//
// Routine to load a PCM file (FITS format)
// to model cross-coupling within the receiver and noise source system
//
void sdhdf_loadPCM(sdhdf_calibration *polCal,int *nPolCalChan,char *observatory, char *rcvr,char *pcmFile,int av,float av1freq,float av2freq,int pcmClear)
{
  int status=0;
  fitsfile *fptr;
  char fname[1024];
  char runtimeDir[1024];
  int nchan_cal_poln;
  int nchan_pcm;
  
  int colnum_data,colnum_freq;
  float *data;
  double *dataFreq;
  
  float n_fval=0;
  double n_dval=0;
  
  int initflag=0;
  int i,j;

  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdfProc_calibration requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  sprintf(fname,"%s/observatory/%s/calibration/%s/%s",runtimeDir,observatory,rcvr,pcmFile);
  printf("Opening: %s\n",fname);
  
  fits_open_file(&fptr,fname,READONLY,&status);
  fits_report_error(stderr,status);

  // Obtain the Stokes parameters of the noise source
  // NOTE THAT THIS DOESN'T HAVE A REQUENCY AXIS - That's probably a mistake in the FITS file definition
  //
  fits_movnam_hdu(fptr,BINARY_TBL,"CAL_POLN",1,&status);
  fits_read_key(fptr,TINT,"NCHAN",&nchan_cal_poln,NULL,&status);
  printf("Number of channels = %d\n",nchan_cal_poln);
  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum_data,&status);

  // Each entry has Q/I, U/I and V/I (i.e., 3 parameters)
  data = (float *)malloc(sizeof(float)*nchan_cal_poln*3);
  fits_read_col(fptr,TFLOAT,colnum_data,1,1,nchan_cal_poln*3,&n_fval,data,&initflag,&status);
  for (i=0;i<nchan_cal_poln;i++)
    {
      if (pcmClear==1)
	{
	  polCal[i].noiseSource_QoverI = 1.0;
	  polCal[i].noiseSource_UoverI = 0.0;
	  polCal[i].noiseSource_VoverI = 0.0;
	}
      else
	{
	  polCal[i].noiseSource_QoverI = data[i*3];
	  polCal[i].noiseSource_UoverI = data[i*3+1];
	  polCal[i].noiseSource_VoverI = data[i*3+2];      
	}
    }
  
  *nPolCalChan = nchan_cal_poln;

  free(data);

  
  // Now load the PCM parameteres
  fits_movnam_hdu(fptr,BINARY_TBL,"FEEDPAR",1,&status);
  fits_report_error(stderr,status);
  
  fits_read_key(fptr,TINT,"NCHAN",&nchan_pcm,NULL,&status);
  if (nchan_pcm != nchan_cal_poln)
    {
      printf("ERROR: In sdhdfProc_calibration: nchan_pcm != nchan_cal_poln\n");
      printf("nchan_pcm = %d\n",nchan_pcm);
      printf("nchan_cal_poln = %d\n",nchan_cal_poln);
      exit(1);      
    }
  //
  // Have 7 parameters to define the feed and noise source system
  //
  data = (float *)malloc(sizeof(float)*nchan_pcm*7);
  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum_data,&status);


  fits_report_error(stderr,status);
  // NOTE SOME OF THESE SEEM TO BE NAN - BE CAREFUL READING IN -- FIX ME
  fits_read_col(fptr,TFLOAT,colnum_data,1,1,nchan_pcm*7,&n_fval,data,&initflag,&status);
  fits_report_error(stderr,status);
  for (i=0;i<nchan_cal_poln;i++)
    {
      if (isnan(data[i*7+6]) || isnan(data[i*7+5]) || isnan(data[i*7+5]) ||
	  isnan(data[i*7+4]) || isnan(data[i*7+3]) || isnan(data[i*7+2]) || isnan(data[i*7+1]) || isnan(data[i*7+0]))
	{
	  polCal[i].constant_gain = 0;
	  polCal[i].constant_diff_gain = 0;
	  polCal[i].constant_diff_phase = 0;
	  polCal[i].constant_b1 = 0;
	  polCal[i].constant_b2 = 0;
	  polCal[i].constant_r1 = 0;
	  polCal[i].constant_r2 = 0;
	  polCal[i].bad = 1;
	}
      else if (pcmClear==1) // Adding in a method to switch off the PCM
	{
	  //	  printf("I am in here\n");
	  polCal[i].constant_gain = data[i*7];
	  polCal[i].constant_diff_gain = data[i*7+1];
	  polCal[i].constant_diff_phase = data[i*7+2];
	  printf("Got %g %g %g\n",polCal[i].constant_gain,polCal[i].constant_diff_gain,polCal[i].constant_diff_phase);
	  polCal[i].constant_b1 = 0;
	  polCal[i].constant_b2 = 0;
	  polCal[i].constant_r1 = 0;
	  polCal[i].constant_r2 = 0;
	  polCal[i].bad = 0;	 
	}
      else
	{
	  polCal[i].constant_gain = data[i*7];
	  polCal[i].constant_diff_gain = data[i*7+1];
	  polCal[i].constant_diff_phase = data[i*7+2];
	  polCal[i].constant_b1 = data[i*7+3];
	  polCal[i].constant_b2 = data[i*7+4];
	  polCal[i].constant_r1 = data[i*7+5];
	  polCal[i].constant_r2 = data[i*7+6];
	  polCal[i].bad = 0;
	}
    }
  free(data);

  
  // Read the frequency axis
  dataFreq = (double *)malloc(sizeof(double)*nchan_pcm);
  fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&colnum_freq,&status);
  fits_read_col(fptr,TDOUBLE,colnum_freq,1,1,nchan_pcm,&n_dval,dataFreq,&initflag,&status);
  fits_report_error(stderr,status);
  for (i=0;i<nchan_cal_poln;i++)
    polCal[i].freq = dataFreq[i];
  
  free(dataFreq);


  // Do we require any averaging?
    if (av==1 || av==2)
    {
      double av1,av2,av3,av4,av5,av6,av7,av8,av9,av10;
      int n=0;
      av1=av2=av3=av4=av5=av6=av7=av8=av9=av10=0;
      for (i=0;i<nchan_cal_poln;i++)
	{
	  if (polCal[i].freq >= av1freq  && polCal[i].freq <= av2freq)
	    {
	      av1 += polCal[i].noiseSource_QoverI;
	      av2 += polCal[i].noiseSource_UoverI;
	      av3 += polCal[i].noiseSource_VoverI;
	      av4 += polCal[i].constant_gain;
	      av5 += polCal[i].constant_diff_gain;
	      av6 += polCal[i].constant_diff_phase;
	      av7 += polCal[i].constant_b1;
	      av8 += polCal[i].constant_b2;
	      av9 += polCal[i].constant_r1;
	      av10 += polCal[i].constant_r2;
	      n++;
	    }
	}
      av1 /= (double)n;
      av2 /= (double)n;
      av3 /= (double)n;
      av4 /= (double)n;       av5 /= (double)n;      av6 /= (double)n;       av7 /= (double)n;
      av8 /= (double)n;       av9 /= (double)n;      av10 /= (double)n;
      printf("Setting PCM noise source to (%g,%g,%g)\n",av1,av2,av3);
      printf("Setting PCM 7-parameters to (%g,%g,%g,%g,%g,%g,%g)\n",av4,av5,av6,av7,av8,av9,av10);
      for (i=0;i<nchan_cal_poln;i++)
	{
	  polCal[i].noiseSource_QoverI = av1;
	  polCal[i].noiseSource_UoverI = av2;
	  polCal[i].noiseSource_VoverI = av3;
	  polCal[i].constant_gain = av4;
	  polCal[i].constant_diff_gain = av5;
	  polCal[i].constant_diff_phase = av6;
	  polCal[i].constant_b1 = av7;
	  polCal[i].constant_b2 = av8;
	  polCal[i].constant_r1 = av9;
	  polCal[i].constant_r2 = av10;	  
	}      
    }

  

  
  
  
  fits_close_file(fptr,&status);
  
}
/*
int sdhdf_loadTcal(sdhdf_tcal_struct *tcalData,char *fname)
{
  FILE *fin;
  int n=0;
  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open filename >%s<\n",fname);
      exit(1);
    }
  while (!feof(fin))
    {
      if (fscanf(fin,"%lf %lf %lf",&(tcalData[n].freq),&(tcalData[n].tcalA),&(tcalData[n].tcalB))==3)
	n++;
    }
  fclose(fin);
  return n;
}
*/


void sdhdf_get_tcal(sdhdf_tcal_struct *tcalData,int n,double f0,double *tcalA,double *tcalB)
{
  int i;
  double m,c;
  for (i=0;i<n-1;i++)
    {
      // Should do an interpolation here
      if (f0 >= tcalData[i].freq && f0 < tcalData[i+1].freq)
	{
	  m = (tcalData[i].tcalA-tcalData[i+1].tcalA)/(tcalData[i].freq-tcalData[i+1].freq);
	  c = tcalData[i].tcalA-m*tcalData[i].freq;
	  *tcalA = m*f0+c;

	  m = (tcalData[i].tcalB-tcalData[i+1].tcalB)/(tcalData[i].freq-tcalData[i+1].freq);
	  c = tcalData[i].tcalB-m*tcalData[i].freq;
	  *tcalB = m*f0+c;

	  break;
	}
    }
}

//
// This agrees with compute_stokes in simpol.C in PSRCHIVE
//
void sdhdf_convertStokes(float p1,float p2,float p3,float p4,float *stokesI,float *stokesQ,float *stokesU,float *stokesV)
{
  *stokesI = p1 + p2;
  *stokesQ = p1 - p2;
  *stokesU = 2*p3;
  *stokesV = 2*p4;
}

void sdhdf_loadTcal(sdhdf_fluxCalibration *fluxCal,int *nFluxCal,char *observatory,char *rcvr,char *tcalFile)
{
  char fname[1024];
  char runtimeDir[1024];
  int nchan;
  int i,j;
  FILE *fin;
  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdfProc_calibration requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  sprintf(fname,"%s/observatory/%s/calibration/%s/%s",runtimeDir,observatory,rcvr,tcalFile);
  printf("Opening: %s\n",fname);
  fin = fopen(fname,"r");
  while (!feof(fin))
    {
      if (fscanf(fin,"%lf %f %f",&(fluxCal[*nFluxCal].freq),&(fluxCal[*nFluxCal].scalAA),&(fluxCal[*nFluxCal].scalBB))==3)
	{
	  // Note that we're assuming that we add AA + BB to get Stokes I
	  fluxCal[*nFluxCal].scalAA /= 2.;  // SHOULD CHECK THIS CAREFULLY
	  fluxCal[*nFluxCal].scalBB /= 2.;
	  fluxCal[*nFluxCal].type = 2;
	  (*nFluxCal)++;
	}
    }
  fclose(fin);
}
//
// Routine to load a flux calibration file
// if "auto" then automatically choose a sensible file
//
void sdhdf_loadFluxCal(sdhdf_fluxCalibration *fluxCal,int *nFluxCalChan,char *observatory, char *rcvr,char *fluxCalFile,float mjdObs)
{
  int status=0;
  fitsfile *fptr;
  char fname[1024];
  char runtimeDir[1024];
  int nchan;
  
  int colnum_data,colnum_freq;
  float *data;
  double *dataFreq;
  
  float n_fval=0;
  double n_dval=0;
  
  int initflag=0;
  int i,j;

  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdfProc_calibration requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  if (strcmp(fluxCalFile,"auto")==0)
    {
      FILE *fin;
      char line[1024];
      char loadName[1024];
      float mjd0,mjd1;
      sprintf(fname,"%s/observatory/%s/calibration/%s/automatic_selection_fluxcal.txt",runtimeDir,observatory,rcvr);
      if (!(fin = fopen(fname,"r")))
	{
	  printf("Error: unable to open file %s\n",fname);
	  exit(1);
	}
      while (!feof(fin))
	{
	  if (fgets(line,1024,fin)!=NULL)
	    {
	      if (line[0]!='#') // Not a comment line
		{
		  sscanf(line,"%s %f %f",loadName,&mjd0,&mjd1);
		  if (mjd0 < mjdObs && mjd1 > mjdObs)
		    {
		      strcpy(fluxCalFile,loadName);
		      break;
		    }
		}
	    }
	}
      fclose(fin);
    }
  sprintf(fname,"%s/observatory/%s/calibration/%s/%s",runtimeDir,observatory,rcvr,fluxCalFile);
  printf("Opening: %s\n",fname);
  
  fits_open_file(&fptr,fname,READONLY,&status);
  fits_report_error(stderr,status);

  fits_movnam_hdu(fptr,BINARY_TBL,"FLUX_CAL",1,&status);
  fits_read_key(fptr,TINT,"NCHAN",&nchan,NULL,&status);
  //  printf("Number of channels = %d\n",nchan);
  fits_get_colnum(fptr,CASEINSEN,"S_CAL",&colnum_data,&status);

  data = (float *)malloc(sizeof(float)*nchan*2); // HAVE AA and BB
  fits_read_col(fptr,TFLOAT,colnum_data,1,1,nchan*2,&n_fval,data,&initflag,&status);
  for (i=0;i<nchan;i++)
    {
      fluxCal[i].scalAA = data[i]/1000.; // /1000. because stored in mJy
      fluxCal[i].scalBB = data[nchan+i]/1000.; // Note stored as all AA then all BB
      //      printf("scal = %d %g %g\n",i,fluxCal[i].scalAA,fluxCal[i].scalBB);
      fluxCal[i].type = 1;
    }

  *nFluxCalChan = nchan;
  free(data);
  
  // Read the frequency axis
  dataFreq = (double *)malloc(sizeof(double)*nchan);
  fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&colnum_freq,&status);
  fits_read_col(fptr,TDOUBLE,colnum_freq,1,1,nchan,&n_dval,dataFreq,&initflag,&status);
  fits_report_error(stderr,status);
  for (i=0;i<nchan;i++)
    fluxCal[i].freq = dataFreq[i];
  
  free(dataFreq);

  fits_close_file(fptr,&status);
}
