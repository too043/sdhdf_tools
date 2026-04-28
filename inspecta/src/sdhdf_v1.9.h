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
// These are definitions for SDHDF v1.9 format
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_STRLEN 512
#define MAX_ATTRIBUTES 64
#define MAX_HISTORY 64

// Structure is as follows:
//
// Each data set and group has attributes
//
// filename/
//          obs_meta/
//                  primary header            
//                  software_versions         
//                  history                   
//          obs_config/                     (this group is currently ignored by INSPECTA)
//                  backend_config            
//                  cal_backend_config        
//          beam_XX/
//                  beam_header               
//                  band_header               
//                  band_YY/
//                          astronomy/
//                                data
//                                frequency
//                                  obs_params  
//                          calibrator/
//                              cal_data_on
//                              cal_data_off
//                              cal_frequency
//                              obs_meta/
//                                       cal_band_header 
//                                       cal_obs_params   
//                              cal_proc/
//                                       cal_proc_tsys
//                                       cal_proc_diff_gain
//                                       cal_proc_diff_phase


//
// Attribute structure
//
typedef struct sdhdf_attributes_struct { //2
  char name[MAX_STRLEN];
  char key[MAX_STRLEN];
  char unit[MAX_STRLEN];
  char value[MAX_STRLEN];

  // not used but included while we get it working
  float fvalue;
  int   ivalue;
  int   attributeType; // 0 = string, 1 = float, 2 = int
  //
} sdhdf_attributes_struct;

/*
typedef struct sdhdf_attributes_struct {
  char  key[MAX_STRLEN];
  char  value[MAX_STRLEN];
  float fvalue;
  int   ivalue;
  int   attributeType; // 0 = string, 1 = float, 2 = int
} sdhdf_attributes_struct;
*/

//
// Parameters relating to spectral dumps
//

typedef struct sdhdf_spectralDumpsStruct {
  int    nchan;                // Number of channels
  int    ndump;                // Number of spectral dumps
  int    npol;                 // Number of polarisations allocated
  
  int    freqAllocatedMemory;      // 0 = no, otherwise 1
  int    flagAllocatedMemory;      // 0 = no, otherwise 1
  int    dataWeightsAllocatedMemory;      // 0 = no, otherwise 1
  int    pol1AllocatedMemory;      // 0 = no, otherwise 1
  int    pol2AllocatedMemory;      // 0 = no, otherwise 1
  int    pol3AllocatedMemory;      // 0 = no, otherwise 1
  int    pol4AllocatedMemory;      // 0 = no, otherwise 1
  
  double  *freq;               // Frequency (frequency axis per spectral dump)
  int    nFreqDumps;           // Number of frequency dumps 
  float  *dataWeights;         // Weight information for channels and spectral dumps
  unsigned char *flag;                // Flagging information only for channels
  
  float  *pol1;
  float  *pol2;
  float  *pol3;
  float  *pol4;

  sdhdf_attributes_struct freq_attr[MAX_ATTRIBUTES];
  int nFreqAttributes;

  sdhdf_attributes_struct data_attr[MAX_ATTRIBUTES];
  int nDataAttributes;

} sdhdf_spectralDumpsStruct;



// ***********************************************
// obs_meta definitions
// ***********************************************

// Primary header
typedef struct sdhdf_primaryHeaderStruct {
  char date[20];
  int  hdr_prior4; // = 1 if before v4.0 or 0 otherwise
  char hdr_defn[20];
  char hdr_defn_version[20];
  char file_format[20];
  char file_format_version[20];
  long int sched_block_id;
  char cal_mode[20];
  char instrument[20];
  char observer[20];
  char pid[20];
  char rcvr[20];
  char telescope[64];
  char utc0[64];
  long nbeam;
} sdhdf_primaryHeaderStruct;

// Software versions
typedef struct sdhdf_softwareVersionsStruct {
  char proc_name[64];
  char software[64];
  char software_descr[64];
  char software_version[64];
} sdhdf_softwareVersionsStruct;

// History group
typedef struct sdhdf_historyStruct {
  char date[20];
  char proc_name[64];
  char proc_descr[64];
  char proc_args[1024];
  char proc_host[64];
} sdhdf_historyStruct;

// ***********************************************
// beam definitions
// ***********************************************

typedef struct sdhdf_beamHeaderStruct {
  char label[MAX_STRLEN];
  int  nBand;
  char source[MAX_STRLEN];
  
} sdhdf_beamHeaderStruct;

// ***********************************************
// Band definitions
// ***********************************************

// Band header group
typedef struct sdhdf_bandHeaderStruct {
  char label[64];
  double fc;
  double f0;
  double f1;
  int    nchan;
  int    npol;
  char pol_type[20]; // MOVED
  double dtime;
  int    ndump;
} sdhdf_bandHeaderStruct;

// backend_config
typedef struct sdhdf_backendConfigStruct {
  char backend_phase[64];
  char bw[64];
  char bytes_per_second[64];
  char data_host[64];
  char data_mcast[64];
  char data_port[64];
  char dsb[64];
  char dpsr_git_hash[64];
  char encoding[64];
  char endian[64];
  char end_channel[64];
  char freq[64];
  char hdr_size[64];
  char hdr_version[64];
  char instrument[64];
  char nbit[64];
  char nchan[64];
  char ndim[64];
  char npol[64];
  char psrchive_git_hash[64];
  char psrdada_git_hash[64];
  char rcvr_hand[64];
  char receiver[64];
  char resolution[64];
  char spip_git_hash[64];
  char start_channel[64];
  char state[64];
  char stream_beam_id[64];
  char stream_host[64];
  char stream_subband_id[64];
  char telescope[64];
  char tsamp[64];
  char udp_nsamp[64];
  char adaptive_filter[64];
  char adaptive_filter_epsilon[64];
  char adaptive_filter_nchan[64];
  char adaptive_filter_nsamp[64];
  char calFreq[64];
  char cal_duty_cycle[64];
  char cal_epoch[64];
  char cal_freq[64];
  char cal_phase[64];
  char cal_signal[64];
  char continuum_outnchan[64];
  char continuum_outstokes[64];
  char continuum_outtsamp[64];
  char continuum_outtsubint[64];
  char dec[64];
  char dec_offset[64];
  char observer[64];
  char obs_offset[64];
  char perform_baseband[64];
  char perform_continuum[64];
  char perform_fold[64];
  char perform_search[64];
  char perform_spectral_line[64];
  char perform_vlbi[64];
  char pid[64];
  char ra[64];
  char ra_offset[64];
  char scan_id[64];
  char sched_block_id[64];
  char source[64];
  char tobs[64];
  char tsys_avg_time[64];
  char tsys_freq_res[64];
  char utc_start[64];
  char utc_stop[64];
  char nant[64];
  char nbin[64];
  char picoseconds[64];
  char osamp_numerator[64];
  char osamp_denominator[64];
  char order[64];
  char file_size[64];
  
} sdhdf_backendConfigStruct;


// OBS_PARAMS group
typedef struct sdhdf_obsParamsStruct {
  double timeElapsed;
  double dtime;
  char timedb[64];
  double mjd;
  double fractional_mjd;
  char utc[64];
  char ut_date[64];
  char local_time[64];
  char raStr[64];
  char decStr[64];
  double raDeg;
  double decDeg;
  double raOffset;
  double decOffset;
  double gl;
  double gb;
  double az;
  double ze;
  double el;
  double az_drive_rate;
  double ze_drive_rate;
  double hourAngle;
  double paraAngle;
  double windDir;
  double windSpd;
} sdhdf_obsParamsStruct;

// ***********************************************
// Parameters relating to the data file/observation
// ***********************************************

// Band information
typedef struct sdhdf_bandStruct {
  int haveCal;                // = 0 if no, = 1 if yes
  
  sdhdf_obsParamsStruct     *astro_obsHeader;
  sdhdf_obsParamsStruct     *cal_obsHeader;
  int astro_obsHeaderAllocatedMemory;
  int cal_obsHeaderAllocatedMemory;
  int nAstro_obsHeader;
  int nCal_obsHeader;

  sdhdf_attributes_struct astro_obsHeaderAttr[MAX_ATTRIBUTES];
  int nAstro_obsHeaderAttributes;
  sdhdf_attributes_struct astro_obsHeaderAttr_freq[MAX_ATTRIBUTES];
  int nAstro_obsHeaderAttributes_freq;
  sdhdf_attributes_struct cal_obsHeaderAttr[MAX_ATTRIBUTES];
  int nCal_obsHeaderAttributes;
  
  sdhdf_spectralDumpsStruct astro_data;
  sdhdf_attributes_struct astro_dataAttr[MAX_ATTRIBUTES];
  int nAstro_dataAttributes;

  sdhdf_spectralDumpsStruct cal_on_data;
  sdhdf_attributes_struct cal_on_dataAttr[MAX_ATTRIBUTES];
  int cal_on_dataAttributes;

  sdhdf_spectralDumpsStruct cal_off_data;
  sdhdf_attributes_struct cal_off_dataAttr[MAX_ATTRIBUTES];
  int cal_off_dataAttributes;

  sdhdf_spectralDumpsStruct cal_proc_tsys;
  sdhdf_spectralDumpsStruct cal_proc_diff_gain;
  sdhdf_spectralDumpsStruct cal_proc_diff_phase;
  
} sdhdf_bandStruct;


// Beam information
typedef struct sdhdf_beamStruct {

  sdhdf_bandHeaderStruct *bandHeader;
  sdhdf_attributes_struct bandHeaderAttr[MAX_ATTRIBUTES];
  int nBandHeaderAttributes;


  sdhdf_bandStruct       *bandData;
  sdhdf_attributes_struct bandDataAttr[MAX_ATTRIBUTES];
  int nBandDataAttributes;

  int nBand;
  int bandAllocatedMemory;   // 0 = no, 1 = yes

  sdhdf_bandHeaderStruct *calBandHeader;
  int nCalBand;
  int calBandAllocatedMemory; // 0 = no, 1 = yes
  
} sdhdf_beamStruct;

//
// File information
//
typedef struct sdhdf_fileStruct {
  char  fname[MAX_STRLEN];      // Input filename
  int   fileOpen;               // 0 = file closed, 1 = file opened
  hid_t fileID;                 // File pointer 

  // Primary information
  int nPrimary;
  int primaryAllocatedMemory;              // 0 = no, 1 = yes
  sdhdf_primaryHeaderStruct *primary;
  sdhdf_attributes_struct primaryAttr[MAX_ATTRIBUTES];
  int nPrimaryAttributes;

  // Config information
  char cal_epoch[1024];
  double cal_freq;
  double cal_duty_cycle;
  double cal_phase;
  //  sdhdf_backendConfigStruct config[512]; // Note currently will just load a single line

  
  // Software versions
  int nSoftware;
  int softwareAllocatedMemory;             // 0 = no, 1 = yes
  sdhdf_softwareVersionsStruct *software; 
  sdhdf_attributes_struct softwareAttr[MAX_ATTRIBUTES];
  int nSoftwareAttributes;

  // File history information
  int nHistory;
  int historyAllocatedMemory;              // 0 = no, 1 = yes
  sdhdf_historyStruct *history;
  sdhdf_attributes_struct historyAttr[MAX_ATTRIBUTES];
  int nHistoryAttributes;

  // Beam information
  int nBeam;
  int beamAllocatedMemory;                 // 0 = no, 1 = yes

  sdhdf_beamHeaderStruct  *beamHeader;
  sdhdf_attributes_struct beamHeaderAttr[MAX_ATTRIBUTES];
  int nBeamHeaderAttributes;

  sdhdf_beamStruct *beam;
  sdhdf_attributes_struct beamAttr[MAX_ATTRIBUTES];
  int nBeamAttributes;

  sdhdf_attributes_struct fileAttr[MAX_ATTRIBUTES];
  int nFileAttributes;

} sdhdf_fileStruct;
