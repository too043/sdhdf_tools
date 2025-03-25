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
// These are definitions for SDHDF v1.8 format
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_STRLEN 1024

// History group
typedef struct sdhdf_historyStruct {
  char date[20];
  char proc_name[64];
  char proc_descr[64];
  char proc_args[64];
  char proc_host[64];
} sdhdf_historyStruct;

// Software versions
typedef struct sdhdf_softwareVersionsStruct {
  char proc_name[64];
  char software[64];
  char software_descr[64];
  char software_version[64];
} sdhdf_softwareVersionsStruct;

// Primary header
typedef struct sdhdf_primaryHeaderStruct {
  char date[20];
  char hdr_defn[20];
  char hdr_defn_version[20];
  char file_format[20];
  char file_format_version[20];
  long int n_sub;
  long int sched_block_id;
  double tsamp;
  double totInt;
  char cal_mode[20];
  char dec[64];
  char instrument[20];
  char observer[20];
  char pid[20];
  char pol_type[20];
  char ra[64];
  char rcvr[20];
  char source[64];
  char telescope[64];
  char utc0[64];
  long nbeam;
} sdhdf_primaryHeaderStruct;



// Band header group
typedef struct sdhdf_bandHeaderStruct {
  char label[12];
  double fc;
  double f0;
  double f1;
  long int nchan;
  long int npol;
  double dtime;
  long int ndumps;
  long int nbeams;
} sdhdf_bandHeaderStruct;

// Attributes
typedef struct sdhdf_frequency_attributes_struct {
  char frame[MAX_STRLEN];
  char unit[MAX_STRLEN];
} sdhdf_frequency_attributes_struct;

typedef struct sdhdf_data_attributes_struct {
  char unit[MAX_STRLEN];
} sdhdf_data_attributes_struct;

