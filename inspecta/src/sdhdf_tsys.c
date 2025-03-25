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

//
// Software to process and view the calibrator
//
// Usage:
// sdhdf_tsys
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"
#include <cpgplot.h>

void help()
{
  printf("sdhdf_tsys\n");
  printf("Provides ways to visualise the system temperature variations\n");
}


int main(int argc,char *argv[])
{
  int i,j,k;
  char **fname; // [MAX_STRLEN];
  sdhdf_fileStruct *inFile;
  int iband=0,nchan,idump;
  int ndump=0;
  int beam=0;
  int totChan=0;
  int av=0;
  int nFiles=0;

  fname = (char **)malloc(sizeof(char *)*MAX_FILES);
  for (i=0;i<MAX_FILES;i++)
    fname[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);
    
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);
    
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-sb")==0 || strcmp(argv[i],"-b")==0)
	sscanf(argv[++i],"%d",&iband);
      else
  	strcpy(fname[nFiles++],argv[++i]);	
    }
  printf("Number of input files = %d\n",nFiles);

  for (i=0;i<nFiles;i++)
    {    
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_loadMetaData(inFile);
      ndump = inFile->beam[beam].calBandHeader[iband].ndump;
      nchan = inFile->beam[beam].calBandHeader[iband].nchan;
      sdhdf_loadBandData(inFile,beam,iband,2);
      sdhdf_loadBandData(inFile,beam,iband,3);
      printf("%s nchan/ndump = %d/%d\n",fname[i],ndump,nchan);
      
      sdhdf_closeFile(inFile);
    }

  
  
  free(inFile);

  for (i=0;i<MAX_FILES;i++)
    free(fname[i]);
  free(fname);

}

  
 
