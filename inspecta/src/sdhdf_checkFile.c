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
// Software to provide a quick-look at the metadata in a SDHDF file.
//
// Usage:
// sdhdf_checkFile <filename.hdf> <filename2.hdf> ...
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"

#define VERSION "v1.0"
#define MAX_HISTORY_LOG 512
#define MAX_COMMANDS 16

void checkVersion(char *ver);
void checkLog(sdhdf_fileStruct *inFile);


typedef struct historyLogStruct {
  int type;
  double mjdStart;
  double mjdEnd;
  char detailedReport[1024];
  char report[1024];
} historyLogStruct;

void help()
{
  printf("sdhdf_checkFile %s (INSPECTA %s)\n",VERSION,SOFTWARE_VER);
  printf("Authors: G. Hobbs\n");

  printf("\n\n");
  printf("-h        This help\n");
  printf("Filenames are given on the command line\n\n");
  printf("Example\n");
  printf("sdhdf_checkFile *.hdf\n");
  exit(1);

}

int main(int argc,char *argv[])
{
  int i,j,beam;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  sdhdf_fileStruct *inFile;
  char cmdList[MAX_COMMANDS][MAX_STRLEN];
  int nCommands=0;
  int iband=0;

  // Display help if no commands given
  if (argc==1)
    help();
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-h")==0)
	help();      
      else if (strcmp(argv[i],"-c")==0)
	strcpy(cmdList[nCommands++],argv[++i]);
      else
	{
	  strcpy(fname[nFiles++],argv[i]);
	  if (nFiles == MAX_FILES)
	    {
	      printf("ERROR: Maximum number of files = %d\n",MAX_FILES);
	      exit(1);
	    }
	}
    }

  inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));
  
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      if (sdhdf_openFile(fname[i],inFile,1)==-1)
	printf("Warning: unable to open file >%s<. Skipping\n",fname[i]);
      else
	{
	  sdhdf_loadMetaData(inFile);
	  if (nCommands==0) // Produce a report
	    {
	      printf("Report for %s\n",inFile->fname);
	      printf(" ... sdhdf version: %s\n",inFile->primary[0].hdr_defn_version);
	      checkVersion(inFile->primary[0].hdr_defn_version);       
	      checkLog(inFile);
	    }
	  else
	    {
	      printf("%s ",inFile->fname);
	      for (j=0;j<nCommands;j++)
		{
		  if (strcasecmp(cmdList[j],"ver")==0)
		    printf("%s ",inFile->primary[0].hdr_defn_version);
		  else if (strcasecmp(cmdList[j],"noisesource")==0)
		    printf("%s ",inFile->primary[0].cal_mode);
		}
	      printf("\n");
	    }
	  sdhdf_closeFile(inFile);
	}
    }

  free(inFile);
}

void checkLog(sdhdf_fileStruct *inFile)
{
  int i;
  int ibeam = 0; // FIX ME
  int iband = 0; // FIX ME
  char telescope[1024];
  char observatoryDir[1024]="NULL";
  FILE *fin;
  char loadLine[1024];
  char fname[1024],fname1[2024];
  char runtimeDir[1024];
  int nLog=0;
  double mjd;
  historyLogStruct *historyLog;
  char *tok;
  char entry[1024];
  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdf_checkFile requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));

  sprintf(fname1,"%s/observatory/observatories.list",runtimeDir);
  fin = fopen(fname1,"r");
  while (!feof(fin))
    {
      if (fgets(loadLine,1024,fin)!=NULL)
	{
	  if (loadLine[0]!='#')
	    {
	      if (strstr(loadLine,inFile->primary[0].telescope)!=NULL)
		{
		  sscanf(loadLine,"%s",observatoryDir);
		  break;
		}
	    }
	}	
    }
  fclose(fin);
  
  
  strcpy(telescope,inFile->primary[0].telescope);
  printf(" ... searching history log for observatory: %s\n",telescope);
  sprintf(fname,"%s/observatory/%s/historyLog/historyLog.dat",runtimeDir,observatoryDir);

  historyLog = (historyLogStruct *)malloc(sizeof(historyLogStruct)*MAX_HISTORY_LOG);

  if (!(fin = fopen(fname,"r")))
    {
      printf("UNABLE to open file %s\n",fname);
      return;
    }
  
  while (!feof(fin))
    {
      fgets(loadLine,1024,fin);
      if (loadLine[0]!='#' && strlen(loadLine) > 1)
	{
	  sscanf(loadLine,"%d",&historyLog[nLog].type);
	  if (historyLog[nLog].type==1)
	    {
	      tok = strtok(loadLine," ");
	      tok = strtok(NULL," ");
	      sscanf(tok,"%lf",&historyLog[nLog].mjdStart);
	      tok = strtok(NULL," ");
	      sscanf(tok,"%lf",&historyLog[nLog].mjdEnd);
	      tok = strtok(NULL," ");
	      strcpy(historyLog[nLog].detailedReport,tok);
	      tok = strtok(NULL,"\n");
	      strcpy(historyLog[nLog].report,tok);
	      nLog++;
	    }
	  else
	    printf("UNKNOWN TYPE %d in history log\n",historyLog[nLog].type);
	}
    }

  fclose(fin);

  for (i=0;i<nLog;i++)
    {
      mjd = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].mjd; // NOTE 0 HERE
      if (mjd >= historyLog[i].mjdStart && mjd  <= historyLog[i].mjdEnd)
	{
	  printf("         WARNING: %s ",historyLog[i].report);
	  if (strcmp(historyLog[i].detailedReport,"NULL")!=0)
	    printf("(for more details see %s/observatory/%s/historyLog/%s)\n",runtimeDir,observatoryDir,historyLog[i].detailedReport);
	  else
	    printf("\n");
	}
    }
  

  free(historyLog);
}

void checkVersion(char *ver)
{
  int primaryVer;
  sscanf(ver,"%d.",&primaryVer);
  if (primaryVer < 2)
    printf("         WARNING: This is an early version of the file format.\n");
}
