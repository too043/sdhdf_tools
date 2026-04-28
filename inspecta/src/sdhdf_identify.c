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
// Software to provide a quick-look at the metadata in a SDHDF file.
//
// Usage:
// sdhdf_identify <filename.hdf> <filename2.hdf> ...
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"
#include "hdf5.h"

#define MAX_SRC 128
#define MAX_BANDS 26

typedef struct infoStruct {
  char fname[MAX_STRLEN];
  char source[MAX_STRLEN];
  char projid[MAX_STRLEN];
  char hdrVersion[MAX_STRLEN];
  char observer[MAX_STRLEN];
  char rcvr[MAX_STRLEN];
  char telescope[MAX_STRLEN];
  int npol;
  char calMode[20];
  int scheduleNumber;
  float ra;
  float dec;
  double mjd;
  int nband;
  int ndump;
  int sb;
  
} infoStruct;

double haversine(double centre_long,double centre_lat,double src_long,double src_lat);


void help()
{
  printf("sdhdf_identify\n\n");
  printf("Command line options\n\n");
  printf("-allBand                        - look for matches amongst all sub-bands in the files (otherwise use band defined by -sb)\n");
  printf("-cal <id>                       - selects based on the noise source - usually ON or OFF\n");
  printf("-coord <ra0> <dec0> <raddist>   - selects file within a particular angular radius of those coordinates (RA = HH:MM:SS.SSS, DEC = DD:MM:SS.SSS, raddist = cone radius in degrees\n");
  printf("-h                              - this help\n");
  printf("-hdrversion <num>               - selects based on SDHDF header version\n");
  printf("-info                           - display more information output lines\n");
  printf("-npol <n>                       - selects based on the number of polarisation channels (1,2 or 4)\n");
  printf("-o <outFile>                    - write the output to an output file\n");
  printf("-observer <name>                - selects based on specified observer\n");
  printf("-pair                           - produce paired 'on' and 'off' output - need -on and -off set\n");
  printf("-partial                        - does not require an exact match when checking source names\n");
  printf("-rcvr <recevier>                - selects based on receiver name\n");
  printf("-sb <subband>                   - select sub-band number \n");
  printf("-sched <id>                     - selects based on the scheduling block identifier\n");
  printf("-src <name>                     - selects based on source name\n");
  printf("-tel <telescope>                - selects based on telescope name\n");
  printf("-nband <num>                    - selects based on the number of sub-bands in the file\n");
  printf("-ndump <num>                    - selects based on the number of spectral dumps in the file\n");
  printf("-on <source name>               - defines 'on' source to be the specified name\n");
  printf("-off <source name>              - defines 'off' source to be the specified name\n");
  printf("\n\n");
  printf("Example: sdhdf_identify -pair -on 0407-658 -off 0407-658_R *.hdf\n");
  printf("Example: sdhdf_identify -src 0407-658 *.hdf\n");
  printf("Example: sdhdf_identify -cal ON *.hdf\n");
  exit(1);
}

int main(int argc,char *argv[])
{
  int i,j,k,l;
  char fname[MAX_FILES][128];
  int nFiles=0;
  int ibeam=0;
  int iband=0;
  int allBand=0;
  infoStruct *info;
  sdhdf_fileStruct *inFile;
  int useCoord=0;
  char ra0Str[1024],dec0Str[1024];
  float ra0,dec0,raddist;
  double dist;
  int useSource=0;
  int infoDisp=0;
  char **src;
  char **onSrc;
  char **offSrc;
  int nSrc=0;
  int nOnSrc=0;
  int nOffSrc=0;
  int exactMatch=1;
  int openFile=0;
  int pair=0;
  char outFile[MAX_STRLEN];
  int useNband=0,selNband;
  int useProjID=0;
  char projID[MAX_STRLEN];
  int useHdrVersion=0;
  char hdrVersion[MAX_STRLEN];
  int useObserver=0;
  char observer[MAX_STRLEN];
  int useRcvr=0;
  char rcvr[MAX_STRLEN];
  int useTelescope=0;
  char telescope[MAX_STRLEN];
  int useCalMode=0;
  char calMode[MAX_STRLEN];
  int useNpol=0;
  int select_npol;
  int useSchedNumber=0;
  int select_sched;
  int useNdump=0;
  int select_ndump;

  
  int writeOut=0;
  FILE *fout;
  int disp;
  int nEntry=0;
  
  src = (char **)malloc(sizeof(char *)*MAX_SRC);
  onSrc = (char **)malloc(sizeof(char *)*MAX_SRC);
  offSrc = (char **)malloc(sizeof(char *)*MAX_SRC);
  for (i=0;i<MAX_SRC;i++)
    {
      src[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);
      onSrc[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);
      offSrc[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);
    }
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-coord")==0)
	{
	  useCoord=1;
	  sscanf(argv[++i],"%s",ra0Str);
	  sscanf(argv[++i],"%s",dec0Str);
	  sscanf(argv[++i],"%f",&raddist);
	}
      else if (strcasecmp(argv[i],"-allband")==0)
	allBand=1;
      else if (strcmp(argv[i],"-sb")==0)
	sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outFile,argv[++i]);
	  writeOut=1;
	}
      else if (strcmp(argv[i],"-pair")==0)
	pair=1;
      else if (strcmp(argv[i],"-partial")==0)
	exactMatch=0;
      else if (strcmp(argv[i],"-info")==0)
	infoDisp=1;
      else if (strcmp(argv[i],"-src")==0)
	{
	  useSource=1;
	  strcpy(src[nSrc++],argv[++i]);
	}
      else if (strcmp(argv[i],"-nband")==0)
	{
	  useNband = 1;
	  sscanf(argv[++i],"%d",&selNband);
	}
      else if (strcmp(argv[i],"-on")==0)
	strcpy(onSrc[nOnSrc++],argv[++i]);
      else if (strcmp(argv[i],"-off")==0)
	strcpy(offSrc[nOffSrc++],argv[++i]);
      else if (strcasecmp(argv[i],"-projid")==0)
	{  useProjID=1; strcpy(projID,argv[++i]); }
      else if (strcasecmp(argv[i],"-hdrversion")==0)
	{  useHdrVersion=1; strcpy(hdrVersion,argv[++i]); }
      else if (strcasecmp(argv[i],"-observer")==0)
	{  useObserver=1; strcpy(observer,argv[++i]); }
      else if (strcasecmp(argv[i],"-rcvr")==0)
	{  useRcvr=1; strcpy(rcvr,argv[++i]); }
      else if (strcasecmp(argv[i],"-tel")==0)
	{  useTelescope=1; strcpy(telescope,argv[++i]); }
      else if (strcasecmp(argv[i],"-cal")==0)
	{  useCalMode=1; strcpy(calMode,argv[++i]); }
      else if (strcasecmp(argv[i],"-npol")==0)
	{  useNpol=1; sscanf(argv[++i],"%d",&select_npol);}
      else if (strcasecmp(argv[i],"-sched")==0)
	{  useSchedNumber=1; sscanf(argv[++i],"%d",&select_sched);}
      else if (strcasecmp(argv[i],"-ndump")==0)
	{  useNdump=1; sscanf(argv[++i],"%d",&select_ndump);}
      else
	{
	  strcpy(fname[nFiles++],argv[i]);
	}
    }

  info = (infoStruct *)malloc(sizeof(infoStruct)*nFiles*MAX_BANDS); // SHOULD FIX THE MAX_BANDS HERE

  nEntry=0;
  for (i=0;i<nFiles;i++)
    {

      sdhdf_initialiseFile(inFile);
      openFile = sdhdf_openFile(fname[i],inFile,1);
      if (openFile==-1)
	{
	  printf("ERROR: Unable to open file %s\n",fname[i]);
	}
      else
	{
	  sdhdf_loadMetaData(inFile);
	  if (allBand == 0)
	    {
	      strcpy(info[nEntry].fname,fname[i]);
	      strcpy(info[nEntry].source,inFile->beamHeader[ibeam].source);
	      info[nEntry].ra = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].raDeg; // Note 0 here
	      info[nEntry].dec = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].decDeg;
	      info[nEntry].mjd = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].mjd + inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].fractional_mjd;
	      info[nEntry].nband = inFile->beam[ibeam].nBand;
	      info[nEntry].sb = iband;
	      strcpy(info[nEntry].projid,inFile->primary[0].pid);
	      strcpy(info[nEntry].hdrVersion,inFile->primary[0].hdr_defn_version);
	      strcpy(info[nEntry].observer,inFile->primary[0].observer);
	      strcpy(info[nEntry].rcvr,inFile->primary[0].rcvr);
	      strcpy(info[nEntry].calMode,inFile->primary[0].cal_mode);
	      info[nEntry].scheduleNumber = inFile->primary[0].sched_block_id;
	      info[nEntry].ndump = inFile->beam[ibeam].bandHeader[iband].ndump;
	      info[nEntry].npol = inFile->beam[ibeam].bandHeader[iband].npol;
	      
	      nEntry++;
	    }
	  else
	    {
	      for (j=0;j<inFile->beam[ibeam].nBand;j++)
		{
		  strcpy(info[nEntry].fname,fname[i]);
		  strcpy(info[nEntry].source,inFile->beamHeader[ibeam].source);
		  info[nEntry].ra = inFile->beam[ibeam].bandData[j].astro_obsHeader[0].raDeg; // Note 0 here
		  info[nEntry].dec = inFile->beam[ibeam].bandData[j].astro_obsHeader[0].decDeg;
		  info[nEntry].mjd = inFile->beam[ibeam].bandData[j].astro_obsHeader[0].mjd + inFile->beam[ibeam].bandData[j].astro_obsHeader[0].fractional_mjd;
		  info[nEntry].nband = inFile->beam[ibeam].nBand;
		  info[nEntry].sb = j;
		  strcpy(info[nEntry].projid,inFile->primary[0].pid);
		  strcpy(info[nEntry].hdrVersion,inFile->primary[0].hdr_defn_version);
		  strcpy(info[nEntry].observer,inFile->primary[0].observer);
		  strcpy(info[nEntry].rcvr,inFile->primary[0].rcvr);
		  strcpy(info[nEntry].calMode,inFile->primary[0].cal_mode);
		  info[nEntry].scheduleNumber = inFile->primary[0].sched_block_id;
		  info[nEntry].ndump = inFile->beam[ibeam].bandHeader[j].ndump;
		  info[nEntry].npol = inFile->beam[ibeam].bandHeader[j].npol;
		  nEntry++;
		}
	    }
	}
      sdhdf_closeFile(inFile);
    }
  //    printf("Loaded data for %d files\n",nFiles);


  if (writeOut==1)
    fout = fopen(outFile,"w");
    
  if (pair==1)
    {
      int foundOn=0;
      int foundOff=0;
      int bestOff;
      double timeDiff;
      for (i=0;i<nEntry;i++)
	{
	  // Find the on source
	  foundOn=0;
	  for (j=0;j<nOnSrc;j++)
	    {
	      if ((exactMatch==1 && strcmp(onSrc[j],info[i].source)==0) || (exactMatch==0 && strstr(info[i].source,onSrc[j])!=NULL))
		{foundOn=1; break;}
	    }
	  //	  printf("Found = %d %d\n",foundOn,exactMatch);
	  if (foundOn==1)
	    {
	      // Now find a matching off source
	      bestOff=-1;
	      for (j=0;j<nEntry;j++)
		{
		  foundOff=0;
		  for (k=0;k<nOffSrc;k++)
		    {
		      if ((exactMatch==1 && strcmp(offSrc[k],info[j].source)==0) || (exactMatch==0 && strstr(info[j].source,offSrc[k])!=NULL))
			{foundOff=1; break;}
		    }
		  if (foundOff==1)
		    {
		      if (bestOff == -1)
			{
			  bestOff = j;
			  timeDiff = fabs(info[j].mjd-info[i].mjd);
			}
		      else if (fabs(info[j].mjd-info[i].mjd) < timeDiff)
			{
			  bestOff = j;
			  timeDiff = fabs(info[j].mjd-info[i].mjd);
			}
		    }
		}
	      printf("[MATCH] %s %s # SRC_ON: %s SRC_OFF: %s TIMEDIFF: %g\n",info[i].fname,info[bestOff].fname,info[i].source,info[bestOff].source,timeDiff);
	      if (writeOut==1)
		fprintf(fout,"%s %s %s.%s.onoff # SRC_ON: %s SRC_OFF: %s TIMEDIFF: %g\n",info[i].fname,info[bestOff].fname,info[i].fname,info[bestOff].fname,info[i].source,info[bestOff].source,timeDiff);
	      //	      printf("Found %d as off source. Time difference = %g\n",bestOff,timeDiff);
	    }
	}
    }
  else
    {
      for (i=0;i<nEntry;i++)
	{
	  // Check source names
	  if (useSource==1)
	    {
	      for (j=0;j<nSrc;j++)
		{
		  disp=0;
		  if (exactMatch==1)
		    {
		      if (strcmp(info[i].source,src[j])==0) disp=1;			
		    }
		  else
		    {
		      if (strstr(info[i].source,src[j])!=NULL) disp=1;
		    }
		  if (disp==1)
		    {
		      if (infoDisp==1) printf("[SRC MATCH] %s %s sub-band %d %.3f\n",info[i].fname,info[i].source,info[i].sb,info[i].mjd);
		      else printf("%s\n",info[i].fname);
		    }
		}
	      
	    }
	  // Check coordinates
	  if (useCoord==1)
	    {	      
	      ra0 = turn_deg(hms_turn(ra0Str));
	      dec0 = turn_deg(dms_turn(dec0Str));
	      dist = haversine(ra0,dec0,info[i].ra,info[i].dec);
	      
	      if (dist < raddist)
		{
		  if (infoDisp==1) printf("[DIST MATCH] %s %g %g %g sub-band %d\n",info[i].fname,info[i].ra,info[i].dec,dist,info[i].sb);
		  else printf("%s\n",info[i].fname);		    
		}
	    }
	  // Check project ID
	  if (useProjID==1)
	    {
	      if (strcmp(info[i].projid,projID)==0)
		{
		  if (infoDisp==1) printf("[PROJID MATCH] %s %s sub-band %d\n",info[i].fname,info[i].projid,info[i].sb);
		  else printf("%s\n",info[i].fname);		    
		}
	    }
	  // Check HdrVersion
	  if (useHdrVersion==1)
	    {
	      if (strcmp(info[i].hdrVersion,hdrVersion)==0)
		{
		  if (infoDisp==1) printf("[HDR_VERSION MATCH] %s %s sub-band %d\n",info[i].fname,info[i].hdrVersion,info[i].sb);
		  else printf("%s\n",info[i].fname);		    
		}
	    }
	  // Check observer
	  if (useObserver==1)
	    {
	      if (strcmp(info[i].observer,observer)==0)
		{
		  if (infoDisp==1) printf("[OBSERVER MATCH] %s %s sub-band %d\n",info[i].fname,info[i].observer,info[i].sb);
		  else printf("%s\n",info[i].fname);		    
		}
	    }
	  // Check receiver
	  if (useRcvr==1)
	    {
	      if (strcmp(info[i].rcvr,rcvr)==0)
		{
		  if (infoDisp==1) printf("[RECEIVER MATCH] %s %s sub-band %d\n",info[i].fname,info[i].rcvr,info[i].sb);
		  else printf("%s\n",info[i].fname);		    
		}
	    }
	  // Check telescope
	  if (useTelescope==1)
	    {
	      if (strcmp(info[i].telescope,telescope)==0)
		{
		  if (infoDisp==1) printf("[TELESCOPE MATCH] %s %s sub-band %d\n",info[i].fname,info[i].telescope,info[i].sb);
		  else printf("%s\n",info[i].fname);		    
		}
	    }
	  // Check cal_mode
	  if (useCalMode==1)
	    {
	      if (strcmp(info[i].calMode,calMode)==0)
		{
		  if (infoDisp==1) printf("[CAL MATCH] %s %s sub-band %d\n",info[i].fname,info[i].calMode,info[i].sb);
		  else printf("%s\n",info[i].fname);		    
		}
	    }
	  // Check nband
	  if (useNband == 1)
	    {
	      if (info[i].nband == selNband)
		{
		  if (infoDisp==1) printf("[NBAND MATCH] %s sub-band %d\n",fname[i],info[i].sb);
		  else
		    printf("%s\n",fname[i]);
		}
	    }
	  // Check npol
	  if (useNpol == 1)
	    {
	      if (info[i].npol  == select_npol)
		{
		  if (infoDisp==1) printf("[NPOL MATCH] %s sub-band %d\n",fname[i],info[i].sb);
		  else
		    printf("%s\n",fname[i]);
		}
	    }
	  if (useSchedNumber == 1)
	    {
	      if (info[i].scheduleNumber  == select_sched)
		{
		  if (infoDisp==1) printf("[NPOL MATCH] %s sub-band %d\n",fname[i],info[i].sb);
		  else
		    printf("%s\n",fname[i]);
		}
	    }
	  if (useNdump == 1)
	    {
	      if (info[i].ndump  == select_ndump)
		{
		  if (infoDisp==1) printf("[NPOL MATCH] %s sub-band %d\n",fname[i],info[i].sb);
		  else
		    printf("%s\n",fname[i]);
		}
	    }
	}
    }
  if (writeOut==1)
    fclose(fout);

  free(inFile);
  free(info);
  for (i=0;i<MAX_SRC;i++)
    {
      free(src[i]);
      free(onSrc[i]);
      free(offSrc[i]);
    }
  free(src);
  free(onSrc);
  free(offSrc);

}



// Input in degrees
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


