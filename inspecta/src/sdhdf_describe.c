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
// sdhdf_describe <filename.hdf> <filename2.hdf> ...
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inspecta.h"

#define VERSION "v1.0"

void help()
{
  printf("sdhdf_describe %s (INSPECTA %s)\n",VERSION,SOFTWARE_VER);
  printf("Authors: G. Hobbs\n");
  printf("Purpose: to present meta-data information for multiple SDHDF files\n");
  printf("Example: sdhdf_describe file1.hdf file2.hdf\n");
  printf("\n");
  printf("Command line arguments:\n\n");
  printf("-atoa            Provide information useful for the ATOA\n");
  printf("-attributes      Show group and data set attributes\n");
  printf("-band            Provide band information\n");
  printf("-cband           Provide band information for the calibrator source\n");
  printf("-cdump           Provide spectral dump information for the calibrator source\n");
  printf("-dump            Provide spectral dump information\n");
  printf("-h               This help\n");
  printf("-history         Provide history information\n");
  printf("-sb <band>       Select band identifier (default = 0)\n");
  printf("-software        Provide software information\n");

  exit(1);

}


int main(int argc,char *argv[])
{
  int i,j,beam;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  float maxTime=0;
  float intTime=0;
  sdhdf_fileStruct *inFile;

  int showBands=0;
  int showDump=0;
  int showHistory=0;
  int showSoftware=0;
  int showCalBands=0;
  int showCalDump=0;
  int showAttributes=0;
  int atoa=0;
  int calInfo=0;
  
  int iband=0;
  int test=0;
  
  // Display help if no commands given
  if (argc==1)
    help();
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-h")==0)
	help();      
      else if (strcmp(argv[i],"-atoa")==0)
	atoa=1;
      else if (strcmp(argv[i],"-band")==0)
	showBands=1;
      else if (strcmp(argv[i],"-test")==0)
	test=1;
      else if (strcmp(argv[i],"-cband")==0)
	showCalBands=1;
      else if (strcmp(argv[i],"-attributes")==0)
	showAttributes=1;
      else if (strcmp(argv[i],"-sb")==0)
	sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-dump")==0)
	showDump=1;
      else if (strcmp(argv[i],"-cal")==0)
	calInfo=1;
      else if (strcmp(argv[i],"-cdump")==0)
	showCalDump=1;
      else if (strcmp(argv[i],"-history")==0)
	showHistory=1;
      else if (strcmp(argv[i],"-software")==0)
	showSoftware=1;
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

  if (atoa==0)
    {
      printf("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
      printf("%-36.36s %-5.5s %-22.22s %-5.5s %-6.6s %-8.8s %-20.20s %-10.10s %-10.10s %-7.7s %-6.6s %-10.10s %-5.5s %-6.6s %-12.12s %-12.12s\n","File","Beam","UTC","SDHDF","PID","Sched_ID","Source","Tel","Observer","RCVR","Bands","POL_TYPE","CAL","M_Time","RA_s","DEC_s");
      printf("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    }
  
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      if (sdhdf_openFile(fname[i],inFile,1)==-1)
	printf("Warning: unable to open file >%s<. Skipping\n",fname[i]);
      else
	{
	  sdhdf_loadMetaData(inFile);
	  if (atoa==1)
	    {
	      maxTime=0;
	      //
	      // Note here that we are hardcoding to beam 0
	      //
	      for (j=0;j<inFile->beam[0].nBand;j++)
		{
		  intTime = inFile->beam[0].bandHeader[j].dtime*inFile->beam[0].bandHeader[j].ndump;
		  if (maxTime < intTime)
		    maxTime = intTime;
		}
	      printf("%s %s %s %s %s %s %8.3f\n",inFile->fname,inFile->primary[0].pid,inFile->beamHeader[0].source,inFile->beam[0].bandData[0].astro_obsHeader[0].raStr,inFile->beam[0].bandData[0].astro_obsHeader[0].decStr,inFile->primary[0].utc0,maxTime);

              for (j=0;j<inFile->beam[0].nBand;j++)
		{
		  intTime = inFile->beam[0].bandHeader[j].dtime*inFile->beam[0].bandHeader[j].ndump;
		  printf("%d %-10.10s %d %d %d %d %d %g\n",j,inFile->beam[0].bandHeader[j].label,inFile->beam[0].bandHeader[j].nchan,(int)inFile->beam[0].bandHeader[j].f0,(int)inFile->beam[0].bandHeader[j].f1,(int)inFile->beam[0].bandHeader[j].fc,(int)(inFile->beam[0].bandHeader[j].f1-inFile->beam[0].bandHeader[j].f0),inFile->beam[0].bandHeader[j].dtime);
		}

	    }
	  else
	    {
	      for (beam=0;beam<inFile->nBeam;beam++)
		{
		  maxTime=0;  
		  for (j=0;j<inFile->beam[beam].nBand;j++)
		    {
		      intTime = inFile->beam[beam].bandHeader[j].dtime*inFile->beam[beam].bandHeader[j].ndump;
		      if (maxTime < intTime)
			maxTime = intTime;
		    }
		  
		  printf("%-36.36s %-5d %-22.22s %-5.5s %-6.6s %-8d %-20.20s %-10.10s %-10.10s %-7.7s %-6d %-10.10s %-5.5s %-6.1f %s %s\n",
			 inFile->fname, beam,inFile->primary[0].utc0,inFile->primary[0].hdr_defn_version,
			 inFile->primary[0].pid,(int)inFile->primary[0].sched_block_id,inFile->beamHeader[iband].source,
			 inFile->primary[0].telescope,inFile->primary[0].observer,inFile->primary[0].rcvr,inFile->beam[beam].nBand,
			 inFile->beam[beam].bandHeader[0].pol_type,inFile->primary[0].cal_mode,maxTime,inFile->beam[beam].bandData[iband].astro_obsHeader[0].raStr,inFile->beam[beam].bandData[iband].astro_obsHeader[0].decStr);

		  if (test==1)
		    {
		      int na,kk,ndims,ii;
		      hsize_t     size[32];
		      hsize_t     maxsize[32];
		      hid_t attr_id;
		      hid_t type, space, type_class,space_type,ptype;
		      char attrName[1024];
		      hsize_t     alloc_size;
		      hsize_t     nelmts = 1;
		      void       *buf = NULL;
		      int nread;
		      printf("IN TESTING ATTRIBUTES\n");
		      na = sdhdf_getNattributes(inFile,"/beam_0/band_SB0/astronomy_data/data");
		      printf("na = %d\n",na);
		      for (kk=0;kk<na;kk++)
			{
			  nelmts = 1;
			  printf("Attribute %d\n",kk);
			  attr_id = H5Aopen_by_idx(inFile->fileID,"/beam_0/band_SB0/astronomy_data/data", H5_INDEX_NAME, H5_ITER_NATIVE,kk,H5P_DEFAULT,H5P_DEFAULT);
			  H5Aget_name(attr_id,1024,attrName);
			  type = H5Aget_type(attr_id);

			  type_class   = H5Tget_class(type);
			  printf("attr_id = %d, type = %d, type_class = %d > %s <\n",(int)attr_id,(int)type,(int)type_class,attrName);
			  if (type_class == H5T_STRING) printf("STRING\n");
			  else if (type_class == H5T_FLOAT) printf("FLOAT\n");
			  else if (type_class == H5T_INTEGER) printf("INTEGER\n");
			  else printf("NO IDEA WHAT TYPE_CLASS THIS IS\n");

			  space = H5Aget_space(attr_id);
			  space_type = H5Sget_simple_extent_type(space);
			  printf("space = %d\n",(int)space);
			  if (space_type == H5S_SCALAR) printf("SPACE: H5S_SCALAR\n");
			  else if (space_type == H5S_SIMPLE) printf("SPACE: SIMPLE\n");
			  else if (space_type == H5S_NULL) printf("SPACE: NULL\n");
			  else printf("SPACE: UNKNOWN\n");
			  
			  ndims = H5Sget_simple_extent_dims(space, size, maxsize);
			  printf("ndims = %d\n",ndims);
			  for (ii=0;ii<ndims;ii++)
			    {
			      nelmts *= size[ii];
			      if (maxsize[ii] == H5S_UNLIMITED)
				printf("dims %d, size = %d, maxsize = %d (UNLIMITED)\n",ii,(int)size[ii],(int)maxsize[ii]);
			      else
				printf("dims %d, size = %d, maxsize = %d\n",ii,(int)size[ii],(int)maxsize[ii]);
			    }
			  printf("nelmts = %d\n",(int)nelmts);
			  alloc_size = nelmts * H5Tget_size(type);
			  printf("alloc_size = %d\n",(int)alloc_size);
			  buf =  malloc((size_t)alloc_size);
			  nread = H5Aread(attr_id, type, buf);
			  printf("READ %d\n",(int)nread);

			  // TO PRINT IT -- SEE h5tools_dump_simple_data
			}
		      
		    }
		  if (calInfo==1)
		    {
		      printf("[CAL] %s %s %d %s %s %f %f %f\n",inFile->fname,inFile->beamHeader[iband].source,(int)inFile->primary[0].sched_block_id,inFile->primary[0].cal_mode,inFile->cal_epoch,inFile->cal_freq,inFile->cal_phase,inFile->cal_duty_cycle);
		    }
		  
		  if (showBands==1)
		    {
		      printf("-----------------------------------------------------------------------------\n");
		      printf("        #   BandID           F0       F1      NCHAN    TDUMP    NPOL NDUMP TOBS\n");
		      printf("-----------------------------------------------------------------------------\n");
		      for (j=0;j<inFile->beam[beam].nBand;j++)
			{
			  intTime = inFile->beam[beam].bandHeader[j].dtime*inFile->beam[beam].bandHeader[j].ndump;
			  printf(" [Band] %3.3d %-15.15s %8.2f %8.2f %-8d %-8.3f %-4d %-5d %-8.3f %s\n",j,inFile->beam[beam].bandHeader[j].label,inFile->beam[beam].bandHeader[j].f0,inFile->beam[beam].bandHeader[j].f1,inFile->beam[beam].bandHeader[j].nchan,inFile->beam[beam].bandHeader[j].dtime,inFile->beam[beam].bandHeader[j].npol,inFile->beam[beam].bandHeader[j].ndump,intTime,inFile->fname);
			}
		    }
		  
		  // Showing dumps for band 0
		  if (showDump==1)
		    {
		      printf("------------------------------------------------------------------------------------------------------------------------------------\n");
		      printf("          ElapsedTime MJD           UTC         LOCAL        RA          DEC          RADEG   DECDEG AZ       EL     GL      GB     PA \n");
		      printf("------------------------------------------------------------------------------------------------------------------------------------\n");
		      for (j=0;j<inFile->beam[beam].bandData[iband].nAstro_obsHeader;j++)
			{
			  printf(" [SDUMP] %12.3f %.7f %-11.11s %-11.11s %-11.11s %-11.11s %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %s %d %d\n",
				 inFile->beam[beam].bandData[iband].astro_obsHeader[j].timeElapsed,inFile->beam[beam].bandData[iband].astro_obsHeader[j].mjd,
				 inFile->beam[beam].bandData[iband].astro_obsHeader[j].utc,inFile->beam[beam].bandData[iband].astro_obsHeader[j].local_time,
				 inFile->beam[beam].bandData[iband].astro_obsHeader[j].raStr,inFile->beam[beam].bandData[iband].astro_obsHeader[j].decStr,
				 inFile->beam[beam].bandData[iband].astro_obsHeader[j].raDeg,inFile->beam[beam].bandData[iband].astro_obsHeader[j].decDeg,
				 inFile->beam[beam].bandData[iband].astro_obsHeader[j].az,inFile->beam[beam].bandData[iband].astro_obsHeader[j].el,
				 inFile->beam[beam].bandData[iband].astro_obsHeader[j].gl,inFile->beam[beam].bandData[iband].astro_obsHeader[j].gb,
				 inFile->beam[beam].bandData[iband].astro_obsHeader[j].paraAngle,inFile->fname,j,beam);
			}
		    }
		
		  if (showCalBands==1)
		    {
		      printf("-----------------------------------------------------------------------------\n");
		      printf("        #   BandID      F0       F1      NCHAN    TDUMP    NPOL NDUMP TOBS\n");
		      printf("-----------------------------------------------------------------------------\n");
		      for (j=0;j<inFile->beam[beam].nBand;j++)
			{
			  intTime = inFile->beam[beam].calBandHeader[j].dtime*inFile->beam[beam].calBandHeader[j].ndump;
			  printf(" [Band] %3.3d %-10.10s %8.2f %8.2f %-8d %-8.3f %-4d %-5d %-8.3f\n",j,inFile->beam[beam].calBandHeader[j].label,inFile->beam[beam].calBandHeader[j].f0,inFile->beam[beam].calBandHeader[j].f1,inFile->beam[beam].calBandHeader[j].nchan,inFile->beam[beam].calBandHeader[j].dtime,inFile->beam[beam].calBandHeader[j].npol,inFile->beam[beam].calBandHeader[j].ndump,intTime);
			}
		    }
		  
		  // Showing dumps for band 0
		  if (showCalDump==1)
		    {
		      printf("------------------------------------------------------------------------------------------------------------------------------------\n");
		      printf("          ElapsedTime MJD           UTC         LOCAL        RA          DEC          RADEG   DECDEG AZ       EL     GL      GB\n");
		      printf("------------------------------------------------------------------------------------------------------------------------------------\n");
		      for (j=0;j<inFile->beam[beam].bandData[iband].nCal_obsHeader;j++)
			{
			  printf(" [SDUMP] %12.3f %.7f %-11.11s %-11.11s %-11.11s %-11.11s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
				 inFile->beam[beam].bandData[iband].cal_obsHeader[j].timeElapsed,inFile->beam[beam].bandData[iband].cal_obsHeader[j].mjd,
				 inFile->beam[beam].bandData[iband].cal_obsHeader[j].utc,inFile->beam[beam].bandData[iband].cal_obsHeader[j].local_time,
				 inFile->beam[beam].bandData[iband].cal_obsHeader[j].raStr,inFile->beam[beam].bandData[iband].cal_obsHeader[j].decStr,
				 inFile->beam[beam].bandData[iband].cal_obsHeader[j].raDeg,inFile->beam[beam].bandData[iband].cal_obsHeader[j].decDeg,
				 inFile->beam[beam].bandData[iband].cal_obsHeader[j].az,inFile->beam[beam].bandData[iband].cal_obsHeader[j].el,
				 inFile->beam[beam].bandData[iband].cal_obsHeader[j].gl,inFile->beam[beam].bandData[iband].cal_obsHeader[j].gb);
			  
			}
		    }



		}
	    }

	  if (showHistory==1)
	    {
	      printf("------------------------------------------------------------------------------------------------------------------------------------\n");
	      printf("DATE                  PROC               PROC_DESCR                                   PROC_ARGS          PROC_HOST\n");
	      printf("------------------------------------------------------------------------------------------------------------------------------------\n");

	      for (j=0;j<inFile->nHistory;j++)
		{
		  printf("%-20.20s | %-16.16s | %-42.42s | %s | %s\n",inFile->history[j].date,inFile->history[j].proc_name,inFile->history[j].proc_descr,inFile->history[j].proc_args,inFile->history[j].proc_host);
		}
	    }

	  if (showSoftware==1)
	    {
	      printf("------------------------------------------------------------------------------------------------------------------------------------\n");
	      printf("PROC                  SOFTWARE             DESCRIPTION                                  VERSION\n");
	      printf("------------------------------------------------------------------------------------------------------------------------------------\n");
	      for (j=0;j<inFile->nSoftware;j++)
		{
		  printf("%-20.20s | %-18.18s | %-42.42s | %-24.24s \n",inFile->software[j].proc_name,inFile->software[j].software,inFile->software[j].software_descr,inFile->software[j].software_version);
		}
	    }
	  if (showAttributes==1)
	    {
	      printf("\n\n");
	      printf("Attributes\n\n");
	      printf("Primary_header\n\n");
	      for (j=0;j<inFile->nPrimaryAttributes;j++)
		printf("%-20.20s %s\n",inFile->primaryAttr[j].key,inFile->primaryAttr[j].value);

	      printf("\n\nBeam_header\n\n");
	      for (j=0;j<inFile->nBeamHeaderAttributes;j++)
		printf("%-20.20s %s\n",inFile->beamHeaderAttr[j].key,inFile->beamHeaderAttr[j].value);
	    }
      
	  sdhdf_closeFile(inFile);
	}
    }

  free(inFile);
}
