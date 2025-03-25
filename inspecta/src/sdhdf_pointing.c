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
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "hdf5.h"
#include "inspecta.h"
#include <cpgplot.h>
#include "TKnonlinearFit.h"

#define MAX_DUMPS 4096
#define MAX_BANDS 26

typedef struct modelFitStruct {
  double amp;
  double offset;
  double width;
  double baseline; 
  double minx;
  double maxx;
} modelFitStruct;

typedef struct sdhdf_scanResult {
  float x[MAX_DUMPS];
  float pol1[MAX_DUMPS];
  float pol2[MAX_DUMPS];
} sdhdf_scanResult;


double haversine(double centre_long,double centre_lat,double src_long,double src_lat);
void processScan(sdhdf_fileStruct *inFile,int ibeam,int iband,float *fy,int *nPlot,int pol);
void sdhdf_calc1MinMax(float *v,int n,float *min,float *max);
void fitBeam(sdhdf_scanResult scan,int npts,int pol,modelFitStruct *model);
void plotModel(modelFitStruct model,float minx,float maxx);
void initialiseModel(modelFitStruct *pol1,modelFitStruct *pol2, sdhdf_fileStruct *file,int ibeam,int iband);

int main(int argc,char *argv[])
{
  int i,j,k,band,pol,scan,zoom;
  sdhdf_fileStruct *ra_f,*ra_b,*dec_f,*dec_b;

  char fname_raf[MAX_STRLEN];
  char fname_rab[MAX_STRLEN];
  char fname_decf[MAX_STRLEN];
  char fname_decb[MAX_STRLEN];

  sdhdf_scanResult result_raf;
  sdhdf_scanResult result_rab;
  sdhdf_scanResult result_decf;
  sdhdf_scanResult result_decb;
  sdhdf_scanResult result_ra;
  sdhdf_scanResult result_dec;
  
  modelFitStruct model_raf_p1;
  modelFitStruct model_rab_p1;
  modelFitStruct model_decf_p1;
  modelFitStruct model_decb_p1;

  modelFitStruct model_raf_p2;
  modelFitStruct model_rab_p2;
  modelFitStruct model_decf_p2;
  modelFitStruct model_decb_p2;
  
  int   nPlot;
  float miny,maxy;
  float mx,my;
  char key;
  float wx1,wx2,wy1,wy2;
  int ibeam = 0;
  char xBox[MAX_STRLEN];
  char yBox[MAX_STRLEN];
  int col;
  
  double ra0,dec0;
  double ra1,dec1;
  double angle;
  float zx[2],zy[2];
  char title[1024];
  char text[1024];

  int interactive=0;

  float freq[MAX_BANDS];
  float raOff_p1_fwd[MAX_BANDS],raOff_p2_fwd[MAX_BANDS];
  float raOff_p1_bwd[MAX_BANDS],raOff_p2_bwd[MAX_BANDS];
  float decOff_p1_fwd[MAX_BANDS],decOff_p2_fwd[MAX_BANDS];
  float decOff_p1_bwd[MAX_BANDS],decOff_p2_bwd[MAX_BANDS];
  float beamWidth_ra_p1_fwd[MAX_BANDS],beamWidth_ra_p2_fwd[MAX_BANDS];
  float beamWidth_dec_p1_fwd[MAX_BANDS],beamWidth_dec_p2_fwd[MAX_BANDS];
  float beamWidth_ra_p1_bwd[MAX_BANDS],beamWidth_ra_p2_bwd[MAX_BANDS];
  float beamWidth_dec_p1_bwd[MAX_BANDS],beamWidth_dec_p2_bwd[MAX_BANDS];

  float raOff_p1_avg[MAX_BANDS],raOff_p2_avg[MAX_BANDS];
  float decOff_p1_avg[MAX_BANDS],decOff_p2_avg[MAX_BANDS];
  float raOff_avg[MAX_BANDS];
  float decOff_avg[MAX_BANDS];

  float avAntEff[MAX_BANDS];
  float gain_dpfu[MAX_BANDS];

  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-beam")==0)
	sscanf(argv[++i],"%d",&ibeam);
      else if (strcasecmp(argv[i],"-raf")==0)
	strcpy(fname_raf,argv[++i]);
      else if (strcasecmp(argv[i],"-rab")==0)
	strcpy(fname_rab,argv[++i]);
      else if (strcasecmp(argv[i],"-decf")==0)
	strcpy(fname_decf,argv[++i]);
      else if (strcasecmp(argv[i],"-decb")==0)
	strcpy(fname_decb,argv[++i]);     
      else if (strcasecmp(argv[i],"-i")==0)
	interactive=1;
    }  

  ra_f = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));
  ra_b = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));
  dec_f = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));
  dec_b = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));

  sdhdf_initialiseFile(ra_f);
  sdhdf_initialiseFile(ra_b);
  sdhdf_initialiseFile(dec_f);
  sdhdf_initialiseFile(dec_b);

  if (sdhdf_openFile(fname_raf,ra_f,1)==-1)
    {printf("Error: unable to open file >%s<.\n",fname_raf); exit(1);}
  if (sdhdf_openFile(fname_rab,ra_b,1)==-1)
    {printf("Error: unable to open file >%s<.\n",fname_rab); exit(1);}
  if (sdhdf_openFile(fname_decf,dec_f,1)==-1)
    {printf("Error: unable to open file >%s<.\n",fname_decf); exit(1);}
  if (sdhdf_openFile(fname_decb,dec_b,1)==-1)
    {printf("Error: unable to open file >%s<.\n",fname_decb); exit(1);}

  sdhdf_loadMetaData(ra_f);
  sdhdf_loadMetaData(ra_b);
  sdhdf_loadMetaData(dec_f);
  sdhdf_loadMetaData(dec_b);

		  
  // FIX ME
  printf("Setting default position to 0407-657\n");
  //      04 08 20.37884 -65 45 09.0806
  ra0 = (4.+8./60.0+20.37884/60./60.)*180./12.;
  dec0 = -65.0-45.0/60.-9.0806/60./60.;
  
  // FIX ME: should check that all the files have the same numbers of bands etc.

  if (interactive==1)
    cpgbeg(0,"/xs",1,1);
  else
    cpgbeg(0,"sdhdf_pointing.ps/cps",1,1);
  cpgask(0);  cpgsch(1.4);  cpgslw(2);  cpgscf(2);
  for (band = 0;band < ra_f->beam[ibeam].nBand;band++)
    {
      freq[band] = ra_f->beam[ibeam].bandHeader[band].fc;

      initialiseModel(&model_raf_p1,&model_raf_p2,ra_f,ibeam,band);
      initialiseModel(&model_rab_p1,&model_rab_p2,ra_b,ibeam,band);
      initialiseModel(&model_decf_p1,&model_decf_p2,dec_f,ibeam,band);
      initialiseModel(&model_decb_p1,&model_decb_p2,dec_b,ibeam,band);

      if (interactive==1)
	cpgeras();
      else
	cpgpage();
      
      cpgsch(1.4);
      cpgsvp(0.0,1.0,0.15,1);      
      cpgswin(0,1,0,1);
      cpgbox("",0,0,"",0,0);

      cpglab("Angular offset (arc seconds)","",title);
      cpgsch(0.8);
      sprintf(text,"Input files:"); cpgtext(0.02,0.97,text);
      sprintf(text,"RA  +: %s",fname_raf); cpgtext(0.02,0.94,text);
      sprintf(text,"RA  -: %s",fname_rab); cpgtext(0.02,0.91,text);
      sprintf(text,"DEC +: %s",fname_decf); cpgtext(0.02,0.88,text);
      sprintf(text,"RA  -: %s",fname_decb); cpgtext(0.02,0.85,text);

      sprintf(text,"%s: %.1f-%.1f MHz",ra_f->beam[ibeam].bandHeader[band].label,ra_f->beam[ibeam].bandHeader[band].f0,ra_f->beam[ibeam].bandHeader[band].f1);
      cpgsci(2); cpgtext(0.5,0.97,text); cpgsci(1);

      cpgsch(1.4);

      
      
      for (scan=0;scan < 4; scan++)
	{
	  if (scan==0)
	    {
	      for (i=0;i<ra_f->beam[ibeam].bandHeader[band].ndump;i++)
		{
		  ra1 = ra_f->beam[ibeam].bandData[band].astro_obsHeader[i].raDeg;
		  result_raf.x[i] = (ra1-ra0)*cos(dec0*M_PI/180.0)*60.*60.;
		}
	      wy1 = 0.15; wy2 = 0.475;
	      col=1;
	    }
	  else if (scan==1)
	    {
	      for (i=0;i<ra_b->beam[ibeam].bandHeader[band].ndump;i++)
		{
		  ra1 = ra_b->beam[ibeam].bandData[band].astro_obsHeader[i].raDeg;
		  angle = (ra1-ra0)*cos(dec0*M_PI/180.0)*60.*60.;
		  result_rab.x[i] = angle;
		}
	      wy1 = 0.15; wy2 = 0.475;
	      col=5;
	    }
	  else if (scan==2)
	    {
	      for (i=0;i<dec_f->beam[ibeam].bandHeader[band].ndump;i++)
		{
		  dec1 = dec_f->beam[ibeam].bandData[band].astro_obsHeader[i].decDeg;
		  
		  angle = (dec1-dec0)*60.*60.;
		  result_decf.x[i] = angle;
		}
	      wy1 = 0.475; wy2 = 0.8;
	      col=1;
	    }
	  else if (scan==3)
	    {
	      for (i=0;i<dec_b->beam[ibeam].bandHeader[band].ndump;i++)
		{
		  dec1 = dec_b->beam[ibeam].bandData[band].astro_obsHeader[i].decDeg;
		  
		  angle = (dec1-dec0)*60.*60.;
		  result_decb.x[i] = angle;
		}
	      wy1 = 0.475; wy2 = 0.8;
	      col=5;
	    }
	  for (pol=0;pol<2;pol++)
	    {
	      if (pol==0)
		{
		  if (scan==0)
		    {
		      processScan(ra_f,ibeam,band,result_raf.pol1,&nPlot,pol+1);
		      sdhdf_calc1MinMax(result_raf.pol1,nPlot,&miny,&maxy);
		      for (k=0;k<nPlot;k++)
			result_ra.pol1[k]  = result_raf.pol1[k];
		      
		      fitBeam(result_raf,nPlot,1,&model_raf_p1);
		      
		      raOff_p1_fwd[band] = model_raf_p1.offset;
		      beamWidth_ra_p1_fwd[band] = model_raf_p1.width*2.355/60.;

		    }
		  else if (scan==1)
		    {
		      processScan(ra_b,ibeam,band,result_rab.pol1,&nPlot,pol+1);
		      sdhdf_calc1MinMax(result_rab.pol1,nPlot,&miny,&maxy);
		      fitBeam(result_rab,nPlot,1,&model_rab_p1);



		      
		      raOff_p1_bwd[band] = model_rab_p1.offset;
		      beamWidth_ra_p1_bwd[band] = model_rab_p1.width*2.355/60.;
		    }
		  else if (scan==2)
		    {
		      processScan(dec_f,ibeam,band,result_decf.pol1,&nPlot,pol+1);
		      sdhdf_calc1MinMax(result_decf.pol1,nPlot,&miny,&maxy);
		      fitBeam(result_decf,nPlot,1,&model_decf_p1);

 		      for (k=0;k<nPlot;k++)
			result_dec.pol1[k]  = result_decf.pol1[k];



		      
		      decOff_p1_fwd[band]        = model_decf_p1.offset;
		      beamWidth_dec_p1_fwd[band] = model_decf_p1.width*2.355/60.;

		    }
		  else if (scan==3)
		    {
		      processScan(dec_b,ibeam,band,result_decb.pol1,&nPlot,pol+1);
		      sdhdf_calc1MinMax(result_decb.pol1,nPlot,&miny,&maxy);
		      fitBeam(result_decb,nPlot,1,&model_decb_p1);


		      
		      decOff_p1_bwd[band]        = model_decb_p1.offset;
		      beamWidth_dec_p1_bwd[band] = model_decb_p1.width*2.355/60.;
		    }

		  
		  wx1 = 0.1;
		  wx2 = 0.3833;
		}
	      else if (pol==1)
		{
		  if (scan==0)
		    {
		      processScan(ra_f,ibeam,band,result_raf.pol2,&nPlot,pol+1);
		      sdhdf_calc1MinMax(result_raf.pol2,nPlot,&miny,&maxy);
		      fitBeam(result_raf,nPlot,2,&model_raf_p2);
		      raOff_p2_fwd[band] = model_raf_p2.offset;
		      beamWidth_ra_p2_fwd[band] = model_raf_p2.width*2.355/60.;
 		      for (k=0;k<nPlot;k++)
			result_ra.pol2[k]  = result_raf.pol2[k];
		      
		    }
		  else if (scan==1)
		    {
		      processScan(ra_b,ibeam,band,result_rab.pol2,&nPlot,pol+1);
		      sdhdf_calc1MinMax(result_rab.pol2,nPlot,&miny,&maxy);
		      fitBeam(result_rab,nPlot,2,&model_rab_p2);
		      raOff_p2_bwd[band] = model_rab_p2.offset;
		      beamWidth_ra_p2_bwd[band] = model_rab_p2.width*2.355/60.;

		    }
		  else if (scan==2)
		    {
		      processScan(dec_f,ibeam,band,result_decf.pol2,&nPlot,pol+1);
		      sdhdf_calc1MinMax(result_decf.pol2,nPlot,&miny,&maxy);
		      fitBeam(result_decf,nPlot,2,&model_decf_p2);
		      decOff_p2_fwd[band] = model_decf_p2.offset;
		      beamWidth_dec_p2_fwd[band] = model_decf_p2.width*2.355/60.;
 		      for (k=0;k<nPlot;k++)
			result_dec.pol2[k]  = result_decf.pol2[k];

		    }
		  else if (scan==3)
		    {
		      processScan(dec_b,ibeam,band,result_decb.pol2,&nPlot,pol+1);
		      sdhdf_calc1MinMax(result_decb.pol2,nPlot,&miny,&maxy);
		      fitBeam(result_decb,nPlot,2,&model_decb_p2);
		      decOff_p2_bwd[band] = model_decb_p2.offset;
		      beamWidth_dec_p2_bwd[band] = model_decb_p2.width*2.355/60.;

		    }

 		  wx1 = 0.3833;
		  wx2 = 0.6666;
		}

	      
	      for (zoom=0;zoom<2;zoom++)
		{
		  if (zoom==0)
		    {
		      cpgsch(1.4);
		      cpgsvp(wx1,wx2,wy1,wy2);      
		      cpgswin(-900,900.,miny,maxy+(maxy-miny)*0.1);
		    }
		  else if (zoom==1)
		    {
		      cpgsch(0.8);
		      cpgsvp(wx1+(wx2-wx1)*0.65,wx2-(wx2-wx1)*0.05,wy1+(wy2-wy1)*0.6,wy2-(wy2-wy1)*0.1);      
		      cpgswin(-1.5*60.,1.5*60.,0.9,1.02);		      
		    }
		  if (zoom==1)
		    strcpy(yBox,"BC");
		  else
		    {
		      if (pol == 0)
			strcpy(yBox,"BCTSN");
		      else
			strcpy(yBox,"BCTS");
		    }
		  if (scan==0 || zoom==1)
		    strcpy(xBox,"BCTSN");
		  else
		    strcpy(xBox,"BCTS");
		  cpgbox(xBox,0,0,yBox,0,0);
		  
		  zx[0] = zx[1] = 0;
		  zy[0] = miny; zy[1] = 1.1;
		  cpgsls(4); cpgline(2,zx,zy); cpgsls(1);
		  cpgsci(col);
		  if (pol==0)
		    {
		      if (scan==0)
			{
			  plotModel(model_raf_p1,-900,900);
			  cpgpt(nPlot,result_raf.x,result_raf.pol1,17);
			  if (zoom==0){
			    cpgsch(0.8); cpgtext(-800,0.9,"RA(+,1)"); cpgsch(1.0);
			  }
			}
		      else if (scan==1)
			{
			  plotModel(model_rab_p1,-900,900);
			  cpgpt(nPlot,result_rab.x,result_rab.pol1,17);
			  if (zoom==0){
			    cpgsch(0.8); cpgtext(-800,0.8,"RA(-,1)"); cpgsch(1.0);
			  }
			}
		      else if (scan==2)
			{
			  plotModel(model_decf_p1,-900,900);
			  cpgpt(nPlot,result_decf.x,result_decf.pol1,17);
			  if (zoom==0){
			    cpgsch(0.8); cpgtext(-800,0.9,"DEC(+,1)"); cpgsch(1.0);
			  }
			}
		      else if (scan==3)
			{
			  plotModel(model_decb_p1,-900,900);
			  cpgpt(nPlot,result_decb.x,result_decb.pol1,17);
			  if (zoom==0){
			    cpgsch(0.8); cpgtext(-800,0.8,"DEC(-,1)"); cpgsch(1.0);
			  }
			}
		    }
		  else if (pol==1)
		    {
		      if (scan==0)
			{
			  plotModel(model_raf_p2,-900,900);
			  cpgpt(nPlot,result_raf.x,result_raf.pol2,17);
			  if (zoom==0){
			    cpgsch(0.8); cpgtext(-800,0.9,"RA(+,2)"); cpgsch(1.0);
			  }
			}
		      else if (scan==1)
			{
			  plotModel(model_rab_p2,-900,900);
			  cpgpt(nPlot,result_rab.x,result_rab.pol2,17);
			  if (zoom==0){
			    cpgsch(0.8); cpgtext(-800,0.8,"RA(-,2)"); cpgsch(1.0);
			  }
			}
		      else if (scan==2)
			{
			  plotModel(model_decf_p2,-900,900);
			  cpgpt(nPlot,result_decf.x,result_decf.pol2,17);
			  if (zoom==0){
			    cpgsch(0.8); cpgtext(-800,0.9,"DEC(+,2)"); cpgsch(1.0);
			  }
			}
		      else if (scan==3)
			{
			  plotModel(model_decb_p2,-900,900);
			  cpgpt(nPlot,result_decb.x,result_decb.pol2,17);
			  if (zoom==0){
			    cpgsch(0.8); cpgtext(-800,0.8,"DEC(-,2)"); cpgsch(1.0);
			  }
			}
			  
		    }
		  cpgsci(1);
		}

	    }
	}
      cpgsvp(0.0,1.0,0.15,1);      
      cpgswin(0,1,0,1);

      sprintf(text,"RA(+,1) offset: %.2f arc sec",model_raf_p1.offset);
      cpgsch(0.8); cpgtext(0.35,0.94,text); cpgsch(1.4);
      sprintf(text,"RA(-,1) offset: %.2f arc sec",model_rab_p1.offset);
      cpgsch(0.8); cpgtext(0.35,0.91,text); cpgsch(1.4);
      sprintf(text,"DEC(+,1) offset: %.2f arc sec",model_decf_p1.offset);
      cpgsch(0.8); cpgtext(0.35,0.88,text); cpgsch(1.4);
      sprintf(text,"DEC(-,1) offset: %.2f arc sec",model_decb_p1.offset);
      cpgsch(0.8); cpgtext(0.35,0.85,text); cpgsch(1.4);

      sprintf(text,"RA(+,2) offset: %.2f arc sec",model_raf_p2.offset);
      cpgsch(0.8); cpgtext(0.65,0.94,text); cpgsch(1.4);
      sprintf(text,"RA(-,2) offset: %.2f arc sec",model_rab_p2.offset);
      cpgsch(0.8); cpgtext(0.65,0.91,text); cpgsch(1.4);
      sprintf(text,"DEC(+,2) offset: %.2f arc sec",model_decf_p2.offset);
      cpgsch(0.8); cpgtext(0.65,0.88,text); cpgsch(1.4);
      sprintf(text,"DEC(-,2) offset: %.2f arc sec",model_decb_p2.offset);
      cpgsch(0.8); cpgtext(0.65,0.85,text); cpgsch(1.4);

      // Now plot beam polarisation
      //	      wy1 = 0.475; wy2 = 0.8;
      cpgsch(1.4);
      cpgsvp(0.66666,0.949933,0.15,0.3125);      
      cpgswin(-900,900.,miny,maxy+(maxy-miny)*0.1);
      cpgbox("BCTSN",0,0,"BCTS",0,0);
      cpgsci(4); cpgline(nPlot,result_raf.x,result_raf.pol1);
      cpgsch(0.8); cpgtext(-800,0.9,"RA(+,1)"); cpgsch(1.4);
      cpgsci(2); cpgline(nPlot,result_raf.x,result_raf.pol2);
      cpgsch(0.8); cpgtext(-800,0.75,"RA(+,2)"); cpgsch(1.4);
      cpgsci(1); cpgsls(4); cpgline(2,zx,zy); cpgsls(1);
      cpgsci(1);

      cpgsvp(0.66666,0.949933,0.3125,0.475);      
      cpgswin(-900,900.,miny,maxy+(maxy-miny)*0.1);
      cpgbox("BCTS",0,0,"BCTS",0,0);
      cpgsci(4); cpgline(nPlot,result_rab.x,result_rab.pol1);
      cpgsch(0.8); cpgtext(-800,0.9,"RA(-,1)"); cpgsch(1.4);
      cpgsci(2); cpgline(nPlot,result_rab.x,result_rab.pol2);
      cpgsch(0.8); cpgtext(-800,0.75,"RA(-,2)"); cpgsch(1.4);
      cpgsci(1); cpgsls(4); cpgline(2,zx,zy); cpgsls(1);
      cpgsci(1);

      cpgsvp(0.66666,0.949933,0.475,0.6375);      
      cpgswin(-900,900.,miny,maxy+(maxy-miny)*0.1);
      cpgbox("BCTS",0,0,"BCTS",0,0);
      cpgsci(4); cpgline(nPlot,result_decf.x,result_decf.pol1);
      cpgsch(0.8); cpgtext(-800,0.9,"DEC(+,1)"); cpgsch(1.4);
      cpgsci(2); cpgline(nPlot,result_decf.x,result_decf.pol2);
      cpgsch(0.8); cpgtext(-800,0.75,"DEC(+,2)"); cpgsch(1.4);
      cpgsci(1); cpgsls(4); cpgline(2,zx,zy); cpgsls(1);
      cpgsci(1);

      cpgsvp(0.66666,0.949933,0.6375,0.8);      
      cpgswin(-900,900.,miny,maxy+(maxy-miny)*0.1);
      cpgbox("BCTS",0,0,"BCTS",0,0);
      cpgsci(4); cpgline(nPlot,result_decb.x,result_decb.pol1);
      cpgsch(0.8); cpgtext(-800,0.9,"DEC(-,1)"); cpgsch(1.4);
      cpgsci(2); cpgline(nPlot,result_decb.x,result_decb.pol2);
      cpgsch(0.8); cpgtext(-800,0.75,"DEC(-,2)"); cpgsch(1.4);
      cpgsci(1); cpgsls(4); cpgline(2,zx,zy); cpgsls(1);
      cpgsci(1);

      if (interactive==1)
	cpgcurs(&mx,&my,&key);
    }
  cpgsch(1.4);

  // Beam width

  cpgpage();
  cpgsvp(0.1,0.95,0.15,0.5);
  cpgswin(700,4100,0,29);
  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
  cpglab("Frequency (MHz)","                      Gaussian FWHM (arc minutes)","");
  cpgline(26,freq,beamWidth_ra_p1_fwd);
  cpgsls(2); cpgline(26,freq,beamWidth_ra_p1_bwd); cpgsls(1);
  cpgsci(5);
  cpgline(26,freq,beamWidth_ra_p2_fwd);
  cpgsls(2); cpgline(26,freq,beamWidth_ra_p2_bwd); cpgsls(1);
  cpgsci(1);
  cpgtext(800,2,"Right ascension scans");

  
  cpgsvp(0.1,0.95,0.5,0.85);
  cpgswin(700,4100,0,29);
  cpgbox("ABCTS",0,0,"ABCTSN",0,0);
  cpglab("","","Beam widths");
  
  cpgline(26,freq,beamWidth_dec_p1_fwd);
  cpgsls(2); cpgline(26,freq,beamWidth_dec_p1_bwd); cpgsls(1);
  cpgsci(5);
  cpgline(26,freq,beamWidth_dec_p2_fwd);
  cpgsls(2); cpgline(26,freq,beamWidth_dec_p2_bwd); cpgsls(1);
  cpgsci(1);
  cpgtext(800,2,"Declination scans");
  zx[0] = 3100; zx[1] = 3450;  zy[0] = 21; zy[1] = 21;
  cpgline(2,zx,zy);
  zx[0] = 3100; zx[1] = 3450;  zy[0] = 19; zy[1] = 19;
  cpgsls(2); cpgline(2,zx,zy); cpgsls(1);
  cpgsci(5);
  zx[0] = 3100; zx[1] = 3450;  zy[0] = 17; zy[1] = 17;
  cpgline(2,zx,zy);
  zx[0] = 3100; zx[1] = 3450;  zy[0] = 15; zy[1] = 15;
  cpgsls(2); cpgline(2,zx,zy); cpgsls(1);
  cpgsci(1);
  cpgsch(0.8);
  cpgtext(3500,20.7,"Forward pol 1");
  cpgtext(3500,18.7,"Backward pol 1");
  cpgtext(3500,16.7,"Forward pol 2");
  cpgtext(3500,14.7,"Backward pol 2");
  cpgsch(1.4);

  printf("Results\n");
  {
    float fwhm_ra_p1,fwhm_ra_p2;
    float fwhm_dec_p1,fwhm_dec_p2;
    double fwhm_predict;
    double dAct=64;      // FIX ME: HARDCODE TO PARKES
    double lambda,f;
    double antEff_ra_p1;
    double antEff_dec_p1;
    double antEff_ra_p2;
    double antEff_dec_p2;
    double avA_eff;
    
    for (i=0;i<26;i++)
      {
	fwhm_ra_p1  = (beamWidth_ra_p1_fwd[i]+beamWidth_ra_p1_bwd[i])/2.;
	fwhm_ra_p2  = (beamWidth_ra_p2_fwd[i]+beamWidth_ra_p2_bwd[i])/2.;
	fwhm_dec_p1 = (beamWidth_dec_p1_fwd[i]+beamWidth_dec_p1_bwd[i])/2.;
	fwhm_dec_p2 = (beamWidth_dec_p2_fwd[i]+beamWidth_dec_p2_bwd[i])/2.;

	// SHOULD SET TO A WEIGHTED MEAN -- FIX ME
	f = ra_f->beam[ibeam].bandHeader[i].fc*1e6; 
	lambda = SPEED_LIGHT/f;
	
	fwhm_predict = 1.02 * lambda/dAct * 180./M_PI*60.; // In arc minutes
	antEff_ra_p1 = pow(fwhm_predict/fwhm_ra_p1,2);  // CHECK EQUATION -- FIX ME
	antEff_ra_p2 = pow(fwhm_predict/fwhm_ra_p2,2);  // CHECK EQUATION -- FIX ME
	antEff_dec_p1 = pow(fwhm_predict/fwhm_dec_p1,2);  // CHECK EQUATION -- FIX ME
	antEff_dec_p2 = pow(fwhm_predict/fwhm_dec_p2,2);  // CHECK EQUATION -- FIX ME
      
	avAntEff[i] = (antEff_ra_p1+antEff_ra_p2+antEff_dec_p1+antEff_dec_p2)/4.;
	avA_eff = avAntEff[i] * M_PI *32.0*32.0; // FIX ME -- HARD CODE TO PARKS
	gain_dpfu[i] = avA_eff/BOLTZMANN/2.*1e-26;
	
	printf("%-3d %-10.10s %-8.1f %-8.1f %-5.2f %-5.2f %-5.2f %-5.2f %-5.2f %-5.2f %-5.2f %-5.2f %-5.2f %-5.2f %-5.2f %-5.2f\n",i,ra_f->beam[ibeam].bandHeader[i].label,ra_f->beam[ibeam].bandHeader[i].f0,ra_f->beam[ibeam].bandHeader[i].f1,fwhm_ra_p1,fwhm_ra_p2,fwhm_dec_p1,fwhm_dec_p2,fwhm_predict,antEff_ra_p1,antEff_ra_p2,antEff_dec_p1,antEff_dec_p2,avAntEff[i],gain_dpfu[i],1./gain_dpfu[i]);
      }
  }


  
  // Aperture efficiency and gain_dpfu
  cpgpage();
  cpgsvp(0.1,0.95,0.15,0.5);
  cpgswin(700,4100,0.3,0.9);
  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
  cpglab("Frequency (MHz)","Antenna efficiency","");
  cpgline(26,freq,avAntEff);
  cpgsvp(0.1,0.95,0.5,0.85);
  cpgswin(700,4100,0.3,0.9);
  cpgbox("ABCTS",0,0,"ABCTSN",0,0);
  cpglab("","Gain\\ddpfu\\u(K/Jy)","Derived from beam widths");
  cpgline(26,freq,gain_dpfu);
  
  // Offsets  
  cpgpage();
  cpgsvp(0.1,0.95,0.15,0.5);
  cpgswin(700,4100,-60,60);
  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
  cpglab("Frequency (MHz)","                     Angular offset (arc seconds)","");
  cpgline(26,freq,raOff_p1_fwd);
  cpgsls(2); cpgline(26,freq,raOff_p1_bwd); cpgsls(1);
  cpgsci(5);
  cpgline(26,freq,raOff_p2_fwd);
  cpgsls(2); cpgline(26,freq,raOff_p2_bwd); cpgsls(1);
  cpgsci(1);
  cpgtext(800,45,"Right ascension scans");

  
  cpgsvp(0.1,0.95,0.5,0.85);
  cpgswin(700,4100,-60,60);
  cpgbox("ABCTS",0,0,"ABCTSN",0,0);
  cpgline(26,freq,decOff_p1_fwd);
  cpgsls(2); cpgline(26,freq,decOff_p1_bwd); cpgsls(1);
  cpgsci(5);
  cpgline(26,freq,decOff_p2_fwd);
  cpgsls(2); cpgline(26,freq,decOff_p2_bwd); cpgsls(1);
  cpgsci(1);
  cpgtext(800,45,"Declination scans");

  // Average values
  for (i=0;i<26;i++)
    {
      raOff_p1_avg[i] = 0.5*(raOff_p1_fwd[i]+raOff_p1_bwd[i]);
      raOff_p2_avg[i] = 0.5*(raOff_p2_fwd[i]+raOff_p2_bwd[i]);
      decOff_p1_avg[i] = 0.5*(decOff_p1_fwd[i]+decOff_p1_bwd[i]);
      decOff_p2_avg[i] = 0.5*(decOff_p2_fwd[i]+decOff_p2_bwd[i]);
      raOff_avg[i] = 0.5*(raOff_p1_avg[i]+raOff_p2_avg[i]);
      decOff_avg[i] = 0.5*(decOff_p1_avg[i]+decOff_p2_avg[i]);
    }
   
  cpgpage();
  cpgsvp(0.1,0.95,0.15,0.5);
  cpgswin(700,4100,-30,30);
  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
  cpglab("Frequency (MHz)","                     Angular offset (arc seconds)","");
  cpgline(26,freq,raOff_p1_avg);
  cpgsci(5); cpgline(26,freq,raOff_p2_avg); cpgsci(1);
  cpgsci(2); cpgsls(2); cpgline(26,freq,raOff_avg); cpgsls(1); cpgsci(1);
  cpgtext(800,20,"Right ascension");

  
  cpgsvp(0.1,0.95,0.5,0.85);
  cpgswin(700,4100,-30,30);
  cpgbox("ABCTS",0,0,"ABCTSN",0,0);
  cpglab("","","Averaging forward and backward scans");
  cpgline(26,freq,decOff_p1_avg);
  cpgsci(5); cpgline(26,freq,decOff_p2_avg); cpgsci(1);
  cpgsci(2); cpgsls(2); cpgline(26,freq,decOff_avg); cpgsls(1); cpgsci(1);
  cpgtext(800,20,"Declination");


  cpgend();
  // Create a "pop" file
  {
    FILE *fout;
    char popName[1024];
    double dRA,dDEC;
    double ha_hh,ha_mm,ha;
    double dec_hh,dec_mm,dec;
    double ha_1,ha_2,ha_3,ha_4;
    double para;
    double para_1,para_2,para_3,para_4;
    double dec_1,dec_2,dec_3,dec_4;
    double axd,eld;
    int central;
    int neg=0;
    
    for (i=0;i<26;i++)
      {
	sprintf(popName,"%s_%.2f.pop",ra_f->primary[0].utc0,ra_f->beam[ibeam].bandHeader[i].fc);
	fout = fopen(popName,"w");
	central = (int)(ra_f->beam[ibeam].bandHeader[i].ndump/2.+0.5);  ha_1 = ra_f->beam[ibeam].bandData[i].astro_obsHeader[central].hourAngle;
	central = (int)(ra_b->beam[ibeam].bandHeader[i].ndump/2.+0.5);  ha_2 = ra_b->beam[ibeam].bandData[i].astro_obsHeader[central].hourAngle;
	central = (int)(dec_f->beam[ibeam].bandHeader[i].ndump/2.+0.5); ha_3 = dec_f->beam[ibeam].bandData[i].astro_obsHeader[central].hourAngle;
	central = (int)(dec_b->beam[ibeam].bandHeader[i].ndump/2.+0.5); ha_4 = dec_b->beam[ibeam].bandData[i].astro_obsHeader[central].hourAngle;

	central = (int)(ra_f->beam[ibeam].bandHeader[i].ndump/2.+0.5);  para_1 = ra_f->beam[ibeam].bandData[i].astro_obsHeader[central].paraAngle;
	central = (int)(ra_b->beam[ibeam].bandHeader[i].ndump/2.+0.5);  para_2 = ra_b->beam[ibeam].bandData[i].astro_obsHeader[central].paraAngle;
	central = (int)(dec_f->beam[ibeam].bandHeader[i].ndump/2.+0.5); para_3 = dec_f->beam[ibeam].bandData[i].astro_obsHeader[central].paraAngle;
	central = (int)(dec_b->beam[ibeam].bandHeader[i].ndump/2.+0.5); para_4 = dec_b->beam[ibeam].bandData[i].astro_obsHeader[central].paraAngle;

	central = (int)(ra_f->beam[ibeam].bandHeader[i].ndump/2.+0.5);  dec_1 = ra_f->beam[ibeam].bandData[i].astro_obsHeader[central].decDeg;
	central = (int)(ra_b->beam[ibeam].bandHeader[i].ndump/2.+0.5);  dec_2 = ra_b->beam[ibeam].bandData[i].astro_obsHeader[central].decDeg;
	central = (int)(dec_f->beam[ibeam].bandHeader[i].ndump/2.+0.5); dec_3 = dec_f->beam[ibeam].bandData[i].astro_obsHeader[central].decDeg;
	central = (int)(dec_b->beam[ibeam].bandHeader[i].ndump/2.+0.5); dec_4 = dec_b->beam[ibeam].bandData[i].astro_obsHeader[central].decDeg;

	ha = (ha_1+ha_2+ha_3+ha_4)/4.0; // Average all scans
	ha_hh = (int)(ha*24.0/360.0);	
	ha_mm = (int)((ha*24.0/360.0-ha_hh)*60.);

	para = (para_1+para_2+para_3+para_4)/4.0; // Average all scans
	
	
	dec = (dec_1+dec_2+dec_3+dec_4)/4.0;
	if (dec < 0)
	  {
	    neg=1;
	    dec*=-1;
	  }
	else
	  neg=0;
	
	dec_hh = (int)(dec);
	dec_mm = (int)((dec-dec_hh)*60.0);
	if (neg==1)
	  dec_hh *= -1;
	
	// Pol1
	dRA = raOff_p1_avg[i];          // In arcseconds **FIX ME IF NEEDED
	dDEC = decOff_p1_avg[i];        // In arcseconds
	//	dRA = 8.9;
	//	dDEC = -4.9;
	
	sdhdf_para(dRA,dDEC,para,&axd,&eld);

	// To seconds
	dRA *= 12.0/180.;
	
	// dRA(s) dDEC('') HA_hh HA_mm DEC_hh DEC_mm Source FirstFile RCVR Para AXD ELD
	fprintf(fout,"%6.2f %6.2f %.0f %.0f %.0f %.0f %s %s %s %.2f %6.2f %6.2f\n",dRA,dDEC,ha_hh,ha_mm,dec_hh,dec_mm,ra_f->beamHeader[ibeam].source,ra_f->fname,ra_f->primary[0].rcvr,para,axd,eld);

	// Pol 2
	dRA = raOff_p2_avg[i];          // In arcseconds ** FIX ME IF NEEDED
	dDEC = decOff_p2_avg[i];        // In arcseconds
	sdhdf_para(dRA,dDEC,para,&axd,&eld);

	dRA *= 12.0/180.;

	fprintf(fout,"%6.2f %6.2f %.0f %.0f %.0f %.0f %s %s %s %.2f %6.2f %6.2f\n",dRA,dDEC,ha_hh,ha_mm,dec_hh,dec_mm,ra_f->beamHeader[ibeam].source,ra_f->fname,ra_f->primary[0].rcvr,para,axd,eld);
	
	fclose(fout);
      }
  }



  
  sdhdf_closeFile(ra_f);
  sdhdf_closeFile(ra_b);
  sdhdf_closeFile(dec_f);
  sdhdf_closeFile(dec_b);
  
  free(ra_f);
  free(ra_b);
  free(dec_f);
  free(dec_b);
} 

void sdhdf_calc1MinMax(float *v,int n,float *min,float *max)
{
  int i;
  *min = *max = v[0];
  for (i=0;i<n;i++)
    {
      if (*min > v[i]) *min = v[i];
      if (*max < v[i]) *max = v[i];
    }
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

void processScan(sdhdf_fileStruct *inFile,int ibeam,int iband,float *fy,int *nPlot,int pol)
{
  int i,j,k;
  int npol,nchan,ndump;
  int nchan_cal,ndump_cal;
  double wsum=0;
  double result=0;
  float min,max;

  
  npol = inFile->beam[ibeam].bandHeader[iband].npol;
  nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
  //  nchan_cal = inFile->beam[ibeam].calBandHeader[iband].nchan;
  ndump = inFile->beam[ibeam].bandHeader[iband].ndump;
  //  ndump_cal = inFile->beam[ibeam].calBandHeader[iband].ndump;
  
  sdhdf_loadBandData(inFile,ibeam,iband,1);

  // FIX ME -- NEED A CAL
  //  sdhdf_loadBandData(inFile,ibeam,iband,2);
  //  sdhdf_loadBandData(inFile,ibeam,iband,3);

  for (i=0;i<ndump;i++)
    {
      wsum = 0;
      result = 0;
      for (j=0;j<nchan;j++)
	{
	  if (inFile->beam[ibeam].bandData[iband].astro_data.flag[j] == 0)
	    {
	      if (pol==1)
		result += inFile->beam[ibeam].bandData[iband].astro_data.pol1[j+i*nchan]; 
	      else if (pol==2)
		result += inFile->beam[ibeam].bandData[iband].astro_data.pol2[j+i*nchan]; 
	      wsum++;
	    }
	}
      fy[i] = result/wsum;      
    }
  sdhdf_calc1MinMax(fy,ndump,&min,&max);
  for (i=0;i<ndump;i++)
    fy[i] = (fy[i]-min)/(max-min);
    
  *nPlot = ndump;

  // Process the cal FIX ME
  /* 
  if (iband==6)
    {
      for (j=0;j<nchan_cal;j++)
	{
	  printf("cal: %g %g\n",inFile->beam[ibeam].bandData[iband].cal_on_data.pol1[j],
		 inFile->beam[ibeam].bandData[iband].cal_off_data.pol1[j]);
	}
    }
  */




  
  
  sdhdf_releaseBandData(inFile,ibeam,iband,1);

}


void fitBeam(sdhdf_scanResult scan,int npts,int pol,modelFitStruct *model)
{
  int i;
  int nfit=4;
  double pval[nfit];
  double px[npts];
  double py[npts];
  int nptsInFit;
  lm_status_struct status;
  lm_control_struct control = lm_control_double;


  nptsInFit=0;
  for (i=0;i<npts;i++)
    {

      if (scan.x[i] > model->minx && scan.x[i] < model->maxx)
	{
	  px[nptsInFit] = scan.x[i];
	  
	  if (pol==1)
	    py[nptsInFit] = scan.pol1[i];
	  else
	    py[nptsInFit] = scan.pol2[i];
	  nptsInFit++;
	}
    }

  pval[0] = model->amp; 
  pval[1] = model->offset;
  pval[2] = model->width; 
  pval[3] = model->baseline;

  lmcurve_fit(nfit,pval,nptsInFit,px,py,nonlinearFunc,&control,&status);  

  model->amp      = pval[0];
  model->offset   = pval[1];
  model->width    = fabs(pval[2]);
  model->baseline = pval[3];
  
  //  printf("Output fit = %g %g %g %g\n",pval[0],pval[1],pval[2],pval[3]);
  //  		      fitBeam(result_decb,2);
}


void plotModel(modelFitStruct model,float minx,float maxx)
{
  float px[1024];
  float py[1024];
  int i;

  minx = model.minx;
  maxx = model.maxx;
  for (i=0;i<1024;i++)
    {
      px[i] = minx + i*(maxx-minx)/1024.;
      py[i] = model.amp*exp(-pow((px[i]-model.offset),2)/2./pow(model.width,2))+model.baseline;
    }

  cpgline(1024,px,py);

}
  double amp;
  double offset;
  double width;
  double baseline; 
  double minx;
  double maxx;

void initialiseModel(modelFitStruct *pol1,modelFitStruct *pol2, sdhdf_fileStruct *file,int ibeam,int iband)
{
  double width;
  
  pol1->amp = pol2->amp = 1;
  pol1->offset = pol2->offset = 0;
  pol1->baseline = pol2->baseline = 0;

  // FWHM (Gaussian) = 2.355 sigma
  width = 1.02*SPEED_LIGHT/(file->beam[ibeam].bandHeader[iband].fc*1e6)*180./M_PI/64.0/2.355; // HARD CODE TO 64m---  FIX ME
  
  pol1->width = pol2->width = width*60.*60.*2; // * 2 ??? FIX ME
  pol1->minx  = pol2->minx = -width*60.*60.*2;
  pol1->maxx  = pol2->maxx = +width*60.*60.*2;
}
