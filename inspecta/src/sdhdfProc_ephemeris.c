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
#include <erfa.h>
#include <calceph.h>

long sdhdf_loadEOP(sdhdf_eopStruct *eop)
{
  long n=0;
  FILE *fin;
  char line[4096];
  int idummy;
  float year;
  double mjd,x,y,dut1;
  char runtimeDir[1024];
  char fname[1024];


    if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdf_convertTo requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  sprintf(fname,"%s/earth/eopc04_IAU2000.62-now",runtimeDir);
  fin = fopen(fname,"r");
  while(!feof(fin))
    {
      if (fgets(line,4096,fin)!=NULL)
	{
	  if (sscanf(line,"%f %d %d %lf %lf %lf %lf",&year,&idummy,&idummy,&mjd,&x,&y,&dut1)==7)
	    {
	      if (mjd > 37665 && mjd < 70000)
		{
		  eop[n].mjd = mjd;
		  eop[n].x = x/3600.*M_PI/180.0;
		  eop[n].y = y/3600.*M_PI/180.0;
		  eop[n].dut1 = dut1;
		  n++;
		}
	    }
	}
    }
  fclose(fin);
  //  printf("Have loaded %d entries from the EOP file\n",n);
  return n;
}

//double sdhdf_calcVoverC(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,sdhdf_eopStruct *eop,int nEOP,int lsr,char *ephemName)

// Load in an array of MJDs for a given observatory and return an array of vOverC values
//double sdhdf_calcVoverC(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,sdhdf_eopStruct *eop,int nEOP,int lsr,char *ephemName)
// bary_lsr = 0 (or 1) = bary, = 2 = lsr
void sdhdf_calcVoverC(double *mjd,double *raDeg,double *decDeg,int nvals,double *vOverC,char *tel,char *ephemName,sdhdf_eopStruct *eop,int nEOP,int bary_lsr)
{
  t_calcephbin *eph;
  int i,j,k,ii;
  long double jd;
  double jd0,jd1;
  double earth_ssb[6],vobs[3],pos[3];
  double r,ctheta,voverc,freqf,pospos;  
  double alpha,delta,ca,sa,cd,sd;
  double observatory_earth[3];
  double zenith[3];
  double siteVel[3];
  double trs[3];
  double long_grs80, lat_grs80,height_grs80;
  double vlsr_ra,vlsr_dec;

  char runtimeDir[1024];
  char fname[1024];
  FILE *fin;
  char fname1[1024];
  char fname2[1024];
  char loadLine[1024];
  char observatoryDir[1024]="NULL";
  char label[1024];
  char entry[1024];
	

  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdf_convertTo requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  //  sprintf(fname,"%s/ephemeris/DE436.1950.2050",runtimeDir);
  sprintf(fname,"%s/ephemeris/%s",runtimeDir,ephemName);

  eph = calceph_open(fname);
  if (eph) {
    //    printf("Successfully opened ephemeris\n");
  } else {
    printf("Error: unable to open ephemeris\n");
    exit(1);
  }

  for (ii=0;ii<nvals;ii++)
    {
      jd = mjd[ii]+2400000.5;
      jd0 = (double)((int)jd);
      jd1 = (double)(jd-(int)jd);
      calceph_compute_unit(eph,jd0,jd1,3,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,earth_ssb);
      for (i=0;i<6;i++)
	earth_ssb[i]*=1000.0; // Convert to m and m/s
      for (i=0;i<3;i++)
	vobs[i] = earth_ssb[i+3];
      
      //      printf("%g %g %g %g %g %g %g\n",mjd[ii],earth_ssb[0],earth_ssb[1],earth_ssb[2],earth_ssb[3],earth_ssb[4],earth_ssb[5]);
      
      
      // Get site velocity
      // Obtain the coordinates of the observatory
      if (ii==0)
	{
	  sprintf(fname1,"%s/observatory/observatories.list",runtimeDir);
	  fin = fopen(fname1,"r");
	  while (!feof(fin))
	    {
	      if (fgets(loadLine,1024,fin)!=NULL)
		{
		  if (loadLine[0]!='#')
		    {
		      if (strstr(loadLine,tel)!=NULL)
			{
			  sscanf(loadLine,"%s",observatoryDir);
			  break;
			}
		    }
		}	
	    }
	  fclose(fin);
	  if (strcmp(observatoryDir,"NULL")==0)
	    {
	      printf("ERROR: in sdhdfProc_ephemeris.c - unable to find observatory %s in %s\n",tel,fname1);
	      exit(1);
	    }
	  sprintf(fname2,"%s/observatory/%s/observatory.properties",runtimeDir,observatoryDir);
	  fin = fopen(fname2,"r");
	  while (!feof(fin))
	    {
	      if (fgets(loadLine,1024,fin)!=NULL)
		{
		  if (loadLine[0]!='#')
		    {
		      sscanf(loadLine,"%s %s",label,entry);
		      if (strcasecmp(label,"antenna_x")==0) sscanf(entry,"%lf",&trs[0]);
		      if (strcasecmp(label,"antenna_y")==0) sscanf(entry,"%lf",&trs[1]);
		      if (strcasecmp(label,"antenna_z")==0) sscanf(entry,"%lf",&trs[2]);
		    }
		}
	    }
	  fclose(fin);
	  printf("Loaded observatory coordinates (x,y,z) = %g %g %g (telescope %s)\n",trs[0],trs[1],trs[2],tel);
	  sdhdf_ITRF_to_GRS80(trs[0],trs[1],trs[2],&long_grs80,&lat_grs80,&height_grs80);
	  
	  zenith[0] = height_grs80 * cos(long_grs80) * cos(lat_grs80);
	  zenith[1] = height_grs80 * sin(long_grs80) * cos(lat_grs80);
	  zenith[2] = height_grs80 * sin(lat_grs80);
	  
	}
      
      
      //  printf("Inputs to IAU2000B: %g %g %g %g\n",trs[0],trs[1],trs[2],(double)inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].mjd);
      
      sdhdf_obsCoord_IAU2000B(trs, zenith,mjd[ii],mjd[ii],
			      observatory_earth, zenith, siteVel,eop,nEOP);
      
      //  printf("Observatory_Earth = %g %g %g\n",observatory_earth[0],observatory_earth[1],observatory_earth[2]);
      //  printf("Site velocity = %g %g %g\n",siteVel[0],siteVel[1],siteVel[2]);
      
      for (i=0;i<3;i++)
	vobs[i] += siteVel[i];
      // Convert to the local standard of rest
      
      // From Carl Heiles' code:
      // THE STANDARD LSR IS DEFINED AS FOLLOWS: THE SUN MOVES AT 20.0 KM/S
      // TOWARD RA=18.0H, DEC=30.0 DEG IN 1900 EPOCH COORDS
      // using PRECESS, this works out to ra=18.063955 dec=30.004661 in 2000 coords.
      if (bary_lsr==2) // Approx LSRK
	{
	  vlsr_ra = 2*M_PI*18./24.;
	  vlsr_dec = 30*M_PI/180.0;
	  vobs[0] += 20e3*cos(vlsr_dec)*cos(vlsr_ra);
	  vobs[1] += 20e3*cos(vlsr_dec)*sin(vlsr_ra);
	  vobs[2] += 20e3*sin(vlsr_dec);
	}
      else if (bary_lsr==3) // More precise LSRK
	{
	  vlsr_ra = 2*M_PI*18.063955/24.;
	  vlsr_dec = 30.004661*M_PI/180.0;
	  vobs[0] += 20e3*cos(vlsr_dec)*cos(vlsr_ra);
	  vobs[1] += 20e3*cos(vlsr_dec)*sin(vlsr_ra);
	  vobs[2] += 20e3*sin(vlsr_dec);
	}
      
      // Calculate Doppler shift
      // (ignoring Einstein delay)
      
      alpha = raDeg[ii]*M_PI/180.;
      delta = decDeg[ii]*M_PI/180.;
      
      //      printf("Source position = %g %g\n",alpha,delta);
      ca = cos(alpha);
      sa = sin(alpha);
      cd = cos(delta);
      sd = sin(delta);
      
      pos[0] = ca*cd;
      pos[1] = sa*cd;
      pos[2] = sd;
      
      pospos = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
      for (k=0;k<3;k++)
	pos[k] /= pospos;
      
      vOverC[ii] = sdhdf_dotproduct(pos,vobs)/SPEED_LIGHT;
    }
  calceph_close(eph);
      
  //  printf("voverc inputs: %g %g %g | %g %g %g\n",pos[0],pos[1],pos[2],vobs[0],vobs[1],vobs[2]);
  //  printf("mjd = %g %g\n",inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].mjd,voverc);
}

void sdhdf_obsCoord_IAU2000B(double *observatory_trs,
			     double *zenith_trs,
			     long double tt_mjd, long double utc_mjd,
			     double *observatory_crs,
			     double *zenith_crs,
			     double *observatory_velocity_crs,sdhdf_eopStruct *eop,int nEOP)
{
  double trs[2][3], crs[2][3], north[3], pole_itrs[3], omega_itrs[3];
  double t2c[3][3], polarmotion[3][3];
  double dut1, dut1dot, eradot, xp, yp;
  long double tt_jd = tt_mjd + (long double)2400000.5, ut1_jd;
  double tt_jd1 = (int)tt_jd, tt_jd2 = tt_jd-tt_jd1;
  double ut1_jd1, ut1_jd2;
  double sprime=0.0;
  int i;
  // Get Earth orientation parameters
  //
  sdhdf_interpolate_EOP((double)utc_mjd, &xp, &yp, &dut1, &dut1dot, eop, nEOP);
  //  printf("EOP input %g\n",(double)utc_mjd);
  //  printf("EOP output %g %g %g %g\n",xp,yp,dut1,dut1dot);
  
  ut1_jd = utc_mjd + dut1/86400.0 + (long double)2400000.5;
  ut1_jd1 = (int)ut1_jd;
  ut1_jd2 = ut1_jd-ut1_jd1;

  // Stick the site position vector in
  // Note that 'trs' is two dimensional
  eraCp(observatory_trs, trs[0]);

  // Work out site velocity in CRS... 
  // first we need to know the angular velocity vector
  eraPom00(xp, yp, sprime, polarmotion); // polar motion matrix
  north[0] = north[1] = 0.0; north[2] = 1.0; // Vector to +ve pole
  eraRxp(polarmotion, north, pole_itrs); // Spin pole in ITRS 
  eradot = 2.0*M_PI*1.00273781191135448*(1.0+dut1dot)/86400.0;
  eraSxp(eradot, pole_itrs, omega_itrs); // Angular velocity in ITRS (rad/s)
  eraPxp(omega_itrs, trs[0], trs[1]); // Tangential velocity (m/s)
  
  // Get the Celestial->Terrestrial matrix
  eraC2t00b(tt_jd1, tt_jd2, ut1_jd1, ut1_jd2, xp, yp, t2c);
  
  // Multiply the itrs position/velocity vector by its transpose (=inverse)
  // to transform trs->crs
  eraTrxpv(t2c, trs, crs);

  // NOTE:::::: DON'T NEED ZENITH STUFF ???
  //iauTrxp(t2c, zenith_trs, zenith_crs);
  for (i=0;i<3;i++)
    {
      observatory_crs[i]          = crs[0][i];
      observatory_velocity_crs[i] = crs[1][i];
      //      printf("crs: %g %g\n",zenith_trs[i],zenith_crs[i],crs[0][i],crs[1][i]);
    }
  
}

void sdhdf_interpolate_EOP(double mjd, double *xp, double *yp, double *dut1, double *dut1dot, sdhdf_eopStruct *eop,int nEOP)
{
  long i;
  long isamp;
  double leap=0.0;
  double f;
  
  // Find first sample after requested time
  for (i=0;i<nEOP;i++)
    {
      if (eop[i].mjd > mjd)
	{
	  isamp = i;
	  break;
	}
    }

  if (mjd > eop[nEOP-1].mjd)
    {
      isamp = nEOP-1;
      mjd = eop[nEOP-1].mjd;
      printf("WARNING: You need to update the EOP file in the SDHDF_RUNTIME/earth directory.\n");
    }
  // cope with leap second ... take it off second point as jump happens
  // right AT second point
  if (eop[isamp].dut1 - eop[isamp-1].dut1 > 0.5)
    leap = 1.0;
  else if (eop[isamp].dut1 - eop[isamp-1].dut1 < -0.5)
    leap = -1.0;
  else
    leap = 0.0;

  *dut1dot = (eop[isamp].dut1 - leap - eop[isamp-1].dut1) / 86400.0;

  /* interpolate */
  f = (mjd - eop[isamp-1].mjd) / (eop[isamp].mjd - eop[isamp-1].mjd);
  *xp = eop[isamp-1].x + f*(eop[isamp].x - eop[isamp-1].x);
  *yp = eop[isamp-1].y + f*(eop[isamp].y - eop[isamp-1].y);
  *dut1 = eop[isamp-1].dut1 + f*(eop[isamp].dut1 - leap - eop[isamp-1].dut1);
  
  //  printf("HERE EOP WITH %d %g %g %g [ %g ]  (%g %g %g %g)\n",isamp,leap,eop[isamp].y,eop[isamp-1].y,f,*xp,*yp,*dut1,*dut1dot);
  
}

// Geocentric to geodetic.
// Uses Vermeille (2004)'s method:
//http://www.springerlink.com/app/home/contribution.asp?wasp=08ea5d2c4c62464789a7961196d84ab5&referrer=parent&backto=issue,11,18;journal,9,85;linkingpublicationresults,1:100435,1
void sdhdf_ITRF_to_GRS80(double x,double y,double z,double *long_grs80,double *lat_grs80,double *height_grs80)
{
    double p = (x*x + y*y)/ (GRS80_A*GRS80_A);
    double esq = GRS80_F*(2.0-GRS80_F);
    double q = (1.0-esq)/(GRS80_A*GRS80_A)*z*z;
    double r = (p+q-esq*esq)/6.0;
    double s = esq*esq*p*q/(4*r*r*r);
    double t = pow(1.0+s+sqrt(s*(2.0+s)), 1.0/3.0);
    double u = r*(1.0+t+1.0/t);
    double v = sqrt(u*u+esq*esq*q);
    double w = esq*(u+v-q)/(2.0*v);
    double k = sqrt(u+v+w*w)-w;
    double D = k*sqrt(x*x+y*y)/(k+esq);

    *height_grs80 = (k+esq-1.0)/k * sqrt(D*D+z*z);
    *lat_grs80 = 2.0*atan2(z, D+sqrt(D*D+z*z));
    if (y >= 0.0)
        *long_grs80 =  0.5*M_PI - 2.0*atan2(x, sqrt(x*x+y*y)+y);
    else
        *long_grs80 = -0.5*M_PI + 2.0*atan2(x, sqrt(x*x+y*y)-y);

    //    printf("grs80 = %g %g %g\n",*height_grs80,*lat_grs80,*long_grs80);
}
