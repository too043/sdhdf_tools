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
// ITU-R P.676 gaseous attenuation -- table-free approximate (Annex 2 Sec 1)
// slant-path method. Direct port of atnf-parkes-opacity's P676.java (itself
// a port of the itur Python package's P.676-10 approximate implementation),
// so values match that reference. Valid for f <= 350 GHz, elevation 5-90 deg.
//

#include <math.h>
#include "sdhdf_opacity.h"

static double phi(double rp,double rt,double a,double b,double c,double d)
{
  return pow(rp,a)*pow(rt,b)*exp(c*(1-rp)+d*(1-rt));
}

static double gfunc(double f,double fi)
{
  double r = (f-fi)/(f+fi);
  return 1 + r*r;
}

static double sq(double x)
{
  return x*x;
}

double sdhdf_opacity_waterVapourDensity(double tempC,double relHumidityPct)
{
  double es = 6.1121*exp((18.678 - tempC/234.5)*tempC/(257.14+tempC));
  double e  = es*(relHumidityPct/100.0);
  double tK = tempC + 273.15;
  double rho = 216.7*e/tK;

  if (rho < 0.1) rho = 0.1;
  if (rho > 30.0) rho = 30.0;
  return rho;
}

// Dry-air (oxygen) specific attenuation gamma_0 (dB/km)
static double gamma0(double fGhz,double pHpa,double rho,double tK)
{
  double f  = fGhz;
  double rp = pHpa/1013.0;
  double rt = 288.0/tK;
  double delta,xi1,xi2,xi3,xi4,xi5,xi6,xi7;
  double gamma54,gamma58,gamma60,gamma62,gamma64,gamma66;

  delta = -0.00306*phi(rp,rt,3.211,-14.94,1.583,-16.37);

  xi1 = phi(rp,rt,0.0717,-1.8132,0.0156,-1.6515);
  xi2 = phi(rp,rt,0.5146,-4.6368,-0.1921,-5.7416);
  xi3 = phi(rp,rt,0.3414,-6.5851,0.2130,-8.5854);
  xi4 = phi(rp,rt,-0.0112,0.0092,-0.1033,-0.0009);
  xi5 = phi(rp,rt,0.2705,-2.7192,-0.3016,-4.1033);
  xi6 = phi(rp,rt,0.2445,-5.9191,0.0422,-8.0719);
  xi7 = phi(rp,rt,-0.1833,6.5589,-0.2402,6.131);

  if (f <= 54.0)
    return ((7.2*pow(rt,2.8))/(f*f+0.34*rp*rp*pow(rt,1.6))
	    + (0.62*xi3)/(pow(54.0-f,1.16*xi1)+0.83*xi2))
      * f*f*rp*rp*1e-3;

  gamma54 = 2.192*phi(rp,rt,1.8286,-1.9487,0.4051,-2.8509);
  gamma58 = 12.59*phi(rp,rt,1.0045,3.5610,0.1588,1.2834);
  gamma60 = 15.00*phi(rp,rt,0.9003,4.1335,0.0427,1.6088);
  gamma62 = 14.28*phi(rp,rt,0.9886,3.4176,0.1827,1.3429);
  gamma64 = 6.819*phi(rp,rt,1.4320,0.6258,0.3177,-0.5914);
  gamma66 = 1.908*phi(rp,rt,2.0717,-4.1404,0.4910,-4.8718);

  if (f <= 60.0)
    return exp(log(gamma54)/24.0*(f-58)*(f-60)
	       - log(gamma58)/8.0*(f-54)*(f-60)
	       + log(gamma60)/12.0*(f-54)*(f-58));
  if (f <= 62.0)
    return gamma60 + (gamma62-gamma60)*(f-60)/2.0;
  if (f <= 66.0)
    return exp(log(gamma62)/8.0*(f-64)*(f-66)
	       - log(gamma64)/4.0*(f-62)*(f-66)
	       + log(gamma66)/8.0*(f-62)*(f-64));
  if (f <= 120.0)
    return (3.02e-4*pow(rt,3.5)
	    + (0.283*pow(rt,3.8))/(pow(f-118.75,2)+2.91*rp*rp*pow(rt,1.6))
	    + (0.502*xi6*(1-0.0163*xi7*(f-66)))
	      /(pow(f-66,1.4346*xi4)+1.15*xi5))
      * f*f*rp*rp*1e-3;

  return ((3.02e-4)/(1+1.9e-5*pow(f,1.5))
	  + (0.283*pow(rt,0.3))/(pow(f-118.75,2)+2.91*rp*rp*pow(rt,1.6)))
    * f*f*rp*rp*pow(rt,3.5)*1e-3 + delta;
}

// Water-vapour specific attenuation gamma_w (dB/km)
static double gammaw(double fGhz,double pHpa,double rho,double tK)
{
  double f  = fGhz;
  double rp = pHpa/1013.0;
  double rt = 288.0/tK;
  double eta1 = 0.955*rp*pow(rt,0.68) + 0.006*rho;
  double eta2 = 0.735*rp*pow(rt,0.50) + 0.0353*pow(rt,4)*rho;
  double sum;

  sum =   (3.98*eta1*exp(2.23*(1-rt)))/(sq(f-22.235)+9.42*eta1*eta1)*gfunc(f,22.0)
        + (11.96*eta1*exp(0.70*(1-rt)))/(sq(f-183.310)+11.14*eta1*eta1)
        + (0.081*eta1*exp(6.44*(1-rt)))/(sq(f-321.226)+6.29*eta1*eta1)
        + (3.660*eta1*exp(1.60*(1-rt)))/(sq(f-325.153)+9.22*eta1*eta1)
        + (25.37*eta1*exp(1.09*(1-rt)))/sq(f-380.000)
        + (17.40*eta1*exp(1.46*(1-rt)))/sq(f-448.000)
        + (844.6*eta1*exp(0.17*(1-rt)))/sq(f-557.000)*gfunc(f,557.0)
        + (290.0*eta1*exp(0.41*(1-rt)))/sq(f-752.000)*gfunc(f,752.0)
        + (8.3328e4*eta2*exp(0.99*(1-rt)))/sq(f-1780.00)*gfunc(f,1780.0);

  return sum*f*f*pow(rt,2.5)*rho*1e-4;
}

// Equivalent heights h0, hw (km). Per P.676 Annex 2 this uses the total
// pressure p = P + e (e = water-vapour partial pressure).
static void equivalentHeights(double fGhz,double pTotalHpa,double *h0Out,double *hwOut)
{
  double f  = fGhz;
  double rp = pTotalHpa/1013.0;
  double t1,t2,t3,h0,sigmaw,hw;

  t1 = 4.64/(1+0.066*pow(rp,-2.3))*exp(-sq((f-59.7)/(2.87+12.4*exp(-7.9*rp))));
  t2 = (0.14*exp(2.21*rp))/(sq(f-118.75)+0.031*exp(2.2*rp));
  t3 = 0.0114/(1+0.14*pow(rp,-2.6))*f
    *(-0.0247+0.0001*f+1.61e-6*f*f)
    /(1-0.0169*f+4.1e-5*f*f+3.2e-7*f*f*f);

  h0 = 6.1/(1+0.17*pow(rp,-1.1))*(1+t1+t2+t3);
  if (f < 70.0)
    {
      double cap = 10.7*pow(rp,0.3);
      if (h0 > cap) h0 = cap;
    }

  sigmaw = 1.013/(1+exp(-8.6*(rp-0.57)));
  hw = 1.66*(1
	     + (1.39*sigmaw)/(sq(f-22.235)+2.56*sigmaw)
	     + (3.37*sigmaw)/(sq(f-183.31)+4.69*sigmaw)
	     + (1.58*sigmaw)/(sq(f-325.1)+2.89*sigmaw));

  *h0Out = h0;
  *hwOut = hw;
}

double sdhdf_opacity_attenuationDb(double freqGHz,double elevationDeg,double rho,double pressureHpa,double tempK)
{
  double g0 = gamma0(freqGHz,pressureHpa,rho,tempK);
  double gw = gammaw(freqGHz,pressureHpa,rho,tempK);
  double e  = rho*tempK/216.7; // water-vapour partial pressure (hPa)
  double h0,hw,a0,aw;

  equivalentHeights(freqGHz,pressureHpa+e,&h0,&hw);
  a0 = g0*h0;
  aw = gw*hw;
  return (a0+aw)/sin(elevationDeg*M_PI/180.0);
}

double sdhdf_opacity_tau(double freqGHz,double elevationDeg,double rho,double pressureHpa,double tempK)
{
  double aDb = sdhdf_opacity_attenuationDb(freqGHz,elevationDeg,rho,pressureHpa,tempK);
  return aDb/(10.0/log(10.0)); // dB -> nepers
}

double sdhdf_opacity_correctionFactor(double freqMHz,double elevationDeg,double tempC,double pressureHpa,double relHumidityPct)
{
  double freqGHz = freqMHz/1000.0;
  double tempK   = tempC + 273.15;
  double rho     = sdhdf_opacity_waterVapourDensity(tempC,relHumidityPct);
  double aDb     = sdhdf_opacity_attenuationDb(freqGHz,elevationDeg,rho,pressureHpa,tempK);
  double tauSlant,correction;

  if (aDb < 0.0)  aDb = 0.0;
  if (aDb > 10.0) aDb = 10.0;

  tauSlant   = aDb/(10.0/log(10.0));
  correction = exp(tauSlant);
  if (correction < 1.0) correction = 1.0;
  if (correction > 3.0) correction = 3.0;

  return correction;
}
