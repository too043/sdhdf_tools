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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "inspecta.h"
#include <complex.h>

double turn_deg(double turn){
 
  /* Converts double turn to string "sddd.ddd" */
  return turn*360.0;
}

/* ************************************************************** *
 * Spline routines                                                *
 * ************************************************************** */
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))


void TKcmonot (int n, float *x, float *y, float **yd)
/*****************************************************************************
compute cubic interpolation coefficients via the Fritsch-Carlson method,
which preserves monotonicity
******************************************************************************
Input:
n		number of samples
x  		array[n] of monotonically increasing or decreasing abscissae
y		array[n] of ordinates

Output:
yd		array[n][4] of cubic interpolation coefficients (see notes)
******************************************************************************
Notes:
The computed cubic spline coefficients are as follows:
yd[i][0] = y(x[i])    (the value of y at x = x[i])
yd[i][1] = y'(x[i])   (the 1st derivative of y at x = x[i])
yd[i][2] = y''(x[i])  (the 2nd derivative of y at x = x[i])
yd[i][3] = y'''(x[i]) (the 3rd derivative of y at x = x[i])

To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
use the computed coefficients as follows:
y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))

The Fritsch-Carlson method yields continuous 1st derivatives, but 2nd
and 3rd derivatives are discontinuous.  The method will yield a
monotonic interpolant for monotonic data.  1st derivatives are set
to zero wherever first divided differences change sign.

For more information, see Fritsch, F. N., and Carlson, R. E., 1980,
Monotone piecewise cubic interpolation:  SIAM J. Numer. Anal., v. 17,
n. 2, p. 238-246.

Also, see the book by Kahaner, D., Moler, C., and Nash, S., 1989,
Numerical Methods and Software, Prentice Hall.  This function was
derived from SUBROUTINE PCHEZ contained on the diskette that comes
with the book.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 09/30/89
Modified:  Dave Hale, Colorado School of Mines, 02/28/91
	changed to work for n=1.
Modified:  Dave Hale, Colorado School of Mines, 08/04/91
	fixed bug in computation of left end derivative
*****************************************************************************/
{
  int i;
  float h1,h2,del1,del2,dmin,dmax,hsum,hsum3,w1,w2,drat1,drat2,divdf3;
  
  /* copy ordinates into output array */
  for (i=0; i<n; i++)
    yd[i][0] = y[i];
  
  /* if n=1, then use constant interpolation */
  if (n==1) {
    yd[0][1] = 0.0;
    yd[0][2] = 0.0;
    yd[0][3] = 0.0;
    return;
    
    /* else, if n=2, then use linear interpolation */
  } else if (n==2) {
    yd[0][1] = yd[1][1] = (y[1]-y[0])/(x[1]-x[0]);
    yd[0][2] = yd[1][2] = 0.0;
    yd[0][3] = yd[1][3] = 0.0;
    return;
  }
  
  /* set left end derivative via shape-preserving 3-point formula */
  h1 = x[1]-x[0];
  h2 = x[2]-x[1];
  hsum = h1+h2;
  del1 = (y[1]-y[0])/h1;
  del2 = (y[2]-y[1])/h2;
  w1 = (h1+hsum)/hsum;
  w2 = -h1/hsum;
  yd[0][1] = w1*del1+w2*del2;
  if (yd[0][1]*del1<=0.0)
    yd[0][1] = 0.0;
  else if (del1*del2<0.0) {
    dmax = 3.0*del1;
    if (ABS(yd[0][1])>ABS(dmax)) yd[0][1] = dmax;
  }
  
  /* loop over interior points */
  for (i=1; i<n-1; i++) {
    
    /* compute intervals and slopes */
    h1 = x[i]-x[i-1];
    h2 = x[i+1]-x[i];
    hsum = h1+h2;
    del1 = (y[i]-y[i-1])/h1;
    del2 = (y[i+1]-y[i])/h2;
    
    /* if not strictly monotonic, zero derivative */
    if (del1*del2<=0.0) {
      yd[i][1] = 0.0;
      
      /*
       * else, if strictly monotonic, use Butland's formula:
       *      3*(h1+h2)*del1*del2
       * -------------------------------
       * ((2*h1+h2)*del1+(h1+2*h2)*del2)
       * computed as follows to avoid roundoff error
       */
    } else {
      hsum3 = hsum+hsum+hsum;
      w1 = (hsum+h1)/hsum3;
      w2 = (hsum+h2)/hsum3;
      dmin = MIN(ABS(del1),ABS(del2));
      dmax = MAX(ABS(del1),ABS(del2));
      drat1 = del1/dmax;
      drat2 = del2/dmax;
      yd[i][1] = dmin/(w1*drat1+w2*drat2);
    }
  }
  
  /* set right end derivative via shape-preserving 3-point formula */
  w1 = -h2/hsum;
  w2 = (h2+hsum)/hsum;
  yd[n-1][1] = w1*del1+w2*del2;
  if (yd[n-1][1]*del2<=0.0)
    yd[n-1][1] = 0.0;
  else if (del1*del2<0.0) {
    dmax = 3.0*del2;
    if (ABS(yd[n-1][1])>ABS(dmax)) yd[n-1][1] = dmax;
  }
  
  /* compute 2nd and 3rd derivatives of cubic polynomials */
  for (i=0; i<n-1; i++) {
    h2 = x[i+1]-x[i];
    del2 = (y[i+1]-y[i])/h2;
    divdf3 = yd[i][1]+yd[i+1][1]-2.0*del2;
    yd[i][2] = 2.0*(del2-yd[i][1]-divdf3)/h2;
    yd[i][3] = (divdf3/h2)*(6.0/h2);
  }
  yd[n-1][2] = yd[n-2][2]+(x[n-1]-x[n-2])*yd[n-2][3];
  yd[n-1][3] = yd[n-2][3];
}

float sdhdf_splineValue(float x,int n,float *xSpline,float **yd)
{
  float val;
  float h;
  float x0,dx;
  int i0;
  
  // To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
  // use the computed coefficients as follows:
  // y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))

  // Assume equal sampling
  // FIX ME -- SHOULDN'T DO THIS AS SHOULD REMOVE RFI FIRST IN THE INTERPOLATION
  //
  x0 = xSpline[0];
  dx = xSpline[1] - x0;

  i0 = (x - x0)/dx;
  if (i0 < 0) i0 = 0;
  if (i0 > xSpline[n-1]) i0 = n-1;

  h = x - xSpline[i0];
  val = yd[i0][0]+h*(yd[i0][1]+h*(yd[i0][2]/2.0+h*yd[i0][3]/6.0));
  return val;
}

// Interpolate the spline fit on to a given set of x-values
void TKspline_interpolate(int n,float *x,float **yd,float *interpX,
			  float *interpY,int nInterp)
{
  float h;
  int jpos;
  int i,j;

  //To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
  //use the computed coefficients as follows:
  //y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))
  
  // Assume the data are sorted so x[i] < x[i+1]
  for (i=0;i<nInterp;i++)
    {
      jpos=-1;
      for (j=0;j<n-1;j++)
	{
	  if (interpX[i]>=x[j] && interpX[i]<x[j+1])
	    {
	      jpos=j;
	      break;
	    }
	}
      if (jpos!=-1)
	{
	  h = interpX[i]-x[jpos];
	  interpY[i] = yd[jpos][0]+h*(yd[jpos][1]+h*(yd[jpos][2]/2.0+h*yd[jpos][3]/6.0));
	}
      else
	interpY[i] = 0.0;
    }

}


int turn_dms(double turn, char *dms){
  
  /* Converts double turn to string "sddd:mm:ss.sss" */
  
  int dd, mm, isec;
  double trn, sec;
  char sign;
  
  sign=' ';
  if (turn < 0.){
    sign = '-';
    trn = -turn;
  }
  else{
    sign = '+';
    trn = turn;
  }
  dd = trn*360.;
  mm = (trn*360.-dd)*60.;
  sec = ((trn*360.-dd)*60.-mm)*60.;
  isec = (sec*1000. +0.5)/1000;
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        dd=dd+1;
      }
    }
  sprintf(dms,"%c%02d:%02d:%010.7f",sign,dd,mm,sec);
  return 0;
}


int turn_hms(double turn, char *hms){
 
  /* Converts double turn to string " hh:mm:ss.ssss" */
  
  int hh, mm, isec;
  double sec;

  hh = turn*24.;
  mm = (turn*24.-hh)*60.;
  sec = ((turn*24.-hh)*60.-mm)*60.;
  isec = (sec*10000. +0.5)/10000;
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        hh=hh+1;
        if(hh==24){
          hh=0;
        }
      }
    }

  sprintf(hms," %02d:%02d:%010.7f",hh,mm,sec);
  return 0;
}


double hms_turn(char *line){

  /* Converts string " hh:mm:ss.ss" or " hh mm ss.ss" to double turn */
  
  int i;int turn_hms(double turn, char *hms);
  double hr, min, sec, turn=0;
  char hold[MAX_STRLEN];

  strcpy(hold,line);

  /* Get rid of ":" */
  for(i=0; *(line+i) != '\0'; i++)if(*(line+i) == ':')*(line+i) = ' ';

  i = sscanf(line,"%lf %lf %lf", &hr, &min, &sec);
  if(i > 0){
    turn = hr/24.;
    if(i > 1)turn += min/1440.;
    if(i > 2)turn += sec/86400.;
  }
  if(i == 0 || i > 3)turn = 1.0;


  strcpy(line,hold);

  return turn;
}

double dms_turn(char *line){

  /* Converts string "-dd:mm:ss.ss" or " -dd mm ss.ss" to double turn */
  
  int i;
  char *ic, ln[40];
  double deg, min, sec, sign, turn=0;

  /* Copy line to internal string */
  strcpy(ln,line);

  /* Get rid of ":" */
  for(i=0; *(ln+i) != '\0'; i++)if(*(ln+i) == ':')*(ln+i) = ' ';

  /* Get sign */
  if((ic = strchr(ln,'-')) == NULL)
     sign = 1.;
  else {
     *ic = ' ';
     sign = -1.;
  }

  /* Get value */
  i = sscanf(ln,"%lf %lf %lf", &deg, &min, &sec);
  if(i > 0){
    turn = deg/360.;
    if(i > 1)turn += min/21600.;
    if(i > 2)turn += sec/1296000.;
    if(turn >= 1.0)turn = turn - 1.0;
    turn *= sign;
  }
  if(i == 0 || i > 3)turn =1.0;

  return turn;
}



// FIX ME: Check this one
// Conjugate and transpose
void sdhdf_complex_matrix_2x2_dagger(double complex R[2][2],double complex Rdag[2][2])
{
  Rdag[0][0] = creal(R[0][0])-I*cimag(R[0][0]);
  Rdag[0][1] = creal(R[1][0])-I*cimag(R[1][0]);
  Rdag[1][0] = creal(R[0][1])-I*cimag(R[0][1]);
  Rdag[1][1] = creal(R[1][1])-I*cimag(R[1][1]);
}

// CHECK THIS ONE VERY CAREFULLY --- FIX ME
// return 0 if okay
// return 1 if bad
int sdhdf_complex_matrix_2x2_inverse(double complex a[2][2],double complex inv[2][2])
{
  double complex det;
  double rdet;
  
  det = a[0][0]*a[1][1] - a[1][0]*a[0][1];
  rdet = creal(det);
  if (isnan(rdet) || rdet == 0)
    return 1;


  inv[0][0] =  1.0/rdet * a[1][1];
  inv[1][0] = -1.0/rdet * a[1][0];
  inv[0][1] = -1.0/rdet * a[0][1];
  inv[1][1] =  1.0/rdet * a[0][0];
  return 0;
}


void sdhdf_complex_matrix_2x2(double complex J[2][2],double complex e00,double complex e10,double complex e01,double complex e11)
{
  //  printf("Entering with %g %g %g %g\n",creal(e00),creal(e10),creal(e01),creal(e11));
  //  printf("Entering with %g %g %g %g\n",cimag(e00),cimag(e10),cimag(e01),cimag(e11));
  J[0][0] = e00;
  J[1][0] = e10;
  J[0][1] = e01;
  J[1][1] = e11;
}

void sdhdf_copy_complex_matrix_2x2(double complex to[2][2],double complex from[2][2])
{
  to[0][0] = from[0][0];
  to[1][0] = from[1][0];
  to[0][1] = from[0][1];
  to[1][1] = from[1][1];
}

void sdhdf_multiply_complex_matrix_2x2(double complex a[2][2],double complex b[2][2])
{
  double a_11r,a_11i;
  double a_12r,a_12i;
  double a_21r,a_21i;
  double a_22r,a_22i;

  double b_11r,b_11i;
  double b_12r,b_12i;
  double b_21r,b_21i;
  double b_22r,b_22i;

  double c_11r,c_11i;
  double c_12r,c_12i;
  double c_21r,c_21i;
  double c_22r,c_22i;


  a_11r = creal(a[0][0]); a_11i = cimag(a[0][0]);
  a_12r = creal(a[1][0]); a_12i = cimag(a[1][0]);
  a_21r = creal(a[0][1]); a_21i = cimag(a[0][1]);
  a_22r = creal(a[1][1]); a_22i = cimag(a[1][1]);

  b_11r = creal(b[0][0]); b_11i = cimag(b[0][0]);
  b_12r = creal(b[1][0]); b_12i = cimag(b[1][0]);
  b_21r = creal(b[0][1]); b_21i = cimag(b[0][1]);
  b_22r = creal(b[1][1]); b_22i = cimag(b[1][1]);
  
  c_11r = a_11r*b_11r - a_11i*b_11i + a_12r*b_21r - a_12i*b_21i;
  c_11i = a_11i*b_11r + a_11r*b_11i + a_12i*b_21r + a_12r*b_21i;

  c_12r = a_11r*b_12r - a_11i*b_12i + a_12r*b_22r - a_12i*b_22i;
  c_12i = a_11i*b_12r + a_11r*b_12i + a_12i*b_22r + a_12r*b_22i;

  c_21r = a_21r*b_11r - a_21i*b_11i + a_22r*b_21r - a_22i*b_21i;
  c_21i = a_21i*b_11r + a_21r*b_11i + a_22i*b_21r + a_22r*b_21i;

  c_22r = a_21r*b_12r - a_21i*b_12i + a_22r*b_22r - a_22i*b_22i;
  c_22i = a_21i*b_12r + a_21r*b_12i + a_22i*b_22r + a_22r*b_22i;

  // Now replace
  a[0][0] = c_11r + I*c_11i;
  a[1][0] = c_12r + I*c_12i;
  a[0][1] = c_21r + I*c_21i;
  a[1][1] = c_22r + I*c_22i;
  
}


void sdhdf_display_complex_matrix_2x2(double complex J[2][2])
{
  /*
  printf("%g + i%g         %g + i%g\n",creal(J[0][0]),cimag(J[0][0]),creal(J[1][0]),cimag(J[1][0]));
  printf("%g + i%g         %g + i%g\n",creal(J[0][1]),cimag(J[0][1]),creal(J[1][1]),cimag(J[1][1]));
  */

  printf("%g + i%g         %g + i%g      %g + i%g      %g + i%g\n",creal(J[0][0]),cimag(J[0][0]),creal(J[1][0]),cimag(J[1][0]),creal(J[0][1]),cimag(J[0][1]),creal(J[1][1]),cimag(J[1][1]));
}

void sdhdf_setIdentity_4x4(float mat[4][4])
{
  int i,j;
  for (i=0;i<4;i++)
    {
      for (j=0;j<4;j++)
	mat[i][j]=0;
    }
  mat[0][0]=mat[1][1]=mat[2][2]=mat[3][3]=1;    
}
void sdhdf_copy_mat4(float in[4][4],float out[4][4])
{
  int i,j;
  for (i=0;i<4;i++)
    {
      for (j=0;j<4;j++)
	out[i][j] = in[i][j];
    }
}


void sdhdf_copy_vec4(float *in,float *out)
{
  int i;
  for (i=0;i<4;i++)
    out[i]=in[i];
}

void sdhdf_display_vec4(float *vec)
{
  printf("%g\n",vec[0]);
  printf("%g\n",vec[1]);
  printf("%g\n",vec[2]);
  printf("%g\n",vec[3]);
}

void sdhdf_multMat_vec_replace(float mat[4][4],float *vec)
{
  float newVec[4];
  int i;
  
  newVec[0] = mat[0][0]*vec[0] + mat[1][0]*vec[1] + mat[2][0]*vec[2] + mat[3][0]*vec[3];
  newVec[1] = mat[0][1]*vec[0] + mat[1][1]*vec[1] + mat[2][1]*vec[2] + mat[3][1]*vec[3];
  newVec[2] = mat[0][2]*vec[0] + mat[1][2]*vec[1] + mat[2][2]*vec[2] + mat[3][2]*vec[3];
  newVec[3] = mat[0][3]*vec[0] + mat[1][3]*vec[1] + mat[2][3]*vec[2] + mat[3][3]*vec[3];
  for (i=0;i<4;i++)
    vec[i] = newVec[i];
}

void sdhdf_mult4x4_replace(float src1[4][4], float src2[4][4])
{
  float temp[4][4];
  int i,j;
  sdhdf_mult4x4(src1,src2,temp);
  for (i=0;i<4;i++)
    {
      for (j=0;j<4;j++)
	src1[i][j] = temp[i][j];
    }
}

// gamma in radians **
void sdhdf_setFeed(float mf[4][4],float gamma)
{
  sdhdf_setIdentity_4x4(mf);

  // This is from Equation 25 of Robishaw, which is different from the Handbook of pulsar astronomy
  mf[1][1] = mf[2][2] = cos(2*gamma);
  mf[2][1] = sin(2*gamma);
  mf[1][2] = -sin(2*gamma);

}



// Delta G and Delta Phi (radians)
void sdhdf_setGainPhase(float ma[4][4],float diffGain,float diffPhase)
{
  sdhdf_setIdentity_4x4(ma);
  
  ma[1][0] = ma[0][1] = diffGain/2;
  ma[2][2] = ma[3][3] = cos(diffPhase);
  ma[3][2] = -sin(diffPhase);
  ma[2][3] = sin(diffPhase);
  
}

// Delta G/2 and Delta Phi (radians)
void sdhdf_setGain2Phase(float ma[4][4],float diffGain2,float diffPhase)
{
  sdhdf_setIdentity_4x4(ma);
  
  ma[1][0] = ma[0][1] = diffGain2;
  ma[2][2] = ma[3][3] = cos(diffPhase);
  ma[3][2] = -sin(diffPhase);
  ma[2][3] = sin(diffPhase);
  
}


// pa in radians **
void sdhdf_setParallacticAngle(float msky[4][4],float pa)
{
  sdhdf_setIdentity_4x4(msky);
  msky[1][1] = msky[2][2] = cos(2*pa);
  msky[2][1] = sin(2*pa);
  msky[1][2] = -sin(2*pa);
       
}

// Multiply 4x4 matrix
void sdhdf_mult4x4(float src1[4][4], float src2[4][4], float dest[4][4])
{

  dest[0][0] = src1[0][0]*src2[0][0] + src1[1][0]*src2[0][1] + src1[2][0]*src2[0][2] + src1[3][0]*src2[0][3];
  dest[0][1] = src1[0][1]*src2[0][0] + src1[1][1]*src2[0][1] + src1[2][1]*src2[0][2] + src1[3][1]*src2[0][3];
  dest[0][2] = src1[0][2]*src2[0][0] + src1[1][2]*src2[0][1] + src1[2][2]*src2[0][2] + src1[3][2]*src2[0][3];
  dest[0][3] = src1[0][3]*src2[0][0] + src1[1][3]*src2[0][1] + src1[2][3]*src2[0][2] + src1[3][3]*src2[0][3];

  dest[1][0] = src1[0][0]*src2[1][0] + src1[1][0]*src2[1][1] + src1[2][0]*src2[1][2] + src1[3][0]*src2[1][3];
  dest[1][1] = src1[0][1]*src2[1][0] + src1[1][1]*src2[1][1] + src1[2][1]*src2[1][2] + src1[3][1]*src2[1][3];
  dest[1][2] = src1[0][2]*src2[1][0] + src1[1][2]*src2[1][1] + src1[2][2]*src2[1][2] + src1[3][2]*src2[1][3];
  dest[1][3] = src1[0][3]*src2[1][0] + src1[1][3]*src2[1][1] + src1[2][3]*src2[1][2] + src1[3][3]*src2[1][3];

  dest[2][0] = src1[0][0]*src2[2][0] + src1[1][0]*src2[2][1] + src1[2][0]*src2[2][2] + src1[3][0]*src2[2][3];
  dest[2][1] = src1[0][1]*src2[2][0] + src1[1][1]*src2[2][1] + src1[2][1]*src2[2][2] + src1[3][1]*src2[2][3];
  dest[2][2] = src1[0][2]*src2[2][0] + src1[1][2]*src2[2][1] + src1[2][2]*src2[2][2] + src1[3][2]*src2[2][3];
  dest[2][3] = src1[0][3]*src2[2][0] + src1[1][3]*src2[2][1] + src1[2][3]*src2[2][2] + src1[3][3]*src2[2][3];

  dest[3][0] = src1[0][0]*src2[3][0] + src1[1][0]*src2[3][1] + src1[2][0]*src2[3][2] + src1[3][0]*src2[3][3];
  dest[3][1] = src1[0][1]*src2[3][0] + src1[1][1]*src2[3][1] + src1[2][1]*src2[3][2] + src1[3][1]*src2[3][3];
  dest[3][2] = src1[0][2]*src2[3][0] + src1[1][2]*src2[3][1] + src1[2][2]*src2[3][2] + src1[3][2]*src2[3][3];
  dest[3][3] = src1[0][3]*src2[3][0] + src1[1][3]*src2[3][1] + src1[2][3]*src2[3][2] + src1[3][3]*src2[3][3];

  /* WRONG BELOW ***

  dest[0][0] = src1[0][0] * src2[0][0] + src1[0][1] * src2[1][0] + src1[0][2] * src2[2][0] + src1[0][3] * src2[3][0];
  dest[0][1] = src1[0][0] * src2[0][1] + src1[0][1] * src2[1][1] + src1[0][2] * src2[2][1] + src1[0][3] * src2[3][1];
  dest[0][2] = src1[0][0] * src2[0][2] + src1[0][1] * src2[1][2] + src1[0][2] * src2[2][2] + src1[0][3] * src2[3][2];
  dest[0][3] = src1[0][0] * src2[0][3] + src1[0][1] * src2[1][3] + src1[0][2] * src2[2][3] + src1[0][3] * src2[3][3];
  dest[1][0] = src1[1][0] * src2[0][0] + src1[1][1] * src2[1][0] + src1[1][2] * src2[2][0] + src1[1][3] * src2[3][0];
  dest[1][1] = src1[1][0] * src2[0][1] + src1[1][1] * src2[1][1] + src1[1][2] * src2[2][1] + src1[1][3] * src2[3][1];
  dest[1][2] = src1[1][0] * src2[0][2] + src1[1][1] * src2[1][2] + src1[1][2] * src2[2][2] + src1[1][3] * src2[3][2];
  dest[1][3] = src1[1][0] * src2[0][3] + src1[1][1] * src2[1][3] + src1[1][2] * src2[2][3] + src1[1][3] * src2[3][3];
  dest[2][0] = src1[2][0] * src2[0][0] + src1[2][1] * src2[1][0] + src1[2][2] * src2[2][0] + src1[2][3] * src2[3][0];
  dest[2][1] = src1[2][0] * src2[0][1] + src1[2][1] * src2[1][1] + src1[2][2] * src2[2][1] + src1[2][3] * src2[3][1];
  dest[2][2] = src1[2][0] * src2[0][2] + src1[2][1] * src2[1][2] + src1[2][2] * src2[2][2] + src1[2][3] * src2[3][2];
  dest[2][3] = src1[2][0] * src2[0][3] + src1[2][1] * src2[1][3] + src1[2][2] * src2[2][3] + src1[2][3] * src2[3][3];
  dest[3][0] = src1[3][0] * src2[0][0] + src1[3][1] * src2[1][0] + src1[3][2] * src2[2][0] + src1[3][3] * src2[3][0];
  dest[3][1] = src1[3][0] * src2[0][1] + src1[3][1] * src2[1][1] + src1[3][2] * src2[2][1] + src1[3][3] * src2[3][1];
  dest[3][2] = src1[3][0] * src2[0][2] + src1[3][1] * src2[1][2] + src1[3][2] * src2[2][2] + src1[3][3] * src2[3][2];
  dest[3][3] = src1[3][0] * src2[0][3] + src1[3][1] * src2[1][3] + src1[3][2] * src2[2][3] + src1[3][3] * src2[3][3];
  */
}


// Invert 4x4 matrix
int sdhdf_inv4x4(float m[4][4],float invOut[4][4])
{
  float det;
  float inv[4][4];
  int i,j;

  inv[0][0] = m[1][1]  * m[2][2] * m[3][3] -
    m[1][1]  * m[2][3] * m[3][2] -
    m[2][1]  * m[1][2]  * m[3][3] +
    m[2][1]  * m[1][3]  * m[3][2] +
    m[3][1]  * m[1][2]  * m[2][3] -
    m[3][1]  * m[1][3]  * m[2][2];
  
  inv[1][0] = -m[1][0]  * m[2][2] * m[3][3] +
    m[1][0]  * m[2][3] * m[3][2] +
    m[2][0]  * m[1][2]  * m[3][3] -
    m[2][0]  * m[1][3]  * m[3][2] -
    m[3][0] * m[1][2]  * m[2][3] +
    m[3][0] * m[1][3]  * m[2][2];
  
  inv[2][0] = m[1][0]  * m[2][1] * m[3][3] -
    m[1][0]  * m[2][3] * m[3][1] -
    m[2][0]  * m[1][1] * m[3][3] +
    m[2][0]  * m[1][3] * m[3][1] +
    m[3][0] * m[1][1] * m[2][3] -
    m[3][0] * m[1][3] * m[2][1];
  
  inv[3][0] = -m[1][0]  * m[2][1] * m[3][2] +
    m[1][0]  * m[2][2] * m[3][1] +
    m[2][0]  * m[1][1] * m[3][2] -
    m[2][0]  * m[1][2] * m[3][1] -
    m[3][0] * m[1][1] * m[2][2] +
    m[3][0] * m[1][2] * m[2][1];
  
  inv[0][1] = -m[0][1]  * m[2][2] * m[3][3] +
    m[0][1]  * m[2][3] * m[3][2] +
    m[2][1]  * m[0][2] * m[3][3] -
    m[2][1]  * m[0][3] * m[3][2] -
    m[3][1] * m[0][2] * m[2][3] +
    m[3][1] * m[0][3] * m[2][2];
  
  inv[1][1] = m[0][0]  * m[2][2] * m[3][3] -
    m[0][0]  * m[2][3] * m[3][2] -
    m[2][0]  * m[0][2] * m[3][3] +
    m[2][0]  * m[0][3] * m[3][2] +
    m[3][0] * m[0][2] * m[2][3] -
    m[3][0] * m[0][3] * m[2][2];
  
  inv[2][1] = -m[0][0]  * m[2][1] * m[3][3] +
    m[0][0]  * m[2][3] * m[3][1] +
    m[2][0]  * m[0][1] * m[3][3] -
    m[2][0]  * m[0][3] * m[3][1] -
    m[3][0] * m[0][1] * m[2][3] +
    m[3][0] * m[0][3] * m[2][1];
  
  inv[3][1] = m[0][0]  * m[2][1] * m[3][2] -
    m[0][0]  * m[2][2] * m[3][1] -
    m[2][0]  * m[0][1] * m[3][2] +
    m[2][0]  * m[0][2] * m[3][1] +
    m[3][0] * m[0][1] * m[2][2] -
    m[3][0] * m[0][2] * m[2][1];
  
  inv[0][2] = m[0][1]  * m[1][2] * m[3][3] -
    m[0][1]  * m[1][3] * m[3][2] -
    m[1][1]  * m[0][2] * m[3][3] +
    m[1][1]  * m[0][3] * m[3][2] +
    m[3][1] * m[0][2] * m[1][3] -
    m[3][1] * m[0][3] * m[1][2];
  
  inv[1][2] = -m[0][0]  * m[1][2] * m[3][3] +
    m[0][0]  * m[1][3] * m[3][2] +
    m[1][0]  * m[0][2] * m[3][3] -
    m[1][0]  * m[0][3] * m[3][2] -
    m[3][0] * m[0][2] * m[1][3] +
    m[3][0] * m[0][3] * m[1][2];
  
  inv[2][2] = m[0][0]  * m[1][1] * m[3][3] -
    m[0][0]  * m[1][3] * m[3][1] -
    m[1][0]  * m[0][1] * m[3][3] +
    m[1][0]  * m[0][3] * m[3][1] +
    m[3][0] * m[0][1] * m[1][3] -
    m[3][0] * m[0][3] * m[1][1];
  
  inv[3][2] = -m[0][0]  * m[1][1] * m[3][2] +
    m[0][0]  * m[1][2] * m[3][1] +
    m[1][0]  * m[0][1] * m[3][2] -
    m[1][0]  * m[0][2] * m[3][1] -
    m[3][0] * m[0][1] * m[1][2] +
    m[3][0] * m[0][2] * m[1][1];
  
  inv[0][3] = -m[0][1] * m[1][2] * m[2][3] +
    m[0][1] * m[1][3] * m[2][2] +
    m[1][1] * m[0][2] * m[2][3] -
    m[1][1] * m[0][3] * m[2][2] -
    m[2][1] * m[0][2] * m[1][3] +
    m[2][1] * m[0][3] * m[1][2];
  
  inv[1][3] = m[0][0] * m[1][2] * m[2][3] -
    m[0][0] * m[1][3] * m[2][2] -
    m[1][0] * m[0][2] * m[2][3] +
    m[1][0] * m[0][3] * m[2][2] +
    m[2][0] * m[0][2] * m[1][3] -
    m[2][0] * m[0][3] * m[1][2];
  
  inv[2][3] = -m[0][0] * m[1][1] * m[2][3] +
    m[0][0] * m[1][3] * m[2][1] +
    m[1][0] * m[0][1] * m[2][3] -
    m[1][0] * m[0][3] * m[2][1] -
    m[2][0] * m[0][1] * m[1][3] +
    m[2][0] * m[0][3] * m[1][1];
  
  inv[3][3] = m[0][0] * m[1][1] * m[2][2] -
    m[0][0] * m[1][2] * m[2][1] -
    m[1][0] * m[0][1] * m[2][2] +
    m[1][0] * m[0][2] * m[2][1] +
    m[2][0] * m[0][1] * m[1][2] -
    m[2][0] * m[0][2] * m[1][1];
  
  det = m[0][0] * inv[0][0] + m[0][1] * inv[1][0] + m[0][2] * inv[2][0] + m[0][3] * inv[3][0];
  //  printf("det = %g\n",det);
  if (det == 0)
    return -1;

  det = 1.0 / det;
  
  for (i = 0; i < 4; i++)
    {
      for (j=0;j<4;j++)
	{
	  invOut[i][j] = inv[i][j] * det;
	  //	  printf("output: %g\n",invOut[i][j]);
	}
    }
  return 0;
}

/* Computes the dot product of two vectors */
double sdhdf_dotproduct(double *v1,double *v2)
{
  double dot;
  dot = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
  return dot;
}

// Function based on spot.py 
// FIX ME: MUST CHECK THIS CAREFULLY
void sdhdf_para(double dxd,double ddc,double q,double *axd,double *eld)
{
  // Determine position of beam in (RA,DEC) frame
  // given its RA (dxd) and DEC (ddc) offsets, plus
  // parallactic angle (q).
  q = q * M_PI/180.0;
  *axd = -dxd*cos(q) + ddc*sin(q);
  *eld =  ddc*cos(q) + dxd*sin(q);
  
}


void displayMatrix_4x4(float matrix[4][4])
{
  printf("%g %g %g %g\n",matrix[0][0],matrix[1][0],matrix[2][0],matrix[3][0]);
  printf("%g %g %g %g\n",matrix[0][1],matrix[1][1],matrix[2][1],matrix[3][1]);
  printf("%g %g %g %g\n",matrix[0][2],matrix[1][2],matrix[2][2],matrix[3][2]);
  printf("%g %g %g %g\n",matrix[0][3],matrix[1][3],matrix[2][3],matrix[3][3]);
}
