#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <gsl/gsl_integration.h>
using namespace std;

//////////////////////////
//
// Submitted to the Elog by Jonathan Mulholland May 6th 2013
//
// g++ compilation line: g++ -o coils.exe coils.cxx -lgsl -lgslcblas -lm
// no dependency files necessary
// must have installed the gnu scientific libraries

///////////////////////////////////////
//
//These are the functions that are the integrands of the BiotSavart law. 
//The variable being integrated over is theta.
// The point of inquiry and the loop radius are passed as a params[] as require by
// the gsl integration routine that uses bx, by, and bz.
//
double bx (double theta, void * params){
  double x = *((double *) params);
  double y = *((double *) params+1);
  double z = *((double *) params+2);
  double R = *((double *) params+3);
  double bx = -R/4.0/3.141593*z*cos(theta)/pow(gsl_hypot3(x-R*cos(theta),y-R*sin(theta),z),3);
  return bx;
}

double by (double theta, void * params){
  double x = *((double *) params);
  double y = *((double *) params+1);
  double z = *((double *) params+2);
  double R = *((double *) params+3);
  double by = -R/4.0/3.141593*z*sin(theta)/pow(gsl_hypot3(x-R*cos(theta),y-R*sin(theta),z),3);
  return by;
}

double bz (double theta, void * params){
  double x = *((double *) params);
  double y = *((double *) params+1);
  double z = *((double *) params+2);
  double R = *((double *) params+3);
  double bz = R/4.0/3.141593 * (-y*sin(theta) - x*cos(theta) + R )/pow(gsl_hypot3(x-R*cos(theta),y-R*sin(theta),z),3);
  return bz;
}
//
//////////////////////////////////////


int main(){                                                         //main 
  
  ///////////////////////////////////////////////////////
  //  Some variables for the integration routine
  //
  //
  double mu0 = 1.25663706e-6;               // permeability of free space
  double d2r = 3.141593/180.0;              // convert degrees to radians
  double position[] = {0.00,0.00,0.00,0.0}; // (x,y,z,r) position of field point and loop radius passed to integrand
  double resultx,resulty,resultz;           // The result of integration over one segment of one coil
  double errorx,errory,errorz;              // error on said integration
  double bxSum,bySum,bzSum;                 // Running sum of the field contributions from each segment
  double bxErrSum,byErrSum,bzErrSum;        // Running sum of errors on those integrations
 
  ///////////////////////////////////////////////////////////////////////////
  //  The size, position, and orientation of each coil tray
  //
  //  the coordinates I use here: the beam travels in the negative z direction. The magnet bends in the positive y direction
  //  the coordinate system is right handed
  
  double angleL[11];                        // tilt of the coil (degrees) (0 degree aligns with the z-axis)
  double xL[11],yL[11],zL[11];              // the middle of the centeral axis of the coil (meters)
  double lengthL[11],thickL[11],inRL[11];   // length of the coil, thickness of the winding, inner radius of the coil (meters)
  double densityL[11];                      // current density: amps per m^2
  //  For breaking each coil into loops for the BiotSavart treatment
  int divL[11] = {0,0,0,0,0,0,0,0,0,0,0};   // how many segments to break the length of the coil into 
  int divR[11] = {0,0,0,0,0,0,0,0,0,0,0};   // how many segments to break the thickness of the winding into
  double current;                           // magnet current (Amps)

  //coil 1 info
  current = 120.0;
  angleL[0] = 0.0*d2r;
  xL[0] = 0.0;
  yL[0] = 0.0;
  zL[0] = 0.03;
  lengthL[0] = 0.06;
  thickL[0] = 0.0128;
  inRL[0] = 0.07;
  divL[0] = 3;
  divR[0] = 2;
  densityL[0] = current*1422/lengthL[0]/thickL[0]; //current density
  
  //coil 2 info
  angleL[1] = 0.0*d2r;
  xL[1] = 0.0;
  yL[1] = 0.0;
  zL[1] = 0.03;
  lengthL[1] = 0.06;
  thickL[1] = 0.01965;
  inRL[1] = 0.0832;
  divL[1] = 3;
  divR[1] = 2;
  densityL[1] = current*1949/lengthL[1]/thickL[1]; //current density

  //coil 3 info
  angleL[2] = 0.0*d2r;
  xL[2] = 0.0;
  yL[2] = 0.0;
  zL[2] = 0.215;
  lengthL[2] = 0.3;
  thickL[2] = 0.00885;
  inRL[2] = 0.07;
  divL[2] = 15;
  divR[2] = 2;
  densityL[2] = current*4787/lengthL[2]/thickL[2]; //current density

  //coil 4 info
  angleL[3] = 0.0*d2r;
  xL[3] = 0.0;
  yL[3] = 0.0;
  zL[3] = 0.215;
  lengthL[3] = 0.3;
  thickL[3] = 0.0105;
  inRL[3] = 0.0791;
  divL[3] = 15;
  divR[3] = 2;
  densityL[3] = current*5252/lengthL[3]/thickL[3]; //current density

  //coil 5 info
  angleL[4] = 0.0*d2r;
  xL[4] = 0.0;
  yL[4] = 0.0;
  zL[4] = 0.3775;
  lengthL[4] = 0.015;
  thickL[4] = 0.0057;
  inRL[4] = 0.07;
  divL[4] = 2;
  divR[4] = 2;
  densityL[4] = current*152/lengthL[4]/thickL[4]; //current density

  //coil 6 info
  angleL[5] = 0.0*d2r;
  xL[5] = 0.0;
  yL[5] = 0.0;
  zL[5] = 0.3775;
  lengthL[5] = 0.015;
  thickL[5] = 0.01672;
  inRL[5] = 0.07605;
  divL[5] = 2;
  divR[5] = 8;
  densityL[5] = current*396/lengthL[5]/thickL[5]; //current density

  //coil 7 info
  angleL[6] = 0.0*d2r;
  xL[6] = 0.0;
  yL[6] = 0.0;
  zL[6] = 0.3775;
  lengthL[6] = 0.015;
  thickL[6] = 0.01143;
  inRL[6] = 0.1046;
  divL[6] = 2;
  divR[6] = 5;
  densityL[6] = current*303/lengthL[6]/thickL[6]; //current density
  
  //coil 8 info
  angleL[7] = -4.75*d2r;
  xL[7] = 0.0;
  yL[7] = 0.0;
  zL[7] = 0.4235;
  lengthL[7] = 0.035;
  thickL[7] = 0.0315;
  inRL[7] = 0.07;
  divL[7] = 4;
  divR[7] = 15;
  densityL[7] = current*2008/lengthL[7]/thickL[7]; //current density

  //coil 9 info
  angleL[8] = -9.5*d2r;
  xL[8] = 0.0;
  yL[8] = 0.0072;
  zL[8] = 0.4699;
  lengthL[8] = 0.015;
  thickL[8] = 0.03756;
  inRL[8] = 0.07;
  divL[8] = 2;
  divR[8] = 19;
  densityL[8] = current*987/lengthL[8]/thickL[8]; //current density

  //coil 10 info
  angleL[9] = -9.5*d2r;
  xL[9] = 0.0;
  yL[9] = 0.0237;
  zL[9] = 0.5685;
  lengthL[9] = 0.1750;
  thickL[9] = 0.0158;
  inRL[9] = 0.07;
  divL[9] = 18;
  divR[9] = 8;
  densityL[9] = current*5125/lengthL[9]/thickL[9]; //current density

  //coil 11 info
  angleL[10] = -9.5*d2r;
  xL[10] = 0.0;
  yL[10] = 0.0445;
  zL[10] = 0.6931;
  lengthL[10] = 0.0675;
  thickL[10] = 0.0252;
  inRL[10] = 0.08;
  divL[10] = 7;
  divR[10] = 12;
  densityL[10] = current*3204/lengthL[10]/thickL[10]; //current density
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////
  //
  // Setting up the array of points for which we want to calculate the field
  // this program just uses nested loops to hit all combinations of the entries
  // in the x,y,z arrays
  //
  double xPoints[451];  
  double yPoints[451];
  double zPoints[451];
  
  //number of cells in the x,y,z directions
  int numX = 100;
  int numY = 150;
  int numZ = 450;

  double widthX = 0.002;  // distance in m between measured points in x 
  double widthY = 0.002;  // distance in m between measured points in y 
  double widthZ = 0.002;  // distance in m between measured points in z
  //filling the position arrays
  xPoints[0] = -0.1;
  for(int i=1; i< numX; i++){
    xPoints[i] = xPoints[i-1]+widthX;
  }
  yPoints[0] = -0.1;
  for(int i=1; i< numY; i++){
    yPoints[i] = yPoints[i-1]+widthY;
  }
  zPoints[0] = -0.1;
  for(int i=1; i< numZ; i++){
    zPoints[i] = zPoints[i-1]+widthZ;
  }
  //
  //
  //////////////////////////////////////
  
  //////////////////////////////////////
  //
  //prepping for integrating... need a workspace and output files
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc(1000);
  FILE * xfile ;
  // xfile = fopen("coilsx.dat","w");  
  FILE * yfile ;
  //  yfile = fopen("coilsy.dat","w");  
  FILE * zfile ;
  zfile = fopen("coilsz.dat","w"); 
  //  FILE * centerfile ;
  //  centerfile = fopen("center.dat","w"); 
  //  FILE * axfile ;
  //  axfile = fopen("coilsax.dat","w"); 
  //  FILE * testfile ;
  //  testfile = fopen("test.dat","w"); 
  //
  //////////////////////

  //////////////////////
  ///
  ///            Looping over measurement points
  for(int iX = 0; iX<numX;iX++ ){	
    for(int iY = 0; iY<numY;iY++ ){
      for(int iZ = 0; iZ<numZ;iZ++ ){

	// initialize fields and errors. 
	double dI;
	bxSum = 0.0;
	bySum = 0.0;
	bzSum = 0.0;
	bxErrSum = 0.0;
	byErrSum = 0.0;
	bzErrSum = 0.0;

	for(int i=0; i<11; i++){ // loop over coils--the code works just fine if you want change this loop to check a particular coil
	  double dL,dxL,dyL,dzL,dR;
	  double x0,y0,z0,R0;
	  // starting position for the end of the coil (first element)
	  x0 = xL[i];
	  y0 = yL[i]-0.5*lengthL[i]*(1-1/(1.0*divL[i]))*sin(-angleL[i]);
	  z0 = zL[i]-0.5*lengthL[i]*(1-1/(1.0*divL[i]))*cos(-angleL[i]);
	  R0 = inRL[i]+0.5*thickL[i]/(1.0*divR[i]);
	  // steps for going through elements
	  dL = lengthL[i]/divL[i];  // length along axis of coil
	  dxL = 0.0; // step along axis of coil in x direction
	  dyL = dL*sin(-angleL[i]); // step along axis of coil in y direction
	  dzL = dL*cos(-angleL[i]); // step along axis of coil in z direction
	  dR = thickL[i]/divR[i];  // step along coil radius
	  dI = densityL[i]*dL*dR;
	  for(int iL=0; iL<divL[i]; iL++){ // loop over elements (lengthwise and radial)
	    for(int iR=0; iR<divR[i]; iR++){
	      double xElem,yElem,zElem; // element position
	      double R,theta;
	      xElem = x0+dxL*iL;
	      yElem = y0+dyL*iL;
	      zElem = z0+dzL*iL;
	      R = R0+dR*iR;
	      //rotating the positions into the element's coordinate system
	      position[0] = xPoints[iX];
	      position[1] = (yPoints[iY] - yElem)*cos(angleL[i])+ (zPoints[iZ]-zElem)*sin(angleL[i]);
	      position[2] = -(yPoints[iY] - yElem)*sin(angleL[i])+ (zPoints[iZ]-zElem)*cos(angleL[i]);
	      position[3] = R;
	      if(! // stopping integration if the point is too close to the loop element
		 (abs( (R-gsl_hypot(position[0],position[1]))/R ) < 0.10   && abs(position[2]) < 0.001)){
		////////Integration///////////////////////////////////////////
		////Integrating x component over theta 0..2pi
		gsl_function Fx;
		Fx.function = &bx;
		Fx.params = &position;
		gsl_integration_qag(&Fx,0.0,2*3.141593,0.001,0.0,1000,6,w,&resultx,&errorx);      
		//Integrating y component over theta 0..2pi
		gsl_function Fy;
		Fy.function = &by;
		Fy.params = &position;
		gsl_integration_qag(&Fy,0.0,2*3.141593,0.001,0.0,1000,6,w,&resulty,&errory);		
		//Integrating z component over theta 0..2pi
		gsl_function Fz;
		Fz.function = &bz;
		Fz.params = &position;
		gsl_integration_qag(&Fz,0.0,2*3.141593,0.001,0.0,1000,6,w,&resultz,&errorz);
		//
		//  done integration
		//
		//////////////////////////////////////
		//   rotating the b vector back into the original coordinate fram
		//   B' -> B : (Bx',By'*cos(-theta) + Bz'*sin(-theta), -By'*sin(-theta) + Bz*cos(-theta))
		//
		resultx = resultx;
		double tempRY = resulty;
		double tempRZ = resultz;
		resulty = tempRY*cos(-angleL[i])+tempRZ*sin(-angleL[i]);
		resultz = -tempRY*sin(-angleL[i])+tempRZ*cos(-angleL[i]);
		errorx = errorx;
		tempRY = errory;
		tempRZ = errorz;
		errory = tempRY*cos(-angleL[i])+tempRZ*sin(-angleL[i]);
		errorz = -tempRY*sin(-angleL[i])+tempRZ*cos(-angleL[i]);
	      }else{   // counting point as having zero contribution from the coil, if the point is in the coil segment being integrated over
		resultx = 0.0;
		resulty = 0.0;
		resultz = 0.0;
		errorx = 0.0;
		errory = 0.0;
		errorz = 0.0;
	      }

	      // summing from different segements

	      bxSum = bxSum + dI*mu0*resultx;
	      bySum = bySum + dI*mu0*resulty;
	      bzSum = bzSum + dI*mu0*resultz;
	      bxErrSum = bxErrSum + errorx;
	      byErrSum = byErrSum + errory;
	      bzErrSum = bzErrSum + errorz;

	    }  // loop over the elements along the radius of the tray
	  }  // loop over elements along the length of the tray
	}  // loop over coil trays

	// Output for this measurement point
	//	fprintf (xfile,"%f %f\n",zPoints[iZ], bxSum);
	//	fprintf (yfile,"%f %f\n",zPoints[iZ], bySum);
	fprintf (zfile,"%f %f %f %f %f %f\n",xPoints[iX], yPoints[iY], zPoints[iZ], bxSum, bySum, bzSum);
	double mag = gsl_hypot3(bxSum,bySum,bzSum);
	double phase = atan(bySum/bzSum);  
	////////
	//
	//  output for an angled path through the magnet
	//	if(zPoints[iZ]>0.427){ //0.427 -> 0.5
	//	  double z = 0.427+gsl_hypot(zPoints[iZ]-0.427,yPoints[iZ]); 
	//	  fprintf (axfile,"%f %f %f %f \n", z , bzSum*cos(9.5*3.14159/180)+bySum*sin(9.5*3.14159/180), mag, phase);
	//	}else{
	//	  fprintf (axfile,"%f %f %f %f \n", zPoints[iZ] , bzSum, mag, phase);
	//	}

      }  // loop over points for field calculations
      fprintf(zfile,"\n");
    }  // loop over points for field calculations
  }  // loop over points for field calculations
  gsl_integration_workspace_free (w); //free workspace
  //  fclose(xfile);                      // close output file
  //  fclose(yfile);                      // close output file
  fclose(zfile);                      // close output file
  //  fclose(centerfile);
  //  fclose(axfile);                     // close output file 

  return 0;
}
