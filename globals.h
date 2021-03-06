//globals.h

#include "Vector3.h"

#ifndef GLOBAL_H // header guards
#define GLOBAL_H
//#include "logfile.h"


// extern tells the compiler the variable is declared elsewhere

const double pi = 3.14159265358979323;
const double el_charge  = 1.60217733e-19;	//el lading [C]
const double amu = 1.6605402E-27;//massa van een ion in kg
const double ke=8.987551787E9; //=1/(4*pi*eps0) [N*(m/C)^2]
const double gravity=6.67428e-11; // [m^3/(kg*s^2)]
const double Joule_to_eV = 6.24150e18; //1J = 6.24150e18 eV
const double eV_to_Joule = 1.60217646e-19 ;
const double mass_Mev = 931.494;
const double eps0 = 8.85419e-12 ;

//const double pi = acos(-1.0);
const double pi2 = 2 * pi;
const double sqrt2 = sqrt(2);
const double sqrt2o2 = sqrt2/2;

//standard vars concerning the gas
const double kb = 1.3806504E-23;//Boltzmann constant  [J/K]
const double epsilon_0 = 8.8541878176E-12; //[Farad/m]=[C/(V*m)]
//const double He_xi_e = 1.0e-7; //dielectric susceptebility of Helium

const double mass_electron_amu = 0.000548579867;
const double mass_proton_amu = 1.00727638;
const double mass_proton_kg = 1.6726217E-27;

const Vector3 x_hat(1.,0.,0.);
const Vector3 y_hat(0.,1.,0.);
const Vector3 z_hat(0.,0.,1.);
const Vector3 NullVector(0.,0.,0.);

//const double e2m_r = 1.672621777e−27;

/* const double massa_He4_amu = 4.00260325415; */
/* const double massa_Ar40_amu = 39.9623831225; */
/* const double massa_Ar35_amu = 34.9752576;   //34.9752576(8) */
/* const double massa_Cl35_amu = 34.96885268;  //34.96885268(4) */
/* const double massa_K39_amu = 38.96370668;   //38.96370668(20)  */
/* const double massa_Cs133_amu = 132.905451932; //132.905451(932) */
/* const double massa_144_Eu = 143.918816823;  //(30000) */
/* const double massa_144_Gd = 143.922963000;  //(30000) */
/* const double massa_H2O_amu=18.0106; */

/* const double atomic_radii_He4 = 0.31E-10; */
/* const double atomic_radii_Ar40 = 0.71E-10;  //is not ionic radius but the atomic */
/* const double atomic_radii_Ar35_m = 0.88E-10; //is not ionic radius but the atomic */
/* const double atomic_radii_Cl35_m = 0.97E-10; //is not ionic radius but the atomic */
/* const double atomic_radii_K39_m = 2.27E-10; */
/* const double atomic_radii_Eu144_m = 1.85E-10; */
/* const double atomic_radii_Gd144_m = 1.80E-10; */
/* const double atomic_radii_H2O=3.7E-10; */

#endif
