#include <iostream>
#include "Vector3.h"
#include "EFieldMap.h"
#include "EFieldPoint.h"
#include "MagFieldMap.h"
#include "globals.h"

#ifndef CHARGEDPARTICLE_H
#define CHARGEDPARTICLE_H

class ChargedParticle {
 private:

  double c2m_r;

  Vector3 position;
  Vector3 oldposition;
  Vector3 oldaccel;
  Vector3 velocity;
  Vector3 oldvelocity;
  Vector3 oldB;
  double charge;
  double mass;
  double time;
  double CTolerance;

  Vector3 vold_half;
  Vector3 xold;

  double c_dp[14];
  double a_dp[14][13];
  double b_dp1[14];
  double b_dp2[14];

  double c_rkf45[7];
  double a_rkf45[7][6];
  double b_rkf451[7];
  double b_rkf452[7];

  double c_rkkc45[7];
  double a_rkkc45[7][6];
  double b_rkkc451[7];
  double b_rkkc452[7];

  double rkf45_stepsize;
  double search_radiuse;
  double search_radiusb;

  double pot_E;
  double kin_E;

 public:
  ChargedParticle();
  ChargedParticle(const Vector3 & pos,
		  const Vector3 & vel,
		  const double & cha,
		  const double & ma,
		  const double & potent,
		  const double & ti,
		  const double & sre,
		  const double & srb);
  //ChargedParticle& operator=(const ChargedParticle& rhs);
  ~ChargedParticle();
  void DefineConstants();
  void SetTolerance(double tol);

  void Prop_DT(const Vector3 & theForce, double timestep, const my_kd_tree_t & index);
  void Prop_DT(const Vector3 & theForce, double timestep,
	       const Vector3 & workForce, int changes, const my_kd_tree_t & index);
  void Prop_DT_RK4(double t, const EFieldMap & efield,
		   const MagFieldMap & bfield,
		   const my_kd_tree_t & index);
  void Prop_DT_RKF45(const EFieldMap & efield,
		     const MagFieldMap & bfield,
		     const my_kd_tree_t & index,
		     double & stepsize);
  bool Prop_DT_RKF451(const EFieldMap & efield,
		     const EFieldMap & bfield,
		      const my_kd_tree_t & index_e,
		     const my_kd_tree_t & index_b,
		     double & stepsize);
  bool Prop_DT_RKF45_GC(const EFieldMap & efield,
		     const EFieldMap & bfield,
		      const my_kd_tree_t & index_e,
		     const my_kd_tree_t & index_b,
		     double & stepsize);
  bool Prop_DT_RKKC45(const EFieldMap & efield,
		     const EFieldMap & bfield,
		      const my_kd_tree_t & index_e,
		     const my_kd_tree_t & index_b,
		     double & stepsize);
  bool Prop_DT_VelVer(const EFieldMap & efield,
		     const EFieldMap & bfield,
		      const my_kd_tree_t & index_e,
		     const my_kd_tree_t & index_b,
		     double & stepsize);
  bool Prop_DT_RKF452(const EFieldMap & efield,
  		     const EFieldMap & bfield,
		      const my_kd_tree_t & index_e,
  		     const my_kd_tree_t & index_b,
  		     double & stepsize);
  void UpdateVelocityBoris(double dt, Vector3 ef, Vector3 bf);
  bool Prop_DT_DOPRI(const EFieldMap & efield,
		     const EFieldMap & bfield, 
		     const my_kd_tree_t & index_e,
		     const my_kd_tree_t & index_b,
		     double & stepsize);
  void Prop_DT_DOPRI2(const EFieldMap & efield,
		     const EFieldMap & bfield, 
		     const my_kd_tree_t & index_e,
		     const my_kd_tree_t & index_b,
		     double & stepsize);
  Vector3 GetForce(Vector3 point, Vector3 vel, double t,
		   const EFieldMap & efield,
		   const MagFieldMap & bfield, const my_kd_tree_t & index);
  bool GetForce(Vector3 point, Vector3 vel, double t,
		const EFieldMap & efield,
		const EFieldMap & bfield, const my_kd_tree_t & index_e,
		const my_kd_tree_t & index_b, Vector3 & force);
  bool GetForce_E(Vector3 point, Vector3 vel, double t,
		  const EFieldMap & efield,
		  const my_kd_tree_t & index_e,
		  Vector3 & force);
  bool GetForce_B(Vector3 point, Vector3 vel, double t,
		  const EFieldMap & bfield,
		  const my_kd_tree_t & index_b,
		  Vector3 & force);
  bool Get_E(Vector3 point, Vector3 vel, double t,
	     const EFieldMap & efield,
	     const my_kd_tree_t & index_e,
	     Vector3 & force);
  bool Get_B(Vector3 point, Vector3 vel, double t,
	     const EFieldMap & bfield,
	     const my_kd_tree_t & index_b,
	     Vector3 & force);
  void SetOldPos();
  Vector3 GetFakeB(Vector3 point);
  Vector3 GetFakeE(Vector3 point);
  Vector3 GetPos();
  Vector3 GetOldPos();
  Vector3 GetVel();
  Vector3 GetOldVel();
  Vector3 GetOldB();
  double GetPE();
  double GetKE();
  double GetCharge();
  double GetTime();
  double GetSRB();
  double GetSRE();


  friend std::ostream& operator<< (std::ostream &out, ChargedParticle &tChaPar);
  /* friend std::istream& operator>> (std::istream &in, ChargedParticle &tVect); */
  /* friend ChargedParticle operator+(const ChargedParticle &v1, const ChargedParticle &v2); */
  /* friend ChargedParticle operator+(const ChargedParticle &v1, const double & scalar); */
  /* friend ChargedParticle operator-(const ChargedParticle &v1, const ChargedParticle &v2); */
  /* friend ChargedParticle operator-(const ChargedParticle &v1, const double & scalar); */
  /* friend ChargedParticle operator*(const ChargedParticle &v1, const double & scalar); */

};


#endif
