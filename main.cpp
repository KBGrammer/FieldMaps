#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "EFieldMap.h"
#include "EFieldPoint.h"
#include "Vector3.h"
#include "MagFieldMap.h"
#include "ChargedParticle.h"
#include "nanoflann.hpp"
#include "rng.h"
#include "globals.h"
#include <iomanip>

Vector3 Permute_BField(Vector3 vec);
Vector3 Get_BField_Point(Vector3 vec);
void printbmap (int start);
void unit_vector(Rng & random, GRng & gauss,
		 Vector3 & vect, Vector3 & origin,
		 int type);
void write_BField(int i );
Vector3 Get_BField_Divergence(Vector3 point,
			      const EFieldMap & bfield,
			      const my_kd_tree_t & index_b,
			      double search_radius,
			      ChargedParticle p);


int main (int argc, char *argv[])
{
  int start = atoi(argv[1]);
  // printbmap(start);
  // return 0;
 
  // int i = atoi(argv[1]);
  // write_BField(i);
  // return 0;

  EFieldMap efield_map(1.0/512.0);
  // MagFieldMap bfield_map_old;

  // for(double tx = -100; tx < 100; tx+=5) {
  //   for(double ty = -100; ty < 100; ty+=5) {
  //     Vector3 vec(tx, ty, 0);
  //     Vector3 b = bfield_map_old.FieldEval(vec);
  //     std::cout << vec << ' ' << b << std::endl;
  //   }
  // }

  // for(double tx = -25; tx < 25; tx+=0.5) {
  //   for(double ty = -25; ty < 25; ty+=0.5) {
  //     Vector3 vec(tx, ty, 0);
  //     Vector3 b = bfield_map_old.FieldEval(vec);
  //     std::cout << vec << ' ' << b << std::endl;
  //   }
  // }

  // for(double tx = -3; tx < 3; tx+=0.05) {
  //   for(double ty = -3; ty < 3; ty+=0.05) {
  //     Vector3 vec(tx, ty, 0);
  //     Vector3 b = bfield_map_old.FieldEval(vec);
  //     std::cout << vec << ' ' << b << std::endl;
  //   }
  // }

  // return 0;
  //efield_map.ReadEFieldMap(3305196, "efield_sort_closed");
  //  efield_map.ReadEFieldMap(3305196, "efield_trap_open2");
  //  efield_map.ReadEFieldMap(1792199, "trap_open_20130906.dat");
  //efield_map.ReadEFieldMap(4098329, "trap_open_ramp_20130906.dat");
  efield_map.ReadEFieldMap(0, "comsol_efield_L9_ramp.dat");

  EFieldMap bfield_map(1.0/512.0);
  //bfield_map.ReadEFieldMap(3305196, "bfield_sort3");
  //bfield_map.ReadEFieldMap(0, "fieldmap_coarse_sort.dat");
  //  bfield_map.ReadEFieldMap(0, "magfield_highprecision.dat");
  //bfield_map.ReadEFieldMap(0, "Bfield_RDK2_Matt2_meters2.txt");
  bfield_map.ReadEFieldMap(0, "b_field_comsol_3312014_sort.dat");

  
  my_kd_tree_t index_efield(3 /*dim*/, efield_map.map, 
		       KDTreeSingleIndexAdaptorParams(5 /* max leaf */));
  index_efield.buildIndex();
  
  my_kd_tree_t index_bfield(3 /*dim*/, bfield_map.map, 
		       KDTreeSingleIndexAdaptorParams(5 /* max leaf */));
  index_bfield.buildIndex();
  //efield_map.GenerateDensity(0.01, index_efield);
  //bfield_map.GenerateDensity(0.01, index_bfield);
  int rank = 0;

  unsigned long long rngseed = 2385556849720357773ULL;
  unsigned long long grngseed = 18359139775082460051ULL;

  Rng random(rngseed);
  GRng gauss(grngseed);

  Vector3 grad(0.65416, 0.0186, 0.0348949);
  Vector3 origin(0.65416, 0.0186, 0.0348949);
  Vector3 nullvector (0, 0, 0);
  Vector3 trap_axis (1, 0, 0);
  

  double delta = 1.0/8.0/20.0;
  double xlim = delta*20;
  

  Vector3 initialpos(0.12435, 0.0, 0.0);
  
  const double mev_2_joule = 1.0/(6.24150647996E+12);

  bool stop = false;
  const Vector3 detcenter(0.699, 0.041, 0.0);
  const Vector3 normal(0.699-0.447, 0.041, 0);
  //const Vector3 detcenter(0.697045, 0.04324, 0.0);
  //const Vector3 normal(.67952-.699308, .041731-.045043, 0);

  // 679.52, 41.731, 0
  // 699.308, 45.043, 0

  double rkf45_stepsize = 1e-11;
  
  char identify = ' ';
  
  for(int i = 0; i < start; i++) {
    rkf45_stepsize = 1e-10;
    identify = ' ';
    //Vector3 tvect(1/sqrt(3), 1/sqrt(3),1/sqrt(3));

    Vector3 tvect;
    Vector3 origin (0.19,0.0,0);
    unit_vector(random,gauss,tvect,origin,2);
    //Vector3 origin (0.33,0.0,0);
    //unit_vector(random,gauss,tvect,origin,2);
    //tvect = Vector3(0,sqrt(2.0)/2.0, sqrt(2.0)/2.0);
    
    tvect.scalefactor(sqrt(751.0 * eV_to_Joule * 2.0 / mass_proton_kg));

    tvect.scalecomponents(-1,1,1);

    Vector3 B;
    double searchrad = 0.01;
    bool goodB = bfield_map.Interpolate_Field(origin, index_bfield, searchrad, B);
    double frequency = (el_charge / mass_proton_kg) * B.magnitude() / (2.0 * pi);
    //rkf45_stepsize = 1.0 / (frequency * 120.0);
    //    rkf45_stepsize = 1.0 / (frequency * 120.0);

    rkf45_stepsize = 1e-10;
    //std::cerr << rkf45_stepsize << std::endl;

    //0.220398 -0.00253615 -0.00608275 -51370.7 -240562 -289640
    //Vector3 tvect(51370.7, -240562, -289640);
    //Vector3 origin (0.220398, -0.00253615, -0.00608275);


    ChargedParticle proton (origin, tvect, el_charge, 
			    mass_proton_kg, 0, 0, 0.01, 0.001);

    proton.SetTolerance(1e-12);
    
    Vector3 diff;

    int q = 0;


    //    if( i == 7) {
      std::cout << i << ' ' << -2 << ' ' << proton << std::endl;
      std::cerr << i << ' ' << -2 << ' ' << proton << std::endl;
      //    }

    bool pass_mid = false;
    bool good_field = true;
    double x, vx;
    Vector3 Forcesum(0,0,0);
    //    if( i == 7) {
    do {
      // broken 0 0.150398 -0.00253615 -0.00608275 -51370.7 -240562 -289640 1.20763e-16 0
      if(true) {
	good_field = proton.Prop_DT_RKF452(efield_map,
					   bfield_map,
					   index_efield,
					   index_bfield,
					   rkf45_stepsize);
      }
      else if(false) {
	good_field = proton.Prop_DT_VelVer(efield_map,
					   bfield_map,
					   index_efield,
					   index_bfield,
					   rkf45_stepsize);
      }
      else if(false) {
	good_field = proton.Prop_DT_RKKC45(efield_map,
					   bfield_map,
					   index_efield,
					   index_bfield,
					   rkf45_stepsize);
      }
      else if(false) {
	good_field = proton.Prop_DT_RKF451(efield_map,
					   bfield_map,
					   index_efield,
					   index_bfield,
					   rkf45_stepsize);
      }
      else {
	good_field = proton.Prop_DT_DOPRI(efield_map,
					  bfield_map,
					  index_efield,
					  index_bfield,
					  rkf45_stepsize);
      }
      diff = proton.GetPos() - detcenter;
      //Vector3 EF = proton.GetForce_E(proton.GetPos(), proton.GetVel(), 0, efield_map, index_efield);
      //Vector3 BF = proton.GetForce_B(proton.GetPos(), proton.GetVel(), 0, bfield_map, index_bfield);
      //std::cout << proton << ' ' << EF << ' ' << BF << std::endl;

      if(proton.GetTime() > 1e-4) break;
      
      if(efield_map.CheckBoundary(proton.GetPos())) {
      	proton.SetOldPos();
      }
      
      if((q % 50 == 0 && false) && good_field) {
	std::cout.precision(10);
	Vector3 ef, bf, E, B;

	proton.GetForce_E(proton.GetPos(),
			  proton.GetVel(), 0,
			  efield_map, index_efield,ef);
	proton.GetForce_B(proton.GetPos(),
			  proton.GetVel(), 0,
			  bfield_map, index_bfield,bf);
	proton.Get_E(proton.GetPos(),
		     proton.GetVel(), 0,
		     efield_map, index_efield, E);
	proton.Get_B(proton.GetPos(),
		     proton.GetVel(), 0,
		     bfield_map, index_bfield, B);
	Vector3 diver = Get_BField_Divergence(proton.GetPos(), bfield_map,
					      index_bfield, proton.GetSRB(), proton);
	
	stop = (proton.GetVel().dot(B) < 0 && pass_mid);
	std::cout << i << ' ' << q << ' ' << proton << ' ' 
		  << ef << ' '<< bf << ' ' << E << ' ' << B
		  << ' ' << proton.GetVel().dot(B) << ' ' << diver << std::endl;
	
      }
      if(!good_field) std::cerr << "bad field " << i << std::endl;
      Vector3 bf;
      proton.GetForce_B(proton.GetPos(),
			proton.GetVel(), 0,
			bfield_map, index_bfield,bf);
      Forcesum = Forcesum + bf;
      
      x = proton.GetPos().xc();
      vx = proton.GetVel().xc();
      if(x > 0.3 && vx > 0.0) {
        pass_mid = true;
      }

      //stop = (x > 0.3 && vx < -1e5);

      if(x > 0.3 && vx < -1e5) {
	identify = 'f';
      }

      stop = stop && (proton.GetPos().xc() == proton.GetOldPos().xc());

      q++; 
      }while (efield_map.CheckBoundary(proton.GetPos()) && good_field);
      //}while (true);
    // ensures points beyond detector aren't used
    //    if(i == 7) {
      std::cerr << i << ' ' << -1 << ' ' << proton 
		<< ' ' << Forcesum << ' ' << identify << std::endl;
      std::cout << i << ' ' << -1 << ' ' << proton 
		<< ' ' << Forcesum << ' ' << identify << std::endl;
      //    }
  }

  return 0;
}

Vector3 Permute_BField(Vector3 vec)
{
  return Vector3(vec.zc(), vec.yc(), vec.xc());
}

void write_BField(int start)
{
  MagFieldMap bfield;
  if (false) {
    std::ifstream fin;
    fin.open("comsol_efield.dat");
    double t_ex, t_ey, t_ez, t_cx, t_cy, t_cz;
    int i = 0;
    start = start * 80000;
    do {
      //fin >> t_ex >> t_ey >> t_ez >> trash >> trash >> trash
      //	>> t_cx >> t_cy >> t_cz;
      fin >> t_ex >> t_ey >> t_ez
	  >> t_cx >> t_cy >> t_cz;
      if(i >= start) {
	Vector3 vect (t_cx, t_cy, t_cz);
	Vector3 vectb (t_cz-0.427+0.447, t_cy, t_cx);
	Vector3 B = bfield.FieldEval(vectb);
	std::cout << B << ' ' << vect << std::endl;
      }
      i++;
    }while(!fin.eof() && i < start + 80000);
  }
  else if(false) {
    const double xybound = 0.02;
    const double xybins = 40;
    const double dxyz = 0.001;
    const double zbins = 750;
    for(int i = start; i < start+1 && i < zbins; i++) {
      double x = 0.0 + double(i) * dxyz;
      if(x < 0.447) {
	for(int k = 0; k < xybins; k++) {
	  double z = -xybound + double(k) * dxyz;
	  for(int j = 0; j < xybins; j++) {
	    double y = -xybound + double(j) * dxyz; 
	    Vector3 vect (x, y, z);
	    Vector3 vectb (z-0.427+0.447, y, x);
	    Vector3 B = bfield.FieldEval(vectb);
	    std::cout << std::setprecision(15) << B << ' ' << vect << std::endl;
	  }
	}
      }
    }

  }
  else if(true) {
    const double xybound = 0.02;
    const double xybins = 40;
    const double dxyz = 0.001;
    const double zbins = 750;
    for(int i = start; i < start+100 && i < zbins; i++) {
      double x = 0.0 + double(i) * dxyz;
      if(x < 0.447) {
	for(int k = 0; k < xybins; k++) {
	  double z = -xybound + double(k) * dxyz;
	  for(int j = 0; j < xybins; j++) {
	    double y = -xybound + double(j) * dxyz; 
	    Vector3 vect (x, y, z);
	    Vector3 vectb (z-0.427+0.447, y, x);
	    Vector3 B = bfield.FieldEval(vectb);
	    std::cout << std::setprecision(15) << B << ' ' << vect << std::endl;
	  }
	}
      }
      else {
	for(int k = 0; k < xybins; k++) {
	  double z = -xybound + double(k) * dxyz;
	  for(int j = 0; j < xybins; j++) {
	    double y = -xybound + double(j) * dxyz + (0.045/0.25)*(x-0.447); 
	    Vector3 vect (x, y, z);
	    Vector3 vectb (z-0.427+0.447, y, x);
	    Vector3 B = bfield.FieldEval(vectb);
	    std::cout << std::setprecision(15) << B << ' ' << vect << std::endl;
	  }
	}
      }
    }
  }
  else {
    const double xybound = 0.011324;
    const double xybins = 50;
    const double zbins = 350;
    for(int i = start; i < start+50 && i < zbins; i++) {
      double x = 0.29764 + double(i) * 0.2/zbins;
      if(x < 0.447) {
	for(int k = 0; k < xybins; k++) {
	  double z = -xybound + double(k) * 2*xybound/xybins;
	  for(int j = 0; j < xybins; j++) {
	    double y = -xybound + double(j) * 2*xybound/xybins; 
	    Vector3 vect (x, y, z);
	    Vector3 vectb (z-0.427+0.447, y, x);
	    Vector3 B = bfield.FieldEval(vectb);
	    std::cout << B << ' ' << vect << std::endl;
	  }
	}
      }
      else {
	for(int k = 0; k < xybins; k++) {
	  double z = -xybound + double(k) * 2*xybound/xybins;
	  for(int j = 0; j < xybins; j++) {
	    double y = -xybound + double(j) * 2*xybound/xybins + (0.045/0.25)*(x-0.447); 
	    Vector3 vect (x, y, z);
	    Vector3 vectb (z-0.427+0.447, y, x);
	    Vector3 B = bfield.FieldEval(vectb);
	    std::cout << B << ' ' << vect << std::endl;
	  }
	}
      }
    }
  }
}

Vector3 Get_BField_Divergence(Vector3 point,
			      const EFieldMap & bfield,
			      const my_kd_tree_t & index_b,
			      double search_radiusb,
			      ChargedParticle p)
{
  Vector3 field1, field2;
  double h = 0.0001; //search_radiusb;
  Vector3 vect;
  double t;
  //  std::cerr << search_radiusb << std::endl;
  p.Get_B(point - Vector3(h,0,0), vect, t, bfield, index_b, field1);
  p.Get_B(point + Vector3(h,0,0), vect, t, bfield, index_b, field2);
  double delta_x = (field2 - field1).dot(x_hat);

  p.Get_B(point - Vector3(0,h,0), vect, t, bfield, index_b, field1);
  p.Get_B(point + Vector3(0,h,0), vect, t, bfield, index_b, field2);;
  double delta_y = (field2 - field1).dot(y_hat);

  p.Get_B(point - Vector3(0,0,h), vect, t, bfield, index_b, field1);
  p.Get_B(point + Vector3(0,0,h), vect, t, bfield, index_b, field2);
  double delta_z = (field2 - field1).dot(z_hat);

  return Vector3(delta_x, delta_y, delta_z) / (2.0 * h);
}

Vector3 Get_BField_Point(Vector3 vec)
{
  return Vector3(vec.zc()-17.559*.0254+0.427, vec.yc(), vec.xc());
}

void printbmap (int start)
{
  //std::cout << "#Generating " << N << " point cloud..." << std::endl;
  MagFieldMap bfield;

  std::ifstream fin;
  fin.open("efield_sort_open2");
  double t_ex, t_ey, t_ez, t_cx, t_cy, t_cz;
  double trash;

  std::vector <Vector3> points;
  
  int i = 0;
  do {
    fin >> t_ex >> t_ey >> t_ez >> t_cx >> t_cy >> t_cz;
    //map.pts[i].x = t_cx;
    //map.pts[i].y = t_cy;
    //map.pts[i].z = t_cz;
    //map.efpoints[i].x = t_ex;
    //map.efpoints[i].y = t_ey;
    //map.efpoints[i].z = t_ez;

    Vector3 point (t_cx, t_cy, t_cz);
    //Vector3 B = bfield.FieldEval(point);

    points.push_back(point);

    //std::cout << B << ' ' << point << std::endl;

    i++;
  }while(!fin.eof());

  fin.close();

  for(int j = start; j < i+1; j += 6) {
    Vector3 B = bfield.FieldEval(points[j]);
    std::cout << B << ' ' << points[j] << std::endl;
  }

}

void unit_vector(Rng & random, GRng & gauss, Vector3 &vect, Vector3 & origin, int type)
{
  double u,theta,temp,x,y,z,ox,oy,oz;
  u = random.rdouble()*2.0-1.0;
  theta = random.rdouble_exc()*pi2;
  temp = sqrt(1-u*u);
  x = temp*cos(theta);
  y = temp*sin(theta);
  z = u;
  vect.set_coords(x,y,z);
    // double mag1 = vect.magnitude();
  // vect.make_unitvector();
  // double mag2 = vect.magnitude();
  // cout << mag2 - mag1 << huzzah;

  //  double slabsize = 9*2.54;

  ox = origin.xc();
  oy = origin.yc();
  oz = origin.zc();

  if(type == 1) { // Point source at origin
    origin.set_coords(ox,oy,oz);
  }
  else if(type == 2) {// cylinder source
    double tz,ty;
    do{
      tz = (random.rdouble()*2.0-1.0);
      ty = (random.rdouble()*2.0-1.0);
    }while(tz*tz+ty*ty>1);
    tz = tz*0.01; // 4.5*2.54/2 = radius in cm
    ty = ty*0.01;
    double tx = (random.rdouble())*0.05; // width
    origin.set_coords(ox+tx,oy+ty,oz+tz);
  }
  else if(type == 3) { // pill source
    origin.set_coords(ox,oy,oz+double(random.rdouble()*60.0-30.0));
  }
}
