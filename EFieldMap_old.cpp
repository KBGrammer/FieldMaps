#include <vector>
#include <fstream>
#include "EFieldMap.h"
#include "EFieldPoint.h"
#include "Vector3.h"
#include "mymath.h"
#include <cmath>
#include <iostream>

// EFieldMap constructor

// EFieldMap member function
void EFieldMap::ReadEFieldMap()
{

  EFieldPoint my_point_temp;

  bin_spacing = 0.0001;
  double current = -1e15;

  std::ifstream fin;
  fin.open("temp_efield");
  double t_ex, t_ey, t_ez, t_cx, t_cy, t_cz;
  double trash;
  
  Vector3 t_c, t_ef;
  int counter = 0;

  do {
    fin >> t_ex >> t_ey >> t_ez >> trash >> trash >> trash
	>> t_cx >> t_cy >> t_cz;
    
    if(current == -1e15 || current + bin_spacing < t_cx) {
      current = t_cx;
      bins.push_back(counter);
      keys.push_back(current);
    }

    // if (counter < 10) {
    //   std::cout << t_ex << ' ' << t_ey << ' ' << t_ez << ' ' << trash << ' ' << trash
    // 		<< ' ' << t_cx << ' ' << t_cy << ' ' << t_cz << std::endl;
    // }

    counter++;

    t_c.set_coords(t_cx, t_cy, t_cz);
    t_ef.set_coords(t_ex, t_ey, t_ez);

    my_point_temp.SetEFieldPoint(t_c, t_ef);
    map.push_back(my_point_temp);
  }while(!fin.eof());
}

void EFieldMap::PrintMap(int count)
{
  //std::cout << map.size() << std::endl;
  if (count ==  0) count = map.size();
  for (int i = 0; i < count; i++) {
    Vector3 tPoint = map[i].GetPoint();
    std::cout << tPoint << ' ' << std::endl;
  }
}

void EFieldMap::PrintKeysBins()
{
  std::cout << keys.size() << std::endl;
  for (int i = 0; i < keys.size(); i++) {
    std::cout << bins[i] << ' ' << keys[i] << std::endl;
  }
}

void EFieldMap::ListPointsinRad(Vector3 vect, double rad)
{
  for (int i= 0; i < map.size(); i++) {
    double dist = vect.distance(map[i].GetPoint());
    // std::cout << dist << ' ' << vect.xc() << ' ' 
    // 	      << vect.yc() << ' ' << vect.zc() << ' ' 
    // 	      << map[i].GetPoint().xc() << ' ' 
    // 	      << map[i].GetPoint().yc() << ' '
    // 	      << map[i].GetPoint().zc() << ' ' << std::endl;
    if(dist < rad) {
      Vector3 goodpoint = map[i].GetPoint();
      std::cout << goodpoint << ' ' << dist << std::endl;
    }
  }
}

// Inverse Distance Weighting
// Uses weight function = squ((R-d))/(R*d))
Vector3 EFieldMap::Interpolate_Field(Vector3 vect) const
{
  double rad_limit = 0.002; 
  double rad_limit2 = squ(rad_limit);
  double denom_sum = 0;
  double num_sum_ex = 0;
  double num_sum_ey = 0;
  double num_sum_ez = 0;
  int start = 0;
  int end = 0;
  int tstart = 0;

  do {
    tstart++;
  }while(keys[tstart] < vect.xc());
  
  start = tstart - int(rad_limit / bin_spacing) - 1;
  end   = tstart + int(rad_limit / bin_spacing) + 1;

  // std::cout << vect.xc() << ' ' << start << ' '
  // 	    << end << ' ' << keys[tstart]
  // 	    << ' ' << rad_limit / bin_spacing << std::endl;

  if(start < 0) start = 0;
  if(end > keys.size()) end = keys.size()-1;

  // std::cout << vect.xc() << ' ' << start << ' '
  // 	    << end << ' ' << keys[tstart]
  // 	    << ' ' << keys[start] << ' ' << keys[end] << std::endl;

  start = bins[start];
  end = bins[end];

  int counter = 0;
  int retry = 0;

  for(int i = start; i < end; i++) {
    Vector3 tpoint = map[i].GetPoint();
    double dx = map[i].GetPoint().xc()- vect.xc();
    double dy = map[i].GetPoint().yc() - vect.yc();
    double dz = map[i].GetPoint().zc() - vect.zc();
    if((dx < rad_limit && dx > -rad_limit)
       && (dy < rad_limit && dy > -rad_limit)
       && (dz < rad_limit && dz > -rad_limit)) {
      double dist = vect.distancesq(map[i].GetPoint());
      if(dist < rad_limit2) {
	dist = sqrt(dist);
	counter++;
	Vector3 t_efield = map[i].GetField();
	Vector3 t_point = map[i].GetPoint();
	double weight = squ((rad_limit - dist) / (rad_limit * dist));
	denom_sum += weight;
	num_sum_ex += weight * t_efield.xc();
	num_sum_ey += weight * t_efield.yc();
	num_sum_ez += weight * t_efield.zc();
	//std::cout << "0 " << t_point << ' ' << t_efield << ' ' << dist << std::endl;
      }
      if(i == end - 1 && counter < 20 && retry == 0) {
	denom_sum = num_sum_ex = num_sum_ey = num_sum_ez = 0;
	rad_limit *= 2;
	rad_limit2 *= 4;
	start = tstart - int(rad_limit / bin_spacing) - 1;
	end   = tstart + int(rad_limit / bin_spacing) + 1;
	i = start;
	if(start < 0) start = 0;
	if(end > keys.size()) end = keys.size()-1;
	
	start = bins[start];
	end = bins[end];
	
	retry = 1;
      }
    }
  }
  //std::cout << counter << std::endl;
  Vector3 result(num_sum_ex / denom_sum, 
	  num_sum_ey / denom_sum,  
	  num_sum_ez / denom_sum);
  //std::cout << "100 " << vect << ' ' << result << std::endl;
  return result;
}

Vector3 EFieldMap::Interpolate_Field_Gradient(Vector3 vect,
					      Vector3 &gradient)
{
  double rad_limit = 0.005; 
  double rad_limit2 = squ(rad_limit);
  double denom_sum = 0;
  double num_sum_ex = 0;
  double num_sum_ey = 0;
  double num_sum_ez = 0;

  double xdivnum[4][2];
  double ydivnum[4][2];
  double zdivnum[4][2];
  for(int i=0; i<4; i++) {
    xdivnum[i][0] = xdivnum[i][1] = 0;
    ydivnum[i][0] = ydivnum[i][1] = 0;
    zdivnum[i][0] = zdivnum[i][1] = 0;
  }
  int start = 0;
  int end = 0;

  do {
    start++;
  }while(keys[start] < vect.xc());
  
  start = start - 1;
  end   = start + 2;

  if(start < 0) start = 0;
  if(end > keys.size()) end = keys.size()-1;

  start = bins[start];
  end = bins[end];

  for(int i = start; i < end; i++) {
    double dist = vect.distancesq(map[i].GetPoint());
    if(dist < rad_limit2) {
      dist = sqrt(dist);
      Vector3 t_efield = map[i].GetField();
      Vector3 t_point = map[i].GetPoint();
      double weight = squ((rad_limit - dist) / (rad_limit * dist));
      denom_sum += weight;
      num_sum_ex += weight * t_efield.xc();
      num_sum_ey += weight * t_efield.yc();
      num_sum_ez += weight * t_efield.zc();
      //std::cout << "0 " << t_point << ' ' << t_efield << ' ' << dist << std::endl;
      // calculate the x derivative
      double tstep = -2.0*stepsize;
      for(int j = 0; j < 4; j++) {
	Vector3 tempvect = vect+Vector3(tstep,0,0);
	double dist2 = tempvect.distance(t_point);
	weight = squ((rad_limit - dist2) / (rad_limit * dist2));
	xdivnum[j][0] += weight;
	xdivnum[j][1] += weight * t_efield.xc();
	//std::cout << "100 " << tempvect <<  ' ' << tstep << std::endl;
	tstep += stepsize;
	if(j == 1) tstep += stepsize;
      }
      // calculate the y derivative
      tstep = -2.0*stepsize;
      for(int j = 0; j < 4; j++) {
	Vector3 tempvect = vect+Vector3(0,tstep,0);
	double dist2 = tempvect.distance(t_point);
	weight = squ((rad_limit - dist2) / (rad_limit * dist2));
	ydivnum[j][0] += weight;
	ydivnum[j][1] += weight * t_efield.yc();
	//std::cout << "100 " << tempvect << ' ' << tstep << std::endl;
	tstep += stepsize;
	if(j == 1) tstep += stepsize;
      }
      // calculate the z derivative
      tstep = -2.0*stepsize;
      for(int j = 0; j < 4; j++) {
	Vector3 tempvect = vect+Vector3(0,0,tstep);
	double dist2 = tempvect.distance(t_point);
	weight = squ((rad_limit - dist2) / (rad_limit * dist2));
	zdivnum[j][0] += weight;
	zdivnum[j][1] += weight * t_efield.zc();
	//std::cout << "100 " << tempvect <<' ' << tstep <<  std::endl;
	tstep += stepsize;
	if(j == 1) tstep += stepsize;
      }
    }
  }
  Vector3 result(num_sum_ex/ denom_sum, 
	  num_sum_ey/ denom_sum,  
	  num_sum_ez/ denom_sum);

  for(int i = 0; i < 4; i++) {
    xdivnum[i][1] = xdivnum[i][1] / xdivnum[i][0];
    ydivnum[i][1] = ydivnum[i][1] / ydivnum[i][0];
    zdivnum[i][1] = zdivnum[i][1] / zdivnum[i][0];
  }
  xdivnum[0][1] = xdivnum[0][1] - 8 * xdivnum[1][1]
    + 8 * xdivnum[2][1] - xdivnum[3][1];
  ydivnum[0][1] = ydivnum[0][1] - 8 * ydivnum[1][1]
    + 8 * ydivnum[2][1] - ydivnum[3][1];
  zdivnum[0][1] = zdivnum[0][1] - 8 * zdivnum[1][1]
    + 8 * zdivnum[2][1] - zdivnum[3][1];

  gradient.set_coords(xdivnum[0][1] / (12 * stepsize),
		      ydivnum[0][1] / (12 * stepsize),
		      zdivnum[0][1] / (12 * stepsize));
  //std::cout << "100 " << vect << ' ' << result << ' ' << gradient << std::endl;
  return result;
}

Vector3 EFieldMap::Interpolate_Gradient(Vector3 vect,
					Vector3 &gradient)
{
  Vector3 tdx = Interpolate_Field(vect+Vector3(2*stepsize,0,0))*-1 +
    Interpolate_Field(vect+Vector3(stepsize,0,0))*8 +
    Interpolate_Field(vect+Vector3(-stepsize,0,0))*-8 +
    Interpolate_Field(vect+Vector3(-2*stepsize,0,0));
  double dx = tdx.xc() / (12 * stepsize);
  tdx = Interpolate_Field(vect+Vector3(0,2*stepsize,0))*-1 +
    Interpolate_Field(vect+Vector3(0,stepsize,0))*8 +
    Interpolate_Field(vect+Vector3(0,-stepsize,0))*-8 +
    Interpolate_Field(vect+Vector3(0,-2*stepsize,0));
  double dy = tdx.yc() / (12 * stepsize);
  tdx = Interpolate_Field(vect+Vector3(0,0,2*stepsize))*-1 +
    Interpolate_Field(vect+Vector3(0,0,stepsize))*8 +
    Interpolate_Field(vect+Vector3(0,0,-stepsize))*-8 +
    Interpolate_Field(vect+Vector3(0,0,-2*stepsize));
  double dz = tdx.zc() / (12 * stepsize);

  gradient.set_coords(dx,dy,dz);

}
