/***********************************************************************
This uses the code package nanoflann in order to search through the
three dimensional grid of electric field points.
 *************************************************************************/
/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright 2011 Jose Luis Blanco (joseluisblancoc@gmail.com).
 *   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/
#include <vector>
#include <fstream>
#include <string>
#include "EFieldMap.h"
#include "EFieldPoint.h"
#include "Vector3.h"
#include "mymath.h"
#include <cmath>
#include <iostream>
#include "nanoflann.hpp"

// EFieldMap constructor

// EFieldMap member function
void EFieldMap::ReadEFieldMap(const size_t N, std::string filename)
{
  std::cerr << "#Generating " << N << " point cloud..." << std::endl;

  std::ifstream fin;
  //  fin.open("efield_sort_closed");
  fin.open(filename.c_str());

  int numLines = 0;
  while ( fin.good() )
    {
      std::string line;
      std::getline(fin, line);
      ++numLines;
    }

  fin.close();

  std::cerr << "#Number of points " << numLines << std::endl;

  fin.open(filename.c_str());


  map.pts.resize(numLines);
  map.efpoints.resize(numLines);
  map.ptdensity.resize(numLines);


  double t_ex, t_ey, t_ez, t_cx, t_cy, t_cz;
  double trash;
  
  int i = 0;
  do {
    //fin >> t_ex >> t_ey >> t_ez >> trash >> trash >> trash
    //	>> t_cx >> t_cy >> t_cz;
    fin >> t_ex >> t_ey >> t_ez
	>> t_cx >> t_cy >> t_cz;
    map.pts[i].x = t_cx;
    map.pts[i].y = t_cy;
    map.pts[i].z = t_cz;
    map.efpoints[i].x = t_ex;
    map.efpoints[i].y = t_ey;
    map.efpoints[i].z = t_ez;
    map.ptdensity[i] = 0;
    i++;
  }while(!fin.eof());

  fin.close();

  mapsize = i;

  std::cerr << "#done" << std::endl;
}

void EFieldMap::PrintMap(int count)
{
  //std::cout << map.size() << std::endl;
  // if (count ==  0) count = map.size();
  // for (int i = 0; i < count; i++) {
  //   Vector3 tPoint = map[i].GetPoint();
  //   std::cout << tPoint << ' ' << std::endl;
  // }
}

void EFieldMap::PrintKeysBins()
{
  // std::cout << keys.size() << std::endl;
  // for (int i = 0; i < keys.size(); i++) {
  //   std::cout << bins[i] << ' ' << keys[i] << std::endl;
  // }
}

void  EFieldMap::GenerateDensity(double rad, const my_kd_tree_t & index)
{
  std::cerr << "#generating density" << std::endl;
  for(int i = 1; i < mapsize; i++) {
    double query_pt[3] = {map.pts[i].x,
			  map.pts[i].y,
			  map.pts[i].z};
    int nc = 0;
    do {
      nanoflann::SearchParams params;
      params.sorted = false;

      const double search_radius = static_cast<double>(rad*rad);
      std::vector<std::pair<size_t,double> > ret_matches;
      const size_t nMatches = index.radiusSearch(&query_pt[0], 
						 search_radius, ret_matches, params);
      nc = nMatches;
      if(nc < 16) rad *= 2;
      if(nc > 128) rad /= 2;
    }while(nc == 0 || nc > 128);
    map.ptdensity[i] = (rad) / double(nc);

    //if(i % 1000 == 0) std::cerr << i << ' ';

  }
  std::cerr << "#done" << std::endl;
}

void EFieldMap::ListPointsinRad(Vector3 vect, double rad, const my_kd_tree_t & index)
{
  const double search_radius = static_cast<double>(rad*rad);
  std::vector<std::pair<size_t,double> >   ret_matches;

  double query_pt[3] = {vect.xc(), vect.yc(), vect.zc()};
  nanoflann::SearchParams params;
  //params.sorted = false;

  const size_t nMatches = index.radiusSearch(&query_pt[0],search_radius, ret_matches, params);

  std::cout << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
  for (size_t i=0;i<nMatches;i++)
    std::cout << "idx["<< i << "]=" << ret_matches[i].first << " dist["
	 << i << "]=" << ret_matches[i].second << std::endl;

}

// Inverse Distance Weighting
// Uses weight function = squ((R-d))/(R*d))
bool EFieldMap::Interpolate_Field(Vector3 vect, const my_kd_tree_t & index,
				  double & rad, Vector3 & field) const
{
  double trad_limit = rad;
  //std::cerr <<  trad_limit << std::endl;
  double query_pt[3] = {vect.xc(), vect.yc(), vect.zc()};
  int count = 0;
  if(true) {
    count++;
    int count2 = 0;
    int nc = 0;
    double denom_sum = 0;
    double num_sum_ex = 0;
    double num_sum_ey = 0;
    double num_sum_ez = 0;
    do {
      count2++;
      denom_sum = 0;
      num_sum_ex = 0;
      num_sum_ey = 0;
      num_sum_ez = 0;
      const double search_radius = static_cast<double>(squ(trad_limit));
      std::vector<std::pair<size_t,double> >   ret_matches;
      nanoflann::SearchParams params;
      params.sorted = false;
      const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
      for (size_t i=0;i<nMatches;i++) {
	double weight = ((trad_limit - sqrt(ret_matches[i].second)) 
			 / (trad_limit * sqrt(ret_matches[i].second)));
	// double weight = (squ(query_pt[0]-map.pts[ret_matches[i].first].x)+
	// 		 squ(query_pt[1]-map.pts[ret_matches[i].first].y)+
	// 		 squ(query_pt[2]-map.pts[ret_matches[i].first].z));
	//double weight = exp(-squ(sqrt(ret_matches[i].second)-trad_limit)/(2*0.00001));
	// double weight = exp(-squ(sqrt(ret_matches[i].second)-trad_limit)/(2*0.0000001))
	//   /(map.ptdensity[ret_matches[i].first]);
	//	weight = weight*weight;
	//weight = weight * map.ptdensity[ret_matches[i].first];
	  denom_sum += weight;
	num_sum_ex += weight * map.efpoints[ret_matches[i].first].x;
	num_sum_ey += weight * map.efpoints[ret_matches[i].first].y;
	num_sum_ez += weight * map.efpoints[ret_matches[i].first].z;
      }
      nc = nMatches;
      if(count > 5) {
	std::cerr << vect << ' ' << nc << ' ' <<  trad_limit << ' ' 
		  << search_radius << std::endl;
      }

      if(trad_limit != trad_limit || trad_limit > 1 || trad_limit < 1e-9) trad_limit = 0.001;

      if(nc < 8) {
	trad_limit *= 3.0;
      }
      if(nc > 32) {
	trad_limit /= 2.0;
      }
      if(count2 > 10) std::cerr << count2 << ' ' << nc << ' ' << trad_limit << std::endl;
      // else if (nc > 128) {
      // 	trad_limit /=2.0;
      // }
    }while((nc < 8 || nc > 32) && (count2 < 10 && nc != 0));
    Vector3 result(num_sum_ex / denom_sum, 
		   num_sum_ey / denom_sum,  
		   num_sum_ez / denom_sum);
    rad = trad_limit;

    if(nc == 0) return false;

    field = result;
    return true;

  }
  else {
    size_t num_results = 16;
    double denom_sum = 0;
    double num_sum_ex = 0;
    double num_sum_ey = 0;
    double num_sum_ez = 0;
    int good_bins;
    //do{
      good_bins = 0;
      denom_sum = 0;
      num_sum_ex = 0;
      num_sum_ey = 0;
      num_sum_ez = 0;
      std::vector<size_t> ret_index(num_results);
      std::vector<double> out_dist_sqr(num_results);
      index.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
      for (size_t i=0;i<num_results;i++) {
	if(CheckBoundary(Vector3(map.pts[ret_index[i]].x,
			      map.pts[ret_index[i]].y,
				 map.pts[ret_index[i]].z))) {
	  // ensures points beyond detector are not used
	  double weight = squ((trad_limit - sqrt(out_dist_sqr[i])) 
			      / (trad_limit * sqrt(out_dist_sqr[i])));
	  denom_sum += weight;
	  num_sum_ex += weight * map.efpoints[ret_index[i]].x;
	  num_sum_ey += weight * map.efpoints[ret_index[i]].y;
	  num_sum_ez += weight * map.efpoints[ret_index[i]].z;
	  good_bins++;
	}
      }
      //num_results *= 2;
      //}while(good_bins < 8);
    Vector3 result(num_sum_ex / denom_sum, 
		   num_sum_ey / denom_sum,  
		   num_sum_ez / denom_sum);
    field = result;

    return true;
  }
}

Vector3 EFieldMap::Interpolate_Field_Gradient(Vector3 vect,
					      Vector3 &gradient)
{
  return vect;
}

Vector3 EFieldMap::Interpolate_Gradient(Vector3 vect,
					Vector3 &gradient)
{
  // Vector3 tdx = Interpolate_Field(vect+Vector3(2*stepsize,0,0))*-1 +
  //   Interpolate_Field(vect+Vector3(stepsize,0,0))*8 +
  //   Interpolate_Field(vect+Vector3(-stepsize,0,0))*-8 +
  //   Interpolate_Field(vect+Vector3(-2*stepsize,0,0));
  // double dx = tdx.xc() / (12 * stepsize);
  // tdx = Interpolate_Field(vect+Vector3(0,2*stepsize,0))*-1 +
  //   Interpolate_Field(vect+Vector3(0,stepsize,0))*8 +
  //   Interpolate_Field(vect+Vector3(0,-stepsize,0))*-8 +
  //   Interpolate_Field(vect+Vector3(0,-2*stepsize,0));
  // double dy = tdx.yc() / (12 * stepsize);
  // tdx = Interpolate_Field(vect+Vector3(0,0,2*stepsize))*-1 +
  //   Interpolate_Field(vect+Vector3(0,0,stepsize))*8 +
  //   Interpolate_Field(vect+Vector3(0,0,-stepsize))*-8 +
  //   Interpolate_Field(vect+Vector3(0,0,-2*stepsize));
  // double dz = tdx.zc() / (12 * stepsize);

  gradient.set_coords(0,0,0);

  return vect;

}

bool EFieldMap::CheckBoundary(Vector3 tvect) const
{
  {
    // const Vector3 detcenter(0.69951, 0.043837, 0.0);
    // const Vector3 normal(0.67165-0.700727, 0.031697-.036563, 0);
    const Vector3 detcenter(0.699, 0.041, 0.0);
    const Vector3 normal(0.9862856015, 0.1650476059, 0);
    
    Vector3 diff = tvect - detcenter;
    if(diff.dot(normal) > 0) {
      return false;
    }
  }
  {
    const Vector3 detcenter(0,0,0);
    const Vector3 normal(1,0,0);
    
    Vector3 diff = tvect - detcenter;
    if(diff.dot(normal) < 0) {
      return false;
    }
  }
  return true;
}
