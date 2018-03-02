
#include "EFieldPoint.h"
#include <vector>
#include <string>
#include "nanoflann.hpp"

using namespace nanoflann;

#ifndef E_FIELD_MAP_H
#define E_FIELD_MAP_H

template <typename T>
struct PointCloud
{
  struct Point
  {
    T  x,y,z;
  };

  std::vector<Point>  pts;
  std::vector<Point> efpoints;
  std::vector<double> ptdensity;

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return pts.size(); }

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t size) const
  {
    const T d0=p1[0]-pts[idx_p2].x;
    const T d1=p1[1]-pts[idx_p2].y;
    const T d2=p1[2]-pts[idx_p2].z;
    return d0*d0+d1*d1+d2*d2;
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline T kdtree_get_pt(const size_t idx, int dim) const
  {
    if (dim==0) return pts[idx].x;
    else if (dim==1) return pts[idx].y;
    else return pts[idx].z;
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX &bb) const { return false; }

};


// This is an exampleof a custom data set class
typedef KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<double, PointCloud<double> > ,
		  PointCloud<double>,
		  3 /* dim */
		  > my_kd_tree_t;
/* typedef KDTreeSingleIndexAdaptor< */
/* 		L1_Adaptor<double, PointCloud<double> > , */
/* 		  PointCloud<double>, */
/* 		  3 /\* dim *\/ */
/* 		  > my_kd_tree_t; */
 
class EFieldMap
{
 private:
  //  std::vector < EFieldPoint > map;
  std::vector < int > bins;
  std::vector < double > keys;
  double bin_spacing;
  double stepsize;
  double rad_limit;
  int nmatch;
  int mapsize;

	
 public:
  //EFieldMap();
  PointCloud<double> map;

  EFieldMap(double tstep) {stepsize = tstep;} // private default constructor
 
  void ReadEFieldMap(const size_t N, std::string filename);

  void PrintMap(int count);

  //  void SetRadius(double radius);

  void PrintKeysBins();

  void ListPointsinRad(Vector3 vect, double rad, const my_kd_tree_t & index);
  void GenerateDensity(double rad, const my_kd_tree_t & index);
  
  bool Interpolate_Field(Vector3 vect, const my_kd_tree_t & index,
			 double & rad, Vector3 & field) const;
  Vector3 Interpolate_Field_Gradient(Vector3 vect, Vector3 &gradient);
  Vector3 Interpolate_Gradient(Vector3 vect, Vector3 &gradient);
  bool CheckBoundary(Vector3 tvect) const;


};

#endif
