#ifndef DEFS_H
#define DEFS_H

// CGAL kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/squared_distance_3.h>

//-- for mark_domain()
struct FaceInfo2
{
  FaceInfo2() {}
  int nesting_level;
  bool in_domain() {
    return nesting_level % 2 == 1;
  }
};  

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2                  Point2;
typedef K::Point_3                  Point3;
typedef CGAL::Polygon_2<K>          Polygon2;
typedef K::Plane_3                  Plane;
typedef K::Vector_3 Vector3;

typedef CGAL::Triangulation_vertex_base_with_id_2 <K>             Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>   Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>       Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_intersections_tag                             Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CT;


#endif
