#include "json.hpp"
using json = nlohmann::json;
#include <cmath>
#include "geomtools.h"


std::vector<std::vector<int>>
construct_ct_one_face(const std::vector<std::vector<int>>& lsRings, 
                      const std::vector<Point3>& lspts)
{
  std::vector<std::vector<int>> re;

  //-- find best fitted plane (only based on oring)
  std::vector<Point3> planepts;
  for (auto& each : lsRings[0]) {
    planepts.push_back(lspts[each]);
  }
  Plane bestfitplane = get_best_fitted_plane(planepts);
  //-- check orientation (for good normals for the output, pointing outwards)
  bool reversed = false;
  Polygon2 pgn;
  for (auto& each : lsRings[0]) {
    Point3 p = lspts[each];
    pgn.push_back(bestfitplane.to_2d(p));
  }

  //-- stop if the oring is not simple (could crash)
  if (pgn.is_simple() == false) {
    return re;
  }
  //-- check orientation, this works all must be ccw
  if (pgn.is_counterclockwise_oriented() == false) {
    reversed = true;
  }
    
  CT ct;
  for (auto& ring : lsRings) {
    //-- make another *closed* ring for simplicity
    std::vector<int> r2 = ring;
    r2.push_back(ring.front());
    std::vector<int>::const_iterator it;
    std::vector<int>::iterator it2 = std::prev(r2.end());  // r2.end() points to the one after the last element, not
                                                            // and it does not exist
    for (it = r2.begin(); it != it2; it++) {
      Point2 p0 = bestfitplane.to_2d(lspts[*it]);
      CT::Vertex_handle v0 = ct.insert(p0);
      v0->id() = *it;
      auto it3 = it;
      it3++;
      Point2 p1 = bestfitplane.to_2d(lspts[*it3]);
      CT::Vertex_handle v1 = ct.insert(p1);
      v1->id() = *it3;
      if (v0 != v1) {
        ct.insert_constraint(v0, v1);
      }
    }
  }
  mark_domains(ct); 
  if (!ct.is_valid()) 
    return re;
  for (CT::Finite_faces_iterator fit = ct.finite_faces_begin();
       fit != ct.finite_faces_end(); 
       ++fit) 
  {
    if (fit->info().in_domain()) {
      std::vector<int> t;
      t.push_back(fit->vertex(0)->id());
      if (reversed) {
        t.push_back(fit->vertex(2)->id());
        t.push_back(fit->vertex(1)->id());
      } else {
        t.push_back(fit->vertex(1)->id());
        t.push_back(fit->vertex(2)->id());
      }
      re.push_back(t);
    }
  }
  return re;
}


Plane get_best_fitted_plane(const std::vector<Point3> &lspts)
{
  Plane p;
  linear_least_squares_fitting_3(lspts.begin(), lspts.end(), p, CGAL::Dimension_tag<0>());  
  K::FT tol = 1e-12;
  K::FT a = p.a();
  K::FT b = p.b();
  K::FT c = p.c();
  bool updated = false;
  if (CGAL::abs(p.a()) < tol) {
    a = 0;
    updated = true;  
  }
  if (CGAL::abs(p.b()) < tol) {
    b = 0;
    updated = true;  
  }
  if (CGAL::abs(p.c()) < tol) {
    c = 0;
    updated = true;  
  }
  if (updated == false)
    return p;
  else {
    return Plane(a, b, c, p.d());
  }
}

void mark_domains(CT& ct,
  CT::Face_handle start,
  int index,
  std::list<CT::Edge>& border) 
{
  if (start->info().nesting_level != -1) {
    return;
  }
  std::list<CT::Face_handle> queue;
  queue.push_back(start);
  while (!queue.empty()) {
    CT::Face_handle fh = queue.front();
    queue.pop_front();
    if (fh->info().nesting_level == -1) {
      fh->info().nesting_level = index;
      for (int i = 0; i < 3; i++) {
        CT::Edge e(fh, i);
        CT::Face_handle n = fh->neighbor(i);
        if (n->info().nesting_level == -1) {
          if (ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}

//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident 
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void mark_domains(CT& ct) {
  for (CT::All_faces_iterator it = ct.all_faces_begin(); it != ct.all_faces_end(); ++it) {
    it->info().nesting_level = -1;
  }
  std::list<CT::Edge> border;
  mark_domains(ct, ct.infinite_face(), 0, border);
  while (!border.empty()) {
    CT::Edge e = border.front();
    border.pop_front();
    CT::Face_handle n = e.first->neighbor(e.second);
    if (n->info().nesting_level == -1) {
      mark_domains(ct, n, e.first->info().nesting_level + 1, border);
    }
  }
}

int cal_orientation(const std::vector<int> tri, const std::vector<Point3>& lspts, const Point3 &poi) {
    Point3 p1 = lspts[tri[0]];
    Point3 p2 = lspts[tri[1]];
    Point3 p3 = lspts[tri[2]];
    auto orient = CGAL::orientation(p1, p2, p3, poi);
    if(orient==CGAL::POSITIVE) {return -1;}
    else if (orient==CGAL::NEGATIVE) {return 1;}
    else {return 0;}
}

double cal_area(const std::vector<int> tri, const std::vector<Point3>& lspts) {
    Point3 p1 = lspts[tri[0]];
    Point3 p2 = lspts[tri[1]];
    Point3 p3 = lspts[tri[2]];

    Vector3 v1 = p2 - p1;
    Vector3 v2 = p3 - p1;

    // Cross product of two vectors is perpendicular to both
    Vector3 cross = CGAL::cross_product(v1, v2);

    // Magnitude of cross product is area of parallelogram defined by v1 and v2
    double area = sqrt(cross.squared_length())/2;

    return area;
}

double cal_distance(const std::vector<int> tri, const std::vector<Point3>& lspts, const Point3 &poi) {
    Point3 p1 = lspts[tri[0]];
    Point3 p2 = lspts[tri[1]];
    Point3 p3 = lspts[tri[2]];

    // Calculate the normal vector of the plane defined by the triangle
    Plane plane(p1, p2, p3);

    Point3 projected = plane.projection(poi);

    double distance = squared_distance(projected, poi);
    double dist = sqrt(distance);
    // Calculate the distance between poi and the closest point on the plane
    return dist;
}

std::vector<double> calculate_volume_area(const std::vector<std::vector<std::vector<int>>> trss,
                                          const std::vector<Point3>& lspts) {
    int bp_index = trss[0][0][0];
    Point3 bp = lspts[bp_index];
    double volume = 0;
    double area = 0;
    for (auto trs:trss) {
        for (auto& tri:trs) {
            double tri_area = cal_area(tri, lspts);
            double distance = cal_distance(tri, lspts, bp);
            int orient = cal_orientation(tri,lspts,bp);
            volume += orient * tri_area * distance/3;
            area += tri_area;
        }
    }
    std::vector<double> result;
    result.push_back(volume);
    result.push_back(area);
    return result;
}