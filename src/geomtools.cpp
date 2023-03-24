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

std::pair<std::vector<double>, std::vector<double>> calculate_volume_area(const std::vector<std::vector<std::vector<int>>> trss,
                                                                          const std::vector<Point3>& lspts) {
    int bp_index = trss[0][0][0];
    Point3 bp = lspts[bp_index];
    double volume = 0;
    double area = 0;
    std::vector<double> area_list;
    for (auto trs:trss) {
        for (auto& tri:trs) {
            double tri_area = cal_area(tri, lspts);
            double distance = cal_distance(tri, lspts, bp);
            int orient = cal_orientation(tri,lspts,bp);
            volume += orient * tri_area * distance/3;
            area += tri_area;
            //save area of every triangle for later
            area_list.push_back(tri_area);
        }
    }
    std::vector<double> result;
    result.push_back(volume);
    result.push_back(area);
    return make_pair(result,area_list);
}

double hemisphericality(double volume,double area){
    double hem = (3 * sqrt(2*M_PI) * volume) / pow(area,1.5);
    return hem;
}

//caculate distance of a point to center point
double pt_to_center(Point3 pt,double x,double y,double z){
    double distance = sqrt(pow((pt.x() - x),2)+ pow((pt.y() - y),2)+pow((pt.z() - z),2));
    return distance;
}
//caculate the distance to center of the building
double distance_to_center(double area,const std::vector<double> area_list,const std::vector<std::vector<std::vector<int>>> trss,
                          const std::vector<Point3>& lspts){
    //caculate the center of the building
    //The weight of each triangle's center based on the its area / all trianles' area
    std::vector<double> weights;
    for (auto each_area:area_list) {
        double weight = each_area / area;
        weights.push_back(weight);
    }
    std::vector<std::vector<double>> tri_centers;//center of each triangle
    std::vector<std::vector<Point3>> randompts;
    for (auto trs:trss) {
        for (auto& tri:trs) {
            std::vector<double> tri_center;
            Point3 p1 = lspts[tri[0]];
            Point3 p2 = lspts[tri[1]];
            Point3 p3 = lspts[tri[2]];
            tri_center.push_back((p1.x() + p2.x() + p3.x()) / 3.0);
            tri_center.push_back((p1.y() + p2.y() + p3.y()) / 3.0);
            tri_center.push_back((p1.z() + p2.z() + p3.z()) / 3.0);
            tri_centers.push_back(tri_center);
            // Define the Random_points_in_triangle_3 generator
            CGAL::Random_points_in_triangle_3<Point3> generator(p1, p2, p3);
            // Generate 2 random points within the triangle
            std::vector<Point3> randompt;
            CGAL::copy_n(generator, 2, std::back_inserter(randompt));
            randompts.push_back(randompt);
        }
    }
    //center of the building
    double x,y,z;
    x = y = z =0;
    //caculate the center of the building with weights
    for (int i =0; i < weights.size(); i++){
        x+=tri_centers[i][0]*weights[i];
        y+=tri_centers[i][1]*weights[i];
        z+=tri_centers[i][2]*weights[i];
    }
    double distance_sum;//the average distance of points sampled to the center
    for (int i =0; i < weights.size(); i++){
        Point3 pt1 = randompts[i][0];
        Point3 pt2 = randompts[i][1];
        //caculate the avg distance of two random points to center
//        std::cout << x<<" " << y <<" "<<z<<std::endl;
//        std::cout << pt1.x()<<" " << pt1.y() <<" "<<pt1.z()<<std::endl;
        double dist = (pt_to_center(pt1,x,y,z) +pt_to_center(pt2,x,y,z))/2;
//        std::cout << dist<<std::endl;
        distance_sum+= dist * weights[i];
    }
    return distance_sum;
}

double Roughness_index(double volume,double area,const std::vector<double> area_list,const std::vector<std::vector<std::vector<int>>> trss,
                       const std::vector<Point3>& lspts){
    double down = pow(area, 1.5) + volume;
    double distance_sum = distance_to_center(area,area_list,trss,lspts);
    double ri = pow(distance_sum,3) * 48.735 /down;
    return ri;
}

double calculate_rectangularity(const std::vector<Point3> exterior_pts, double vol){
    // Compute the optimal bounding box, and return the 8 coordinates of the cube
    std::array<Point3, 8> obb_points;
    CGAL::oriented_bounding_box(exterior_pts, obb_points, CGAL::parameters::use_convex_hull(true));

    // Calculate the volume of the oobb
    double length = sqrt(squared_distance(obb_points[1], obb_points[0]));
    double width = sqrt(squared_distance(obb_points[2], obb_points[1]));
    double height = sqrt(squared_distance(obb_points[7], obb_points[2]));
    double volume_oobb = length * width * height;

    // compute the rectangularity
    double rectangularity = vol / volume_oobb;

    // varify whether the rectangularity is in the interval of [0, 1]
    if (rectangularity > 1){
        std::cout << "Rectangularity is greater than 1" << std::endl;
    }
    else if(rectangularity < 0){
        std::cout << "Rectangularity is negative" << std::endl;
    }
    else{
        return rectangularity;
    }
}

std::pair<int, std::string> roof_orientation(const std::vector<int> tri, const std::vector<Point3>& lspts){
    std::pair <int, std::string> orientation;
    // compute the normal vector of the selected triangle
    Point3 p1 = lspts[tri[0]];
    Point3 p2 = lspts[tri[1]];
    Point3 p3 = lspts[tri[2]];

    Vector3 v1 = p2 - p1;
    Vector3 v2 = p3 - p1;

    Vector3 cross =  CGAL::cross_product(v1, v2);

    K::FT x = cross.x();
    K::FT y = cross.y();

    double tan = y / x;

    if (abs(x) <= 0.0001 && abs(y) <= 0.0001) {
        orientation= std::make_pair(8, "horizontal");
    }
    else if (abs(x) <= 0.0001 && y > 0){ // North
        if (tan > 1) {
            orientation = std::make_pair(0, "NE");
        }
        else{
            orientation = std::make_pair(7, "NW");
        }
    }
    else if (abs(x) <= 0.0001 && y < 0){ // South
        if (tan > 1) {
            orientation = std::make_pair(1, "SW");
        }
        else{
            orientation = std::make_pair(6, "SE");
        }
    }
    else if ( x > 0 && y >= 0){
        if (tan >= 1 ){
            orientation = std::make_pair(0, "NE");
        }
        else{
            orientation = std::make_pair(2, "EN");
        }
    }
    else if ( x > 0 && y < 0){
        if (tan >= -1 ){
            orientation = std::make_pair(4, "ES");
        }
        else{
            orientation = std::make_pair(6, "SE");
        }
    }
    else if ( x < 0 && y <= 0){
        if (tan >= 1 ){
            orientation = std::make_pair(1, "SW");
        }
        else{
            orientation = std::make_pair(3, "WS");
        }
    }
    else if ( x < 0 && y > 0){
        if (tan >= -1 ){
            orientation = std::make_pair(5, "WN");
        }
        else{
            orientation = std::make_pair(7, "NW");
        }
    }
    return orientation;
}

std::pair<int, std::string> roof_orientation(Plane best_face){
    std::pair <int, std::string> orientation;

    K::FT x = best_face.a();
    K::FT y = best_face.b();

    double tan = y / x;

    if (abs(x) <= 0.0001 && abs(y) <= 0.0001) {
        orientation= std::make_pair(8, "horizontal");
    }
    else if (abs(x) <= 0.0001 && y > 0){ // North
        if (tan > 1) {
            orientation = std::make_pair(0, "NE");
        }
        else{
            orientation = std::make_pair(7, "NW");
        }
    }
    else if (abs(x) <= 0.0001 && y < 0){ // South
        if (tan > 1) {
            orientation = std::make_pair(1, "SW");
        }
        else{
            orientation = std::make_pair(6, "SE");
        }
    }
    else if ( x > 0 && y >= 0){
        if (tan >= 1 ){
            orientation = std::make_pair(0, "NE");
        }
        else{
            orientation = std::make_pair(2, "EN");
        }
    }
    else if ( x > 0 && y < 0){
        if (tan >= -1 ){
            orientation = std::make_pair(4, "ES");
        }
        else{
            orientation = std::make_pair(6, "SE");
        }
    }
    else if ( x < 0 && y <= 0){
        if (tan >= 1 ){
            orientation = std::make_pair(1, "SW");
        }
        else{
            orientation = std::make_pair(3, "WS");
        }
    }
    else if ( x < 0 && y > 0){
        if (tan >= -1 ){
            orientation = std::make_pair(5, "WN");
        }
        else{
            orientation = std::make_pair(7, "NW");
        }
    }
    return orientation;
}