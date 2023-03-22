#ifndef __geomtools__
#define __geomtools__

#include "definitions.h"


Plane                 get_best_fitted_plane(const std::vector<Point3> &lspts);
void                  mark_domains(CT& ct);
void                  mark_domains(CT& ct, CT::Face_handle start, int index, std::list<CT::Edge>& border);

std::vector<std::vector<int>>
construct_ct_one_face(const std::vector<std::vector<int>>& lsRings, 
                      const std::vector<Point3>& lspts);
std::pair<std::vector<double>, std::vector<double>> calculate_volume_area(const std::vector<std::vector<std::vector<int>>> trss,
                                          const std::vector<Point3>& lspts);
double hemisphericality(double volume,double area);
double Roughness_index(double volume,double area,const std::vector<double> area_list,const std::vector<std::vector<std::vector<int>>> trss,
                       const std::vector<Point3>& lspts);
#endif 