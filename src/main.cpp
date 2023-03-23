/*
  geo1004.2023
  hw02 help code
  Hugo Ledoux <h.ledoux@tudelft.nl>
  2023-03-01
*/

#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include "definitions.h"
#include "geomtools.h"

//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder
using json = nlohmann::json;

std::vector<Point3> get_coordinates(const json& j, bool translate = true);
void                save2obj(std::string filename, const json& j);
void                enrich_and_save(std::string filename, json& j);



int main(int argc, const char * argv[]) {
  //-- will read the file passed as argument or 2b.city.json if nothing is passed
  const char* filename = (argc > 1) ? argv[1] : "../data/myfile.city.json";
  std::cout << "Processing: " << filename << std::endl;
  std::ifstream input(filename);
  json j;
  input >> j; //-- store the content of the file in a nlohmann::json object
  input.close();
  
  //-- convert each City Object in the file to OBJ and save to a file
  save2obj("out.obj", j);

  //-- enrich with some attributes and save to a new CityJSON 
  enrich_and_save("out.city.json", j);

  return 0;
}


//-- write the OBJ file
void save2obj(std::string filename, const json& j) {
  std::ofstream ofile(filename);
  //-- fetch all the vertices in real-world coordinates (so "transform" is applied)
  std::vector<Point3> lspts = get_coordinates(j, true);
  for (auto& p : lspts) {
    ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
  }
  //-- iterate over each object in the file and output the CDT
  for (auto& co : j["CityObjects"].items()) {
    for (auto& g : co.value()["geometry"]) {
      if ( (g["type"] == "Solid") && (g["lod"] == "2.2") ) {   //-- LoD2.2 only!!!!!
        ofile << "o " << co.key() << std::endl;
        for (int i = 0; i < g["boundaries"].size(); i++) { //-- iterate over each shell
          for (int j = 0; j < g["boundaries"][i].size(); j++) {  // iterate over each face
            std::vector<std::vector<int>> gb = g["boundaries"][i][j];
            std::vector<std::vector<int>> trs = construct_ct_one_face(gb, lspts);
            for (auto& tr : trs) {
              ofile << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
            }
          }
        }
      }
    }
  }
  ofile.close();
  std::cout << "OBJ file written to disk: " << filename << std::endl;
}

//-- add a new attribute "volume" to each City Object and assign a random value
void enrich_and_save(std::string filename, json& j) {
    //-- seed to generate a random number
    //-- https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution

    std::vector<Point3> lspts = get_coordinates(j, true);
    std::map<std::string, int> faces_with_orientation = {{"SW",         1},
                                                         {"WS",         2},
                                                         {"WN",         3},
                                                         {"NW",         4},
                                                         {"NE",         5},
                                                         {"EN",         6},
                                                         {"ES",         7},
                                                         {"SE",         8},
                                                         {"horizontal", 9}};
    for (auto &co: j["CityObjects"].items()) {    //.items() return the key-value pairs belong to the CityObjects
        if (co.value()["type"] == "BuildingPart") {
            std::vector<std::vector<std::vector<int>>> trss;
            std::vector<Point3> exterior_pts;
            for (auto &g: co.value()["geometry"]) {


                auto origin_surface_num = g["semantics"]["surfaces"].size();
                for (const auto &ori: faces_with_orientation) {
                    json new_surface = {
                            {"type",        "RoofSurface"},
                            {"Orientation", ori.first}
                    };
                    g["semantics"]["surfaces"].push_back(new_surface);
                }

                for (int i = 0; i < g["boundaries"].size(); i++) { //-- iterate over each shell
                    for (int j = 0; j < g["boundaries"][i].size(); j++) {  // iterate over each face
                        std::vector<std::vector<int>> gb = g["boundaries"][i][j];
                        std::vector<std::vector<int>> trs = construct_ct_one_face(gb, lspts);
                        trss.push_back(trs);
                        if (g["semantics"]["values"][0][j] == 1) { // extract roof surface ct
                            if (trs.size() < 1) {
                                continue;
                            } else {
                                std::string orient = roof_orientation(trs[0], lspts);
                                g["semantics"]["values"][0][j] =
                                        faces_with_orientation[orient] + origin_surface_num - 1;
                            }
                        }
                    }
                }

                for (auto &ex_faces: g["boundaries"][0]) { // extract the exterior shell only
                    for (auto &ex_face: ex_faces) {  // iterate over each face of exterior shell
                        for (int i = 0; i < ex_face.size(); i++) { // iterate over the pt_index of each face
                            Point3 exterior_pt = lspts[ex_face[i]];
                            exterior_pts.push_back(exterior_pt);
                        }
                    }
                }



                std::pair<std::vector<double>, std::vector<double>> pair = calculate_volume_area(trss, lspts);
                std::vector<double> vol = pair.first;
                std::vector<double> area_list = pair.second;
                double rec = calculate_rectangularity(exterior_pts, vol[0]);
                double hem = hemisphericality(vol[0], vol[1]);
                double ri = Roughness_index(vol[0], vol[1], area_list, trss, lspts);
                co.value()["attributes"]["volume"] = vol[0];
                co.value()["attributes"]["area"] = vol[1];
                co.value()["attributes"]["rectangularity"] = rec;
                co.value()["attributes"]["hemisphericality"] = hem;
                co.value()["attributes"]["roughness"] = ri;
            }
        }
    }

    for (auto &co: j["CityObjects"].items()) {    //.items() return the key-value pairs belong to the CityObjects
        if (co.value()["type"] == "Building") {
            std::cout <<"1" <<std::endl;
            if (co.value()["children"].size()<1) {
                continue;
            }
            else {
                std::string child_name = co.value()["children"][0];
                auto child = j["CityObjects"].find(child_name);
                co.value()["attributes"]["volume"] = child.value()["attributes"]["volume"];
                co.value()["attributes"]["rectangularity"] = child.value()["attributes"]["rectangularity"];
                co.value()["attributes"]["hemisphericality"] = child.value()["attributes"]["hemisphericality"];
                co.value()["attributes"]["roughness"] = child.value()["attributes"]["roughness"];
                child.value()["attributes"].erase("volume");
                child.value()["attributes"].erase("rectangularity");
                child.value()["attributes"].erase("hemisphericality");
                child.value()["attributes"].erase("roughness");
            }
        }
    }

        //-- write to disk the modified city model (myfile.city.json)
        std::ofstream o(filename);
        o << j.dump(2) << std::endl;
        o.close();
        std::cout << "Enriched CityJSON file written to disk: " << filename << std::endl;

}


//-- get real-world coordinates from the "vertices" property
//-- https://www.cityjson.org/specs/#transform-object
//-- param translate is to use the translation in the "transform",
//-- it can be put to false to make the coords smaller (and better for computations)
std::vector<Point3> get_coordinates(const json& j, bool translate) {
  std::vector<Point3> lspts;
  std::vector<std::vector<int>> lvertices = j["vertices"];
  if (translate) {
    for (auto& vi : lvertices) {
      double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
      double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
      double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
      lspts.push_back(Point3(x, y, z));
    } 
  } else {
    //-- do not translate, useful to keep the values low for downstream processing of data
    for (auto& vi : lvertices) {
      double x = (vi[0] * j["transform"]["scale"][0].get<double>());
      double y = (vi[1] * j["transform"]["scale"][1].get<double>());
      double z = (vi[2] * j["transform"]["scale"][2].get<double>());
      lspts.push_back(Point3(x, y, z));
    }
  }
  return lspts;
}
