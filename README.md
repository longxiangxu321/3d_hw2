# Enriching the 3D BAG with building metrics

The code included processes LoD2.2 semantic building models in stored in .json following [cityjson](https://www.cityjson.org/) specification.

It enriches following attributes for the CityObjects of type "Building":

- Volume
- Rectangularity
- Hemisphericality
- Roughness index
- Orientation of the "RoofSurface" surfaces

# Team Member

- *Bingshiuan Tsai*
- *Qiuxian Wei*
- *Longxiang Xu*



# Requirements

The code has been tested on WSL2 with .

Required packages:

- CGAL (with Eigen3 support)
- Eigen 3.3.7



# Folder structure

```
  |-- src/
      |-- main.cpp
      |-- definitions.h
      |-- geomtools.h
      |-- geomtools.cpp
  |-- include/
      |-- json.hpp
  |-- report/
      |-- report.pdf
  |-- CMakeList.txt
```



# Run

To compile and run:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ ./hw02 myfile.city.json





