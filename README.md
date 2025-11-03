![CMake Windows](https://github.com/filthynobleman/disk-segmentation/actions/workflows/cmake-windows.yml/badge.svg)

# Topological Disk Segmentation
This repository implements algorithm presented in the paper *UV Parametrization via Topological Disk Segmentation of Surfaces*.

## Building instructions
The building process is entirely carried out with CMake. If you have not already cloned the repository recursively, or if you have not updated the submodules, please run
```
git submodule update --init --recursive --remote
```

For building the project, you need [CMake](https://cmake.org/) and a C++ compiler compliant with [C++17 standard](https://en.cppreference.com/w/cpp/compiler_support/17).  

From the root directory of the project, execute the following commands:
```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX="../install"
cmake --build . --config release
cmake --install .
```
This will produce the application executables and will install the header and library files.


## Usage
The building process should produce the executable `OptimAtlas`.  

The program computes a geodesic Voronoi partitioning of the surface into the given number of samples. Then it subsamples the regions until each region is a topological disk, and tries to merge as many regions as possible, collapsing pairs of adjacent regions if their union is still a topological disk. Then, it computes the UV parametrization of each region with the given algorithm.  
The program supports the following syntax
```
OptimAtlas path/to/mesh.ext [options]
```
|Option Flag|Option Full Name|Value Types|Meaning|Default|
|-----------|----------------|-----------|-------|-------|
|-n|--num-regions|Integer|The number of Voronoi regions into which the mesh must be subdivided.|10|
|-s|--sub-regions|Integer|The number of subregions into which non-disk regions must be subdivided at each iteration.|5|
|-t|--threshold|Float|The minimum size ratio of the intersection between regions in order to consider them for merging.|0.5|
|-a|--algorithm|harmonic, conformal, arap, tutte|The UV unwrapping algorithm that will be applied to each region.|tutte|
|-o|--output|String|The path to the output file. Must be in OBJ file format.|\*See below|
|-v|--verbosity|Integer|The amount of logging information that will be printed on the standard output.|1|
|-p|--packing|No args|Forces the algorithm to produce islands that do not overlap with each other.|No|

In order to avoid the creation of regions that are too difficult to unwrap, the algorithm tries to avoid merging regions that are touching across a very small boundary. In order for two regions to be mergeable, their number of boundary edges must be at least a fraction of the maximum number of triangles. That is, given two regions `Ri`, `Rj`, their union `Ri + Rj` can be considered only if the intersection `Ri * Rj` satisfies the following relation `boundary(Ri * Rj) >= theta * sqrt(area(Ri + Rj))`. Here, `theta` is the parameter fixed with the option `-t|--threshold`.  

By default, each region is unwrapped to a UV island that is rescaled to fit inside the [0, 1] square. However, for some applications it may be required that islands do not overlap, and thus the option `-p|--packing` can be used place the islands in such a way that they are guaranteed not to overlap.  

Regarding the output filename, the default value is obtained from the input mesh. Unless specified with the option `-o|--output`, the output file will be in the same folder as the input mesh and will have the same name, with the `-uv` suffix and the `.obj` extension. For instance, if the input path is `path/to/mesh.off`, the default output path will be `path/to/mesh-uv.obj`. Even when specified with the `-o|--output` option, the output filename must have the `.obj` extension. If not, the program will ask to automatically change it.  

The verbosity levels are specified in the table below
|Verbosity Level|Information Printed|
|---------------|-------------------|
|0|No output. Errors and warnings are still printed.|
|1|Total runtime for computing the unwrapping. Time for loading and exporting the mesh are not included.|
|2|Runtime for each step of the algorithm, including loading and exporting times.|

The program also supports the help command as
```
OptimAtlas -h|--help