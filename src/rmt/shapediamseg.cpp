/**
 * @file        shapediamseg.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2025-08-29
 */
#include <rmt/segmentation.hpp>
#include <rmt/utils.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
 
#include <CGAL/mesh_segmentation.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/property_map.h>
 
#include <iostream>
#include <fstream>
 
typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;
 
typedef CGAL::Surface_mesh<Point_3>                           CGALMesh;
 
typedef boost::graph_traits<CGALMesh>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<CGALMesh>::face_descriptor        face_descriptor;


Eigen::VectorXi rmt::ShapeDiameterSegmentation(const rmt::Mesh& M)
{
    CGALMesh mesh;
    for (int i = 0; i < M.NumVertices(); ++i)
        mesh.add_vertex({ M.GetVertices()(i, 0), M.GetVertices()(i, 1), M.GetVertices()(i, 2) });
    for (int i = 0; i < M.NumTriangles(); ++i)
    {
        CGAL::SM_Vertex_index v0(M.GetTriangles()(i, 0));
        CGAL::SM_Vertex_index v1(M.GetTriangles()(i, 1));
        CGAL::SM_Vertex_index v2(M.GetTriangles()(i, 2));
        mesh.add_face(v0, v1, v2);
    }
    
    typedef CGALMesh::Property_map<face_descriptor,double> Facet_double_map;
    Facet_double_map sdf_property_map;
    
    sdf_property_map = mesh.add_property_map<face_descriptor,double>("f:sdf").first;
    
    // compute SDF values
    // We can't use default parameters for number of rays, and cone angle
    // and the postprocessing
    CGAL::sdf_values(mesh, sdf_property_map, 2.0 / 3.0 * CGAL_PI, 25, true);
    
    // create a property-map for segment-ids
    typedef CGALMesh::Property_map<face_descriptor, std::size_t> Facet_int_map;
    Facet_int_map segment_property_map = mesh.add_property_map<face_descriptor,std::size_t>("f:sid").first;;
    
    // segment the mesh using default parameters for number of levels, and smoothing lambda
    // Any other scalar values can be used instead of using SDF values computed using the CGAL function
    std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);
    
    Eigen::VectorXi Partition;
    Partition.resize(M.NumTriangles());
    
    // print segment-ids
    int i = 0;
    for(face_descriptor f : faces(mesh) ) {
        Partition[i++] =  segment_property_map[f];
    }

    
    return Partition;
}