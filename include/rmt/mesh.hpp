/**
 * @file        mesh.hpp
 * 
 * @brief       A data structure representing a triangular mesh.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#pragma once

#include <Eigen/Dense>
#include <string>

#include <rmt/graph.hpp>

namespace rmt
{

enum SubMeshAccessType
{
    BY_VERTICES = 1,
    BY_VERTICES_ALL = BY_VERTICES | (BY_VERTICES << 1),
    BY_VERTICES_ANY = BY_VERTICES,
    BY_TRIANGLES = BY_VERTICES << 2
};

enum UVMappingAlgorithm
{
    HARMONIC,
    CONFORMAL,
    ARAP,
    TUTTE,
    NONE
};
    
class Mesh
{
private:
    // Mesh vertices #V x 3
    Eigen::MatrixXd m_V;

    // Mesh triangles #F x 3
    Eigen::MatrixXi m_F;

    // Mesh (unique) edges #E x 2
    Eigen::MatrixXi m_E;

    // Incidence map edges on triangles
    Eigen::MatrixXi m_T2E;


public:
    Mesh(const Eigen::MatrixXd& V,
         const Eigen::MatrixXi& F);
    Mesh(const std::string& Filename);
    Mesh(const rmt::Mesh& M);
    Mesh(rmt::Mesh&& M);
    rmt::Mesh& operator=(const rmt::Mesh& M);
    rmt::Mesh& operator=(rmt::Mesh&& M);
    ~Mesh();

    int NumVertices() const;
    int NumEdges() const;
    int NumTriangles() const;

    const Eigen::MatrixXd& GetVertices() const;
    const Eigen::MatrixXi& GetTriangles() const;
    const Eigen::MatrixXi& GetEdges() const;
    const Eigen::MatrixXi& GetT2E() const;

    void ComputeEdgesAndIncidence();

    void Scale(double Alpha);
    void Translate(const Eigen::Vector3d& Movement);
    void Translate(const Eigen::VectorXd& Movement);

    void CenterAtOrigin();
    void RescaleInsideUnitBox();
    void RescaleInsideUnitSphere();

    void Resample(int OutputSize);
    void MakeManifold();
    void RemoveSmallComponents(double AreaFraction = 1e-2);
    void RemoveDegenaracies(double DistanceThreshold = 1e-4);
    void CleanUp(double AreaFraction = 1e-2, double DistanceThreshold = 1e-4);

    rmt::Mesh SubMesh(const std::vector<int>& Indices,
                      rmt::SubMeshAccessType AccType,
                      std::vector<int>& VIdx,
                      std::vector<int>& TIdx) const;

    Eigen::MatrixXd UVMapping(rmt::UVMappingAlgorithm Algo,
                              bool OptimalRotation = false) const;


    rmt::Graph AsGraph() const;
    rmt::Graph AsDualGraph(const std::vector<std::pair<int, int>>& Ignore = { })const;


    void Export(const std::string& Filename) const;
};


} // namespace rmt
