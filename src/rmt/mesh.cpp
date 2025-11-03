/**
 * @file        mesh.cpp
 * 
 * @brief       Implementation of class rmt::Mesh.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#include <rmt/mesh.hpp>
#include <rmt/io.hpp>
#include <rmt/preprocess.hpp>
#include <rmt/clean.hpp>
#include <rmt/utils.hpp>
#include <rmt/embed.hpp>
#include <cut/cut.hpp>
#include <cassert>
#define NOMINMAX
#include <igl/unique_edge_map.h>
#include <igl/barycenter.h>
#include <igl/doublearea.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <unordered_set>

rmt::Mesh::Mesh(const Eigen::MatrixXd& V,
                const Eigen::MatrixXi& F)
    : m_V(V), m_F(F)
{
    m_E.resize(0, 2);
    m_T2E.resize(0, 3);
}

rmt::Mesh::Mesh(const std::string& Filename)
{
    CUTAssert(rmt::LoadMesh(Filename, m_V, m_F));
    m_E.resize(0, 2);
    m_T2E.resize(0, 3);
}

rmt::Mesh::Mesh(const rmt::Mesh& M)
    : m_V(M.m_V), m_F(M.m_F), m_E(M.m_E), m_T2E(M.m_T2E)
{ }

rmt::Mesh& rmt::Mesh::operator=(const rmt::Mesh& M)
{
    m_V = M.m_V;
    m_F = M.m_F;
    m_E = M.m_E;
    m_T2E = M.m_T2E;
    return *this;
}

rmt::Mesh::Mesh(rmt::Mesh&& M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    m_E = std::move(M.m_E);
    m_T2E = std::move(M.m_T2E);
}

rmt::Mesh& rmt::Mesh::operator=(rmt::Mesh&& M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    m_E = std::move(M.m_E);
    m_T2E = std::move(M.m_T2E);
    return *this;
}

rmt::Mesh::~Mesh() { }


int rmt::Mesh::NumVertices() const          { return m_V.rows(); }
int rmt::Mesh::NumEdges() const             { return m_E.rows(); }
int rmt::Mesh::NumTriangles() const         { return m_F.rows(); }


const Eigen::MatrixXd& rmt::Mesh::GetVertices() const       { return m_V; }
const Eigen::MatrixXi& rmt::Mesh::GetTriangles() const      { return m_F; }
const Eigen::MatrixXi& rmt::Mesh::GetEdges() const
{ 
    CUTAssert(m_E.size() > 0);
    return m_E;
}
const Eigen::MatrixXi& rmt::Mesh::GetT2E() const
{
    CUTAssert(m_T2E.size() > 0);
    return m_T2E;
}


void rmt::Mesh::ComputeEdgesAndIncidence()
{
    // Compute all the edges and their incidence on triangles
    std::unordered_map<std::pair<int, int>, int, rmt::PairHash<int>> m_Ecount;
    m_Ecount.reserve(2 * NumTriangles());
    m_T2E.resize(NumTriangles(), 3);
    for (int i = 0; i < m_F.rows(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int j1 = j + 1;
            if (j1 > 2)
                j1 = 0;
            std::pair<int, int> e{m_F(i, j), m_F(i, j1)};
            if (e.first > e.second)
                std::swap(e.first, e.second);
            if (m_Ecount.find(e) == m_Ecount.end())
                m_Ecount.emplace(e, (int)m_Ecount.size());
            m_T2E(i, j) = m_Ecount[e];
        }
    }

    // Edge matrix
    m_E.resize(m_Ecount.size(), 2);
    for (auto e : m_Ecount)
        m_E.row(e.second) = Eigen::RowVector2i{ e.first.first, e.first.second };
}


void rmt::Mesh::Scale(double Alpha)
{
    m_V *= Alpha;
}


void rmt::Mesh::Translate(const Eigen::VectorXd& Movement)
{
    Eigen::Vector3d Mov3 = Movement.segment<3>(0);
    Translate(Mov3);
}

void rmt::Mesh::Translate(const Eigen::Vector3d& Movement)
{
    m_V.rowwise() += Movement.transpose();
}



void rmt::Mesh::CenterAtOrigin()
{
    Eigen::MatrixXd BC;
    igl::barycenter(m_V, m_F, BC);
    Eigen::VectorXd Center = BC.colwise().mean();
    Translate(Center);
}


void rmt::Mesh::RescaleInsideUnitBox()
{
    double L = m_V.cwiseAbs().maxCoeff();
    m_V /= L;
}

void rmt::Mesh::RescaleInsideUnitSphere()
{
    double L = std::sqrt(m_V.rowwise().squaredNorm().maxCoeff());
    m_V /= L;
}


void rmt::Mesh::Resample(int OutputSize)
{
    double MEL = rmt::MaxEdgeLength(m_V, m_F, OutputSize);
    rmt::ResampleMesh(m_V, m_F, MEL);
}


void rmt::Mesh::MakeManifold()
{
    rmt::MakeManifold(m_V, m_F);
}

void rmt::Mesh::RemoveSmallComponents(double AreaFraction)
{
    rmt::RemoveSmallComponents(m_V, m_F, AreaFraction);
}

void rmt::Mesh::RemoveDegenaracies(double DistanceThreshold)
{
    rmt::RemoveDegeneracies(m_V, m_F, DistanceThreshold);
}

void rmt::Mesh::CleanUp(double AreaFraction, double DistanceThreshold)
{
    rmt::CleanUp(m_V, m_F, AreaFraction, DistanceThreshold);
}





rmt::Graph rmt::Mesh::AsGraph() const
{
    return rmt::Graph(m_V, m_F);
}

rmt::Graph rmt::Mesh::AsDualGraph(const std::vector<std::pair<int, int>>& Ignore) const
{
    // Ignore map
    std::unordered_set<std::pair<int, int>, rmt::PairHash<int>> IgnoreSet;
    for (auto e : Ignore)
    {
        auto ee = e;
        if (ee.first > ee.second)
            std::swap(ee.first, ee.second);
        IgnoreSet.emplace(ee);
    }

    // Get edge 2 triangle map
    std::unordered_map<std::pair<int, int>, std::pair<int, int>, rmt::PairHash<int>> E2T;
    E2T.reserve(2 * NumTriangles());
    for (int i = 0; i < NumTriangles(); ++i)
    {
        for (int j1 = 0; j1 < 3; ++j1)
        {
            int j2 = j1 + 1;
            if (j2 == 3)
                j2 = 0;
            std::pair<int, int> e{ m_F(i, j1), m_F(i, j2) };
            if (e.first > e.second)
                std::swap(e.first, e.second);
            if (IgnoreSet.find(e) != IgnoreSet.end())
                continue;
            if (E2T.find(e) == E2T.end())
            {
                std::pair<int, int> newt{ i, -1 };
                E2T.emplace(e, newt);
            }
            else
                E2T[e].second = i;
        }
    }

    // Remove boundary edges, they are unnecessary
    std::vector<std::pair<int, int>> ToRemove;
    for (const auto it : E2T)
    {
        if (it.second.second == -1)
            ToRemove.emplace_back(it.first);
    }
    for (const auto& e : ToRemove)
        E2T.erase(e);

    // Get triangle to triangle adjacency list
    std::vector<std::pair<int, int>> T2T;
    T2T.reserve(E2T.size());
    for (const auto it : E2T)
        T2T.emplace_back(it.second);
    
    // Compute barycenters of triangles
    Eigen::MatrixXd Barycs;
    igl::barycenter(m_V, m_F, Barycs);

    // Return dual graph
    return rmt::Graph(Barycs, T2T);
}




rmt::Mesh rmt::Mesh::SubMesh(const std::vector<int>& Indices,
                             rmt::SubMeshAccessType AccType,
                             std::vector<int>& VIdx,
                             std::vector<int>& TIdx) const
{
    Eigen::VectorXi Potential;
    Eigen::VectorXi TriPotential;
    // If access type is by vertices, we must obtain the triangles
    if ((AccType & rmt::SubMeshAccessType::BY_VERTICES) != 0)
    {
        Potential.setZero(NumVertices());
        Potential(Indices).setConstant(1);
        TriPotential.setZero(NumTriangles());
        for (int i = 0; i < NumTriangles(); ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (Potential(GetTriangles()(i, j)) == 0)
                    continue;
                TriPotential[i]++;
            }
        }
        int Thresh = 1;
        if (AccType == rmt::SubMeshAccessType::BY_VERTICES_ALL)
            Thresh = 3;
        for (int i = 0; i < NumTriangles(); ++i)
            TriPotential[i] = TriPotential[i] >= Thresh ? 1 : 0;
    }
    // Otherwise, we already have the triangles
    else
    {
        TriPotential.setZero(NumTriangles());
        TriPotential(Indices).setConstant(1);
        Potential.setZero(NumVertices());
    }

    // Get the right vertex potential
    for (int i = 0; i < NumTriangles(); ++i)
    {
        if (TriPotential[i] == 0)
            continue;
        for (int j = 0; j < 3; ++j)
            Potential[GetTriangles()(i, j)] = 1;
    }

    // Extract vertex indices
    VIdx.clear();
    VIdx.reserve(Potential.sum());
    std::map<int, int> VMap;
    for (int i = 0; i < NumVertices(); ++i)
    {
        if (Potential[i] == 0)
            continue;
        
        VMap.emplace(i, (int)VIdx.size());
        VIdx.emplace_back(i);
    }

    // Make triangles
    TIdx.clear();
    TIdx.reserve(TriPotential.sum());
    for (int i = 0; i < NumTriangles(); ++i)
    {
        if (TriPotential[i] == 0)
            continue;
        
        TIdx.emplace_back(i);
    }

    Eigen::MatrixXi Tris = m_F(TIdx, Eigen::all);
    for (int i = 0; i < Tris.rows(); ++i)
    {
        for (int j = 0; j < 3; ++j)
            Tris(i, j) = VMap[Tris(i, j)];
    }

    // Return submesh
    return rmt::Mesh(m_V(VIdx, Eigen::all), Tris);
}



void rmt::Mesh::Export(const std::string& Filename) const
{
    rmt::ExportMesh(Filename, m_V, m_F);
}




Eigen::MatrixXd rmt::Mesh::UVMapping(rmt::UVMappingAlgorithm Algo,
                                     bool OptimalRotation) const
{
    Eigen::MatrixXd UV;
    if (Algo == rmt::UVMappingAlgorithm::NONE)
    {
        UV.setZero(NumVertices(), 2);
        return UV;
    }
    if (Algo == rmt::UVMappingAlgorithm::CONFORMAL)
    {
        if (!igl::lscm(GetVertices(), GetTriangles(), UV))
            throw std::runtime_error("Failed to find a conformal mapping.");
    }
    else if (Algo == rmt::UVMappingAlgorithm::TUTTE)
    {
        UV = rmt::TutteEmbedding(*this);
    }
    else
    {
        Eigen::VectorXi Bnd;
        igl::boundary_loop(GetTriangles(), Bnd);

        Eigen::MatrixXd BndUV(0, 2);
        if (Bnd.rows() > 0)
            igl::map_vertices_to_circle(GetVertices(), Bnd, BndUV);

        if (!igl::harmonic(GetVertices(), GetTriangles(), Bnd, BndUV, 1, UV))
            throw std::runtime_error("Failed to find a harmonic mapping.");

        if (Algo == rmt::UVMappingAlgorithm::ARAP)
        {
            igl::ARAPData ARAPData;
            ARAPData.with_dynamics = true;
            ARAPData.max_iter = 100;

            Eigen::VectorXi b = Eigen::VectorXi::Zero(0);
            Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0, 0);

            if (!igl::arap_precomputation(GetVertices(), GetTriangles(), 2, b, ARAPData))
                throw std::runtime_error("Failed ARAP precomputation.");
            if (!igl::arap_solve(bc, ARAPData, UV))
                throw std::runtime_error("Failed to find an ARAP mapping.");
        }
    }

    // Count flipped triangle and check if it's better to flip the entire unwrapping
    int NFlips = 0;
    for (int i = 0; i < NumTriangles(); ++i)
    {
        Eigen::Vector2d A, B, C;
        A = UV.row(m_F(i, 0));
        B = UV.row(m_F(i, 1));
        C = UV.row(m_F(i, 2));
        if ((B - A).cross(C - B) < 0)
            NFlips++;
    }
    if (NFlips > NumTriangles() - NFlips)
        UV.col(0) = -UV.col(0);

    // Rescale into [0, 1]^2
    UV = UV.array() - UV.minCoeff();
    UV = UV.array() / UV.maxCoeff();

    // If output it's fine as it is, we are done
    if (!OptimalRotation)
        return UV;

    // We have to search for the rotation that maximizes the total area
    // The rescaled version fit to the unit square has a scaling factor of 1 / (maxCoeff - minCoeff)
    // Iterate over a bunch of rotations and see what happens
    double MaxCoeff = 1.0;
    double BestTheta = 0.0;
    Eigen::Matrix2d R;
    R.setZero();
    UV = UV.array() - 0.5;
    double Theta = M_PI / 180.0; // The angle is always one degree, because we are iterating the rotations on UV
    // The matrix must be the trasposed rotation, because vectors are rowwise
    R(0, 0) = std::cos(Theta);
    R(0, 1) = std::sin(Theta);
    R(1, 0) = -R(0, 1);
    R(1, 1) = R(0, 0);
    for (int i = 0; i < 360; ++i)
    {
        UV *= R;
        double Coeff = UV.maxCoeff() - UV.minCoeff();
        if (Coeff < MaxCoeff)   // Check if less, beucase we are comparing 1 / Coeff
        {
            MaxCoeff = Coeff;
            BestTheta = (i - 1) * Theta;
        }
    }

    // Last rotation is 360, back to origin, so we can directly apply the rotation
    R(0, 0) = std::cos(BestTheta);
    R(0, 1) = std::sin(BestTheta);
    R(1, 0) = -R(0, 1);
    R(1, 1) = R(0, 0);
    UV *= R;
    
    // Rescale into [0, 1]^2
    UV = UV.array() - UV.minCoeff();
    UV = UV.array() / UV.maxCoeff();

    return UV;
}