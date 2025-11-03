/**
 * @file        spectralseg.cpp
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
#include <rmt/voronoifps.hpp>

#include <Eigen/Sparse>

#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/edge_lengths.h>
#include <igl/squared_edge_lengths.h>
#include <igl/eigs.h>

void stiffmatrix(const Eigen::MatrixXd V,
                 const Eigen::MatrixXi T,
                 Eigen::SparseMatrix<double>& S)
{
    S.resize(V.rows(), V.rows());
    S.reserve(10 * V.rows());

    std::vector<Eigen::Triplet<double>> Angles;
    Eigen::VectorXd Diag;
    Diag.setZero(V.rows());
    for (int i0 = 0; i0 < 3; ++i0)
    {
        int i1 = (i0 + 1) % 3;
        int i2 = (i1 + 1) % 3;

        Eigen::MatrixXd E01 = V(T(Eigen::all, i1), Eigen::all) - V(T(Eigen::all, i0), Eigen::all);
        Eigen::MatrixXd E02 = V(T(Eigen::all, i2), Eigen::all) - V(T(Eigen::all, i0), Eigen::all);
        E01.rowwise().normalize();
        E02.rowwise().normalize();

        for (int j = 0; j < T.rows(); ++j)
        {
            // double A = -0.5 / std::tan(std::acos(E01.row(j).dot(E02.row(j))));
            // This should be equivalent and more efficient
            // A = cos(alpha) --> sqrt(1 - A*A) = sin(alpha) --> A / sqrt(1 - A*A) = cos(alpha) / sin(alpha) = cot(alpha)
            double A = E01.row(j).dot(E02.row(j));
            A = -0.5 * A / std::sqrt(1 - A * A);
            Angles.emplace_back(T(j, i1), T(j, i2), A);
            Angles.emplace_back(T(j, i2), T(j, i1), A);
            Diag[i1] -= A;
            Diag[i2] -= A;
        }
    }

    for (int i = 0; i < V.rows(); ++i)
        Angles.emplace_back(i, i, Diag[i]);
    
    S.setZero();
    S.setFromTriplets(Angles.begin(), Angles.end());
}

Eigen::VectorXi rmt::SpectralSegmentation(const rmt::Mesh& M,
                                          int k,
                                          double t)
{
    Eigen::SparseMatrix<double> S, A;
    Eigen::MatrixXd V = M.GetVertices();
    igl::cotmatrix(V, M.GetTriangles(), S);
    igl::massmatrix(V, M.GetTriangles(), igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, A);

    Eigen::VectorXd EL2 = (V(M.GetEdges()(Eigen::all, 0), Eigen::all) - V(M.GetEdges()(Eigen::all, 1), Eigen::all)).rowwise().squaredNorm();
    t = t * EL2.mean();

    Eigen::MatrixXd U;
    U.resize(M.NumVertices(), 3);
    for (int i = 0; i < k; ++i)
    {
        double s = (V.colwise().maxCoeff() - V.rowwise().minCoeff()).norm();
        V = V.array() / s;
        Eigen::Vector3d c = 0.5 * (V.colwise().maxCoeff() + V.colwise().minCoeff());
        V -= c.transpose().replicate(V.rows(), 1);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> Solver;
        U = A * V;
        A -= t * S;
        Solver.analyzePattern(A);
        Solver.factorize(A);
        V = Solver.solve(U);

        igl::massmatrix(V, M.GetTriangles(), igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, A);
    }

    igl::cotmatrix(V, M.GetTriangles(), S);

    Eigen::MatrixXd EVecs;
    Eigen::VectorXd EVals;
    // igl::eigs(S, A, 2, igl::EigsType::EIGS_TYPE_SM, EVecs, EVals);
    if (!igl::eigs(S, A, 2, igl::EigsType::EIGS_TYPE_SM, EVecs, EVals))
        throw std::runtime_error("Cannot solve eigendecomposition.");
    Eigen::VectorXd Phi = EVecs.col(1);
    
    std::vector<bool> IsSmaller;
    IsSmaller.resize(M.NumVertices(), true);
    std::vector<bool> IsLarger;
    IsLarger.resize(M.NumVertices(), true);
    for (int i = 0; i < M.NumEdges(); ++i)
    {
        int v0 = M.GetEdges()(i, 0);
        int v1 = M.GetEdges()(i, 1);
        if (Phi(v0) > Phi(v1))
        {
            IsLarger[v1] = false;
            IsSmaller[v0] = false;
        }
        else if (Phi(v0) < Phi(v1))
        {
            IsLarger[v0] = false;
            IsSmaller[v1] = false;
        }
        else
        {
            IsLarger[v0] = false;
            IsSmaller[v0] = false;
            IsLarger[v1] = false;
            IsSmaller[v1] = false;
        }
    }

    Eigen::VectorXi SaddleTris;
    SaddleTris.setZero(M.NumVertices());
    for (int i = 0; i < M.NumTriangles(); ++i)
    {
        for (int j0 = 0; j0 < 3; ++j0)
        {
            int j1 = j0 + 1;
            if (j1 == 3)
                j1 = 0;
            int j2 = j1 + 1;
            if (j2 == 3)
                j2 = 0;
            
            int v0 = M.GetTriangles()(i, j0);
            int v1 = M.GetTriangles()(i, j1);
            int v2 = M.GetTriangles()(i, j2);
            
            if (Phi(v0) <= std::min(Phi(v1), Phi(v2)))
                continue;
            if (Phi(v0) >= std::max(Phi(v1), Phi(v2)))
                continue;
            
            SaddleTris[v0] += 1;
        }
    }


    std::vector<int> CriticalPoints;
    for (int i = 0; i < M.NumVertices(); ++i)
    {
        if (IsSmaller[i])
        {
            CriticalPoints.push_back(i);
            continue;
        }
        if (IsLarger[i])
        {
            CriticalPoints.push_back(i);
            continue;
        }
        if (SaddleTris[i] >= 4)
        {
            CriticalPoints.push_back(i);
            continue;
        }
    }


    rmt::VoronoiPartitioning VPart(M);
    for (int cp : CriticalPoints)
        VPart.AddSample(cp);


    Eigen::VectorXi Partition;
    Partition.resize(M.NumTriangles());
    for (int i = 0; i < M.NumTriangles(); ++i)
    {
        int p0 = VPart.GetPartition(M.GetTriangles()(i, 0));
        int p1 = VPart.GetPartition(M.GetTriangles()(i, 1));
        int p2 = VPart.GetPartition(M.GetTriangles()(i, 2));
        // // Same partition for all vertices --> triangle is in that partition
        // if (p0 == p1 && p1 == p2)
        // {
        //     Partition[i] = p0;
        //     continue;
        // }
        // // Two vertices have the same partition --> triangle goes there
        // if (p0 == p1 || p0 == p2)
        // {
        //     Partition[i] = p0;
        //     continue;
        // }
        // if (p1 == p2)
        // {
        //     Partition[i] = p1;
        //     continue;
        // }
        // // All different paritions --> just choose one
        // Partition[i] = p0;
        // All this stuff can be simplified as
        if (p1 == p2)
            Partition[i] = p1;
        else
            Partition[i] = p0;
    }

    return Partition;
}