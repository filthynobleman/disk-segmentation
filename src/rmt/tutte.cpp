/**
 * @file        tutte.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2025-09-18
 */
#include <rmt/embed.hpp>

#include <Eigen/Sparse>

#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/adjacency_matrix.h>
#include <igl/harmonic.h>


#define _USE_MATH_DEFINES
#include <math.h>

void TutteSystem(const rmt::Mesh& M,
                 const std::vector<int>& BLoop,
                 Eigen::SparseMatrix<double>& S)
{
    S.resize(M.NumVertices(), M.NumVertices());
    S.setZero();
    S.reserve(10 * M.NumVertices());

    igl::adjacency_matrix(M.GetTriangles(), S);


    // Eigen::VectorXi Degree;
    // Degree.setZero(M.NumVertices());
    // for (int i = 0; i < M.NumEdges(); ++i)
    // {
    //     Degree[M.GetEdges()(i, 0)] += 1;
    //     Degree[M.GetEdges()(i, 1)] += 1;
    // }

    // for (int v : BLoop)
    //     Degree[v] = 0;

    // std::vector<Eigen::Triplet<double>> SEntries;
    // SEntries.reserve(10 * M.NumVertices());
    // for (int i = 0; i < M.NumEdges(); ++i)
    // {
    //     if (Degree[M.GetEdges()(i, 0)] != 0)
    //         SEntries.emplace_back(M.GetEdges()(i, 0), M.GetEdges()(i, 1), -1.0 / Degree[M.GetEdges()(i, 0)]);
    //     if (Degree[M.GetEdges()(i, 1)] != 0)
    //         SEntries.emplace_back(M.GetEdges()(i, 1), M.GetEdges()(i, 0), -1.0 / Degree[M.GetEdges()(i, 1)]);
    // }
    // for (int i = 0; i < M.NumVertices(); ++i)
    // {
    //     SEntries.emplace_back(i, i, 1);
    // }

}

Eigen::MatrixXd rmt::TutteEmbedding(const rmt::Mesh& M,
                                    rmt::TutteBoundary BndCondition)
{
    std::vector<int> BLoop;
    igl::boundary_loop(M.GetTriangles(), BLoop);
    Eigen::MatrixXd Bnd;
    Eigen::VectorXi BndVert = Eigen::VectorXi::Map(BLoop.data(), BLoop.size());
    igl::map_vertices_to_circle(M.GetVertices(), BndVert, Bnd);

    if (BndCondition == rmt::TutteBoundary::SQUARE)
    {
        for (int i = 0; i < Bnd.rows(); ++i)
        {
            double t = std::atan2(Bnd(i, 1), Bnd(i, 0));
            Bnd(i, 0) = std::cos(t);
            Bnd(i, 1) = std::sin(t);
            Bnd.row(i) /= std::max(std::abs(Bnd(i, 0)), std::abs(Bnd(i, 1)));
        }
    }

    Eigen::MatrixXd UV;
    igl::harmonic(M.GetTriangles(), BndVert, Bnd, 1, UV);
    return UV;

    // Eigen::SparseMatrix<double> TS;
    // TutteSystem(M, BLoop, TS);
    // Eigen::MatrixXd rhs;
    // rhs.setZero(M.NumVertices(), 2);
    // rhs(BLoop, Eigen::all) = Bnd;

    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> Solver;
    // Solver.compute(TS);
    // return Solver.solve(rhs);
}