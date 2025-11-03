/**
 * @file        io.hpp
 * 
 * @brief       Functions for loading and exporting data.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-26
 */
#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <rmt/graph.hpp>

namespace rmt
{

bool LoadMesh(const std::string& Filename,
              Eigen::MatrixXd& V,
              Eigen::MatrixXi& F);

bool ExportMesh(const std::string& Filename,
                const Eigen::MatrixXd& V,
                const Eigen::MatrixXi& F);


bool ExportWeightmap(const std::string& Filename,
                     const Eigen::SparseMatrix<double>& WM);

bool ExportGraph(const std::string& Filename,
                 const rmt::Graph& G,
                 const Eigen::MatrixXd& V);


template<typename T>
bool ExportList(const std::string& Filename,
                const std::vector<T>& List)
{
    std::ofstream Stream;
    Stream.open(Filename, std::ios::out);
    if (!Stream.is_open())
        return false;

    for (const auto& v : List)
        Stream << v << '\n';

    Stream.close();
    return true;
}

template<typename T>
bool LoadList(const std::string& Filename,
              std::vector<T>& List)
{
    std::ifstream Stream;
    Stream.open(Filename, std::ios::in);
    if (!Stream.is_open())
        return false;

    std::string Line;
    while (!Stream.eof())
    {
        std::getline(Stream, Line);
        if (Line.empty())
            continue;
        std::istringstream ss(Line);
        T v;
        ss >> v;
        List.emplace_back(v);
    }

    Stream.close();
    return true;
}
    
} // namespace rmt
