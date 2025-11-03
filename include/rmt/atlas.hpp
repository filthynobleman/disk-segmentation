/**
 * @file        atlas.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2025-08-01
 */
#pragma once

#include <rmt/mesh.hpp>
#include <rmt/region.hpp>


namespace rmt
{
    
class Atlas
{
private:
    const rmt::Mesh& m_Mesh;
    std::vector<rmt::SurfaceRegion> m_Regions;
    Eigen::VectorXi m_Part;

public:
    Atlas(const rmt::Mesh& M);
    ~Atlas();

    void MakeRegions(const Eigen::VectorXi& TriPart);
    void OptimizeRegions(double Threshold);
    
    int NumRegions() const;
    const rmt::SurfaceRegion& GetRegion(int i) const;
    bool IsValid() const;

    const Eigen::VectorXi& GetPartitions() const;

    void GenerateUVMap(rmt::UVMappingAlgorithm UVAlgo,
                       Eigen::MatrixXd& UV,
                       Eigen::MatrixXi& TUV,
                       bool Packing = false) const;
};


} // namespace rmt
