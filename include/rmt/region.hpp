/**
 * @file        region.hpp
 * 
 * @brief       A class representing a region over a triangulated surface.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <rmt/utils.hpp>

#include <set>


namespace rmt
{
    
class SurfaceRegion
{
private:
    std::set<int> m_Verts;
    std::set<int> m_Edges;
    std::set<int> m_Faces;

public:
    SurfaceRegion();
    SurfaceRegion(const rmt::SurfaceRegion& SR);
    SurfaceRegion(rmt::SurfaceRegion&& SR);
    rmt::SurfaceRegion& operator=(const rmt::SurfaceRegion& SR);
    rmt::SurfaceRegion& operator=(rmt::SurfaceRegion&& SR);
    ~SurfaceRegion();

    int EulerCharacteristic() const;

    void AddFace(int fid);
    void AddEdge(int eid);
    void AddVertex(int vid);

    int NumVertices() const;
    int NumEdges() const;
    int NumFaces() const;
    void GetVertices(std::vector<int>& Vertices) const;
    void GetEdges(std::vector<int>& Edges) const;
    void GetFaces(std::vector<int>& Faces) const;

    friend rmt::SurfaceRegion operator+(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2);
    friend rmt::SurfaceRegion operator*(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2);
};

rmt::SurfaceRegion operator+(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2);
rmt::SurfaceRegion operator*(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2);



} // namespace rmt
