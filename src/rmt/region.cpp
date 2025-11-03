/**
 * @file        region.cpp
 * 
 * @brief       Implementation of class rmt::SurfaceRegion.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#include <rmt/region.hpp>


rmt::SurfaceRegion::SurfaceRegion() { }

rmt::SurfaceRegion::SurfaceRegion(const rmt::SurfaceRegion& SR) 
{ 
    m_Verts = SR.m_Verts;
    m_Edges = SR.m_Edges;
    m_Faces = SR.m_Faces;
}

rmt::SurfaceRegion& rmt::SurfaceRegion::operator=(const rmt::SurfaceRegion& SR) 
{
    m_Verts = SR.m_Verts;
    m_Edges = SR.m_Edges;
    m_Faces = SR.m_Faces;
    return *this;
}

rmt::SurfaceRegion::SurfaceRegion(rmt::SurfaceRegion&& SR) 
{ 
    m_Verts = std::move(SR.m_Verts);
    m_Edges = std::move(SR.m_Edges);
    m_Faces = std::move(SR.m_Faces);
}

rmt::SurfaceRegion& rmt::SurfaceRegion::operator=(rmt::SurfaceRegion&& SR) 
{
    m_Verts = std::move(SR.m_Verts);
    m_Edges = std::move(SR.m_Edges);
    m_Faces = std::move(SR.m_Faces);
    return *this;
}

rmt::SurfaceRegion::~SurfaceRegion() { }


int rmt::SurfaceRegion::EulerCharacteristic() const
{
    // return NumVertices() + NumFaces() - NumEdges();
    return m_Verts.size() + m_Faces.size() - m_Edges.size();
}


void rmt::SurfaceRegion::AddFace(int fid)          { m_Faces.insert(fid); }
void rmt::SurfaceRegion::AddVertex(int vid)        { m_Verts.insert(vid); }
void rmt::SurfaceRegion::AddEdge(int eid)          { m_Edges.insert(eid); }

int rmt::SurfaceRegion::NumFaces() const          { return m_Faces.size(); }
int rmt::SurfaceRegion::NumVertices() const       { return m_Verts.size(); }
int rmt::SurfaceRegion::NumEdges() const          { return m_Edges.size(); }


void rmt::SurfaceRegion::GetVertices(std::vector<int>& Vertices) const
{
    Vertices.clear();
    Vertices.insert(Vertices.begin(), m_Verts.begin(), m_Verts.end());
}

void rmt::SurfaceRegion::GetEdges(std::vector<int>& Edges) const
{
    Edges.clear();
    Edges.insert(Edges.begin(), m_Edges.begin(), m_Edges.end());
}
void rmt::SurfaceRegion::GetFaces(std::vector<int>& Faces) const
{
    Faces.clear();
    Faces.insert(Faces.begin(), m_Faces.begin(), m_Faces.end());
}


rmt::SurfaceRegion rmt::operator+(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2)
{
    rmt::SurfaceRegion SR;
    std::set_union(SR1.m_Verts.begin(), SR1.m_Verts.end(),
                   SR2.m_Verts.begin(), SR2.m_Verts.end(),
                   std::inserter(SR.m_Verts, SR.m_Verts.begin()));
    std::set_union(SR1.m_Edges.begin(), SR1.m_Edges.end(),
                   SR2.m_Edges.begin(), SR2.m_Edges.end(),
                   std::inserter(SR.m_Edges, SR.m_Edges.begin()));
    std::set_union(SR1.m_Faces.begin(), SR1.m_Faces.end(),
                   SR2.m_Faces.begin(), SR2.m_Faces.end(),
                   std::inserter(SR.m_Faces, SR.m_Faces.begin()));

    return SR;
}

rmt::SurfaceRegion rmt::operator*(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2)
{
    rmt::SurfaceRegion SR;
    std::set_intersection(SR1.m_Verts.begin(), SR1.m_Verts.end(),
                          SR2.m_Verts.begin(), SR2.m_Verts.end(),
                          std::inserter(SR.m_Verts, SR.m_Verts.begin()));
    std::set_intersection(SR1.m_Edges.begin(), SR1.m_Edges.end(),
                          SR2.m_Edges.begin(), SR2.m_Edges.end(),
                          std::inserter(SR.m_Edges, SR.m_Edges.begin()));
    // std::set_intersection(SR1.m_Faces.begin(), SR1.m_Faces.end(),
    //                       SR2.m_Faces.begin(), SR2.m_Faces.end(),
    //                       std::inserter(SR.m_Faces, SR.m_Faces.begin()));

    return SR;
}