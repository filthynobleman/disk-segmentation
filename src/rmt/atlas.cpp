/**
 * @file        atlas.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2025-08-01
 */
#include <rmt/atlas.hpp>

#include <igl/doublearea.h>

#include <queue>


rmt::Atlas::Atlas(const rmt::Mesh& Mesh) : m_Mesh(Mesh) { }
rmt::Atlas::~Atlas() { }


void rmt::Atlas::MakeRegions(const Eigen::VectorXi& Regions)
{
    m_Part = Regions;
    int NRegs = Regions.maxCoeff() + 1;
    m_Regions.clear();
    m_Regions.resize(NRegs);

    for (int i = 0; i < m_Mesh.NumTriangles(); ++i)
    {
        // Add the triangle to the corresponding region
        m_Regions[Regions[i]].AddFace(i);

        // Add the triangle vertices to the region
        for (int j = 0; j < 3; ++j)
            m_Regions[Regions[i]].AddVertex(m_Mesh.GetTriangles()(i, j));
        
        // Add the triangle's edges to the region
        for (int j = 0; j < 3; ++j)
            m_Regions[Regions[i]].AddEdge(m_Mesh.GetT2E()(i, j));
    }
}


void rmt::Atlas::OptimizeRegions(double Threshold)
{
    // Build adjacency graph of regions
    rmt::Graph DM = m_Mesh.AsDualGraph();
    std::set<std::pair<int, int>> Edges;
    for (int i = 0; i < DM.NumVertices(); ++i)
    {
        int nadj = DM.NumAdjacents(i);
        for (int jj = 0; jj < nadj; ++jj)
        {
            int j = DM.GetAdjacent(i, jj).first;
            if (m_Part[i] == m_Part[j])
                continue;
            std::pair<int, int> e{ m_Part[i], m_Part[j] };
            if (e.first > e.second)
                std::swap(e.first, e.second);
            Edges.insert(e);
        }
    }
    Eigen::MatrixXd V = Eigen::MatrixXd::Random(NumRegions(), 3);
    rmt::Graph G(V, Edges);

    Eigen::VectorXd Areas;
    igl::doublearea(m_Mesh.GetVertices(), m_Mesh.GetTriangles(), Areas);
    Areas *= 0.5;
    Eigen::VectorXd ELens = (m_Mesh.GetVertices()(m_Mesh.GetEdges()(Eigen::all, 1), Eigen::all) - m_Mesh.GetVertices()(m_Mesh.GetEdges()(Eigen::all, 0), Eigen::all)).rowwise().norm();

    bool SomethingMerged = false;
    // std::vector<std::pair<int, int>> Mergeable;
    do
    {
        // Greedy approach:
        // For each region, try to merge it with as many regions as possible
        // After the first mergeable region has been merged, skip everything else
        // If no unions could be found, we are in a local optimum
        SomethingMerged = false;
        // Scored approach:
        // Assign a score to each edge, and select the best scoring edge
        // Score must not be less than a threshold
        // Score = Boundary length / maximum area
        std::priority_queue<std::tuple<double, int, int>> Q;
        for (int i = 0; i < G.NumVertices(); ++i)
        {
            int nadj = G.NumAdjacents(i);
            std::vector<int> Fi;
            m_Regions[i].GetFaces(Fi);
            double Ai = Areas(Fi).sum();
            for (int jj = 0; jj < nadj; ++jj)
            {
                int j = G.GetAdjacent(i, jj).first;
                if (j < i)
                    continue;
                rmt::SurfaceRegion Intersection = m_Regions[i] * m_Regions[j];
                if (Intersection.EulerCharacteristic() != 1)
                    continue;
                m_Regions[j].GetFaces(Fi);
                double Aj = Areas(Fi).sum();
                Intersection.GetEdges(Fi);
                double Eij = ELens(Fi).sum();
                double Score = Eij / std::sqrt(std::max(Ai, Aj));
                if (Score < Threshold)
                    continue;
                Q.emplace(Score, i, j);
                // if (Intersection.NumEdges() < Threshold * std::sqrt(std::max(m_Regions[i].NumFaces(), m_Regions[j].NumFaces())))
                //     continue;
                // if (Eij < Threshold * std::sqrt(std::max(Ai, Aj)))
                //     continue;
                // rmt::SurfaceRegion Union = m_Regions[i] + m_Regions[j];
                // if (Union.EulerCharacteristic() != 1)
                //     continue;

                // m_Regions[i] = Union;
                // m_Regions.erase(m_Regions.begin() + j);
                // SomethingMerged = true;

                // auto E2 = Edges;
                // Edges.clear();
                // for (auto e : E2)
                // {
                //     auto ee = e;
                //     if (ee.first == j)
                //         ee.first = i;
                //     if (ee.second == j)
                //         ee.second = i;
                //     if (ee.first > j)
                //         ee.first--;
                //     if (ee.second > j)
                //         ee.second--;
                //     if (ee.first > ee.second)
                //         std::swap(ee.first, ee.second);
                //     if (ee.first == ee.second)
                //         continue;
                //     Edges.insert(ee);
                // }
                // V = V(Eigen::seq(0, V.rows() - 2), Eigen::all);
                // G = rmt::Graph(V, Edges);

                // break;
            }
            // if (SomethingMerged)
            //     break;
        }
        if (Q.empty())
            break;
        
        SomethingMerged = true;
        double Score;
        int i, j;
        std::tie(Score, i, j) = Q.top();
        rmt::SurfaceRegion Union = m_Regions[i] + m_Regions[j];
        if (Union.EulerCharacteristic() != 1)
            continue;

        m_Regions[i] = Union;
        m_Regions.erase(m_Regions.begin() + j);
        SomethingMerged = true;

        auto E2 = Edges;
        Edges.clear();
        int MaxE = 0;
        int MinE = std::numeric_limits<int>::max();
        for (auto e : E2)
        {
            auto ee = e;
            if (ee.first == j)
                ee.first = i;
            if (ee.second == j)
                ee.second = i;
            if (ee.first > j)
                ee.first--;
            if (ee.second > j)
                ee.second--;
            if (ee.first > ee.second)
                std::swap(ee.first, ee.second);
            if (ee.first == ee.second)
                continue;
            MaxE = std::max(MaxE, std::max(ee.first, ee.second));
            MinE = std::min(MinE, std::min(ee.first, ee.second));
            Edges.insert(ee);
        }
        V.conservativeResize(V.rows() - 1, V.cols());
        G = rmt::Graph(V, Edges);
    } while(SomethingMerged);

    std::vector<int> FID;
    for (int i = 0; i < NumRegions(); ++i)
    {
        GetRegion(i).GetFaces(FID);
        m_Part(FID).setConstant(i);
    }
}



int rmt::Atlas::NumRegions() const { return m_Regions.size(); }
const rmt::SurfaceRegion& rmt::Atlas::GetRegion(int i) const { return m_Regions[i]; }
bool rmt::Atlas::IsValid() const
{
    for (int i = 0; i < NumRegions(); ++i)
    {
        if (GetRegion(i).EulerCharacteristic() != 1)
            return false;
    }

    // Atlas cannot be valid if it does not map anything
    return NumRegions() > 0;
}

const Eigen::VectorXi& rmt::Atlas::GetPartitions() const { return m_Part; }



void rmt::Atlas::GenerateUVMap(rmt::UVMappingAlgorithm UVAlgo,
                               Eigen::MatrixXd& UV,
                               Eigen::MatrixXi& TUV,
                               bool Packing) const
{
    UV.resize(0, 2);
    TUV.resize(m_Mesh.NumTriangles(), 3);
    std::vector<int> Region;
    Region.reserve(m_Mesh.NumTriangles());
    std::vector<int> Tmp1, Tmp2;
    int SideSize = std::ceil(std::sqrt(NumRegions()));
    for (int i = 0; i < NumRegions(); ++i)
    {
        GetRegion(i).GetFaces(Region);
        rmt::Mesh SM = m_Mesh.SubMesh(Region, rmt::SubMeshAccessType::BY_TRIANGLES, Tmp1, Tmp2);

        int PrevUV = UV.rows();
        UV.conservativeResize(PrevUV + SM.NumVertices(), 2);
        Eigen::MatrixXd VertUV = SM.UVMapping(UVAlgo, UVAlgo != rmt::UVMappingAlgorithm::HARMONIC);
        VertUV = VertUV.array() - VertUV.minCoeff();
        if (Packing)
        {
            VertUV = VertUV.array() / VertUV.maxCoeff() * 0.975 + 0.0125;
            VertUV.col(0) = VertUV.col(0).array() + ((double)(i % SideSize));
            VertUV.col(1) = VertUV.col(1).array() + ((double)(i / SideSize));
            VertUV = VertUV.array() / SideSize;
        }
        else
            VertUV = VertUV.array() / VertUV.maxCoeff() * 0.975 + 0.0125;
        UV(Eigen::seq(PrevUV, UV.rows() - 1), Eigen::all) = VertUV;
        TUV(Region, Eigen::all) = SM.GetTriangles().array() + PrevUV;
    }
}