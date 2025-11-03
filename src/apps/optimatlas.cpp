/**
 * @file        remap.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2025-07-29
 */
#include <rmt/rmt.hpp>
#include <iostream>
#include <chrono>
#include <filesystem>
#include <cctype>
#include <algorithm>
#include <cut/time/time.hpp>

#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/writeOBJ.h>


struct
{
    std::string InputMesh;
    std::string OutputMesh;
    int NumRegions;
    int SubRegions;
    double IntersectThreshold;
    rmt::UVMappingAlgorithm Algorithm;
    bool UVPacking;
    bool WithNormals;
    bool DebugOutput;

    int Verbosity;
} Args;

void ParseArgs(int argc, const char* const argv[]);
void Usage(const std::string& argv0);


int main(int argc, const char* const argv[])
{
    ParseArgs(argc, argv);

    if (Args.Verbosity > 1)
        cut::Timer::AttachTimer("mesh_loader");
    rmt::Mesh M(Args.InputMesh);
    M.ComputeEdgesAndIncidence();
    Eigen::MatrixXd N;
    Eigen::MatrixXi FN;
    Eigen::MatrixXd UV;
    Eigen::MatrixXi FTC;
    if (Args.WithNormals)
    {
        igl::per_vertex_normals(M.GetVertices(), M.GetTriangles(), N);
        FN = M.GetTriangles();
    }
    else
    {
        N.setZero(0, 3);
        FN.setZero(0, 3);
    }
    if (Args.Verbosity > 1)
    {
        std::cout << "Mesh loaded in ";
        std::cout << cut::Timer::GetTimer("mesh_loader").GetTime();
        std::cout << " seconds." << std::endl;
    }

    if (Args.Verbosity > 0);
        cut::Timer::AttachTimer("total_runtime");

    if (Args.Verbosity > 1)
        cut::Timer::AttachTimer("dual_graph");
    rmt::Graph DG = M.AsDualGraph();
    if (Args.Verbosity > 1)
    {
        std::cout << "Dual graph construction in ";
        std::cout << cut::Timer::GetTimer("dual_graph").GetTime();
        std::cout << " seconds." << std::endl;
    }

    if (Args.Verbosity > 1)
        cut::Timer::AttachTimer("voronoi");
    rmt::VoronoiPartitioning VPart(DG);
    Args.NumRegions = std::min(Args.NumRegions, DG.NumVertices());
    while (VPart.NumSamples() < Args.NumRegions)
        VPart.AddSample(VPart.FarthestVertex());
    if (Args.Verbosity > 1)
    {
        std::cout << "Voronoi partitioning in ";
        std::cout << cut::Timer::GetTimer("voronoi").GetTime();
        std::cout << " seconds." << std::endl;
    }
        
    if (Args.Verbosity > 1)
        cut::Timer::AttachTimer("optimize_regions");
    rmt::Atlas Atlas(M);
    Atlas.MakeRegions(VPart.GetPartitions());
    if (Args.DebugOutput)
    {
        Atlas.GenerateUVMap(rmt::UVMappingAlgorithm::TUTTE, UV, FTC);
        std::string vorout = Args.OutputMesh;
        vorout = vorout.substr(0, vorout.rfind(".obj")) + "-voronoi.obj";
        igl::writeOBJ(vorout, M.GetVertices(), M.GetTriangles(), N, FN, UV, FTC);
    }
    std::vector<int> Faces;
    Faces.reserve(M.NumTriangles());
    std::vector<int> Tmp1, Tmp2;
    int NParts = VPart.NumSamples();
    Eigen::VectorXi Partitions = VPart.GetPartitions();
    while (!Atlas.IsValid())
    {
        // Identify non-disk regions
        for (int i = 0; i < Atlas.NumRegions(); ++i)
        {
            if (Atlas.GetRegion(i).EulerCharacteristic() == 1)
                continue;
            
            // Try subsampling these regions
            Atlas.GetRegion(i).GetFaces(Faces);
            rmt::Mesh SM = M.SubMesh(Faces, rmt::SubMeshAccessType::BY_TRIANGLES, Tmp1, Tmp2);
            rmt::Graph SDG = SM.AsDualGraph();
            rmt::VoronoiPartitioning SVP(SDG);
            while (SVP.NumSamples() < std::min(Args.SubRegions, SDG.NumVertices()))
                SVP.AddSample(SVP.FarthestVertex());
            Eigen::VectorXi SubParts = SVP.GetPartitions();
            for (int j = 0; j < Faces.size(); ++j)
            {
                if (SubParts[j] != 0)
                    Partitions[Faces[j]] = SubParts[j] + NParts - 1;
            }
            NParts += SVP.NumSamples() - 1;
        }

        Atlas.MakeRegions(Partitions);
    }
    if (Args.DebugOutput)
    {
        Atlas.GenerateUVMap(rmt::UVMappingAlgorithm::TUTTE, UV, FTC);
        std::string disks = Args.OutputMesh;
        disks = disks.substr(0, disks.rfind(".obj")) + "-disks.obj";
        igl::writeOBJ(disks, M.GetVertices(), M.GetTriangles(), N, FN, UV, FTC);
    }
    Atlas.OptimizeRegions(Args.IntersectThreshold);
    if (Args.Verbosity > 1)
    {
        std::cout << "Optimized regions in ";
        std::cout << cut::Timer::GetTimer("optimize_regions").GetTime();
        std::cout << " seconds." << std::endl;
    }


    if (Args.Verbosity > 1)
        cut::Timer::AttachTimer("unwrapping");
    try
    {
        Atlas.GenerateUVMap(Args.Algorithm, UV, FTC, Args.UVPacking);
    }
    catch(const std::exception& e)
    {
        return EXIT_FAILURE;
    }
    // Sometimes, NaN values are generated. Dirty solution: replace with zeros
    for (int i = 0; i < UV.rows(); ++i)
    {
        for (int j = 0; j < UV.cols(); ++j)
        {
            if (std::isnan(UV(i, j)))
                UV(i, j) = 0;
        }
    }
    if (Args.Verbosity > 1)
    {
        std::cout << "Unwrapped " << Atlas.NumRegions() << " regions in ";
        std::cout << cut::Timer::GetTimer("unwrapping").GetTime();
        std::cout << " seconds." << std::endl;
    }

    double Runtime = 0.0;
    if (Args.Verbosity > 0)
        Runtime = cut::Timer::GetTimer("total_runtime").GetTime();
    
    if (Args.Verbosity > 1)
        cut::Timer::AttachTimer("write_mesh");
    igl::writeOBJ(Args.OutputMesh, M.GetVertices(), M.GetTriangles(), N, FN, UV, FTC);
    if (Args.Verbosity > 1)
    {
        std::cout << "Mesh exported in ";
        std::cout << cut::Timer::GetTimer("write_mesh").GetTime();
        std::cout << " seconds." << std::endl;;
    }

    if (Args.Verbosity > 0)
    {
        std::cout << "Total runtime is ";
        std::cout << Runtime;
        std::cout << " seconds." << std::endl;
    }

    return 0;
}

void ParseArgs(int argc, const char* const argv[])
{
    // Set default parameters
    Args.Algorithm = rmt::UVMappingAlgorithm::TUTTE;
    Args.Verbosity = 1;
    Args.NumRegions = 10;
    Args.SubRegions = 5;
    Args.IntersectThreshold = 0.76; // About the ratio between pentagon side and sqrt(pentagon area)
    Args.UVPacking = false;
    Args.WithNormals = false;
    Args.DebugOutput = false;

    // Unset parameters
    Args.InputMesh = "";
    Args.OutputMesh = "";

    for (int i = 0; i < argc; ++i)
    {
        // Only required input is the input mesh
        if (argv[i][0] != '-')
        {
            Args.InputMesh = argv[i];
            continue;
        }
        std::string Argvi = argv[i];
        // Number of samples
        if (Argvi == "-n" || Argvi == "--num-regions")
        {
            Args.NumRegions = std::stoi(argv[++i]);
            continue;
        }
        // Number of subsamples
        if (Argvi == "-s" || Argvi == "--sub-regions")
        {
            Args.SubRegions = std::stoi(argv[++i]);
            continue;
        }
        // Threshold
        if (Argvi == "-t" || Argvi == "--threshold")
        {
            Args.IntersectThreshold = std::stod(argv[++i]);
            continue;
        }
        // Output file
        if (Argvi == "-o" || Argvi == "--output")
        {
            Args.OutputMesh = argv[++i];
            continue;
        }
        // Algorithm
        if (Argvi == "-a" || Argvi == "--algorithm")
        {
            Argvi = argv[++i];
            std::transform(Argvi.begin(), Argvi.end(), Argvi.begin(),
                           [](unsigned char c) { return std::tolower(c); });

            if (Argvi == "harmonic")
                Args.Algorithm = rmt::UVMappingAlgorithm::HARMONIC;
            else if (Argvi == "conformal")
                Args.Algorithm = rmt::UVMappingAlgorithm::CONFORMAL;
            else if (Argvi == "arap")
                Args.Algorithm = rmt::UVMappingAlgorithm::ARAP;
            else if (Argvi == "tutte")
                Args.Algorithm = rmt::UVMappingAlgorithm::TUTTE;
            continue;
        }
        // Verbosity
        if (Argvi == "-v" || Argvi == "--verbosity")
        {
            Args.Verbosity = std::stoi(argv[++i]);
            continue;
        }
        // Packing
        if (Argvi == "-p" || Argvi == "--packing")
        {
            Args.UVPacking = true;
            continue;
        }
        // Normals in output file
        if (Argvi == "--normals")
        {
            Args.WithNormals = true;
            continue;
        }
        // Debug mode
        if (Argvi == "--debug")
        {
            Args.DebugOutput = true;
            continue;
        }
        // Help
        if (Argvi == "-h" || Argvi == "--help")
        {
            Usage(argv[0]);
            exit(EXIT_SUCCESS);
        }
        // Any other option is not acceptable
        std::cerr << "Unrecongnized option " << Argvi << '.' << std::endl;
        std::cerr << "The correct syntax is as follows." << std::endl;
        Usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // If input mesh is unset, invalid syntax
    if (Args.InputMesh.empty())
    {
        std::cerr << "No input mesh provided." << std::endl;
        std::cerr << "The correct syntax is as follows." << std::endl;
        Usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // If no output mesh is given, it's the same as input mesh, but with "-uv" at the end of the filename
    if (Args.OutputMesh.empty())
    {
        Args.OutputMesh = Args.InputMesh.substr(0, Args.InputMesh.rfind('.'));
        Args.OutputMesh += "-uv.obj";
    }

    // Output mesh must be in OBJ format
    std::string Ext = Args.OutputMesh.substr(Args.OutputMesh.rfind('.'));
    std::transform(Ext.begin(), Ext.end(), Ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (Ext != ".obj")
    {
        std::cerr << "WARNING: output mesh is not in OBJ format. Do you want to replace the extension with \'.obj\' and continue? y/n" << std::endl;
        char choice;
        do
        {
            std::cin >> choice;
            switch (choice)
            {
            case 'y':
            case 'Y':
                break;
            
            case 'n':
            case 'N':
                std::cerr << "Terminating the application." << std::endl;
                exit(EXIT_FAILURE);

            default:
                std::cerr << "Invalid choice." << std::endl;
            }
        } while (choice != 'y' && choice != 'Y');
        Args.OutputMesh = Args.OutputMesh.substr(0, Args.OutputMesh.rfind('.')) + ".obj";
    }
}


void Usage(const std::string& argv0)
{
    std::cout << "Usage:" << std::endl;
    std::cout << "  " << argv0 << " input_mesh [-n num_regions] [-s sub_regions] [-t threshold] [-o output_mesh] [-a algorithm] [-v verbosity] [-p]" << std::endl;
    std::cout << "  " << argv0 << " input_mesh [--num-regions num_regions] [--sub-regions sub_regions] [--threshold threshold] [--output output_mesh] [--algorithm algorithm] [--verbosity verbosity] [--packing]" << std::endl;
    std::cout << "  " << argv0 << "-h" << std::endl;
    std::cout << "  " << argv0 << "--help" << std::endl;
    std::cout << std::endl;
    std::cout << "    input_mesh is the path to a triangulated mesh." << std::endl;
    std::cout << "    num_regions is the number of regions that must subdivide the mesh. Default is 10." << std::endl;
    std::cout << "    sub_regions is the number of regions that must subdivide non topological disks. Default is 5." << std::endl;
    std::cout << "    threshold is the threshold parameter that determines if two regions can be merged. Default is 0.5" << std::endl;
    std::cout << "    algorithm is the UV unwrapping algorithm for each region. Acceptable values are \'harmonic\', \'conformal\', \'arap\'. Default is \'harmonic\'." << std::endl;
    std::cout << "    verbosity is the verbosity level of the output, ranging from 0 (no output) to 2 (runtime of each step). Default is 1 (total runtime)." << std::endl;
    std::cout << "    -p|--packing option forces each UV island to not overlap with the others. By default, island are all rescaled to [0, 1]^2." << std::endl;
    std::cout << "    -h|--help option prints this help message." << std::endl;
}