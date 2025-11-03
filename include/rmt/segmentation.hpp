/**
 * @file        segmentation.hpp
 * 
 * @brief       Segmentation methods.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2025-08-28
 */
#pragma once


#include <rmt/mesh.hpp>


namespace rmt
{

Eigen::VectorXi ShapeDiameterSegmentation(const rmt::Mesh& M);
Eigen::VectorXi SpectralSegmentation(const rmt::Mesh& M,
                                     int k = 10,
                                     double t = 1e-3);

} // namespace rmt
