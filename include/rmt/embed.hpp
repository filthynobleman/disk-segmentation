/**
 * @file        embed.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2025-09-18
 */
#pragma once


#include <rmt/mesh.hpp>


namespace rmt
{

enum TutteBoundary
{
    CIRCLE,
    SQUARE
};

    
Eigen::MatrixXd TutteEmbedding(const rmt::Mesh& M, 
                               rmt::TutteBoundary BndCondition = rmt::TutteBoundary::CIRCLE);
    
} // namespace rmt
