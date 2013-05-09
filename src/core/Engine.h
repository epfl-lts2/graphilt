/****************************************************************************
**
** Copyright (C) 2013 EPFL-LTS2
** Contact: Kirell Benzi (first.last@epfl.ch)
**
** This file is part of Graphilt.
**
**
** GNU General Public License Usage
** This file may be used under the terms of the GNU
** General Public License version 3.0 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPLv3 included in the
** packaging of this file.  Please review the following information to
** ensure the GNU General Public License version 3.0 requirements
** will be met: http://www.gnu.org/licenses/
**
****************************************************************************/

#ifndef GHT_ENGINE_H
#define GHT_ENGINE_H

#include <vector>
#include <viennacl/ocl/platform.hpp>
#include <viennacl/meta/result_of.hpp>
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/prod.hpp>

#include "util/log.h"

namespace ght {
namespace core {

class Engine
{
public:
    Engine() {}

    template<typename MatrixType, typename VectorType>
    static bool runCPU( const MatrixType& laplacian, const VectorType& signal,
              const std::vector<std::vector<float> >& coeff )
    {
        typedef typename viennacl::result_of::value_type<VectorType>::type VType;
        typedef typename viennacl::result_of::value_type<MatrixType>::type MType;


        return true;
    }

    template<typename MatrixType, typename VectorType, typename ScalarType>
    static bool runGPU( const MatrixType& laplacian, const VectorType& signal,
                 const std::vector<std::vector<ScalarType> >& coeff,
                 std::vector<std::vector<ScalarType> >& result )
    {
        // Init GPU values
        size_t nbScales = coeff.size();
        viennacl::vector<ScalarType> gSignal(signal.size());
        viennacl::vector<ScalarType> gTemp(signal.size());
        viennacl::compressed_matrix<ScalarType> gLaplacian(signal.size(), signal.size());
        viennacl::matrix<ScalarType> gCoeff(nbScales, coeff.at(0).size());
        result.resize(nbScales);

        // Copy
        viennacl::copy(signal, gSignal);
        viennacl::copy(laplacian, gLaplacian);
        viennacl::copy(coeff, gCoeff);

        // For earch scale
        for( size_t i = 0; i < nbScales; ++i ) {
            // For each scale coefficients
            size_t scaleCoeffs = coeff.at(i).size();
            result[i].resize(scaleCoeffs);
            for( size_t j = 0; j < scaleCoeffs; ++j ) {
                if( j == 0 ) {
                    gTemp = viennacl::linalg::prod(gLaplacian, gSignal);
                    gTemp *= gCoeff(i,j);
                }
                else {
                    gTemp += viennacl::linalg::prod(gLaplacian, gTemp);
                    gTemp *= gCoeff(i,j);
                }
            }
            // Copy back results
            viennacl::copy(gTemp, result[i]);
        }
        return true;
    }

    bool checkFitInGPUMem( size_t matrixSize, size_t signalSize )
    {
//        uint64_t inputSize = matrixSize + signalSize;
//        LOG(logDEBUG) << "Input size: " << inputSize / (1024 * 1024) << " MB";
//        LOG(logDEBUG) << "Max allocable: " << viennacl::ocl::current_device().max_allocable_memory() / (1024 * 1024) << " MB";

//        if( inputSize > viennacl::ocl::current_device().max_allocable_memory() ) {
//            return false;
//        }
        return true;
    }

};

} // end namespace core
} // end namespace ght

#endif // GHT_ENGINE_H
