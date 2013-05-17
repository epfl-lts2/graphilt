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

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/lu.hpp>

// Must be set if you want to use ViennaCL algorithms on ublas objects
#define VIENNACL_WITH_UBLAS 1

#include <viennacl/ocl/platform.hpp>
#include <viennacl/meta/result_of.hpp>
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/prod.hpp>
#include <viennacl/matrix_proxy.hpp>
#include <viennacl/linalg/power_iter.hpp>
#include <viennacl/linalg/lanczos.hpp>

#include "util/log.h"
#include "util/benchmark.h"
#include "io/matrix-io.h"

namespace bn = boost::numeric;
static const double kERROR_FACTOR = 1.01;

namespace ght {
namespace core {

// Template typedefs
// boost
template <typename ScalarType>
using UMatrix = bn::ublas::compressed_matrix<ScalarType>;
template <typename ScalarType>
using UVector = bn::ublas::vector<ScalarType>;
// Std
template <typename ScalarType>
using SVecvec = std::vector<std::vector<ScalarType> >;
template <typename ScalarType>
using SVec = std::vector<ScalarType>;
template <typename ScalarType>
using SCMatrix = std::vector<std::map<uint32_t, ScalarType> >;
template <typename ScalarType>
using VCMatrix = viennacl::compressed_matrix<ScalarType>;
template <typename ScalarType>
using VCMatrixPtr = std::shared_ptr<VCMatrix<ScalarType> >;

class Engine
{
public:
    Engine() {}

    template <typename ScalarType>
    static bool runNaiveCPU( const UMatrix<ScalarType>& laplacian,
                             const UVector<ScalarType>& signal,
                 const SVecvec<ScalarType>& coeff,
                 std::vector<std::vector<ScalarType> >& result )
    {
        size_t nbScales = coeff.size();
        bn::ublas::vector<ScalarType> temp(signal.size());
        bn::ublas::vector<ScalarType> res(signal.size());
        result.resize(nbScales);
        // For earch scale
        for( size_t i = 0; i < nbScales; ++i ) {
            result[i].resize(temp.size());
            // For each scale coefficients
            size_t scaleCoeffs = coeff.at(i).size();
            res = signal * coeff[i][0]; // Init bias
            for( size_t j = 1; j < scaleCoeffs; ++j ) {
                if( j == 1 ) {
                    temp = bn::ublas::prod(laplacian, signal);
                }
                else { // Recurrence relation
                    temp = bn::ublas::prod(laplacian, temp);
                }
                res += temp * coeff[i][j];
            }
            std::copy(res.begin(), res.end(), result[i].begin());
        }
        return true;
    }

    template<typename MatrixType, typename VectorType, typename ScalarType>
    static bool runNaiveGPU( const MatrixType& laplacian, const VectorType& signal,
                 const SVecvec<ScalarType>& coeff, SVecvec<ScalarType>& result )
    {
        // Init GPU values
        size_t nbScales = coeff.size();
        size_t signalSize = signal.size();
        viennacl::vector<ScalarType> gSignal(signalSize);
        viennacl::vector<ScalarType> gCombLaplacian(signalSize);
        viennacl::vector<ScalarType> gResult(signalSize);
        viennacl::compressed_matrix<ScalarType> gLaplacian(signalSize, signalSize);
        viennacl::matrix<ScalarType> gCoeff(nbScales, coeff.at(0).size());
        result.resize(nbScales);
        // Copy
        try {
            viennacl::copy(signal, gSignal);
            viennacl::copy(laplacian, gLaplacian);
            viennacl::copy(coeff, gCoeff);
        } catch( std::exception& e ) {
            LOG(logERROR) << "Engine::runNaiveGPU: " << e.what();
            return false;
        }
        // For earch scale
        for( size_t i = 0; i < nbScales; ++i ) {
            // For each scale coefficients
            result[i].resize(signalSize);
            size_t scaleCoeffs = coeff.at(i).size();
            gResult = gSignal * gCoeff(i,0); // first term bias
            for( size_t j = 1; j < scaleCoeffs; ++j ) {
                if( j == 1 ) {
                    gCombLaplacian = viennacl::linalg::prod(gLaplacian, gSignal);
                }
                else { // Recurrence relation
                    gCombLaplacian = viennacl::linalg::prod(gLaplacian, gCombLaplacian);
                }
                gResult += gCombLaplacian * gCoeff(i,j);
            }

            // Copy back results
            viennacl::copy(gResult, result[i]);
        }
        return true;
    }

    template<typename MatrixType, typename ScalarType>
    static bool runGPU( const MatrixType& laplacian, const SVec<ScalarType>& signal,
                 const SVecvec<ScalarType>& coeff, SVecvec<ScalarType>& result )
    {
        using namespace viennacl;
        typedef vector<ScalarType> VCLVector;
        typedef matrix<ScalarType> VCLMatrix;
        // Init GPU values
        size_t nbScales = coeff.size();
        size_t signalSize = signal.size();
        size_t filterOrder = coeff.at(0).size(); // tmp

        VCLVector gCombLaplacian(signalSize);
        VCMatrix<ScalarType> gLaplacian(signalSize, signalSize);
        VCLMatrix gCoeff(nbScales, filterOrder);
        VCLMatrix gResult(nbScales, signalSize);
        VCLVector gSignalVector(signalSize);

        // Copy to GPU
        fast_copy(signal, gSignalVector);
        copy(laplacian, gLaplacian);
        copy(coeff, gCoeff);

        // First coeff
        range coeffRow(0, coeff.size());
        range coeffCol(0, 1); // 1 column only
        matrix_range<VCLMatrix> currentCoeffs = project(gCoeff, coeffRow, coeffCol);
        // Init result matrix
        gResult = linalg::prod(currentCoeffs, VCLMatrix(gSignalVector.handle().opencl_handle(), signalSize, 1));

        // 2nd coeff, first Laplacian
        gCombLaplacian = linalg::prod(gLaplacian, gSignalVector);
        VCLMatrix gCombLaplacianMat(gCombLaplacian.handle().opencl_handle(), 1, signalSize);
        gResult += linalg::prod(currentCoeffs, gCombLaplacianMat);

        // 3rd coeff and above
        for( size_t i = 2; i < filterOrder; ++i ) {
            // Get current coefficients for filter order
            range coeffCol(i, i+1); // 1 column only
            matrix_range<VCLMatrix> currentCoeffs = project(gCoeff, coeffRow, coeffCol);

            gCombLaplacian = linalg::prod(gLaplacian, gCombLaplacian);
            // Promote to matrix
            VCLMatrix gCombLaplacianMat(gCombLaplacian.handle().opencl_handle(), 1, signalSize);
            gResult += linalg::prod(currentCoeffs, gCombLaplacianMat);
        }
        copy(gResult, result);
        return true;
    }

    template<typename ScalarType>
    static VCMatrixPtr<ScalarType> copyLaplacian2GPU( const SCMatrix<ScalarType>& laplacian )
    {
        VCMatrixPtr<ScalarType> gLaplacian(new VCMatrix<ScalarType>(laplacian.size(), laplacian.size()));
        viennacl::copy(laplacian, *gLaplacian);
        return gLaplacian;
    }

    template<typename ScalarType>
    static ScalarType getLargestEigenValue( const VCMatrixPtr<ScalarType>& gLaplacianPtr, double precision = 1e-8, int maxIter = 400 );

    template<typename ScalarType>
    static bool loadLaplacian( const std::string& filename, SCMatrix<ScalarType>& A )
    {
        LOG(logINFO) << "Loading graph: " << filename;
        if( io::readMatrixMarketFile(A, filename) <= 1 ) {
            return false;
        }
        return true;
    }

};

// Template specialization
template<>
double Engine::getLargestEigenValue<double>( const VCMatrixPtr<double>& gLaplacianPtr, double precision, int /*maxIter*/ )
{
    viennacl::linalg::lanczos_tag ltag(precision, 1, 0, 20);
    double lmax = viennacl::linalg::eig((*gLaplacianPtr), ltag).at(0) * kERROR_FACTOR; // first
    return lmax;
}

template<>
float Engine::getLargestEigenValue<float>( const VCMatrixPtr<float>& gLaplacianPtr, double precision, int maxIter )
{
    viennacl::linalg::power_iter_tag ptag(precision, maxIter);
    float lmax = viennacl::linalg::eig(*gLaplacianPtr, ptag) * kERROR_FACTOR;
    return lmax;
}

} // end namespace core
} // end namespace ght

#endif // GHT_ENGINE_H
