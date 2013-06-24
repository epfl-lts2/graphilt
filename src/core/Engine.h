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
#include <viennacl/linalg/cg.hpp>

#include "util/log.h"
#include "util/benchmark.h"
#include "io/matrix-io.h"
#include "util/benchmark.h"
#include "util/maths.h"

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
    static bool chebyNaiveForwardCPU( const UMatrix<ScalarType>& laplacian,
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
    static bool chebyNaiveForwardGPU( const MatrixType& laplacian, const VectorType& signal,
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

    template<typename ScalarType>
    static bool chebyForwardGPU( SVecvec<ScalarType>& result,
                                 const SVec<ScalarType>& signal,
                                 const VCMatrixPtr<ScalarType>& laplacian,
                                 const SVecvec<ScalarType>& coeff,
                                 const std::pair<ScalarType, ScalarType>& arange )
    {
        using namespace viennacl;
        typedef matrix<ScalarType> VCLMatrix;
        typedef vector<ScalarType> VCLVector;

        size_t signalSize = signal.size();
        size_t nbScales = coeff.size();

        VCLVector gSignalVector(signalSize);
        fast_copy(signal, gSignalVector);

        VCLMatrix gResult = chebyForwardImpl(gSignalVector, laplacian, coeff, arange);

        // Copy back results from GPU
        result.resize(nbScales);
        for( size_t i = 0; i < result.size(); ++i )
            result[i].resize(signalSize);
        copy(gResult, result);
        return true;
    }

    template<typename ScalarType>
    static bool chebyInverseTransform( SVec<ScalarType>& result,
                                       const SVecvec<ScalarType>& input,
                                       const VCMatrixPtr<ScalarType>& laplacian,
                                       const SVecvec<ScalarType>& coeffs,
                                       const std::pair<ScalarType, ScalarType>& arange,
                                       double tolerance = 1e-6,
                                       int maxIter = 200 )
    {

        using namespace viennacl;
        typedef matrix<ScalarType> VCLMatrix;
        typedef vector<ScalarType> VCLVector;
        size_t nbScales = input.size();
        size_t signalSize = input.at(0).size();
        // Copy input in GPU memory
        VCLMatrix in(nbScales, signalSize);
        copy(input, in);

        VCLVector adjoint = adjointImpl(in, laplacian, coeffs, arange);

        SVecvec<ScalarType> padCoeffs = normalizeCoeff(coeffs);
        SVec<ScalarType> mergeCoeffs;
        for( size_t i = 0; i < padCoeffs.size(); ++i ) {
            if( i == 0 )
                mergeCoeffs = chebySquarePolynomial(padCoeffs.at(i));
            else
                util::addInPlace(mergeCoeffs, chebySquarePolynomial(padCoeffs.at(i)));
        }

//        VCLMatrix gWstarW = chebyForwardImpl(adjoint, laplacian, SVecvec<ScalarType>{ mergeCoeffs }, arange);

//        linalg::cg_tag custom_cg(tolerance, maxIter);
//        VCLVector res = viennacl::linalg::solve(gWstarW, adjoint, custom_cg);
//        //print number of iterations taken and estimated error:
//        std::cout << "No. of iters: " << custom_cg.iters() << std::endl;
//        std::cout << "Est. error: " << custom_cg.error() << std::endl;

//        std::cout << res << std::endl;

        return true;
    }

    template<typename ScalarType>
    static bool chebyAdjoint( SVec<ScalarType>& result,
                              const SVecvec<ScalarType>& input,
                              const VCMatrixPtr<ScalarType>& laplacian,
                              const SVecvec<ScalarType>& coeffs,
                              const std::pair<ScalarType, ScalarType>& arange )
    {
        using namespace viennacl;
        typedef matrix<ScalarType> VCLMatrix;
        typedef vector<ScalarType> VCLVector;
        size_t nbScales = input.size();
        size_t signalSize = input.at(0).size();
        // Copy input in GPU memory
        VCLMatrix in(nbScales, signalSize);
        copy(input, in);

        VCLVector res = adjointImpl(in, laplacian, coeffs, arange);
        // Copy back result
        result.resize(res.size());
        copy(res, result);
        return true;
    }

    /**
     * @brief Chebyshev coefficients for square of polynomial
     * @param input Chebyshev coefficients for p(x) = sum c(1+k) T_k(x); 0<=K<=M
     * @return Chebyshev coefficients for p(x)^2 = sum d(1+k) T_k(x); 0<=k<=2*M
     */
    template<typename ScalarType>
    static SVec<ScalarType> chebySquarePolynomial( const SVec<ScalarType>& input )
    {
        SVec<ScalarType> cp = input;
        cp[0] *= 0.5;

        size_t M = input.size();
        SVec<ScalarType> out(2 * M - 1);
        for( size_t i = 0; i < out.size(); ++i )
            out[i] = 0;

        // Calc each coeff
        for( size_t k = 0; k < out.size(); ++k ) {
            if( k == 0 ) {
                out[k] += 0.5 * cp[k] * cp[k];
                for( size_t i = 0; i < M; ++i )
                    out[k] += 0.5 * cp[i] * cp[i];
            }
            else if( k < M ) {
                for( size_t i = 0; i <= k; ++i )
                    out[k] += 0.5 * cp[i] * cp[k-i];
                for( size_t i = 0; i < M - k; ++i )
                    out[k] += 0.5 * cp[i] * cp[i+k];
                for( size_t i = k; i < M; ++i )
                    out[k] += 0.5 * cp[i] * cp[i-k];
            }
            else { // M <= k <= 2 * M
                for( size_t i = k-M+1; i < M; ++i )
                    out[k] += 0.5 * cp.at(i) * cp.at(k-i);
            }
        }
        out[0] *= 2;
        return out;
    }


    /**
     * @brief Normalize coeff to have a square matrix, fill with 0 if coeffs < MaxCoeffs
     * @param in un-normalized coeffs
     * @return normalized coeffs
     */
    template<typename ScalarType>
    static SVecvec<ScalarType> normalizeCoeff( const SVecvec<ScalarType>& in )
    {
        SVecvec<ScalarType> result;
        size_t maxOrder = in.at(0).size();
        bool expand = false;
        for( size_t i = 0; i < in.size(); ++i ) {
            size_t t = in.at(i).size();
            if( t > maxOrder ) {
                expand = true;
                maxOrder = t;
            }
        }
        // Needs expanding
        if( expand ) {
            result = in;
            for( size_t i = 0; i < in.size(); ++i ) {
                size_t t = in.at(i).size();
                // If those coeffs needs expansion
                if( t < maxOrder ) {
                    result[i].resize(maxOrder); // use default ScalarType constructor which is 0;
                }
            }
        }
        else {
            return in;
        }
        return result;
    }

    template<typename ScalarType>
    static VCMatrixPtr<ScalarType> copyLaplacian2GPU( const SCMatrix<ScalarType>& laplacian )
    {
        VCMatrixPtr<ScalarType> gLaplacian(new VCMatrix<ScalarType>(laplacian.size(), laplacian.size()));
        viennacl::copy(laplacian, *gLaplacian);
        return gLaplacian;
    }

    template<typename ScalarType>
    static ScalarType getLargestEigenValue( const VCMatrixPtr<ScalarType>& gLaplacianPtr,
                                            double precision = 1e-8, int maxIter = 400 );

    template<typename ScalarType>
    static bool loadGraph( const std::string& filename, SCMatrix<ScalarType>& A )
    {
        LOG(logINFO) << "Loading graph: " << filename;
        if( io::readMatrixMarketFile(A, filename) <= 1 ) {
            return false;
        }
        return true;
    }



private:

    template<typename ScalarType>
    static viennacl::matrix<ScalarType> chebyForwardImpl( viennacl::vector<ScalarType>& signal,
                                  const VCMatrixPtr<ScalarType>& laplacian,
                                  const SVecvec<ScalarType>& inCoeff,
                                  const std::pair<ScalarType, ScalarType>& arange )
    {
        using namespace viennacl;
        typedef vector<ScalarType> VCLVector;
        typedef matrix<ScalarType> VCLMatrix;
        // Init GPU values
        SVecvec<ScalarType> coeff = normalizeCoeff(inCoeff);
        size_t nbScales = coeff.size();
        size_t signalSize = signal.size();
        size_t filterOrder = coeff.at(0).size(); // coeff have been normalized

        VCLVector gTwfNew(signalSize);
        VCLMatrix gCoeff(nbScales, filterOrder);
        VCLMatrix gResult(nbScales, signalSize);

        // Copy to GPU
        // Signal and Laplacian are already on gpu memory
        VCLVector gSignalVector = signal;
//        fast_copy(signal, gSignalVector);
        VCMatrix<ScalarType> gLaplacian = *laplacian;
        copy(coeff, gCoeff);

//        util::Timer timer;
//        timer.start();

        scalar<ScalarType> a1 = ScalarType((arange.second - arange.first) / 2.0);
        scalar<ScalarType> a2 = ScalarType((arange.second + arange.first) / 2.0);
        scalar<ScalarType> a3 = ScalarType(4.0 / (arange.second + arange.first)); // 2 / a1

        // Starts here
        VCLVector gTwfOld = gSignalVector;
        VCLVector gTwfCur = (linalg::prod(gLaplacian, gSignalVector) - a2*gSignalVector) / a1;

        // Bias coeff
        range coeffRow(0, coeff.size());
        range coeffCol(0, 1); // 1 column only
        matrix_range<VCLMatrix> coeff0 = project(gCoeff, coeffRow, coeffCol);

        // First coeff
        coeffCol = range(1, 2); // 1 column only
        matrix_range<VCLMatrix> coeff1 = project(gCoeff, coeffRow, coeffCol);

        // Init result matrix
        // r{j}=.5*c{j}(1)*Twf_old + c{j}(2)*Twf_cur;
        VCLVector gTemp = gTwfOld * 0.5;
        gResult = linalg::prod(coeff0, VCLMatrix(gTemp.handle().opencl_handle(), 1, signalSize));
        gResult += linalg::prod(coeff1, VCLMatrix(gTwfCur.handle().opencl_handle(), 1, signalSize));

        // 2nd coeff and above
        for( size_t i = 2; i < filterOrder; ++i ) {
            // Get current coefficients for filter order
            range coeffCol(i, i+1); // 1 column only
            matrix_range<VCLMatrix> currentCoeffs = project(gCoeff, coeffRow, coeffCol);
            // Twf_new = (2/a1)*(L*Twf_cur-a2*Twf_cur)-Twf_old;
            gTwfNew = a3 * (linalg::prod(gLaplacian, gTwfCur) - a2*gTwfCur) - gTwfOld;
            // Promote to matrix
            gResult += linalg::prod(currentCoeffs, VCLMatrix(gTwfNew.handle().opencl_handle(), 1, signalSize));
            gTwfOld = gTwfCur;
            gTwfCur = gTwfNew;
        }
//        double res = timer.get();
//        LOG(logINFO) << " - GPU Execution time gpu only: " << res;
        return gResult;
    }

    template<typename ScalarType>
    static viennacl::vector<ScalarType> adjointImpl( viennacl::matrix<ScalarType>& input,
                                  const VCMatrixPtr<ScalarType>& laplacian,
                                  const SVecvec<ScalarType>& coeff,
                                  const std::pair<ScalarType, ScalarType>& arange )
    {
        using namespace viennacl;
        typedef vector<ScalarType> VCLVector;
        typedef matrix<ScalarType> VCLMatrix;

        size_t signalSize = input.size2();
        VCLVector adjoint;

        range coeffCol(0, signalSize);
        for( size_t i = 0; i < input.size1(); ++i ) {
            range coeffRow(i, i+1); // 1 row only
            matrix_range<VCLMatrix> row = project(input, coeffRow, coeffCol);
            // Tmp variable to convert to Vector
            VCLMatrix line = row;
            // Create vector from Matrix(signalsize, 1)
            VCLVector signal(line.handle().opencl_handle(), signalSize);
            VCLMatrix result = chebyForwardImpl(signal, laplacian, SVecvec<ScalarType>{ coeff.at(i) }, arange);
            // Result is a Matrix(signalsize, 1)
            VCLVector resVec(result.handle().opencl_handle(), signalSize);

            if( i == 0 )
                adjoint = resVec;
            else
                adjoint += resVec;
        }
        return adjoint;
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
