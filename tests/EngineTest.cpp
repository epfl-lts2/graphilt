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

#include <gtest/gtest.h>
#ifndef NDEBUG
 #define NDEBUG
#endif
//
// ublas includes
//
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

// Must be set if you want to use ViennaCL algorithms on ublas objects
#define VIENNACL_WITH_UBLAS 1

#include <viennacl/compressed_matrix.hpp>

#include "core/Engine.h"
#include "io/matrix-io.h"
#include "io/vector-io.h"

using namespace ght;
using namespace ght::core;

using namespace boost::numeric;

TEST( EngineTest, SmallGraph )
{
//    typedef float       ScalarType;
//    ublas::compressed_matrix<ScalarType> A;
//    if( !io::readMatrixMarketFile(A, "../resources/testdata/mat65k.mtx") ) {
//      std::cout << "Error reading Matrix file" << std::endl;
//    }
//    size_t mSize = A.nnz() * sizeof(ScalarType);
//    Engine engine;
//    std::vector<int> signal(10);
//    EXPECT_EQ(engine.checkFitInGPUMem(mSize, signal.size() * sizeof(ScalarType)), true);
}

TEST( EngineTest, TooBigGraph )
{
//    typedef float       ScalarType;
//    ublas::compressed_matrix<ScalarType> A;
//    std::string filename = "/Volumes/Storage/datasets/graphs/undirected/com-lj.ungraph.txt.mtx";
//    if( !io::readMatrixMarketFile(A, filename) ) {
//      std::cout << "Error reading Matrix file" << std::endl;
//    }
//    size_t mSize = A.nnz() * sizeof(ScalarType);
//    Engine engine;
//    std::vector<int> signal(10);
//    EXPECT_EQ(engine.checkFitInGPUMem(mSize , signal.size() * sizeof(ScalarType)), false);
}

TEST( EngineTest, GPU1 )
{
    typedef float       ScalarType;
    std::size_t size = 5;
    ublas::vector<ScalarType> rhs(size, ScalarType(size));
    ublas::compressed_matrix<ScalarType> ublas_matrix(size, size);

    // Laplacian
    ublas_matrix(0,0) =  2.0f; ublas_matrix(0,1) = -1.0f;
    ublas_matrix(1,0) = -1.0f; ublas_matrix(1,1) =  2.0f; ublas_matrix(1,2) = -1.0f;
    ublas_matrix(2,1) = -1.0f; ublas_matrix(2,2) =  2.0f; ublas_matrix(2,3) = -1.0f;
    ublas_matrix(3,2) = -1.0f; ublas_matrix(3,3) =  2.0f; ublas_matrix(3,4) = -1.0f;
    ublas_matrix(4,3) = -1.0f; ublas_matrix(4,4) =  2.0f;

    // Signal
    for( size_t i = 0; i < size; ++i )
        rhs[i] = i;

    // Coefficients
    std::vector<std::vector<ScalarType> > coeff;
    std::vector<ScalarType> coeff0;
    for( size_t i = 0; i < size; ++i )
        coeff0.push_back(i);
    coeff.push_back(coeff0);

    // Final
    std::vector<std::vector<ScalarType> > result;
    bool ok = Engine::runGPU(ublas_matrix, rhs, coeff, result);
}

