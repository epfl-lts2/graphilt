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

#include <iostream>
//
// ublas includes
//
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/lu.hpp>

// Must be set if you want to use ViennaCL algorithms on ublas objects
#define VIENNACL_WITH_UBLAS 1

#include <viennacl/compressed_matrix.hpp>

#include "core/Engine.h"
#include "io/matrix-io.h"
#include "io/vector-io.h"
#include "util/benchmark.h"
#include "util/log.h"

using namespace ght;
using namespace ght::util;
using namespace ght::core;

using namespace boost::numeric;

namespace  {

template<typename ScalarType>
std::vector<std::vector<ScalarType> > genCoeff( int scales, int numCoeff )
{
    std::vector<std::vector<ScalarType> > coeff(scales);
    for( int i = 0; i < scales; ++i ) {
        std::vector<ScalarType> scale(numCoeff);
        for( int j = 0; j < numCoeff; ++j ) {
            scale[j] = j+2;
        }
        coeff[i]= scale;
    }
    return coeff;
}

template<typename ScalarType>
std::vector<ScalarType> genSignal( int size )
{
    std::vector<ScalarType> signal(size);
    for( int i = 0; i < size; ++i ) {
        signal[i] = i+1;
    }
    return signal;
}

template<typename VectorType>
VectorType genSignal2( int size )
{
    VectorType signal(size);
    for( int i = 0; i < size; ++i ) {
        signal(i) = i+1;
    }
    return signal;
}

template<typename ScalarType>
bool loadGraph( const std::string& filename, std::vector<std::map<uint32_t, ScalarType> >& A )
{
    LOG(logINFO) << "Loading graph: " << filename;
    if( io::readMatrixMarketFile(A, filename) <= 1 ) {
        return false;
    }
    return true;
}

} // end namespace anonymous

static const int kNBSCALES = 5;
static const int kFILTERORDER = 8;
//static const std::string kGRAPHPATH = "../resources/as-skitter.txt.mtx";
//static const std::string kGRAPHPATH = "../resources/com-lj.ungraph.txt.mtx";
//static const std::string kGRAPHPATH = "../resources/randomregular-10000-50.mtx";
static const std::string kGRAPHPATH = "../resources/randomregular-2000-30.mtx";
//static const std::string kGRAPHPATH = "../resources/comet-10000-50.mtx";
//static const std::string kGRAPHPATH = "../resources/minnesota.mtx";

double gCPUTime = 0;
double gGPUNAIVE = 0;
typedef double ScalarType;
std::vector<std::map<uint32_t, ScalarType> > kA;
bool ok = loadGraph(kGRAPHPATH, kA);

TEST( EngineTest, SmallGraph )
{
//    typedef float       ScalarType;
//    ublas::compressed_matrix<ScalarType> A;
//    if( io::readMatrixMarketFile(A, "../resources/testdata/mat65k.mtx") <= 1 ) {
//      std::cout << "Error reading Matrix file" << std::endl;
//      return;
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
//    if( io::readMatrixMarketFile(A, filename) <= 1 ) {
//        ASSERT_TRUE(false);
//        return;
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
    std::vector<std::vector<ScalarType> > coeff = genCoeff<ScalarType>(3,3);

    // Final
    std::vector<std::vector<ScalarType> > result;
    bool ok = Engine::runNaiveGPU(ublas_matrix, rhs, coeff, result);
    EXPECT_EQ(ok, true);
}

TEST( EngineTest, CPU )
{
    Timer timer;
    double exec_time;
    ASSERT_TRUE(ok);

    ublas::compressed_matrix<ScalarType> A(kA.size(), kA.size());
    for( size_t i = 0; i < kA.size(); ++i ) {
        for( auto it = kA.at(i).begin(); it != kA.at(i).end(); ++it ) {
            A(i, it->first) = it->second;
        }
    }

    auto signal = genSignal2<ublas::vector<ScalarType> >(A.size1());
    auto coeff = genCoeff<ScalarType>(kNBSCALES, kFILTERORDER);

    std::vector<std::vector<ScalarType> > result;
    timer.start();
    bool ok = Engine::runNaiveCPU(A, signal, coeff, result);
    exec_time = timer.get();
    gCPUTime = exec_time;
    LOG(logINFO) << " - CPU Execution time: " << exec_time;
    EXPECT_EQ(ok, true);
}

TEST( EngineTest, NAIVEGPU )
{
    Timer timer;
    double exec_time;

    ASSERT_TRUE(ok);
    auto signal = genSignal<ScalarType>(kA.size());
    auto coeff = genCoeff<ScalarType>(kNBSCALES, kFILTERORDER);

    std::vector<std::vector<ScalarType> > result;
    timer.start();
    bool ok = Engine::runNaiveGPU(kA, signal, coeff, result);
    exec_time = timer.get();
    gGPUNAIVE = exec_time;
    LOG(logINFO) << " - GPU Execution time: " << exec_time;
    LOG(logINFO) << " - GPU speedup x" << gCPUTime / exec_time;
    EXPECT_EQ(ok, true);
}

TEST( EngineTest, GPU2 )
{
    Timer timer;
    double exec_time;

    ASSERT_TRUE(ok);
    auto signal = genSignal<ScalarType>(kA.size());
    auto coeff = genCoeff<ScalarType>(kNBSCALES, kFILTERORDER);

    std::vector<std::vector<ScalarType> > result;
    timer.start();
    bool ok = Engine::runGPU(kA, signal, coeff, result);
    exec_time = timer.get();
    LOG(logINFO) << " - GPU Execution time: " << exec_time;
    LOG(logINFO) << " - GPU speedup x" << gCPUTime / exec_time;
    EXPECT_EQ(ok, true);
}









