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

#ifndef NDEBUG
 #define NDEBUG
#endif

#include <iostream>

#include "core/Engine.h"
#include "core/FilterFactory.h"
#include "util/benchmark.h"
#include "util/log.h"
#include "io/vector-io.h"

using namespace ght;
using namespace ght::util;
using namespace ght::core;

namespace {

template<typename ScalarType>
std::vector<ScalarType> genSignal( int size, int zone = 0 )
{
    std::vector<ScalarType> signal(size);
    for( int i = 0; i < size; ++i ) {
        if( zone == 0 ) {
            signal[i] = i+1 / std::sqrt(2);
        }
        else { // partion signal in zones
//            signal[i] = (i+1) % zone;
            signal[i] = i % 2;
        }
    }
    return signal;
}

} // end namespace anonymous

static const int kNBSCALES = 2;
static const int kFILTERORDER = 5;
//static const std::string kGRAPHPATH = "../resources/as-skitter.txt.mtx";
//static const std::string kGRAPHPATH = "../resources/com-lj.ungraph.txt.mtx";
//static const std::string kGRAPHPATH = "../resources/randomregular-10000-50.mtx";
//static const std::string kGRAPHPATH = "../resources/randomregular-2000-30.mtx";
//static const std::string kGRAPHPATH = "../resources/comet-10000-50.mtx";
static const std::string kGRAPHPATH = "../resources/minnesota.mtx";
//static const std::string kGRAPHPATH = "../resources/paques.mtx";
//static const std::string kGRAPH_SIG_PATH = "../resources/paques_sig.txt";

int main( int argc, char *argv[] )
{
    typedef float ScalarType;
    LOG(logINFO) << "Welcome to Graphilt demo 1";
    SCMatrix<ScalarType> laplacian;
    bool ok = Engine::loadGraph(kGRAPHPATH, laplacian);
    if( !ok ) {
        LOG(logERROR) << "Error reading laplacian: " << kGRAPHPATH;
    }
    VCMatrixPtr<ScalarType> gLaplacian = Engine::copyLaplacian2GPU(laplacian);

    auto signal = genSignal<ScalarType>(laplacian.size(), 10);
//    std::vector<ScalarType> signal;
//    ok = io::readVectorFromFile(kGRAPH_SIG_PATH, signal);
//    if( !ok ) {
//        LOG(logERROR) << "Cannot read signal from file " << kGRAPH_SIG_PATH;
//    }

    ScalarType lmax = Engine::getLargestEigenValue(gLaplacian);
    std::pair<ScalarType, ScalarType> range(0.0, lmax);
    auto res = FilterFactory::createFilter<ScalarType>(FilterFactory::MEXICAN_HAT, lmax, kNBSCALES);
    auto coeff = FilterFactory::createChebyCoeff(res, kFILTERORDER, 0, range);
    std::vector<std::vector<ScalarType> > result;

    Timer timer;
    double exec_time;
    timer.start();
    ok = Engine::chebyForwardGPU(result, signal,gLaplacian, coeff, range);
    exec_time = timer.get();
    if( ok )
        LOG(logINFO) << " - GPU Execution time: " << exec_time;
    else
        LOG(logERROR) << " Forward transform failed";


//    std::cout << result[0][0] << std::endl;
    for( auto& v: result ) {
        std::cout << "SCALE: \n";
        for( auto val: v ) {
            std::cout << val << " ";
        }
        std::cout << "\n" << std::endl;
    }

    return 0;
}

