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

using namespace ght;
using namespace ght::util;
using namespace ght::core;

namespace {

template<typename ScalarType>
std::vector<ScalarType> genSignal( int size )
{
    std::vector<ScalarType> signal(size);
    for( int i = 0; i < size; ++i ) {
        signal[i] = i+1 / std::sqrt(2);
    }
    return signal;
}

} // end namespace anonymous

static const int kNBSCALES = 1;
static const int kFILTERORDER = 8;
//static const std::string kGRAPHPATH = "../resources/as-skitter.txt.mtx";
//static const std::string kGRAPHPATH = "../resources/com-lj.ungraph.txt.mtx";
//static const std::string kGRAPHPATH = "../resources/randomregular-10000-50.mtx";
//static const std::string kGRAPHPATH = "../resources/randomregular-2000-30.mtx";
static const std::string kGRAPHPATH = "../resources/comet-10000-50.mtx";
//static const std::string kGRAPHPATH = "../resources/paques.mtx";

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

    auto signal = genSignal<ScalarType>(laplacian.size());
    ScalarType lmax = Engine::getLargestEigenValue(gLaplacian);
    std::pair<ScalarType, ScalarType> range(0.0, lmax);
    auto res = FilterFactory::createFilter<ScalarType>(FilterFactory::MEXICAN_HAT, lmax, kNBSCALES);
    auto coeff = FilterFactory::createChebyCoeff(res, kFILTERORDER, 0, range);
    std::vector<std::vector<ScalarType> > result;
    Timer timer;
    double exec_time;
    timer.start();
    ok = Engine::chebyForwardGPU(result, signal, gLaplacian, coeff, range);
    exec_time = timer.get();
    if( ok )
        LOG(logINFO) << " - GPU Execution time: " << exec_time;
    else
        LOG(logERROR) << " Forward transform failed";

    return 0;
}

