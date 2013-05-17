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
#include "util/benchmark.h"
#include "util/log.h"

using namespace ght;
using namespace ght::util;
using namespace ght::core;

static const int kNBSCALES = 5;
static const int kFILTERORDER = 8;
//static const std::string kGRAPHPATH = "../resources/as-skitter.txt.mtx";
//static const std::string kGRAPHPATH = "../resources/com-lj.ungraph.txt.mtx";
//static const std::string kGRAPHPATH = "../resources/randomregular-10000-50.mtx";
//static const std::string kGRAPHPATH = "../resources/randomregular-2000-30.mtx";
//static const std::string kGRAPHPATH = "../resources/comet-10000-50.mtx";
static const std::string kGRAPHPATH = "../resources/paques.mtx";

int main( int argc, char *argv[] )
{
    typedef float ScalarType;
    LOG(logINFO) << "Welcome to Graphilt demo 1";
    SCMatrix<ScalarType> laplacian;
    bool ok = Engine::loadLaplacian(kGRAPHPATH, laplacian);
    if( !ok ) {
        LOG(logERROR) << "Error reading laplacian: " << kGRAPHPATH;
    }
    VCMatrixPtr<ScalarType> gLaplacian = Engine::copyLaplacian2GPU(laplacian);
    LOG(logINFO) << "Largest eigenvalue: " << Engine::getLargestEigenValue(gLaplacian);

    return 0;
}

