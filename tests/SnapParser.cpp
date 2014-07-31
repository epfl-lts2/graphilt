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
#include "io/parsers.h"
#include "io/matrix-io.h"
#include "util/log.h"

#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace ght;
using namespace ght::io;
using namespace boost::numeric;

typedef int ScalarType;

TEST( ParseSnapFormat, UndirectedGraphVector )
{
//    std::string filename = "../resources/as-skitter.txt";
    std::string filename = "../resources/com-lj.ungraph.txt";
    auto* data = new std::vector<std::map<uint32_t, ScalarType> >();

    try {
        LOG(logDEBUG) << "Start parsing " << filename;
        bool ok = parseSnapFormat(filename, *data);
        LOG(logDEBUG) << "End parsing";
        EXPECT_EQ(ok, true);
        writeMatrixMarketFile(*data, filename + ".mtx");
    } catch( std::exception& e ) {
        std::cout << e.what() << std::endl;
    }

    delete data;
}

//TEST( ParseSnapFormat, UndirectedGraphMatrix )
//{
////    std::string filename = "/Volumes/Storage/datasets/graphs/undirected/as-skitter.txt";
//    std::string filename = "/Volumes/Storage/datasets/graphs/undirected/com-lj.ungraph.txt";
//    auto* data = new ublas::compressed_matrix<ScalarType>();

//    try {
//        LOG(logDEBUG) << "Start parsing " << filename;
//        bool ok = parseSnapFormat(filename, *data);
//        LOG(logDEBUG) << "End parsing";
//        EXPECT_EQ(ok, true);
////        writeMatrixMarketFile(*data, filename + ".mtx");
//    } catch( std::exception& e ) {
//        std::cout << e.what() << std::endl;
//    }

//    delete data;
//}

