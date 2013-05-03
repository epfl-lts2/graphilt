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

using namespace ght::io;

TEST( ParseSnapFormat, UndirectedGraph )
{
    std::string filename = "/Volumes/Storage/datasets/graphs/undirected/as-skitter.txt";
    std::vector<std::map<uint32_t, int> > data;
    EXPECT_EQ(parseSnapFormat(filename, data), true);
    writeMatrixMarketFile(data, filename + ".mtx");
}
