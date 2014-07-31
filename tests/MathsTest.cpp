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
#include <iostream>

#include "util/maths.h"

using namespace ght;
using namespace ght::util;

TEST( MathsTest, Linespace )
{
    auto res = linspace(-5.56, 8.0, 3);
    for( auto v: res ) {
        std::cout << v << std::endl;
    }
}
