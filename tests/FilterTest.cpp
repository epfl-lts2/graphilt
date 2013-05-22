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

#include "core/FilterFactory.h"
#include "util/maths.h"

using namespace ght;
using namespace ght::core;

TEST( FilterTest, waveletScales )
{
    auto res = FilterFactory::waveletScales(1.0, 20.0, 5);
    EXPECT_EQ(res[0], 2.0);
    EXPECT_EQ(float(res[4]), float(0.05) );
    for( auto v: res ) {
        std::cout << v << std::endl;
    }
}

TEST( FilterTest,  createFilt )
{
    typedef float ScalarType;
    auto res = FilterFactory::createFilter<ScalarType>(FilterFactory::MEXICAN_HAT, 1.0, 5, 20);
//    for( auto v: res ) {
//        float q = v->apply(2.0);
//        std::cout << q << std::endl;
//    }
//    Filter<ScalarType> filt;
//    typename Func<ScalarType>::FuncPtr scale( new ScaleFunc<ScalarType>(2.0) );
//    typename Func<ScalarType>::FuncPtr gb( new XExpMinusFunc<ScalarType>(scale) );
//    filt.push_back(gb);
//    filt.push_back(gb);
//    filt.push_back(gb);

    std::cout << res[0]->apply(0.2) << std::endl;
    std::cout << res << std::endl;

}


