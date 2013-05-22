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

#ifndef GHT_MATHS_H
#define GHT_MATHS_H

#include <iostream>
#include <vector>
#include <cmath>

namespace ght {
namespace util {

/**
 * @brief linspace(a, b, N) generates N points between a and b.
 * @param a first bound
 * @param b lowest bound
 * @param N number of points
 * @return vector of equally spaced points
 */
template <typename ScalarType>
std::vector<ScalarType> linspace( ScalarType a, ScalarType b, int N )
{
    std::vector<ScalarType> array;
    if( N <= 0 ) {
        return array;
    }

    array.reserve(N);
    if( N == 1 ) {
        array.push_back(b);
    }
    else if( N == 2 ) {
        array.push_back(a);
        array.push_back(b);

    }
    else {
        double step = std::fabs(b-a) / (N-1);
        if( a <= b ) {
            while( a <= b ) {
                array.push_back(a);
                a += step;
            }
            if( a < b )
                array.push_back(b);
        }
        else {
            while( a >= b ) {
                array.push_back(a);
                a -= step;
            }
            if( a > b )
                array.push_back(b);
        }
    }
    return array;
}

template <typename ScalarType>
std::vector<ScalarType> exp( const std::vector<ScalarType>& in )
{
    std::vector<ScalarType> res;
    res.reserve(in.size());
    for( auto& v: in ) {
        res.push_back(std::exp(v));
    }
    return res;
}

} // end namespace util
} // end namespace ght

#endif // GHT_MATHS_H
