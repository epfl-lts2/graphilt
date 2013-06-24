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
#include <boost/math/constants/constants.hpp>

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
    // Upcast
    double A = a;
    double B = b;
    std::vector<ScalarType> array;
    if( N <= 0 ) {
        return array;
    }

    array.reserve(N);
    if( N == 1 ) {
        array.push_back(B);
    }
    else if( N == 2 ) {
        array.push_back(A);
        array.push_back(B);

    }
    else {
        double step = std::fabs(B-A) / (N-1);
        if( A <= B ) {
            while( A <= B ) {
                array.push_back(A);
                A += step;
            }
            if( A < B )
                array.push_back(B);
        }
        else {
            while( A >= B ) {
                array.push_back(A);
                A -= step;
            }
            if( A > B )
                array.push_back(B);
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

template <typename ScalarType>
void addInPlace( std::vector<ScalarType>& inout, const std::vector<ScalarType>& rhs )
{
    for( size_t i = 0; i < inout.size(); ++i ) {
        inout[i] += rhs.at(i);
    }
}


} // end namespace util
} // end namespace ght

#endif // GHT_MATHS_H
