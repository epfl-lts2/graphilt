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

#ifndef GHT_MATRIXIO_H
#define GHT_MATRIXIO_H

#include <viennacl/io/matrix_market.hpp>

namespace ght {
namespace io {

/**
* @brief Reads a sparse matrix from a file (MatrixMarket format), wrapper around viennacl
* @param mat The matrix that is to be read (ublas-types and std::vector< std::map <unsigned int, ScalarType> > are supported)
* @param file The filename
* @param index_base The index base, typically 1
* @tparam MatrixType A generic matrix type. Type requirements: size1() returns number of rows, size2() returns number columns, operator() writes array entries, resize() allows resizing the matrix.
* @return Returns nonzero if file is read correctly
*/
template<typename ScalarType>
long readMatrixMarketFile( std::vector< std::map<unsigned int, ScalarType> >& mat, const char* file,
                             long index_base = 1 )
{
    return viennacl::io::read_matrix_market_file(mat, file, index_base);
}

template<typename ScalarType>
long readMatrixMarketFile( std::vector< std::map<unsigned int, ScalarType> >& mat, const std::string& file,
                             long index_base = 1 )
{
    return viennacl::io::read_matrix_market_file(mat, file, index_base);
}


template<typename MatrixType>
long readMatrixMarketFile( MatrixType& mat, const std::string& file,
                             long index_base = 1 )
{
    return viennacl::io::read_matrix_market_file(mat, file, index_base);
}

template<typename MatrixType>
long readMatrixMarketFile( MatrixType& mat, const char* file,
                             long index_base = 1 )
{
    return viennacl::io::read_matrix_market_file(mat, file, index_base);
}

/**
 *
 * @brief Writes a sparse matrix to a file (MatrixMarket format)
 *
 * @param mat The matrix that is to be read (ublas-types and std::vector< std::map <unsigned int, ScalarType> > are supported)
 * @param file The filename
 * @param index_base The index base, typically 1
 * @tparam MatrixType A generic matrix type. Type requirements: size1() returns number of rows, size2() returns number columns, operator() writes array entries, resize() allows resizing the matrix.
 * @return Returns nonzero if file is read correctly
**/
template<typename MatrixType>
void writeMatrixMarketFile( const MatrixType& mat, const std::string& file, long index_base = 1 )
{
    viennacl::io::write_matrix_market_file(mat, file, index_base);
}

template <typename MatrixType>
void writeMatrixMarketFile( const MatrixType& mat, const char* file, long index_base = 1 )
{
    viennacl::io::write_matrix_market_file(mat, file, index_base);
}

template<typename ScalarType>
void writeMatrixMarketFile( const std::vector<std::map<unsigned int, ScalarType> >& mat, const std::string& file,
                              long index_base = 1 )
{
    viennacl::io::write_matrix_market_file(mat, file, index_base);
}

template<typename ScalarType>
void writeMatrixMarketFile( const std::vector<std::map<unsigned int, ScalarType> >& mat, const char* file,
                              long index_base = 1 )
{
    viennacl::io::write_matrix_market_file(mat, file, index_base);
}

} // end namespace io
} // end namespace ght

#endif // GHT_MATRIXIO_H
