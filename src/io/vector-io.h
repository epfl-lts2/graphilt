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

#ifndef GHT_VECTORIO_H
#define GHT_VECTORIO_H

#include <string>
#include <iostream>
#include <fstream>

#include <viennacl/meta/result_of.hpp>

namespace ght {
namespace io {

template <typename VectorType>
void resize_vector( VectorType & vec, unsigned int size )
{
  vec.resize(size);
}

template <typename VectorType>
bool readVectorFromFile( const std::string& filename, VectorType& vec )
{
  typedef typename viennacl::result_of::value_type<VectorType>::type ScalarType;
  std::ifstream file(filename.c_str());
  if( !file ) return false;

  unsigned int size;
  file >> size;
  
  resize_vector(vec, size);

  for( unsigned int i = 0; i < size; ++i ) {
    ScalarType element;
    file >> element;
    vec[i] = element;
  }

  return true;
}

} // end namespace io
} // end namespace ght

#endif // GHT_VECTORIO_H
