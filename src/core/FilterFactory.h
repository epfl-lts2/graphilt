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

#ifndef GHT_FILTERFACTORY_H
#define GHT_FILTERFACTORY_H

namespace ght {
namespace core {

template <typename ScalarType>
class FilterFactory {

public:
    static void createFilter();
};

} // end namespace core
} // end namespace ght

#endif // GHT_FILTERFACTORY_H
