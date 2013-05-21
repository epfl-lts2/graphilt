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

#include <vector>
#include <functional>
#include "util/maths.h"

namespace ght {
namespace core {

template <typename ScalarType>
class Func {
public:
    typedef std::shared_ptr<Func<ScalarType>> FuncPtr;
    explicit Func() {}
    virtual ~Func() {}

    Func( FuncPtr f ): m_f(f) {}
    ScalarType apply( ScalarType x )
    {
        if( m_f ) {
            x = m_f->apply(x);
        }
        return op(x);
    }
    virtual ScalarType op( ScalarType x ) = 0;

protected:
    FuncPtr m_f;
};

template <typename ScalarType>
ScalarType Func<ScalarType>::op( ScalarType x ) { return x; }

template <typename ScalarType>
class ExpFunc: public Func<ScalarType>
{
public:
    ExpFunc(): Func<ScalarType>() {}
    explicit ExpFunc( typename Func<ScalarType>::FuncPtr f ): Func<ScalarType>(f) {}
    virtual ScalarType op( ScalarType x ) override
    {
        return static_cast<ScalarType>(std::exp(x));
    }
};

template <typename ScalarType>
class Filter {
    typedef typename Func<ScalarType>::FuncPtr FuncPtr;
    typedef std::vector<FuncPtr> FuncVec;

public:
    typedef typename FuncVec::iterator iterator;
    typedef const typename FuncVec::const_iterator const_iterator;
    Filter() {}
    ~Filter() {}
    FuncPtr& operator[]( size_t i ) { return m_data[i]; }
    const FuncPtr& operator[]( size_t i ) const { return m_data[i]; }
    void push_back( FuncPtr&& val ) { m_data.push_back(val); }
    void push_back( const FuncPtr& f ) { m_data.push_back(f); }
    iterator begin() { return m_data.begin(); }
    const_iterator begin() const { return m_data.begin(); }
    iterator end() { return m_data.end(); }
    const_iterator end() const { return m_data.end(); }

private:
    FuncVec m_data;
};

class FilterFactory {

public:
    template <typename ScalarType>
    static Filter<ScalarType> createFilter()
    {
        typedef std::shared_ptr<ExpFunc<ScalarType>> ExpPtr;
        Filter<ScalarType> filt;
        ExpPtr expFunc(new ExpFunc<ScalarType>());
        ExpPtr expFunc2(new ExpFunc<ScalarType>(expFunc));
        filt.push_back(expFunc2);

        return filt;
    }


    /**
     * @brief Compute a set of wavelet scales adapted to spectrum bounds. Scales logarithmicaly spaced between minimum and maximum
     * "effective" scales : i.e. scales below minumum or above maximum
     * will yield the same shape wavelet (due to homogoneity of sgwt kernel :
     * currently assuming sgwt kernel g given as abspline with t1=1, t2=2)
     * @param lmin minimum nonzero eigenvalue of Laplacian
     * @param lmax maximum eigenvalue of Laplacian
     * @param N number of wavelet scales
     * @return a (possibly good) set of wavelet scales given minimum nonzero and
     maximum eigenvalues of laplacian.
     */
    template <typename ScalarType>
    static std::vector<ScalarType> waveletScales( ScalarType lmin, ScalarType lmax, int N )
    {
        std::vector<ScalarType> res;
        if( lmin == 0 || lmax == 0 || N <= 0 )
            return res;

        ScalarType t1 = 1.0;
        ScalarType t2 = 2.0;

        ScalarType smin = t1 / lmax;
        ScalarType smax = t2 / lmin;
        // scales should be decreasing ... higher j should give larger s
        res = util::exp(util::linspace(std::log(smax), std::log(smin), N));
        return res;
    }
};

} // end namespace core
} // end namespace ght

#endif // GHT_FILTERFACTORY_H
