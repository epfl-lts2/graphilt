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
#include "Func.h"
#include "util/maths.h"

namespace ght {
namespace core {

template <typename ScalarType>
class Filter {
    typedef typename Func<ScalarType>::FuncPtr FuncPtr;
    typedef std::vector<FuncPtr> FuncVec;

public:
    typedef typename FuncVec::iterator iterator;
    typedef const typename FuncVec::const_iterator const_iterator;
    Filter() {}
    ~Filter() {}
    Filter( const Filter& rhs )
        : m_data(rhs.m_data)
    {}

    void swap( Filter& lhs, Filter& rhs )
    {
        using std::swap;
        swap(lhs.m_data, rhs.m_data);
    }

    Filter& operator =( Filter rhs )
    {
        swap(*this, rhs);
        return *this;
    }

    FuncPtr& operator[]( size_t i ) { return m_data[i]; }
    const FuncPtr& operator[]( size_t i ) const { return m_data[i]; }
    void push_back( FuncPtr&& val ) { m_data.push_back(val); }
    void push_back( const FuncPtr& f ) { m_data.push_back(f); }
    iterator begin() { return m_data.begin(); }
    const_iterator begin() const { return m_data.begin(); }
    iterator end() { return m_data.end(); }
    const_iterator end() const { return m_data.end(); }

    std::string print() const
    {
        std::string os("");
        for( size_t i = 0; i < m_data.size(); ++i ) {
            os += "[" + boost::lexical_cast<std::string>(i) + "]" + " " + m_data.at(i)->print() + " \n";
        }
        return os;
    }


private:
    FuncVec m_data;
};

template <typename ScalarType>
std::ostream& operator<<( std::ostream& os, const Filter<ScalarType>& filt )
{
    return os << filt.print();
}

class FilterFactory {
public:

    enum FilterClass {
        MEXICAN_HAT,
        MEYER,
        ABSPLINE3,
        UNDEFINIED
    };

    template <typename ScalarType>
    static Filter<ScalarType> createFilter( FilterClass type, ScalarType lmax, int Nscales, ScalarType lpFactor )
    {
        Filter<ScalarType> filt;
        switch( type ) {
            case MEXICAN_HAT:
                filt = buildMexicanHat(lmax, Nscales, lpFactor);
                break;
            default:
                break;
        }

        //        case 'meyer'
//            t=(4/(3*lmax)) * 2.^(Nscales-1:-1:0);
//            g{1}= @(x) sgwt_kernel_meyer(t(1)*x,'sf');
//            for j=1:Nscales
//                g{j+1}= @(x) sgwt_kernel_meyer(t(j)*x,'wavelet');
//            end
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

private:
    template <typename ScalarType>
    static Filter<ScalarType> buildMexicanHat( ScalarType lmax, int Nscales, ScalarType lpFactor )
    {
        Filter<ScalarType> filt;

        // Bias term
        // 1.2*exp(-1) * exp((-x/lminfac)^4);
        ScalarType lmin = lmax / lpFactor;
        ScalarType lminFac =  0.4 * lmin;
        typename Func<ScalarType>::FuncPtr scale( new ScaleFunc<ScalarType>(-1.0/lminFac) );
        typename Func<ScalarType>::FuncPtr powFunc(new PowFunc<ScalarType>(scale, 4));
        typename Func<ScalarType>::FuncPtr expFunc(new ExpFunc<ScalarType>(powFunc));
        typename Func<ScalarType>::FuncPtr gb( new ScaleFunc<ScalarType>(expFunc, 1.2*std::exp(-1)) );
        filt.push_back(gb);

        // Scales
        std::vector<ScalarType> t = waveletScales(lmin, lmax, Nscales);

        // scale[i]*x * exp(-scale[i]*x)
        for( int i = 0; i < Nscales; ++i ) {
            typename Func<ScalarType>::FuncPtr scale( new ScaleFunc<ScalarType>(t.at(i)) );
            typename Func<ScalarType>::FuncPtr gb( new XExpMinusFunc<ScalarType>(scale) );
            filt.push_back(gb);
        }
        return filt;
    }
};

} // end namespace core
} // end namespace ght

#endif // GHT_FILTERFACTORY_H
