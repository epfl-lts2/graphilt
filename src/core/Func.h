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

#ifndef GHT_FUNC_H
#define GHT_FUNC_H

#include <functional>
#include <cmath>
#include <boost/lexical_cast.hpp>

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
    virtual std::string print() const = 0;

protected:
    FuncPtr m_f;
};
template <typename ScalarType>
std::ostream& operator<<( std::ostream& os, const Func<ScalarType>& f )
{
    return os << f.print();
}

template <typename ScalarType>
ScalarType Func<ScalarType>::op( ScalarType x ) { return x; }

template <typename ScalarType>
class ExpFunc: public Func<ScalarType>
{
public:
    typedef std::shared_ptr<ExpFunc<ScalarType> > Ptr;
    ExpFunc(): Func<ScalarType>() {}
    explicit ExpFunc( typename Func<ScalarType>::FuncPtr f ): Func<ScalarType>(f) {}
    virtual ScalarType op( ScalarType x ) override
    {
        return static_cast<ScalarType>(std::exp(x));
    }
    virtual std::string print() const override
    {
        if( Func<ScalarType>::m_f )
            return "exp(" + Func<ScalarType>::m_f->print() + ")";
        else
            return "exp(x)";
    }
};

template <typename ScalarType>
class XExpMinusFunc: public Func<ScalarType>
{
public:
    typedef std::shared_ptr<XExpMinusFunc<ScalarType>> Ptr;
    XExpMinusFunc(): Func<ScalarType>() {}
    explicit XExpMinusFunc( typename Func<ScalarType>::FuncPtr f ): Func<ScalarType>(f) {}
    virtual ScalarType op( ScalarType x ) override
    {
        return static_cast<ScalarType>(x * std::exp(-x));
    }
    virtual std::string print() const override
    {
        if( Func<ScalarType>::m_f )
            return Func<ScalarType>::m_f->print() + " * exp(-" + Func<ScalarType>::m_f->print() + ")";
        else
            return "x * exp(x)";
    }
};

template <typename ScalarType>
class PowFunc: public Func<ScalarType>
{
    ScalarType m_fac;
public:
    typedef std::shared_ptr<PowFunc<ScalarType>> Ptr;
    PowFunc(): Func<ScalarType>() {}
    explicit PowFunc( typename Func<ScalarType>::FuncPtr f, ScalarType fac )
        : Func<ScalarType>(f)
        , m_fac(fac)
    {}
    virtual ScalarType op( ScalarType x ) override
    {
        return static_cast<ScalarType>(std::pow(x, m_fac));
    }
    virtual std::string print() const override
    {
        if( Func<ScalarType>::m_f )
            return "(" + Func<ScalarType>::m_f->print() + ")^" + boost::lexical_cast<std::string>(m_fac);
        else
            return "x^" + boost::lexical_cast<std::string>(m_fac);
    }
};

template <typename ScalarType>
class ScaleFunc: public Func<ScalarType>
{
    ScalarType m_fac;
public:
    typedef std::shared_ptr<ScaleFunc<ScalarType>> Ptr;
    ScaleFunc(): Func<ScalarType>() {}
    explicit ScaleFunc( typename Func<ScalarType>::FuncPtr f, ScalarType fac )
        : Func<ScalarType>(f)
        , m_fac(fac)
    {}
    explicit ScaleFunc( ScalarType fac )
        :  m_fac(fac)
    {}

    virtual ScalarType op( ScalarType x ) override
    {
        return static_cast<ScalarType>(x * m_fac);
    }

    virtual std::string print() const override
    {
        if( Func<ScalarType>::m_f )
            return boost::lexical_cast<std::string>(m_fac) + " * " + Func<ScalarType>::m_f->print();
        else
            return boost::lexical_cast<std::string>(m_fac) + "x";
    }
};

template <typename ScalarType>
class MinusFunc: public Func<ScalarType>
{
public:
    typedef std::shared_ptr<MinusFunc<ScalarType>> Ptr;
    MinusFunc(): Func<ScalarType>() {}
    explicit MinusFunc( typename Func<ScalarType>::FuncPtr f )
        : Func<ScalarType>(f)
    {}
    virtual ScalarType op( ScalarType x ) override
    {
        return static_cast<ScalarType>(-x);
    }
    virtual std::string print() const override
    {
        if( Func<ScalarType>::m_f )
            return "-" + Func<ScalarType>::m_f->print();
        else
            return "-x";
    }
};

template <typename ScalarType>
class NoOpFunc: public Func<ScalarType>
{
public:
    typedef std::shared_ptr<NoOpFunc<ScalarType>> Ptr;
    NoOpFunc(): Func<ScalarType>() {}
    virtual ScalarType op( ScalarType x ) override
    {
        return static_cast<ScalarType>(x);
    }
    virtual std::string print() const override
    {
        if( Func<ScalarType>::m_f )
            return Func<ScalarType>::m_f->print();
        else
            return "x";
    }
};

template <typename ScalarType>
class SumFunc: public Func<ScalarType>
{
    int m_low;
    int m_high;
public:
    typedef std::shared_ptr<SumFunc<ScalarType> > Ptr;
    SumFunc()
        : Func<ScalarType>()
        , m_low(0)
        , m_high(0)
    {}

    explicit SumFunc( typename Func<ScalarType>::FuncPtr f, int low, int high )
        : Func<ScalarType>(f)
        , m_low(low)
        , m_high(high)
    {}

    virtual ScalarType op( ScalarType x ) override
    {
        ScalarType res = 0.0;
        for( int i = m_low; i < m_high; ++i ) {
            res += x;
        }
        return static_cast<ScalarType>(std::exp(x));
    }
    virtual std::string print() const override
    {
        if( Func<ScalarType>::m_f )
            return "sum(" + Func<ScalarType>::m_f->print() + ")";
        else
            return "sum(x)";
    }
};






} // end namespace core
} // end namespace ght

#endif // GHT_FUNC_H
