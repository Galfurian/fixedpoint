/// @file fixedPoint.hpp
/// @brief Class for fixed-point operations.
/// @date 7/4/19
/// @author Enrico Fraccaroli
/// @copyright
/// Copyright (c) 2019 Enrico Fraccaroli <enrico.fraccaroli@gmail.com>
/// Permission is hereby granted, free of charge, to any person obtaining a
/// copy of this software and associated documentation files (the "Software"),
/// to deal in the Software without restriction, including without limitation
/// the rights to use, copy, modify, merge, publish, distribute, sublicense,
/// and/or sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following conditions:
///     The above copyright notice and this permission notice shall be included
///     in all copies or substantial portions of the Software.
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
/// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
/// DEALINGS IN THE SOFTWARE.

#pragma once

#include "bitvector.hpp"
#include "support.hpp"

#include <fstream>
#include <sstream>

/// @brief The fixed point class.
/// @tparam WHOLE the size of the whole part.
/// @tparam FRACT the size of the fractional part.
template <long unsigned int WHOLE, long unsigned int FRACT>
class FixedPoint {
private:
    /// The original value.
    double original;
    /// Whole part.
    BitVector<WHOLE> whole;
    /// Fractional part.
    BitVector<FRACT> fractional;

public:
    /// @brief Empty constructor.
    FixedPoint()
        : original(),
          whole(),
          fractional()
    {
        // Nothing to do.
    }

    /// @brief Empty constructor.
    explicit FixedPoint(double val)
        : original(val),
          whole(),
          fractional()
    {
        this->_set_from_double(val);
    }

    /// @brief Constructor from string, e.g., "0011","1000" = "3.5".
    FixedPoint(const std::string &_whole, const std::string &_fractional)
        : original(),
          whole(lpad(_whole, WHOLE, '0')),
          fractional(rpad(_fractional, FRACT, '0'))
    {
        // Nothing to do.
    }

    /// @brief Constructor from bit-vectors, e.g., <0011>,<1000> = 3.5.
    FixedPoint(const BitVector<WHOLE + FRACT> &_full)
        : original(),
          whole(),
          fractional()
    {
        split<WHOLE, FRACT>(_full, whole, fractional);
    }

    /// @brief Constructor from bit-vectors, e.g., <0011>,<1000> = 3.5.
    FixedPoint(const BitVector<WHOLE> &_whole, const BitVector<FRACT> &_fractional)
        : original(),
          whole(_whole),
          fractional(_fractional)
    {
        // Nothing to do.
    }

    /// @brief Constructor from another fixed point.
    FixedPoint(const FixedPoint &other)
        : original(other.original),
          whole(other.whole),
          fractional(other.fractional)
    {
        // Nothing to do.
    }

    /// @brief Constructor from another fixed point.
    template <long unsigned int WHOLE2, long unsigned int FRACT2>
    FixedPoint(const FixedPoint<WHOLE2, FRACT2> &other)
        : original(),
          whole(),
          fractional()
    {
        (*this) = other;
    }

    inline double get_original() const
    {
        return original;
    }

    inline const BitVector<WHOLE> &get_whole() const
    {
        return whole;
    }

    inline bool get_whole(long unsigned int position) const
    {
        return whole[position];
    }

    inline const BitVector<FRACT> &get_fractional() const
    {
        return fractional;
    }

    inline bool get_fractional(long unsigned int position) const
    {
        return fractional[position];
    }

    inline void set_whole(const BitVector<WHOLE> &_whole)
    {
        whole = _whole;
    }

    inline void set_whole(long unsigned int _whole)
    {
        whole = _whole;
    }

    inline void set_whole(const std::string &_whole)
    {
        whole = lpad(_whole, WHOLE, '0');
    }

    inline void set_fractional(const BitVector<FRACT> &_fractional)
    {
        fractional = _fractional;
    }

    inline void set_fractional(long unsigned int _fractional)
    {
        fractional = _fractional;
    }

    inline void set_fractional(const std::string &_fractional)
    {
        fractional = rpad(_fractional, FRACT, '0');
    }

    static inline constexpr long unsigned int get_whole_max()
    {
        return BitVector<WHOLE>::ones().reset(WHOLE - 1).to_number();
    }

    std::string to_string() const
    {
        return whole.to_string() + "." + fractional.to_string();
    }

    double to_number() const
    {
        if (!whole.sign())
            return whole.to_number() + fractional_to_float<FRACT>(fractional);

        // Support variables.
        BitVector<WHOLE + FRACT> _recombined;
        BitVector<WHOLE> _whole;
        BitVector<FRACT> _fractional;

        // 1. Recombine whole and fract.
        recombine(whole, fractional, _recombined);

        // 2. Perform 2's complement.
        _recombined.two_complement();

        // 3. Split the complemented value back into whole and fractional.
        split(_recombined, _whole, _fractional);
        return -(_whole.to_number() + fractional_to_float<FRACT>(_fractional));
    }

    // ========================================================================
    // MATHEMATICAL OPERATIONS
    // ========================================================================
    inline FixedPoint<WHOLE, FRACT> &sum(const FixedPoint<WHOLE, FRACT> &rhs, FixedPoint<WHOLE, FRACT> &solution) const
    {
        // Create bitvector for support.
        BitVector<WHOLE + FRACT> op1;
        BitVector<WHOLE + FRACT> op2;

        // Recombine the operators.
        recombine(whole, fractional, op1);
        recombine(rhs.get_whole(), rhs.get_fractional(), op2);

        // Compute subtraction and return result.
        solution = op1 + op2;
        return solution;
    }

    inline FixedPoint<WHOLE, FRACT> &sub(const FixedPoint<WHOLE, FRACT> &rhs, FixedPoint<WHOLE, FRACT> &solution) const
    {
        // Create bitvector for support.
        BitVector<WHOLE + FRACT> op1;
        BitVector<WHOLE + FRACT> op2;

        // Recombine the operators.
        recombine(whole, fractional, op1);
        recombine(rhs.get_whole(), rhs.get_fractional(), op2);

        // Compute subtraction and return result.
        solution = op1 - op2;
        return solution;
    }

    inline FixedPoint<WHOLE, FRACT> &mul(const FixedPoint<WHOLE, FRACT> &rhs, FixedPoint<WHOLE, FRACT> &solution) const
    {
        // Create bitvector for support.
        BitVector<WHOLE * 2 + WHOLE * 2> support;
        BitVector<WHOLE + FRACT> op1;
        BitVector<WHOLE + FRACT> op2;

        // Recombine whole and fractional.
        recombine(whole, fractional, op1);
        recombine(rhs.get_whole(), rhs.get_fractional(), op2);

        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();

        // Change signes.
        if (op1_neg)
            op1.two_complement();
        if (op2_neg)
            op2.two_complement();

        // Compute multiplication with full precision.
        support = op1 * op2;

        // Apply sign.
        if (op1_neg != op2_neg)
            support.two_complement();

        // Change to string.
        const std::string str = support.to_string();

        // Split whole and fractional.
        solution.set_whole(str.substr(WHOLE, WHOLE));
        solution.set_fractional(str.substr(WHOLE * 2, FRACT));
        return solution;
    }

    inline FixedPoint<WHOLE, FRACT> &div(const FixedPoint<WHOLE, FRACT> &rhs, FixedPoint<WHOLE, FRACT> &solution) const
    {
        // Create bitvector for support.
        BitVector<WHOLE + FRACT * 2> support;
        BitVector<WHOLE + FRACT> _op1;
        BitVector<WHOLE + FRACT> _op2;

        // Recombine whole and fractional.
        recombine(whole, fractional, _op1);
        recombine(rhs.get_whole(), rhs.get_fractional(), _op2);

        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();

        // Change signes.
        if (op1_neg)
            _op1.two_complement();
        if (op2_neg)
            _op2.two_complement();

        // Resize both bit-vectors.
        BitVector<WHOLE + FRACT * 2> op1(rpad(_op1.to_string(), WHOLE + FRACT * 2, '0'));
        BitVector<WHOLE + FRACT * 2> op2(_op2.to_string());

        // Perform division.
        support = op1 / op2;

        // Apply sign.
        if (op1_neg != op2_neg)
            support.two_complement();

        // Change to string.
        const std::string str = support.to_string();

        // Split whole and fractional.
        solution.set_whole(str.substr(FRACT, WHOLE));
        solution.set_fractional(str.substr(WHOLE + FRACT, FRACT));
        return solution;
    }

    // ========================================================================
    // ASSIGN
    // ========================================================================
    template <long unsigned int WHOLE2, long unsigned int FRACT2>
    FixedPoint<WHOLE, FRACT> &operator=(const FixedPoint<WHOLE2, FRACT2> &rhs)
    {
        original = rhs.get_original();

        // Copy whole part.
        whole.assign(rhs.get_whole());

        // Copy fractional part.
        fractional.rassign(rhs.get_fractional());
        if (rhs.get_whole().sign() && (WHOLE > WHOLE2))
            for (unsigned int it = WHOLE2; it < WHOLE; ++it)
                whole[it] = 1;
        return (*this);
    }

    FixedPoint<WHOLE, FRACT> &operator=(double rhs)
    {
        original = rhs;
        this->_set_from_double(rhs);
        return (*this);
    }

    FixedPoint<WHOLE, FRACT> &operator=(const BitVector<WHOLE + FRACT> &rhs)
    {
        split<WHOLE, FRACT>(rhs, whole, fractional);
        return *this;
    }

    // ========================================================================
    // SUM
    // ========================================================================
    inline FixedPoint<WHOLE, FRACT> operator+(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        FixedPoint<WHOLE, FRACT> result;
        return this->sum(rhs, result);
    }

    inline FixedPoint<WHOLE, FRACT> operator+(double rhs) const
    {
        FixedPoint<WHOLE, FRACT> result;
        return this->sum(FixedPoint<WHOLE, FRACT>(rhs), result);
    }

    inline FixedPoint<WHOLE, FRACT> operator+=(const FixedPoint<WHOLE, FRACT> &rhs)
    {
        return this->sum(rhs, *this);
    }

    inline FixedPoint<WHOLE, FRACT> operator+=(double rhs) const
    {
        return this->sum(FixedPoint<WHOLE, FRACT>(rhs), *this);
    }

    // ========================================================================
    // SUB
    // ========================================================================
    inline FixedPoint<WHOLE, FRACT> operator-(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        FixedPoint<WHOLE, FRACT> result;
        return this->sub(rhs, result);
    }

    inline FixedPoint<WHOLE, FRACT> operator-(double rhs) const
    {
        FixedPoint<WHOLE, FRACT> result;
        return this->sub(FixedPoint<WHOLE, FRACT>(rhs), result);
    }

    inline FixedPoint<WHOLE, FRACT> &operator-=(const FixedPoint<WHOLE, FRACT> &rhs)
    {
        return this->sub(rhs, *this);
    }

    inline FixedPoint<WHOLE, FRACT> &operator-=(double rhs)
    {
        return this->sub(FixedPoint<WHOLE, FRACT>(rhs), *this);
    }

    // ========================================================================
    // MUL
    // ========================================================================
    inline FixedPoint<WHOLE, FRACT> operator*(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        FixedPoint<WHOLE, FRACT> result;
        return this->mul(rhs, result);
    }

    inline FixedPoint<WHOLE, FRACT> operator*(double rhs) const
    {
        FixedPoint<WHOLE, FRACT> result;
        return this->mul(FixedPoint<WHOLE, FRACT>(rhs), result);
    }

    inline FixedPoint<WHOLE, FRACT> operator*=(const FixedPoint<WHOLE, FRACT> &rhs)
    {
        return this->mul(rhs, *this);
    }

    inline FixedPoint<WHOLE, FRACT> &operator*=(double rhs)
    {
        return this->mul(FixedPoint<WHOLE, FRACT>(rhs), *this);
    }

    // ========================================================================
    // DIV
    // ========================================================================
    inline FixedPoint<WHOLE, FRACT> operator/(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        FixedPoint<WHOLE, FRACT> result;
        return this->div(rhs, result);
    }

    inline FixedPoint<WHOLE, FRACT> operator/(double rhs) const
    {
        FixedPoint<WHOLE, FRACT> result;
        return this->div(FixedPoint<WHOLE, FRACT>(rhs), result);
    }

    inline FixedPoint<WHOLE, FRACT> operator/=(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        return this->div(rhs, *this);
    }

    inline FixedPoint<WHOLE, FRACT> &operator/=(double rhs)
    {
        return this->div(FixedPoint<WHOLE, FRACT>(rhs), *this);
    }

    // ========================================================================
    // EQUALITY OPs
    // ========================================================================
    inline bool operator==(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        return (whole == rhs.whole) && (fractional == rhs.fractional);
    }

    inline bool operator==(double rhs) const
    {
        return (*this) == FixedPoint<WHOLE, FRACT>(rhs);
    }

    inline bool operator!=(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        return (whole != rhs.whole) || (fractional != rhs.fractional);
    }

    inline bool operator!=(double rhs) const
    {
        return (*this) != FixedPoint<WHOLE, FRACT>(rhs);
    }

    inline bool operator<(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();
        if (op1_neg && !op2_neg)
            return true;
        if (!op1_neg && op2_neg)
            return false;
        // Recombine whole and fractional.
        return recombine(whole, fractional) <
               recombine(rhs.whole, rhs.fractional);
    }

    inline bool operator<(double rhs) const
    {
        return (*this) < FixedPoint<WHOLE, FRACT>(rhs);
    }

    inline bool operator<=(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();
        if (op1_neg && !op2_neg)
            return true;
        if (!op1_neg && op2_neg)
            return false;
        // Recombine whole and fractional.
        return recombine(whole, fractional) <=
               recombine(rhs.whole, rhs.fractional);
    }

    inline bool operator<=(double rhs) const
    {
        return (*this) <= FixedPoint<WHOLE, FRACT>(rhs);
    }

    inline bool operator>(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();
        if (!op1_neg && op2_neg)
            return true;
        if (op1_neg && !op2_neg)
            return false;
        // Recombine whole and fractional.
        return recombine(whole, fractional) >
               recombine(rhs.whole, rhs.fractional);
    }

    inline bool operator>(double rhs) const
    {
        return (*this) > FixedPoint<WHOLE, FRACT>(rhs);
    }

    inline bool operator>=(const FixedPoint<WHOLE, FRACT> &rhs) const
    {
        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();
        if (!op1_neg && op2_neg)
            return true;
        if (op1_neg && !op2_neg)
            return false;
        // Recombine whole and fractional.
        return recombine(whole, fractional) >=
               recombine(rhs.whole, rhs.fractional);
    }

    inline bool operator>=(double rhs) const
    {
        return (*this) >= FixedPoint<WHOLE, FRACT>(rhs);
    }

private:
    void _set_from_double(double val)
    {
        bool sign = (val < 0);
        if (sign)
            val *= -1;
        // Get fractional and whole part.
        double fractpart, intpart;
        fractpart = modf(val, &intpart);
        // Set whole part.
        whole = static_cast<long unsigned int>(intpart);
        // Set fractional part.
        fractional = float_to_fractional<FRACT>(fractpart);
        // Apply the sign.
        if (sign) {
            BitVector<WHOLE + FRACT> full = recombine(whole, fractional);
            full.two_complement();
            split(full, whole, fractional);
        }
    }
};

// ========================================================================
// ARITHMETIC OPERATIONS
// ========================================================================
template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint<WHOLE, FRACT> operator+(double lhs, const FixedPoint<WHOLE, FRACT> &rhs)
{
    return FixedPoint<WHOLE, FRACT>(lhs) + rhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint<WHOLE, FRACT> operator-(double lhs, const FixedPoint<WHOLE, FRACT> &rhs)
{
    return FixedPoint<WHOLE, FRACT>(lhs) - rhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint<WHOLE, FRACT> operator*(double lhs, const FixedPoint<WHOLE, FRACT> &rhs)
{
    return FixedPoint<WHOLE, FRACT>(lhs) * rhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint<WHOLE, FRACT> operator/(double lhs, const FixedPoint<WHOLE, FRACT> &rhs)
{
    static_assert(FRACT != 0, "You are trying to divide with a fractional precision of zero!");
    return FixedPoint<WHOLE, FRACT>(lhs) / rhs;
}

// ========================================================================
// OPERATORS FOR OSTREAM
// ========================================================================
template <long unsigned int WHOLE, long unsigned int FRACT>
std::ostream &operator<<(std::ostream &lhs, const FixedPoint<WHOLE, FRACT> &rhs)
{
    lhs << rhs.to_number();
    return lhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
std::stringstream &operator<<(std::stringstream &lhs, const FixedPoint<WHOLE, FRACT> &rhs)
{
    lhs << rhs.to_number();
    return lhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
std::ifstream &operator>>(std::ifstream &lhs, FixedPoint<WHOLE, FRACT> &rhs)
{
    std::string s;
    lhs >> s;
    size_t size;
    rhs = std::stof(s, &size);
    return lhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
std::stringstream &operator>>(std::stringstream &lhs, FixedPoint<WHOLE, FRACT> &rhs)
{
    std::string s;
    lhs >> s;
    size_t size;
    rhs = std::stof(s, &size);
    return lhs;
}

// ========================================================================
// FUNCTIONS
// ========================================================================
template <long unsigned int WHOLE, long unsigned int FRACT>
FixedPoint<WHOLE, FRACT> round(const FixedPoint<WHOLE, FRACT> &fp)
{
    FixedPoint<WHOLE, FRACT> rounded = fp;
    if ((FRACT >= 1) && rounded.get_fractional().at(FRACT - 1))
        rounded = rounded + 1;
    rounded.set_fractional(0);
    return rounded;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
FixedPoint<WHOLE, FRACT> abs(const FixedPoint<WHOLE, FRACT> &fp)
{
    if (!fp.get_whole().sign())
        return fp;
    return recombine(fp.get_whole(), fp.get_fractional()).two_complement();
}

template <long unsigned int WHOLE, long unsigned int FRACT>
FixedPoint<WHOLE, FRACT> sqrt(
    const FixedPoint<WHOLE, FRACT> &value,
    long max_iterations            = 20,
    FixedPoint<WHOLE, FRACT> error = FixedPoint<WHOLE, FRACT>(0.001))
{
    FixedPoint<WHOLE, FRACT> result(1.0);
    int iteration = 0;
    while ((iteration++ < max_iterations) &&
           (abs((result * result) - value) >= error))
        result = 0.5f * (result + (value / result));
    return result;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
FixedPoint<WHOLE, FRACT> pow(const FixedPoint<WHOLE, FRACT> &value, unsigned exponent)
{
    if (exponent == 0)
        return FixedPoint<WHOLE, FRACT>(1);
    FixedPoint<WHOLE, FRACT> result = value;
    for (unsigned it = 1; it < exponent; ++it)
        result *= value;
    return result;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint<WHOLE, FRACT> exp(FixedPoint<WHOLE, FRACT> x, unsigned prec = 12)
{
    x = 1.0 + (x / std::pow(2, prec));
    for (unsigned i = 0; i < prec; ++i)
        x *= x;
    return x;
}