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

#include <fstream>

#include "bitvector.hpp"

/// @brief Adds character c to the left of s, until s has a length of n.
inline std::string lpad(std::string const &s, size_t n, char c)
{
    std::string _s = s;
    if (n > _s.length())
        _s.insert(0, std::string(n - _s.length(), c));
    return _s;
}

/// @brief Adds character c to the right of s, until s has a length of n.
inline std::string rpad(std::string const &s, size_t n, char c)
{
    std::string _s = s;
    if (n > _s.length())
        _s.append(std::string(n - _s.length(), c));
    return _s;
}

/// @brief
/// @tparam WHOLE
/// @tparam FRACT
template <long unsigned int WHOLE, long unsigned int FRACT>
class FixedPoint_t {
private:
    /// The original value.
    double original;
    /// Whole part.
    BitVector<WHOLE> whole;
    /// Fractional part.
    BitVector<FRACT> fractional;

    void set_from_double(double val)
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

public:
    /// @brief Empty constructor.
    FixedPoint_t()
        : original(),
          whole(),
          fractional()
    {
        // Nothing to do.
    }

    /// @brief Empty constructor.
    explicit FixedPoint_t(double val)
        : original(val),
          whole(),
          fractional()
    {
        this->set_from_double(val);
    }

    /// @brief Constructor from string, e.g., "0011","1000" = "3.5".
    FixedPoint_t(std::string const &_whole,
                 std::string const &_fractional)
        : original(),
          whole(lpad(_whole, WHOLE, '0')),
          fractional(rpad(_fractional, FRACT, '0'))
    {
        // Nothing to do.
    }

    /// @brief Constructor from bit-vectors, e.g., <0011>,<1000> = 3.5.
    FixedPoint_t(BitVector<WHOLE + FRACT> const &_full)
        : original(),
          whole(),
          fractional()
    {
        split<WHOLE, FRACT>(_full, whole, fractional);
    }

    /// @brief Constructor from bit-vectors, e.g., <0011>,<1000> = 3.5.
    FixedPoint_t(BitVector<WHOLE> const &_whole,
                 BitVector<FRACT> const &_fractional)
        : original(),
          whole(_whole),
          fractional(_fractional)
    {
        // Nothing to do.
    }

    /// @brief Constructor from another fixed point.
    FixedPoint_t(FixedPoint_t const &fixedPoint)
        : original(fixedPoint.original),
          whole(fixedPoint.whole),
          fractional(fixedPoint.fractional)
    {
        // Nothing to do.
    }

    /// @brief Constructor from another fixed point.
    template <long unsigned int WHOLE2, long unsigned int FRACT2>
    FixedPoint_t(FixedPoint_t<WHOLE2, FRACT2> const &fixedPoint)
        : original(),
          whole(),
          fractional()
    {
        (*this) = fixedPoint;
    }

    inline double get_original() const
    {
        return original;
    }

    inline BitVector<WHOLE> const &get_whole() const
    {
        return whole;
    }

    inline bool get_whole(long unsigned int position) const
    {
        return whole[position];
    }

    inline BitVector<FRACT> const &get_fractional() const
    {
        return fractional;
    }

    inline bool get_fractional(long unsigned int position) const
    {
        return fractional[position];
    }

    inline void set_whole(BitVector<WHOLE> const &_whole)
    {
        whole = _whole;
    }

    inline void set_whole(long unsigned int _whole)
    {
        whole = _whole;
    }

    inline void set_fractional(BitVector<FRACT> const &_fractional)
    {
        fractional = _fractional;
    }

    inline void set_fractional(long unsigned int _fractional)
    {
        fractional = _fractional;
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
    // ASSIGN
    // ========================================================================
    template <long unsigned int RHS_WHOLE, long unsigned int RHS_FRACT>
    FixedPoint_t<WHOLE, FRACT> &operator=(
        FixedPoint_t<RHS_WHOLE, RHS_FRACT> const &rhs)
    {
        original = rhs.get_original();
        // Copy whole part.
        whole.assign(rhs.get_whole());
        // Copy fractional part.
        fractional.rassign(rhs.get_fractional());
        if (rhs.get_whole().sign() && (WHOLE > RHS_WHOLE))
            for (unsigned int it = RHS_WHOLE; it < WHOLE; ++it)
                whole[it] = 1;
        return (*this);
    }

    FixedPoint_t<WHOLE, FRACT> &operator=(double rhs)
    {
        original = rhs;
        this->set_from_double(rhs);
        return (*this);
    }

    FixedPoint_t<WHOLE, FRACT> &operator=(BitVector<WHOLE + FRACT> const &rhs)
    {
        split<WHOLE, FRACT>(rhs, whole, fractional);
    }

    // ========================================================================
    // SUM
    // ========================================================================
    inline FixedPoint_t<WHOLE, FRACT> sum(FixedPoint_t<WHOLE, FRACT> const &rhs,
                                          BitVector<WHOLE + FRACT> &op1,
                                          BitVector<WHOLE + FRACT> &op2) const
    {
        // Recombine the operators.
        recombine(whole, fractional, op1);
        recombine(rhs.get_whole(), rhs.get_fractional(), op2);
        // Compute summation and return result.
        return FixedPoint_t<WHOLE, FRACT>(op1 + op2);
    }

    inline FixedPoint_t<WHOLE, FRACT>
    operator+(FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        BitVector<WHOLE + FRACT> op1, op2;
        return this->sum(rhs, op1, op2);
    }

    inline FixedPoint_t<WHOLE, FRACT> operator+(double rhs) const
    {
        return (*this) + FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline FixedPoint_t<WHOLE, FRACT>
    operator+=(FixedPoint_t<WHOLE, FRACT> const &rhs)
    {
        // Compute summation and return result.
        return (*this) = (*this) + rhs;
    }

    inline FixedPoint_t<WHOLE, FRACT> operator+=(double rhs) const
    {
        return (*this) = (*this) + rhs;
    }

    // ========================================================================
    // SUB
    // ========================================================================
    inline FixedPoint_t<WHOLE, FRACT> sub(FixedPoint_t<WHOLE, FRACT> const &rhs,
                                          BitVector<WHOLE + FRACT> &op1,
                                          BitVector<WHOLE + FRACT> &op2) const
    {
        // Recombine the operators.
        recombine(whole, fractional, op1);
        recombine(rhs.get_whole(), rhs.get_fractional(), op2);
        // Compute subtraction and return result.
        return FixedPoint_t<WHOLE, FRACT>(op1 - op2);
    }

    inline FixedPoint_t<WHOLE, FRACT>
    operator-(FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        BitVector<WHOLE + FRACT> op1, op2;
        return this->sub(rhs, op1, op2);
    }

    inline FixedPoint_t<WHOLE, FRACT> operator-(double rhs) const
    {
        return (*this) - FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline FixedPoint_t<WHOLE, FRACT>
    operator-=(FixedPoint_t<WHOLE, FRACT> const &rhs)
    {
        // Compute subtraction and return result.
        return (*this) = (*this) - rhs;
    }

    inline FixedPoint_t<WHOLE, FRACT> &operator-=(double rhs)
    {
        return (*this) = (*this) - rhs;
    }

    // ========================================================================
    // MUL
    // ========================================================================
    inline FixedPoint_t<WHOLE, FRACT> mul(
        FixedPoint_t<WHOLE, FRACT> const &rhs,
        BitVector<WHOLE * 2 + WHOLE * 2> &full_precision,
        BitVector<WHOLE + FRACT> &op1,
        BitVector<WHOLE + FRACT> &op2) const
    {
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
        full_precision = op1 * op2;
        // Apply sign.
        if (op1_neg != op2_neg)
            full_precision.two_complement();
        // Change to string.
        const std::string str = full_precision.to_string();
        // Split whole and fractional.
        return FixedPoint_t<WHOLE, FRACT>(
            str.substr(WHOLE, WHOLE),
            str.substr(WHOLE * 2, FRACT));
    }

    inline FixedPoint_t<WHOLE, FRACT>
    operator*(FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        BitVector<WHOLE * 2 + WHOLE * 2> full_precision;
        BitVector<WHOLE + FRACT> op1, op2;
        return this->mul(rhs, full_precision, op1, op2);
    }

    inline FixedPoint_t<WHOLE, FRACT> operator*(double rhs) const
    {
        return (*this) * FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline FixedPoint_t<WHOLE, FRACT>
    operator*=(FixedPoint_t<WHOLE, FRACT> const &rhs)
    {
        return (*this) = (*this) * rhs;
    }

    inline FixedPoint_t<WHOLE, FRACT> &operator*=(double rhs)
    {
        return (*this) = (*this) * rhs;
    }

    // ========================================================================
    // DIV
    // ========================================================================
    inline FixedPoint_t<WHOLE, FRACT> div(
        FixedPoint_t<WHOLE, FRACT> const &rhs,
        BitVector<WHOLE + FRACT * 2> &op1,
        BitVector<WHOLE + FRACT * 2> &op2,
        BitVector<WHOLE + FRACT * 2> &support) const
    {
        // Recombine whole and fractional.
        op1 = recombine(whole, fractional);
        op2 = recombine(rhs.get_whole(), rhs.get_fractional());

        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();

        // Change signes.
        if (op1_neg)
            op1.two_complement();
        if (op2_neg)
            op2.two_complement();

        // Resize both bit-vectors.
        op1 = (rpad(op1.to_string(), WHOLE + FRACT * 2, '0'));

        // Perform multiplication.
        support = op1 / op2;

        // Apply sign.
        if (op1_neg != op2_neg)
            support.two_complement();

        // Change to string.
        const std::string str = support.to_string();

        // Split whole and fractional.
        return FixedPoint_t<WHOLE, FRACT>(
            str.substr(FRACT, WHOLE),
            str.substr(WHOLE + FRACT, FRACT));
    }

    inline FixedPoint_t<WHOLE, FRACT> operator/(
        FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        BitVector<WHOLE + FRACT * 2> op1, op2, support;
        return this->div(rhs, op1, op2, support);
    }

    inline FixedPoint_t<WHOLE, FRACT> operator/(double rhs) const
    {
        return (*this) / FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline FixedPoint_t<WHOLE, FRACT> operator/=(
        FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        return (*this) = (*this) / rhs;
    }

    inline FixedPoint_t<WHOLE, FRACT> &operator/=(double rhs)
    {
        return (*this) = (*this) / rhs;
    }

    // ========================================================================
    // EQUALITY OPs
    // ========================================================================
    inline bool operator==(FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        return (whole == rhs.whole) && (fractional == rhs.fractional);
    }

    inline bool operator==(double rhs) const
    {
        return (*this) == FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline bool operator!=(FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        return (whole != rhs.whole) || (fractional != rhs.fractional);
    }

    inline bool operator!=(double rhs) const
    {
        return (*this) != FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline bool operator<(FixedPoint_t<WHOLE, FRACT> const &rhs) const
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
        return (*this) < FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline bool operator<=(FixedPoint_t<WHOLE, FRACT> const &rhs) const
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
        return (*this) <= FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline bool operator>(FixedPoint_t<WHOLE, FRACT> const &rhs) const
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
        return (*this) > FixedPoint_t<WHOLE, FRACT>(rhs);
    }

    inline bool operator>=(FixedPoint_t<WHOLE, FRACT> const &rhs) const
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
        return (*this) >= FixedPoint_t<WHOLE, FRACT>(rhs);
    }
};

// ========================================================================
// ARITHMETIC OPERATIONS
// ========================================================================
template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint_t<WHOLE, FRACT> operator+(
    double lhs,
    FixedPoint_t<WHOLE, FRACT> const &rhs)
{
    return FixedPoint_t<WHOLE, FRACT>(lhs) + rhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint_t<WHOLE, FRACT> operator-(
    double lhs,
    FixedPoint_t<WHOLE, FRACT> const &rhs)
{
    return FixedPoint_t<WHOLE, FRACT>(lhs) - rhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint_t<WHOLE, FRACT> operator*(
    double lhs,
    FixedPoint_t<WHOLE, FRACT> const &rhs)
{
    return FixedPoint_t<WHOLE, FRACT>(lhs) * rhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint_t<WHOLE, FRACT> operator/(
    double lhs,
    FixedPoint_t<WHOLE, FRACT> const &rhs)
{
    static_assert(FRACT != 0, "You are trying to divide with a fractional "
                              "precision of zero!");
    return FixedPoint_t<WHOLE, FRACT>(lhs) / rhs;
}

// ========================================================================
// OPERATORS FOR OSTREAM
// ========================================================================
template <long unsigned int WHOLE, long unsigned int FRACT>
std::ostream &operator<<(std::ostream &lhs,
                         FixedPoint_t<WHOLE, FRACT> const &rhs)
{
    lhs << rhs.to_number();
    return lhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
std::stringstream &operator<<(std::stringstream &lhs,
                              FixedPoint_t<WHOLE, FRACT> const &rhs)
{
    lhs << rhs.to_number();
    return lhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
std::ifstream &operator>>(std::ifstream &lhs,
                          FixedPoint_t<WHOLE, FRACT> &rhs)
{
    std::string s;
    lhs >> s;
    size_t size;
    rhs = std::stof(s, &size);
    return lhs;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
std::stringstream &operator>>(std::stringstream &lhs,
                              FixedPoint_t<WHOLE, FRACT> &rhs)
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
FixedPoint_t<WHOLE, FRACT> round(FixedPoint_t<WHOLE, FRACT> const &fp)
{
    FixedPoint_t<WHOLE, FRACT> rounded = fp;
    if ((FRACT >= 1) && rounded.get_fractional().at(FRACT - 1))
        rounded = rounded + 1;
    rounded.set_fractional(0);
    return rounded;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
FixedPoint_t<WHOLE, FRACT> abs(FixedPoint_t<WHOLE, FRACT> const &fp)
{
    if (!fp.get_whole().sign())
        return fp;
    return recombine(fp.get_whole(), fp.get_fractional()).two_complement();
}

template <long unsigned int WHOLE, long unsigned int FRACT>
FixedPoint_t<WHOLE, FRACT> sqrt(
    FixedPoint_t<WHOLE, FRACT> const &value,
    long max_iterations              = 20,
    FixedPoint_t<WHOLE, FRACT> error = FixedPoint_t<WHOLE, FRACT>(0.001))
{
    FixedPoint_t<WHOLE, FRACT> result(1.0);
    int iteration = 0;
    while ((iteration++ < max_iterations) &&
           (abs((result * result) - value) >= error))
        result = 0.5f * (result + (value / result));
    return result;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
FixedPoint_t<WHOLE, FRACT> pow(
    FixedPoint_t<WHOLE, FRACT> const &value,
    unsigned exponent)
{
    if (exponent == 0)
        return FixedPoint_t<WHOLE, FRACT>(1);
    FixedPoint_t<WHOLE, FRACT> result = value;
    for (unsigned it = 1; it < exponent; ++it)
        result *= value;
    return result;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline FixedPoint_t<WHOLE, FRACT> exp(FixedPoint_t<WHOLE, FRACT> x,
                                      unsigned prec = 12)
{
    x = 1.0 + (x / std::pow(2, prec));
    for (unsigned i = 0; i < prec; ++i)
        x *= x;
    return x;
}