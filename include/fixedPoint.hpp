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

template <long unsigned int FRACT>
inline double fractional_to_float(BitVector<FRACT> const &_fract)
{
    double result = 0;
    for (long i = FRACT - 1; i >= 0; --i)
        if (_fract[i])
            result += 1 / std::pow(2, FRACT - i);
    return result;
}

template <long unsigned int FRACT>
inline BitVector<FRACT> float_to_fractional(double value)
{
    BitVector<FRACT> result;
    double acc = 0;
    for (long it = (FRACT - 1); it >= 0; --it) {
        double element = 1.0 / (1U << (FRACT - it));
        if ((acc + element) <= value) {
            result.flip(it);
            acc += element;
        }
    }
    return result;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline BitVector<WHOLE + FRACT> recombine(
    BitVector<WHOLE> const &whole,
    BitVector<FRACT> const &fractional)
{
    BitVector<WHOLE + FRACT> full;
    for (size_t it = 0; it < WHOLE; ++it)
        full[it + FRACT] = whole[it];
    for (size_t it = 0; it < FRACT; ++it)
        full[it] = fractional[it];
    return full;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline void split(
    BitVector<WHOLE + FRACT> const &full,
    BitVector<WHOLE> &whole,
    BitVector<FRACT> &fractional)
{
    for (size_t it = 0; it < WHOLE; ++it)
        whole[it] = full[it + FRACT];
    for (size_t it = 0; it < FRACT; ++it)
        fractional[it] = full[it];
}

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
            BitVector<WHOLE + FRACT> full(
                whole.to_string() + fractional.to_string());
            full.two_complement();
            whole      = BitVector<WHOLE>(full.to_string().substr(0, WHOLE));
            fractional = BitVector<FRACT>(
                full.to_string().substr(WHOLE, FRACT));
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
        BitVector<WHOLE + FRACT> full = recombine(whole, fractional);
        full.two_complement();
        BitVector<WHOLE> _whole;
        BitVector<FRACT> _fract;
        split(full, _whole, _fract);
        return -(_whole.to_number() + fractional_to_float<FRACT>(_fract));
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
    inline FixedPoint_t<WHOLE, FRACT>
    operator+(FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        // Compute summation and return result.
        return FixedPoint_t<WHOLE, FRACT>(
            recombine(whole, fractional) +
            recombine(rhs.get_whole(), rhs.get_fractional()));
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
    inline FixedPoint_t<WHOLE, FRACT>
    operator-(FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        // Compute subtraction and return result.
        return FixedPoint_t<WHOLE, FRACT>(
            recombine(whole, fractional) -
            recombine(rhs.get_whole(), rhs.get_fractional()));
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
    inline FixedPoint_t<WHOLE, FRACT>
    operator*(FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        // Create bitvector for solution.
        BitVector<WHOLE * 2 + FRACT * 2> mul;
        // Recombine whole and fractional.
        auto op1 = recombine(whole, fractional);
        auto op2 = recombine(rhs.get_whole(), rhs.get_fractional());
        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();
        // Change signes.
        if (op1_neg)
            op1.two_complement();
        if (op2_neg)
            op2.two_complement();
        // Compute multiplication.
        mul = op1 * op2;
        // Apply sign.
        if (op1_neg != op2_neg)
            mul.two_complement();
        // Change to string.
        const std::string str = mul.to_string();
        // Split whole and fractional.
        return FixedPoint_t<WHOLE, FRACT>(
            str.substr(WHOLE, WHOLE),
            str.substr(WHOLE * 2, FRACT));
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
    inline FixedPoint_t<WHOLE, FRACT> operator/(
        FixedPoint_t<WHOLE, FRACT> const &rhs) const
    {
        // Create bitvector for solution.
        BitVector<WHOLE + FRACT * 2> div;

        // Recombine whole and fractional.
        auto _op1 = recombine(whole, fractional);
        auto _op2 = recombine(rhs.get_whole(), rhs.get_fractional());

        // Get the sign.
        bool op1_neg = whole.sign();
        bool op2_neg = rhs.get_whole().sign();

        // Change signes.
        if (op1_neg)
            _op1.two_complement();
        if (op2_neg)
            _op2.two_complement();

        // Resize both bit-vectors.
        BitVector<WHOLE + FRACT * 2> op1(
            rpad(_op1.to_string(), WHOLE + FRACT * 2, '0'));
        BitVector<WHOLE + FRACT * 2> op2(_op2.to_string());

        // Perform multiplication.
        div = op1 / op2;

        // Apply sign.
        if (op1_neg != op2_neg)
            div.two_complement();

        // Change to string.
        const std::string str = div.to_string();

        // Split whole and fractional.
        return FixedPoint_t<WHOLE, FRACT>(
            str.substr(FRACT, WHOLE),
            str.substr(WHOLE + FRACT, FRACT));
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
    FixedPoint_t<WHOLE, FRACT> positive = fp;
    if (fp.get_whole().sign()) {
        auto full = recombine(fp.get_whole(), fp.get_fractional());
        full.two_complement();
        positive = full;
    }
    return positive;
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