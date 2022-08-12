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

#include "bvlib/bitvector.hpp"
#include "bvlib/math.hpp"
#include "support.hpp"

#include <fstream>
#include <sstream>
#include <iostream>

namespace fplib
{

/// @brief The fixed point class.
/// @tparam WHOLE the size of the whole part.
/// @tparam FRACT the size of the fractional part.
template <std::size_t WHOLE, std::size_t FRACT>
class FixedPoint {
private:
    /// The original value.
    double original;
    /// Whole part.
    bvlib::BitVector<WHOLE> whole;
    /// Fractional part.
    bvlib::BitVector<FRACT> fractional;

public:
    const std::size_t whole_size = WHOLE;
    const std::size_t fract_size = FRACT;

    /// @brief Empty constructor.
    FixedPoint()
        : original(),
          whole(),
          fractional()
    {
        // Nothing to do.
    }

    /// @brief Empty constructor.
    FixedPoint(double val)
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
    FixedPoint(const bvlib::BitVector<WHOLE + FRACT> &_full)
        : original(),
          whole(),
          fractional()
    {
        split<WHOLE, FRACT>(_full, whole, fractional);
    }

    /// @brief Constructor from bit-vectors, e.g., <0011>,<1000> = 3.5.
    FixedPoint(const bvlib::BitVector<WHOLE> &_whole, const bvlib::BitVector<FRACT> &_fractional)
        : original(),
          whole(_whole),
          fractional(_fractional)
    {
        // Nothing to do.
    }

    /// @brief Constructor from bit-vectors, e.g., <0011>,<1000> = 3.5.
    template <std::size_t WHOLE2, std::size_t FRACT2>
    FixedPoint(const bvlib::BitVector<WHOLE2> &_whole, const bvlib::BitVector<FRACT2> &_fractional)
        : original(),
          whole(_whole),
          fractional(_fractional)
    {
        // Nothing to do.
    }

    /// @brief Constructor from another fixed point.
    FixedPoint(const FixedPoint<WHOLE, FRACT> &other)
        : original(other.original),
          whole(other.whole),
          fractional(other.fractional)
    {
        // Nothing to do.
    }

    /// @brief Constructor from another fixed point.
    template <std::size_t WHOLE2, std::size_t FRACT2>
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

    inline const bvlib::BitVector<WHOLE> &get_whole() const
    {
        return whole;
    }

    inline bool get_whole(std::size_t position) const
    {
        return whole[position];
    }

    inline const bvlib::BitVector<FRACT> &get_fractional() const
    {
        return fractional;
    }

    inline bool get_fractional(std::size_t position) const
    {
        return fractional[position];
    }

    inline void set_whole(const bvlib::BitVector<WHOLE> &_whole)
    {
        whole = _whole;
    }

    inline void set_whole(std::size_t _whole)
    {
        whole = _whole;
    }

    inline void set_whole(const std::string &_whole)
    {
        whole = lpad(_whole, WHOLE, '0');
    }

    inline void set_fractional(const bvlib::BitVector<FRACT> &_fractional)
    {
        fractional = _fractional;
    }

    inline void set_fractional(std::size_t _fractional)
    {
        fractional = _fractional;
    }

    inline void set_fractional(const std::string &_fractional)
    {
        fractional = rpad(_fractional, FRACT, '0');
    }

    static inline constexpr std::size_t get_whole_max()
    {
        return bvlib::BitVector<WHOLE>::ones().reset(WHOLE - 1).to_number();
    }

    std::string to_string() const
    {
        return whole.to_string() + "." + fractional.to_string();
    }

    double to_number() const
    {
        if (!whole.sign())
            return whole.template to_number<double>() + fractional_to_float<FRACT>(fractional);

        // Support variables.
        bvlib::BitVector<WHOLE + FRACT> _recombined;
        bvlib::BitVector<WHOLE> _whole;
        bvlib::BitVector<FRACT> _fractional;

        // 1. Recombine whole and fract.
        recombine(whole, fractional, _recombined);

        // 2. Perform 2's complement.
        _recombined.two_complement();

        // 3. Split the complemented value back into whole and fractional.
        split(_recombined, _whole, _fractional);
        return -(_whole.template to_number<double>() + fractional_to_float<FRACT>(_fractional));
    }

    // ========================================================================
    // ASSIGN
    // ========================================================================
    inline auto &operator=(const FixedPoint<WHOLE, FRACT> &rhs)
    {
        // Copy the original value.
        original = rhs.get_original();
        // Copy whole part.
        whole.assign(rhs.get_whole());
        // Copy fractional part.
        fractional.assign(rhs.get_fractional());
        return (*this);
    }

    template <std::size_t WHOLE2, std::size_t FRACT2>
    inline auto &operator=(const FixedPoint<WHOLE2, FRACT2> &rhs)
    {
        // Copy the original value.
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

    inline auto &operator=(double rhs)
    {
        original = rhs;
        this->_set_from_double(rhs);
        return (*this);
    }

    inline auto &operator=(const bvlib::BitVector<WHOLE + FRACT> &rhs)
    {
        split<WHOLE, FRACT>(rhs, whole, fractional);
        return *this;
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
        whole = static_cast<std::size_t>(intpart);
        // Set fractional part.
        fractional = float_to_fractional<FRACT>(fractpart);
        // Apply the sign.
        if (sign) {
            bvlib::BitVector<WHOLE + FRACT> full = recombine(whole, fractional);
            full.two_complement();
            split(full, whole, fractional);
        }
    }
};

} // namespace fplib
