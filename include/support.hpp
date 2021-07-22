/// @file support.hpp
/// @brief Support functions.
/// @date 21/07/2021
/// @author Enrico Fraccaroli
/// @copyright
/// Copyright (c) 2019-2021 Enrico Fraccaroli <enrico.fraccaroli@gmail.com>
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

#include <string>

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

template <long unsigned int FRACT>
inline double fractional_to_float(BitVector<FRACT> const &fractional)
{
    double result = 0;
    for (long i = FRACT - 1; i >= 0; --i)
        if (fractional[i])
            result += 1 / std::pow(2, FRACT - i);
    return result;
}

template <long unsigned int FRACT>
inline BitVector<FRACT> &float_to_fractional(double value, BitVector<FRACT> &fractional)
{
    double acc = 0;
    for (long it = (FRACT - 1); it >= 0; --it) {
        double element = 1.0 / (1U << (FRACT - it));
        if ((acc + element) <= value) {
            fractional.flip(it);
            acc += element;
        }
    }
    return fractional;
}

template <long unsigned int FRACT>
inline BitVector<FRACT> float_to_fractional(double value)
{
    BitVector<FRACT> fractional;
    return float_to_fractional<FRACT>(value, fractional);
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline BitVector<WHOLE + FRACT> &recombine(
    BitVector<WHOLE> const &whole,
    BitVector<FRACT> const &fractional,
    BitVector<WHOLE + FRACT> &recombined)
{
    size_t it;
    for (it = 0; it < WHOLE; ++it)
        recombined[it + FRACT] = whole[it];
    for (it = 0; it < FRACT; ++it)
        recombined[it] = fractional[it];
    return recombined;
}

template <long unsigned int WHOLE, long unsigned int FRACT>
inline BitVector<WHOLE + FRACT> recombine(
    BitVector<WHOLE> const &whole,
    BitVector<FRACT> const &fractional)
{
    BitVector<WHOLE + FRACT> recombined;
    return recombine(whole, fractional, recombined);
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
