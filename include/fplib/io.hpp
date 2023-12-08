/// @file io.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Input and output stream operators.

#pragma once

#include <iostream>
#include <sstream>
#include <fstream>

template <std::size_t WHOLE, std::size_t FRACT>
std::ostream &operator<<(std::ostream &lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    lhs << rhs.to_number();
    return lhs;
}

template <std::size_t WHOLE, std::size_t FRACT>
std::stringstream &operator<<(std::stringstream &lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    lhs << rhs.to_number();
    return lhs;
}

template <std::size_t WHOLE, std::size_t FRACT>
std::ifstream &operator>>(std::ifstream &lhs, fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    double value;
    lhs >> value;
    rhs = value;
    return lhs;
}

template <std::size_t WHOLE, std::size_t FRACT>
std::stringstream &operator>>(std::stringstream &lhs, fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    double value;
    lhs >> value;
    rhs = value;
    return lhs;
}
