/// @file math.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include "fixedpoint.hpp"
#include "support.hpp"

namespace fplib
{

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline auto sum(const FixedPoint<WHOLE1, FRACT1> &lhs, const FixedPoint<WHOLE2, FRACT2> &rhs)
{
    constexpr std::size_t whole_size = std::max(WHOLE1, WHOLE2);
    constexpr std::size_t fract_size = std::max(FRACT1, FRACT2);
    // Create bvlib::BitVector for support.
    bvlib::BitVector<WHOLE1 + FRACT1> op1;
    bvlib::BitVector<WHOLE2 + FRACT2> op2;
    // Recombine the operators.
    fplib::recombine(lhs.get_whole(), lhs.get_fractional(), op1);
    fplib::recombine(rhs.get_whole(), rhs.get_fractional(), op2);
    // Compute the sum and return result.
    return FixedPoint<whole_size, fract_size>(op1 + op2);
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline auto sub(const FixedPoint<WHOLE1, FRACT1> &lhs, const FixedPoint<WHOLE2, FRACT2> &rhs)
{
    constexpr std::size_t whole_size = std::max(WHOLE1, WHOLE2);
    constexpr std::size_t fract_size = std::max(FRACT1, FRACT2);
    // Create bvlib::BitVector for support.
    bvlib::BitVector<WHOLE1 + FRACT1> op1;
    bvlib::BitVector<WHOLE2 + FRACT2> op2;
    // Recombine the operators.
    fplib::recombine(lhs.get_whole(), lhs.get_fractional(), op1);
    fplib::recombine(rhs.get_whole(), rhs.get_fractional(), op2);
    // Compute the sum and return result.
    return FixedPoint<whole_size, fract_size>(op1 - op2);
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline auto mul(const FixedPoint<WHOLE1, FRACT1> &lhs, const FixedPoint<WHOLE2, FRACT2> &rhs)
{
    constexpr std::size_t whole_size = std::max(WHOLE1, WHOLE2);
    constexpr std::size_t fract_size = std::max(FRACT1, FRACT2);

    // Create bvlib::BitVector for support.
    bvlib::BitVector<WHOLE1 + FRACT1> op1;
    bvlib::BitVector<WHOLE2 + FRACT2> op2;
    // Recombine the operators.
    fplib::recombine(lhs.get_whole(), lhs.get_fractional(), op1);
    fplib::recombine(rhs.get_whole(), rhs.get_fractional(), op2);
    // Get the sign.
    bool op1_neg = lhs.get_whole().sign();
    bool op2_neg = rhs.get_whole().sign();
    // Change signes.
    if (op1_neg)
        op1.two_complement();
    if (op2_neg)
        op2.two_complement();
    // Compute multiplication with full precision.
    auto support = bvlib::mul(op1, op2);
    // Apply sign.
    if (op1_neg != op2_neg)
        support.two_complement();
    // Change to string.
    const std::string str = support.to_string();
    // Return the fixed point value by splitting the string.
    return fplib::FixedPoint<whole_size, fract_size>(
        str.substr((support.size() / 2) - whole_size, whole_size),
        str.substr((support.size() / 2), fract_size));
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline auto div(const FixedPoint<WHOLE1, FRACT1> &lhs, const FixedPoint<WHOLE2, FRACT2> &rhs)
{
    constexpr std::size_t WHOLE = std::max(WHOLE1, WHOLE2);
    constexpr std::size_t FRACT = std::max(FRACT1, FRACT2);
    // Create bit vectors for support.
    bvlib::BitVector<2 * WHOLE + FRACT> op1, op2;
    bvlib::BitVector<WHOLE + FRACT> result;
    // Recombine the operators.
    fplib::recombine<WHOLE1, FRACT1, 2 * WHOLE, FRACT>(lhs.get_whole(), lhs.get_fractional(), op1);
    fplib::recombine<WHOLE2, FRACT2, 2 * WHOLE, FRACT>(rhs.get_whole(), rhs.get_fractional(), op2);
    // Get the sign.
    bool op1_neg = lhs.get_whole().sign();
    bool op2_neg = rhs.get_whole().sign();
    // Change signes.
    if (op1_neg) {
        for (std::size_t i = 0; i < (2 * WHOLE-WHOLE1);++i)
            op1.bits[i] = 1;
        op1.two_complement();
    }
    if (op2_neg) {
        for (std::size_t i = 0; i < (2 * WHOLE-WHOLE2);++i)
            op2.bits[i] = 1;
        op2.two_complement();
    }
    // Perform multiplication.
    auto [quotient, remainder] = bvlib::div(op1 << FRACT, op2);
    (void) remainder;
    // Apply sign.
    if (op1_neg != op2_neg)
        quotient.two_complement();
    // Change to string.
    for (std::size_t it = 0; it < (WHOLE + FRACT); ++it)
        result[it] = quotient[it];
    // Split whole and fractional.
    return FixedPoint<WHOLE, FRACT>(result);
}

} // namespace fplib

// ========================================================================
// ARITHMETIC OPERATIONS (fixedpoint & fixedpoint)
// ========================================================================
template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline auto operator+(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    return fplib::sum(lhs, rhs);
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline auto operator-(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    return fplib::sub(lhs, rhs);
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline auto operator*(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    return fplib::mul(lhs, rhs);
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline auto operator/(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    return fplib::div(lhs, rhs);
}

// ========================================================================
// ARITHMETIC OPERATIONS (scalar & fixedpoint)
// ========================================================================
template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline auto operator+(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::sum(fplib::FixedPoint<WHOLE, FRACT>(lhs), rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline auto operator+(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return fplib::sum(lhs, fplib::FixedPoint<WHOLE, FRACT>(rhs));
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline auto operator-(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::sub(fplib::FixedPoint<WHOLE, FRACT>(lhs), rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline auto operator-(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return fplib::sub(lhs, fplib::FixedPoint<WHOLE, FRACT>(rhs));
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline auto operator*(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::mul(fplib::FixedPoint<WHOLE, FRACT>(lhs), rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline auto operator*(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return fplib::mul(lhs, fplib::FixedPoint<WHOLE, FRACT>(rhs));
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline auto operator/(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    static_assert(FRACT != 0, "You are trying to divide with a fractional precision of zero!");
    return fplib::div(fplib::FixedPoint<WHOLE, FRACT>(lhs), rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline auto operator/(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    static_assert(FRACT != 0, "You are trying to divide with a fractional precision of zero!");
    return fplib::div(lhs, fplib::FixedPoint<WHOLE, FRACT>(rhs));
}

// ========================================================================
// EQUALITY OPERATIONS (fixedpoint & fixedpoint)
// ========================================================================
template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline bool operator==(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    return (WHOLE1 == WHOLE2) && (FRACT1 == FRACT2) && (lhs.get_whole() == rhs.get_whole()) && (lhs.get_fractional() == rhs.get_fractional());
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline bool operator!=(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    return (WHOLE1 != WHOLE2) || (FRACT1 != FRACT2) || (lhs.get_whole() != rhs.get_whole()) || (lhs.get_fractional() != rhs.get_fractional());
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline bool operator<(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    // Get the sign.
    bool op1_neg = lhs.get_whole().sign(), op2_neg = rhs.get_whole().sign();
    // Check the sign first, it's faster.
    if (op1_neg && !op2_neg)
        return true;
    if (!op1_neg && op2_neg)
        return false;
    // Recombine whole and fractional.
    return fplib::recombine(lhs.get_whole(), lhs.get_fractional()) < fplib::recombine(rhs.get_whole(), rhs.get_fractional());
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline bool operator>(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    // Get the sign.
    bool op1_neg = lhs.get_whole().sign(), op2_neg = rhs.get_whole().sign();
    // Check the sign first, it's faster.
    if (!op1_neg && op2_neg)
        return true;
    if (op1_neg && !op2_neg)
        return false;
    // Recombine whole and fractional.
    return fplib::recombine(lhs.get_whole(), lhs.get_fractional()) > fplib::recombine(rhs.get_whole(), rhs.get_fractional());
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline bool operator<=(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    // Get the sign.
    bool op1_neg = lhs.get_whole().sign(), op2_neg = rhs.get_whole().sign();
    // Check the sign first, it's faster.
    if (op1_neg && !op2_neg)
        return true;
    if (!op1_neg && op2_neg)
        return false;
    // Recombine whole and fractional.
    return fplib::recombine(lhs.get_whole(), lhs.get_fractional()) <= fplib::recombine(rhs.get_whole(), rhs.get_fractional());
}

template <std::size_t WHOLE1, std::size_t FRACT1, std::size_t WHOLE2, std::size_t FRACT2>
inline bool operator>=(const fplib::FixedPoint<WHOLE1, FRACT1> &lhs, const fplib::FixedPoint<WHOLE2, FRACT2> &rhs)
{
    // Get the sign.
    bool op1_neg = lhs.get_whole().sign(), op2_neg = rhs.get_whole().sign();
    // Check the sign first, it's faster.
    if (!op1_neg && op2_neg)
        return true;
    if (op1_neg && !op2_neg)
        return false;
    // Recombine whole and fractional.
    return fplib::recombine(lhs.get_whole(), lhs.get_fractional()) >= fplib::recombine(rhs.get_whole(), rhs.get_fractional());
}

// ========================================================================
// EQUALITY OPERATIONS (scalar & fixedpoint)
// ========================================================================
template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator==(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return lhs == fplib::FixedPoint<WHOLE, FRACT>(rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator==(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::FixedPoint<WHOLE, FRACT>(lhs) == rhs;
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator!=(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return lhs != fplib::FixedPoint<WHOLE, FRACT>(rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator!=(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::FixedPoint<WHOLE, FRACT>(lhs) != rhs;
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator>(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return lhs > fplib::FixedPoint<WHOLE, FRACT>(rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator>(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::FixedPoint<WHOLE, FRACT>(lhs) > rhs;
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator<(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return lhs < fplib::FixedPoint<WHOLE, FRACT>(rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator<(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::FixedPoint<WHOLE, FRACT>(lhs) < rhs;
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator>=(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return lhs >= fplib::FixedPoint<WHOLE, FRACT>(rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator>=(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::FixedPoint<WHOLE, FRACT>(lhs) >= rhs;
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator<=(const fplib::FixedPoint<WHOLE, FRACT> &lhs, T rhs)
{
    return lhs <= fplib::FixedPoint<WHOLE, FRACT>(rhs);
}

template <typename T, std::size_t WHOLE, std::size_t FRACT>
inline bool operator<=(T lhs, const fplib::FixedPoint<WHOLE, FRACT> &rhs)
{
    return fplib::FixedPoint<WHOLE, FRACT>(lhs) <= rhs;
}

// ========================================================================
// FUNCTIONS
// ========================================================================
template <std::size_t WHOLE, std::size_t FRACT>
inline auto round(const fplib::FixedPoint<WHOLE, FRACT> &fp)
{
    fplib::FixedPoint<WHOLE, FRACT> rounded = fp;
    if ((FRACT >= 1) && rounded.get_fractional().at(FRACT - 1))
        rounded = rounded + 1;
    rounded.set_fractional(0);
    return rounded;
}

template <std::size_t WHOLE, std::size_t FRACT>
inline auto abs(const fplib::FixedPoint<WHOLE, FRACT> &fp)
{
    if (!fp.get_whole().sign())
        return fp;
    return recombine(fp.get_whole(), fp.get_fractional()).two_complement();
}

template <std::size_t WHOLE, std::size_t FRACT>
inline auto sqrt(
    const fplib::FixedPoint<WHOLE, FRACT> &value,
    long max_iterations                   = 20,
    fplib::FixedPoint<WHOLE, FRACT> error = fplib::FixedPoint<WHOLE, FRACT>(0.001))
{
    fplib::FixedPoint<WHOLE, FRACT> result(1.0);
    int iteration = 0;
    while ((iteration++ < max_iterations) &&
           (abs((result * result) - value) >= error))
        result = 0.5f * (result + (value / result));
    return result;
}

template <std::size_t WHOLE, std::size_t FRACT>
inline auto pow(const fplib::FixedPoint<WHOLE, FRACT> &value, unsigned exponent)
{
    if (exponent == 0)
        return fplib::FixedPoint<WHOLE, FRACT>(1);
    fplib::FixedPoint<WHOLE, FRACT> result = value;
    for (unsigned it = 1; it < exponent; ++it)
        result *= value;
    return result;
}

template <std::size_t WHOLE, std::size_t FRACT>
inline auto exp(fplib::FixedPoint<WHOLE, FRACT> x, unsigned prec = 12)
{
    x = 1.0 + (x / std::pow(2, prec));
    for (unsigned i = 0; i < prec; ++i)
        x *= x;
    return x;
}