#include "fplib/fixedpoint.hpp"
#include "fplib/math.hpp"

#include <iostream>
#include <iomanip>
#include <random>
#include <vector>

int main(int, char *[])
{
    double __a = 0.125, __b = 2;
    fplib::FixedPoint<5, 4> a = __a;
    fplib::FixedPoint<3, 6> b = __b;

    std::cout << a << " (" << a.to_string() << ")\n";
    std::cout << b << " (" << b.to_string() << ")\n";
    std::cout << " " << __a << " / " << __b << " = " << (__a / __b) << "\n";
    std::cout << " " << a << " / " << b << " = " << (a / b) << "\n";

#if 0
    double a = 3, b = 7;
    fplib::FixedPoint<4, 4> l_a   = a;
    fplib::FixedPoint<4, 4> b   = b;
    //fplib::FixedPoint<16, 16> m_a = a;
    //fplib::FixedPoint<16, 16> m_b = b;
    //fplib::FixedPoint<32, 32> h_a = a;
    //fplib::FixedPoint<32, 32> h_b = b;

    std::cout << " " << a << " + " << b << " = " << (a + b) << "\n";
    std::cout << "    Low    # bits : " << (l_a + b) << "\n";
    //std::cout << "    Medium # bits : " << (m_a + m_b) << "\n";
    //std::cout << "    High   # bits : " << (h_a + h_b) << "\n";

    std::cout << " " << a << " - " << b << " = " << (a - b) << "\n";
    std::cout << "    Low    # bits : " << (l_a - b) << "\n";
    //std::cout << "    Medium # bits : " << (m_a - m_b) << "\n";
    //std::cout << "    High   # bits : " << (h_a - h_b) << "\n";

    std::cout << " " << a << " * " << b << " = " << (a * b) << "\n";
    std::cout << "    Low    # bits : " << (l_a * b) << "\n";
    //std::cout << "    Medium # bits : " << (m_a * m_b) << "\n";
    //std::cout << "    High   # bits : " << (h_a * h_b) << "\n";

    std::cout << " " << a << " / " << b << " = " << (a / b) << "\n";
    auto r0 = __div(l_a, b);
    //auto r1 = __div(m_a, m_b);
    //auto r2 = __div(h_a, h_b);
    std::cout << "\n";
    std::cout << "    Low    # bits : " << r0 << "\n";
    //std::cout << "    Medium # bits : " << r1 << "\n";
    //std::cout << "    High   # bits : " << r2 << "\n";
#endif
    return 0;
}