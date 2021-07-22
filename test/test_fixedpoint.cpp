#include <iostream>
#include <iomanip>
#include <vector>

#include "fixedPoint.hpp"

#define TEST_OP(op)                                                     \
    std::cout << "== Testing " #op " " << std::string(25, '=') << "\n"; \
    for (size_t it = 0; it < inputs.size(); ++it) {                     \
        for (size_t it2 = 0; it2 < inputs.size(); ++it2) {              \
            auto a = inputs[it], b = inputs[it2];                       \
            auto c = a.get_original(), d = b.get_original();            \
            std::cout << std::setw(12) << a << " " #op " "              \
                      << std::setw(12) << b << " = "                    \
                      << std::setw(12) << (a op b) << " == "            \
                      << std::setw(12) << (c op d) << "\n";             \
        }                                                               \
    }

int main(int argc, char *argv[])
{
    std::vector<FixedPoint<32, 32>> inputs;
    inputs.emplace_back(FixedPoint<32, 32>(+6.5625));
    inputs.emplace_back(FixedPoint<32, 32>(+4.2500));
    inputs.emplace_back(FixedPoint<32, 32>(-2.2500));
    inputs.emplace_back(FixedPoint<32, 32>(-4.2500));
    inputs.emplace_back(FixedPoint<32, 32>(321.6543));
    inputs.emplace_back(FixedPoint<32, 32>(76.7657));
    inputs.emplace_back(FixedPoint<32, 32>(1023.2156346));
    inputs.emplace_back(FixedPoint<32, 32>(-653.864786786));

    std::cout << std::setprecision(4) << std::fixed;
    std::cout << "== INPUTS " << std::string(25, '=') << "\n";
    for (size_t it = 0; it < inputs.size(); ++it) {
        std::cout << "op" << it << " ["
                  << inputs[it].to_string() << "] "
                  << inputs[it].to_number()
                  << " original(" << inputs[it].get_original() << ")"
                  << "\n";
    }
    std::cout << std::string(35, '=') << "\n";
    TEST_OP(+)
    TEST_OP(-)
    TEST_OP(*)
    TEST_OP(/)
    std::cout << std::string(35, '=') << "\n";
    TEST_OP(<)
    TEST_OP(<=)
    TEST_OP(>)
    TEST_OP(>=)
    TEST_OP(==)
    TEST_OP(!=)
    std::cout << std::string(35, '=') << "\n";
    return 0;
}