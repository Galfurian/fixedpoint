#include "fplib/fixedpoint.hpp"

#include <iostream>
#include <iomanip>
#include <random>
#include <vector>

template <std::size_t WHOLE, std::size_t FRACT>
struct test_fp_pair_t {
    fplib::FixedPoint<WHOLE, FRACT> fixedpoint;
    double original;
    test_fp_pair_t(double _original = 0)
        : fixedpoint(_original),
          original(_original)
    {
    }
};

#if 0

template <
    std::size_t WHOLE1,
    std::size_t FRACT1,
    std::size_t WHOLE2,
    std::size_t FRACT2,
    unsigned long QT>
int test_operators()
{
    test_fp_pair_t<WHOLE1, FRACT1> inputs1[QT];
    test_fp_pair_t<WHOLE2, FRACT2> inputs2[QT];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distr(5, 5);

    for (unsigned long i = 0; i < QT; ++i) {
        inputs1[i] = test_fp_pair_t<WHOLE1, FRACT1>(distr(gen));
        inputs2[i] = test_fp_pair_t<WHOLE2, FRACT2>(distr(gen));
    }

    std::stringstream ss;
    ss << "<" << WHOLE1 << ", " << FRACT1 << "><" << WHOLE2 << ", " << FRACT2 << ">";
    auto sizes = ss.str();

    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint + inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original + inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " + " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " + " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            if (inputs1[i].original < inputs2[j].original)
                continue;
            auto bv_result = inputs1[i].fixedpoint - inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original - inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " - " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " - " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint * inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original * inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " * " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits) != \n";
                std::cerr << inputs1[i].original << " * " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint / inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original / inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " / " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " / " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint < inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original < inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " < " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " < " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint <= inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original <= inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " <= " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " <= " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint > inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original > inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " > " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " > " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint >= inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original >= inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " >= " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " >= " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint == inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original == inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " == " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " == " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    for (unsigned long i = 0; i < QT; ++i) {
        for (unsigned long j = 0; j < QT; ++j) {
            auto bv_result = inputs1[i].fixedpoint != inputs2[j].fixedpoint;
            auto nm_result = inputs1[i].original != inputs2[j].original;
            if (bv_result != nm_result) {
                std::cerr << inputs1[i].fixedpoint << " != " << inputs2[j].fixedpoint << " = " << bv_result << " (" << sizes << " bits)!= \n";
                std::cerr << inputs1[i].original << " != " << inputs2[j].original << " = " << nm_result << "\n\n";
                return 1;
            }
        }
    }
    return 0;
}

int main(int, char *[])
{
#if 0
    test_operators<8, 8, 8, 8, 50>();
    test_operators<16, 16, 16, 16, 50>();
    test_operators<32, 32, 32, 32, 50>();
    test_operators<64, 64, 64, 64, 50>();
#endif
    return 0;
}

#else

int main(int, char *[])
{
    return 0;
}

#endif