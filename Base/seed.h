#ifndef SEED_H
#define SEED_H

#include <string>

struct Seed {
    std::string str_;
    int ocnFg_;
    int ocnBg_;
    double FET_;

    Seed(std::string str, int ocnFg, int ocnBg, double FET):
        str_(str), ocnFg_(ocnFg),ocnBg_(ocnBg), FET_(FET) {}
};

#endif /* SEED_H */
