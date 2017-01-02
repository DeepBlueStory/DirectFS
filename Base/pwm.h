#ifndef PWM_H
#define PWM_H

#include <deque>
#include "seed.h"
#include <cmath>

struct PWM {
    std::deque<double> pwm_;
    double b_;
    int ocnFg_;
    int ocnBg_;
    double FET_;

    PWM();

    PWM(const Seed seed, const double b)
        : pwm_(std::deque<double>(seed.str_.size() * 4, log(0 + minconst))),
        b_(b), ocnFg_(seed.ocnFg_), ocnBg_(seed.ocnBg_), FET_(seed.FET_)
    {
        double log1 = log(1 + minconst);
        const std::string& str = seed.str_;
        for(int i = 0; i < (int)seed.str_.size(); ++i)
        {
            switch(str[i])
            {
                case 'A' :
                    pwm_[i * 4 + 0] = log1;
                    break;
                case 'C' :
                    pwm_[i * 4 + 1] = log1;
                    break;
                case 'G' : 
                    pwm_[i * 4 + 2] = log1;
                    break;
                case 'T' :
                    pwm_[i * 4 + 3] = log1;
                    break;
            }
        }
    }

private:
    constexpr static double minconst = 0.001;
};

#endif /* PWM_H */
