#ifndef DIRECTFS_H
#define DIRECTFS_H

#include <vector>
#include <limits>
#include "motifbase.h"

const double MIN_DOUBLE = std::numeric_limits<double>::min();
const double MAX_DOUBLE = std::numeric_limits<double>::max();
const std::vector<char> DNAalphabet = {'A', 'C', 'G', 'T'};

typedef struct PWM {
    std::vector<double> w;  //the pwm vector which lenght is 4l+1
    int f;                  //the appearance time of w in positive set
    int b;                  //the appearance time of w in negative set
    double fetsw;           //the fet of w
} PWM;

class DirectFS : private MotifBase
{
private:
    //when to stop?
    double m_sshd;

    //all the motif pwm
    std::vector<PWM> m_pwms;

    static constexpr double CONVERGENCE = 0.001;

public:
    //read the positive & negative set.
    //
    //ls means seed length
    //lm means motif length
    //fp means filename of positive set
    //fn means filename of negative set
    //shreshold_b means the shreshold that can't make a pwm to be a motif
    DirectFS(std::string fp, std::string fn, int ls, int lm, double shreshold);
    
    virtual ~DirectFS();

    //find all the pwm
    void findAllPwm();

private:
    //optimize pwm
    //return the fet of pwm.
    double optimizePWM();

    //extend pwm from length(seed) to length(motif)
    void extendPWM();
 

};

#endif /* DIRECTFS_H */
