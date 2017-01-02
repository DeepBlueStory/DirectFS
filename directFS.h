#ifndef DIRECTFS_H
#define DIRECTFS_H

#include <map>
#include <vector>
#include <unordered_map>
#include "./Base/motifbase.h"
#include "./Base/pwm.h"
#include "./Base/seed.h"

class DirectFS : private MotifBase
{

public:
    DirectFS(std::string filePathPos, std::string filePathNeg, 
            int lenSeed, int lenMotif, double b);
    
    void findAllPwm();

private:
    double optimizePWM(PWM &pwm);

    double extendPWM(PWM &pwm);

    PWM initANewPwm();

    int whichTypeOfSeq(const PWM& wk, const std::string& seq, const int ei,
            double& maxIs, double& maxOs);

    int whichScenarioOfSeq(const int& typeOfSeq, const double& maxOs);

    void getAllBreakPoint(std::map<double, int>& breakPointFg, 
            std::map<double, int>& breakPointBg, 
            const PWM& wk, const int& ei);

    double getOptimalStepSize(const PWM& wk, const int ei, int& ocnFg, int& ocnBg);
    
    double getMinFETofIDim(PWM& wk, int i);

    void iteratedOptimition(PWM& wk);

private:
    double b_;
    //all the motif pwm
    std::vector<PWM> pwms_;

    //all the lmers and its occurence number in both foreground and background
    std::unordered_map<std::string, int> lmerInFg_;
    std::unordered_map<std::string, int> lmerInBg_;
};


#endif /* DIRECTFS_H */
