#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <cfloat>
#include <algorithm>
#include "./Base/pwm.h"
#include "./FET/FisherExactTest.h"
#include "directFS.h"
#include "Toolkit/motiftoolkit.cpp"
#include "./Base/seed.h"

const double CONVERGENCE = 0.001;
const std::vector<char> DNAalphabet = {'A', 'C', 'G', 'T'};
const double INTERVALSIZE = 0.1;

DirectFS::DirectFS(std::string filePathPos, std::string filePathNeg, 
        int lenSeed, int lenMotif, double b)
    : MotifBase(filePathPos, filePathNeg, lenSeed, lenMotif), b_(b)
{
    lmerInFg_ = slideMapSeq(seqsInFg_, lenSeed);
    lmerInBg_ = slideMapSeq(seqsInBg_, lenSeed);
}

PWM DirectFS::initANewPwm()
{
    Seed seed = findSeed(lmerInFg_, lmerInBg_, nFg_, nBg_);
    PWM pwm(seed, b_);
    return pwm;
}

void DirectFS::findAllPwm()
{
    while(true)
    {
        PWM pwm = initANewPwm();
        //stop condition 1
        //if(pwm.FET_ < 0.05) break;
        
        //stop condition 2
        //if(optimizePWM(pwm) > shreshold) break;

        if(lenSeed_ < lenMotif_)
            extendPWM(pwm);

        pwms_.push_back(pwm);

        //deletePWMInPositiveSet;
    }
}
int DirectFS::whichTypeOfSeq(const PWM& wk, const std::string& seq, const int ei,
        double& maxIs, double& maxOs)
{
    if(-1 == ei)
    {
        double tempbMaxIs = -DBL_MAX;
        for(int i = 0; i < (int) seq.size(); ++i)
        {
            tempbMaxIs = wTx(wk, seq, i);
            maxIs = tempbMaxIs > maxIs ? tempbMaxIs : maxIs;
        }
        return 3;
    }
    bool flagIs = false;
    bool flagOs = false;
    maxIs = -DBL_MAX;
    maxOs = -DBL_MAX;
    //the position and the character
    int eiPosition = ei % 4;
    char eiAlphabet = DNAalphabet[eiPosition];

    double tempMaxIs;
    double tempMaxOs;
    for(int i = 0; i < (int) seq.size(); ++i)
    {
        if(seq[i + eiPosition] == eiAlphabet)
        {
            flagIs = true;
            tempMaxIs = wTx(wk, seq, i);
            maxIs = tempMaxIs > maxIs ? tempMaxIs : maxIs;
        }
        else
        {
            flagOs = true;
            tempMaxOs = wTx(wk, seq, i);
            maxOs = tempMaxOs > maxOs ? tempMaxOs : maxOs;
        }
    }
    
    if(flagIs && flagOs)
        return 1;
    if(flagOs)
        return 2;
    else
        return 3;
}
int DirectFS::whichScenarioOfSeq(const int& typeOfSeq, const double& maxOs)
{
    if((1 == typeOfSeq && 0 >= maxOs) || 3 == typeOfSeq)
    {
        return 1;
    }
    if((1 == typeOfSeq || 2 == typeOfSeq) && 0 < maxOs)
    {
        return 2;
    }
    if(2 == typeOfSeq && 0 >= maxOs)
    {
        return 3;
    }
    std::cout << "error in scenario jugement" << std::endl;
    return -1;
}
void DirectFS::getAllBreakPoint(std::map<double, int>& breakPointFg, 
        std::map<double, int>& breakPointBg, 
        const PWM& wk, const int& ei)
{
    //Foreground
    //parallel for
    for(auto& seq : seqsInFg_)
    {
        double tempMaxIs = -DBL_MAX;
        double tempMaxOs = -DBL_MAX;
        int tempType = whichTypeOfSeq(wk, seq, ei, tempMaxIs, tempMaxOs);
        int tempScenario = whichScenarioOfSeq(tempType, tempMaxOs);
        if(1 == tempScenario)
        {
            double ts = -tempMaxIs;
            if(breakPointFg.find(ts) == breakPointFg.end())
            {
                breakPointFg[ts] = 1;
                breakPointBg[ts] = 0;
            }
            else
            {
                breakPointFg[ts] ++;
            }
        }
    }
    //Background
    //parallel for
    for(auto& seq : seqsInBg_)
    {
        double tempMaxIs = -DBL_MAX;
        double tempMaxOs = -DBL_MAX;
        int tempType = whichTypeOfSeq(wk, seq, ei, tempMaxIs, tempMaxOs);
        int tempScenario = whichScenarioOfSeq(tempType, tempMaxOs);
        if(1 == tempScenario)
        {
            double ts = -tempMaxIs;
            if(breakPointFg.find(ts) == breakPointFg.end())
            {
                breakPointBg[ts] = 1;
                breakPointFg[ts] = 0;
            }
            else
            {
                breakPointBg[ts] ++;
            }
        }
    }
}
double DirectFS::getOptimalStepSize(const PWM& wk, const int ei, int& ocnFg, int& ocnBg)
{
    std::map<double, int> breakPointFg;
    std::map<double, int> breakPointBg;
    getAllBreakPoint(breakPointFg, breakPointBg, wk, ei);

    size_t tempSizeOfInterval = breakPointFg.size() + 1;
    std::vector<double> interval(tempSizeOfInterval);
    std::vector<int> numMotifInInterFg(tempSizeOfInterval);
    std::vector<int> numMotifInInterBg(tempSizeOfInterval);

    auto iFg = breakPointFg.begin();
    double tlFg = iFg->first - INTERVALSIZE;
    double trFg = iFg->first;
    int diffBetweenFg = 0;
    int accFg = 0;

    auto iBg = breakPointBg.begin();
    double tlBg = iBg->first - INTERVALSIZE;
    double trBg = iBg->first;
    int diffBetweenBg = 0;
    int accBg = 0;

    int indexInterval = 0;

    while(iFg != breakPointFg.end() || iBg != breakPointBg.end())
    {
        trFg = iFg->first;
        trBg = iBg->first;
        if(0 <= trFg && tlFg < 0)
        {
            diffBetweenFg = accFg;
        }
        if( 0 <= trBg && tlBg < 0)
        {
            diffBetweenBg = accBg;
        }
        numMotifInInterFg[indexInterval] = accFg;
        numMotifInInterBg[indexInterval] = accBg;
        interval[indexInterval] = (tlFg + trFg) / 2;
        accFg += iFg->second;
        accBg += iBg->second;
        
        tlFg = trFg;
        tlBg = trBg;
        iFg++;
        iBg++;
        indexInterval++;
    }
    if(iFg != breakPointFg.end() || iBg != breakPointBg.end())
    {
        std::cout << "intervals in foreground and background not equal!" << std::endl;
        exit(0);
    }

    numMotifInInterFg[indexInterval] = accFg;
    numMotifInInterBg[indexInterval] = accBg;
    interval[indexInterval] = tlFg + INTERVALSIZE / 2;

    int optimalIndex = getIndexOfLeastPvalue(numMotifInInterFg, numMotifInInterBg, nFg_, nBg_);
    ocnFg = numMotifInInterFg[optimalIndex];
    ocnBg = numMotifInInterBg[optimalIndex];
    return interval[optimalIndex];
}
        
double DirectFS::getMinFETofIDim(PWM& wk, int i)
{
    //modify the wk, let wk become the best pwm at ei == i;
    int ocnFg = 0;
    int ocnBg = 0;
    double t = getOptimalStepSize(wk, i, ocnFg, ocnBg);
    wk.ocnFg_ = ocnFg;
    wk.ocnBg_ = ocnBg;
    if(-1 == i)
    {
        wk.b_ += t;
    }
    else
    {
        wk.pwm_[i] += t;
    }
    wk.FET_ = getLogPvalue(ocnFg, ocnBg, nFg_, nBg_);
    return wk.FET_;
}

void DirectFS::iteratedOptimition(PWM& wk)
{
    //pseudo code
    PWM wk_1;
    PWM tempWk_1;
    
    //parallel for
    int tempLenWk = wk.pwm_.size();
    for(int i = 0; i < tempLenWk; ++i)
    {
        tempWk_1 = wk;
        if(getMinFETofIDim(tempWk_1, i) < wk_1.FET_)
        {
            wk_1 = tempWk_1;
        }
    }
    
    //optimize b
    tempWk_1 = wk;
    if(getMinFETofIDim(tempWk_1, -1) < wk_1.FET_)
    {
        wk_1 = tempWk_1;
    }
    wk = wk_1;
}

double DirectFS::optimizePWM(PWM& wk)
{
    double tempFET0 = wk.FET_;
    double tempFET1 = 0;
    while(CONVERGENCE > tempFET0 - tempFET1)
    {
        iteratedOptimition(wk);
        tempFET1 = tempFET0;
        tempFET0 = wk.FET_;
    }
    return tempFET0;
}

double DirectFS::extendPWM(PWM& wk)
{
    int extendLen = lenMotif_ - lenSeed_;
    PWM wkLeft;
    PWM wkRight;
    for(int i = 0; i < extendLen; ++i)
    {
        wkLeft = wk;
        for(int j = 0; j < 4; ++j)
            wkLeft.pwm_.push_front(0.0);
        wkRight = wk;
        for(int j = 0; j < 4; ++j)
            wkRight.pwm_.push_back(0.0);

        if(optimizePWM(wkLeft) < optimizePWM(wkRight))
            wk = wkLeft;
        else
            wk = wkRight;
    }
    return wk.FET_;
}
