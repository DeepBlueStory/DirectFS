#include <iostream>
#include "directFS.h"
#include "motiftoolkit.h"

DirectFS::DirectFS(std::string fp, std::string fn, int ls, int lm, double shreshold):
    MotifBase(fp, fn, ls, lm),
    m_sshd(shreshold) {}

void DirectFS::findAllPwm() {
    while(true) {
        if(!getSeed()) break;

        //convert the seed to pwm
        auto w = seedToPWM(m_seeds.back().str, true);
        PWM pwm;
        pwm.w = w;
        pwm.f = m_seeds.back().nf;
        pwm.b = m_seeds.back().nb;
        pwm.fetsw = m_seeds.back().fet;
        m_pwms.push_back(pwm);

        //optimize the last pwm in m_pwms
        double fet = optimizePWM();

        if(fet < m_sshd) {
            m_pwms.pop_back();
            break;
        }
        
        if(m_lseed < m_lmotif) {
            extendPWM();
        }
        
        //delete the motif in input positive set
        deletePWMInPositiveSet(m_pwms.back().w, m_seqsP, m_mapPSeq, []
                (const int iseq, const std::string motif, const bool isPosStrand) {
                std::cout << iseq << '\t';
                std::cout << motif << '\t';
                std::cout << isPosStrand << std::endl;
                });
    }
}

double DirectFS::optimizePWM() {
    return 0;
}

void DirectFS::extendPWM() {
}
