#include <vector>
#include <map>
#include <string>
#include <limits>
#include <algorithm>
typedef struct PWM {
    std::vector<double> w;
    int F;
    int B;
    double FETSw;
} PWM;

const double convergence = 0.001;
int m_lmotif = 8;

PWM optimizePWM(PWM w0) {
    PWM wk = w0;
    PWM wk1;
    do {
        PWM wk1i0 = wk;
        PWM wk1i;
        //4l
        //parallel for
        for(int i = 0; i < 4*m_lmotif; ++i) {
            wk1i = minFETSw(wk, i);
            if(wk1i0.FETSw > wk1i.FETSw) {
                wk1i0 = wk1i;
            }
        }
        //b
        wk1i = minFETSwb(wk);
        PWM wk1 = wk1i0.FETSw > wk1i.FETSw ? wk1i0 : wk1i;
    }while((wk1.FETSw-wk.FETSw) > convergence);
    return wk1;
}

std::vector<char> DNAalphabet = {
    'A', 'C', 'G', 'T'};
std::vector<std::string> Psequence;
std::vector<std::string> Nsequence;
PWM minFETSw(PWM wk, int i4l1) {
    int ic = i4l1 / 4;
    char cc = DNAalphabet[ic];

    std::map<double, int> nf;
    std::map<double, int> nb;
    for(auto ps : Psequence) {
        //
        //
    }
    for(auto ns : Nsequence) {
        //
        //
    }
}

const double mindouble = std::numeric_limits<double>::min();
bool teOfS(std::string &S, PWM wk, double &te, int ic, char cc) {
    double tI = mindouble;
    double tO = mindouble;
    bool hasI = false;
    bool hasO = false;
    for(int istart = 0; istart <S.size()-m_lmotif+1; ++istart) {
        if(S[istart+ic] == cc) {
            hasI = true;
            tI = max(tI, xTwk(S, istart, wk));
        }else {
            hasO = true;
            tO = max(tO, xTwk(S, istart, wk));
        }
    }
    if((hasI && !hasO) || ((hasI && hasO) && tO <=0)) {
        te = -tI;
        return true;
    }else {
        return false;
    }
}

double xTwk(std::string &S, int istart, PWM wk) {
    for(int i = istart; i < istart + m_lmotif; ++i) {

