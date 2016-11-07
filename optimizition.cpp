#include <vector>
#include <map>
#include <string>
#include <limits>
#include <algorithm>
#include "optimizition.h"


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

PWM minFETSw(PWM wk, int i4l1) {
    int ic = i4l1 / 4;
    char cc = DNAalphabet[ic];

    std::map<double, int> nf;
    std::map<double, int> nb;
    for(auto ps : Psequence) {
        double te;
        if(teOfS(ps, wk, te, ic, cc)) {
            if(nf.find(te) == nf.end()) {
                nf[te] = 1;
            }else {
                nf[te] ++;
            }
        }
    }
    for(auto ns : Nsequence) {
        double te;
        if(teOfS(ns, wk, te, ic, cc)) {
            if(nb.find(te) == nb.end()) {
                nb[te] = 1;
            }else {
                nb[te] ++;
            }
        }
    }
    accadjMap(nf, wk.F);
    accadjMap(nb, wk.B);

    std::map<double, int>::iterator inf = nf.begin();
    std::map<double, int>::iterator inb = nb.begin();

    while(inf != nf.end() && inb != nb.end()) {


    

    //调整
    //累加
    //求最大。

}

void accadjMap(std::map<double, int> &nm, int reference) {
    auto itnm = std::lower_bound(nm.begin(), nm.end(), reference, 
            [](std::pair<double, int> a, std::pair<double, int> b) {
                return a.first < b.first;
            });
    std::map<double, int>::iterator right;
    std::map<double, int>::iterator left;
    if(itnm->first == reference) {
        left = itnm;
        right = ++itnm;
    }else {
        right = itnm;
        left = --itnm;
    }
    int curr;
    curr = reference;
    for(auto i = right; i != nm.end(); ++i) {
        i->second += curr;
        curr = i->second;
    }
    curr = reference;
    for(auto i = left; i != nm.begin(); --i) {
        left->second = curr;
        curr -= left->second;
    }
    nm.begin()->second = curr;
}

bool teOfS(std::string &S, PWM wk, double &te, int ic, char cc) {
    double tI = mindouble;
    double tO = mindouble;
    bool hasI = false;
    bool hasO = false;
    for(unsigned long istart = 0; istart <S.size()-m_lmotif+1; ++istart) {
        if(S[istart+ic] == cc) {
            hasI = true;
            tI = max(tI, xTwk(S, istart, wk, m_lmotif));
        }else {
            hasO = true;
            tO = max(tO, xTwk(S, istart, wk, m_lmotif));
        }
    }
    if((hasI && !hasO) || ((hasI && hasO) && tO <=0)) {
        te = -tI;
        return true;
    }else {
        return false;
    }
}

double xTwk(std::string &S, int istart, PWM wk, int length) {
    double res = 0;
    auto w = wk.w;
    for(int i = 0; i <length; ++i) {
        char c = S[istart + i];
        switch(c) {
            case 'A' :
                res += w[4*i+0];
                break;
            case 'C' :
                res += w[4*i+1];
                break;
            case 'G' :
                res +=w[4*i+2];
                break;
            case 'T' :
                res += w[4*i+3];
                break;
        }
    }
    res += w[4*length];
    return res;
}
