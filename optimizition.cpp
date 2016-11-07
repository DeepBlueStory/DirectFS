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

PWM minFETSwb(PWM wk) {
    
}

PWM minFETSw(PWM wk, int i4l1) {
    int ic = i4l1 / 4;
    char cc = DNAalphabet[i4l1 % 4];
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

    std::map<double, int>::iterator inf = --nf.end();
    std::map<double, int>::iterator inb = --nb.end();
    if(inf->first > inb->first) {
        std::pair<double, int> p(inf->first, inb->second);
        nb.erase(inb);
        nb.insert(p);
    }else if(inf->first < inb->first) {
        std::pair<double, int> p(inb->first, inf->second);
        nf.erase(inf);
        nf.insert(p);
    }

    double minfetsw = maxdouble;
    int mf;
    int mb;
    double t;
    std::vector<double> w;
    while(inf != nf.end() && inb != nb.end()) {
        if(inf == nf.end()) --inf;
        if(inb == nb.end()) --inb;
        double fet = FEST(inf->second, inb->second);
        if(minfetsw > fet) {
            minfetsw = fet;
            mf = inf->second;
            mb = inb->second;
            t = inf->first < inb->first ? inf->first : inb->first;
        }
        
        if(inf->first > inb->first) ++inb;
        else if(inf->first < inb->first) ++inf;
        else {++inb; ++inf;}
    }
    for(auto i = wk.w.begin(); i != wk.w.end(); ++i) {
        w.push_back(*i);
    }
    w[i4l1] += t;
    
    PWM p;
    p.w = w;
    p.F = mf;
    p.B = mb;
    p.FETSw = minfetsw;

    return p;

}

void accadjMap(std::map<double, int> &nm, int reference) {
    auto itnm = std::lower_bound(nm.begin(), nm.end(), reference, 
            [](std::pair<double, int> a, std::pair<double, int> b) {
                return a.first < b.first;
            });
    std::map<double, int>::iterator right;
    std::map<double, int>::iterator left;
    left = itnm;
    right = ++itnm;
    int curr;

    curr = reference + left->second;
    for(auto i = left; i != nm.begin(); --i) {
        i->second = curr - i->second;
        curr -= left->second;
    }
    nm.begin()->second = curr - nm.begin()->second;

    curr = reference + left->second;
    for(auto i = right; i != nm.end(); ++i) {
        i->second = curr;
        curr += i->second;
    }
    //这里我为了对最后一个分段做处理，多加入一个te。
    //这个te的取值为最后一个te+1
    //如果太大，可以改
    //也有可能太小，需要取一个合适的值
    double te = (--nm.end())->first;
    te += 1;
    nm[te] = curr;
}

bool teOfS(std::string &S, PWM wk, double &te, int ic, char cc) {
    double tI = mindouble;
    double tO = mindouble;
    bool hasI = false;
    bool hasO = false;
    for(unsigned long istart = 0; istart <S.size()-m_lmotif+1; ++istart) {
        if(S[istart+ic] == cc) {
            hasI = true;
            tI = std::max(tI, xTwk(S, istart, wk, m_lmotif));
        }else {
            hasO = true;
            tO = std::max(tO, xTwk(S, istart, wk, m_lmotif));
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
