#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include "optimizition.h"
#include "motiftoolkit.h"

bool hasTeIn4l(const std::string &S, const std::vector<double>& wk, int i4l1, double &te);
void accadjMap(std::map<double, int> &nm, int reference);
double FET(int mf, int mb, int nf, int nb);
PWM minFETSwIn4l(const PWM wk, int i4l);
PWM minFETSwb(PWM wk);
double teForB(const std::string& S, const std::vector<double>& wk);

double optimizePWM(PWM& w0){
    PWM wk = w0;
    PWM wk1 = w0;
    do {
        wk = wk1;
        //wk1i original value is w0
        PWM wk1i0 = w0;
        PWM wk1i;
        //4l
        //parallel for
        for(unsigned long int i = 0; i < w0.w.size()-1; ++i) {
            wk1i = minFETSwIn4l(wk, i);
            if(wk1i0.fetsw > wk1i.fetsw) {
                wk1i0 = wk1i;
            }
        }
        //b
        wk1i = minFETSwb(wk);
        wk1 = wk1i0.fetsw < wk1i.fetsw ? wk1i0 : wk1i;
    }while((wk1.fetsw-wk.fetsw) > CONVERGENCE);
    //optimal pwm in the last of pwms vector
    w0 = wk1;
    return w0.fetsw;
}

PWM minFETSwb(PWM wk) {
    //for optimizing b, every x's (4l+1)th is 1
    //so every S in foreground and background have to calculate the te.
    
    std::map<double, int> nf;
    std::map<double, int> nb;
    double te;

    for(auto ps : m_seqsP) {
        te = teForB(ps, wk.w);
        if(nf.find(te) == nf.end()) {
            nf[te] = 1;
        }else {
            nf[te] ++;
        }
    }
    for(auto ns : m_seqsN) {
        te = teForB(ns, wk.w);
        if(nb.find(te) == nb.end()) {
            nb[te] = 1;
        }else {
            nb[te] ++;
        }
    }
    accadjMap(nf, wk.f);
    accadjMap(nb, wk.b);

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

    PWM p;
    p.fetsw = MAX_DOUBLE;
    double t;
    inf = nf.begin();
    inb = nb.begin();
    while(inf != nf.end() && inb != nb.end()) {
        if(inf == nf.end()) --inf;
        if(inb == nb.end()) --inb;
        double fet = FET(inf->second, inb->second, m_seqsP.size(), m_seqsN.size());
        if(p.fetsw > fet) {
            p.fetsw = fet;
            p.f= inf->second;
            p.b= inb->second;
            t = inf->first < inb->first ? inf->first : inb->first;
        }
        
        if(inf->first > inb->first) ++inb;
        else if(inf->first < inb->first) ++inf;
        else {++inb; ++inf;}
    }
    for(auto i = wk.w.begin(); i != wk.w.end(); ++i) {
        p.w.push_back(*i);
    }
    p.w.back() += t;
    
    return p;
}


//moveto
PWM minFETSwIn4l(const PWM wk, int i4l) {
    std::map<double, int> nf;
    std::map<double, int> nb;
    double te;
    for(auto ps : m_seqsP) {
        if(hasTeIn4l(ps, wk.w, i4l, te)) {
            if(nf.find(te) == nf.end()) {
                nf[te] = 1;
            }else {
                nf[te] ++;
            }
        }
    }
    for(auto ns : m_seqsN) {
        if(hasTeIn4l(ns, wk.w, i4l, te)) {
            if(nb.find(te) == nb.end()) {
                nb[te] = 1;
            }else {
                nb[te] ++;
            }
        }
    }
    accadjMap(nf, wk.f);
    accadjMap(nb, wk.b);

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

    PWM p;
    p.fetsw = MAX_DOUBLE;
    double t;
    inf = nf.begin();
    inb = nb.begin();
    while(inf != nf.end() && inb != nb.end()) {
        if(inf == nf.end()) --inf;
        if(inb == nb.end()) --inb;
        double fet = FET(inf->second, inb->second, m_seqsP.size(), m_seqsN.size());
        if(p.fetsw > fet) {
            p.fetsw = fet;
            p.f= inf->second;
            p.b= inb->second;
            t = inf->first < inb->first ? inf->first : inb->first;
        }
        
        if(inf->first > inb->first) ++inb;
        else if(inf->first < inb->first) ++inf;
        else {++inb; ++inf;}
    }
    for(auto i = wk.w.begin(); i != wk.w.end(); ++i) {
        p.w.push_back(*i);
    }
    p.w[i4l] += t;
    
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

bool hasTeIn4l(const std::string &S, const std::vector<double>& wk, int i4l, double &te){
    int ic = i4l / 4;
    char cc = DNAalphabet[i4l % 4];
    double tI = MIN_DOUBLE;
    double tO = MIN_DOUBLE;
    bool hasI = false;
    bool hasO = false;
    for(unsigned long istart = 0; istart <S.size()-(wk.size()-1)/4+1; ++istart) {
        if(S[istart+ic] == cc) {
            hasI = true;
            tI = std::max(tI, xTwk(S, wk, istart));
        }else {
            hasO = true;
            tO = std::max(tO, xTwk(S, wk, istart));
        }
    }
    if((hasI && !hasO) || ((hasI && hasO) && tO <=0)) {
        te = -tI;
        return true;
    }else {
        return false;
    }
}

double teForB(const std::string& S, const std::vector<double>& wk) {
    double tI = MIN_DOUBLE;
    for(unsigned long istart = 0; istart < S.size()-(wk.size()-1)/4+1; ++istart) {
        tI = std::max(tI, xTwk(S, wk, istart));
    }
    return -tI;
}

void extendPWM() {
    PWM& wk = m_pwms.back();
    PWM wkl, wkr;
    do {
        wkl.f = wk.f;
        wkl.b = wk.b;
        wkl.fetsw = wk.fetsw;
        for(int i = 0; i < 4; ++i) wkl.w.push_back(0);
        for(auto i = wk.w.begin(); i != wk.w.end(); ++i) 
            wkl.w.push_back(*i);
        
        wkr.f = wk.f;
        wkr.b = wk.b;
        wkr.fetsw = wk.fetsw;
        for(auto i = wk.w.begin(); i != wk.w.end(); ++i) 
            wkr.w.push_back(*i);
        double te = wkr.w.back();
        for(int i = 0; i < 4; ++i) wkr.w.push_back(0);
        wkr.w.push_back(te);

        wk = optimizePWM(wkl) < optimizePWM(wkr) ? wkl : wkr;
    }while(wk.w.size() < (unsigned long)m_lmotif);
}
