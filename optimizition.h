#ifndef OPTIMIZITION_H
#define OPTIMIZITION_H

#include <map>
#include <vector>
#include <string>
#include <limits>
typedef struct PWM {
    std::vector<double> w;
    int F;
    int B;
    double FETSw;
} PWM;

const double mindouble = std::numeric_limits<double>::min();
const double maxdouble = std::numeric_limits<double>::max();
const double convergence = 0.001;
int m_lmotif = 8;
std::vector<char> DNAalphabet = {
    'A', 'C', 'G', 'T'};
std::vector<std::string> Psequence;
std::vector<std::string> Nsequence;

PWM optimizePWM(PWM w0);

PWM minFETSw(PWM wk, int i4l1);

void accadjMap(std::map<double, int> &nm, int reference);
bool teOfS(std::string &S, PWM wk, double &te, int ic, char cc);
double xTwk(std::string &S, int istart, PWM wk, int length);
#endif //OPTIMIZITION_H
