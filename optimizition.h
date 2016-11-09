#ifndef DFSOPTIMIZITION_H
#define DFSOPTIMIZITION_H

#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <limits>

typedef struct Seed {
    std::string str;
    int nf;
    int nb;
    int fet;
}Seed;

typedef struct PWM {
    std::vector<double> w;  //the pwm vector which lenght is 4l+1
    int f;                  //the appearance time of w in positive set
    int b;                  //the appearance time of w in negative set
    double fetsw;           //the fet of w
} PWM;

const double MIN_DOUBLE = std::numeric_limits<double>::min();
const double MAX_DOUBLE = std::numeric_limits<double>::max();
const double CONVERGENCE = 0.001;
std::vector<char> DNAalphabet = {'A', 'C', 'G', 'T'};

//length of seed
int m_lseed;

//length of motif
int m_lmotif;

//positive set
std::vector<std::string> m_seqsP;

//negative set
std::vector<std::string> m_seqsN;

//map of positive set
std::unordered_map<std::string, int> m_mapPSeq;

//map of negative set
std::unordered_map<std::string, int> m_mapNSeq;

//all the seed
std::vector<Seed> m_seeds;

//when to stop?
double m_sshd;

//all the motif pwm
std::vector<PWM> m_pwms;
double optimizePWM();

PWM minFETSw(PWM wk, int i4l1);

void accadjMap(std::map<double, int> &nm, int reference);

bool teOfS(std::string &S, PWM wk, double &te, int ic, char cc);


#endif /* DFSOPTIMIZITION_H */
