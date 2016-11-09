#include "motiftoolkit.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "vague_FisherExactTest.h"
#include "FisherExactTest.h"


//read fasta file
//return the vector of sequence of positive set
//for c++11 just call this function with no need to consider the efficiency of copy
//
//filepath : the path of fasta file

std::vector<std::string> readFasta(const std::string& filepath) {
    std::vector<std::string> res;
    std::ifstream infile;
    infile.open(filepath);

    std::string sequence;
    if(infile) {
        while(getline(infile, sequence)) {
            //skip the blank line or the '>' line
            if(sequence.size() != 0 && sequence[0] != '>') {
                res.push_back(sequence);
            }
        }
        infile.close();
    }
    //for the vector, c++11 have make replace the copy as move.
    return res;
}


//get the reverse complementary DNA trand
//return the complementary trand string
//
//sequence : the positive DNA trand

std::string getReverseComplementaryTrand(const std::string& sequence) {
    std::string res;
    for(auto c : sequence) {
        switch(c) {
            case 'A' :
                res.push_back('T');
                break;
            case 'T' :
                res.push_back('A');
                break;
            case 'C' :
                res.push_back('G');
                break;
            case 'G' :
                res.push_back('C');
                break;
        }
    }
    std::reverse(res.begin(), res.end());
    return res;
}


//find all the string in seqs whoes size is length and its apparance time
//in both postive strand and supplement strand
//return the map from string to int
//for c++11, just call this function with no need to consider the efficiency of copy
//
//seqs : const reference of sequence set
//length : const length of string

std::unordered_map<std::string, int> 
slideMapSeq(const std::vector<std::string> seqs, const int length) {
    std::unordered_map<std::string, int> res;
    for(auto s : seqs) {
        //get the reverse complementary trand of s
        std::string rcoms = getReverseComplementaryTrand(s);
        for(unsigned long i = 0; i < s.size()-length+1; ++i) {
            std::string pos = s.substr(i, length);
            std::string rcs = rcoms.substr(i, length);
            if(res.find(pos) == res.end()) {
                res[pos] = 1;
            }else {
                res[pos] ++;
            }
            if(res.find(rcs) == res.end()) {
                res[rcs] = 1;
            } else {
                res[rcs] ++;
            }
        }
    }
    return res;
}


//find the most seed
//return the pair:
//  pair.first.first : seed(string)
//  pair.first.second : fet(double)
//  pair.second.first : nf(int)
//  pair.second.second : nb(int)
//
//mapPSeq : map of positive set <string, int>
//mapNSeq : map of negative set <string, int>

std::pair<std::pair<std::string, double>, std::pair<int, int>> 
findSeed(const std::unordered_map<std::string, int>& mapPSeq, 
        const std::unordered_map<std::string, int>& mapNSeq, 
        int nF, int nB) {
    //get string's apparance in both positive set and negative set which is a key of positive map
    std::vector<int> vpos;
    std::vector<int> vneg;
    for(auto i = mapPSeq.begin(); i != mapPSeq.end(); ++i) {
        vpos.push_back(i->second);
        if(mapNSeq.find(i->first) == mapNSeq.end()) {
            vneg.push_back(0);
        } else {
            vneg.push_back(mapNSeq.find(i->first)->second);
        }
    }
    auto ite = mapPSeq.begin();
    int index = getIndexOfMostFET(vpos, vneg, nF, nB);
    std::advance(ite, index);
    std::pair<std::string, double> pf(ite->first, 
            FET(vpos[index], vneg[index], nF, nB));
    std::pair<int, int> ps(vpos[index], vneg[index]);
    
    return 
    std::pair<std::pair<std::string, double>, std::pair<int, int>>(pf, ps);
}


//convert seed(string) to pwm(vector<double>)
//return the pwm of seed
//
//seed : const reference of seed.
//minconst : min const for log pwm
//b : const double value of pwm means the threshold of PWM

std::vector<double> seedToPWM(const std::string& seed, const double minconst, const double b) {
    std::vector<double> res;
    //according to the DNAalphabet
    for(auto c : seed) {
        switch(c) {
            case 'A' :
                res.push_back(1+minconst);
                res.push_back(0+minconst);
                res.push_back(0+minconst);
                res.push_back(0+minconst);
                break;
            case 'C' :
                res.push_back(0+minconst);
                res.push_back(1+minconst);
                res.push_back(0+minconst);
                res.push_back(0+minconst);
                break;
            case 'G' :
                res.push_back(0+minconst);
                res.push_back(0+minconst);
                res.push_back(1+minconst);
                res.push_back(0+minconst);
                break;
            case 'T' :
                res.push_back(0+minconst);
                res.push_back(0+minconst);
                res.push_back(0+minconst);
                res.push_back(1+minconst);
                break;
        }
    }
    
    //b != 0 means need to add -b after pwm
    if(b != 0) {
        res.push_back(-b);
    }

    //log pwm
    for(auto& d : res) {
        d = std::log(d);
    }

    return res;
}


//get the result of x(string) multiplied w(vector)
//return the result 
//
//S : const reference of x
//wk : const reference of w
//istart: const value indicate the begin index of x in S, default is 0

double xTwk(const std::string &S, 
        const std::vector<double>& wk, 
        const int istart) {
    double res = 0;
    //the length of pwm
    int length = wk.size() / 4;
    //where there has a b
    bool hasb = wk.size() % 4;
    for(int i = 0; i < length; ++i) {
        switch(S[istart+i]) {
            case 'A' :
                res += wk[4*i+0];
                break;
            case 'C' :
                res += wk[4*i+1];
                break;
            case 'G' :
                res += wk[4*i+2];
                break;
            case 'T' :
                res += wk[4*i+3];
                break;
        }
    }
    
    if(hasb) res += wk.back();
    return res;
}


//delete motif in every positive sequence that have highest macthing score of w.
//return void 
//
//w : const reference of w
//seqP : reference of positive sequence set
//mapPSeq: reference of map of positive set, that need to get change for next findseed.
//pfprint : output function for print the delete information of motif.
//  iseq : const value of index of motif in positive set
//  motif : const string of motif to be delete
//  isPosStrand : motif in the positive strand or supplementary strand

void deletePWMInPositiveSet(const std::vector<double>& w, 
        std::vector<std::string>& seqP, 
        std::unordered_map<std::string, int>& mapPSeq, 
        double shreshold,
        void (*pfprint) (const int iseq, 
            const std::string motif, 
            const bool isPosStrand)) {
    //when delete, I need to think the reverse complementary strand
}
