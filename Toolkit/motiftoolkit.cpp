#include <vector>
#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include "../Base/seed.h"
#include "../FET/FisherExactTest.h"
#include "../Base/pwm.h"

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
//in both strand and supplement strand
//return the map from string to int
//for c++11, just call this function with no need to consider the efficiency of copy
//
//seqs : const reference of sequence set
//length : const length of string

std::unordered_map<std::string, int>  
slideMapSeq(const std::vector<std::string>& seqs, const int length) {
    std::unordered_map<std::string, int> ret;
    for(auto s : seqs) {
        //get the reverse complementary trand of s
        std::string rcoms = getReverseComplementaryTrand(s);
        for(unsigned long i = 0; i < s.size()-length+1; ++i) {
            std::string pos = s.substr(i, length);
            std::string rcs = rcoms.substr(i, length);
            if(ret.find(pos) == ret.end()) {
                ret[pos] = 1;
            }else {
                ret[pos] ++;
            }
            if(ret.find(rcs) == ret.end()) {
                ret[rcs] = 1;
            } else {
                ret[rcs] ++;
            }
        }
    }
    return ret;
}


//find the most seed
//
//mapPSeq : map of positive set <string, int>
//mapNSeq : map of negative set <string, int>

Seed findSeed(const std::unordered_map<std::string, int>& mapPSeq, 
                const std::unordered_map<std::string, int>& mapNSeq, 
                int nFg, int nBg)
{
    //get string's apparance in both positive set and negative set which is a key of positive map
    std::vector<int> vpos;
    std::vector<int> vneg;
    for(auto i = mapPSeq.begin(); i != mapPSeq.end(); ++i) {
        vpos.push_back(i->second);
        if(mapNSeq.find(i->first) == mapNSeq.end()) {
            vneg.push_back(0);
        } else {
            vneg.push_back((mapNSeq.find(i->first))->second);
        }
    }
    FETOneTail fet;
    auto ite = mapPSeq.begin();
    int index = fet.getIndexOfLeastPvalue(vpos, vneg, nFg, nBg);
    int ocnFg= vpos[index];
    int ocnBg= vneg[index];
    std::advance(ite, index);

    return Seed(ite->first, fet.getLogPvalue(ocnFg, ocnBg, nFg, nBg), ocnFg,ocnBg);
}

double wTx(const PWM& wk, const std::string &seq, const int iBegin)
{
    double ret = 0;
    auto& pwm = wk.pwm_;
    int length = pwm.size() / 4;
    for(int i = 0; i < length; ++i)
    {
        switch(seq[iBegin+i])
        {
            case 'A' :
                ret += pwm[4*i+0];
                break;
            case 'C' :
                ret += pwm[4*i+1];
                break;
            case 'G' :
                ret += pwm[4*i+2];
                break;
            case 'T' :
                ret += pwm[4*i+3];
                break;
        }
    }
    ret += wk.b_;
    return ret;
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

void deletePWMInSeqs(PWM wk, std::vector<std::string>& seqs, 
        std::unordered_map<std::string, int>& numLmer, double shreshold)
{
    //when delete, I need to think the reverse complementary strand
    wk.b_ = 0;
    int lenMotif = wk.pwm_.size() / 4;
    //parallel for
    for(int i = 0; i < (int) seqs.size(); i++)
    {
        std::string& exactTrand = seqs[i];
        double minFET = 1.0;
        int minIndex = -1;
        auto complTrand = getReverseComplementaryTrand(exactTrand);
        for(int j = 0; j < (int) exactTrand.size(); ++j)
        {
            double tempExactFET = wTx(wk, exactTrand, j);
            double tempComplFET = wTx(wk, complTrand, j);
            if(tempExactFET < minFET)
            {
                minFET = tempExactFET;
                minIndex = j;
            }
            if(tempComplFET < minFET)
            {
                minFET = tempComplFET;
                minIndex = j;
            }
        }
        if(minFET < shreshold)
        {
            int pos = minIndex - lenMotif / 2;
            int len = lenMotif % 2 == 0 ? 2 * lenMotif : 2 * lenMotif - 1;
            len = pos > 0 ? len : len + pos;
            pos = pos > 0 ? pos : 0;
            std::string tempDelete(exactTrand, pos, len);

            exactTrand.erase(pos, len);
            for(int k = 0; k < len - lenMotif + 1; ++k)
            {
                std::string tempLmer = tempDelete.substr(k, lenMotif);
                numLmer[tempLmer]--;
            }

            std::cout << "delete in " << i << "sequence"<< std::endl;
            std::cout << '\t' << "position is " << minIndex << std::endl;
            std::cout << '\t' << "motif is " << exactTrand.substr(minIndex, lenMotif) << std::endl;
        }
    }
}
