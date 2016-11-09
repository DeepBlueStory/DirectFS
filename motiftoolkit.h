#ifndef MOTIF_H
#define MOTIF_H

#include <vector>
#include <string>
#include <unordered_map>

//here to mark the encoding method
//A 1000
//C 0100
//G 0010
//T 0001

std::vector<char> DNAalphabet = {'A', 'C', 'G', 'T'};


//read fasta file
//return the vector of sequence of positive set
//for c++11 just call this function with no need to consider the efficiency of copy
//
//filepath : the path of fasta file

std::vector<std::string> readFasta(const std::string& filepath);


//get the reverse complementary DNA trand
//return the complementary trand string
//
//sequence : the positive DNA trand

std::string getReverseComplementaryTrand(const std::string& sequence);


//find all the string in seqs whoes size is length
//return the map from string to int
//
//seqs : const reference of sequence set
//length : const length of string

std::unordered_map<std::string, int> 
slideMapSeq(const std::vector<std::string> seqs, const int length);


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
        const std::unordered_map<std::string, int>& mapNSeq);


//convert seed(string) to pwm(vector<double>)
//return the pwm of seed
//
//seed : const reference of seed.
//minconst : min const for log pwm
//b : const double value of pwm means the threshold of PWM

std::vector<double> seedToPWM(const std::string& seed, const double minconst, const double b);


//get the result of x(string) multiplied w(vector)
//return the result 
//
//S : const reference of x
//wk : const reference of w
//istart: const value indicate the begin index of x in S, default is 0

double xTwk(const std::string &S, 
        const std::vector<double>& wk, 
        const int istart = 0);


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
            const bool isPosStrand) = NULL);

#endif /* MOTIF_H  */
