#ifndef MOTIFBASE_H
#define MOTIFBASE_H

#include <vector>
#include <fstream>
#include "typedef.h"


class MotifBase {
public:
    MotifBase(std::string positiveFilePath, std::string negativeFilePath, int lenSeed, int lenMotif);

private:
    void readFasta(const std::string filePath, std::vector<sequence>& seqs);

protected:
    //length of seed
    const number lenSeed_;

    //length of motif
    const number lenMotif_;

    //sequences in foreground set;
    std::vector<sequence>  seqsInFg_;

    //sequences in background set;
    std::vector<sequence>  seqsInBg_;

    //number of sequences in foreground set;
    unsigned int nFg_;

    //number of sequences in background set;
    unsigned int nBg_;
};

MotifBase::MotifBase(std::string filePathFg, 
        std::string filePathBg, 
        int lenSeed, int lenMotif):
    lenSeed_(lenSeed), lenMotif_(lenMotif), 
    seqsInFg_(std::vector<sequence>()), 
    seqsInBg_(std::vector<sequence>()) {

    //read fasta file
    readFasta(filePathFg, seqsInFg_);
    readFasta(filePathBg, seqsInBg_);

    nFg_ = seqsInFg_.size();
    nBg_ = seqsInBg_.size();
}
    

void MotifBase::readFasta(const std::string filePath, std::vector<sequence>& seqs) {
    sequence str;
    std::ifstream infile;
    infile.open(filePath);
    if(infile) {
        while(getline(infile, str)) {
            //skip the blank line or the '>' line
            if(str.size() != 0 && str[0] != '>') {
                seqs.push_back(str);
            }
        }
        infile.close();
    }
}

#endif /* MOTIFBASE_H */
