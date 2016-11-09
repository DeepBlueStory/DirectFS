#ifndef MOTIFBASE_H
#define MOTIFBASE_H

#include <vector>
#include <string>
#include <unordered_map>

typedef struct Seed {
    std::string str;
    int nf;
    int nb;
    int fet;
}Seed;

class MotifBase {
protected:
    //length of seed
    const int m_lseed;

    //length of motif
    const int m_lmotif;

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

public:
    MotifBase(std::string fp, std::string fn, int ls, int lm);

    virtual ~MotifBase();

protected:
    //get the most seed
    //if successed, return true and put the seed into m_seeds
    //if not, return false
    
    bool getSeed();
};

#endif /* MOTIFBASE_H */
