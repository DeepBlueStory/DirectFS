#include "motifbase.h"
#include "motiftoolkit.h"

MotifBase::MotifBase(std::string fp, std::string fn, int ls, int lm):
    m_lseed(ls), m_lmotif(lm) {
    //readfasta
    m_seqsP = readFasta(fp);
    m_seqsN = readFasta(fn);

    //get <string, int> map
    m_mapPSeq = slideMapSeq(m_seqsP, ls);
    m_mapNSeq = slideMapSeq(m_seqsN, ls);
}

bool MotifBase::getSeed() {
    auto pairseed = findSeed(m_mapPSeq, m_mapNSeq);
    if(pairseed.first.first == "") return false;
    Seed seed;
    seed.str = pairseed.first.first;
    seed.nf = pairseed.second.first;
    seed.nb = pairseed.second.second;
    seed.fet = pairseed.first.second;
    m_seeds.push_back(seed);
    return true;
}
