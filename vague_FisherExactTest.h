#ifndef VAGUEFET_H
#define VAGUEFET_H

//get idea from fet_v7

#include <vector>

int getIndexOfMostFET(const std::vector<int>& vpos,
        const std::vector<int>& vneg, 
        int nF,
        int nB);

#endif /* VAGUEFET_H */
