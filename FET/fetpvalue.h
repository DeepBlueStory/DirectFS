//FileName: fetpvalue.h
//Date: 2017-01-02
//Author: Li Ning
//E-mail: jclining@126.com
//
//Brief: 
//  Calculate the p value by fisher exact test(FET).
//  
//There are 3 type of p values.
//1. Left: Use this when the alternative to independence is that there is negative association between the variables.
//   That is, the observations tend to lie in lower left and upper right.
//2. Right: Use this when the alternative to independence is that there is positive association between the variables.
//   That is, the observations tend to lie in upper left and lower right.
//3. 2-Tail: Use this when there is no prior alternative. /  
//
//Actually, almost all paper use 'right' p values not the 2-tail.
#ifndef FETPVALUE_H
#define FETPVALUE_H

#include <cmath>
#include <unordered_map>
#include <vector>
#include <iostream>

constexpr double LOG_ZERO = -1e10;
constexpr double LOG_SMALL = -0.5e10;
const double LOG0_99999999 = log(0.99999999);
const double LOG1_00000001 = log(1.00000001);
const double MM_NATS = log(-LOG_ZERO);
class FETPvalue
{
private:
    double my_exp(double x)
    {
        if(x < LOG_SMALL)
            return 0.0;
        else
            return exp(x);
    }
    double log_sum(double logx, double logy)
    {
        if(logy > logx)
        {
            double temp = logx;
            logx = logy;
            logy = temp;
        }
        if(logx - logy > MM_NATS)
            return logx;
        else
            return logx + log(1 + my_exp(logy - logx));
    }
    double lngamm(int z)
    {
        double x;
        x = 0.0;
        x = x + 0.1659470187408462e-06/(z+7.0);
        x = x + 0.9934937113930748e-05/(z+6.0);
        x = x - 0.1385710331296526    /(z+5.0);
        x = x + 12.50734324009056     /(z+4.0);
        x = x - 176.6150291498386     /(z+3.0);
        x = x + 771.3234287757674     /(z+2.0);
        x = x - 1259.139216722289     /(z+1.0);
        x = x + 676.5203681218835     /(z);
        x = x + 0.9999999999995183;
        return log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5);
    }
    double lnfact(int n)
    {
        if(n <= 1)
            return 0.0;
        if(lnfact_hash.find(n) != lnfact_hash.end())
            return lnfact_hash[n];
        double result = lngamm(n + 1);
        lnfact_hash[n] = result;
        return result;
    }
    double lnbico(int n, int k)
    {
        //log binomial coefficient n choose k
        return lnfact(n) - lnfact(k) - lnfact(n - k);
    }
    double log_hyper_323(int n11, int n1_, int n_1, int n)
    {
        return lnbico(n1_, n11) + lnbico(n-n1_, n_1-n11)-lnbico(n, n_1);
    }
    double log_hyper0(int n11i, int n1_i, int n_1i, int ni)
    {
        if((n1_i | n_1i | ni) == 0)
        {
            if(n11i % 10 != 0)
            {
                if(n11i == _sn11 + 1)
                {
                    _log_sprob += log(((_sn1_-_sn11) * 1.0 / n11i) * 
                            ((_sn_1 - _sn11) * 1.0 / (n11i + _sn - _sn1_ - _sn_1)));
                    _sn11 = n11i;
                    return _log_sprob;
                }
                if(n11i == _sn11 - 1)
                {
                    _log_sprob += log((_sn11 * 1.0 / (_sn1_ - n11i)) * 
                            ((_sn11 + _sn - _sn1_ - _sn_1) * 1.0 / (_sn_1 - n11i)));
                    _sn11 = n11i;
                    return _log_sprob;
                }
            }
            _sn11 = n11i;
        }
        else
        {
            _sn11 = n11i;
            _sn1_ = n1_i;
            _sn_1 = n_1i;
            _sn = ni;
        }
        _log_sprob = log_hyper_323(_sn11, _sn1_, _sn_1, _sn);
        return _log_sprob;
    }
    double log_hyper(int n11)
    {
        return log_hyper0(n11, 0, 0, 0);
    }
    double log_getFETprob(int a1, int a2, int b1, int b2)
    {
        int n = a1 + a2 + b1 + b2;
        int row1 = a1 + a2;
        int col1 = a1 + b1;
        int max = row1 < col1 ? row1 : col1;
        int min = row1 + col1 - n;
        if(min < 0)
            min = 0;
        if(min == max)
            //rt = (prob, sless, sright, sleft, slarg) = (1.0,1.0,1.0,1.0,1.0)
            //rt = (log_prob, log_sless, log_sright, log_sleft, log_slarg) = (0,0,0,0,0)
            return 0;
        double log_prob = log_hyper0(a1, row1, col1, n);
        double log_p;

        //double log_sleft = LOG_ZERO;
        //log_p = log_hyper(min);
        //int i = min + 1;
        //while(log_p < LOG0_99999999 + log_prob)
        //{
            //log_sleft = log_sum(log_sleft, log_p);
            //log_p = log_hyper(i);
            //++i;
        //}
        //--i;
        //if(log_p < LOG1_00000001 + log_prob)
            //log_sleft = log_sum(log_sleft, log_p);
        //else
            //--i;

        double log_sright = LOG_ZERO;
        log_p = log_hyper(max);
        int j = max - 1;
        while(log_p < LOG0_99999999 + log_prob)
        {
            log_sright = log_sum(log_sright, log_p);
            log_p = log_hyper(j);
            --j;
        }
        ++j;
        if(log_p < LOG1_00000001 + log_prob)
            log_sright = log_sum(log_sright, log_p);
        else
            ++j;
          
        //double log_sless;
        //double log_slarg;
        //if(abs(i - a1) < abs(j - a1))
        //{
            //log_sless = log_sleft;
            //log_slarg = log(1.0 - exp(log_sleft));
            //log_slarg = log_sum(log_slarg, log_prob);
        //}
        //else
        //{
            //log_sless = log(1.0 - exp(log_sright));
            //log_sless = log_sum(log_sless, log_prob);
            //log_slarg = log_sright;
        //}
        //return log_slarg;
        return log_sright;
    } 
    int getIndexOfThreshold(std::vector<int>& ocnPos, std::vector<int>& ocnNeg)
    {
        int intervalMax = 0;
        int index = 0;
        int size = (int) ocnPos.size();
        for(int i = 0; i < size; ++i)
        {
            int interval = ocnPos[i] - ocnNeg[i];
            if(interval < 0)
                continue;
            if(interval > intervalMax)
            {
                intervalMax = interval;
                index = i;
            }
                
        }
        return index;
    }

public:
    double getLogFETPvalue(int p, int P, int n, int N)
    {
        return log_getFETprob(N - n, n, P - p, p);
    }
    double getFETPvalue(int p, int P, int n, int N)
    {
        return exp(log_getFETprob(N - n, n, P - p, p));
    }
    int getIndexOfLeastPvalue(std::vector<int>& ocnPos, std::vector<int>& ocnNeg,
            int nPos, int nNeg)
    {
        int index = getIndexOfThreshold(ocnPos, ocnNeg);
        double threshold = getLogFETPvalue(ocnPos[index], nPos, ocnNeg[index], nNeg);

        int size = (int) ocnPos.size();
        int n = nPos + nNeg;
        for(int i = 0; i < size; ++i)
        {
            int row1 = nNeg;
            int col1 = n - ocnPos[i] - ocnNeg[i];
            int max = row1 < col1 ? row1 : col1;
            double log_prob = log_hyper0(nNeg - ocnNeg[i], row1, col1, n);
            double log_sright = LOG_ZERO;
            double log_p = log_hyper(max);
            int j = max - 1;
            while(log_p < LOG0_99999999 + log_prob)
            {
                log_sright = log_sum(log_sright, log_p);
                if(log_sright > threshold) break;
                log_p = log_hyper(j);
                --j;
            }
            if(log_p < LOG1_00000001 + log_prob)
                log_sright = log_sum(log_sright, log_p);
            if(log_sright < threshold)
            {
                threshold = log_sright;
                index = i;
            }
        }
        return index;
    }

private:
    static std::unordered_map<int, double> lnfact_hash;
    int _sn11, _sn1_, _sn_1, _sn;
    double _log_sprob;
};
std::unordered_map<int, double> FETPvalue::lnfact_hash;
#endif /* FETPVALUE_H */
