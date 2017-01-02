# FET 菲舍精确检验

## FisherExactTest.h

Caculate the precise value of FisherExactTest
```
double FET(int ocnFg, int ocnBg, int nFg, int nBg);
```
## FisherExactTest_vague.h

Find the least FET from a vector of ocnFg and a vector of ocnBg
```
unsigned int FET_vague(const std::vector<int>& ocnFg, 
        const std::vector<int>& ocnBg, 
        int nFg, int nBg);
```


# Base

## motifbase.h
motif的父类，目的是管理一些基本资源。
其中存储着：
1. seed的长度
2. motif的长度
3. 前景序列
4. 背景序列
5. 前景序列的个数
6. 背景序列的个数

## pwm.h
关于pwm矩阵的向量式定义
其中存储着：
1. motif的pwm矩阵的向量形式
2. 在前景序列中出现的次数
3. 在背景序列中出现的次数
4. 精确的FET值

## seed.h
关于seed的定义
其中存储着
1. seed序列
2. 前景序列中出现的次数
3. 背景序列中出现的次数
4. 精确的FET值

# Toolkit
一些常用的工具


```
//由正例负例集合求得map
std::unordered_map<std::string, int> 
slideMapSeq(const std::vector<std::string> seqs, const int length);

//找出seed
Seed findSeed(const std::unordered_map<std::string, int>& mapPSeq, 
        const std::unordered_map<std::string, int>& mapNSeq);
//返回值的first是seed字符串
//second中first是seed在正例中出现的次数，second是负例中出现的次数。
//mapPSeq为正例序列
//mapNSeq为负例序列

//将字符串seed转化为pwm矩阵
std::vector<double> seedToPWM(const std::string& seed, bool hasb = false);
//hasb代表是否需要在pwm之后加一位的b
//返回值采用move的形式


//删除pwm
//既然要删除就一定要继续计算，所以mapPSeq是必须要的。
void deletePWMInPositiveSet(const std::vector<double>& w, 
        std::vector<std::string>& seqP, 
        std::unordered_map<std::string, int>& mapPSeq, 
        void (*pfprint) (const int iseq, const string motif, const bool isPosStrand) = NULL);
//w 是pwm矩阵
//seqP是正例序列
//mapPSeq是正例中string到int的map
//pfprint是输出函数， iseq是正例序列的序号，motif为被删除的序列，isPosStrand为是否为正链。

```
