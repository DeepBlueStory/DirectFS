### FET 菲舍精确检验
菲舍精确检验可以写在一个外部文件中。
其中.h文件中就只有一个函数：
```
double FET(int nf, int nb, int nF, int nB);
```

### motif 父类
位置权重矩阵包里面有的函数应该有
```
//由正例负例集合求得map
std::unordered_map<std::string, int> slideMapSeq(const std::vector<std::string> seqs,
        const int length);

//找出seed
std::pair<string, std::pair<int, int>> findSeed(
        const std::unordered_map<std::string, int>& mapPSeq, 
        const std::unordered_map<std::string, int>& mapNSeq);
//返回值的first是seed字符串
//second中first是seed在正例中出现的次数，second是负例中出现的次数。
//mapPSeq为正例序列
//mapNSeq为负例序列

//将字符串seed转化为pwm矩阵
std::vector<double> seedToPWM(const std::string& seed, bool hasb = false);
//hasb代表是否需要在pwm之后加一位的b
//返回值采用move的形式

//xTwk表示向量相乘
//重载的两个函数分别表示不同情况。
//同时，如果wk的长度为奇数时，说明pwm后跟一位b。这时在最后加一位b。

//当x为字符串时，需要表明字符串的起始位置。默认为0
double xTwk(const std::string &S, const std::vector<double>& wk, const int istart = 0);

//当x为向量时，直接相乘。
double xTwk(const std::vector<bool>& S, const std::vector<double>& wk);

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


### slidewindow 滑动窗口
各种形式的滑动窗口

