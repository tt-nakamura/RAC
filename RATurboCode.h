// uses NTL
//   http://www.shoup.net/ntl

#ifndef __RATurboCode_h__
#define __RATurboCode_h__

#include<NTL/vec_GF2.h>
using namespace NTL;

struct RATurboCode {
    long K,r,s;
    long rK,M,N;
    Vec<long> P;// permutation index
    RATurboCode(long, long=3, long=0);
    void encode_(Vec<GF2>&, const Vec<GF2>&);
    void encode(Vec<GF2>&, const Vec<GF2>&);
    long decode(Vec<GF2>&, const double*, long, double);
    long decode(Vec<GF2>&, const Vec<double>&, long=10, double=0);
};

void conv(Vec<GF2>&, const std::string&);
void conv(std::string&, const Vec<GF2>&);
void AddNoise(Vec<double>&, const Vec<GF2>&, double);

#endif // __RATurboCode_h__
