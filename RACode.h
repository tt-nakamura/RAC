// uses NTL
//   http://www.shoup.net/ntl

#ifndef __RACode_h__
#define __RACode_h__

#include<NTL/vec_GF2.h>
using namespace NTL;

struct RACode {
    long K,r,s;
    long rK,M,N;
    Vec<long> P;// permutation index
    Vec<Vec<long> > I,J;// check matrix (sparse representation)
    RACode(long, long=3, long=0);
    void encode_(Vec<GF2>&, const Vec<GF2>&);
    void encode(Vec<GF2>&, const Vec<GF2>&);
    long decode(Vec<GF2>&, const double*, long);
    long decode(Vec<GF2>&, const Vec<double>&, long=32);
};

void conv(Vec<GF2>&, const std::string&);
void conv(std::string&, const Vec<GF2>&);
void AddNoise(Vec<double>&, const Vec<GF2>&, double);

#endif // __RACode_h__
