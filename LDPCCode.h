// uses NTL
//   http://www.shoup.net/ntl

#ifndef __LDPCCode_h__
#define __LDPCCode_h__

#include<NTL/mat_GF2.h>
using namespace NTL;

struct LDPCCode {// Low Density Parity Check Code
    long N;// code length
    long K;// message length
    Mat<GF2> H;// parity check matrix (transposed)
    Mat<GF2> D;// some part of H
    Vec<Vec<long> > I,J;// sparse representation of H
    Vec<Vec<long> > A,B,C,E,T;// some part of I
    LDPCCode(long, long, long=3);
    void encode_(Vec<GF2>&, const Vec<GF2>&);
    void encode(Vec<GF2>&, const Vec<GF2>&);
    long decode_(Vec<GF2>&, const double[], long);
    long decode(Vec<GF2>&, const Vec<double>&, long=20);
};

void indexr(const int[], int[], int);
void conv(Vec<GF2>&, const std::string&);
void conv(std::string&, const Vec<GF2>&);
void AddNoise(Vec<double>&, const Vec<GF2>&, double);

#endif // __LDPCCode_h__