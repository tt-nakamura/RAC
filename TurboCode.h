// uses NTL
//   http://www.shoup.net/ntl

#ifndef __TurboCode_h__
#define __TurboCode_h__

#include<NTL/GF2X.h>
#include<NTL/mat_GF2.h>
using namespace NTL;

struct TurboCode {
    static const long coeff[2];
    long m;// memory order
    long K;// message length
    long M;// number of transitions = 1<<(m+1)
    GF2X G[2];// denominator and numerator of generator function
    Mat<GF2> X;// table of output bits
    Vec<long> I;// permutation index
    TurboCode(long=1024, const long* =coeff);
    void encode_(Vec<GF2>&, const Vec<GF2>&);
    void encode(Vec<GF2>&, const Vec<GF2>&);
    long decode(Vec<GF2>&, const double[], long, double);
    long decode(Vec<GF2>&, const Vec<double>&, long=10, double=0);
    void BCJR(Vec<double>&, const Mat<double>&);
};

void GF2XFromLong(GF2X&, long);
void DivRemLSB(GF2X&, GF2X&, const GF2X&, const GF2X&, long=-1);
void AddNoise(Vec<double>&, const Vec<GF2>&, double);
void conv(std::string&, const Vec<GF2>&);
void conv(Vec<GF2>&, const std::string&);
void conv(Vec<GF2>&, const Vec<double>&);

#endif // __TurboCode_h__
