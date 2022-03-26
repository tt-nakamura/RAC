// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/GF2X.h>
using namespace NTL;

void GF2XFromLong(GF2X& x, long a) {
    clear(x);
    x.SetLength(NumBits(a));
    for(long i=0; a; a>>=1, i++) if(a&1) set(x[i]);
}

long LongFromGF2X(const GF2X& a) {
    if(IsZero(a)) return 0;
    long x(1),i;
    for(i=deg(a)-1; i>=0; i--) {
        x<<=1;
        if(IsOne(a[i])) x|=1;
    }
    return x;
}

void DivRemLSB(GF2X& q, GF2X& r, const GF2X& a, const GF2X& b, long n)
// long division from Least Significant Bit
// q,r = quotient and remainder of a/b such that a=bq+r
// n = upper bound of deg(q); deg(q)<=n; r[i]=0 for 0<=i<=n;
// if n<0, n is set to deg(a)
// assume ConstTerm(b) != 0; &r==&a is allowed
{
    long i;
    GF2X s(b);
    if(n<0) n = deg(a);
    if(&r!=&a) r=a;
    clear(q);
    for(i=0; i<=n && !IsZero(r); i++, s<<=1)
        if(IsOne(r[i])) { SetCoeff(q,i); r-=s; }
}

void conv(Vec<GF2>& x, const std::string& a) {
    GF2X f;
    GF2XFromBytes(f, (const unsigned char*)a.c_str(), a.length());
    conv(x,f);
}

void conv(std::string& x, const Vec<GF2>& a) {
    GF2X f;
    conv(f,a);
    x.resize(NumBytes(f));
    BytesFromGF2X((unsigned char*)x.data(), f, x.length());
}