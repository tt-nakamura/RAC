// uses NTL
//   http://www.shoup.net/ntl

#include "RATurboCode.h"
#include "perm.h"

RATurboCode::RATurboCode(long K_, long r_, long s_) 
: K(K_), r(r_), s(s_>0 ? s_:r),
  rK(r*K), M((rK-1)/s+1), N(K+M)
// K = message length before encoding
// r = repeat length
// s = shrink factor
// M = number of parity check equations
// N = code length after encoding
{ RandomPerm(P,rK); }

void RATurboCode::encode_(Vec<GF2>& x, const Vec<GF2>& u)
// x = encode of u
// assume length(u) == K
{
    long i,j,k;
    Vec<GF2> v;

    VectorCopy(x,u,N);
    v.SetLength(rK);
    for(i=k=0; i<K; i++)
        for(j=0; j<r; j++) v[k++] = u[i];// repeat
    permute(v,v,P);
    for(i=K, k=0; i<N; i++)
        for(j=0; j<s && k<rK; j++) x[i] += v[k++];
    for(i=K+1; i<N; i++) x[i] += x[i-1];// accumulate
}

long RATurboCode::decode(Vec<GF2>& u, const double *y,
                         long IMAX, double EPS)
// output: u = decoding of y
// input:
//   y = received code with gaussian noise (prob(u==1))
//   IMAX = max number of iterations for decoding
//   EPS = termination criterion: terminate if
//         min(|prob-0.5|) > 0.5-EPS
// return number of iterations if successful
// return IMAX+1 if failed after IMAX iterations
// assume y.length == N
// reference: S. J. Johnson
//            "Iterative Error Correction" section 6.2
{
    long i,j,k,l,m;
    double w[M],p[K],a,b,T(0.5-EPS);
    double A0[M],A1[M],B0[M],B1[M];
    double C0[M],C1[M],C2[M],C3[M];
    const double *t(y+K);
    Vec<double> v;

    for(i=0; i<K; i++) p[i] = y[i];

    u.SetLength(K);
    v.SetLength(rK);
    A0[0] = 1; A1[0] = 0;
    B0[M-1] = B1[M-1] = 0.5;

    for(m=1; m<=IMAX; m++) {
        for(i=k=0; i<K; i++) {
            v[k++] = 1 - 2*p[i];
            for(j=1; j<r; j++, k++) v[k] = v[k-1];// repeat
        }
        permute(v,v,P);
        for(i=k=0; i<M; i++) {
            w[i] = v[k++];
            for(j=1; j<s && k<rK; j++)
                w[i] *= v[k++];
        }
        for(i=0; i<M; i++) {
            a = (1 - w[i])/2;
            C0[i] = (1-t[i])*(1-a);
            C1[i] = t[i]*a;
            C2[i] = (1-t[i])*a;
            C3[i] = t[i]*(1-a);
        }
        for(i=0; i<M-1; i++) {// BCJR
            A0[i+1] = A0[i]*C0[i] + A1[i]*C2[i];
            A1[i+1] = A0[i]*C1[i] + A1[i]*C3[i];
            a = A0[i+1] + A1[i+1];
            if(a) {A0[i+1] /= a; A1[i+1] /= a; }
            else A0[i+1] = A1[i+1] = 0.5;
        }
        for(i=M-1; i>0; i--) {// BCJR
            B0[i-1] = B0[i]*C0[i] + B1[i]*C1[i];
            B1[i-1] = B0[i]*C2[i] + B1[i]*C3[i];
            b = B0[i-1] + B1[i-1];
            if(b) { B0[i-1] /= b; B1[i-1] /= b; }
            else B0[i-1] = B1[i-1] = 0.5;
        }
        for(i=k=0; i<M; i++) {
            a = A0[i]*C1[i]*B1[i] + A1[i]*C2[i]*B0[i];
            b = A0[i]*C0[i]*B0[i] + A1[i]*C3[i]*B1[i];
            if(a||b) a = (b-a)/(a+b);
            else a = 0;
            for(j=0; j<s && k<rK; j++, k++) {
                if(v[k]==0) {// exceptional case
                    b = a;
                    for(l=0; l<s && i+l<rK; l++)
                        if(j!=l) b *= v[i+l];
                    v[k] = (1-b)/2;
                }
                else v[k] = (1 - a*w[i]/v[k])/2;
            }
        }
        InvPermute(v,v,P);
        clear(u);
        for(i=k=0; i<K; i++) {
            a = y[i];
            b = 1-y[i];
            for(j=0; j<s; j++, k++) {// accumulate
                a *= v[k];
                b *= 1-v[k];
            }
            if(b+=a) p[i] = a/b;
            else p[i] = 0.5;
            if(p[i]>=0.5) set(u[i]);
        }
        for(i=0; i<K; i++)
            if(fabs(p[i]-0.5) <= T) break;
        if(i==K) break;
    }
    return m;
}

void RATurboCode::encode(Vec<GF2>& c, const Vec<GF2>& b)
// c = encoding of b
// c.length = multiple of N
{
    Vec<GF2> v,u(b);
    c.SetLength(0);
    while(!IsZero(u)) {
        VectorCopy(v,u,K);
        encode_(v,v);
        c.append(v);
        shift(u,u,-K);
    }
}

long RATurboCode::decode(Vec<GF2>& b, const Vec<double>& d,
                         long IMAX, double EPS)
// b = decoding of d
// assume d.length == multiple of N
// return 0 if successful, -1 otherwise
{
    long i,k(0);
    Vec<GF2> v;
    b.SetLength(0);
    for(i=0; i<d.length(); i+=N) {
        if(decode(v, &d[i], IMAX, EPS) > IMAX) k=-1;
        b.append(v);
    }
    return k;
}