// uses NTL
//   http://www.shoup.net/ntl

#include "TurboCode.h"
#include "perm.h"

const long TurboCode::coeff[2] = {037, 021};

TurboCode::TurboCode(long K1, const long *c)
// K1 = message length
// c = array of coefficients of generator polynomials
: K(K1)
{
    long i,j;
    for(i=0; i<2; i++)
        GF2XFromLong(G[i], c[i]);
    m = max(deg(G[0]), deg(G[1]));
    M = 1<<(m+1);
    X.SetDims(2,M);
    for(i=0; i<2; i++)
        for(j=1; j<M; j++)
            if(weight(j&c[i])&1) set(X[i][j]);
    RandomPerm(I,K);
}

void TurboCode::encode_(Vec<GF2>& t, const Vec<GF2>& s)
// t = encoding of s
// assume length(s)==K; assume &t!=&s
{
    long i,j,L(K+m);
    GF2X f[2],q,r;
    Vec<GF2> v;
    permute(v,s,I);
    conv(f[0],s);
    conv(f[1],v);
    for(i=1; i>=0; i--) {
        DivRemLSB(q,r,f[i],G[0],K-1);
        mul(f[i],q,G[1]);
        f[i].SetLength(L);
    }
    r.SetLength(L);
    t.SetLength(L<<1);
    for(i=j=0; i<K; i++, j+=2) {
        t[j] = s[i];
        t[j+1] = f[i&1][i];// puncture
    }
    for(; i<L; i++, j+=2) {// tail
        t[j] = r[i];
        t[j+1] = f[0][i];
    }
}

void TurboCode::encode(Vec<GF2>& t, const Vec<GF2>& s)
// t = encoding of s; &t==&s is allowed
{
    Vec<GF2> c,d,u(s);
    t.SetLength(0);
    while(!IsZero(u)) {
        VectorCopy(c,u,K);
        encode_(d,c);
        t.append(d);
        shift(u,u,-K);
    }
}

void TurboCode::BCJR(Vec<double>& p, const Mat<double>& P)
// output:
//   p[i] = prob(message[i]==1) after decoding
// input:
//   p[i] = prob(message[i]==1) before decoding
//   P = table of log(transition probability)
{
    long i,j,k,L(p.length()),N(1<<m);
    double a,b,c;
    Mat<double> A,B,C(P);
    A.SetDims(L,N);
    B.SetDims(L,N);
    for(i=0; i<L; i++)
        for(j=0; j<M; j++)
            C[i][j] *= (IsOne(X[0][j]) ? p[i] : 1-p[i]);

    A[0][0] = B[L-1][0] = 1;
    for(j=1; j<N; j++) A[0][j] = 0;
    for(j=1; j<N; j++) B[L-1][j] = (L==K);
    for(i=0; i<L-1; i++) {
        for(a=j=0; j<N; j++) {
            k = j|N;
            c = A[i][j>>1]*C[i][j] + A[i][k>>1]*C[i][k];
            A[i+1][j] = c;
            a += c;
        }
        if(a) for(j=0; j<N; j++) A[i+1][j] /= a;
        else  for(j=0; j<N; j++) A[i+1][j] = 1;        
    }
    for(i=L-1; i>0; i--) {
        for(b=j=0; j<M; j+=2) {
            k = j&(N-1);
            c = B[i][k]*C[i][j] + B[i][k+1]*C[i][j+1];
            B[i-1][j>>1] = c;
            b += c;
        }
        if(b) for(j=0; j<N; j++) B[i-1][j] /= b;
        else  for(j=0; j<N; j++) B[i-1][j] = 1;
    }
    for(i=0; i<L; i++) {
        for(a=b=j=0; j<M; j++) {
            c = A[i][j>>1]*C[i][j]*B[i][j&(N-1)];
            if(IsOne(X[0][j])) a += c;
            else b += c;
        }
        if(b+=a) p[i] = a/b;
        else p[i] = 0.5;
    }
}

long TurboCode::decode(Vec<GF2>& s, const double t[],
                       long IMAX, double EPS)
// output: s = decoding of t
// input:
//   t = received code with gaussian noise (prob(s==1))
//   IMAX = max number of iterations for decoding
//   EPS = termination criterion: terminate if
//         min(|prob-0.5|) > 0.5-EPS
// return number of iterations
// assume t.length == 2*L
{
    long i,j,k,L(K+m);
    double T(0.5-EPS);
    Vec<double> p,q;
    Mat<double> P[2];// table of log(transition probablity)
    p.SetLength(L);
    P[0].SetDims(L,M);
    P[1].SetDims(K,M);
    for(i=k=0; i<L; i++, k+=2) p[i] = t[k];
    for(j=0; j<M; j++) {
        for(i=0, k=1; i<K; i++, k+=2) {
            P[i&1][i][j] = (IsOne(X[1][j]) ? t[k] : 1-t[k]);
            P[(i&1)^1][i][j] = 0.5;// puncture
        }
        for(; i<L; i++, k+=2)// tail
            P[0][i][j] = (IsOne(X[1][j]) ? t[k] : 1-t[k]);
    }
    for(j=1; j<=IMAX; j++) {
        BCJR(p,P[0]);
        p.SetLength(K);
        permute(q,p,I);
        BCJR(q,P[1]);
        InvPermute(p,q,I);
        p.SetLength(L);
        for(i=0; i<L; i++)
            if(fabs(p[i]) <= T) break;
        if(i==L) break;
    }
    s.SetLength(K);
    clear(s);
    for(i=0; i<K; i++)
        if(p[i]>=0.5) set(s[i]);
    return j;
}

long TurboCode::decode(Vec<GF2>& s, const Vec<double>& t,
                       long IMAX, double EPS)
// output: s = decoding of t
// input:
//   t = received code with gaussian noise
//   IMAX = max number of iterations for decoding
//   EPS = termination criterion: terminate if
//         min(|prob-0.5|) > 0.5-EPS
// return 0 if number of iterations <= IMAX for all blocks
//        -1 otherwise
// assume t.length == multiple of 2*L
{
    long i,k(0),N((K+m)<<1);
    Vec<GF2> c;
    s.SetLength(0);
    for(i=0; i<t.length(); i+=N) {
        if(decode(c, &t[i], IMAX, EPS) > IMAX) k=-1;
        s.append(c);
    }
    return k;
}