// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/mat_GF2.h>
#include "RACode.h"
#include "perm.h"

#define RA_NTRY 20 // max number of trials to make check matrix

RACode::RACode(long K_, long r_, long s_) 
: K(K_), r(r_), s(s_>0 ? s_:r),
  rK(r*K), M((rK-1)/s+1), N(K+M)
// K = message length before encoding
// r = repeat length
// s = shrink factor
// M = number of parity check equation
// N = code length after encoding
{
    long i,j,k,l,m,c(0);
    Mat<GF2> H;// parity check matrix (transpose)
    Vec<GF2> d;

    H.SetDims(N,M);
    I.SetLength(M);
    J.SetLength(N);
    d.SetLength(M);
a:
    RandomPerm(P,rK);
    for(k=0; k<rK; k++) {
        j = k/s;
        for(m=k; m<rK; m++) {
            i = P[m]/r;
            if(IsOne(H[i][j])) continue;
            if(j && IsOne(H[i][j-1])) continue;
            for(l=0; l<j; l++)
                if(IsOne(H[i][l]) && IsOne(d[l]))
                    break;// detect cycle
            if(l==j) break;
        }
        if(m==rK) {
            if(c++ == RA_NTRY)
                Error("too large r");
            clear(H);
            goto a;// retry
        }
        set(H[i][j]);
        if(k%s == s-1) clear(d);
        else for(l=0; l<j; l++)
            if(IsOne(H[i][l])) set(d[l]);
        if(m>k) swap(P[m], P[k]);
    }
    for(j=0; j<M; j++) set(H[K+j][j]);// accumulate
    for(j=1; j<M; j++) set(H[K+j-1][j]);
    for(i=0; i<N; i++)// sparse representation
        for(j=0; j<M; j++)
            if(IsOne(H[i][j]))
            { I[j].append(i); J[i].append(j); }
}

void RACode::encode_(Vec<GF2>& x, const Vec<GF2>& u)
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

template<class T>
void mul_sparse(Vec<T>& y, const Vec<T>& x, const Vec<Vec<long> >& A)
// y=xA (A is sparse); assume &y!=&x
{
    long i,j;
    y.SetLength(A.length());
    clear(y);
    for(i=0; i<A.length(); i++)
        for(j=0; j<A[i].length(); j++)
            y[i] += x[A[i][j]];
}

long RACode::decode(Vec<GF2>& u, const double *y, long IMAX)
// output:
//   u = decoding of y by sum product algorithm
// input:
//   y = received code with gaussian noise
//       y[i] = prob(code[i]==1)
//   IMAX = maximum number of iteration
// return number of iterations if successful
// return IMAX+1 if failed after IMAX iterations
// reference: S. J. Johnson
//            "Iterative Error Correction" section 6.2
{
    long i,j,k,l(1);
    double p[N],q[N],a,b;
    Mat<double> P,Q;
    Vec<GF2> v;

    P.SetDims(N,M);
    Q.SetDims(N,M);
    
    for(i=0; i<N; i++)
        for(j=0; j<J[i].length(); j++)
            P[i][J[i][j]] = y[i];
    u.SetLength(N);
    for(;;) {
        for(j=0; j<M; j++) {
            a = 1;
            for(i=0; i<I[j].length(); i++) {
                p[i] = 1 - 2*P[I[j][i]][j];
                if(p[i]==0) break;
                a *= p[i];
            }
            if(i==I[j].length()) {
                for(i=0; i<I[j].length(); i++)
                    Q[I[j][i]][j] = (1 - a/p[i])/2;
            }
            else {// in case y[i]==0.5
                for(k=0; k<I[j].length(); k++)
                    if(k!=i) Q[I[j][k]][j] = 0.5;
                for(k=i+1; k<I[j].length(); k++)
                    a *= 1 - 2*P[I[j][k]][j];
                Q[I[j][i]][j] = (1-a)/2;
            }
        }
        clear(u);
        for(i=0; i<N; i++) {
            p[i] = y[i];
            q[i] = 1 - y[i];
            for(j=0; j<J[i].length(); j++) {
                p[i] *= Q[i][J[i][j]];
                q[i] *= 1 - Q[i][J[i][j]];
            }
            if(p[i] >= q[i]) set(u[i]);
        }
        mul_sparse(v,u,I);// syndrome
        if(IsZero(v) || ++l > IMAX) break;
        for(i=0; i<N; i++) {
            for(j=0; j<J[i].length(); j++) {
                if(Q[i][J[i][j]]==0) {// exceptional case
                    a = y[i];
                    for(k=0; k<J[i].length(); k++)
                        if(j!=k) a *= Q[i][J[i][k]];
                }
                else a = p[i]/Q[i][J[i][j]];
                if(Q[i][J[i][j]]==1) {// exceptional case
                    b = 1 - y[i];
                    for(k=0; k<J[i].length(); k++)
                        if(j!=k) b *= 1 - Q[i][J[i][k]];
                }
                else b = q[i]/(1 - Q[i][J[i][j]]);
                if(b+=a) P[i][J[i][j]] = a/b;
                else P[i][J[i][j]] = 0.5;
            }
        }
    }
    u.SetLength(K);
    return l;
}

void RACode::encode(Vec<GF2>& c, const Vec<GF2>& b)
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

long RACode::decode(Vec<GF2>& b, const Vec<double>& d, long IMAX)
// b = decoding of d
// assume d.length == multiple of N
// return 0 if successful, -1 otherwise
{
    long i,k(0);
    Vec<GF2> v;
    b.SetLength(0);
    for(i=0; i<d.length(); i+=N) {
        if(decode(v, &d[i], IMAX) > IMAX) k=-1;
        b.append(v);
    }
    return k;
}