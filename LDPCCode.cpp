// uses NTL
//   http://www.shoup.net/ntl

#include "LDPCCode.h"

#define LDPC_NTRY 20

LDPCCode::LDPCCode(long N1, long K1, long W1)
: N(N1), K(K1)
// N1 = code length; K1 = message length
// W1 = column weight (must be odd)
// reference: T. J. Richardson and R. L. Urbanke
//   IEEE Transaction on Information Theory 47 (2001) 638
{
    long M(N-K);// number of parity checks
    long W(W1|1);// make W odd
    long i,j,k,l,c,g; 
    int count[M],ind[M];
    Vec<GF2> v;
    Mat<GF2> HB,HD,HE,HT,U,V;

    if(W<3) Error("too low density");
    H.SetDims(N,M);
    I.SetLength(M);
    J.SetLength(N);
    v.SetLength(N);
a:
    for(j=0; j<M; j++) count[j] = 0;
    for(i=0; i<N; i++) {
        clear(v);
        indexr(count, ind, M);
        for(j=c=0; j<M && c<W; j++) {
            l = ind[j];
            for(k=0; k<i; k++)
                if(IsOne(H[k][l]) && IsOne(v[k]))
                    break;// detect cycle
            if(k<i) continue;
            set(H[i][l]);
            count[l]++;
            for(k=0; k<i; k++)
                if(IsOne(H[k][l])) set(v[k]);
            c++;
        }
        if(c<W && j==M) Error("too high density");
    }

    for(i=l=0; i<N; i++) {// approximate triangulation
        for(j=0; j<M; j++)
            if(IsOne(H[i][j])) break;
        if(j<l) continue;
        if(i>l) swap(H[i], H[l]);
        if(j>l) for(k=0; k<N; k++)
            swap(H[k][j], H[k][l]);
        l++;
    }
    g = M-l;// gap
    HB.SetDims(g,l);
    HD.SetDims(g,g);
    HE.SetDims(l,g);
    HT.SetDims(l,l);

    for(i=0; i<l; i++) {
        for(j=i; j<l; j++) HT[i][j] = H[i][j];
        for(j=0; j<g; j++) HE[i][j] = H[i][j+l];
    }
    for(k=l-1, i=N-1; k>=0; k--, i--)
        swap(H[i], H[k]);
    inv(U,HT);

    for(k=0; k<LDPC_NTRY; k++) {
        for(i=0; i<g; i++) {
            for(j=0; j<l; j++) HB[i][j] = H[i+K][j];
            for(j=0; j<g; j++) HD[i][j] = H[i+K][j+l];
        }
        mul(V,HB,U);
        V *= HE;
        HD -= V;
        if(gauss(V=HD) == g) break;
        for(i=K+g-1; i>=0; i--)
            swap(H[i], H[RandomBnd(i+1)]);
    }
    if(k==LDPC_NTRY) { clear(H); goto a; }

    inv(D,HD);

    A.SetLength(l);
    B.SetLength(l);
    T.SetLength(l);
    C.SetLength(g);
    E.SetLength(g);

    for(i=0; i<N; i++)// sparse representation
        for(j=0; j<M; j++)
            if(IsOne(H[i][j]))
            { I[j].append(i); J[i].append(j); }

    for(j=0; j<l; j++) {
        for(i=0; i<K; i++)
            if(IsOne(H[i][j])) A[j].append(i);
        for(i=0; i<g; i++)
            if(IsOne(HB[i][j])) B[j].append(i);
        for(i=0; i<j; i++)
            if(IsOne(HT[i][j])) T[j].append(i);
    }
    for(j=0; j<g; j++) {
        for(i=0; i<K; i++)
            if(IsOne(H[i][j+l])) C[j].append(i);
        for(i=0; i<l; i++)
            if(IsOne(HE[i][j])) E[j].append(i);
    }
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

template<class T>
void back_subst(Vec<T>& y, const Vec<T>& x, const Vec<Vec<long> >& A)
// solve yA=x for y where A is upper triangular and sparse
// &y==&x is allowed
{
    long i,j;
    if(&y!=&x) y=x;
    for(i=0; i<A.length(); i++)
        for(j=0; j<A[i].length(); j++)
            y[i] -= y[A[i][j]];
}

void LDPCCode::encode_(Vec<GF2>& x, const Vec<GF2>& a)
// x = encoding of a; &x==&a is allowed
{
    Vec<GF2> u,v,b,c;

    mul_sparse(u,a,A);
    back_subst(v,u,T);
    mul_sparse(b,v,E);
    mul_sparse(v,a,C);
    b -= v;
    b *= D;

    mul_sparse(c,b,B);
    c += u;
    back_subst(c,c,T);

    VectorCopy(x,a,K);
    x.append(b);
    x.append(c);
}

long LDPCCode::decode_(Vec<GF2>& x, const double y[], long IMAX)
// output: x = decoding of y by sum product algorithm
// input: y = received code with gaussian noise (prob(x==1))
//        IMAX = max number of iteration for decoding
// return number of iterations if successful
// return IMAX+1 if failed after IMAX iterations
// reference: S. J. Johnson
//   "Iterative Error Correction" algorithm 2.3
{
    long i,j,k,l(1),M(N-K);
    double p[N],q[N],s,t;
    Mat<double> P,Q;
    Vec<GF2> v;
    
    P.SetDims(N,M);
    Q.SetDims(N,M);
    
    for(i=0; i<N; i++)
        for(j=0; j<J[i].length(); j++)
            P[i][J[i][j]] = y[i];
    x.SetLength(N);
    for(;;) {
        for(j=0; j<M; j++) {
            s = 1;
            for(i=0; i<I[j].length(); i++) {
                p[i] = 1 - 2*P[I[j][i]][j];
                if(p[i] == 0) break;
                s *= p[i];
            }
            if(i==I[j].length()) {
                for(i=0; i<I[j].length(); i++)
                    Q[I[j][i]][j] = (1 - s/p[i])/2;
            }
            else {// in case y[i]==0
                for(k=0; k<I[j].length(); k++)
                    if(k!=i) Q[I[j][k]][j] = 0.5;
                for(k=i+1; k<I[j].length(); k++)
                    s *= 1 - 2*P[I[j][k]][j];
                Q[I[j][i]][j] = (1-s)/2;
            }
        }
        clear(x);
        for(i=0; i<N; i++) {
            p[i] = y[i];
            q[i] = 1 - y[i];
            for(j=0; j<J[i].length(); j++) {
                p[i] *= Q[i][J[i][j]];
                q[i] *= 1 - Q[i][J[i][j]];
            }
            if(p[i] >= q[i]) set(x[i]);
        }
        mul_sparse(v,x,I);// syndrome
        if(IsZero(v) || ++l > IMAX) break;
        for(i=0; i<N; i++) {
            for(j=0; j<J[i].length(); j++) {
                if(Q[i][J[i][j]]==0) {// exceptional case
                    s = y[i];
                    for(k=0; k<J[i].length(); k++)
                        if(j!=k) s *= Q[i][J[i][k]];
                }
                else s = p[i]/Q[i][J[i][j]];
                if(Q[i][J[i][j]]==1) {// exceptional case
                    t = 1 - y[i];
                    for(k=0; k<J[i].length(); k++)
                        if(j!=k) t *= 1 - Q[i][J[i][k]];
                }
                else t = q[i]/(1 - Q[i][J[i][j]]);
                if(t+=s) P[i][J[i][j]] = s/t;
                else P[i][J[i][j]] = 0.5;
            }
        }
    }
    x.SetLength(K);
    return l;
}

void LDPCCode::encode(Vec<GF2>& x, const Vec<GF2>& a)
// x = encoding of a; x.length = multiple of N
{
    Vec<GF2> v,u(a);
    x.SetLength(0);
    while(!IsZero(u)) {
        VectorCopy(v,u,K);
        encode_(v,v);
        x.append(v);
        shift(u,u,-K);
    }
}

long LDPCCode::decode(Vec<GF2>& x, const Vec<double>& y, long IMAX)
// x = decoding of y; assume y.length == multiple of N
// IMAX = max number of iteration for decoding
// return 0 if successful, -1 if any error exist
{
    long i,k(0);
    Vec<GF2> v;
    x.SetLength(0);
    for(i=0; i<y.length(); i+=N) {
        if(decode_(v, &y[i], IMAX) > IMAX) k=-1;
        x.append(v);
    }
    return k;
}