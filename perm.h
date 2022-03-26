// uses NTL
//   http://www.shoup.net/ntl

#ifndef __perm_h__
#define __perm_h__

#include<NTL/vector.h>

void RandomPerm(Vec<long>& p, long n);

template<class T>
void permute(Vec<T>& v, const Vec<T>& u, const Vec<long>& p)
// v = permuted elements of u according to p
// assume length(u)==length(p); &v==&u is allowed
{
    if(&v==&u) { Vec<T> w(u); permute(v,w,p); return; }
    long i;
    v.SetLength(u.length());
    for(i=0; i<p.length(); i++) v[i] = u[p[i]];
}

template<class T>
void InvPermute(Vec<T>& v, const Vec<T>& u, const Vec<long>& p)
// v = permuted elements of u according to inverse of p
// assume length(u)==length(p); &v==&u is allowed
{
    if(&v==&u) { Vec<T> w(u); InvPermute(v,w,p); return; }
    long i;
    v.SetLength(u.length());
    for(i=0; i<p.length(); i++) v[p[i]] = u[i];
}

#endif // __perm_h__