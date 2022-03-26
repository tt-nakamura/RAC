// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/vec_ZZ.h>
using namespace NTL;

void RandomPerm(Vec<long>& p, long n)
// p = random permutation of {0,1,2,...,n-1}
{
    long i;
    p.SetLength(n);
    for(i=0; i<n; i++) p[i] = i;
    for(i=n-1; i>0; i--)
        swap(p[i], p[RandomBnd(i+1)]);
}