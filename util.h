#ifndef __util_h__
#define __util_h__

#include <cstdlib>
#include <string>
#include <iostream>
using namespace std;

template<class T>
inline void SWAP(T& a, T& b) { T c(a); a=b; b=c; }

template<class T>
inline const T SQR(const T a) {return a*a;}

template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

template<class T>
inline const T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline void error(const string error_text)
{
    cerr << error_text << endl;
    exit(1);
}

#endif // __util_h__