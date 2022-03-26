#include "RACode.h"
#include "TurboCode.h"
#include "LDPCCode.h"

main() {
    double sigma(0.5),t,t1[3],t2[3],EPS(0.01);
    long K(256),N(512),IMAX(32),M(256);
    long i;
    Vec<GF2> u,v;
    Vec<double> y;
    GF2X g;

    RACode c(K);
    TurboCode b(K);
    LDPCCode l(N,K);
    for(i=0; i<3; i++) t1[i] = t2[i] = 0;
    for(i=0; i<M; i++) {
        random(u,K);
        // RA
        t = GetTime();
        c.encode(v,u);
        t1[0] += GetTime() - t;
        AddNoise(y,v,sigma);
        t = GetTime();
        c.decode(v,y,IMAX);
        t2[0] += GetTime() - t;
        // Turbo
        t = GetTime();
        b.encode(v,u);
        t1[1] += GetTime() - t;
        AddNoise(y,v,sigma);
        t = GetTime();
        b.decode(v,y,IMAX,EPS);
        t2[1] += GetTime() - t;
        // LDPC
        t = GetTime();
        l.encode(v,u);
        t1[2] += GetTime() - t;
        AddNoise(y,v,sigma);
        t = GetTime();
        l.decode(v,y);
        t2[2] += GetTime() - t;
    }
    for(i=0; i<3; i++) {
        std::cout << t1[i]/M << ' ';
        std::cout << t2[i]/M << ' ';
    }
    std::cout << std::endl;
}