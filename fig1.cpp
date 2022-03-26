#include "RACode.h"
#include "RATurboCode.h"
#include "TurboCode.h"
#include "LDPCCode.h"
#include<fstream>

main() {
    long N(512),K(256),IMAX(32);
    long n(10),M1(400),M;
    double s1(0.65), s2(0.9), EPS(0.01);
    double sigma, e[4], w, ds((s2-s1)/n);
    long i,j,k;
    Vec<GF2> u,v;
    Vec<double> y;
    std::ofstream f("fig1.txt");

    RACode c(K);
    RATurboCode r(K);
    TurboCode b(K);
    LDPCCode l(N,K);

    for(i=0; i<=n; i++) {
        sigma = s1 + i*ds;
        std::cout << sigma;
        f << sigma;
        M = long(M1/sigma);
        for(k=0; k<4; k++) e[k] = 0;
        for(j=0; j<M; j++) {
            random(u,K);
            // RACode
            c.encode(v,u);
            AddNoise(y,v,sigma);
            c.decode(v,y,IMAX);
            w = weight(v-=u);
            e[0] += w/K;
            // RATurboCode
            r.encode(v,u);
            AddNoise(y,v,sigma);
            r.decode(v,y,IMAX,EPS);
            w = weight(v-=u);
            e[1] += w/K;
            // Turbo
            b.encode(v,u);
            AddNoise(y,v,sigma);
            b.decode(v,y,IMAX,EPS);
            w = weight(v-=u);
            e[2] += w/K;
            // LDPC
            l.encode(v,u);
            AddNoise(y,v,sigma);
            l.decode(v,y,IMAX);
            w = weight(v-=u);
            e[3] += w/K;
        }
        for(k=0; k<4; k++) {
            e[k] /= M;
            std::cout << ' ' << e[k];
            f << ' ' << e[k];
        }
        std::cout << std::endl;
        f << std::endl;
    }
}