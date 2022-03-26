#include "RACode.h"
#include "TurboCode.h"
#include "LDPCCode.h"

main() {
    double sigma(0.6);
    Vec<double> d;
    Vec<GF2> u,v;
    std::string s;
    
    s = "Hello world of error correcting codes.";
    conv(v,s);

    std::cout << "RA code:" << std::endl;
    RACode c(512,3);
    c.encode(u,v);
    AddNoise(d,u,sigma);
    c.decode(u,d);
    conv(s,u);
    std::cout << s << std::endl;
    
    std::cout << "Turbo code:" << std::endl;
    TurboCode b(512);
    b.encode(u,v);
    AddNoise(d,u,sigma);
    b.decode(u,d);
    conv(s,u);
    std::cout << s << std::endl;    

    std::cout << "LDPC code:" << std::endl;
    LDPCCode l(512,3);
    l.encode(u,v);
    AddNoise(d,u,sigma);
    l.decode(u,d);
    conv(s,u);
    std::cout << s << std::endl;    
}