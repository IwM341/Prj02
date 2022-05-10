#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"

#if defined(TEST1)
int main (void){
    const double Vdisp = 0.1*U0;
    const int Nmk =10000000;
    srand (time(NULL));
    auto G = [](){return (double(rand() + 1))/(RAND_MAX+1);};
    std::default_random_engine stdGen;
    std::uniform_real_distribution<double> stdUniform(1.0,0.0);
    //auto G = [&stdGen,&stdUniform]()->double{return stdUniform(stdGen);};
    auto res = MC::MCIntegrate([&G,Vdisp](){
        auto x = Velocity(G,U0,Vdisp);
        return /*x.Result.norm() **/x.Result.norm()/Vdisp*x.RemainDensity;
    },Nmk);
    auto res1 = MC::MCIntegrate([&G,Vdisp](){
        auto ksiMK = Gauss2(G,Vdisp);
        auto v  = (vec3::PolarCos(ksiMK.norm(),RandomCos(G),RandomPhi(G))+vec3::PolarCos(U0,RandomCos(G),0)).norm();
        auto ksi  = ksiMK.norm();
        return sqrt(v*v+U0*U0)/Vdisp*ksi/Vdisp*2/sqrt(2*M_PI);
    },Nmk);
    std::cout << res.Result << "\t" << res.Sigma <<"\t" << res.Sigma/sqrt(Nmk) << std::endl;
    //std::cout << res1.Result << "\t" << res1.Sigma <<"\t" << res1.Sigma/sqrt(Nmk) << std::endl;
    auto tres = integrateABP([Vdisp](double U){
        auto D = Vdisp;
        auto D2 = D*D;
        auto y =  2*U*U0/D2;
        auto V = sqrt(U*U+U0*U0);
        return exp(-fast_pow((U-U0)/D,2)/2)*(1-exp(-y))/y/fast_pow(sqrt(2*M_PI*D2),3)*4*M_PI*U*V*V/Vdisp;
    },U0*0.000001,20*U0,100000);
    std::cout << tres << std::endl;
    std::cout << res.Result - tres << std::endl;

}
#elif defined(TEST2)
int main (void){
    const int Nmk =100000000;
    srand (time(NULL));
    auto G = [](){return (double(rand()%RAND_MAX + 1))/RAND_MAX;};
    std::default_random_engine stdGen;
    std::uniform_real_distribution<double> stdUniform(0.0,1.0);
    //auto G = [&stdGen,&stdUniform]()->double{return stdUniform(stdGen);};
    auto res = MC::MCIntegrate([&G](){
        return sqrt(G()/4);
    },Nmk);;
    std::cout << res.Result << "\t" << res.Sigma <<"\t" << res.Sigma/sqrt(Nmk) << std::endl;
    std::cout << res.Result - 1.0/3 << "\t" << std::abs(res.Result - 2.0/3) <<std::endl;
    return 0;
}
#elif defined(TEST3)
int main (void){
    const int Nmk =100000000;
    srand (time(NULL));
    auto G = [](){return (double(rand()%RAND_MAX + 1))/RAND_MAX;};
    std::cout << RAND_MAX <<std::endl;

    auto V = Vector(100000,[](size_t i){return rand();});
    std::sort(V.begin(),V.end());
    std::cout << V.front() << "\t" << V.back() <<std::endl;
    return 0;
}
#endif
