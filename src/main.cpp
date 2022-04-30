#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"


/*
int main(void){
    double mp = 0.938;
    double me = 0.51e-3;
    auto PhiFactor55p = [mp](double mk,double q)->double{
            return fast_pow(q*q/(mp*mk),2);
        };
    auto PhiFactor55e = [me,mp](double mk,double q)->double{
            return fast_pow(q*q/(me*mk),2);
        };
    auto PhiFactorSS = [mp](double mk,double q)->double{
            return 1.0;
        };

    const auto& BM = BodyModel(2.03e-4,"D:/Important/articles/DMFramework/celstial_models/jupiter_model.dat");


    //std::vector<double> Vmk = {0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5};
    std::vector<double> Vmk = {0.4};
    size_t N = Vmk.size();

    std::vector<std::pair<double,double>> Fel(N);
    std::vector<std::pair<double,double>> FelE(N);
    std::vector<std::pair<double,double>> Fione(N);
    std::vector<std::pair<double,double>> Fion(N);


    #pragma omp parallel for
    for(size_t i=0;i<N;++i){
        Fel[i] = SupressFactor(Vmk[i],0,ELASTIC,PROTON,PhiFactor55p,BM,"H",100000000,U0*1.1,U0);
        //FelE[i] = SupressFactor(Vmk[i],0,ELASTIC,ELECTRON,PhiFactor55e,BM,"H",1000000,U0*1.1,U0);
        //Fion[i] = SupressFactor(Vmk[i],0,IONIZATION,PROTON,PhiFactor55p,BM,"H",1000000,U0*1.1,U0);
        //Fione[i] = SupressFactor(Vmk[i],0,IONIZATION,ELECTRON,PhiFactor55e,BM,"H",1000000,U0*1.1,U0);
    }

    std::ofstream ofsEl("el55.dat");
    std::ofstream ofsIon("ion55.dat");
    std::ofstream ofsIonE("ione55.dat");
    std::ofstream ofsElE("elE55.dat");

    ofsEl << "mk\tF\n"<< Function::GridFunction1(Vmk,Fel) << std::endl;
    ofsIon << "mk\tF\n"<< Function::GridFunction1(Vmk,Fion) << std::endl;
    ofsIonE << "mk\tF\n"<< Function::GridFunction1(Vmk,Fione) << std::endl;
    ofsElE << "mk\tF\n"<< Function::GridFunction1(Vmk,FelE) << std::endl;

    return 0;

}
*/

int main (void){
    auto G = [](){return (double(rand()%RAND_MAX + 1))/RAND_MAX;};
    auto res = MC::MCIntegrate([&G](){
        auto x = Velocity(G,U0);
        return x.Result.norm() *x.RemainDensity;
    },10000000);
    std::cout << res.Result << "\t" << res.Sigma << std::endl;
    std::cout << integrateAB5P([](double U){
        auto D = 1.1*U0;
        auto D2 = D*D;
        auto y =  2*U*U0/D2;
        auto V = sqrt(U*U+U0*U0);
        return exp(-fast_pow((U-U0)/D,2)/2)*(1-exp(-y))/y/fast_pow(sqrt(2*M_PI*D2),3)*4*M_PI*V*U*V;
    },U0*0.0001,30*U0,1000) << std::endl;

}
