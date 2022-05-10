#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>

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


    std::vector<double> Vmk = {0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.20};
    //std::vector<double> Vmk = {0.4};
    size_t N = Vmk.size();

    std::vector<std::pair<double,double>> Fel(N);
    std::vector<std::pair<double,double>> FelE(N);
    std::vector<std::pair<double,double>> Fione(N);
    std::vector<std::pair<double,double>> Fion(N);
    std::vector<std::pair<double,double>> Fmgd(N);

    auto t_start = clock();
    #pragma omp parallel for
    for(size_t i=0;i<N;++i){
        Fel[i] = SupressFactor(Vmk[i],0,ELASTIC,PROTON,PhiFactor55p,BM,"H",1000000,U0*1.1,U0);
        FelE[i] = SupressFactor(Vmk[i],0,ELASTIC,ELECTRON,PhiFactor55e,BM,"H",1000000,U0*1.1,U0);
        Fion[i] = SupressFactor(Vmk[i],0,IONIZATION,PROTON,PhiFactor55p,BM,"H",1000000,U0*1.1,U0);
        Fione[i] = SupressFactor(Vmk[i],0,IONIZATION,ELECTRON,PhiFactor55e,BM,"H",1000000,U0*1.1,U0);
        Fmgd[i] = SupressFactor(Vmk[i],0,MIGDAL,PROTON,PhiFactor55p,BM,"H",1000000,U0*1.1,U0);
    }
    std::cout << "time = " << (clock()- t_start) <<std::endl;
    std::ofstream ofsEl("el55.dat");
    std::ofstream ofsIon("ion55.dat");
    std::ofstream ofsIonE("ione55.dat");
    std::ofstream ofsElE("elE55.dat");
    std::ofstream ofsMgd("mgd55.dat");

    ofsEl << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,Fel) << std::endl;
    ofsIon << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,Fion) << std::endl;
    ofsIonE << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,Fione) << std::endl;
    ofsElE << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,FelE) << std::endl;

    ofsMgd << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,Fmgd) << std::endl;
    return 0;

}



