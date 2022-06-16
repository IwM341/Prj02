#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>
const double mp = 0.938;
const double me = 0.51e-3;
auto PhiFactor55p = [](double mk,double q)->double{
    return fast_pow(q*q/(mp*mk),2);
};
auto PhiFactor55e = [](double mk,double q)->double{
    return fast_pow(q*q/(me*mk),2);
};
auto PhiFactorSS = [](double mk,double q)->double{
    return 1.0;
};

#if defined(MASS)
int main(void){


    const auto& BM = BodyModel(2.03e-4,"D:/Important/articles/DMFramework/celstial_models/jupiter_model.dat");


    std::vector<double> Vmk = {0.5,1,2,10};
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
#elif defined(DELTAMASS)
int main(void){


    const auto& BM = BodyModel(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");
    const double mk = 10.0;

    std::vector<double> Vdmk = Vector(10,[](size_t i){return i*1e-4;});
    //PVAR(Vdmk);
    //PVAR(BM["Vesc"]);
    //std::vector<double> Vmk = {0.4};
    size_t N = Vdmk.size();

    std::vector<std::pair<double,double>> Fel(N);
    std::vector<std::pair<double,double>> FelE(N);
    std::vector<std::pair<double,double>> Fione(N);
    std::vector<std::pair<double,double>> Fion(N);
    std::vector<std::pair<double,double>> Fmgd(N);

    auto t_start = clock();
    #pragma omp parallel for
    for(size_t i=0;i<N;++i){
        Fel[i] = SupressFactor(mk,Vdmk[i],ELASTIC,PROTON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
        FelE[i] = SupressFactor(mk,Vdmk[i],ELASTIC,ELECTRON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
        Fion[i] = SupressFactor(mk,Vdmk[i],IONIZATION,PROTON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
        Fione[i] = SupressFactor(mk,Vdmk[i],IONIZATION,ELECTRON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
        Fmgd[i] = SupressFactor(mk,Vdmk[i],MIGDAL,PROTON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
    }
    std::cout << "time = " << (clock()- t_start) <<std::endl;
    std::ofstream ofsEl("el55dm.dat");
    std::ofstream ofsIon("ion55dm.dat");
    std::ofstream ofsIonE("ione55dm.dat");
    std::ofstream ofsElE("elE55dm.dat");
    std::ofstream ofsMgd("mgd55dm.dat");

    ofsEl << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,Fel) << std::endl;
    ofsIon << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,Fion) << std::endl;
    ofsIonE << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,Fione) << std::endl;
    ofsElE << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,FelE) << std::endl;

    ofsMgd << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,Fmgd) << std::endl;
    return 0;

}
#endif


