#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <utils>
#include <random>
#include <cmath>

#define U0 0.7667e-3
#define alpha 0.0073

template <class Generator> 
inline double RandomCos(Generator G){
	return G()*2.0-1;
}
template <class Generator> 
inline double RandomPhi(Generator G){
	return G()*2*M_PI;
}

template <class Generator>
inline double RandomPhi(Generator G,double phi0,double phi1){
    return phi0 + G()*(phi1-phi0);
}

template <class Generator>
inline vec3 RandomN(Generator G){
    return vec3::PolarCos(1.0,RandomCos(G),RandomPhi(G));
}

template <class Generator> 
inline MC::MCResult<vec3> Boltsman(Generator G, double Vdisp,double Vmin = 0){
	double V = sqrt(Vmin*Vmin-2*Vdisp*Vdisp*log(G()));
    return MC::MCResult<vec3>(vec3::PolarCos(V,RandomCos(G),RandomPhi(G)),
				exp(-Vmin*Vmin/(2*Vdisp*Vdisp))*sqrt(2/M_PI)*V/Vdisp);
}

template <class Generator>
inline double Gauss2Norm(Generator G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(G()));
    return V;
}

template <class Generator>
inline vec3 Gauss2(Generator G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(G()));
    return vec3::PolarCos(V,0,RandomPhi(G));
}

template <class Generator>
inline double Gauss(Generator G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(G()));
    return V*cos(RandomPhi(G));
}
template <class Generator>
inline vec3 Gauss3(Generator G, double Vdisp){
    double V1 = Vdisp*sqrt(-2*log(G()));
    double V2 = Vdisp*sqrt(-2*log(G()));
    double phi1 = RandomPhi(G);
    double phi2 = RandomPhi(G);
    return vec3(V1*cos(phi1),V1*sin(phi1),V2*cos(phi2));
}

template <class Generator>
/*MK generator of input velocity*/
inline MC::MCResult<vec3> Velocity(Generator G,double VescTmp,
						double Vdisp = 1.1*U0,double mU0 = U0){
    auto ksi = sqrt(-2*log(G()));
    auto phi = RandomPhi(G,0,M_PI);


    auto sinPhi = sin(phi);
    auto cosPhi = cos(phi);

    auto u0 = mU0/Vdisp;
    auto ve = VescTmp/Vdisp;

    auto u = sqrt(u0*u0+ksi*ksi+2*u0*ksi*cosPhi);


    auto sinTheta = (u != 0.0 ? ksi*sinPhi/u : 0);
    auto v = sqrt(u*u+ve*ve);
    auto n = RandomN(G);
    return MC::MCResult<vec3>(n*(v*Vdisp),sinTheta*v*sqrt(M_PI/2));
}

template <class Generator> 
/*MK generator of output nu'*/
inline MC::MCResult<vec3> NuOut(Generator G,const vec3& Vcm,const vec3&Nu,
						double Vesc,double mp,double mk,double deltaE = 0){
	double VcmN = Vcm.norm();
	
    vec3  n_v = Vcm/VcmN;

    double cosThetaVcm = n_v.z;
    double sinThetaVcm = sqrt(n_v.x*n_v.x+n_v.y*n_v.y);

    double cosPhiVcm = 1.0;
    double sinPhiVcm = 0.0;

    if(sinThetaVcm > 1e-10){
        cosPhiVcm = n_v.x/sinThetaVcm;
        sinPhiVcm =  n_v.y/sinThetaVcm;
    }

    vec3 n_1(cosThetaVcm*cosPhiVcm,cosThetaVcm*sinPhiVcm,-sinThetaVcm);
    vec3 n_2(-sinPhiVcm,cosPhiVcm,0);
	
    /*
    std::cout << n_1*n_1 << "\t" << n_1*n_2 << "\t" << n_1*n_v << std::endl;
    std::cout << n_2*n_1 << "\t" << n_2*n_2 << "\t" << n_2*n_v << std::endl;
    std::cout << n_v*n_1 << "\t" << n_v*n_2 << "\t" << n_v*n_v << std::endl<< std::endl;
    */

	double Nu1_squared = 	Nu.quad()-deltaE*2*mp/(mk*(mp+mk));
	if(Nu1_squared<=0.0)
		return MC::MCResult<vec3>(vec3(0,0,0),0);
	
	double Nu1 = sqrt(Nu1_squared);
	
	double cosTh1max = (Vesc*Vesc-Nu1_squared-VcmN*VcmN)/(2*VcmN*Nu1);
	
	if(cosTh1max <= -1)
		return MC::MCResult<vec3>(vec3(0,0,0),0);
	else if(cosTh1max >= 1){
		cosTh1max = 1;
	}
	
	double cosTh1 = (1+cosTh1max)*G()-1;
    double sinTh1 = sqrt(1.0-cosTh1*cosTh1);
    double phi1 = RandomPhi(G);


    return MC::MCResult<vec3>(Nu1*(n_v*cosTh1+n_1*sinTh1*cos(phi1)+n_2*sinTh1*sin(phi1)),
								0.5*(1.0+cosTh1max)*Nu1);
}

enum ScatteringType{
	ELASTIC,IONIZATION,MIGDAL
};
enum Target{
    ELECTRON,PROTON
};

template <typename ScatterFuncType,typename VescRFuncType,typename TempRFuncType,typename Generator>
inline MC::MCResult<vec3> Vout(double mk,double delta_mk,ScatteringType ST,Target T,
                              ScatterFuncType PhiFactor, VescRFuncType VescR,
                               TempRFuncType TempR, Generator G,
                               double Vdisp, double mU0){

    double mp = 0.938;
    double me = 0.52e-3;
    double factor = 1.0;

    //generate radius
    double r_nd = pow(G(),1.0/3.0);

    //gain escape velocity from redius
    double Vesc = VescR(r_nd);

    //random input velocity
    auto VelocityMk = Velocity(G,Vesc,Vdisp,mU0);
    auto V_wimp = VelocityMk.Result;
    factor *= VelocityMk.RemainDensity;


    double n_nd = 1.0;//TODO n_nd as a function of radius
    vec3 V1(0,0,0);//TODO: add thermal distribution of nuclei velocity

    //Vcm - is a vrlocity of momentum center
    vec3 Vcm = (V_wimp*mk + V1*mp)/(mp+mk);
    //Nu is input velocity of WIMP in cm coordinatesd
    vec3 Nu = mp/(mp+mk)*(V_wimp-V1);

    //Ecm - is kinetic energy in cm
    double E_cm = mk*(mp+mk)/mp*Nu.quad()/2;

    // this factor considers inelastic scaterring
    double Inelastic_rd = 1.0;
    double deltaE = delta_mk;
    int nMigdal;

    if(ST == IONIZATION){
        //generating ionization deltaE
        if(E_cm+delta_mk > Rd){
            //dE = max energy loss by ionization - min energyloss od ionization
            double dE = E_cm+delta_mk-Rd;
            deltaE = Rd + G()*dE;
            Inelastic_rd *= dE/Rd;
        }
        else{
            Inelastic_rd = 0;
        }
    }
    else if(ST == MIGDAL){
        double dnMx = nMax(E_cm+delta_mk)-2;
        if(dnMx > 0){
            nMigdal = 2 + G()*(int)(dnMx+1);
            deltaE = deltaEMgd(nMigdal);
        }
        else{
            deltaE = 0;
            nMigdal = -1;
            Inelastic_rd = 0;
        }
        Inelastic_rd *= (dnMx+1);
    }

    factor *= Inelastic_rd;

    // Generating out velocity
    auto Numk = NuOut(G,Vcm,Nu,Vesc,mp,mk,deltaE-delta_mk);
    vec3 Nu1 = Numk.Result;
    factor*=Numk.RemainDensity;

    // q - exchange momentum
    double q = mk*(Nu-Nu1).norm();

    // s - appeared in exponent in scalar product of states
    double s = q/(alpha*me);
    if(T == PROTON){
        s *= me/mp;
    }

    //Integrating factor from scalar product of states
    if(ST == IONIZATION && deltaE > Rd){
        factor *= IonFactor(s,phiMax(deltaE),dE_Rd);
    }
    else if(ST == MIGDAL && nMigdal >= 2){
        factor *= MigdalFactor(s,nMigdal);
    }

    //factor from matrix element
    factor *= PhiFactor(mk,q);

    /*
    if( (Nu1+Vcm).norm() > Vesc && factor > 1e-40){
        PVAR(Vcm.norm());
        PVAR(Nu1.norm());
        PVAR(Vesc);
        PVAR(Nu1+Vcm);
    }
    */

    return MC::MCResult<vec3>(Nu1+Vcm,factor);
}


template <typename FuncType>
MC::MCIntegral<double> SupressFactor(double mk,double delta_mk,ScatteringType ST,Target T,
					FuncType PhiFactor, 
					const BodyModel& BM,const std::string &element,
					size_t Nmk,
                    double Vdisp = U0,double mU0 = U0){

    //std::default_random_engine stdGen;
    //std::uniform_real_distribution<double> stdUniform(1.0,0.0);
    //auto G = [&stdGen,&stdUniform]()->double{return stdUniform(stdGen);};

    srand (time(NULL));
    auto G = [](){return  (1.0 +rand())/(RAND_MAX+1);};

	double sum = 0;
	double sum2 = 0;
	

	auto VescR = BM.UniformRadFunc("Vesc");
    auto TempR = 0;//BM.UniformRadFunc("Temp");

	for(size_t i=0;i<Nmk;++i){
		
        double factor = Vout(mk,delta_mk,ST,T,PhiFactor,VescR,TempR,G,Vdisp,mU0).RemainDensity;
		sum += factor/Nmk;
		sum2 += factor*factor/Nmk;
		
		
	}
    return MC::MCIntegral(sum,sqrt(sum2-sum*sum));
	
}





#endif
