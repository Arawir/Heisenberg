#ifndef MODEL
#define MODEL

#include "itensor/all.h"
#include "interface.h"

using namespace itensor;


Sweeps prepareSweepClass()
{
    auto sweeps = Sweeps(Args::global().getInt("sweeps"));
    sweeps.maxdim() = Args::global().getInt("maxDim");
    sweeps.mindim() = Args::global().getInt("minDim");
    sweeps.cutoff() = Args::global().getReal("cutoff");
    sweeps.niter() = Args::global().getReal("niter");
    return sweeps;
}

MPO heisenbergHamiltonian(SpinHalf &sites,
                       int L, double K)
{
    auto ampo = AutoMPO(sites);
    for(int j=1; j<L; j++){
        ampo += K,"Sp",j,"Sm",j+1;
        ampo += K,"Sm",j,"Sp",j+1;
        ampo += K,"Sz",j,"Sz",j+1;

    }

    if(Args::global().getBool("PBC")){
        ampo += K,"Sp",L,"Sm",1;
        ampo += K,"Sm",L,"Sp",1;
        ampo += K,"Sz",L,"Sz",1;
    }

    return toMPO(ampo);
}
/////////////////////////////////////////////////
void prepareObservables()
{
    ExpCon("Sz1:L") = [](const SpinHalf &sites){
        std::vector<MPO> out;

        for(int i=1; i<=sites.length(); i++){
            auto ampo = AutoMPO(sites);
            ampo += 1.0,"Sz",i;
            out.push_back( toMPO(ampo) );
        }

        return out;
    };
}

std::tuple<SpinHalf,MPS,MPO,Sweeps> prepareExpBasic()
{
    seedRNG(1);
    auto sites = SpinHalf( getI("L") );
    auto psi = prepareInitState(sites);
    auto H = heisenbergHamiltonian(sites,getI("L"),getD("K"));
    auto sweeps = prepareSweepClass();

    return std::make_tuple( sites,psi,H,sweeps );
}



#endif // MODEL

