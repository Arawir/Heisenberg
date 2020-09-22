#include "model.h"
#include "interface.h"
#include "tdvp.h"
#include "basisextension.h"

void tdvpStepWithBasisExtensionIfNeeded(MPS &psi, MPO &H, double dTime, Sweeps &sweeps)
{
    if(maxLinkDim(psi)<10){
       std::vector<Real> epsilonK = {getD("cutoff"),getD("cutoff")};
       addBasis(psi,H,epsilonK,{"Cutoff",getD("cutoff"),"Method","DensityMatrix","KrylovOrd",2,"DoNormalize",true,"Quiet",true});
    }
    tdvp(psi,H,im*dTime,sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});
}

int main(int argc, char *argv[])
{

    Experiments("dmrg") = [](){
        auto [sites, psi, H, sweeps] = prepareExpBasic();
        ExpCon.setSites(sites); ExpCon("E") = H;
        ExpCon.calc(psi, oMode::b,"mem","dim","E");

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        ExpCon.calc(psi, oMode::b,"mem","dim","E","Sz1:L");
    };

    Experiments("timeEv") = [](){
        auto [sites,psi,H,sweeps] = prepareExpBasic();
        ExpCon.setSites(sites);ExpCon("E") = H;

        ExpCon.addPoint("Starting TDVP");

        for(double time=0; time<=getD("maxtime")+getD("dtime")+0.001; time+=getD("dtime")){
            ExpCon.calc(psi,oMode::b,"t:",time,"rtime","mem","dim","E","Sz1:L");
            tdvpStepWithBasisExtensionIfNeeded(psi,H,getD("dtime"),sweeps);
        }
    };


    Params.add("K","double","0.042857143");
    Params.add("L","int","24");
    Params.add("PBC","bool","0");

    Params.add("dtime","double","0.042857143");
    Params.add("maxtime","double","96.0");

    Params.add("Silent","bool","1");
    Params.add("cutoff","double","1E-8");
    Params.add("sweeps","int","4");
    Params.add("minDim","int","10");
    Params.add("maxDim","int","200");
    Params.add("niter","int","30");
    Params.add("state","string","L/2*Up-L/2*Dn");
    Params.add("ConserveSz","bool","0");
    Params.add("exp","string","timeEv");

    Params.add("PBSenable","bool","0");
    Params.add("PBSjobid","int","0");

    Params.set(argc,argv);
    prepareObservables();
    Experiments.run();

    return 0;
}

