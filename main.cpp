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
    tdvp(psi,H,im*dTime/sweeps.nsweep(),sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});
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
    Experiments("timeEvGs") = [](){
        auto [sites,psi,H,sweeps] = prepareExpBasic();
        ExpCon.setSites(sites);ExpCon("E") = H;

        ExpCon.addPoint("Prepare operator");
        auto ampo = AutoMPO(sites);
        ampo += 1,"S+",getI("L")/2;
        MPO OpSpinPlus = toMPO(ampo);

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);
        ExpCon.calc(psi,oMode::b,"rtime","mem","dim","E","Sz","Sz1:L");

        ExpCon.addPoint("Apply operator");
        auto psi2 = removeQNs(psi);
        auto psi3 = applyMPO(OpSpinPlus,psi2,{"Method=","DensityMatrix","MaxDim=",100,"Cutoff=",1E-20,"Normalize=",true});
        psi3.noPrime();

        ExpCon.addPoint("Starting TDVP");

        for(double time=0; time<=getD("maxtime")+getD("dtime")+0.001; time+=getD("dtime")){
            ExpCon.calc(psi3,oMode::b,"t:",time,"rtime","mem","dim","E","Sz","Sz1:L");
            tdvpStepWithBasisExtensionIfNeeded(psi3,H,getD("dtime"),sweeps);
        }
    };


    Params.add("K","double","1.0");
    Params.add("L","int","24");
    Params.add("PBC","bool","0");

    Params.add("dtime","double","0.1");
    Params.add("maxtime","double","196.0");

    Params.add("Silent","bool","1");
    Params.add("cutoff","double","1E-8");
    Params.add("sweeps","int","4");
    Params.add("minDim","int","1");
    Params.add("maxDim","int","100");
    Params.add("niter","int","30");
    Params.add("state","string","Up-Dn");
    Params.add("ConserveSz","bool","0");
    Params.add("ConserveQNs","bool","0");
    Params.add("exp","string","timeEvGs");

    Params.add("PBSenable","bool","0");
    Params.add("PBSjobid","int","0");

    Params.set(argc,argv);
    prepareObservables();
    Experiments.run();

    return 0;
}

