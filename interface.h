#ifndef INTERFACE
#define INTERFACE

#include "itensor/all.h"
#include <ctime>
#include <functional>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

using namespace itensor;

//////SET THIS////////////////////
using SiteSetType=SpinHalf;
//////////////////////////////////

typedef std::complex<double> cpx;
#define im std::complex<double>{0.0,1.0}


inline void massert(std::string info)
{
    std::cerr << info << std::endl;
    assert(false);
}
inline double getD(std::string data)
{
    return Args::global().getReal(data);
}
inline int getI(std::string data)
{
    return Args::global().getInt(data);
}
inline bool getB(std::string data)
{
    return Args::global().getBool(data);
}
inline std::string getS(std::string data)
{
    return Args::global().getString(data);
}

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    auto pos = result.find("\n");
    return result.substr(0,pos-1);
}


template <class T>
class BasicObjectList{
protected:
    std::vector<T*> objects;

public:
    BasicObjectList() : objects{} { }
    ~BasicObjectList()
    {
        for(auto& object : objects){
            delete object;
        }
    }

    T& operator () (std::string name)
    {
        if( exists(name) ){ return *get(name); }
        objects.push_back( new T{name} );
        return *objects.back();
    }

protected:
    T* get(std::string name)
    {
        for(auto &object : objects){
            if(object->name == name) return object;
        }
        massert("ERROR: Wrong object name ("+name+")!");
        return nullptr;
    }
    bool exists(std::string name)
    {
        for(auto &object : objects){
            if(object->name == name) return true;
        }
        return false;
    }


};
/////////////////////////////
struct Param{
    std::string name;
    std::string type;
    std::string value;
    bool wasSet=false;

//    std::vector<std::string> values;
//    std::vector<std::string> operators;

//    void setValue(std::string text)
//    {
//        values.push_back("");
//        for(uint i=0; i<text.length(); i++){
//            if( (text[i]=='+')||(text[i]=='-')||(text[i]=='*')||(text[i]=='/')||(text[i]=='^') ){
//                operators.push_back( std::string(text[i]) );
//            } else {
//                values.back().append(text[i]);
//            }
//        }
//    }



    int getInt(){
        if(type!="int"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return atoi(value.c_str());
    }

    double getDouble(){
        if(type!="double"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return atof(value.c_str());
    }

    bool getBool(){
        if(type!="bool"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return (bool)atoi(value.c_str());
    }
};

class Parameters : public BasicObjectList<Param>{
public:
    Parameters(){ }
    ~Parameters() = default;

    void add(std::string name, std::string type, std::string value){
        objects.push_back( new Param{name, type, value} );
    }

    void set(int argc, char *argv[]){
        for(int i=1; i<argc; i++){
            set(argv[i]);
        }
        addArgs();
        write();
    }

    bool getBool(std::string name){
        return get(name)->getBool();
    }
    int getInt(std::string name){
        return get(name)->getInt();
    }
    double getDouble(std::string name){
        return get(name)->getDouble();
    }

    void write(){
        std::cout << "---------------------------- System parameters ---------------------------" << std::endl;
        for(auto& param : objects){
            std::cout << "  " << param->name << " -> ";
            std::cout << "  " << param->value << " ";
            if(!param->wasSet){ std::cout << "(default)"; }
            std::cout << std::endl;
        }
    }

private:
    void set(std::string data){
        size_t pos = data.find("=");
        get( data.substr(0,pos) )->value = data.substr(pos+1);
        get( data.substr(0,pos) )->wasSet = true;
    }

    void addArgs(){
        for(auto &param : objects){
            if(param->type=="bool"){
                Args::global().add(param->name, param->getBool());
            } else if(param->type=="string"){
                Args::global().add(param->name, param->value);
            } else if(param->type=="int"){
                Args::global().add(param->name, param->getInt());
            } else if(param->type=="double"){
                Args::global().add(param->name, param->getDouble());
            }
        }
    }
} Params;
//////////////////////////////

std::vector<std::string> separateSubBlocs(std::string str)
{
    std::vector<std::string> out{};
    int i=0;
    out.push_back("");
    while(str[i]!='\0'){
        if(str[i]=='-'){
            out.push_back("");
        } else {
            out.back().append(std::string{str[i]});
        }
        i++;
    }
    return out;
}
int identifyMultiplier(std::string data)
{
    if(data.find("L/")!=std::string::npos){
        return getI("L")/atoi( data.substr(data.find("L/")+2).c_str() );
    }
    return atoi(data.c_str());
}
std::tuple<int,std::string> separateMultiplier(std::string data)
{
    int multiplier=1;
    std::string state=data;
    if(data.find("*")!=std::string::npos){
        multiplier = identifyMultiplier( data.substr(0,data.find("*")) );
        state = data.substr(data.find("*")+1);
    }

    return std::make_tuple( multiplier,state );
}
std::vector<std::string> expandSubBlocs(std::vector<std::string> subBlocs)
{
    std::vector<std::string> out{};

    for(auto subBlock : subBlocs){
        auto [n,state] = separateMultiplier(subBlock);
        for(int i=0; i<n; i++){
            out.push_back(state);
        }
    }
    return out;
}
std::vector<std::string> parseInitState(std::string str)
{
    int L = Params.getInt("L");

    std::vector<std::string> subBlocs = separateSubBlocs(str);
    std::vector<std::string> state = expandSubBlocs(subBlocs);


    if(state.size()==0){
        std::cerr << "  Init state is not delcared!" << std::endl;
    }else if(state.size()<(uint)L){
        uint templateLength = state.size()-1;
        std::cout << "  Init state deduced" << std::endl;
        for(int i=state.size()-1; i<L; i++){
            state.push_back( state.at(i-templateLength) );
        }
    }
    return state;
}
itensor::MPS prepareInitState(auto &sites)
{
    if( (Args::global().getString("state")=="random")||(Args::global().getString("state")=="r")){
        auto state = itensor::InitState(sites);
        auto s2 = itensor::randomMPS(state);
        s2.orthogonalize();
        s2.normalize();
        return s2;
    }

    std::vector<std::string> initState = parseInitState(Args::global().getString("state"));
    int L = Params.getInt("L");

    auto state = itensor::InitState(sites);
    for(int i=1;i<=L;i++){
        state.set(i,initState[i-1]);
        std::cout << initState[i-1] << "-";
    }
    std::cout << std::endl;

    return itensor::MPS(state);
}

///////////////////////////////////////////
struct Observable
{
    std::string name;
    std::vector<MPO> matrix;
    bool wasGenerated=false;
    bool useMultMpos = false;
    std::function<MPO(const SiteSetType &sites)> generateMPO;
    std::function<std::vector<MPO>(const SiteSetType &sites)> generateMPOs;
    SiteSetType *sites = nullptr;

    void operator = (std::function<MPO(const SiteSetType &sites)> generate)
    {
        generateMPO = generate;
    }

    void operator = (std::function<std::vector<MPO>(const SiteSetType &sites)> generate)
    {
        useMultMpos=true;
        generateMPOs = generate;
    }

    void operator =  (MPO &nMatrix)
    {
        matrix.push_back( nMatrix );
        wasGenerated=true;
    }

    void generateIfNeeded()
    {
        if(wasGenerated==false){
            wasGenerated=true;
            if(!useMultMpos){
                matrix.push_back( generateMPO(*sites) );
            } else {
                matrix = generateMPOs(*sites);
            }
        }
    }

    std::vector<cpx> expectedValue(MPS psi)
    {
        std::vector<cpx> out;
        for(auto &m : matrix){
            out.push_back( innerC(psi,m,psi) );
        }
        return out;
    }
};

enum class oMode{
    a, b, c
};

class ObservableContainer : public BasicObjectList<Observable>
{
protected:
    oMode mode = oMode::a;
    clock_t  time0, timeLast;

public:
    ObservableContainer() :
        time0{ clock() }
        , timeLast{ clock() }
    { }
    ~ObservableContainer() = default;

    template<typename T,typename... Args>
    void calc(MPS &psi, T val, Args... args)
    {
        std::cout << " ";
        writeObservableValue(val,psi);
        if(mode==oMode::c){ std::cout << std::endl; }
        calc(psi, args...) ;
    }

    template<typename... Args>
    void calc(MPS &psi, oMode val, Args... args)
    {
        mode = val;
        calc(psi, args...) ;
    }

    void calc(MPS &psi)
    {
        std::cout << std::endl;
    }

    void setSites(SiteSetType &sites)
    {
        for(auto& obs : objects){
            obs->sites = &sites;
        }
    }


private:
    void initObservablesIfNeeded()
    {
        for(auto& obs : objects){
            obs->generateIfNeeded();
        }
    }

    void writeObservableValue(std::string name, MPS &psi)
    {
        if( exists(name) == true){
            get(name)->generateIfNeeded();
            if(mode==oMode::b || mode==oMode::c){ std::cout << name << "= "; }
            for(auto val : get(name)->expectedValue(psi)){
                std::cout << val.real() << " ";
            }
        } else if(name=="dim"){
            std::cout << "dim= " << maxLinkDim(psi);
        } else if(name=="mem"){
            if(getB("PBSenable")==true){
                std::string command = "qstat -fx " + std::to_string(getI("PBSjobid")) + " | grep  'vmem' | awk '{print $3}'";
                std::cout << "mem= " << exec(command.c_str());
            } else {
                std::cout << "mem= " << "unknown";
            }
        }  else if(name=="rtime"){
            std::cout << "rtime= " << (clock()-time0)/(double)CLOCKS_PER_SEC;
        } else {
            std::cout << name;
        }
    }
    void writeObservableValue(double val, MPS &psi){ std::cout << val; }
    void writeObservableValue(cpx val, MPS &psi){ std::cout << val; }
    void writeObservableValue(int val, MPS &psi){ std::cout << val; }
};

/////////////////////////////////////////////////
class Controller : public ObservableContainer
{
public:
    Controller() : ObservableContainer{ } { }
    ~Controller() = default;

    void addPoint(std::string text)
    {
        std::cout << "--------------------";
        std::cout << text;
        for(int i=0; i<60-(int)text.size(); i++){
           std::cout << "-";
        }

        std::cout << "(" << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s],"
                  << (clock()-timeLast)/(double)CLOCKS_PER_SEC << " [s])";
        std::cout << std::endl;
        timeLast = clock();
    }
} ExpCon;

////////////////////////////////////////////////
struct Experiment
{
    std::string name;
    std::function<void(void)> experiment;

    Experiment(std::string nName) : name{nName} { }

    void operator = (std::function<void(void)> exp)
    {
        experiment = exp;
    }
};

class ExperimentsClass : public BasicObjectList<Experiment>
{
public:
    ExperimentsClass() { }
    ~ExperimentsClass() = default;

    void run()
    {
        ExpCon.addPoint("Start");
        get(getS("exp"))->experiment();
        ExpCon.addPoint("Finish");
    }

} Experiments;

#endif // INTERFACE



