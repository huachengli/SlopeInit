//
// Created by huacheng on 4/21/22.
//

#include "eos_cpp.h"

double ANEOS::InterpolateTD(double tTem,double tDen,int DataId) const
{
    auto data_ptr = reinterpret_cast<tabletype>(data);
    return ANEOSInterpolateTD(data_ptr,tTem,tDen,DataId);
}

ANEOS::ANEOS(const char * fname): BaseEOS()
{
    this->table = new ANEOSTable;
    this->data = reinterpret_cast<uintptr_t>(table);

    FILE * fp = fopen(fname,"r");
    assert(nullptr != fp);
    LoadANEOS(table,fp);
}

ANEOS::~ANEOS()
{
    UnAllocateANEOS(this->table);
    delete table;
}

double TillEOS::InterpolateTD(double tTem, double tDen, int DataId) const
{
    auto data_ptr = reinterpret_cast<tabletype>(data);
    return TillEOSInterpolateTD(data_ptr,tTem,tDen,DataId);
}

TillEOS::TillEOS(const char * fname)
{
    this->table = new TillotsonTable;
    this->data = reinterpret_cast<uintptr_t>(this->table);
    FILE * fp = fopen(fname,"r");
    assert(nullptr != fp);
    LoadTillEOS(this->table,fp);
}

TillEOS::~TillEOS()
{
    UnAllocateTillEOS(this->table);
    delete table;
}

double ANEOS::InterpolateTP(double tTem, double tPre, int DataId) const
{
    return ANEOSInterpolateTP(this->table,tTem,tPre,DataId);
}

double TillEOS::InterpolateTP(double tTem, double tPre, int DataId) const
{
    return TillEOSInterpolateTP(this->table,tTem,tPre,DataId);
}

TillEOS::TillEOS(uintptr_t _data): BaseEOS(_data){}
ANEOS::ANEOS(uintptr_t _data): BaseEOS(_data){}
BaseEOS::BaseEOS(uintptr_t _data): data(_data){}

double BaseEOS::RhoGrav(double z, double Pres, const double *Grav, const double *Temp) const
{
    double lGrav = (1-z)*Grav[0] + z*Grav[1];
    double lTemp = (1-z)*Temp[0] + z*Temp[1];
    lGrav = fabs(lGrav);
    double lRho  = InterpolateTP(lTemp, Pres, -1);
    return lRho*lGrav;
}

double BaseEOS::PresProfRK3(double *Pres, double *Grav, double *Temp, int steps, double dh) const
{
    double K1,K2,K3;
    for(int k=1;k<steps;k++)
    {
        K1 = RhoGrav(0.0 , Pres[k - 1]                 , Grav + k - 1, Temp + k - 1);
        K2 = RhoGrav(0.50, Pres[k - 1] + 0.50 * K1 * dh, Grav + k - 1, Temp + k - 1);
        K3 = RhoGrav(0.75, Pres[k - 1] + 0.75 * K2 * dh, Grav + k - 1, Temp + k - 1);
        Pres[k] = Pres[k-1] + (2.0*K1 + 3.0*K2 + 4.0*K3)*dh/9.0;
    }
    return Pres[steps-1];
}

BaseEOS::~BaseEOS() noexcept
{}