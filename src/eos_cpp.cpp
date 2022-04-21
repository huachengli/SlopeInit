//
// Created by huacheng on 4/21/22.
//

#include "eos_cpp.h"

double ANEOS::InterpolateTD(double tTem,double tDen,int DataId) const
{
    auto data_ptr = reinterpret_cast<tabletype>(data);
    return ANEOSInterpolateTD(data_ptr,tTem,tDen,DataId);
}

ANEOS::ANEOS(const char * fname): EosTable()
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

TillEOS::TillEOS(uintptr_t _data): EosTable(_data){};
ANEOS::ANEOS(uintptr_t _data): EosTable(_data){};
EosTable::EosTable(uintptr_t _data):data(_data){};
