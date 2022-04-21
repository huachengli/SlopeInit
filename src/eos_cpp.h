//
// Created by huacheng on 4/20/22.
//

#ifndef SLOPEINIT_EOS_CPP_H
#define SLOPEINIT_EOS_CPP_H

#include <cstdint>
#include "eos_state.h"

// c++ wrapper for eos table
class EosTable
{
public:
    virtual double InterpolateTD(double tTem,double tDen,int DataId) const = 0;
    virtual double InterpolateTP(double tTem,double tPre,int DataId) const = 0;
    EosTable() = default;
    EosTable(uintptr_t _data);
protected:
    uintptr_t data;
};

class ANEOS : public EosTable
{
public:
    ANEOS() = delete;
    ANEOS(const char * fname);
    ANEOS(uintptr_t _data);
    double InterpolateTD(double tTem,double tDen,int DataId) const override;
    double InterpolateTP(double tTem,double tPre,int DataId) const override;
    using tabletype = struct ANEOSTable *;
    ~ANEOS();

private:
    tabletype table;
};

class TillEOS : public EosTable
{
public:
    TillEOS() = delete;
    TillEOS(uintptr_t _data);
    TillEOS(const char * fname);
    double InterpolateTD(double tTem,double tDen,int DataId) const override;
    double InterpolateTP(double tTem,double tPre,int DataId) const override;
    using tabletype = struct TillotsonTable *;
    ~TillEOS();
private:
    tabletype table;
};




#endif //SLOPEINIT_EOS_CPP_H
