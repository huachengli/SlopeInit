//
// Created by huacheng on 4/21/22.
//

#include "SlopeSlove.h"

void SlopeEquationData::InitSlopeEos(SlopeInfo * pSlopeInfo)
{
    for(int k=0;k<pSlopeInfo->nmat;++k)
    {
        char eos_file[SLOPENAMELEN*2];
        sprintf(eos_file,"../eos/%s.%s",pSlopeInfo->mname[k],pSlopeInfo->mpostfix[k]);
        if(0== strcasecmp("aneos",pSlopeInfo->mpostfix[k]))
        {
            auto ANEOS_ptr = new ANEOS(eos_file);
            pSlopeInfo->mdata[k] = reinterpret_cast<uintptr_t>(ANEOS_ptr);

        } else if(0 == strcasecmp("tillotson",pSlopeInfo->mpostfix[k]))
        {
            auto TillEOS_ptr = new TillEOS(eos_file);
            pSlopeInfo->mdata[k] = reinterpret_cast<uintptr_t>(TillEOS_ptr);
        } else if(0 == strcasecmp("vacuum_",pSlopeInfo->mname[k]))
        {
            pSlopeInfo->mdata[k] = 0;
        } else
        {
            fprintf(stdout,"unknown eos type [%s]\n",eos_file);
            exit(0);
        }
    }
}

void SlopeEquationData::InitSlopeProfile(SlopeInfo * _s)
{
    _s->dh = 0.2*_s->dx[2];
    _s->nstep = 1 + (int)(_s->sum_depth[_s->nlayers-1]/_s->dh);

    _s->Grav  = (double *) malloc(sizeof(double)*_s->nstep);
    _s->Pre   = (double *) malloc(sizeof(double)*_s->nstep);
    _s->Den   = (double *) malloc(sizeof(double)*_s->nstep);
    _s->Tem   = (double *) malloc(sizeof(double)*_s->nstep);
    _s->Cs    = (double *) malloc(sizeof(double)*_s->nstep);
    _s->Nu    = (double *) malloc(sizeof(double)*_s->nstep);

    // set the initial values for Grav and Tem
    // cosnt grav and temperature
    for(int k=0;k<_s->nstep;++k)
    {
        _s->Grav[k] = _s->gravity;
        _s->Tem[k] = _s->temperature;
    }

    for(int k=0;k<_s->nlayers;++k)
    {
        int layer_step = (int)(_s->depth[k]/_s->dh + 1);
        int layer_init = (int)((_s->sum_depth[k]-_s->depth[k])/_s->dh);

        auto layer_eos = reinterpret_cast<BaseEOS *>(_s->mdata[_s->layerId[k]]);
        layer_eos->PresProfRK3(_s->Pre+layer_init,_s->Grav+layer_init,_s->Tem+layer_init,layer_step,_s->dh);
        for(int j=0;j<layer_step;j++)
        {
            _s->Den[j+layer_init] = layer_eos->InterpolateTP(_s->Tem[j+layer_init],_s->Pre[j+layer_init],-1);
            _s->Cs[j+layer_init] = layer_eos->InterpolateTP(_s->Tem[j+layer_init],_s->Pre[j+layer_init],2);
            _s->Nu[j+layer_init] = _s->nu[_s->layerId[k]];
        }
    }
}
