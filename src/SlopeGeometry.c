//
// Created by huacheng on 4/19/22.
//

#include "SlopeGeometry.h"
#include "InputParser.h"

void LoadSlopeInfo(SlopeInfo * _s,const char * _fname)
{
    InputFile * ifp = OpenInputFile(_fname);
    if(NULL == ifp)
    {
        fprintf(stdout,"cannot open %s\n",_fname);
        goto ERROR;
    }

    char MeshExtOpt[MaxStrLen];
    GetValueS(ifp,"mesh.ext",MeshExtOpt,"off");
    if(0== strcasecmp("on",MeshExtOpt))
    {
        _s->ext[0] = GetValueDk(ifp,"mesh.ex",0,"1.0");
        _s->ext[1] = GetValueDk(ifp,"mesh.ey",0,"1.0");
        _s->ext[2] = GetValueDk(ifp,"mesh.ez",0,"1.0");

        _s->el[0] = GetValueIk(ifp,"mesh.ex",1,"20");
        _s->el[1] = GetValueIk(ifp,"mesh.ey",1,"20");
        _s->el[2] = GetValueIk(ifp,"mesh.ez",1,"20");

        _s->er[0] = GetValueIk(ifp,"mesh.ex",2,"20");
        _s->er[1] = GetValueIk(ifp,"mesh.ey",2,"20");
        _s->er[2] = GetValueIk(ifp,"mesh.ez",2,"20");

    }
    else if(0 == strcasecmp("off",MeshExtOpt))
    {
        _s->ext[0] = _s->ext[1] = _s->ext[2] = 1.0;
        _s->el[0] = _s->el[1] = _s->el[2] = 0;
        _s->er[0] = _s->er[1] = _s->er[2] = 0;
    }
    else
    {
        fprintf(stdout,"unknown mesh.ext [%s]\n",MeshExtOpt);
        goto ERROR;
    }

    _s->nproc[0] = GetValueI(ifp,"processor.npgx","2");
    _s->nproc[1] = GetValueI(ifp,"processor.npgy","2");
    _s->nproc[2] = GetValueI(ifp,"processor.npgz","2");

    _s->npx[0] = GetValueI(ifp,"mesh.npx","16");
    _s->npx[1] = GetValueI(ifp,"mesh.npy","16");
    _s->npx[2] = GetValueI(ifp,"mesh.npz","16");

    _s->dx[0] = GetValueD(ifp,"mesh.dx","1.");
    _s->dx[1] = GetValueD(ifp,"mesh.dy","1.");
    _s->dx[2] = GetValueD(ifp,"mesh.dz","1.");

    _s->noffset = GetValueI(ifp,"mesh.noffset","2");

    _s->x0[0] = GetValueDk(ifp,"mesh.O",0,"0.5");
    _s->x0[1] = GetValueDk(ifp,"mesh.O",1,"0.5");
    _s->x0[2] = GetValueDk(ifp,"mesh.O",2,"0.5");

    GetValueSk(ifp,"target.toptype",_s->type,0,"slope");
    _s->norm_deg[0] = M_PI/180.0*GetValueDk(ifp,"target.toptype",1,"0.");
    _s->norm_deg[1] = M_PI/180.0*GetValueDk(ifp,"target.toptype",2,"0.");

    _s->equation[0] = sin(_s->norm_deg[0])*cos(_s->norm_deg[1]);
    _s->equation[1] = sin(_s->norm_deg[0])*sin(_s->norm_deg[1]);
    _s->equation[2] = cos(_s->norm_deg[0]);
    _s->equation[3] = GetValueD(ifp,"target.toplevel","0.0");

    _s->nmat = GetValueI(ifp,"material.nm","1");
    for(int k=0;k<_s->nmat;++k)
    {
        GetValueSk(ifp,"material.postfix",_s->mpostfix[k],k,"null");
        GetValueSk(ifp,"material.name",_s->mname[k],k,"null");
        _s->nu[k] = GetValueDk(ifp,"material.poisson",k,"0.25");
    }

    _s->nlayers = GetValueI(ifp,"target.number","1");
    for(int k=0;k<_s->nlayers;++k)
    {
        char LayerName[MaxStrLen];
        GetValueSk(ifp,"target.material",LayerName,k,"null");
        _s->layerId[k] = -1;
        for(int j=0;j<_s->nmat;++j)
        {
            if(0== strcasecmp(LayerName,_s->mname[j]))
            {
                _s->layerId[k] = j;
                break;
            }
        }
    }

    for(int k=0;k<SLOPEDIM;++k)
    {
        int Nxx = 2*_s->noffset + _s->nproc[k]*_s->npx[k];
        _s->length[k] = _s->dx[k]*Spacing(Nxx,_s->ext[k],_s->el[k],_s->er[k],Nxx);
        _s->px0[k] = -_s->x0[k]*Spacing(Nxx-1-_s->noffset,_s->ext[k],_s->el[k],_s->er[k],Nxx)
                + (_s->x0[k]-1.0)*Spacing(_s->noffset,_s->ext[k],_s->el[k],_s->er[k],Nxx);
        _s->px0[k] *= _s->dx[k];
    }

    _s->gravity = GetValueDk(ifp,"condition.gravity",1,"-1.0");

    EXIT:
    if(ifp) CloseInputFile(ifp);
    return;

    ERROR:
    if(ifp) CloseInputFile(ifp);
    exit(0);
}

double Sk(int k,double ext)
{
    return ext*(pow(ext,k) - 1.0)/(ext - 1.0);
}

double Spacing(int k,double ext, int eL, int eR, int N)
{
    double Xk = 0.0;
    if(k<eL)
        Xk = Sk(eL,ext) - Sk(eL-k,ext);
    else if(k<N-eR)
        Xk = Sk(eL,ext) + (k - eL);
    else
        Xk = Sk(eL,ext) + (N - eL - eR) + Sk(k-N+eR,ext);
    return Xk;
}

void GetTransformCorner(SlopeInfo * _p,double (*Xi)[SLOPEDIM])
{
    Xi[0][2] = Xi[1][2] = Xi[2][2] = Xi[3][2] = _p->px0[2];

    Xi[0][0]  = _p->px0[0] + _p->length[0];
    Xi[0][1]  = _p->px0[1] + _p->length[1];

    Xi[1][0]  = _p->px0[0];
    Xi[1][1]  = _p->px0[1] + _p->length[1];

    Xi[2][0]  = _p->px0[0];
    Xi[2][1]  = _p->px0[1];

    Xi[3][0]  = _p->px0[0] + _p->length[0];
    Xi[3][1]  = _p->px0[1];

    for(int k=0;k<4;++k)
    {
        Xi[k+4][0] = Xi[k][0];
        Xi[k+4][1] = Xi[k][1];
        Xi[k+4][2] = -(_p->equation[0]*Xi[k][0] + _p->equation[1]*Xi[k][1] + _p->equation[3])/_p->equation[2];
    }
}

// functions used in 8 points interpolate
// see the function array IpV and IpV_B
double Ip3d0(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d1(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d2(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d3(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d4(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d5(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d6(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d7(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[1]) * (1 + _x[2]);
}

double (*IpV[8])(const double *)  = {Ip3d0, Ip3d1, Ip3d2, Ip3d3, Ip3d4, Ip3d5, Ip3d6, Ip3d7};

void TransformInterpolate(double Xi[SLOPEDIM],double xl[SLOPEDIM],double (*Corner)[SLOPEDIM])
{
    Xi[0] = Xi[1] = Xi[2] = 0.;
    for(int j=0;j<7;j++)
    {
        for(int k=0;k<SLOPEDIM;++k)
        {
            Xi[k] += IpV[j](xl)*Corner[j][k];
        }
    }
}