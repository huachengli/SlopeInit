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

    _s->x0[0] = GetValueDk(ifp,"mesh.O",0,"0.5");
    _s->x0[1] = GetValueDk(ifp,"mesh.O",1,"0.5");
    _s->x0[2] = GetValueDk(ifp,"mesh.O",2,"0.5");

    GetValueSk(ifp,"target.toptype",_s->type,0,"slope");
    _s->norm_deg[0] = GetValueDk(ifp,"target.toptype",1,"0.");
    _s->norm_deg[1] = GetValueDk(ifp,"target.toptype",2,"0.");


    EXIT:
    if(ifp) CloseInputFile(ifp);
    return;

    ERROR:
    if(ifp) CloseInputFile(ifp);
    exit(0);

}