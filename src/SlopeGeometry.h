//
// Created by huacheng on 4/19/22.
//

#ifndef SLOPEINIT_SLOPEGEOMETRY_H
#define SLOPEINIT_SLOPEGEOMETRY_H

// cpp warpper for reader
#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define SLOPEDIM 3
#define SLOPENAMELEN 200
#define SLOPELAYER 5
typedef struct _SlopeInfo
{
    // calculation domain
    double length[SLOPEDIM];
    double ext[SLOPEDIM];
    int el[SLOPEDIM];
    int er[SLOPEDIM];
    double x0[SLOPEDIM];
    double dx[SLOPEDIM];
    double px0[SLOPEDIM];
    int nproc[SLOPEDIM];
    int npx[SLOPEDIM];
    int noffset;

    // position of slope
    char type[SLOPENAMELEN];
    double norm_deg[SLOPEDIM];
    double equation[SLOPEDIM+1];
    int nlayers;
    int layerId[SLOPELAYER];
    uintptr_t ldata[SLOPELAYER];

    // materials
    int nmat;
    char mname[SLOPELAYER][SLOPENAMELEN];
    char mpostfix[SLOPELAYER][SLOPENAMELEN];
    uintptr_t mdata[SLOPELAYER];

    // other condition
    double gravity;
    double nu[SLOPELAYER]; // poisson ratio
} SlopeInfo;

void LoadSlopeInfo(SlopeInfo * _s,const char * _fname);
double Sk(int k,double ext);
double Spacing(int k,double ext, int eL, int eR, int N);
void GetTransformCorner(SlopeInfo * _p,double (*Xi)[SLOPEDIM]);
void TransformInterpolate(double Xi[SLOPEDIM],double xl[SLOPEDIM],double (*Corner)[SLOPEDIM]);
#ifdef __cplusplus
}
#endif

#endif //SLOPEINIT_SLOPEGEOMETRY_H
