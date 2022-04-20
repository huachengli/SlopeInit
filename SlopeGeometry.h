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
typedef struct _SlopeInfo
{
    // calculation domain
    double length[SLOPEDIM];
    double ext[SLOPEDIM];
    int el[SLOPEDIM];
    int er[SLOPEDIM];
    double x0[SLOPEDIM];
    double nproc[SLOPEDIM];
    double npx[SLOPEDIM];
    double dx[SLOPEDIM];

    // position of slope
    char type[SLOPENAMELEN];
    double norm_deg[SLOPEDIM];
    double equation[SLOPEDIM+1];
    int nlayers;
    uintptr_t ldata;

    // materials
    int nmat;
    uintptr_t mdata;
} SlopeInfo;

void LoadSlopeInfo(SlopeInfo * _s,const char * _fname);

#ifdef __cplusplus
}
#endif
#endif //SLOPEINIT_SLOPEGEOMETRY_H
