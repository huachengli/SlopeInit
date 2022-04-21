#include <iostream>
#include "src/InputParser.h"
#include "src/SlopeGeometry.h"
#include "src/eos_cpp.h"

void InitSlopeEos(SlopeInfo * pSlopeInfo)
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

int main()
{
    SlopeInfo pSlopeInfo;
    LoadSlopeInfo(&pSlopeInfo,"../SALEc.inp");
    InitSlopeEos(&pSlopeInfo);

    double Xc[8][3];
    GetTransformCorner(&pSlopeInfo,Xc);
    printf("%d\n",pSlopeInfo.noffset);
    return 0;
}
