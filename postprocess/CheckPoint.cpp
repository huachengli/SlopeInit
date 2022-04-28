//
// Created by huacheng on 4/27/22.
//

#include "CheckPoint.h"

void ScalerMove(double *_a, const double *_b, double _p)
{
    // _a = p * _b
    for (int k = 0; k < DIM; k++)
        _a[k] = _p * _b[k];
}

void ChkWrite(const Mesh * _mesh,int _xrank, const char _chkPrefix[])
{
    char fname[50];
    sprintf(fname,"%s.proc%04d.chk",_chkPrefix,_xrank);
    FILE *fp = fopen(fname,"wb");
    fwrite(&(_mesh->ne),sizeof(int),1,fp);
    size_t dsize = sizeof(double);
    char ChkNum = 'K';
    for(int k=0;k<_mesh->ne;k++)
    {
        struct Element * ie = _mesh->Elements + k;
        double * VibVelocity = ie->PsiL3 + 8;
        fwrite(&k,sizeof(int),1,fp);
        fwrite(&(ie->PsiL3[0]),dsize,NPSI,fp);
        fwrite(&(ie->Momentum[0]),dsize,DIM,fp);
        fwrite(&(ie->VOF[0]),dsize,NMT,fp);
        fwrite(&(ie->PhiL3[0][0]),dsize,NMT*NPHI,fp);
        fwrite(&(ie->VibPressure),dsize,1,fp);
        fwrite(&(*VibVelocity),dsize,1,fp);
        fwrite(&(ie->Center[0]),dsize,DIM,fp);
        fwrite(&(ChkNum),sizeof(char),1,fp);
    }
    fclose(fp);
}

void ChkLoad(Mesh * _mesh,int _xrank, const char _chkPrefix[])
{
    char fname[50];
    sprintf(fname,"%s.proc%04d.chk",_chkPrefix,_xrank);
    FILE *fp = fopen(fname,"rb");
    if(NULL == fp)
    {
        fprintf(stdout,"can not open chkfile: %s\n",fname);
        exit(0);
    }

    size_t dsize = sizeof(double);
    char ChkNum = 'K';
    int ChkElements = 0;

    if(1!=fread(&(ChkElements),sizeof(int),1,fp))
    {
        fprintf(stdout,"Error in read number of elements in check point files");
        exit(0);
    }

    int ChkIndex = 0;
    for(int k=0;k<_mesh->ne;k++)
    {
        struct Element * ie = _mesh->Elements + k;
        double * VibVelocity = ie->PsiL3 + 8;
        size_t rChunksSize = 0;
        rChunksSize += fread(&ChkIndex,sizeof(int),1,fp);
        rChunksSize += fread(&(ie->PsiL3[0]),dsize,NPSI,fp);
        rChunksSize += fread(&(ie->Momentum[0]),dsize,DIM,fp);
        rChunksSize += fread(&(ie->VOF[0]),dsize,NMT,fp);
        rChunksSize += fread(&(ie->PhiL3[0][0]),dsize,NMT*NPHI,fp);
        rChunksSize += fread(&(ie->VibPressure),dsize,1,fp);
        rChunksSize += fread(&(*VibVelocity),dsize,1,fp);
        rChunksSize += fread(&(ie->Center[0]),dsize,DIM,fp);
        rChunksSize += fread(&ChkNum,sizeof(char),1,fp);

        const size_t rChunkSizeExpected = (4+NPSI+DIM+NMT*(NPHI+1)) + DIM;
        if(ChkNum!='K' || ChkIndex!=k || rChunksSize != rChunkSizeExpected)
        {
            fprintf(stdout,"CheNum or ChkIndex unexpected in Chkfile\n "
                           "ChunkSize=%ld(%ld expected),ChkNum=%c(K expected)\n",
                    rChunksSize,rChunkSizeExpected,ChkNum);
            exit(0);
        }

        ie->Density = 0.0;
        for(int i=0;i<NMT;i++)
        {
            ScalerMove(ie->sPhiL3[i],ie->PhiL3[i],ie->VOF[i]);
            ie->Density += ie->sPhiL3[i][1];
        }
        ie->Mass = ie->Density*ie->Volume;

        ScalerMove(ie->sPsiL3,ie->PsiL3,ie->Mass);
        ie->Damage = ie->PsiL3[0];
    }
    fclose(fp);
}

Mesh::Mesh(int _ne,int _rank):ne(_ne),rank(_rank)
{
    this->Elements = new Element[ne];
}

Mesh::~Mesh()
{
    delete [] this->Elements;
}

void Mesh::LoadChk()
{
    ChkLoad(this,this->rank,this->prefix);
}

void Mesh::WriteChk(const char *out_prefix) const
{
    ChkWrite(this,this->rank,out_prefix);
}

void Mesh::SetPrefix(const char *_prefix)
{
    strcpy(this->prefix,_prefix);
}