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
    char fname[200];
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
    char fname[200];
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

    _mesh->Resize(ChkElements);

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
    if(ne > 0)
        this->Elements = new Element[ne];
    else
        this->Elements = nullptr;
}

Mesh::~Mesh()
{
    delete [] this->Elements;
}

void Mesh::LoadChk()
{
    ChkLoad(this,this->rank,this->prefix);
}

void Mesh::ExportChk(const char *out_prefix) const
{
    ChkWrite(this,this->rank,out_prefix);
}

void Mesh::SetPrefix(const char *_prefix)
{
    strcpy(this->prefix,_prefix);
}

void Mesh::Resize(int _ne)
{
    if(_ne != this->ne)
    {
        delete [] this->Elements;
        this->ne = _ne;
        this->Elements = new Element[_ne];
    }
}

void Mesh::MergeVtu(vtkUnstructuredGrid *vtu_data)
{

    vtkNew<vtkGenericCell> local_gencell;

#ifdef MULTITHREAD_MERGE
    vtkNew<vtkUnstructuredGrid> local_vtu_data;
    local_vtu_data->ShallowCopy(vtu_data);
#else
    vtkUnstructuredGrid * local_vtu_data = vtu_data;
#endif

    vtkAbstractArray * stress_abstract_array = local_vtu_data->GetPointData()->GetAbstractArray("stress");
    auto * stress_array = vtkArrayDownCast<vtkFloatArray>(stress_abstract_array);
    local_vtu_data->BuildLocator();
    for(size_t k_element=0;k_element< this->ne;++k_element)
    {
        Element * element_ptr = this->Elements + k_element;

        double vacuum_fraction = element_ptr->VOF[0];
        if(vacuum_fraction > 1.0e-5) continue;

        double * element_center = element_ptr->Center;
        int subId;
        double cell_pcoord[24],cell_weight[24], cell_stress[9]={0.};

        vtkIdType element_cell_id = local_vtu_data->FindCell(element_center, nullptr, local_gencell,
                                                             -1, 1.0e-6, subId, cell_pcoord, cell_weight);

        if(element_cell_id < 0)
        {
            /*
             * element in projectile,
             */
            continue;
        }

        local_vtu_data->GetCell(element_cell_id, local_gencell);

        int n_cell_points = local_gencell->GetNumberOfPoints();
        int n_stress_compoents = stress_array->GetNumberOfComponents();
        assert(n_stress_compoents == 9);
        vtkNew<vtkFloatArray> cell_stress_array;
        cell_stress_array->SetNumberOfComponents(n_stress_compoents);
        cell_stress_array->SetNumberOfTuples(n_cell_points);
        stress_array->GetTuples(local_gencell->GetPointIds(), cell_stress_array);

        {
            for(vtkIdType i_point=0;i_point<n_cell_points;++i_point)
                for(vtkIdType k_component=0;k_component<n_stress_compoents;++k_component)
                {
                    cell_stress[k_component] = cell_stress_array->GetValue(i_point*n_stress_compoents + k_component)
                            * cell_weight[i_point];
                }
        }

        double cell_pressure = -(1./3.)*(cell_stress[0] + cell_stress[3] + cell_stress[8]);
        cell_stress[0] += cell_pressure;
        cell_stress[3] += cell_pressure;
        cell_stress[8] += cell_pressure;


        element_ptr->Pressure = cell_pressure;
        double * elememt_stress = element_ptr->PsiL3 + 1; // first compoent in PsiL3 is damage
        for(vtkIdType k_component=0;k_component<n_stress_compoents;++k_component)
        {
            elememt_stress[k_component] = cell_stress[k_component];
        }
    }
}

