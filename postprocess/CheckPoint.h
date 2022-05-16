//
// Created by huacheng on 4/27/22.
//

#ifndef SLOPEINIT_CHECKPOINT_H
#define SLOPEINIT_CHECKPOINT_H
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "vtkMathUtilities.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkTesting.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"

#define DIM  3  // mesh dimension
#define DIM3 27
#define VPB  4  // nodes/vertexes per boundary
#define EPB  2  // neighbor elements per boundary
#define EPV  8  // neighbor elements per node/vertex
#define VPE  8  // neighbor nodes per element
#define BPE  6  // neighbor boundaries per element
#define NMT    4   // the max number of material types in model
#define NPHI   3   // the number of conservation variables in L3
#define NTHETA 3   // the number of variables determined by EOS
#define NPSI   10
#define NSBPB 4
#define NSEPE 8

struct Vertex;
struct Boundary;
struct Vnode{
    int  vId;
    double d;
};

struct Element
{
    double Center[DIM];
    double Scale[DIM];
    double Volume;
    double rVolume;
    double rdx_min;

    double Mass;
    double Density;
    double Damage;
    double Momentum[DIM];

    struct Vertex   *NeV[VPE];
    struct Boundary *NeB[BPE];
    struct Vnode    dCut[VPE];
    double CutRatio[VPE];

    unsigned short CutTag;
    double SubFace[VPE][DIM];
    double ABVArea;

    double GradFactor[DIM][VPE];
    double GradVel[DIM][DIM];
    double ArtificialPressure;
    double Strain[DIM][DIM];
    double DivVel;
    double InvSta2;
    double InvSta3;
    double sInvSte2; // sqrt of 2rd invariant of stress
    double sInvSte2_old; // used in plastic .
    double sInvSta2; // sqrt of 2rd invariant of strain
    // the invariants of strain

    // used in VOF related functions
    double GradVOF[NMT][DIM];
    double VOF[NMT];
    double Courant;

    /*
     * the variables used for L3 calculation
     * they conservation in calculation domain
     */

    double PhiL3[NMT][NPHI];
    double sPhiL3[NMT][NPHI];
    double PsiL3[NPSI];
    double sPsiL3[NPSI];

    // ThetaL3 is some variables updated by EOS
    double ThetaL3[NMT][NTHETA];
    double Temperature;
    double Pressure;
    double Csound;
    double ShearModule;

    double SubVolumeRatio[NSEPE];

    unsigned short State;
    double ShearStrength;
    double TensileStrength;
    unsigned short FailureState;

    double Viscosity;
    double VibPressure;
};

// just for load/write check point (different from SALEc)
class Mesh
{
private:
    int rank;
    char prefix[200];
public:
    Mesh() = delete;
    explicit Mesh(int _ne,int _rank);
    void LoadChk();
    void ExportChk(const char * out_prefix) const;
    void SetPrefix(const char *_prefix);
    void Resize(int _ne);
    void MergeVtu(vtkUnstructuredGrid * vtu_data);
    ~Mesh();

    struct Element * Elements;
    int ne;
};


#define MULTITHREAD_MERGE 1


#endif //SLOPEINIT_CHECKPOINT_H
