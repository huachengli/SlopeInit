//
// Created by huacheng on 1/11/22.
//

#ifndef EOSTOOL_EOS_STATE_H
#define EOSTOOL_EOS_STATE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "InputParser.h"
#include <assert.h>

struct ANEOSTable
{
    /*
     * structured table, interpolate the value linearly
     * between data points
     */
    int nTem;
    int nDen;
    double Pnorm;
    double Dnorm;
    double Tnorm;
    double *yDen;
    double *xTem;
    double ***Data;
    // data is (density*energy,pressure,csound)
    // not same with iSALE2d
};

#define ANEOSDEN -1
#define ANEOSENG 0
#define ANEOSPRE 1
#define ANEOSCSD 2

#define WITHOUTCORRECT 0
#define CONSTSHIFT 1

struct StateReference
{
    /*
     * The parameters shared in one materials
     * which can be used in EOS or strength model.
     */
    int MaterialId;
    double VaporDen;
    double MeltDen;
    double GaCoef; // Constant to convert bulk modulus to shear modulus

    // parameters for simple ROCK strength model
    double Yint0;
    double Yfricint;
    double Ylimint;
    double Ydam0;
    double Yfricdam;
    double Ylimdam;
    double Yten0;

    double IvanA;
    double IvanB;
    double IvanC;

    double minPint;
    double minPdam;


    // parameters for Johnson-Cook strength model
    double JCA;
    double JCB;
    double JCN;
    double JCC;
    double JCM;
    double JCTREF;
    double JCMaxPlastic;

    // soft parameters
    double Asimon;
    double Csimon;
    double Tmelt0;
    double Tfrac;
    double Viscosity;

    // Acoustic Fluidization parameters
    double Toff;
    double Cvib;
    double VibMax;
    double Tdamp;
    double Pvlim;
    double Acvis;
    double GammaEta;
    double GammaBeta;
};

void AllocateANEOS(struct ANEOSTable *_t);

void UnAllocateANEOS(struct ANEOSTable *_t);

void LoadANEOS(struct ANEOSTable *_t, FILE *fp);

void ANEOSInitStateRef(struct ANEOSTable *_t, struct StateReference *_s);

double ANEOSInterpolateTD(struct ANEOSTable *_t, double tTem, double tDen, int DataId);

double ANEOSInterpolateTP(struct ANEOSTable *_t, double tTem, double tPre, int DataId);

void ANEOSWrite(struct ANEOSTable *_t, const char fname[], const char comment[]);

struct AirEOSTable
{
    int nEng;
    int nDen;
    double *xEng;
    double *yDen;
    double ***Data;
    // (All the variables is log10. form)
    // Data = log.10 (Pressure,Temperature,Speed of sound)
};

int BilinearInverse(const double xi[], const double yi[], double xl[]);

void AllocateAirEOS(struct AirEOSTable *_air);

void UnAllocateAirEOS(struct AirEOSTable *_air);

void LoadAirEOS(struct AirEOSTable *_t, FILE *fp);

struct TillotsonTable
{
    double TLRho0;
    double TLCv;
    double TLA;
    double TLB;
    double TLE0;
    double TLa;
    double TLb;
    double TLAlpha;
    double TLBeta;
    double TLEiv;
    double TLEcv;
    double TLTref;

    double *xDen;
    double *yEng;
    double StepDen;
    int nDen;
};

double TillPres(const struct TillotsonTable *_t, double EngIn, double RhoIn, double *Pres, double *Cs);

double TillTemp(const struct TillotsonTable *_t, double EngIn, double RhoIn);

void TillColdEnergy(struct TillotsonTable *_t);

double LoadTillEOS(struct TillotsonTable *_t, FILE *fp);

void TillInitStateRef(struct TillotsonTable *_t, struct StateReference *_s);

double TillEOSInterpolateTP(struct TillotsonTable *_t, double tTem, double tPre, int DataId);

double TillEOSInterpolateTD(struct TillotsonTable *_t, double tTem, double tDen, int DataId);

void UnAllocateTillEOS(struct TillotsonTable *_t);


// some addition function
double Wind(double x, double limitL, double limitR);

void Over(FILE *fp, int n);

double Max(double a, double b);

double Min(double a, double b);

void ANEOSLowDenCorrect(struct ANEOSTable *_t, double plimit);

#ifdef __cplusplus
}
#endif
#endif //EOSTOOL_EOS_STATE_H
