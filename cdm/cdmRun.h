#ifndef RUN_H
#define RUN_H

#include <QtSql>
#include <QtGui>
#include "QsLog.h"

#include <complex>

#ifdef OS
#if OS == LINUX
#include <unistd.h>
#elif OS == WIN32
#include <windows.h>
#endif
#endif

using namespace std;

typedef complex<double> dcmplx;

class Run
{
   public:
   Run(QString _runname);
   ~Run();

   int checkAvailableMemorySize();
   int allocateVectorsMemory(int nmax, int ntheta, int nphi, int nfft2d, int obj_num, int n1m, int nplanm, int nmatim, int nbs, int nxm, int nym, int nzm);
   void cleanVectorsMemory();

   double* getIncidentField();
   double* getLocalField();
   double* getMacroscopicField();
   double* getXc();
   double* getYc();
   double* getZc();
   double* getXcWF();
   double* getYcWF();
   double* getZcWF();
   double* getThetaField();
   double* getPhiField();
   double* getPoyntingField();
   double* getForceX();
   double* getForceY();
   double* getForceZ();
   double* getForceXMulti();
   double* getForceYMulti();
   double* getForceZMulti();
   double* getTorqueX();
   double* getTorqueY();
   double* getTorqueZ();
   double* getTorqueXMulti();
   double* getTorqueYMulti();
   double* getTorqueZMulti();
   dcmplx* getIncidentFieldX();
   dcmplx* getLocalFieldX();
   dcmplx* getMacroscopicFieldX();
   dcmplx* getIncidentFieldY();
   dcmplx* getLocalFieldY();
   dcmplx* getMacroscopicFieldY();
   dcmplx* getIncidentFieldZ();
   dcmplx* getLocalFieldZ();
   dcmplx* getMacroscopicFieldZ();
   dcmplx* getPolarisaField();
   dcmplx* getEpsilonField();
   double* getXY();
   double* getKXY();
   dcmplx* getEimageX();
   dcmplx* getEimageY();
   dcmplx* getEimageZ();
   dcmplx* getEfourierX();
   dcmplx* getEfourierY();
   dcmplx* getEfourierZ();
   dcmplx* getEfourierincX();
   dcmplx* getEfourierincY();
   dcmplx* getEfourierincZ();
   dcmplx* getEimageincX();
   dcmplx* getEimageincY();
   dcmplx* getEimageincZ();
   dcmplx* getEimageXneg();
   dcmplx* getEimageYneg();
   dcmplx* getEimageZneg();
   dcmplx* getEfourierXneg();
   dcmplx* getEfourierYneg();
   dcmplx* getEfourierZneg();
   dcmplx* getEfourierincXneg();
   dcmplx* getEfourierincYneg();
   dcmplx* getEfourierincZneg();
   dcmplx* getEimageincXneg();
   dcmplx* getEimageincYneg();
   dcmplx* getEimageincZneg();
   dcmplx* getMasque();
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
   dcmplx* getFF();
   dcmplx* getFF0();
   dcmplx* getFFloc();
   dcmplx* getxr();
   dcmplx* getxi();
   dcmplx* getwrk();
   dcmplx* getmatrange();
   dcmplx* geta11();
   dcmplx* geta12();
   dcmplx* geta13();
   dcmplx* geta22();
   dcmplx* geta23();
   dcmplx* geta31();
   dcmplx* geta32();
   dcmplx* geta33();
   dcmplx* getb11();
   dcmplx* getb12();
   dcmplx* getb13();
   dcmplx* getb22();
   dcmplx* getb23();
   dcmplx* getb31();
   dcmplx* getb32();
   dcmplx* getb33();
   int *getmatindice();
   int *getmatind();
   int *getmatindplan();
   dcmplx* getEdiffkzpos();
   dcmplx* getEdiffkzneg();
   int* getTabdip();
   int* getTabmulti();
   int* getTabzn();
   int* getTabfft2();

   QString  getName();
   int      getObjectSubunits();
   int      getMeshSubunits();
   int      getNmaxpp();
   double   getMeshSize();
   double   getLambda10n();
   double   getK0();
   double   getToleranceObtained();

   int      getNumberofAx1();
   int      getNumberofAx2();
   double   getReflectivity();
   double   getTransmittivity();
   double   getAbsorptivity();
   double   getExtinctionCrossection();
   double   getAbsorbingCrossection();
   double   getScatteringCrossection();
   double   getScatteringCrossectionWithIntegration();
   double   getScatteringAssymetricParam();

   double   getIrra();
   dcmplx   getE0();

   double   getOpticalForcex();
   double   getOpticalForcey();
   double   getOpticalForcez();
   double   getOpticalForceModulus();
   double   getOpticalTorquex();
   double   getOpticalTorquey();
   double   getOpticalTorquez();
   double   getOpticalTorqueModulus();

   void setName(QString _runname);
   void setObjectSubunits(int _objectsubunits);
   void setMeshSubunits(int _meshsubunits);
   void setNmaxpp(int _nmaxpp);
   void setMeshSize(double _meshsize);
   void setLambda10n(double _lambda10n);
   void setK0(double _k0);
   void setToleranceObtained(double _toleranceobtained);

   void setNumberofAx1(int _numberofax1);
   void setNumberofAx2(int _numberofax2);
   void setReflectivity(double _reflectivity);
   void setTransmittivity(double _transmittivity);
   void setAbsorptivity(double _absorptivity);
   void setExtinctionCrossection(double _extinctioncrosssection);
   void setAbsorbingCrossection(double _absorbingcrosssection);
   void setScatteringCrossection(double _scatteringcrosssection);
   void setScatteringCrossectionWithIntegration(double _scatteringcrosssectionwithintegration);
   void setScatteringAssymetricParam(double _scatteringassymetricparam);

   void setIrra( double _irra);
   void setE0( dcmplx E0 );

   void setOpticalForcex(double _opticalforcex);
   void setOpticalForcey(double _opticalforcey);
   void setOpticalForcez(double _opticalforcez);
   void setOpticalForceModulus(double _opticalforcemodulus);
   void setOpticalTorquex(double _opticaltorquex);
   void setOpticalTorquey(double _opticaltorquey);
   void setOpticalTorquez(double _opticaltorquez);
   void setOpticalTorqueModulus(double _opticaltorquemodulus);

  

   QStringList beamList;
   QStringList objectList;
   QStringList anisotropyList;
   QStringList materialList;

  private:

    QString dbpath;

    double *incidentfield;
    double *localfield;
    double *macroscopicfield;
    double *xc, *yc, *zc, *xcwf, *ycwf, *zcwf, *kxy, *xy;
    double *thetafield, *phifield, *poyntingfield;
    double *forcex, *forcey, *forcez;  
    double *forcexmulti, *forceymulti, *forcezmulti;  
    double *torquex, *torquey, *torquez;
    double *torquexmulti, *torqueymulti, *torquezmulti;
    dcmplx *incidentfieldx, *incidentfieldy, *incidentfieldz;
    dcmplx *localfieldx, *localfieldy, *localfieldz;
    dcmplx *macroscopicfieldx, *macroscopicfieldy, *macroscopicfieldz;
    dcmplx *polarisafield, *epsilonfield;
    dcmplx *eimagex, *eimagey, *eimagez, *efourierx, *efouriery, *efourierz;
    dcmplx *efourierincx, *efourierincy, *efourierincz;
    dcmplx *eimageincx, *eimageincy, *eimageincz;
    dcmplx *eimagexneg, *eimageyneg, *eimagezneg, *efourierxneg, *efourieryneg, *efourierzneg;
    dcmplx *efourierincxneg, *efourierincyneg, *efourierinczneg;
    dcmplx *eimageincxneg, *eimageincyneg, *eimageinczneg;
    dcmplx *masque;
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
    dcmplx *FF, *FF0, *FFloc, *xr, *xi;
    dcmplx *wrk;
    int *matindice, *matind, *matindplan;
    dcmplx *matrange;
    dcmplx *a11, *a12, *a13, *a22, *a23, *a31, *a32, *a33;
    dcmplx *b11, *b12, *b13, *b22, *b23, *b31, *b32, *b33;
    dcmplx *Ediffkzpos, *Ediffkzneg;
    int *Tabdip, *Tabmulti, *Tabzn;
    int *Tabfft2;

    QString runname;

    int    objectsubunits;
    int    meshsubunits;
    int nmaxpp;
    double meshsize;
    double lambda10n;
    double k0;
    double toleranceobtained;

    int    numberofax1, numberofax2;
    double reflectivity;
    double transmittivity;
    double absorptivity;
    double extinctioncrosssection;
    double absorbingcrosssection;
    double scatteringcrosssection;
    double scatteringcrosssectionwithintegration;
    double scatteringassymetricparam;

    double irra;
    dcmplx E0;

    double opticalforcex;
    double opticalforcey;
    double opticalforcez;
    double opticalforcemodulus;
    double opticaltorquex;
    double opticaltorquey;
    double opticaltorquez;
    double opticaltorquemodulus;

    double  wavelength;
    QString beam;
    QString object;
    QString anisotropy;
    QString material;
    double  discretization;
    double  tolerance;

    double incidenceangle_theta_z;
    double incidenceangle_phi_x;
    double polarizationTM;
    double polarizationTE;
    double polarizationRL;

    double sphereradius;
    double cubeside;
    double cubesidex;
    double cubesidey;
    double cubesidez;
    double positionx;
    double positiony;
    double positionz;
    int nxx;
    int nyy;
    int nzz;
    
    dcmplx epsilon;
    dcmplx epsilon11;
    dcmplx epsilon12;
    dcmplx epsilon13;
    dcmplx epsilon21;
    dcmplx epsilon22;
    dcmplx epsilon23;
    dcmplx epsilon31;
    dcmplx epsilon32;
    dcmplx epsilon33;

    int    localfieldCheck;
    int    macroscopicfieldCheck;
    int    nenergieCheck;
    int    crosssectionCheck;
    int    crosssectionpoyntingCheck;
    int    opticalforceCheck;
    int    opticalforcedensityCheck;
    int    opticaltorqueCheck;
    int    opticaltorquedensityCheck;

};
#endif
