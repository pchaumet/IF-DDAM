#include "cdmRun.h"

Run::Run(QString _runname)
{
  runname = _runname;

  xc = NULL;
  yc = NULL;
  zc = NULL;
  kxy = NULL;
  xy = NULL;
  xcwf = NULL;
  ycwf = NULL;
  zcwf = NULL;
  incidentfield = NULL;
  localfield = NULL;
  macroscopicfield = NULL;
  thetafield = NULL;
  phifield = NULL;
  poyntingfield = NULL;  
  forcex = NULL;
  forcey = NULL;
  forcez = NULL;
  forcexmulti = NULL;
  forceymulti = NULL;
  forcezmulti = NULL;
  torquex = NULL;
  torquey = NULL;
  torquez = NULL;
  torquexmulti = NULL;
  torqueymulti = NULL;
  torquezmulti = NULL;
  incidentfieldx = NULL;
  localfieldx = NULL;
  macroscopicfieldx = NULL;
  incidentfieldy = NULL;
  localfieldy = NULL;
  macroscopicfieldy = NULL;
  incidentfieldz = NULL;
  localfieldz = NULL;
  macroscopicfieldz = NULL;
  polarisafield = NULL;
  epsilonfield = NULL;
  eimagex = NULL;
  eimagey = NULL;
  eimagez = NULL;
  efourierx = NULL;
  efouriery = NULL;
  efourierz = NULL;
  efourierincx = NULL;
  efourierincy = NULL;
  efourierincz = NULL;
  eimageincx = NULL;
  eimageincy = NULL;
  eimageincz = NULL;
  eimagexneg = NULL;
  eimageyneg = NULL;
  eimagezneg = NULL;
  efourierxneg = NULL;
  efourieryneg = NULL;
  efourierzneg = NULL;
  efourierincxneg = NULL;
  efourierincyneg = NULL;
  efourierinczneg = NULL;
  eimageincxneg = NULL;
  eimageincyneg = NULL;
  eimageinczneg = NULL;
  masque = NULL;
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
  FF = NULL;
  FF0 = NULL;
  FFloc = NULL;
  xr = NULL;
  xi = NULL;
  wrk = NULL;
  matindice = NULL;
  matindplan = NULL;
  matind = NULL;
  matrange = NULL;
  a11 = NULL;
  a12 = NULL;
  a13 = NULL;
  a22 = NULL;
  a23 = NULL;
  a31 = NULL;
  a32 = NULL;
  a33 = NULL;
  b11 = NULL;
  b12 = NULL;
  b13 = NULL;
  b22 = NULL;
  b23 = NULL;
  b31 = NULL;
  b32 = NULL;
  b33 = NULL;
  Ediffkzpos = NULL;
  Ediffkzneg = NULL;
  Tabdip  = NULL;
  Tabmulti = NULL;
  Tabzn = NULL;
}
Run::~Run() {
  cleanVectorsMemory();
}
int
Run::checkAvailableMemorySize() {
#ifdef OS
#if OS == LINUX
    long long pages = sysconf(_SC_AVPHYS_PAGES);
    long long page_size = sysconf(_SC_PAGE_SIZE);
    int mem_size = (int)((pages * page_size)/1000000L);
    return mem_size;
#elif OS == WIN32
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return (unsigned long long) (status.ullAvailPhys);
#endif
#endif
   return 0;
}
int
Run::allocateVectorsMemory(int nmax, int ntheta, int nphi, int nfft2d, int obj_num, int n1m, int nplanm, int nmatim, int nbs, int nxm, int nym, int nzm) {
  
  QLOG_INFO() << "Run::allocateVectorsMemory> "; 

  int mem_used = 0;

  QLOG_INFO() << "Run::allocateVectorsMemory::Memory allocation dcmplx:" << (int)sizeof(dcmplx);
  QLOG_INFO() << "Run::allocateVectorsMemory::Memory allocation double:" << (int)sizeof(double);
  xc = (double*) malloc(sizeof(double)*nmax);
  if (xc == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xc,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  yc = (double*) malloc(sizeof(double)*nmax);
  if (yc == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(yc,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  zc = (double*) malloc(sizeof(double)*nmax);
  if (zc == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(zc,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  xcwf = (double*) malloc(sizeof(double)*nmax);
  if (xcwf == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xcwf,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  ycwf = (double*) malloc(sizeof(double)*nmax);
  if (ycwf == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(ycwf,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  zcwf = (double*) malloc(sizeof(double)*nmax);
  if (zcwf == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(zcwf,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  incidentfield = (double*) malloc(sizeof(double)*nmax);
  if (incidentfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(incidentfield,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  localfield = (double*) malloc(sizeof(double)*nmax);
  if (localfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(localfield,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  macroscopicfield = (double*) malloc(sizeof(double)*nmax);
  if (macroscopicfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  } 
  memset(macroscopicfield,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  thetafield = (double*) malloc(sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  if (thetafield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(thetafield,0,sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  mem_used+=max((ntheta+1)*nphi,nfft2d*nfft2d)*sizeof(double) / 1000000L;
  phifield = (double*) malloc(sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  if (phifield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(phifield,0,sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  mem_used+=max((ntheta+1)*nphi,nfft2d*nfft2d)*sizeof(double) / 1000000L;
  poyntingfield = (double*) malloc(sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  if (poyntingfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(poyntingfield,0,sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  mem_used+=max((ntheta+1)*nphi,nfft2d*nfft2d)*sizeof(double) / 1000000L;
  forcex = (double*) malloc(sizeof(double)*nmax);
  if (forcex == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcex,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  forcey = (double*) malloc(sizeof(double)*nmax);
  if (forcey == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcey,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  forcez = (double*) malloc(sizeof(double)*nmax);
  if (forcez == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcez,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  
  forcexmulti = (double*) malloc(sizeof(double)*obj_num);
  if (forcexmulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcexmulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  forceymulti = (double*) malloc(sizeof(double)*obj_num);
  if (forceymulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forceymulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  forcezmulti = (double*) malloc(sizeof(double)*obj_num);
  if (forcezmulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcezmulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  torquex = (double*) malloc(sizeof(double)*nmax);
  if (torquex == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }

  memset(torquex,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  torquey = (double*) malloc(sizeof(double)*nmax);
  if (torquey == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torquey,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  torquez = (double*) malloc(sizeof(double)*nmax);
  if (torquez == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torquez,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;

  torquexmulti = (double*) malloc(sizeof(double)*obj_num);
  if (torquexmulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torquexmulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  torqueymulti = (double*) malloc(sizeof(double)*obj_num);
  if (torqueymulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torqueymulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  torquezmulti = (double*) malloc(sizeof(double)*obj_num);
  if (torquezmulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torquezmulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;



  macroscopicfieldx = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (macroscopicfieldx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(macroscopicfieldx,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  macroscopicfieldy = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (macroscopicfieldy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(macroscopicfieldy,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  macroscopicfieldz = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (macroscopicfieldz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(macroscopicfieldz,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  localfieldx = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (localfieldx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(localfieldx,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  localfieldy = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (localfieldy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(localfieldy,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  localfieldz = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (localfieldz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(localfieldz,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  incidentfieldx = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (incidentfieldx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(incidentfieldx,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  incidentfieldy = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (incidentfieldy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(incidentfieldy,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  incidentfieldz = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (incidentfieldz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(incidentfieldz,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  polarisafield = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3*3);
  if (polarisafield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(polarisafield,0,sizeof(dcmplx)*nmax*3*3);
  mem_used+= sizeof(dcmplx)*nmax*3*3/ 1000000L;
  epsilonfield = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3*3);
  if (epsilonfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(epsilonfield,0,sizeof(dcmplx)*nmax*3*3);
  mem_used+= sizeof(dcmplx)*nmax*3*3/ 1000000L;
  kxy = (double*) malloc(sizeof(double)*nfft2d);
  if (kxy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(kxy,0,sizeof(double)*nfft2d);
  mem_used+= sizeof(double)*nfft2d/ 1000000L;
  xy = (double*) malloc(sizeof(double)*nfft2d);
  if (xy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xy,0,sizeof(double)*nfft2d);
  mem_used+= sizeof(double)*nfft2d/ 1000000L;
  eimagex = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimagex == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimagex,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimagey = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimagey == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimagey,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimagez = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimagez == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimagez,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierx = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierx,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efouriery = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efouriery == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efouriery,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierz = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierz,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierincx = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierincx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierincx,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierincy = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierincy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierincy,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierincz = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierincz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierincz,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageincx = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageincx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageincx,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageincy = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageincy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageincy,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageincz = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageincz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageincz,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  // Negatif
  eimagexneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimagexneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }  
  memset(eimagexneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageyneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageyneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageyneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimagezneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimagezneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimagezneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierxneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierxneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierxneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourieryneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourieryneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourieryneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierzneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierzneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierzneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierincxneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierincxneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierincxneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierincyneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierincyneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierincyneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierinczneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierinczneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierinczneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageincxneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageincxneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageincxneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageincyneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageincyneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageincyneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageinczneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageinczneg == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageinczneg,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  masque = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (masque == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(masque,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
  FF = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (FF == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FF,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  FF0 = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (FF0 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FF0,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  FFloc = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (FFloc == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FFloc,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  xr = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (xr == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xr,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  xi = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (xi == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xi,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  wrk = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3*12);
  if (wrk == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(wrk,0,sizeof(dcmplx)*nmax*3*12);
  mem_used+= sizeof(dcmplx)*nmax*3*12/ 1000000L;
//************TABLEAU QUE SURFACE************************************
 
  matindplan = (int*) malloc(sizeof(int)*nzm*nzm);
  if (matindplan == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(matindplan,0,sizeof(int)*nzm*nzm);
  mem_used+= sizeof(int)*nzm*nzm/ 1000000L;
	
  matind = (int*) malloc(sizeof(int)*(2*n1m*n1m+1));
  if (matind == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(matind,0,sizeof(int)*(2*n1m*n1m+1));
  mem_used+= sizeof(int)*(2*n1m*n1m+1)/ 1000000L;	
	
  matindice = (int*) malloc(sizeof(int)*nzm*(nzm+1)*nmatim/2);
  if (matindice == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(matindice,0,sizeof(int)*nzm*(nzm+1)*nmatim/2);
  mem_used+= sizeof(int)*nzm*(nzm+1)*nmatim/2/ 1000000L;

  
  matrange = (dcmplx*) malloc(sizeof(dcmplx)*5*nbs);
  if (matrange == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(matrange,0,sizeof(dcmplx)*5*nbs);
  mem_used+= sizeof(dcmplx)*5*nbs/ 1000000L;
  
  a11 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym*nplanm);
  if (a11 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(a11,0,sizeof(dcmplx)*4*nxm*nym*nplanm);
  mem_used+= sizeof(dcmplx)*4*nxm*nym*nplanm/ 1000000L;
  a12 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym*nplanm);
  if (a12 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(a12,0,sizeof(dcmplx)*4*nxm*nym*nplanm);
  mem_used+= sizeof(dcmplx)*4*nxm*nym*nplanm/ 1000000L;
  a13 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym*nplanm);
  if (a13 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(a13,0,sizeof(dcmplx)*4*nxm*nym*nplanm);
  mem_used+= sizeof(dcmplx)*4*nxm*nym*nplanm/ 1000000L;
  a22 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym*nplanm);
  if (a22 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(a22,0,sizeof(dcmplx)*4*nxm*nym*nplanm);
  mem_used+= sizeof(dcmplx)*4*nxm*nym*nplanm/ 1000000L;
  a23 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym*nplanm);
  if (a23 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(a23,0,sizeof(dcmplx)*4*nxm*nym*nplanm);
  mem_used+= sizeof(dcmplx)*4*nxm*nym*nplanm/ 1000000L;
  a31 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym*nplanm);
  if (a31 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(a31,0,sizeof(dcmplx)*4*nxm*nym*nplanm);
  mem_used+= sizeof(dcmplx)*4*nxm*nym*nplanm/ 1000000L;
  a32 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym*nplanm);
  if (a32 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(a32,0,sizeof(dcmplx)*4*nxm*nym*nplanm);
  mem_used+= sizeof(dcmplx)*4*nxm*nym*nplanm/ 1000000L;
  a33 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym*nplanm);
  if (a33 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(a33,0,sizeof(dcmplx)*4*nxm*nym*nplanm);
  mem_used+= sizeof(dcmplx)*4*nxm*nym*nplanm/ 1000000L;


  b11 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym);
  if (b11 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(b11,0,sizeof(dcmplx)*4*nxm*nym);
  mem_used+= sizeof(dcmplx)*4*nxm*nym/ 1000000L;
  b12 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym);
  if (b12 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(b12,0,sizeof(dcmplx)*4*nxm*nym);
  mem_used+= sizeof(dcmplx)*4*nxm*nym/ 1000000L;
  b13 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym);
  if (b13 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(b13,0,sizeof(dcmplx)*4*nxm*nym);
  mem_used+= sizeof(dcmplx)*4*nxm*nym/ 1000000L;
  b22 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym);
  if (b22 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(b22,0,sizeof(dcmplx)*4*nxm*nym);
  mem_used+= sizeof(dcmplx)*4*nxm*nym/ 1000000L;
  b23 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym);
  if (b23 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(b23,0,sizeof(dcmplx)*4*nxm*nym);
  mem_used+= sizeof(dcmplx)*4*nxm*nym/ 1000000L;
  b31 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym);
  if (b31 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(b31,0,sizeof(dcmplx)*4*nxm*nym);
  mem_used+= sizeof(dcmplx)*4*nxm*nym/ 1000000L;
  b32 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym);
  if (b32 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(b32,0,sizeof(dcmplx)*4*nxm*nym);
  mem_used+= sizeof(dcmplx)*4*nxm*nym/ 1000000L;
  b33 = (dcmplx*) malloc(sizeof(dcmplx)*4*nxm*nym);
  if (b33 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(b33,0,sizeof(dcmplx)*4*nxm*nym);
  mem_used+= sizeof(dcmplx)*4*nxm*nym/ 1000000L;
  
//************TABLEAU QUE SURFACE************************************
 
  Ediffkzpos = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d*3);
  if (Ediffkzpos == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(Ediffkzpos,0,sizeof(dcmplx)*nfft2d*nfft2d*3);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d*3/ 1000000L;
  
  Ediffkzneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d*3);
  if (Ediffkzneg == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  
  memset(Ediffkzneg,0,sizeof(dcmplx)*nfft2d*nfft2d*3);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d*3/ 1000000L;

  Tabdip = (int*) malloc(sizeof(int)*nmax);
  if (Tabdip == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(Tabdip,0,sizeof(int)*nmax);
  mem_used+= sizeof(int)*nmax/ 1000000L;

  Tabmulti = (int*) malloc(sizeof(int)*nmax);
  if (Tabmulti == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(Tabmulti,0,sizeof(int)*nmax);
  mem_used+= sizeof(int)*nmax/ 1000000L;

  Tabzn = (int*) malloc(sizeof(int)*nmax);
  if (Tabzn == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(Tabzn,0,sizeof(int)*nmax);
  mem_used+= sizeof(int)*nmax/ 1000000L;
  
  return mem_used;
}
void 
Run::cleanVectorsMemory() {
  if (xc) {free(xc); xc = NULL;}
  if (yc) {free(yc); yc = NULL;}
  if (zc) {free(zc); zc = NULL;}
  if (xcwf) {free(xcwf); xcwf = NULL;}
  if (ycwf) {free(ycwf); ycwf = NULL;}
  if (zcwf) {free(zcwf); zcwf = NULL;}
  if (xy) {free(xy); xy = NULL;}
  if (kxy) {free(kxy); kxy = NULL;}
  if (incidentfield) {free(incidentfield); incidentfield = NULL;}
  if (localfield) {free(localfield); localfield = NULL;}
  if (macroscopicfield) {free(macroscopicfield); macroscopicfield = NULL;}
  if (thetafield) {free(thetafield); thetafield = NULL;}
  if (phifield) {free(phifield); phifield = NULL;}
  if (poyntingfield) {free(poyntingfield); poyntingfield = NULL;}
  if (forcex) {free(forcex); forcex = NULL;}
  if (forcey) {free(forcey); forcey = NULL;}
  if (forcez) {free(forcez); forcez = NULL;}
  if (forcexmulti) {free(forcexmulti); forcexmulti = NULL;}
  if (forceymulti) {free(forceymulti); forceymulti = NULL;}
  if (forcezmulti) {free(forcezmulti); forcezmulti = NULL;}
  if (torquex) {free(torquex); torquex = NULL;}
  if (torquey) {free(torquey); torquey = NULL;}
  if (torquez) {free(torquez); torquez = NULL;}
  if (torquexmulti) {free(torquexmulti); torquexmulti = NULL;}
  if (torqueymulti) {free(torqueymulti); torqueymulti = NULL;}
  if (torquezmulti) {free(torquezmulti); torquezmulti = NULL;}
  if (incidentfieldx) {free(incidentfieldx); incidentfieldx = NULL;}
  if (localfieldx) {free(localfieldx); localfieldx = NULL;}
  if (macroscopicfieldx) {free(macroscopicfieldx); macroscopicfieldx = NULL;}
  if (incidentfieldy) {free(incidentfieldy); incidentfieldy = NULL;}
  if (localfieldy) {free(localfieldy); localfieldy = NULL;}
  if (macroscopicfieldy) {free(macroscopicfieldy); macroscopicfieldy = NULL;}
  if (incidentfieldz) {free(incidentfieldz); incidentfieldz = NULL;}
  if (localfieldz) {free(localfieldz); localfieldz = NULL;}
  if (macroscopicfieldz) {free(macroscopicfieldz); macroscopicfieldz = NULL;}
  if (polarisafield) {free(polarisafield); polarisafield = NULL;}
  if (epsilonfield) {free(epsilonfield); epsilonfield = NULL;}
  if (eimagex) {free(eimagex); eimagex = NULL;}
  if (eimagey) {free(eimagey); eimagey = NULL;}
  if (eimagez) {free(eimagez); eimagez = NULL;}
  if (efourierx) {free(efourierx); efourierx = NULL;}
  if (efouriery) {free(efouriery); efouriery = NULL;}
  if (efourierz) {free(efourierz); efourierz = NULL;}
  if (efourierincx) {free(efourierincx); efourierincx = NULL;}
  if (efourierincy) {free(efourierincy); efourierincy = NULL;}
  if (efourierincz) {free(efourierincz); efourierincz = NULL;}
  if (eimageincx) {free(eimageincx); eimageincx = NULL;}
  if (eimageincy) {free(eimageincy); eimageincy = NULL;}
  if (eimageincz) {free(eimageincz); eimageincz = NULL;}
  if (eimagexneg) {free(eimagexneg); eimagexneg = NULL;}
  if (eimageyneg) {free(eimageyneg); eimageyneg = NULL;}
  if (eimagezneg) {free(eimagezneg); eimagezneg = NULL;}
  if (efourierxneg) {free(efourierxneg); efourierxneg = NULL;}
  if (efourieryneg) {free(efourieryneg); efourieryneg = NULL;}
  if (efourierzneg) {free(efourierzneg); efourierzneg = NULL;}
  if (efourierincxneg) {free(efourierincxneg); efourierincxneg = NULL;}
  if (efourierincyneg) {free(efourierincyneg); efourierincyneg = NULL;}
  if (efourierinczneg) {free(efourierinczneg); efourierinczneg = NULL;}
  if (eimageincxneg) {free(eimageincxneg); eimageincxneg = NULL;}
  if (eimageincyneg) {free(eimageincyneg); eimageincyneg = NULL;}
  if (eimageinczneg) {free(eimageinczneg); eimageinczneg = NULL;}
  if (masque) {free(masque); masque = NULL;}
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
  if (FF) {free(FF); FF = NULL;}
  if (FF0) {free(FF0); FF0 = NULL;}
  if (FFloc) {free(FFloc); FFloc = NULL;}
  if (xr) {free(xr); xr = NULL;}
  if (xi) {free(xi);  xi = NULL;}
  if (wrk) {free(wrk); wrk = NULL;}

  
  if (a11) {free(a11); a11 = NULL;}
  if (a12) {free(a12); a12 = NULL;}
  if (a13) {free(a13); a13 = NULL;}
  if (a22) {free(a22); a22 = NULL;}
  if (a23) {free(a23); a23 = NULL;}
  if (a31) {free(a31); a31 = NULL;}
  if (a32) {free(a32); a32 = NULL;}
  if (a33) {free(a33); a33 = NULL;}
  if (b11) {free(b11); b11 = NULL;}
  if (b12) {free(b12); b12 = NULL;}
  if (b13) {free(b13); b13 = NULL;}
  if (b22) {free(b22); b22 = NULL;}
  if (b23) {free(b23); b23 = NULL;}
  if (b31) {free(b31); b31 = NULL;}
  if (b32) {free(b32); b32 = NULL;}
  if (b33) {free(b33); b33 = NULL;}
  
  if (matrange) {free(matrange); matrange = NULL;}
  if (matind) {free(matind); matind = NULL;}
  if (matindice) {free(matindice); matindice = NULL;}
  if (matindplan) {free(matindplan); matindplan = NULL;}

  if (Ediffkzpos) {free(Ediffkzpos); Ediffkzpos = NULL;}
  if (Ediffkzneg) {free(Ediffkzneg); Ediffkzneg = NULL;}
  if (Tabdip) {free(Tabdip); Tabdip = NULL;}
  if (Tabmulti) {free(Tabmulti); Tabmulti = NULL;}
  if (Tabzn) {free(Tabzn); Tabzn = NULL;}
}
double* 
Run::getIncidentField(){
  return incidentfield;
}
double* 
Run::getLocalField(){
  return localfield;
}
double* 
Run::getMacroscopicField(){
  return macroscopicfield;
}
double* 
Run::getXc(){
  return xc;
}
double* 
Run::getYc(){
  return yc;
}
double* 
Run::getZc(){
  return zc;
}
double* 
Run::getXcWF(){
  return xcwf;
}
double* 
Run::getYcWF(){
  return ycwf;
}
double* 
Run::getZcWF(){
  return zcwf;
}
double* 
Run::getThetaField(){
  return thetafield;
}
double* 
Run::getPhiField(){
  return phifield;
}
double* 
Run::getPoyntingField(){
  return poyntingfield;
}
double* 
Run::getForceX(){
  return forcex;
}
double* 
Run::getForceY(){
  return forcey;
}
double* 
Run::getForceZ(){
  return forcez;
}
double* 
Run::getForceXMulti(){
  return forcexmulti;
}
double* 
Run::getForceYMulti(){
  return forceymulti;
}
double* 
Run::getForceZMulti(){
  return forcezmulti;
}
double* 
Run::getTorqueX(){
  return torquex;
}
double* 
Run::getTorqueY(){
  return torquey;
}
double* 
Run::getTorqueZ(){
  return torquez;
}
double* 
Run::getTorqueXMulti(){
  return torquexmulti;
}
double* 
Run::getTorqueYMulti(){
  return torqueymulti;
}
double* 
Run::getTorqueZMulti(){
  return torquezmulti;
}
dcmplx* 
Run::getIncidentFieldX(){
  return incidentfieldx;
}
dcmplx* 
Run::getLocalFieldX(){
  return localfieldx;
}
dcmplx* 
Run::getMacroscopicFieldX(){
  return macroscopicfieldx;
}
dcmplx* 
Run::getIncidentFieldY(){
  return incidentfieldy;
}
dcmplx* 
Run::getLocalFieldY(){
  return localfieldy;
}
dcmplx* 
Run::getMacroscopicFieldY(){
  return macroscopicfieldy;
}
dcmplx* 
Run::getIncidentFieldZ(){
  return incidentfieldz;
}
dcmplx* 
Run::getLocalFieldZ(){
  return localfieldz;
}
dcmplx* 
Run::getMacroscopicFieldZ(){
  return macroscopicfieldz;
}
dcmplx*
Run::getPolarisaField() {
  return (dcmplx*)polarisafield;
}
dcmplx*
Run::getEpsilonField() {
  return (dcmplx*)epsilonfield;
}
double*
Run::getXY() {
  return xy;
}
double*
Run::getKXY() {
  return kxy;
}
dcmplx*
Run::getEimageX() {
  return (dcmplx*)eimagex;
}
dcmplx*
Run::getEimageY() {
  return (dcmplx*)eimagey;
}
dcmplx*
Run::getEimageZ() {
  return (dcmplx*)eimagez;
}
dcmplx*
Run::getEfourierX() {
  return (dcmplx*)efourierx;
}
dcmplx*
Run::getEfourierY() {
  return (dcmplx*)efouriery;
}
dcmplx*
Run::getEfourierZ() {
  return (dcmplx*)efourierz;
}
dcmplx*
Run::getEfourierincX() {
  return (dcmplx*)efourierincx;
}
dcmplx*
Run::getEfourierincY() {
  return (dcmplx*)efourierincy;
}
dcmplx*
Run::getEfourierincZ() {
  return (dcmplx*)efourierincz;
}
dcmplx*
Run::getEimageincX() {
  return (dcmplx*)eimageincx;
}
dcmplx*
Run::getEimageincY() {
  return (dcmplx*)eimageincy;
}
dcmplx*
Run::getEimageincZ() {
  return (dcmplx*)eimageincz;
}


dcmplx*
Run::getEimageXneg() {
  return (dcmplx*)eimagexneg;
}
dcmplx*
Run::getEimageYneg() {
  return (dcmplx*)eimageyneg;
}
dcmplx*
Run::getEimageZneg() {
  return (dcmplx*)eimagezneg;
}
dcmplx*
Run::getEfourierXneg() {
  return (dcmplx*)efourierxneg;
}
dcmplx*
Run::getEfourierYneg() {
  return (dcmplx*)efourieryneg;
}
dcmplx*
Run::getEfourierZneg() {
  return (dcmplx*)efourierzneg;
}
dcmplx*
Run::getEfourierincXneg() {
  return (dcmplx*)efourierincxneg;
}
dcmplx*
Run::getEfourierincYneg() {
  return (dcmplx*)efourierincyneg;
}
dcmplx*
Run::getEfourierincZneg() {
  return (dcmplx*)efourierinczneg;
}
dcmplx*
Run::getEimageincXneg() {
  return (dcmplx*)eimageincxneg;
}
dcmplx*
Run::getEimageincYneg() {
  return (dcmplx*)eimageincyneg;
}
dcmplx*
Run::getEimageincZneg() {
  return (dcmplx*)eimageinczneg;
}
dcmplx*
Run::getMasque() {
  return (dcmplx*)masque;
}
dcmplx*
Run::getFF() {
  return (dcmplx*)FF;
}
dcmplx*
Run::getFF0() {
  return (dcmplx*)FF0;
}
dcmplx*
Run::getFFloc() {
  return (dcmplx*)FFloc;
}
dcmplx*
Run::getxr() {
  return (dcmplx*)xr;
}
dcmplx*
Run::getxi() {
  return (dcmplx*)xi;
}
dcmplx*
Run::getwrk() {
  return (dcmplx*)wrk;
}
int*
Run::getmatindplan() {
  return (int*)matindplan;
}
int*
Run::getmatind() {
  return (int*)matind;
}
int*
Run::getmatindice() {
  return (int*)matindice;
}
dcmplx*
Run::getmatrange() {
  return (dcmplx*)matrange;
}
dcmplx*
Run::geta11() {
  return (dcmplx*)a11;
}
dcmplx*
Run::geta12() {
  return (dcmplx*)a12;
}
dcmplx*
Run::geta13() {
  return (dcmplx*)a13;
}
dcmplx*
Run::geta22() {
  return (dcmplx*)a22;
}
dcmplx*
Run::geta23() {
  return (dcmplx*)a23;
}
dcmplx*
Run::geta31() {
  return (dcmplx*)a31;
}
dcmplx*
Run::geta32() {
  return (dcmplx*)a32;
}
dcmplx*
Run::geta33() {
  return (dcmplx*)a33;
}

dcmplx*
Run::getb11() {
  return (dcmplx*)b11;
}
dcmplx*
Run::getb12() {
  return (dcmplx*)b12;
}
dcmplx*
Run::getb13() {
  return (dcmplx*)b13;
}
dcmplx*
Run::getb22() {
  return (dcmplx*)b22;
}
dcmplx*
Run::getb23() {
  return (dcmplx*)b23;
}
dcmplx*
Run::getb31() {
  return (dcmplx*)b31;
}
dcmplx*
Run::getb32() {
  return (dcmplx*)b32;
}
dcmplx*
Run::getb33() {
  return (dcmplx*)b33;
}


dcmplx*
Run::getEdiffkzpos() {
  return (dcmplx*)Ediffkzpos;
}
dcmplx*
Run::getEdiffkzneg() {
  return (dcmplx*)Ediffkzneg;
}
int*
Run::getTabdip() {
  return (int*)Tabdip;
}
int*
Run::getTabmulti() {
  return (int*)Tabmulti;
}

int*
Run::getTabzn() {
  return (int*)Tabzn;
}

void 
Run::setName(QString _runname){
  runname = _runname;
}
void 
Run::setObjectSubunits(int _objectsubunits){
  objectsubunits = _objectsubunits;
}
void 
Run::setMeshSubunits(int _meshsubunits){
  meshsubunits = _meshsubunits;
}
void 
Run::setNmaxpp(int _nmaxpp){
  nmaxpp = _nmaxpp;
}
void 
Run::setMeshSize(double _meshsize){
  meshsize = _meshsize;
}
void 
Run::setLambda10n(double _lambda10n){
  lambda10n = _lambda10n;
}
void 
Run::setK0(double _k0){
  k0 = _k0;
}
void 
Run::setToleranceObtained(double _toleranceobtained){
  toleranceobtained = _toleranceobtained;
}
void 
Run::setNumberofAx1(int _numberofax1){
  numberofax1 = _numberofax1;
}
void 
Run::setNumberofAx2(int _numberofax2){
  numberofax2 = _numberofax2;
}
void 
Run::setReflectivity(double _reflectivity){
  reflectivity = _reflectivity;
}
void 
Run::setTransmittivity(double _transmittivity){
  transmittivity = _transmittivity;
}
void 
Run::setAbsorptivity(double _absorptivity){
  absorptivity = _absorptivity;
}
void 
Run::setExtinctionCrossection(double _extinctioncrosssection){
  extinctioncrosssection = _extinctioncrosssection;
}
void 
Run::setAbsorbingCrossection(double _absorbingcrosssection){
  absorbingcrosssection = _absorbingcrosssection;
}
void 
Run::setScatteringCrossection(double _scatteringcrosssection){
  scatteringcrosssection = _scatteringcrosssection;
}
void 
Run::setScatteringCrossectionWithIntegration(double _scatteringcrosssectionwithintegration){
  scatteringcrosssectionwithintegration = _scatteringcrosssectionwithintegration;
}
void 
Run::setScatteringAssymetricParam(double _scatteringassymetricparam){
  scatteringassymetricparam = _scatteringassymetricparam;
}
void 
Run::setIrra(double _irra){
  irra = _irra;
}
void 
Run::setE0(dcmplx _E0){
  E0 = _E0;
}
void 
Run::setOpticalForcex(double _opticalforcex){
  opticalforcex = _opticalforcex;
}
void 
Run::setOpticalForcey(double _opticalforcey){
  opticalforcey = _opticalforcey;
}
void 
Run::setOpticalForcez(double _opticalforcez){
  opticalforcez = _opticalforcez;
}
void 
Run::setOpticalForceModulus(double _opticalforcemodulus){
  opticalforcemodulus = _opticalforcemodulus;
}
void 
Run::setOpticalTorquex(double _opticaltorquex){
  opticaltorquex = _opticaltorquex;
}
void 
Run::setOpticalTorquey(double _opticaltorquey){
  opticaltorquey = _opticaltorquey;
}
void 
Run::setOpticalTorquez(double _opticaltorquez){
  opticaltorquez = _opticaltorquez;
}
void 
Run::setOpticalTorqueModulus(double _opticaltorquemodulus){
  opticaltorquemodulus = _opticaltorquemodulus;
}

QString
Run::getName(){
  return runname;
}
int
Run::getObjectSubunits(){
  return objectsubunits;
}
int
Run::getMeshSubunits(){
  return meshsubunits;
}
int
Run::getNmaxpp(){
  return nmaxpp;
}
double
Run::getMeshSize(){
  return meshsize;
}
double
Run::getLambda10n(){
  return lambda10n;
}
double
Run::getK0(){
  return k0;
}
double
Run::getToleranceObtained(){
  return toleranceobtained;
}
int
Run::getNumberofAx1(){
  return numberofax1;
}
int
Run::getNumberofAx2(){
  return numberofax2;
}
double
Run::getReflectivity(){
  return reflectivity;
}
double
Run::getTransmittivity(){
  return transmittivity;
}
double
Run::getAbsorptivity(){
  return absorptivity;
}
double
Run::getExtinctionCrossection(){
  return extinctioncrosssection;
}
double
Run::getAbsorbingCrossection(){
  return absorbingcrosssection;
}
double
Run::getScatteringCrossection(){
  return scatteringcrosssection;
}
double
Run::getScatteringCrossectionWithIntegration(){
  return scatteringcrosssectionwithintegration;
}
double
Run::getScatteringAssymetricParam(){
  return scatteringassymetricparam;
}
double
Run::getIrra(){
  return irra;
}
dcmplx
Run::getE0(){
  return E0;
}
double
Run::getOpticalForcex(){
  return opticalforcex;
}
double
Run::getOpticalForcey(){
  return opticalforcey;
}
double
Run::getOpticalForcez(){
  return opticalforcez;
}
double
Run::getOpticalForceModulus(){
  return opticalforcemodulus;
}
double
Run::getOpticalTorquex(){
  return opticaltorquex;
}
double
Run::getOpticalTorquey(){
  return opticaltorquey;
}
double
Run::getOpticalTorquez(){
  return opticaltorquez;
}
double
Run::getOpticalTorqueModulus(){
  return opticaltorquemodulus;
}


