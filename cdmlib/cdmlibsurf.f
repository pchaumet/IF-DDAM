c     Fortran library subroutine entry point
      subroutine cdmlibsurf(
c     input file cdm.in
     $     lambda,beam,object,trope,
     $     materiaumulti,nnnr,tolinit,methodeit,polarizability,
     $     nlecture,filereread,nmatf,h5file,
C     Definition du multicouche
     $     neps,zcouche,epscouche,materiaucouche,
c     output file cdm.out
     $     nlocal,nmacro,ncote,nsection,ndiffracte,nquickdiffracte,nrig,
     $     ninterp,nforce,nforced,ntorque,ntorqued,nproche,
     $     nlentille,nquicklens,nenergie,nobjet,
c     cube, sphere (includes multiple)
     $     density,side, sidex, sidey, sidez, hauteur,
     $     numberobjet, rayonmulti, xgmulti, ygmulti, zgmulti,
     $     epsmulti, epsanimulti,lc,hc,ng,
c     ellipsoid+arbitrary
     $     demiaxea,demiaxeb,demiaxec,thetaobj,phiobj,psiobj,
     $     namefileobj,
c     planewavecircular.in / planewavelinear.in files
     $     theta, phi,pp,ss,ir,P0,w0,xgaus,ygaus,zgaus,namefileinc,
     $     thetam,phim,ppm,ssm,E0m,nbinc,
c     return info stringf
     $     infostr, nstop,
c     return scalar results
     $     nbsphere, ndipole, aretecube,
     $     lambda10n, k0, tol1, ncompte, nloop,
     $     efficacite,efficaciteref,efficacitetrans,
     $     Cext,Cabs,Csca,Cscai,gasym,irra, E0,
     $     forcet, forcem,
     $     couplet, couplem,
     $     nxm, nym, nzm, nxmp, nymp, nzmp, nmaxpp,
     $     incidentfield, localfield, macroscopicfield,
     $     xs, ys, zs, xswf, yswf, zswf,
     $     ntheta, nphi, thetafield,phifield,poyntingfield,
     $     forcex,forcey,forcez,forcexmulti,forceymulti,forcezmulti,
     $     torquex,torquey,torquez,torquexmulti,torqueymulti,
     $     torquezmulti,
     $     incidentfieldx, incidentfieldy,incidentfieldz,
     $     localfieldx, localfieldy, localfieldz,
     $     macroscopicfieldx, macroscopicfieldy, macroscopicfieldz,
     $     polarisa,epsilon,
     $     nfft2d,Eimagexpos,Eimageypos,Eimagezpos,
     $     Eimageincxpos,Eimageincypos,Eimageinczpos,
     $     Efourierxpos, Efourierypos,Efourierzpos,
     $     Efourierincxpos,Efourierincypos, Efourierinczpos,
     $     Eimagexneg,Eimageyneg,Eimagezneg,
     $     Eimageincxneg,Eimageincyneg,Eimageinczneg,
     $     Efourierxneg,Efourieryneg,Efourierzneg,
     $     Efourierincxneg,Efourierincyneg,Efourierinczneg,masque,
     $     kxy,xy,numaperref,numapertra,numaperinc
     $     ,gross,zlensr,zlenst,ntypemic,
c     passe certains arguments pour le dimensionnement
     $     n1m,nmatim,nplanm,nbs,
c     fonction de green de la surface
c     matrange(nbs,5),matindice(nplanm,nmatim),matind(0:2*n1m*n1m)
c     ,matindplan(nzm,nzm)
     $     matrange,matindice,matindplan,matind,
c     a11,a12,a13,a22,a23,a31,a32,a33(2*nxm,2*nym,nplanm)nplanm=nzm*(nzm+1)/2
     $     a11,a12,a13,a22,a23,a31,a32,a33,
     $     b11,b12,b13,b22,b23,b31,b32,b33,
c     taille double complex (3*nxm*nym*nzm)
     $     FF,FF0,FFloc,xr,xi,
c     taille double complex (3*nxm*nym*nzm,12)
     $     wrk,
c     taille double complex (nfft2d,nfft2d,3)
     $     Ediffkzpos,Ediffkzneg,      
c     taille entier (nxm*nym*nzm)
     $     Tabdip,Tabmulti,Tabzn)

      use HDF5

      implicit none
      
c     integer
      integer ii,jj,kk,ll,i,j,k,l,i2,j2,ii2,jj2,cnt,nstop,cntwf,kkm,jjm
     $     ,iim,nrig,nhomo,nmatf
      integer  nlocal,nmacro,nsection,nforce ,nforced
     $     ,nquickdiffracte,nquicklens,ndiffracte,ncote,npolainc
     $     ,ntorque,ntorqued ,nsens ,nproche,nprochefft ,nlecture
     $     ,nlecture1,long,long1 ,ierror,ir,nenergie,ntypemic

c     variables for the object
      integer nbsphere3,nbsphere,ndipole,IP(3),test,numberobjetmax
     $     ,numberobjet,is,ng
      parameter (numberobjetmax=20)
      integer nx,ny,nz,nx2,ny2,nxy2,nz2,nxm,nym,nzm,nxmp,nymp,nzmp
     $     ,ntotal,nxm2,nym2,nzm2,nxym2,nxmpp,nympp,nzmpp,nmaxpp
      integer subunit,nsubunit,comparaison
c     definition of the size for the code
      INTEGER nmax, ntotalm

c     variables for the positions
      double precision x,y,z,xmin,xmax,ymin,ymax ,zmin,zmax,rayon
     $     ,density,side,sidex,sidey,sidez,hauteur
     $     ,xgmulti(numberobjetmax),ygmulti(numberobjetmax)
     $     ,zgmulti(numberobjetmax),rayonmulti(numberobjetmax),demiaxea
     $     ,demiaxeb,demiaxec,thetaobj,phiobj,psiobj,t0,t1,t2,ti,tf,lc
     $     ,hc ,sidemic
      double precision aretecube
      integer iphi,itheta,nphi,ntheta,nthetamod
      DOUBLE PRECISION,DIMENSION(nxm*nym*nzm)::xs,ys,zs
      DOUBLE PRECISION,DIMENSION(nxm*nym*nzm)::xswf,yswf,zswf
      DOUBLE PRECISION,DIMENSION((ntheta+1)*nphi)::thetafield,phifield
     $     ,poyntingfield
      double precision pi,lambda,lambda10n,k0,k03,epi,epr,c

c     variables for the material
      double precision quatpieps0
      double complex eps0
      double complex, dimension(nxm*nym*nzm,3,3) :: polarisa,epsilon
      double complex eps,epsani(3,3),epsmulti(numberobjetmax)
     $     ,epsanimulti(numberobjetmax,3,3)
      character(2) polarizability
      character (64), DIMENSION(numberobjetmax) :: materiaumulti
      character(64) materiau,object,beam,namefileobj,namefileinc
     $     ,filereread,filereread1,message
      character(3) trope,file1
c     variables for the incident field and local field
      DOUBLE PRECISION, DIMENSION(nxm*nym*nzm) :: incidentfield,
     $     localfield,macroscopicfield,forcex,forcey,forcez, torquex
     $     ,torquey,torquez
      double precision forcexmulti(numberobjetmax)
     $     ,forceymulti(numberobjetmax),forcezmulti(numberobjetmax)
     $     ,torquexmulti(numberobjetmax),torqueymulti(numberobjetmax)
     $     ,torquezmulti(numberobjetmax)
      integer nbinc
      double precision ss,pp,theta,phi,psi,I0,Emod,tmp,thetad,phim(10)
     $     ,thetam(10),ssm(10),ppm(10),Emod11,Emod22,Emod12 ,Emod21
      double complex Eloc(3),Em(3),E0,uncomp,icomp,zzero,E0m(10),Emx,Emy
     $     ,Emz
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldx
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldy
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldz
      double complex, dimension(nxm*nym*nzm) :: localfieldx
      double complex, dimension(nxm*nym*nzm) :: localfieldy
      double complex, dimension(nxm*nym*nzm) :: localfieldz
      double complex, dimension(nxm*nym*nzm) :: incidentfieldx
      double complex, dimension(nxm*nym*nzm) :: incidentfieldy
      double complex, dimension(nxm*nym*nzm) :: incidentfieldz
      double complex Stenseur(3,3)
      double complex, dimension(3*nxm*nym*nzm) :: FF,FF0,FFloc

c     Green function
      integer, dimension(nxm*nym*nzm) :: Tabdip,Tabmulti,Tabzn
      integer indice,indicex,indicey,nmat,nbs,n1m,nmatim,n1max,nplanm
     $     ,ninterp ,ninter

c     parameter(n1m=min(nxm,nym)**2+max(nxm,nym)**2,nmatim =min(nxm,nym)
c     $     *(2*max(nxm,nym)-min(nxm,nym)+1))
c     parameter (nplanm=nzm*(nzm+1)/2)
c     parameter(ntotalm=4*nxm*nym,nbs=nzm*nzm *nmatim)
      double precision a(0:2*n1m*n1m)
      double complex matrange(nbs,5)
      integer matindice(nplanm,nmatim),matind(0:2*n1m*n1m)
     $     ,matindplan(nzm,nzm)

      double complex a11(2*nxm,2*nym,nplanm),a12(2*nxm,2*nym,nplanm),
     $     a13(2*nxm,2*nym,nplanm),a22(2*nxm,2*nym,nplanm),a23(2*nxm,2
     $     *nym,nplanm),a31(2*nxm,2*nym,nplanm),a32(2*nxm,2*nym,nplanm)
     $     , a33(2*nxm,2*nym,nplanm),b11(4*nxm*nym),b12(4*nxm*nym),b13(4
     $     *nxm*nym),b22(4*nxm*nym),b23(4*nxm*nym),b31(4 *nxm*nym),b32(4
     $     *nxm*nym),b33(4*nxm*nym)
c     Ftmp(3*nmax)



c     Variables for the computation of the optical force and torque
      double precision forcet(3),forcem,forcemie
      double precision couplet(3),xg,yg,zg,couplem
      double complex Eder(3,3)
      
c     computation of the cross section
      integer imaxk0,nfft2d,nfft2dtmp
      double precision normal(3) ,deltatheta,deltaphi,Csca,Cscai,Cabs
     $     ,Cext,gasym,thetas,phis,MIECEXT,MIECABS,MIECSCA,GSCA
      double complex ctmp,masque(nfft2d,nfft2d),Ediffkzpos(nfft2d
     $     ,nfft2d,3),Ediffkzneg(nfft2d,nfft2d,3)
      
c     variables for the iterative method
      INTEGER ldabi, nlar
      integer nnnr,ncompte,nt
      integer NLIM,ndim,nou,maxit,nstat,nloop,STEPERR
      DOUBLE PRECISION  NORM,TOL,DZNRM2,norm1,norm2,tolinit,tol1,tole

c     COMMON /ONTHEHEAP/ b,xr,xi,wrk
      double complex, dimension(3*nxm*nym*nzm) :: xr,xi
      double complex, dimension(3*nxm*nym*nzm,12) :: wrk
      
c     double complex wrk(*), xi(*), xr(*)
c     POINTER ( xr_p, xr ), ( b_p, b )
c     POINTER ( wrk_p, wrk ), ( xi_p, xi)

c     Poynting vector
      double precision Poyntinginc

c     variable for the multilayer
      integer neps,nepsmax,numerocouche
      parameter (nepsmax=8)
      double precision dcouche(nepsmax),zcouche(0:nepsmax),hcc,epsabs
     $     ,rloin,indicen,indice0,indicem
      double complex epscouche(0:nepsmax+1)
      character (64), DIMENSION(0:nepsmax+1) :: materiaucouche

c     Info string
      character(64) infostr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     nouvelle variable a passer en argument d'entree
c     power et diametre
      double precision P0,w0,xgaus,ygaus,zgaus
      character(12) methodeit

c     nouvelle variable de sortie Irra plus variable de flux pour le
c     faisceau gaussien
      double precision irra,fluxinc ,fluxref ,fluxtrans,fluxreftot
     $     ,fluxtratot,efficacite,efficaciteref,efficacitetrans
      
c     variable pour avoir l'image a travers la lentille
      integer nlentille,nobjet,nfft2d2,nmasque
      double precision kx,ky,kz,kp2,deltakx,deltaky,numaperinc ,deltax
     $     ,signe,zlenst,zlensr,numaperref,numapertra
      double precision kxy(nfft2d),xy(nfft2d),gross
      double complex Eimagexpos(nfft2d*nfft2d),Eimageypos(nfft2d*nfft2d)
     $     ,Eimagezpos(nfft2d*nfft2d),Eimageincxpos(nfft2d*nfft2d)
     $     ,Eimageincypos(nfft2d*nfft2d),Eimageinczpos(nfft2d*nfft2d)
     $     ,Efourierxpos(nfft2d*nfft2d),Efourierypos(nfft2d*nfft2d)
     $     ,Efourierzpos(nfft2d*nfft2d),Efourierincxpos(nfft2d*nfft2d)
     $     ,Efourierincypos(nfft2d*nfft2d),Efourierinczpos(nfft2d
     $     *nfft2d),Eimagexneg(nfft2d*nfft2d),Eimageyneg(nfft2d*nfft2d)
     $     ,Eimagezneg(nfft2d*nfft2d),Eimageincxneg(nfft2d*nfft2d)
     $     ,Eimageincyneg(nfft2d*nfft2d),Eimageinczneg(nfft2d*nfft2d)
     $     ,Efourierxneg(nfft2d*nfft2d),Efourieryneg(nfft2d*nfft2d)
     $     ,Efourierzneg(nfft2d*nfft2d),Efourierincxneg(nfft2d*nfft2d)
     $     ,Efourierincyneg(nfft2d*nfft2d),Efourierinczneg(nfft2d
     $     *nfft2d)
      integer iret,omp_get_max_threads
c     declaration pour FFT
      integer*8 planf,planb,plan2f,plan2b,planfn,planbn
      integer FFTW_FORWARD,FFTW_ESTIMATE,FFTW_BACKWARD
c     declaration pour HDF5
      character(40) :: name
      character(LEN=100) :: h5file
      character(LEN=100) :: datasetname
      integer debug
      integer error
      
      integer(hid_t) :: file_id
      integer(hid_t) :: group_idopt,group_idmic,group_idnf,group_idof
     $     ,group_idff,group_iddip
      integer :: dim(4)
c     declaration calcul temps
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer values(8),values2(8),valuesi(8),valuesf(8)
      
      call cpu_time(ti)
      call date_and_time(date,time,zone,valuesi)
      message=' for the execution of the code '
      call calculatedate(valuesf,valuesi,tf,ti,message)

      call dfftw_init_threads(iret)
      if (iret.eq.0) then
         write(*,*) 'iret',iret
         infostr='Unlikely error during thread initialization'
         nstop=1
         return
      endif
      CALL dfftw_plan_with_nthreads(omp_get_max_threads())
      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64

      

      
      ncote=ncote-1

c     FF0 : champ incident
c     FF  : dipole
c     FFloc : champ local
c     xi: champ local: identique FFloc. Sert dans la methode iterative.
c     xr=(I-A polarisa)xi. A matrice des fonctions de green et polarisa
c     est la polarisabilite=FF0 quand methode iterative d'inversion
c     convergee.
      if (nquickdiffracte.eq.1) nquicklens=1
      if (nquicklens.eq.1) nquickdiffracte=1
      write(*,*) '*************************************************'
      write(*,*) '******************* INPUT DATA ******************'
      write(*,*) '*************************************************'
      write(*,*) 'Wavelength       :',lambda,'nm'
      write(*,*) 'Beam             : ',trim(beam)
      write(*,*) 'Object           : ',trim(object)
      write(*,*) 'Isotropy         : ',trope
      write(*,*) 'Discretization   : ',nnnr
      write(*,*) 'Iterative method : ',trim(methodeit),'tolerance asked'
     $     ,tolinit
      write(*,*) 'Born or rigorous : ',nrig,':0 rigourous'
      write(*,*) 'Green s function : ',ninterp,':0 rigourous'
      write(*,*) 'Write mat file   : ',nmatf,':0 write mat file'
      write(*,*) 'Local field      : ',nlocal,':1 compute local field'
      write(*,*) 'Macroscopic field: ',nmacro,':1 compute mac. field'
      write(*,*) 'Cross section    : ',nsection,'Csca',ndiffracte
      write(*,*) 'Quick cross sec. : ',nquickdiffracte,'FFT used'
      write(*,*) 'Emissivity       : ',nenergie,':1 compute energy'
      write(*,*) 'Lens             : ',nlentille,':1 compute microscopy'
      write(*,*) 'Quick Lens       : ',nquicklens,'FFT used'
      write(*,*) 'Pos. focal refl. : ',zlensr,'nm'
      write(*,*) 'Pos. focal tran. : ',zlenst,'nm'
      write(*,*) 'NA reflexion     : ',numaperref
      write(*,*) 'NA transmission  : ',numapertra
      write(*,*) 'NA Condenser     : ',numaperinc
      write(*,*) 'Side computation : ',ncote,':0 both side'
c      write(*,*) 'Optical force    : ',nforce ,'Density',nforced
c      write(*,*) 'Optical torque   : ',ntorque,'Density',ntorqued
      write(*,*) 'Near field       : ',nproche
      write(*,*) 'Polarizability   : ',polarizability
      write(*,*) 'Waist            : ',w0,'nm'
      write(*,*) 'Power            : ',P0,'W'
      write(*,*) 'Object           : ',nobjet,':1 compute only dipole'
      write(*,*) 'Box size         : ',nxm,nym,nzm
      write(*,*) 'Meshsize         : ',aretecube,'nm'
      write(*,*) 'Number of layer  : ',neps
      write(*,*) 'Write  file      : ',nmatf
     $     ,':0 ascii file: 1 no file: 2 hdf5 file'
      write(*,*) 'Use FFTW'

      if (nmatf.eq.2) then
         debug=1
         call hdf5create(h5file, file_id)
         write(*,*) 'h5 file created :',h5file
         write(*,*) 'file_id         :', file_id
         call h5gcreate_f(file_id,"Option", group_idopt, error)
         call h5gcreate_f(file_id,"Object", group_iddip, error)
         call h5gcreate_f(file_id,"Far Field", group_idff, error)
         call h5gcreate_f(file_id,"Microscopy", group_idmic, error)
c         call h5gcreate_f(file_id,"Optical Force", group_idof, error)
c         write(*,*) 'error',error
         call h5gcreate_f(file_id,"Near Field", group_idnf, error)
      endif


      
c     Initialise nstop
      nstop=0

c     arret de suite si nnnr trop petit par rapport a n*m
      if (nnnr.gt.min(nxm,nym
     $     ,nzm).and.object(1:9).ne.'arbitrary'.and.Object(1:7).ne
     $     .'cuboid2' .and. Object(1:13).ne.'inhomocuboid2') then
         infostr='nxm nym or nzm smaller than discretization'
         nstop = 1
         return         
      endif
c     arret de suite si pas assez de place pour propa
      if (ninterp.ne.0) then
         i=min(nxm,nym)*(2*max(nxm,nym)-min(nxm,nym)+1)/2
         tmp=dble(i)/(aint(dsqrt(dble(nxm*nym+nym*nym))+1.d0)
     $        *dble(ninterp))
         write(*,*) '************* Green function ***************'
         write(*,*) 'Rigorous: Number of Green function :',i
         write(*,*) 'Interpolation Level                :',ninterp
         write(*,*) 'Number of Green function           :'
     $        ,nint((aint(dsqrt(dble(nxm*nym+nym*nym))+1.d0))*ninterp)
         if (tmp.le.1.d0) then
            infostr='Not enough space: decrease ninterp'
            nstop=-1
            return
         endif         
      endif

      nhomo=0
      ctmp=epscouche(0)
      do i=0,neps+1         
         if (epscouche(i).ne.ctmp) nhomo=1
      enddo

      if (nsection.eq.1.and.nhomo.ne.0) then 
         infostr='Can not compute cross section: media no homogeneous'
         nstop=-1
         return
      endif
      
      if (nsection.eq.1.and.beam(1:11).ne.'pwavelinear'.and.beam(1:13)
     $     .ne.'pwavecircular') then 
         infostr='Can not compute cross section: Not a plane wave'
         nstop=-1
         return
      endif

      if (nenergie.eq.1.and.ncote.ne.0) then
         ncote=0 
         write(*,*) 'Due to the computation of the energy the'
         write(*,*) 'code compute the diffracted field for both side'
      endif
      
      if (nlentille.eq.1) nenergie=1

c      if (nenergie.eq.1.and.nquickdiffracte.eq.1) then
c         ndiffracte=1
c      endif
    



      materiau = materiaumulti(1)
c     ne fait rien
      if (nobjet.eq.0.and.nlocal.eq.0.and.nmacro.eq.0
     $     .and.nforce.eq.0.and.nforced.eq.0
     $     .and.ntorque.eq.0.and.ntorqued.eq.0
     $     .and.nlentille.eq.0.and.ndiffracte.eq.0.and.nquickdiffracte
     $     .eq. 0.and.nsection.eq.0.and.nenergie.eq.0) then 
         infostr='No calculation requested!'
         nstop = -1;
         return
      endif

      
c     calculation size parameter initialization
      ntotalm=4*nxm*nym
      nmax = nxm*nym*nzm
      nfft2d2=nfft2d/2
      ldabi = 3*nxm*nym*nzm
      nlar = 12
      nstop=0
      nfft2dtmp=nfft2d

      xg = xgmulti(1)
      yg = ygmulti(1)
      zg = zgmulti(1)
      rayon = rayonmulti(1)
      eps = epsmulti(1)
      epsani = epsanimulti(1,:,:)
c      write(*,*) 'Relative permittivity of the object',epsmulti(1)
c     pas assez discrétisé
      if (nmax.lt.8) then
         nstop=1
         infostr='Check Nxm Nym and Nzm!'
         return
      endif
      if (nnnr.lt.2.and.object(1:9).ne.'arbitrary') then
         nstop=1
         write(*,*) 'There is no discretization!'
         infostr='There is no discretization!'
         return
      endif
c     arret pour tolerance dans la méthode itérative
      if (tolinit.lt.1.d-12) then 
         infostr='Tolerance too small!'
         nstop = 1;
         return
      endif
      if (tolinit.gt.0.1d0) then 
         infostr='Tolerance too large!'
         nstop = 1;
         return
      endif
      write(*,*) '************* Object ***************'
      do k=1, numberobjet
         write (*,*) 'Material object',k,':',materiaumulti(k)
      enddo
      if (nmax.lt.8) then
         nstop=1
         infostr='Check Nxm Nym and Nzm'
         return
      endif
c     if (nx.lt.2.or.ny.lt.2) then
c     nstop=1
c     write(*,*) 'There is no discretization!',nx,ny,nz
c     infostr='There is no discretization!'
c     return
c     endif
  
c     open the output file:
      open(99,file='output')
      write(99,*) '************* OUTPUT FILE ***************'

c     Intensity of the incident field
      open(36,file='incidentfieldx.mat')
      open(37,file='incidentfieldy.mat')
      open(38,file='incidentfieldz.mat')
      open(39,file='incidentfield.mat')
c     Intensity of the local field
      open(40,file='localfieldx.mat')
      open(41,file='localfieldy.mat')
      open(42,file='localfieldz.mat')
      open(43,file='localfield.mat')
c     Intensity of the macroscopic field
      open(44,file='macroscopicfieldx.mat')
      open(45,file='macroscopicfieldy.mat')
      open(46,file='macroscopicfieldz.mat')
      open(47,file='macroscopicfield.mat')
c     Intensity of the incident field
      open(136,file='incidentfieldxwf.mat')
      open(137,file='incidentfieldywf.mat')
      open(138,file='incidentfieldzwf.mat')
      open(139,file='incidentfieldwf.mat')
c     Intensity of the local field
      open(140,file='localfieldxwf.mat')
      open(141,file='localfieldywf.mat')
      open(142,file='localfieldzwf.mat')
      open(143,file='localfieldwf.mat')
c     Intensity of the macroscopic field
      open(144,file='macroscopicfieldxwf.mat')
      open(145,file='macroscopicfieldywf.mat')
      open(146,file='macroscopicfieldzwf.mat')
      open(147,file='macroscopicfieldwf.mat')
c     Far field discretization
      open(67,file='xwf.mat')
      open(68,file='ywf.mat')
      open(69,file='zwf.mat')
c     save the Poynting vecteur
      open(50,file='poynting.mat')
      open(51,file='theta.mat')
      open(52,file='phi.mat')

      
c     save the diffracted field in kx,ky and the total field
      open(500,file='Ediffkposx.mat')
      open(501,file='Ediffkposy.mat')
      open(502,file='Ediffkposz.mat')
      
      open(506,file='Ediffknegx.mat')
      open(507,file='Ediffknegy.mat')
      open(508,file='Ediffknegz.mat')
      
      open(53,file='poyntingpos.mat')
      open(54,file='poyntingneg.mat')

      open(600,file='fourierposinc.mat')
      open(601,file='fourierposincx.mat')
      open(602,file='fourierposincy.mat')
      open(603,file='fourierposincz.mat')
      
      open(604,file='fourierneginc.mat')   
      open(605,file='fouriernegincx.mat')
      open(606,file='fouriernegincy.mat')
      open(607,file='fouriernegincz.mat')

      open(608,file='fourierpos.mat')
      open(609,file='fourierposx.mat')
      open(610,file='fourierposy.mat')
      open(611,file='fourierposz.mat')
      
      open(612,file='fourierneg.mat')   
      open(613,file='fouriernegx.mat')
      open(614,file='fouriernegy.mat')
      open(615,file='fouriernegz.mat')

      open(616,file='imageposinc.mat')
      open(617,file='imageposincx.mat')
      open(618,file='imageposincy.mat')
      open(619,file='imageposincz.mat')
      
      open(620,file='imageneginc.mat')   
      open(621,file='imagenegincx.mat')
      open(622,file='imagenegincy.mat')
      open(623,file='imagenegincz.mat')

      open(624,file='imagepos.mat')
      open(625,file='imageposx.mat')
      open(626,file='imageposy.mat')
      open(627,file='imageposz.mat')
      
      open(628,file='imageneg.mat')   
      open(629,file='imagenegx.mat')
      open(630,file='imagenegy.mat')
      open(631,file='imagenegz.mat')

      
      open(650,file='kxfourier.mat')
      open(651,file='kyfourier.mat')

      open(520,file='kx.mat')
      open(521,file='ky.mat')
c     save the density of the optical force
      open(60,file='forcex.mat')
      open(61,file='forcey.mat')
      open(62,file='forcez.mat')
c     save the density of optical torque
      open(63,file='torquex.mat')
      open(64,file='torquey.mat')
      open(65,file='torquez.mat')
c     save epsilon
      open(66,file='epsilon.mat')
c     initialization of the data
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      pi=dacos(-1.d0)
      zzero=(0.d0,0.d0)
      c=299792458.d0
      eps0=(1.d0,0.d0)
      quatpieps0=1.d0/(c*c*1.d-7)

      infostr='STARTED'

c     transform the input data in meter
      lambda=lambda*1.d-9
      w0=w0*1.d-9
      do i=0,neps
         zcouche(i)=zcouche(i)*1.d-9
      enddo
c     check the potision of the layers
      do i=1,neps
         if (zcouche(i).le.zcouche(i-1)) then
            nstop=1
            infostr='problem with the position of the layers'
            write(99,*) 'problem with the position of the layers'
            write(99,*) 'zcouche',zcouche(i-1),i-1
            write(99,*) 'zcouche',zcouche(i),i
         endif
      enddo
      if (w0.le.0.d0) then
         nstop=1
         infostr='waist=0'
         return
      endif
      if (P0.le.0.d0) then
         nstop=1
         infostr='Power=0'
         return
      endif

c     compute the relative permittivity fron a database versus the
c     wavelength of illumination
c     write(*,*) 'Relative permittivity',eps,materiau(1:2),lambda
      if (object(1:9).eq.'arbitrary') goto 111
      do k=1,numberobjet
         materiau = materiaumulti(k)
         if (materiau(1:2).ne.'xx') then
            
            call interpdielec(lambda,materiau,epr,epi,infostr,nstop)
            if (nstop.eq.1) return
            epsmulti(k)=(epr*uncomp+icomp*epi)
            write(99,*) 'Object relative permittivity',eps
            write(99,*) 'Object relative permittivity',eps,materiau(1:2)
     $           ,lambda
         else
            if (trope(1:3).eq.'iso') then
               write(*,*) 'Object relative permittivity',eps
               write(99,*) 'Object relative permittivity',eps
            else 
               do i=1,3
                  do j=1,3
                     write(*,*) 'Object relative permittivity',epsani(i
     $                    ,j),i,j
                     write(99,*) 'Object relative permittivity',epsani(i
     $                    ,j),i,j
                  enddo
               enddo
            endif
         endif
      enddo
      eps = epsmulti(1)
      materiau = materiaumulti(1)
      write(*,*) '*********** Multilayer ****************'
      do k=0,neps+1
         materiau = materiaucouche(k)
         if (materiau(1:2).ne.'xx') then
            
            call interpdielec(lambda,materiau,epr,epi,infostr,nstop)
            if (nstop.eq.1) return
            epscouche(k)=(epr*uncomp+icomp*epi)
            write(99,*) 'Relative permittivity of the layer',k
     $           ,epscouche(k)
            write(99,*) 'Relative permittivity of the layer'
     $           ,epscouche(k),materiau(1:2),lambda
         else
            write(*,*) 'Relative permittivity of the layer',k
     $           ,epscouche(k)
            write(99,*) 'Relative permittivity of the layer',k,eps
     $           couche(k)          
         endif
      enddo

      
c     epsilon of the layers
 111  if (dimag(epscouche(0)).gt.0.d0) then
         infostr='absorbing incident layer'
         nstop=-1
         return
      endif

      if (ndiffracte.eq.1) then
         ncote=0
         if (dimag(epscouche(neps+1)).gt.0.d0) ncote=-1
      endif

      indice0=dsqrt(dreal(epscouche(0)))
      indicen=dsqrt(dreal(epscouche(neps+1)))

      if (ncote.eq.0.or.ncote.eq.1) then
         if (dreal(epscouche(neps+1)).le.0.d0) then
            infostr='epsilon <0 higher layer computation far field'
            nstop=-1
            return
         endif
         if (dimag(epscouche(neps+1)).gt.0.d0) then
            infostr='absorbing superstrate to computate far field'
            nstop=-1
            return
         endif
         indicem=indicen
         numapertra=numapertra/indicen
      endif

      if (ncote.eq.0.or.ncote.eq.-1) then     
         if (dreal(epscouche(0)).le.0.d0) then
            infostr='epsilon <0 higher layer computation far field'
            nstop=-1
            return
         endif
         if (dimag(epscouche(0)).gt.0.d0) then
            infostr='absorbing substrate to computate far field'
            nstop=-1
            return
         endif
         indicem=indice0
         numaperinc=numaperinc/indice0
         numaperref=numaperref/indice0
      endif
      
      if (ncote.eq.0) then
         indicem=max(indice0,indicen)
      endif

c     change ouverture numerique
      
      if (nlentille.eq.1) then
         
         if (beam(1:9).eq.'arbitrary') then
            infostr='Can not compute total field with this beam'
            nstop=1
            return
         endif
         if (numaperref.le.0.d0.or.numaperref.gt.1.d0) then
            nstop=1
            infostr='problem with numerical aperture in reflexion!'
            return
         endif
         if (numapertra.le.0.d0.or.numapertra.gt.1.d0) then
            nstop=1
            infostr='problem with numerical aperture in transmission!'
            return
         endif
         if (numaperinc.le.0.d0.or.numaperinc.gt.1.d0) then
            nstop=1
            infostr='problem with condenser numerical aperture!'
            return
         endif
         zlensr=zlensr*1.d-9
         zlenst=zlenst*1.d-9
         
      endif

      
 

c     wavenumber
      k0=2.d0*pi/lambda
      k03=k0*k0*k0
c     look  for compute near field with FFT
c     décrément de 1 de nproche pour faciliter le code C
      nproche=nproche-1
      write(*,*) 'number maximum of subunit :nmax = ',nmax     
      write(*,*) 'number of layer for the object :nnnr = ',nnnr

      write(*,*) '************** END INPUT DATA *******************'
      write(*,*) ' '

      write(*,*) '*************************************************'      
      write(*,*) '**************** BEGIN OBJECT *******************'
      write(*,*) '*************************************************'
c     look  for compute near field with FFT
      if (nlecture.eq.1.and.nproche.eq.-1) nproche=0
      
      
      if (nquickdiffracte.eq.1.and.nproche.eq.-1) nproche=0
      if (beam(1:5).eq.'gwave'.and.nproche.eq.-1) nproche=0
      nprochefft=0

c     test si wide field demandé
      if (nproche.ge.1) then
         if (nx.ne.nxm.or.ny.ne.nym.or.nz.ne.nzm) then
            nprochefft=nproche
            nproche=0
         else
            nproche=0
         endif
      endif

      if (nobjet.eq.-1.and.nproche.ge.1) nobjet=0
      
      if (nstop.eq.1) return

c     Built the object
      if (object(1:6).eq.'sphere') then
         numberobjet=1
         call objetspheresurf(trope,eps,epsani,xs,ys,zs,xswf,yswf,zswf
     $        ,k0 ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz
     $        ,polarizability ,nproche,epsilon,polarisa,rayon,xg,yg,zg
     $        ,neps,nepsmax ,dcouche ,zcouche,epscouche,tabzn,nmatf
     $        ,file_id,group_iddip,infostr ,nstop)
         write(99,*) 'sphere',rayon

      elseif (object(1:12).eq.'inhomosphere') then
         numberobjet=1
         if (trope.ne.'iso') then
            nstop=1
            infostr='Permittivity not scalar for inhomogenous sphere'
            write(99,*)
     $           'Permittivity not scalar for inhomogenous sphere'
            return
         endif
         call objetsphereinhomosurf(eps,xs,ys,zs,xswf,yswf,zswf ,k0
     $        ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz
     $        ,polarizability,nproche,epsilon,polarisa,rayon,lc,hc,ng
     $        ,localfieldx,neps,nepsmax,dcouche,zcouche,epscouche,tabzn
     $        ,nmatf,file_id,group_iddip,infostr,nstop)
         localfieldx=0.d0
         write(99,*) 'sphere',rayon
      elseif (object(1:13).eq.'inhomocuboid1') then
         numberobjet=1
         if (trope.ne.'iso') then
            nstop=1
            infostr='Permittivity not scalar for inhomogenous cuboid'
            write(99,*)
     $           'Permittivity not scalar for inhomogenous cuboid'
            return
         endif
         call objetparainhomosurf(eps,xs,ys,zs,xswf,yswf,zswf ,k0
     $        ,aretecube ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm
     $        ,nym ,nzm,polarizability ,nproche ,epsilon,polarisa,sidex
     $        ,sidey ,sidez ,xg ,yg ,zg,lc ,hc,ng,localfieldx,neps
     $        ,nepsmax ,dcouche ,zcouche ,epscouche ,tabzn ,nmatf
     $        ,file_id,group_iddip,infostr,nstop)

         localfieldx=0.d0
         write(99,*) 'cuboid',sidex,sidey,sidez
      elseif (object(1:13).eq.'inhomocuboid2') then
         numberobjet=1
         if (trope.ne.'iso') then
            nstop=1
            infostr='Permittivity not scalar for inhomogenous cuboid'
            write(99,*)
     $           'Permittivity not scalar for inhomogenous cuboid'
            return
         endif
         call objetparanxnynzinhomosurf(eps,xs,ys,zs,xswf,yswf,zswf ,k0
     $        ,aretecube ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm
     $        ,nym,nzm,nxmp,nymp,nzmp,polarizability ,nproche ,epsilon
     $        ,polarisa ,sidex ,sidey,sidez ,xg ,yg,zg,lc ,hc,ng
     $        ,localfieldx,neps ,nepsmax ,dcouche ,zcouche ,epscouche
     $        ,tabzn ,nmatf,file_id,group_iddip,infostr ,nstop)
         localfieldx=0.d0
         write(99,*) 'cuboid',sidex,sidey,sidez
      elseif (object(1:4).eq.'cube') then
         write(99,*) 'cube:side',side
         numberobjet=1
         call objetcubesurf(trope,eps,epsani ,eps0,xs,ys,zs,xswf,yswf
     $        ,zswf,k0,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny
     $        ,nz,polarizability,epsilon,polarisa,side,xg,yg,zg,neps
     $        ,nepsmax ,dcouche,zcouche,epscouche,tabzn,nmatf,file_id
     $        ,group_iddip,infostr ,nstop)

      elseif(object(1:7).eq.'cuboid1') then
         write(99,*) 'cuboid:sidex,sidey,sizez ',sidex,sidey,sidez
         numberobjet=1
         call objetparasurf(trope,eps,epsani,eps0,xs,ys,zs,xswf,yswf
     $        ,zswf,k0,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny
     $        ,nz,nxm,nym,nzm,polarizability,epsilon,polarisa,sidex
     $        ,sidey,sidez,xg ,yg,zg ,phiobj,thetaobj,psiobj,nproche
     $        ,neps,nepsmax ,dcouche ,zcouche,epscouche,tabzn,nmatf
     $        ,file_id,group_iddip,infostr ,nstop)
         write(99,*) 'side',sidex,sidey,sidez,nx,ny,nz
       
       
      elseif(object(1:7).eq.'cuboid2') then
         numberobjet=1

         call objetparanxnynzsurf(trope,eps,epsani,eps0,xs,ys,zs,xswf
     $        ,yswf,zswf,k0 ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole
     $        ,nx ,ny,nz,nxm,nym,nzm,nxmp,nymp,nzmp, polarizability
     $        ,epsilon ,polarisa ,sidex,sidey,sidez,xg,yg,zg,nproche
     $        ,neps ,nepsmax ,dcouche,zcouche ,epscouche,tabzn,nmatf
     $        ,file_id,group_iddip,infostr ,nstop)
         write(99,*) 'side',sidex,sidey,sidez,nx,ny,nz
      elseif(object(1:9).eq.'ellipsoid') then
         numberobjet=1
         call objetellipsesurf(trope,eps,epsani,eps0,xs,ys,zs,xswf,yswf
     $        ,zswf,k0 ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx
     $        ,ny,nz,nxm,nym,nzm,polarizability,nproche,epsilon,polarisa
     $        ,demiaxea,demiaxeb ,demiaxec,xg,yg,zg,phiobj,thetaobj
     $        ,psiobj,neps,nepsmax ,dcouche ,zcouche,epscouche,tabzn
     $        ,nmatf,file_id,group_iddip,infostr,nstop)
         write(99,*) 'ellipse',nbsphere,ndipole,nx,ny,nz
      elseif(object(1:8).eq.'nspheres') then
         call objetnspheressurf(trope,epsmulti,epsanimulti,numberobjet
     $        ,numberobjetmax,xgmulti,ygmulti,zgmulti,rayonmulti,eps0,xs
     $        ,ys,zs,xswf,yswf,zswf,k0 ,aretecube,tabdip,tabmulti,nnnr
     $        ,nmax,nbsphere ,ndipole,nx,ny,nz,polarizability,nproche
     $        ,epsilon ,polarisa,neps,nepsmax,dcouche ,zcouche,epscouche
     $        ,tabzn,nmatf,file_id,group_iddip,infostr ,nstop)
         write(99,*) 'multisphere',nbsphere,ndipole,nx,ny,nz
      elseif(object(1:8).eq.'cylinder') then
         numberobjet=1

         call objetcylindresurf(trope,eps,epsani,eps0,xs,ys,zs,xswf,yswf
     $        ,zswf,k0 ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx
     $        ,ny,nz,nxm,nym,nzm,polarizability,nproche,epsilon,polarisa
     $        ,rayon ,hauteur,xg ,yg,zg,phiobj,thetaobj,psiobj,neps
     $        ,nepsmax ,dcouche ,zcouche ,epscouche,tabzn,nmatf,file_id
     $        ,group_iddip,infostr ,nstop)
         write(99,*) 'cylindre',nbsphere,ndipole,nx,ny,nz
      elseif(object(1:16).eq.'concentricsphere') then
         write(*,*) epsmulti,numberobjet
         call objetsphereconcentricsurf(trope,epsmulti,epsanimulti
     $        ,numberobjet,numberobjetmax,xg,yg,zg,rayonmulti,eps0,xs,ys
     $        ,zs,xswf,yswf,zswf,k0,aretecube,tabdip,tabmulti,nnnr,nmax
     $        ,nbsphere ,ndipole,nx,ny,nz,polarizability,nproche,epsilon
     $        ,polarisa,neps ,nepsmax,dcouche ,zcouche ,epscouche,tabzn
     $        ,nmatf,file_id,group_iddip,infostr,nstop)
         write(99,*) 'concentricsphere',nbsphere,ndipole,nx,ny,nz
      elseif(object(1:9).eq.'arbitrary') then
         numberobjet=1
         call objetarbitrarysurf(trope,eps,epsani,eps0,xs,ys,zs,xswf
     $        ,yswf,zswf,k0,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole
     $        ,nx,ny,nz,polarizability,namefileobj,nproche,epsilon
     $        ,polarisa,neps,nepsmax ,dcouche,zcouche,epscouche,tabzn
     $        ,nmatf,file_id,group_iddip,infostr,nstop)
         write(99,*) 'arbitrary',nbsphere,ndipole,nx,ny,nz
      else
         write(99,*) 'Object unknown'
         write(*,*) 'object',trim(object)
         infostr='Object unknown'
         nstop=1
         return
      endif
      write(*,*) '************ END OBJECT *************************'
      write(*,*) ' '
      if (nstop.eq.1) return

      write(*,*) '*************************************************'      
      write(*,*) '************ NEW MULTILAYER SYSTEM **************'
      write(*,*) '*************************************************'

      do i=0,neps
         write(*,*) 'Epsilon layer',epscouche(i)
         write(*,*) '---------- z=',zcouche(i)         
      enddo
      write(*,*) 'Epsilon layer',epscouche(neps+1)

      write(*,*) '************ END MULTILAYER *********************'
      write(*,*) ' '
      
      if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
         nstop=1
         infostr='Dimension Problem of the Box'
         write(99,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         return
      endif
      if (nstop.eq.1) return
      if (nstop == -1) then
         infostr = 'Calculation cancelled after object created'
         return
      endif

c     ecriture dans fichiers du epsilon
      if (nmatf.eq.0) then
         if (trope.eq.'iso') then
            do i=1,ndipole
               k=tabdip(i)
               if (k.ne.0) then
                  write(66,*) dreal(epsilon(k,2,2)),dimag(epsilon(k,2
     $                 ,2))
               else
                  write(66,*) 1.d0 , 0.d0
               endif
            enddo
         else
            do i=1,ndipole
               k=tabdip(i)
               if (k.ne.0) then
                  do ii=1,3
                     do jj=1,3
                        write(66,*) dreal(epsilon(k,ii,jj))
     $                       ,dimag(epsilon(k,ii,jj))
                     enddo
                  enddo
               else
                  do ii=1,3
                     do jj=1,3
                        write(66,*) 1.d0 , 0.d0
                     enddo
                  enddo
               endif
            enddo          
         endif
         close(66)
      elseif (nmatf.eq.2) then
         dim(1)=1
         dim(2)=1
         datasetname='nx'
         call hdf5write1d_int(group_iddip,datasetname,nx,dim)
         datasetname='ny'
         call hdf5write1d_int(group_iddip,datasetname,ny,dim)
         datasetname='nz'
         call hdf5write1d_int(group_iddip,datasetname,nz,dim)

         if (trope.eq.'iso') then
            do i=1,ndipole
               k=tabdip(i)
               if (k.ne.0) then
                  wrk(i,1)=epsilon(k,2,2)
               else
                  wrk(i,1)=(1.d0,0.d0)
               endif
            enddo
            dim(1)=ndipole
            dim(2)=nmax*3
            datasetname='Epsilon real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,1))
     $           ,dim)
            datasetname='Epsilon imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,1))
     $           ,dim)
         else
            do i=1,ndipole
               k=tabdip(i)
               if (k.ne.0) then
                  wrk(i,1)=epsilon(k,1,1)
                  wrk(i,2)=epsilon(k,2,2)
                  wrk(i,3)=epsilon(k,3,3)
                  wrk(i,4)=epsilon(k,1,2)
                  wrk(i,5)=epsilon(k,1,3)
                  wrk(i,6)=epsilon(k,2,1)
                  wrk(i,7)=epsilon(k,2,3)
                  wrk(i,8)=epsilon(k,3,1)
                  wrk(i,9)=epsilon(k,3,2)
               else
                  wrk(i,1)=(1.d0,0.d0)
                  wrk(i,2)=(1.d0,0.d0)
                  wrk(i,3)=(1.d0,0.d0)
                  wrk(i,4)=(0.d0,0.d0)
                  wrk(i,5)=(0.d0,0.d0)
                  wrk(i,6)=(0.d0,0.d0)
                  wrk(i,7)=(0.d0,0.d0)
                  wrk(i,8)=(0.d0,0.d0)
                  wrk(i,9)=(0.d0,0.d0)
               endif
            enddo   

            dim(1)=ndipole
            dim(2)=nmax*3
            datasetname='epsilon xx real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,1))
     $           ,dim)
            datasetname='epsilon xx imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,1))
     $           ,dim)
            datasetname='epsilon yy real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,2))
     $           ,dim)
            datasetname='epsilon yy imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,2))
     $           ,dim)
            datasetname='epsilon zz real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,3))
     $           ,dim)
            datasetname='epsilon zz imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,3))
     $           ,dim)
            datasetname='epsilon xy real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,4))
     $           ,dim)
            datasetname='epsilon xy imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,4))
     $           ,dim)
            datasetname='epsilon xz real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,5))
     $           ,dim)
            datasetname='epsilon xz imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,5))
     $           ,dim)
            datasetname='epsilon yx real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,6))
     $           ,dim)
            datasetname='epsilon yx imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,6))
     $           ,dim)
            datasetname='epsilon yz real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,7))
     $           ,dim)
            datasetname='epsilon yz imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,7))
     $           ,dim)
            datasetname='epsilon zx real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,8))
     $           ,dim)
            datasetname='epsilon zx imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,8))
     $           ,dim)
            datasetname='epsilon zy real part'
            call hdf5write1d(group_iddip,datasetname,dreal(wrk(:,9))
     $           ,dim)
            datasetname='epsilon zy imaginary part'
            call hdf5write1d(group_iddip,datasetname,dimag(wrk(:,9 ))
     $           ,dim)
         endif
      endif
c     fin ecriture du epsilon
c     epsilon
      
      lambda10n = lambda/10.d0/cdabs(cdsqrt(eps))
      if (aretecube.ge.2.5d0*lambda10n) then
         nstop=1
         infostr='meshsize larger than lambda/4'
         return
      endif
      

c     calcul le imaxk0 et le deltakx et deltax pour energy, microscopy etc...
      if (nenergie+nlentille.ge.1) then
         write(*,*) '**************************************************'
         write(*,*) '*********** INITIALIZE NB OF PT IN NA ************'
         write(*,*) '**************************************************'
         if (nquickdiffracte.eq.0) then
            write(*,*) 'Computation of delta k and delta x'
            write(*,*) 'for the diffracted field with slow method'
            write(*,*) 'Initial step size delta k : ',2.d0*pi
     $           /(dble(nfft2d)*aretecube),'m-1'
            k=0
 222        deltakx=2.d0*pi/(dble(nfft2d)*aretecube)/dble(2**k)
            if (ncote.eq.0) then 
               imaxk0=max(nint(k0*indicen/deltakx)+1,nint(k0*indice0
     $              /deltakx)+1)
            elseif (ncote.eq.1) then
               imaxk0=nint(numapertra*k0*indicen/deltakx)+1
            elseif (ncote.eq.-1) then
               imaxk0=nint(numaperref*k0*indice0/deltakx)+1
            endif
            
            if (imaxk0.le.20) then
               k=k+1
               write(*,*) 'Change delta k :',k,dble(nfft2d*(2**k))
     $              ,nfft2d,imaxk0
               goto 222
            endif
            write(*,*) 'Final delta k',deltakx,'m-1'
            write(*,*) 'Number of point in NA',2*imaxk0+1
            deltaky=deltakx
            deltax=2.d0*pi/dble(nfft2d)/deltakx
            
         else
            deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
            deltaky=deltakx
            if (ncote.eq.0) then 
               imaxk0=max(nint(k0*indicen/deltakx)+1,nint(k0*indice0
     $              /deltakx)+1)
            elseif (ncote.eq.1) then
               imaxk0=nint(numapertra*k0*indicen/deltakx)+1
            elseif (ncote.eq.-1) then
               imaxk0=nint(numaperref*k0*indice0/deltakx)+1
            endif
            write(*,*) 'Delta k',deltakx,'m-1'
            write(*,*) 'Number of point in NA',2*imaxk0+1
            if (2*imaxk0+1.gt.nfft2d) then
               infostr
     $       ='Size of FFT too small to compute the diffracted field'
               nstop=-1
               return
            endif
            if (2*imaxk0+1.lt.7) then
               write(99,*) '2*imax+1',imaxk0,2*imaxk0+1,nfft2d
               infostr='In FFT diffract nfft2d too small'
               nstop = 1
               return
            endif
         endif
         write(*,*) ' *************************************************'
      endif

c     Changement angle si microscopy
      if (nlentille.eq.1) then
         write(*,*) '**************************************************'
         write(*,*) '********** INITIALIZE ANGLE OF INCIDENCE *********'
         write(*,*) '**************************************************'
         if (beam(1:11).eq.'pwavelinear' .or. beam(1:13).eq
     $        .'pwavecircular') then
            write(*,*) 'Angle of incidence asked:', theta,phi
            call anglecalculmic(theta,phi,indice0,k0,deltakx,deltaky)
            write(*,*) 'Angle of incidence taken:', theta,phi

            write(*,*) 'If the new values are too far from the previous'
            write(*,*) 'increase the size of the FFT'
         elseif (beam(1:15).eq.'wavelinearmulti') then
            do i=1,nbinc
               write(*,*) 'Numero of the incidence :',i,'/',nbinc
               write(*,*) 'Angle of incidence asked:', thetam(i),phim(i)
               call anglecalculmic(thetam(i),phim(i),indice0,k0,deltakx
     $              ,deltaky)
               write(*,*) 'Angle of incidence taken:', thetam(i),phim(i)
             
            enddo
            write(*,*) 'If the new values are too far from the previous'
            write(*,*) 'increase the size of the FFT'
         endif
         write(*,*) ' *************************************************'
      endif





      write(99,*) 'number of subunit for the object',nbsphere
      write(99,*) 'number of subunit for the mesh ',ndipole
      write(99,*) 'mesh size',aretecube
      write(99,*) 'lambda/(10n)',lambda/10.d0/cdabs(cdsqrt(eps))


      if (beam(1:11).eq.'pwavelinear') then
         write(99,*) 'Beam pwavelinear'
      elseif (beam(1:13).eq.'pwavecircular') then
         write(99,*) 'Beam pwavecircular'
      elseif (beam(1:15).eq.'wavelinearmulti') then
         write(99,*) 'Beam wavelinearmulti' 
      elseif (beam(1:11).eq.'gwavelinear') then
         write(99,*) 'Beam gwavelinear'    
      elseif (beam(1:13).eq.'gwavecircular') then
         write(99,*) 'Beam gwavecircular'
      elseif (beam(1:8).eq.'gwaveiso') then
         write(99,*) 'Beam isofocus'    
      elseif (beam(1:7).eq.'speckle') then
         write(99,*) 'Beam speckle'    
      elseif (beam(1:9).eq.'arbitrary') then
         write(99,*) 'Beam arbitrary'    
      else
         write(99,*) 'Beam unknown'
         write(*,*) 'beam',trim(beam)
         write(*,*) 'object',trim(object)
         infostr='Beam unknown'
         nstop=1
         return
      endif

c     cré le fichier de data pour connaitre les options pour matlab

      open(900,file='inputmatlab.mat')
      write(900,*) nproche
      write(900,*) nlocal
      write(900,*) nmacro
      write(900,*) nsection
      write(900,*) ndiffracte    
      write(900,*) nquickdiffracte
      write(900,*) nforce
      write(900,*) nforced
      write(900,*) ntorque
      write(900,*) ntorqued
      write(900,*) nlentille
      write(900,*) nquicklens
      write(900,*) nphi
      write(900,*) ntheta+1
      if (trope.eq.'iso') write(900,*) 0
      if (trope.eq.'ani') write(900,*) 1
      write(900,*) nfft2d
      write(900,*) k0
      write(900,*) numaperref*indice0
      write(900,*) numapertra*indicen
      write(900,*) nprochefft
      write(900,*) nobjet
      write(900,*) ncote
      write(900,*) indicen
      write(900,*) indice0
      write(900,*) ntypemic
      write(900,*) nmatf
      close(900)
      if (nmatf.eq.2) then
         open(901,file='filenameh5')
         write(901,*) h5file
         close(901)
         dim(1)=1
         dim(2)=1
         
         datasetname='nproche'
         call hdf5write1d_int(group_idopt,datasetname,nproche,dim)
         datasetname='nlocal'
         call hdf5write1d_int(group_idopt,datasetname,nlocal,dim)
         datasetname='nmacro'
         call hdf5write1d_int(group_idopt,datasetname,nmacro,dim)
         datasetname='nsection'
         call hdf5write1d_int(group_idopt,datasetname,nsection,dim)
         datasetname='ndiffracte'
         call hdf5write1d_int(group_idopt,datasetname,ndiffracte,dim)
         datasetname='nquickdiffracte'
         call hdf5write1d_int(group_idopt,datasetname,nquickdiffracte
     $        ,dim)
         datasetname='nforce'
         call hdf5write1d_int(group_idopt,datasetname,nforce,dim)
         datasetname='nforced'
         call hdf5write1d_int(group_idopt,datasetname,nforced,dim)
         datasetname='ntorque'
         call hdf5write1d_int(group_idopt,datasetname,ntorque,dim)
         datasetname='ntorqued'
         call hdf5write1d_int(group_idopt,datasetname,ntorqued,dim)
         datasetname='nlentille'
         call hdf5write1d_int(group_idopt,datasetname,nlentille,dim)
         datasetname='nquicklens'
         call hdf5write1d_int(group_idopt,datasetname,nquicklens,dim)
         datasetname='nphi'
         call hdf5write1d_int(group_idopt,datasetname,nphi,dim)
         datasetname='ntheta'
         call hdf5write1d_int(group_idopt,datasetname,ntheta+1,dim)
         
         if (trope.eq.'iso') then
            datasetname='iso'
            i=0
            call hdf5write1d_int(group_idopt,datasetname,i,dim)
         endif
         
         if (trope.eq.'ani') then
            datasetname='iso'
            i=1
            call hdf5write1d_int(group_idopt,datasetname,i,dim)
         endif
         datasetname='nfft2d'
         call hdf5write1d_int(group_idopt,datasetname,nfft2d,dim)
         datasetname='k0'
         call hdf5write1d(group_idopt,datasetname,k0,dim)
         datasetname='numaper reflexion'
         tmp=numaperref*indice0
         call hdf5write1d(group_idopt,datasetname,tmp,dim)
         datasetname='numaper transmission'
         tmp=numapertra*indicen
         call hdf5write1d(group_idopt,datasetname,tmp,dim)
         datasetname='nprochefft'
         call hdf5write1d_int(group_idopt,datasetname,nprochefft,dim)
         datasetname='nobjet'
         call hdf5write1d_int(group_idopt,datasetname,nobjet,dim)
         datasetname='nside'
         call hdf5write1d_int(group_idopt,datasetname,ncote,dim)
         datasetname='index upper'
         call hdf5write1d(group_idopt,datasetname,indicen,dim)
         datasetname='index lower'
         call hdf5write1d(group_idopt,datasetname,indice0,dim)
         datasetname='ntypemic'
         call hdf5write1d_int(group_idopt,datasetname,ntypemic,dim)
         datasetname='nmatf'
         call hdf5write1d_int(group_idopt,datasetname,nmatf,dim)
      endif

c     ne fait que l'objet
      if (nobjet.eq.1) then 
         infostr='Dipole calculation completed'
         goto 999
      endif

      write(*,*) '*************************************************'      
      write(*,*) '************** BEGIN INCIDENT FIELD *************'
      write(*,*) '*************************************************'
      write(*,*) 'Beam used  : ',trim(beam)
      write(99,*) '******* Compute the incident field *******'
      write(99,*) 'Beam used',trim(beam)
      write(99,*) 'k0=',k0     
      subunit=0
      write(*,*) 'Initialize plan for FFT'

      call dfftw_plan_dft_2d(plan2b,nfft2d,nfft2d,Eimagexpos,Eimagexpos
     $     ,FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_plan_dft_2d(plan2f,nfft2d,nfft2d,Eimagexpos,Eimagexpos
     $     ,FFTW_FORWARD,FFTW_ESTIMATE)

         
      if (beam(1:11).eq.'pwavelinear') then
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol1=',pp
         write(99,*) 'pol2=',ss
c     compute E0
         call irradiancesurf(P0,w0,E0,irra,epscouche(0))
         write(*,*) 'Power     :',P0
         write(*,*) 'Waist     :',w0
         write(*,*) 'Field     :',E0
         write(*,*) 'Irradiance:',irra
         I0=cdabs(E0)**2
c     write(*,*) 'champ',epscouche,zcouche,neps,nepsmax,k0,E0,ss,pp
c     $        ,theta,phi
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)         
         do i=1,nbsphere
            call champlineaire(epscouche,zcouche,neps,nepsmax,xs(i)
     $           ,ys(i),zs(i),k0,E0,ss,pp,theta,phi,infostr,nstop,FF0(3
     $           *i-2),FF0(3*i-1),FF0(3*i))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         if (nstop.eq.1) return
         
      elseif (beam(1:13).eq.'pwavecircular') then
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol=',ss
c     compute E0
         call irradiancesurf(P0,w0,E0,irra,epscouche(0))
         write(*,*) 'Power     :',P0
         write(*,*) 'Waist     :',w0
         write(*,*) 'Field     :',E0
         write(*,*) 'Irradiance:',irra

         I0=cdabs(E0)**2

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)          
         do i=1,nbsphere
            call champcirculaire(epscouche,zcouche,neps,nepsmax,xs(i)
     $           ,ys(i),zs(i),k0,E0,ss,theta,phi,infostr,nstop,FF0(3*i
     $           -2),FF0(3*i-1) ,FF0(3*i))
c            write(*,*) 'FF0',FF0(3*i -2),FF0(3*i-1) ,FF0(3*i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL           
         if (nstop.eq.1) return
         
      elseif (beam(1:15).eq.'wavelinearmulti') then
         do i=1,nbinc
            write(99,*) 'theta=',thetam(i),i,nbinc
            write(99,*) 'phi=',phim(i)
            write(99,*) 'pol1=',ppm(i)
            write(99,*) 'pol2=',ssm(i)
         enddo

c     compute E0
c     call irradiancesurf(P0,w0,E0,irra,epscouche(0))
c     write(*,*) 'irra',P0,w0,E0,irra
c     I0=cdabs(E0)**2
c     write(*,*) 'champ',epscouche,zcouche,neps,nepsmax,k0,E0,ss,pp
c     $        ,theta,phi
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)        
         do i=1,nbsphere
            call champlineairemulti(epscouche,zcouche,neps,nepsmax,xs(i)
     $           ,ys(i),zs(i),k0,E0m,ssm,ppm,thetam,phim,nbinc,infostr
     $           ,nstop,FF0(3 *i-2),FF0(3*i-1),FF0(3*i))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

         if (nstop.eq.1) return
         
      elseif  (beam(1:11).eq.'gwavelinear') then
         E0=1.d0
         psi=pp
         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
c         write(*,*) 'routine gauss',psi,ss,pp
c     write(*,*) 'data',epscouche,zcouche,neps,nepsmax,k0,w0,x0,y0,z0
c     $        ,theta,phi,psi,E0,xs,ys,zs,aretecube,nx,ny,nz,ndipole
c     $        ,nmax,nfft2d,nfft2d
         call gaussiansurf(epscouche,zcouche,neps,nepsmax,k0,w0,xgaus
     $        ,ygaus,zgaus,theta,phi,psi,E0,xs,ys,zs,FF0,aretecube,nx,ny
     $        ,nz,ndipole,nmax,nfft2d,nfft2d,Efourierincxneg
     $        ,Efourierincyneg ,Efourierinczneg,Efourierincxpos
     $        ,Efourierincypos,Efourierinczpos,fluxinc ,fluxref
     $        ,fluxtrans,irra ,nstop ,infostr,plan2b)
         if (nstop.eq.1) return
c     write(*,*) 'champ',FF0
c     remet tout a l'echelle pour respecter la consigne de la puissance
c     incidente
         write(*,*) 'Power',P0
         tmp=P0/fluxinc
c     tmp=1.d0
         fluxref=fluxref*tmp
         fluxtrans= fluxtrans*tmp
         irra=irra*tmp
         tmp=dsqrt(tmp)
         E0=E0*tmp
         write(*,*) 'Irradiance',irra,'field',E0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)         
         do i=1,nfft2d*nfft2d
            Efourierincxneg(i)=Efourierincxneg(i)*tmp
            Efourierincyneg(i)=Efourierincyneg(i)*tmp
            Efourierinczneg(i)=Efourierinczneg(i)*tmp
            Efourierincxpos(i)=Efourierincxpos(i)*tmp
            Efourierincypos(i)=Efourierincypos(i)*tmp
            Efourierinczpos(i)=Efourierinczpos(i)*tmp
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)        
         do i=1,nbsphere
            FF0(3*i-2)=FF0(3*i-2)*tmp
            FF0(3*i-1)=FF0(3*i-1)*tmp
            FF0(3*i)=FF0(3*i)*tmp
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
         fluxinc=P0
         write(*,*) 'fluxinc',fluxinc,fluxref,fluxtrans,(fluxref
     $        +fluxtrans)/fluxinc
         write(*,*) 'irra',P0,w0,dsqrt(cdabs(FF0(1))**2+ cdabs(FF0(2))
     $        **2+cdabs(FF0(3))**2),irra
         write(*,*) 'irra',irra,'I0',cdabs(FF0(1))**2+ cdabs(FF0(2))**2
     $        +cdabs(FF0(3))**2,'P02mu0c/S',8.d0*pi*1.d-7*299792458.d0
     $        *P0/(w0*w0*pi)
         I0=(cdabs(FF0(1))**2+ cdabs(FF0(2)) **2+cdabs(FF0(3))**2)
      elseif  (beam(1:13).eq.'gwavecircular') then
         E0=1.d0
         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
         call gaussiansurfcirc(epscouche,zcouche,neps,nepsmax,k0,w0
     $        ,xgaus,ygaus,zgaus,theta,phi,ss,E0,xs,ys,zs,FF0,aretecube
     $        ,nx,ny,nz,ndipole,nmax,nfft2d,nfft2d,Efourierincxneg
     $        ,Efourierincyneg ,Efourierinczneg,Efourierincxpos
     $        ,Efourierincypos,Efourierinczpos,fluxinc ,fluxref
     $        ,fluxtrans,irra,nstop ,infostr,plan2b)
         if (nstop.eq.1) return
c     remet tout a l'echelle pour respecter la consigne de la puissance
c     incidente
         tmp=P0/fluxinc
         fluxref=fluxref*tmp
         fluxtrans= fluxtrans*tmp
         irra=irra*tmp
         tmp=dsqrt(tmp)
         E0=E0*tmp

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)        
         do i=1,nfft2d*nfft2d
            Efourierincxneg(i)=Efourierincxneg(i)*tmp
            Efourierincyneg(i)=Efourierincyneg(i)*tmp
            Efourierinczneg(i)=Efourierinczneg(i)*tmp
            Efourierincxpos(i)=Efourierincxpos(i)*tmp
            Efourierincypos(i)=Efourierincypos(i)*tmp
            Efourierinczpos(i)=Efourierinczpos(i)*tmp
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)             
         do i=1,nbsphere
            FF0(3*i-2)=FF0(3*i-2)*tmp
            FF0(3*i-1)=FF0(3*i-1)*tmp
            FF0(3*i)=FF0(3*i)*tmp
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
         fluxinc=P0
         write(*,*) 'fluxinc',fluxinc,fluxref,fluxtrans,(fluxref
     $        +fluxtrans)/fluxinc
         write(*,*) 'irra',P0,w0,dsqrt(cdabs(FF0(1))**2+ cdabs(FF0(2))
     $        **2+cdabs(FF0(3))**2),irra
         write(*,*) 'irra',irra,'I0',cdabs(FF0(1))**2+ cdabs(FF0(2))**2
     $        +cdabs(FF0(3))**2,'P02mu0c/S',8.d0*pi*1.d-7*299792458.d0
     $        *P0/(w0*w0*pi)
         I0=(cdabs(FF0(1))**2+ cdabs(FF0(2)) **2+cdabs(FF0(3))**2)

      elseif  (beam(1:8).eq.'gwaveiso') then
         E0=1.d0
         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
         psi=pp
         write(*,*) 'routine iso',psi,ss,pp
c     write(*,*) 'data',epscouche,zcouche,neps,nepsmax,k0,w0,x0,y0,z0
c     $        ,theta,phi,psi,E0,xs,ys,zs,aretecube,nx,ny,nz,ndipole
c     $        ,nmax,nfft2d,nfft2d

         nmasque=1
         
         call gaussianiso(epscouche,zcouche,neps,nepsmax,k0,w0,xgaus
     $        ,ygaus,zgaus,psi,E0,xs,ys,zs,FF0,aretecube,nx,ny,nz
     $        ,ndipole ,nmax ,nfft2d,nfft2d,Efourierincxneg
     $        ,Efourierincyneg ,Efourierinczneg,Efourierincxpos
     $        ,Efourierincypos,Efourierinczpos,fluxinc ,fluxref
     $        ,fluxtrans,irra,numaperref,masque,nmasque,nstop,infostr
     $        ,plan2b)
         write(*,*) infostr
         if (nstop.eq.1) return
c     write(*,*) 'champ',FF0
c     write(*,*) '000a',Eimageincx(524801) ,Eimageincy(524801)
c     $        ,Eimageincz(524801)
c     remet tout a l'echelle pour respecter la consigne de la puissance
c     incidente
         write(*,*) 'fluxinc apres routine',fluxinc
         tmp=P0/fluxinc
c     tmp=1.d0
         fluxref=fluxref*tmp
         fluxtrans= fluxtrans*tmp
         irra=irra*tmp
         tmp=dsqrt(tmp)
         E0=E0*tmp

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)           
         do i=1,nfft2d*nfft2d
            Efourierincxneg(i)=Efourierincxneg(i)*tmp
            Efourierincyneg(i)=Efourierincyneg(i)*tmp
            Efourierinczneg(i)=Efourierinczneg(i)*tmp
            Efourierincxpos(i)=Efourierincxpos(i)*tmp
            Efourierincypos(i)=Efourierincypos(i)*tmp
            Efourierinczpos(i)=Efourierinczpos(i)*tmp
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            FF0(3*i-2)=FF0(3*i-2)*tmp
            FF0(3*i-1)=FF0(3*i-1)*tmp
            FF0(3*i)=FF0(3*i)*tmp
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         
         fluxinc=P0
         write(*,*) 'fluxinc',fluxinc,fluxref,fluxtrans,(fluxref
     $        +fluxtrans)/fluxinc
         write(*,*) 'irra',P0,w0,dsqrt(cdabs(FF0(1))**2+ cdabs(FF0(2))
     $        **2+cdabs(FF0(3))**2),irra
         write(*,*) 'irra',irra,'I0',cdabs(FF0(1))**2+ cdabs(FF0(2))**2
     $        +cdabs(FF0(3))**2,'P02mu0c/S',8.d0*pi*1.d-7*299792458.d0
     $        *P0/(w0*w0*pi)
         I0=(cdabs(FF0(1))**2+ cdabs(FF0(2)) **2+cdabs(FF0(3))**2)


      elseif  (beam(1:7).eq.'speckle') then
         E0=1.d0
         psi=pp
         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
         call specklesurf(epscouche,zcouche,neps,nepsmax,k0,E0
     $        ,numaperref,IR,xs,ys,zs,xgaus ,ygaus,zgaus,psi,FF0
     $        ,aretecube,nx,ny,nz,nxm,nym,nzm ,ndipole ,nmax ,nfft2d
     $        ,nfft2d,Efourierincxneg,Efourierincyneg ,Efourierinczneg
     $        ,Efourierincxpos,Efourierincypos ,Efourierinczpos,fluxinc
     $        ,fluxref,fluxtrans,irra ,nstop ,infostr,plan2b)
         if (nstop.eq.1) return
         
      elseif  (beam(1:9).eq.'arbitrary') then
         call incidentarbitrary(xs,ys,zs,aretecube,FF0,nxm,nym,nzm
     $        ,nbsphere,nstop,namefileinc,infostr)
         if (nstop.eq.1) return
         
      endif
      write(*,*) '**************** END INCIDENT *******************'
      write(*,*) ' '
      
c     ecriture dans .mat et creation de incidentifeld et composantes.
      subunit=0
      if (nmatf.eq.0) then
         do i=1,ndipole         
            k=tabdip(i)
c     write(*,*) 'iii',i,ndipole,k
            if (k.ne.0) then
               subunit=subunit+1
               incidentfieldx(subunit)=FF0(3*k-2)
               incidentfieldy(subunit)=FF0(3*k-1)
               incidentfieldz(subunit)=FF0(3*k)
               incidentfield(subunit)=dsqrt(dreal(FF0(3*k-2)
     $              *dconjg(FF0(3*k-2))+FF0(3*k-1)*dconjg(FF0(3*k-1))
     $              +FF0(3*k)*dconjg(FF0(3*k))))
               write(36,*)dreal(incidentfieldx(subunit))
     $              ,dimag(incidentfieldx(subunit))
               write(37,*)dreal(incidentfieldy(subunit))
     $              ,dimag(incidentfieldy(subunit))
               write(38,*)dreal(incidentfieldz(subunit))
     $              ,dimag(incidentfieldz(subunit))
               write(39,*) incidentfield(subunit)
            else
               write(36,*) 0.d0,0.d0
               write(37,*) 0.d0,0.d0
               write(38,*) 0.d0,0.d0
               write(39,*) 0.d0
            endif            
         enddo
      elseif (nmatf.eq.2) then
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC)         
         do i=1,ndipole
            k=tabdip(i)
            if (k.ne.0) then
               incidentfieldx(k)=FF0(3*k-2)
               incidentfieldy(k)=FF0(3*k-1)
               incidentfieldz(k)=FF0(3*k)
               incidentfield(k)=dsqrt(dreal(FF0(3*k-2) *dconjg(FF0(3*k
     $              -2))+FF0(3*k-1)*dconjg(FF0(3*k-1)) +FF0(3*k)
     $              *dconjg(FF0(3*k))))
               wrk(i,1)=FF0(3*k-2)
               wrk(i,2)=FF0(3*k-1)
               wrk(i,3)=FF0(3*k)
               wrk(i,4)=incidentfield(k)
            else
               wrk(i,1)=0.d0
               wrk(i,2)=0.d0
               wrk(i,3)=0.d0
               wrk(i,4)=0.d0
            endif            
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL         
         dim(1)=ndipole
         dim(2)=nmax*3
         datasetname='Incident field modulus'
         call hdf5write1d(group_idnf,datasetname,dreal(wrk(:,4)), dim)
         datasetname='Incident field x component real part'
         call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,1)),dim)
         datasetname='Incident field x component imaginary part'
         call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,1)),dim)
         datasetname='Incident field y component real part'
         call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,2)),dim)
         datasetname='Incident field y component imaginary part'
         call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,2)),dim)
         datasetname='Incident field z component real part'
         call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,3)),dim)
         datasetname='Incident field z component imaginary part'
         call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,3)),dim)

      else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC)         
         do i=1,ndipole
            k=tabdip(i)
            if (k.ne.0) then
               incidentfieldx(k)=FF0(3*k-2)
               incidentfieldy(k)=FF0(3*k-1)
               incidentfieldz(k)=FF0(3*k)
               incidentfield(k)=dsqrt(dreal(FF0(3*k-2) *dconjg(FF0(3*k
     $              -2))+FF0(3*k-1)*dconjg(FF0(3*k-1)) +FF0(3*k)
     $              *dconjg(FF0(3*k))))            
            endif            
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
      endif


      
c     close Intensity of the incident field
      close(36)
      close(37)
      close(38)
      close(39)

      write(99,*) 'Field modulus',cdabs(E0)
      if (nstop.eq.1) return
      if (nstop == -1) then
         infostr = 'Calculation cancelled after incident field'
         return
      endif

c     multiplication by a factor 2: Toeplitz matrix transformed in a
c     circulant matrix with a doble size.
      nbsphere3=3*nbsphere
      nx2=2*nx
      ny2=2*ny
      nz2=2*nz
      nxy2=nx2*ny2
      ntotal=8*nx*ny*nz      

      if (nlecture.eq.1) then
c         write(*,*) 'relecture',filereread
         call relecturesurf(lambda,beam,object,trope,nnnr,tolinit, side,
     $        sidex, sidey, sidez, hauteur, numberobjet, rayonmulti,
     $        xgmulti, ygmulti, zgmulti, epsmulti, epsanimulti, demiaxea
     $        ,demiaxeb,demiaxec,thetaobj,phiobj,psiobj, namefileobj,
     $        theta, phi, pp, ss, thetam ,phim, ppm, ssm,E0m, nbinc, P0,
     $        w0,nrig, ninterp,ir,xgaus,ygaus,zgaus,namefileinc,
     $        numberobjetmax ,filereread,nlecture1,neps,nepsmax,zcouche
     $        ,epscouche,aretecube ,nx,ny,nz,hc,lc,ng,nstop ,infostr)
         if (nstop.eq.1) return
c     reread the local field         
         if (nlecture1.eq.1) then
            file1='.lf'
            long = len( trim(filereread  ) )
            long1 = len( trim( file1 ) )
            filereread1=filereread(1:long)//file1(1:long1)
            open(1000,file=filereread1,status='old',iostat=ierror,form
     $           ='unformatted')
            if (ierror.ne.0) then
               nstop=1
               infostr='file for read local field do not exist'
               return
            else
               do i=1,nbsphere3
                  read(1000) FFloc(i)
               enddo

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC)
            do i=1,nbsphere
               k=3*(i-1)
               FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1,2)
     $              *FFloc(k+2)+polarisa(i,1,3)*FFloc(k+3)
               FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2,2)
     $              *FFloc(k+2)+polarisa(i,2,3)*FFloc(k+3)
               FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3,2)
     $              *FFloc(k+2)+polarisa(i,3,3)*FFloc(k+3)            
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL

            endif
            close(1000)
            goto 1000
         endif
      endif

      

c     if the computation asked is rigourous then compute the Green
c     function and its FFT
      if (nrig.eq.0.or.nrig.eq.3) then

         write(*,*) '*************************************************'      
         write(*,*) '****** INITIALIZE PLAN FOR FFT ******************'
         write(*,*) '*************************************************'
         call dfftw_plan_dft_2d(planb,nx2,ny2,b11,b11,FFTW_BACKWARD
     $        ,FFTW_ESTIMATE)
         call dfftw_plan_dft_2d(planf,nx2,ny2,b11,b11,FFTW_FORWARD
     $        ,FFTW_ESTIMATE)

         write(*,*) '************* END PLAN **************************'
         write(*,*) ' '

         hcc=0.3d0
         epsabs=0.d0
         nt=1
c     write(*,*) 'gree surf',hcc,tolinit,epsabs,aretecube,k0,neps
c     $        ,nepsmax,dcouche,zcouche,epscouche,nbsphere ,nmax,n1m
c     $        ,nplanm,nz,nbs,nmat,nmatim,nt 
c     compute Green function (local field)
         if (nobjet.eq.-1.or.nlecture1.eq.1) then
            write(*,*) '**********************************************'      
            write(*,*) '**** DO NOT COMPUTE GREEN FUNCTION ***********'
            write(*,*) '**********************************************'
            goto 200
         endif
         write(*,*) '*************************************************'      
         write(*,*) '**************** BEGIN GREEN FUNCTION ***********'
         write(*,*) '*************************************************'
         write(*,*) 'Interpolation',ninterp  
         call cpu_time(t1)
         call date_and_time(date,time,zone,values)
         if (ninterp.eq.0) then

c     write(*,*) 'fonction interp',hcc,tolinit,epsabs,aretecube,k0
c     $           ,neps,nepsmax,nbsphere ,nmax ,n1m,nzm,nz,nbs,nmat
c     $           ,nmatim,nplanm
            call fonctiongreensurfcomp(hcc,tolinit,epsabs,xswf,yswf,zswf
     $           ,aretecube,k0,neps,nepsmax,dcouche,zcouche,epscouche
     $           ,ndipole ,nmax,n1m,nzm,nz,nbs,nmat,nmatim,nplanm,Tabzn
     $           ,a ,matind ,matindplan,matindice,matrange,nt)
c     write(*,*) 'xyz',xs(1),ys(1),zs(1)
c     compte FFT of the Green function
            write(*,*) '******* BEGIN FFT of GREEN FUNCTION **********'
            call fonctiongreensurffft(nx,ny,nz,nx2,ny2,nxm,nym,n1m,nzm
     $           ,nplanm,nmatim,nbs,ntotalm,aretecube,a,matind
     $           ,matindplan,matindice,matrange,b11,b12,b13,b22 ,b23
     $           ,b31,b32,b33,a11,a12,a13,a22,a23,a31,a32,a33,planb)
         else
c            write(*,*) 'fonction interp',hcc,tolinit,epsabs,nx,ny ,nz
c     $           ,aretecube,k0,neps,nepsmax,dcouche,zcouche ,epscouche
c     $           ,ndipole,nmax,n1m,nzm,nz,nbs,nmat,nmatim ,nplanm

            call  fonctiongreensurfcompinterp(hcc,tolinit,epsabs,nx,ny
     $           ,nz,zswf,aretecube,k0,neps,nepsmax,dcouche,zcouche
     $           ,epscouche,ndipole,nmax,n1m,nzm,nz,nbs,nmat,nmatim
     $           ,nplanm,Tabzn,a,matind,matindplan,matindice,matrange
     $           ,ninter ,ninterp,nt)
            write(*,*) '******* BEGIN FFT of GREEN FUNCTION **********'
            call fonctiongreensurfinterpfft(nx,ny,nz,nx2,ny2,nxm,nym,n1m
     $           ,nzm,nplanm,nmatim,nbs,ntotalm,ninter,ninterp,aretecube
     $           ,a ,matind ,matindplan ,matindice ,matrange,b11,b12,b13
     $           ,b22,b23 ,b31 ,b32,b33,a11 ,a12,a13 ,a22,a23,a31 ,a32
     $           ,a33,planb)
c     write(*,*) 'a11',a11
         endif
         call cpu_time(t2)
         call date_and_time(date,time,zone,values2)
         message='to compute Green function'
         call calculatedate(values2,values,t2,t1,message)
         write(*,*) '**************** END GREEN FUNCTION *************'
         write(*,*) ' '
         
c     Compute the incident field at each subunit of the object
         if (nstop == -1) then
            infostr = 'Calculation cancelled after FFT Green function'
            return
         endif
      endif

      
 200  if (nrig.eq.0) then

         write(*,*) '*************************************************'      
         write(*,*) '************* SOLVE LINEAR SYSTEM ***************'
         write(*,*) '*************************************************'
         
         write(99,*) '***** Solve the linear system *****'
c     Compute the local field at each subunit position by solving Ax=b
         
c     as an initial guess the incident field
         if (FFloc(1).eq.0.d0) then
            do i=1,nbsphere3
               xi(i)=FF0(i)
            enddo
         else
            do i=1,nbsphere3
               xi(i) = FFloc(i)
            enddo
         endif
c     initilization for solve Ax=b
         tol=tolinit

         if (nproche.eq.-1) then
c            write(*,*) 'no opt'
            call inverserigsurf(xi,xr,nbsphere,ndipole,nx,ny,nz,nx2,ny2
     $           ,nxm,nym,nzm,nplanm,ntotalm,nmax,matindplan,Tabdip,b31
     $           ,b32 ,b33,FF,FF0,FFloc,b11,b12,b13,a11,a12,a13,a22,a23
     $           ,a31 ,a32 ,a33 ,WRK,nlar,ldabi,polarisa,methodeit,tol
     $           ,tol1 ,nloop ,ncompte ,planf,planb,nstop ,infostr)
         else
c            write(*,*) 'opt'
            call inverserigsurfopt(xi,xr,nbsphere,ndipole,nx,ny,nz,nx2
     $           ,ny2,nxm,nym,nzm,nplanm,ntotalm,nmax,matindplan,b31
     $           ,b32 ,b33,FF,FF0,FFloc,b11,b12,b13,a11,a12,a13,a22,a23
     $           ,a31 ,a32,a33 ,WRK,nlar,ldabi,polarisa,methodeit,tol
     $           ,tol1 ,nloop,ncompte ,planf,planb,nstop ,infostr)
         endif      
         if (nstop.eq.1) return
         
         if (nlecture.eq.1.and.nlecture1.eq.0) then
            
            file1='.lf'
            long = len( trim(filereread  ) )
            long1 = len( trim( file1 ) )
            filereread1=filereread(1:long)//file1(1:long1)
            open(1000,file=filereread1,status='new',form='unformatted')
            do i=1,nbsphere3
               write(1000) FFloc(i)
            enddo
            close(1000)
         endif

         write(99,*) 'methode',methodeit
         write(99,*) 'Tolerance asked for the iterative method',tolinit
         write(99,*) 'Tolerance obtained for the iterative method',tol1
         write(99,*) 'Number of product Ax for the iterative method'
     $        ,ncompte,nloop

         if (tol1.ge.tolinit) then
            nstop=1
            infostr='Converge do not reach'
            write(*,*) tol1,tolinit
            write(99,*) tol1,tolinit
            return
         endif

         if (nstop == -1) then
            infostr = 'Calculation cancelled after iterative method'
            return
         endif
         write(*,*) '************* END LINEAR SYSTEM *****************'
         write(*,*) ' '
      elseif (nrig.eq.1) then
         write(*,*) '*************************************************'      
         write(*,*) '*** BEGIN RENOMALIZED BORN APPROXIMATION ********'
         write(*,*) '*************************************************'
         
c     Born approximation field
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)          
         do i=1,nbsphere3
            FFloc(i)=FF0(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL           

         write(*,*) '** END  RENOMALIZED BORN APPROXIMATION *********'
         write(*,*) ' '

      elseif (nrig.eq.2) then 
         write(*,*) '*************************************************'      
         write(*,*) '************* BEGIN BORN APPROXIMATION **********'
         write(*,*) '*************************************************'
         nsens=-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,kk,Em,Eloc,ii,jj,epsani,eps0,z)   
!$OMP DO SCHEDULE(STATIC) 
         do k=1,nbsphere
            z=zs(k)
            kk=3*(k-1)
            Em(1)=FF0(kk+1)
            Em(2)=FF0(kk+2)
            Em(3)=FF0(kk+3)
            do ii=1,3
               do jj=1,3
                  epsani(ii,jj)=epsilon(k,ii,jj)
               enddo
            enddo
            eps0=epscouche(numerocouche(z,neps,nepsmax,zcouche))
            call local_macro_surf(Eloc,Em,epsani,eps0,aretecube,k0
     $           ,nsens)
            FFloc(kk+1)=Eloc(1)
            FFloc(kk+2)=Eloc(2)
            FFloc(kk+3)=Eloc(3)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
         
         write(*,*) '************** END BORN APPROXIMATION ***********'
         write(*,*) ' '
      elseif (nrig.eq.3) then
         
         write(*,*) '*************************************************'      
         write(*,*) '************** BEGIN BORN ORDER 1 ***************'
         write(*,*) '*************************************************'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=FF0(i)           
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)  
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo        
!$OMP ENDDO 
!$OMP END PARALLEL

         call produitfftmatvectsurplus(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax,matindplan
     $        ,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12,a13,a22,a23,a31
     $        ,a32 ,a33,planb,planf)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)       
         do i=1,nbsphere
            k=3*(i-1)
            FFloc(k+1)=xr(k+1)
            FFloc(k+2)=xr(k+2)
            FFloc(k+3)=xr(k+3)       
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL      
         
      endif
c     dipole a partir champ local
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC)         
      do i=1,nbsphere
         k=3*(i-1)
         FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1,2)*FFloc(k
     $        +2)+polarisa(i,1,3)*FFloc(k+3)
         FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2,2)*FFloc(k
     $        +2)+polarisa(i,2,3)*FFloc(k+3)
         FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3,2)*FFloc(k
     $        +2)+polarisa(i,3,3)*FFloc(k+3)            
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
      
    

      write(99,*) 'End of the iterative method'
c     ******************************************************
c     compute the near field with FFT
c     ******************************************************

 1000 if (nprochefft.ge.1) then

         write(*,*) '*************************************************'      
         write(*,*) '************* STUDY LARGE NEAR  FIELD ***********'
         write(*,*) '*************************************************'
         
c     compte the size and edge of the box
         xmax=-1.d300
         xmin=1.d300
         ymax=-1.d300
         ymin=1.d300
         zmax=-1.d300
         zmin=1.d300      
!$OMP PARALLEL 
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:xmax,ymax,zmax)
!$OMP& REDUCTION(min:xmin,ymin,zmin)         
         do i=1,nbsphere
            xmax=max(xmax,xs(i))
            xmin=min(xmin,xs(i))
            ymax=max(ymax,ys(i))
            ymin=min(ymin,ys(i))
            zmax=max(zmax,zs(i))
            zmin=min(zmin,zs(i))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
         write(*,*) 'xmin     =',xmin
         write(*,*) 'ymin     =',ymin
         write(*,*) 'zmin     =',zmin
         write(*,*) 'meshsize = ',aretecube
c     compute new position of computation, dipole and incident field at
c     the new grid
       
         nxmpp=nx+2*nxmp
         nympp=ny+2*nymp
         nzmpp=nz+2*nzmp
         nmaxpp=nxmpp*nympp*nzmpp
         write(*,*) 'nx =',nx,'additional',nxmp,'total',nxmpp
         write(*,*) 'ny =',ny,'additional',nymp,'total',nympp
         write(*,*) 'nz =',nz,'additional',nzmp,'total',nzmpp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
         do i=1,3*nmaxpp
            xi(i)=0.d0
            xr(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL        

         
c     initialize 
         cntwf = 0
         l=1
         subunit=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
         do i=1,nmaxpp            
            polarisa(i,1,1)=1.d0
            polarisa(i,2,2)=1.d0
            polarisa(i,3,3)=1.d0
            polarisa(i,1,2)=0.d0
            polarisa(i,1,3)=0.d0
            polarisa(i,2,1)=0.d0
            polarisa(i,2,3)=0.d0
            polarisa(i,3,1)=0.d0
            polarisa(i,3,2)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

c     verification que aucun z n'est sur l'interface
         do i=1,nzmpp
            z=zmin+dble(i-nzmp-1)*aretecube
            do l=0,neps
               call comparaisonreel(z,zcouche(l),test)
               if (test.eq.0) then
                  nstop=1
                  infostr='Dipole on a interface for wide field'
               endif
            enddo
         enddo

c     write xyz wf for hdf5
         if (nmatf.eq.2) then
            do i=1,nzmpp
               zswf(i)=zmin+dble(i-nzmp-1)*aretecube
            enddo
            do j=1,nympp
               yswf(j)=ymin+dble(j-nymp-1)*aretecube
            enddo
            do k=1,nxmpp
               xswf(k)=xmin+dble(k-nxmp-1)*aretecube
            enddo
            dim(1)=nzmpp
            dim(2)=nmax
            datasetname='zwf'
            call hdf5write1d(group_idnf,datasetname,zswf,dim)
            dim(1)=nympp
            datasetname='ywf'
            call hdf5write1d(group_idnf,datasetname,yswf,dim)
            dim(1)=nxmpp
            datasetname='xwf'
            call hdf5write1d(group_idnf,datasetname,xswf,dim)
            
         endif
c     create the new vector position and memorize teh local field
         do i=1,nzmpp
            do j=1,nympp
               do k=1,nxmpp
                  cntwf = cntwf + 1
                  x=xmin+dble(k-nxmp-1)*aretecube
                  y=ymin+dble(j-nymp-1)*aretecube
                  z=zmin+dble(i-nzmp-1)*aretecube
                  
                  zswf(cntwf) = z
                  yswf(cntwf) = y
                  xswf(cntwf) = x
                  if (k.eq.1.and.j.eq.1.and.nmatf.eq.0) write(69,*) z
                  if (i.eq.1.and.k.eq.1.and.nmatf.eq.0) write(68,*) y
                  if (i.eq.1.and.j.eq.1.and.nmatf.eq.0) write(67,*) x
c     if (j.eq.1.and.k.eq.1) then
c     write(69,*) z
c     endif
c     if (i.eq.1.and.k.eq.1) then
c     write(68,*) y
c     endif
c     if (j.eq.1.and.i.eq.1) then
c     write(67,*) x
c     endif
                  subunit=subunit+1
                  Tabzn(subunit)=i
                  nsubunit=3*(subunit-1)
                  if (comparaison(x,y,z,xs(l),ys(l),zs(l)
     $                 ,lambda).eq.1.and.l.le.nbsphere)  then
                     ll=3*(l-1)        
c     memorize the dipole
                     xi(nsubunit+1)=FF(ll+1)
                     xi(nsubunit+2)=FF(ll+2)
                     xi(nsubunit+3)=FF(ll+3)

c     memorize the incident field
                     xr(nsubunit+1)=FF0(ll+1)
                     xr(nsubunit+2)=FF0(ll+2)
                     xr(nsubunit+3)=FF0(ll+3)
c     rewrite the epsilon
                     do ii=1,3
                        do jj=1,3
                           polarisa(subunit,ii,jj)=epsilon(l,ii,jj)
                        enddo
                     enddo
                     l=l+1
                  else
c     reuse pola; put the value of the relative permittivity of the layer
                     eps0=epscouche(numerocouche(z,neps,nepsmax
     $                    ,zcouche))
                     do ii=1,3                       
                        polarisa(subunit,ii,ii)=eps0
                     enddo
                     
c     comptute the incident field if is not computed yet
                     if (beam(1:11).eq.'pwavelinear') then
                        call champlineaire(epscouche,zcouche,neps
     $                       ,nepsmax,x,y,z,k0,E0,ss,pp,theta,phi,nstop
     $                       ,infostr,xr(nsubunit+1),xr(nsubunit+2)
     $                       ,xr(nsubunit +3))
                     elseif (beam(1:13).eq.'pwavecircular') then
                        call champcirculaire(epscouche,zcouche,neps
     $                       ,nepsmax,x,y,z,k0,E0,ss,theta,phi,infostr
     $                       ,nstop ,xr(nsubunit+1),xr(nsubunit+2)
     $                       ,xr(nsubunit +3))
                     elseif (beam(1:15).eq.'wavelinearmulti') then
                        call champlineairemulti(epscouche,zcouche,neps
     $                       ,nepsmax,x,y,z,k0,E0m,ssm,ppm,thetam,phim
     $                       ,nbinc,infostr,nstop,xr(nsubunit+1)
     $                       ,xr(nsubunit+2) ,xr(nsubunit +3))

                     elseif (beam(1:13).eq.'arbitrary') then
                        call incidentarbitrarypos(x,y,z,aretecube
     $                       ,xr(nsubunit+1),xr(nsubunit+2),xr(nsubunit
     $                       +3),nstop,namefileinc,infostr)
                     endif
                   
                  endif
                
               enddo
            enddo
         enddo
c     Far field discretization
         close(67)
         close(68)
         close(69)
         write(*,*) 'sub',subunit,nmaxpp
         
c     compute the incident field in the larger box in the case of the
c     gaussian beam as all the points are computed in one shot.

         if  (beam(1:11).eq.'gwavelinear') then
            call gaussiansurf(epscouche,zcouche,neps,nepsmax,k0,w0,
     $           xgaus,ygaus,zgaus,theta,phi,psi,E0,xswf,yswf,zswf,xr
     $           ,aretecube,nxmpp,nympp,nzmpp,subunit,nmax,nfft2d,nfft2d
     $           ,Eimagexneg,Eimageyneg ,Eimagezneg,Eimagexpos
     $           ,Eimageypos,Eimagezpos,fluxinc ,fluxref ,fluxtrans,irra
     $           ,nstop,infostr,plan2b)
c            write(*,*) 'grosse boite',fluxinc ,fluxref ,fluxtrans,irra
         elseif  (beam(1:13).eq.'gwavecircular') then
            call gaussiansurfcirc(epscouche,zcouche,neps,nepsmax,k0,w0
     $           ,xgaus ,ygaus,zgaus,theta,phi,ss,E0,xswf,yswf,zswf,xr
     $           ,aretecube ,nxmpp,nympp,nzmpp,subunit,nmax,nfft2d
     $           ,nfft2d ,Eimagexneg ,Eimageyneg ,Eimagezneg,Eimagexpos
     $           ,Eimageypos ,Eimagezpos ,fluxinc ,fluxref ,fluxtrans
     $           ,irra,nstop ,infostr,plan2b)
c            write(*,*) 'grosse boite',fluxinc ,fluxref ,fluxtrans,irra
         elseif  (beam(1:7).eq.'speckle') then
            call specklesurf(epscouche,zcouche,neps,nepsmax,k0,E0
     $           ,numaperref,IR,xswf,yswf,zswf,xgaus ,ygaus,zgaus,psi,xr
     $           ,aretecube,nxmpp,nympp,nzmpp,nxm,nym,nzm,subunit,nmax
     $           ,nfft2d ,nfft2d,Eimagexneg ,Eimageyneg ,Eimagezneg
     $           ,Eimagexpos ,Eimageypos ,Eimagezpos,fluxinc ,fluxref
     $           ,fluxtrans ,irra,nstop ,infostr,plan2b)
         endif
         
         
c     save the incident field for wide field
         if (nmatf.eq.0) then
            do i=1,subunit
               incidentfieldx(i) = xr(3*i-2)
               incidentfieldy(i) = xr(3*i-1)
               incidentfieldz(i) = xr(3*i)
               incidentfield(i) = dsqrt(cdabs(xr(3*i-2))**2
     $              +cdabs(xr(3*i-1))**2+cdabs(xr(3*i))**2)
               write(136,*) dreal(xr(3*i-2)),dimag(xr(3*i-2))
               write(137,*) dreal(xr(3*i-1)),dimag(xr(3*i-1))
               write(138,*) dreal(xr(3*i)),dimag(xr(3*i))
               write(139,*)  dsqrt(cdabs(xr(3*i-2))**2+cdabs(xr(3*i-1))
     $              **2+cdabs(xr(3*i))**2)
            enddo
         else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)            
            do i=1,subunit
               incidentfieldx(i) = xr(3*i-2)
               incidentfieldy(i) = xr(3*i-1)
               incidentfieldz(i) = xr(3*i)
               incidentfield(i) = dsqrt(dreal(xr(3*i-2)*dconjg(xr(3*i
     $              -2))+xr(3*i-1)*dconjg(xr(3*i-1))+xr(3*i)*dconjg(xr(3
     $              *i))))          
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL
            if (nmatf.eq.2) then
               dim(1)=subunit
               dim(2)=nmax
               datasetname='Incident field modulus wf'
               call hdf5write1d(group_idnf,datasetname,incidentfield,
     $              dim)
               datasetname='Incident field x component real part wf'
               call hdf5write1d(group_idnf,datasetname
     $              ,dreal(incidentfieldx),dim)
               datasetname
     $              ='Incident field x component imaginary part wf'
               call hdf5write1d(group_idnf,datasetname
     $              ,dimag(incidentfieldx),dim)
               datasetname='Incident field y component real part wf'
               call hdf5write1d(group_idnf,datasetname
     $              ,dreal(incidentfieldy),dim)
               datasetname
     $              ='Incident field y component imaginary part wf'
               call hdf5write1d(group_idnf,datasetname
     $              ,dimag(incidentfieldy),dim)
               datasetname='Incident field z component real part wf'
               call hdf5write1d(group_idnf,datasetname
     $              ,dreal(incidentfieldz),dim)
               datasetname
     $              ='Incident field z component imaginary part wf'
               call hdf5write1d(group_idnf,datasetname
     $              ,dimag(incidentfieldz),dim)
            endif
            
         endif
c     close Intensity of the incident field wide field
         close(136)
         close(137)
         close(138)
         close(139)
         write(*,*) 'Initialize plan for FFT for large field region'

         nxm2=nxmpp*2
         nym2=nympp*2
         
         call dfftw_plan_dft_2d(planbn, nxm2,nym2,b11,b11,FFTW_BACKWARD
     $        ,FFTW_ESTIMATE)
         call dfftw_plan_dft_2d(planfn, nxm2,nym2,b11,b11,FFTW_FORWARD
     $        ,FFTW_ESTIMATE)


c     compute Green function (local field)
         hcc=0.3d0
         epsabs=0.d0
         nt=1
         if (ninterp.eq.0) then

            
c     write(*,*) 'fonction interp2',hcc,tolinit,epsabs,aretecube
c     $           ,k0,neps,nepsmax,dcouche,zcouche,epscouche,subunit,nmax
c     $           ,n1m,nzm,nzm
            
c     write(*,*) 'retour',hcc,tolinit,epsabs,aretecube,k0 ,neps
c     $           ,nepsmax,dcouche,zcouche,epscouche,subunit,nmax ,n1m
c     $           ,nzm,nbs,nmat,nmatim,Tabzn
            call fonctiongreensurfcomp(hcc,tolinit,epsabs,xswf,yswf,zswf
     $           ,aretecube,k0,neps,nepsmax,dcouche,zcouche,epscouche
     $           ,subunit,nmax,n1m,nzm,nzmpp,nbs,nmat,nmatim,nplanm
     $           ,Tabzn ,a ,matind ,matindplan,matindice,matrange,nt)

c     compte FFT of the Green function
            call fonctiongreensurffft(nxmpp,nympp,nzmpp,nxm2,nym2,nxm
     $           ,nym,n1m,nzm,nplanm,nmatim,nbs,ntotalm,aretecube,a
     $           ,matind,matindplan,matindice,matrange,b11,b12,b13,b22
     $           ,b23,b31,b32,b33,a11,a12,a13,a22,a23,a31,a32,a33
     $           ,planbn)
         else
            
c     write(*,*) 'fonction interp2',hcc,tolinit,epsabs,nxm,nym
c     $           ,nzm,aretecube,k0,neps,nepsmax,dcouche,zcouche
c     $           ,epscouche,subunit,nmax,n1m,nzm,nzm,nbs,nmat,nmatim
            call  fonctiongreensurfcompinterp(hcc,tolinit,epsabs,nxmpp
     $           ,nympp,nzmpp,zswf,aretecube,k0,neps,nepsmax,dcouche
     $           ,zcouche,epscouche,subunit,nmax,n1m,nzm,nzmpp,nbs,nmat
     $           ,nmatim,nplanm,Tabzn,a,matind,matindplan,matindice
     $           ,matrange,ninter ,ninterp,nt)

            call fonctiongreensurfinterpfft(nxmpp,nympp,nzmpp,nxm2,nym2
     $           ,nxm,nym,n1m,nzm,nplanm,nmatim,nbs,ntotalm,ninter
     $           ,ninterp,aretecube,a,matind,matindplan,matindice
     $           ,matrange,b11,b12,b13,b22,b23,b31,b32,b33,a11,a12,a13
     $           ,a22,a23 ,a31,a32,a33,planbn)
         endif
c     produit ici de la matrice avec le vecteur
c     write(*,*) 'FF fait',FF(1),FF(2),FF(3)

         call produitfftmatvectsurplusboite(xi,xr,subunit,subunit,nxmpp
     $        ,nympp,nzmpp,nxm2,nym2,nxm,nym,nzm,nzm,nplanm,ntotalm
     $        ,nmax ,matindplan,b31,b32,b33,b11,b12,b13,a11,a12,a13 ,a22
     $        ,a23 ,a31 ,a32 ,a33,planbn,planfn)
        
c     remet la valeur du champ calculé précédemment dans l'objet quand
c     on n'est pas en rigoureux.
         if (nrig.ne.0) then
            subunit=0
            l=1
            do i=1,nzmpp
               do j=1,nympp
                  do k=1,nxmpp
                     x=xmin+dble(k-nxmp-1)*aretecube
                     y=ymin+dble(j-nymp-1)*aretecube
                     z=zmin+dble(i-nzmp-1)*aretecube
                     subunit=subunit+1              
                     nsubunit=3*(subunit-1)
                     if (comparaison(x,y,z,xs(l),ys(l),zs(l)
     $                    ,lambda).eq.1.and.l.le.nbsphere)  then
                        ll=3*(l-1)
                        eps0=epscouche(numerocouche(z,neps,nepsmax
     $                       ,zcouche))
                        
                        tmp=cdabs(epsilon(l,1,1)-eps0) +cdabs(epsilon(l
     $                       ,2,2)-eps0)+cdabs(epsilon(l ,3,3)-eps0)
     $                       +cdabs(epsilon(l,1,2)) +cdabs(epsilon(l,1
     $                       ,3))+cdabs(epsilon(l,2 ,1))+cdabs(epsilon(l
     $                       ,2,3))+cdabs(epsilon(l ,3,1))
     $                       +cdabs(epsilon(l,3,2))

                        if (tmp.ge.1.d-8) then
                           xr(nsubunit+1)=FFloc(ll+1)
                           xr(nsubunit+2)=FFloc(ll+2)
                           xr(nsubunit+3)=FFloc(ll+3)
                        endif
                        l=l+1
                     endif
                  enddo
               enddo
            enddo
            
         endif


         if (nlocal.eq.1) then
            write(*,*) 'compute local field large field'
            if (nmatf.eq.0) then
               do i=1,subunit
                  ii=3*(i-1)
                  localfieldx(i) = xr(ii+1)
                  localfieldy(i) = xr(ii+2)
                  localfieldz(i) = xr(ii+3)
                  localfield(i) = dsqrt(cdabs(xr(ii+1))**2
     $                 +cdabs(xr(ii+2))**2+cdabs(xr(ii+3))**2)
                  write(140,*) dreal(xr(ii+1)),dimag(xr(ii+1))
                  write(141,*) dreal(xr(ii+2)),dimag(xr(ii+2))
                  write(142,*) dreal(xr(ii+3)),dimag(xr(ii+3))
                  write(143,*) dsqrt(cdabs(xr(ii+1))**2+cdabs(xr(ii+2))
     $                 **2+cdabs(xr(ii+3))**2)
               enddo
            else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii)   
!$OMP DO SCHEDULE(STATIC)              
               do i=1,subunit
                  ii=3*(i-1)
                  localfieldx(i) = xr(ii+1)
                  localfieldy(i) = xr(ii+2)
                  localfieldz(i) = xr(ii+3)
                  localfield(i) = dsqrt(cdabs(xr(ii+1))**2 +cdabs(xr(ii
     $                 +2))**2+cdabs(xr(ii+3))**2)                  
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL              

               if (nmatf.eq.2) then
                 dim(1)=subunit
                 dim(2)=nmax
                 datasetname='Local field modulus wf'
                 call hdf5write1d(group_idnf,datasetname,localfield,dim)
                 datasetname='Local field x component real part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dreal(localfieldx),dim)
                 datasetname='Local field x component imaginary part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dimag(localfieldx),dim)
                 datasetname='Local field y component real part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dreal(localfieldy),dim)
                 datasetname='Local field y component imaginary part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dimag(localfieldy),dim)
                 datasetname='Local field z component real part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dreal(localfieldz),dim)
                 datasetname='Local field z component imaginary part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dimag(localfieldz),dim)
               endif
            endif
c     Intensity of the local field wide field
            close(140)
            close(141)
            close(142)
            close(143)            
         endif
         if (nmacro.eq.1) then
            write(*,*) 'compute macroscopic field'
            nsens=1
            if (nmatf.eq.0) then
               do i=1,subunit              
                  ii=3*(i-1)
                  Eloc(1)= xr(ii+1)
                  Eloc(2)= xr(ii+2)
                  Eloc(3)= xr(ii+3)


c     pour l'instant fait au niveau du epsilon pour le champ macro
                  do ii=1,3
                     do jj=1,3
                        epsani(ii,jj)=polarisa(i,ii,jj)
                     enddo
                  enddo 
                  eps0=epscouche(numerocouche(zswf(i),neps,nepsmax
     $                 ,zcouche))
                  call local_macro_surf(Eloc,Em,epsani,eps0,aretecube,k0
     $                 ,nsens)
c     write(*,*) 'local',eps0,zs(i),Eloc,Em,'ani',epsani
                  macroscopicfieldx(i) = Em(1)
                  macroscopicfieldy(i) = Em(2)
                  macroscopicfieldz(i) = Em(3)
                  macroscopicfield(i) = dsqrt(cdabs(Em(1))**2
     $                 +cdabs(Em(2))**2+cdabs(Em(3))**2)
c     write(*,*) 'macro',Em(1),Em(2),Em(3)
                  write(144,*) dreal(Em(1)),dimag(Em(1))
                  write(145,*) dreal(Em(2)),dimag(Em(2))
                  write(146,*) dreal(Em(3)),dimag(Em(3)) 
                  write(147,*) dsqrt(cdabs(Em(1))**2+cdabs(Em(2))**2
     $                 +cdabs(Em(3))**2)
               enddo
            else
               nsens=1
               
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii,jj,Eloc,epsani,Em,eps0)   
!$OMP DO SCHEDULE(STATIC)               
               do i=1,subunit              
                  ii=3*(i-1)
                  Eloc(1)= xr(ii+1)
                  Eloc(2)= xr(ii+2)
                  Eloc(3)= xr(ii+3)
c     pour l'instant faut au niveau du epsilonpour le champ macro
                  do ii=1,3
                     do jj=1,3
                        epsani(ii,jj)=polarisa(i,ii,jj)
                     enddo
                  enddo
                  eps0=epscouche(numerocouche(zswf(i),neps,nepsmax
     $                 ,zcouche))
                  call local_macro_surf(Eloc,Em,epsani,eps0,aretecube,k0
     $                 ,nsens)

                  macroscopicfieldx(i) = Em(1)
                  macroscopicfieldy(i) = Em(2)
                  macroscopicfieldz(i) = Em(3)
                  macroscopicfield(i) = dsqrt(dreal(Em(1)*dconjg(Em(1))
     $                 +Em(2)*dconjg(Em(2))+Em(3)*dconjg(Em(3))))
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
               if (nmatf.eq.2) then
                 dim(1)=subunit
                 dim(2)=nmax
                 datasetname='Macroscopic field modulus wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,macroscopicfield,dim)
                 datasetname
     $                ='Macroscopic field x component real part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dreal(macroscopicfieldx),dim)
                 datasetname
     $                ='Macroscopic field x component imaginary part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dimag(macroscopicfieldx),dim)
                 datasetname
     $                ='Macroscopic field y component real part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dreal(macroscopicfieldy),dim)
                 datasetname
     $                ='Macroscopic field y component imaginary part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dimag(macroscopicfieldy),dim)
                 datasetname
     $                ='Macroscopic field z component real part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dreal(macroscopicfieldz),dim)
                 datasetname
     $                ='Macroscopic field z component imaginary part wf'
                 call hdf5write1d(group_idnf,datasetname
     $                ,dimag(macroscopicfieldz),dim)
               endif
            endif            
c     close Intensity of the macroscopic field wide field
            close(144)
            close(145)
            close(146)
            close(147)
         endif
         write(*,*) '*************** END LARGE NEAR  FIELD ***********'
         write(*,*) ' '
      elseif (nlocal+nmacro.ge.1) then
         write(*,*) '*************************************************'      
         write(*,*) '************* COMPUTE NEAR  FIELD ***************'
         write(*,*) '*************************************************'
         subunit=0
         if (nlocal.eq.1) then
            write(*,*) 'compute local field'
            if (nmatf.eq.0) then
               do i=1,ndipole
                  k=tabdip(i)
                  if (k.ne.0) then
                     ii=3*(k-1)
                     localfieldx(k)=FFloc(ii+1)
                     localfieldy(k)=FFloc(ii+2)
                     localfieldz(k)=FFloc(ii+3)
                     localfield(k) = dsqrt(dreal(FFloc(ii+1)
     $                    *dconjg(FFloc(ii+1))+FFloc(ii+2)
     $                    *dconjg(FFloc(ii +2))+FFloc(ii+3)
     $                    *dconjg(FFloc(ii+3))))
                     write(40,*) dreal(FFloc(ii+1)),dimag(FFloc(ii+1))
                     write(41,*) dreal(FFloc(ii+2)),dimag(FFloc(ii+2))
                     write(42,*) dreal(FFloc(ii+3)),dimag(FFloc(ii+3))             
                     write(43,*) localfield(k)
                  else
                     write(40,*) 0.d0,0.d0
                     write(41,*) 0.d0,0.d0
                     write(42,*) 0.d0,0.d0
                     write(43,*) 0.d0
                  endif
               enddo
            else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,ii)   
!$OMP DO SCHEDULE(STATIC)
               do i=1,ndipole
                  k=tabdip(i)
                  if (k.ne.0) then
                     ii=3*(k-1)
                     localfieldx(k)=FFloc(ii+1)
                     localfieldy(k)=FFloc(ii+2)
                     localfieldz(k)=FFloc(ii+3)
                     localfield(k) = dsqrt(dreal(FFloc(ii+1)
     $                    *dconjg(FFloc(ii+1))+FFloc(ii+2)
     $                    *dconjg(FFloc(ii +2))+FFloc(ii+3)
     $                    *dconjg(FFloc(ii+3))))
                  endif
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL                 

               if (nmatf.eq.2) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC)
                  do i=1,ndipole
                     k=tabdip(i)
                     if (k.ne.0) then
                        wrk(i,1)=localfieldx(k)
                        wrk(i,2)=localfieldy(k)
                        wrk(i,3)=localfieldz(k)
                        wrk(i,4)= dsqrt(dreal(localfieldx(k)
     $                       *dconjg(localfieldx(k))+localfieldy(k)
     $                       *dconjg(localfieldy(k))+localfieldz(k)
     $                       *dconjg(localfieldz(k))))
                     else
                        wrk(i,1)=0.d0
                        wrk(i,2)=0.d0
                        wrk(i,3)=0.d0
                        wrk(i,4)=0.d0 
                     endif
                  enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
                  dim(1)=ndipole
                  dim(2)=nmax*3
                  datasetname='Local field modulus'
                  call hdf5write1d(group_idnf,datasetname,dreal(wrk(:,4)
     $                 ),dim)
                  datasetname='Local field x component real part'
                  call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,1)
     $                 ),dim)
                  datasetname='Local field x component imaginary part'
                  call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,1)
     $                 ),dim)
                  datasetname='Local field y component real part'
                  call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,2)
     $                 ),dim)
                  datasetname='Local field y component imaginary part'
                  call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,2)
     $                 ),dim)
                  datasetname='Local field z component real part'
                  call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,3)
     $                 ),dim)
                  datasetname='Local field z component imaginary part'
                  call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,3)
     $                 ),dim)
               endif

            endif        
c     close Intensity of the local field
            close(40)
            close(41)
            close(42)
            close(43)
         endif
c     compute and save the macroscopic field
         subunit=0
         if (nmacro.eq.1) then
            write(*,*) 'compute macrosocpic field'
            nsens=1
            if (nmatf.eq.0) then
               do i=1,ndipole
                  k=tabdip(i)
                  if (k.ne.0) then
                     ii=3*(k-1)
                     Eloc(1)= FFloc(ii+1)
                     Eloc(2)= FFloc(ii+2)
                     Eloc(3)= FFloc(ii+3)
                     do ii=1,3
                        do jj=1,3
                           epsani(ii,jj)=epsilon(k,ii,jj)
                        enddo
                     enddo 
                     eps0=epscouche(numerocouche(zs(k),neps,nepsmax
     $                    ,zcouche))
                     call local_macro_surf(Eloc,Em,epsani,eps0,aretecube
     $                    ,k0,nsens)
c     write(*,*) 'local',eps0,zs(i),Eloc,Em,'ani',epsani
                     macroscopicfieldx(k)=Em(1)
                     macroscopicfieldy(k)=Em(2)
                     macroscopicfieldz(k)=Em(3)
                     write(44,*) dreal(Em(1)),dimag(Em(1))
                     write(45,*) dreal(Em(2)),dimag(Em(2))
                     write(46,*) dreal(Em(3)),dimag(Em(3)) 
                     macroscopicfield(k)= dsqrt(cdabs(Em(1))**2+
     $                    cdabs(Em(2))**2+cdabs(Em(3))**2)
                     write(47,*) macroscopicfield(k)
                  else
                     write(44,*) 0.d0,0.d0
                     write(45,*) 0.d0,0.d0
                     write(46,*) 0.d0,0.d0
                     write(47,*) 0.d0
                  endif
               enddo
            else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,ii,jj,Eloc,epsani,Em,eps0)   
!$OMP DO SCHEDULE(STATIC)
               do i=1,ndipole
                  k=tabdip(i)
                  if (k.ne.0) then
                     ii=3*(k-1)
                     Eloc(1)= FFloc(ii+1)
                     Eloc(2)= FFloc(ii+2)
                     Eloc(3)= FFloc(ii+3)
                     do ii=1,3
                        do jj=1,3
                           epsani(ii,jj)=epsilon(k,ii,jj)
                        enddo
                     enddo
                     eps0=epscouche(numerocouche(zs(k),neps,nepsmax
     $                    ,zcouche))
                     call local_macro_surf(Eloc,Em,epsani,eps0
     $                    ,aretecube,k0,nsens)
                     macroscopicfieldx(k)=Em(1)
                     macroscopicfieldy(k)=Em(2)
                     macroscopicfieldz(k)=Em(3)
                     macroscopicfield(k)= dsqrt(dreal(Em(1)
     $                    *dconjg(Em(1))+Em(2)*dconjg(Em(2))+Em(3)
     $                    *dconjg(Em(3))))
                  endif
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL          


               if (nmatf.eq.2) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC)
                  do i=1,ndipole
                     k=tabdip(i)
                     if (k.ne.0) then
                        wrk(i,1)=macroscopicfieldx(k)
                        wrk(i,2)=macroscopicfieldy(k)
                        wrk(i,3)=macroscopicfieldz(k)
                        wrk(i,4)=dsqrt(dreal(macroscopicfieldx(k)
     $                       *dconjg(macroscopicfieldx(k))
     $                       +macroscopicfieldy(k)
     $                       *dconjg(macroscopicfieldy(k))
     $                       +macroscopicfieldz(k)
     $                       *dconjg(macroscopicfieldz(k))))
                     else
                        wrk(i,1)=0.d0
                        wrk(i,2)=0.d0
                        wrk(i,3)=0.d0
                        wrk(i,4)=0.d0
                     endif
                  enddo
!$OMP ENDDO 
!$OMP END PARALLEL    

                  dim(1)=ndipole
                  dim(2)=nmax*3
                  datasetname='Macroscopic field modulus'
                  call hdf5write1d(group_idnf,datasetname,dreal(wrk(:,4)
     $                 ),dim)
                  datasetname='Macroscopic field x component real part'
                  call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,1)
     $                 ),dim)
                  datasetname
     $                 ='Macroscopic field x component imaginary part'
                  call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,1)
     $                 ),dim)
                  datasetname='Macroscopic field y component real part'
                  call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,2)
     $                 ),dim)
                  datasetname
     $                 ='Macroscopic field y component imaginary part'
                  call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,2)
     $                 ),dim)
                  datasetname='Macroscopic field z component real part'
                  call hdf5write1d(group_idnf,datasetname,dreal(Wrk(:,3)
     $                 ),dim)
                  datasetname
     $                 ='Macroscopic field z component imaginary part'
                  call hdf5write1d(group_idnf,datasetname,dimag(Wrk(:,3)
     $                 ),dim)
               endif

            endif
c     close Intensity of the macroscopic field
            close(44)
            close(45)
            close(46)
            close(47)
            
         endif
         write(*,*) '***************** END NEAR  FIELD ***************'
         write(*,*) ' '
      endif

c     *******************************************************
c     End of the computation of the local and macrosocpic field
c     *******************************************************


      
c     partie sans sens pour la surface a virer plus tard
      if (nsection.eq.1) then
         write(*,*) '*************************************************'      
         write(*,*) '************* COMPUTE CROS SECTION **************'
         write(*,*) '*************************************************'
c     Compute the extinction cross section and absorbing cross section
         Cext=0.d0   
         Cabs=0.d0
         tmp=dsqrt(dreal(epscouche(0)))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,kk)   
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:Cext,Cabs)  
         do i=1,nbsphere
            kk=3*(i-1)
            Cext=Cext+dimag(dconjg(FF0(kk+1))*FF(kk+1)+dconjg(FF0(kk
     $           +2))*FF(kk+2)+dconjg(FF0(kk+3))*FF(kk+3))
            Cabs=Cabs+dimag((FF(kk+1)*dconjg(FFloc(kk+1))+FF(kk+2)
     $           *dconjg(FFloc(kk+2))+FF(kk+3)*dconjg(FFloc(kk+3))))
     $           -2.d0/3.d0*k03*tmp*tmp*tmp*dreal(FF(kk+1) *dconjg(FF(kk
     $           +1))+FF(kk+2)*dconjg(FF(kk+2))+FF(kk+3) *dconjg(FF(kk
     $           +3)))/tmp/tmp
c     write(*,*) 'rr',dimag((FF(kk+1)*dconjg(FFloc(kk+1))+FF(kk+2)
c     $           *dconjg(FFloc(kk+2))+FF(kk+3)*dconjg(FFloc(kk+3))))
c     $           ,1.d0/3.d0*k03*tmp*tmp*tmp*dreal(FF(kk+1) *dconjg(FF(kk
c     $           +1))+FF(kk+2)*dconjg(FF(kk+2))+FF(kk+3) *dconjg(FF(kk
c     $           +3))),kk
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  

c     Cext=4.d0*pi*k0*Cext/I0*dsqrt(dreal(epscouche(0)))
c     Cabs=4.d0*pi*k0*Cabs/I0*dsqrt(dreal(epscouche(0)))
c     Pdissapeted=W/2Im(p.E^*) et Cext=Pdissapeted/I avec I=0.5 c eps E^2
         Cext=4.d0*pi*k0*Cext/(8.d0*pi*1.d-7*c*P0)*pi*w0*w0
         Cabs=4.d0*pi*k0*Cabs/(8.d0*pi*1.d-7*c*P0)*pi*w0*w0

c     compute the scattering cross section
         Csca=Cext-Cabs
         write(*,*) 'Cext = ',Cext
         if (Cabs.le.Cext*1.d-12) then
            Cabs=0.d0
            write(*,*) 'Cabs = ',Cabs
         else
            write(*,*) 'Cabs = ',Cabs
         endif
         
         
         write(99,*) 'extinction cross section',Cext
         write(99,*) 'absorbing cross section ',Cabs
         if (Cabs.le.Cext*1.d-12) then
            Cabs=0.d0
            write(99,*) 'absorbing cross section ',Cabs
         else
            write(99,*) 'absorbing cross section ',Cabs
         endif
         if (nstop == -1) then
            infostr = 'Calculation cancelled after cross section'
            return
         endif
   

         if (object(1:6).eq.'sphere') then
            tmp=dsqrt(dreal(epscouche(0)))
            CALL CALLBHMIE(tmp,eps,rayon,lambda,MIECEXT ,MIECABS,MIECSCA
     $           ,GSCA)
            write(*,*) 'Comparison with Mie s series'
            write(*,*) 'Cext error in % :',(Cext-MIECEXT)/MIECEXT*100.d0
            write(*,*) 'Csca error in % :',(Csca-MIECSCA)/MIECSCA*100.d0
            if (Cabs.ne.0.d0) then
               write(*,*) 'Cabs error in % :',(Cabs-MIECABS) / MIECABS
     $              *100.d0
            endif
         endif
         write(*,*) '************* END CROSS SECTION *****************'
         write(*,*) ' '
      endif
     
     
c******************************************************************
c     prepare the computation of the far field in the FFT is choosen
c******************************************************************
      if (nenergie+nlentille+ndiffracte.ge.1.and.nquickdiffracte.eq.1)
     $     then
         write(*,*) '*************************************************'      
         write(*,*) '******* COMPUTE DIFFRACTED FIELD WITH FFT *******'
         write(*,*) '*************************************************'

      
         
         tmp=1.d0
         rloin=1.d0
         call diffractefft2dsurf2(nbsphere,nx,ny,nz,nxm,nym,nzm
     $        ,nfft2dtmp,nfft2d,k0,xs,ys,zs,aretecube,Efourierxpos
     $        ,Efourierypos,Efourierzpos,FF,imaxk0,deltakx,deltaky
     $        ,Ediffkzpos,Ediffkzneg,rloin,rloin,tmp,tmp,nepsmax ,neps
     $        ,dcouche,zcouche,epscouche,ncote ,nstop ,infostr,plan2f)

         write(*,*) '******* END DIFFRACTED FIELD WITH FFT ***********'
         write(*,*) ' '
      endif

      
c************************************************************
c     compute the the Poynting vector along the normal in CGS n.P=c/(8
c     pi) Re(ExH^*)=c/(8 pi)|E|^2
c************************************************************
      if (ndiffracte.eq.1) then
         
         write(*,*) '*************************************************'      
         write(*,*) '********** COMPUTE Csca g AND POYNTING ***********'
         write(*,*) '*************************************************'
         if (nsection.eq.0) then
            tmp=dsqrt(dreal(epscouche(0)))
            CALL CALLBHMIE(tmp,eps,rayon,lambda,MIECEXT ,MIECABS,MIECSCA
     $           ,GSCA)
         endif

            
         rloin=1.d0
         if (nquickdiffracte.eq.0) then
            write(*,*) 'Slow method with discretization:',ntheta,nphi,pi
c     compute the diffracted field in summing the dipole
            if (ntheta.le.10.or.nphi.le.20) then
               infostr='ntheta and nphi too small'
               nstop=1
               return
            endif

            thetad=0.d0
            deltatheta=pi/dble(ntheta)
            deltaphi=2.d0*pi/dble(nphi)
            Poyntinginc=0.d0
            cnt = 0
            Cscai=0.d0
            gasym=0.d0

            if (ncote.eq.-1) then
               deltatheta=pi/dble(ntheta)/2.d0
               thetad=pi/2.d0
            endif

            
c     write(*,*) 'coucou',ncote,numaper,deltatheta,deltaphi,thetad
c     write(*,*) 'thetaphi',ntheta,nphi
            z=rloin*dcos(pi/2.d0-deltatheta/100.d0)
c     write(*,*) 'zz',z,zcouche
            if (z.le.zcouche(neps+1)) then              
               rloin=max(rloin,z/dcos(pi/2.d0-deltatheta/100.d0)*2.d0)
c     write(*,*) 'change rloin1',rloin
            endif
            z=rloin*dcos(pi/2.d0+deltatheta/100.d0)
c     write(*,*) 'zz',z,zcouche
            if (z.ge.zcouche(0)) then              
               rloin=max(rloin,z/dcos(pi/2.d0+deltatheta/100.d0)*2.d0)
c     write(*,*) 'change rloin2',rloin
            endif
c     write(*,*)'rr',nbsphere,object,  deltatheta,   deltaphi    
            
            do itheta=0,ntheta
               thetas=thetad+deltatheta*dble(itheta)
c     evite le calcul sur le multicouche
               if (dabs(thetas-pi/2.d0).le.deltatheta/100.d0) then
                  thetas=pi/2.d0+sign(deltatheta/100.d0,thetas-pi/2.d0)
     $                 *deltatheta/100.d0
               endif
               do iphi=0,nphi-1
                  cnt = cnt + 1
                  phis=deltaphi*dble(iphi)
                  normal(1)=dsin(thetas)*dcos(phis)
                  normal(2)=dsin(thetas)*dsin(phis)
                  normal(3)=dcos(thetas)
                  x=rloin*normal(1)
                  y=rloin*normal(2)
                  z=rloin*normal(3)
c     write(*,*) 'obs',x,y,z,iphi,itheta
                  do ii=1,3
                     Em(ii)=0.d0
                  enddo
c     somme du champ rayonne par tous les dipoles.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,kk,ii,jj,Stenseur)   
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:Em) 
                  do i=1,nbsphere
                     kk=3*(i-1)
                     call tenseurmulticoucheloin(x,y,z,xs(i),ys(i),zs(i)
     $                    ,k0,nepsmax,neps,dcouche,zcouche,epscouche
     $                    ,Stenseur)
c     call propa_espace_libre(x,y,z,xs(i),ys(i),zs(i),k0,
c     $                    Stenseur)
c     write(*,*) 'coucounan',x,y,z,xs(i),ys(i),zs(i)
c     $                    ,Stenseur
                     do ii=1,3
                        do jj=1,3
                           Em(ii)=Em(ii)+Stenseur(ii,jj)*FF(kk+jj)
                        enddo
                     enddo
                  enddo   
!$OMP ENDDO 
!$OMP END PARALLEL 
                 
                  Emod=dreal(Em(1)*dconjg(Em(1))+Em(2)*dconjg(Em(2))
     $                 +Em(3)*dconjg(Em(3)))
c     write(*,*) 'Emod',Emod,iphi,itheta,c,pi,quatpieps0
c     $                 ,indicen,indice0
                  Emod=Emod*rloin*rloin
                  
                  if (thetas.le.pi/2.d0) then                    
                     poyntingfield(cnt)=Emod*c/8.d0/pi*quatpieps0
     $                    *indicen                    
                     thetafield(cnt) = thetas*180.d0/pi                    
                     phifield(cnt) = phis*180.d0/pi
                     if (nmatf.eq.0) then
                        write(50,*) Emod*c/8.d0/pi*quatpieps0 *indicen
                        write(51,*) thetas*180.d0/pi
                        write(52,*) phis*180.d0/pi
                     endif
                  else
                     poyntingfield(cnt)=Emod*c/8.d0/pi*quatpieps0
     $                    *indice0
                     thetafield(cnt) = thetas*180.d0/pi
                     phifield(cnt) = phis*180.d0/pi
                     if (nmatf.eq.0) then
                        write(50,*) Emod*c/8.d0/pi*quatpieps0 *indicen
                        write(51,*) thetas*180.d0/pi
                        write(52,*) phis*180.d0/pi
                     endif
                  endif
                  Cscai=Cscai+deltaphi*deltatheta*dsin(thetas)*Emod
                  gasym=gasym+deltaphi*deltatheta*dsin(thetas)*Emod
     $                 *(normal(1)*dsin(theta*pi/180.d0)*dcos(phi*pi
     $                 /180.d0)+normal(2)*dsin(theta*pi/180.d0)*
     $                 dsin(phi *pi/180.d0)+normal(3)*dcos(theta*pi
     $                 /180.d0))
               enddo
            enddo
            close(50)
            close(51)
            close(52)

            if (nmatf.eq.2) then
               dim(1)=cnt
               dim(2)=max((ntheta+1)*nphi,nfft2d*nfft2d)
               datasetname='Poynting'
               call hdf5write1d(group_idff,datasetname,poyntingfield
     $              ,dim)
               datasetname='theta for Poynting'
               call hdf5write1d(group_idff,datasetname,thetafield,dim)
               datasetname='phi for Poynting'
               call hdf5write1d(group_idff,datasetname,phifield,dim)
            endif
            
            gasym=gasym/cscai          
            Cscai=Cscai/I0
            
            if (beam(1:11).ne.'pwavelinear'.and.beam(1:13).ne
     $           .'pwavecircular'.or.ncote.ne.0.or.nhomo.ne.0) then
               gasym=0.d0
               Cscai=0.d0
            else
               write(*,*) 'Csca =',Cscai
               write(*,*) 'g    =',gasym
               if (object(1:6).eq.'sphere') then             
                  write(*,*) 'Comparison with Mie s series'
                  write(*,*) 'Csca error in %',(Cscai-MIECSCA)/MIECSCA
     $                 *100.d0
                  write(*,*) 'g error in %',(gasym-GSCA)/GSCA*100.d0
               endif
            endif
         else
c     compute the diffracted field with FFT
            write(*,*) 'Compute with FFT'
            if (nstop.eq.1) return
            
            call computegcfft2dsurf(imaxk0,deltakx,deltaky,k0,P0,theta
     $           ,phi,nfft2d,Ediffkzpos,Ediffkzneg,gasym,Cscai,w0
     $           ,nepsmax ,neps ,dcouche,zcouche ,epscouche,ncote)
            
c            write(*,*) 'Csca',Cscai,nstop,'NA',numaper,'deltakx'
c     $           ,deltakx,'imaxk0',imaxk0,ncote,beam,object

            if (beam(1:11).ne.'pwavelinear'.and.beam(1:13).ne
     $           .'pwavecircular'.or.ncote.ne.0.or.nhomo.ne.0) then
               gasym=0.d0
               Cscai=0.d0
            else
               write(*,*) 'Csca =',Cscai
               write(*,*) 'g    =',gasym
               if (object(1:6).eq.'sphere') then             
                  write(*,*) 'Comparison with Mie s series'
                  write(*,*) 'Csca error in %',(Cscai-MIECSCA)/MIECSCA
     $                 *100.d0
                  write(*,*) 'g error in %',(gasym-GSCA)/GSCA*100.d0
               endif
            endif

            
            if (nstop.eq.-1) return
            deltatheta=pi/dble(ntheta)
            deltaphi=2.d0*pi/dble(nphi)
c     put the field in memory with the right angles
            cnt=0
            do itheta=0,ntheta
               thetas=deltatheta*dble(itheta)
               do iphi=0,nphi-1
                  cnt = cnt + 1
                  phis=deltaphi*dble(iphi)
                  
                  thetafield(cnt) = thetas*180.d0/pi
                  phifield(cnt) = phis*180.d0/pi

                  if (nmatf.eq.0) then
                     write(52,*) phis*180.d0/pi
                     write(51,*) thetas*180.d0/pi
                  endif
                  
                  kx=k0*indicem*dsin(thetas)*dcos(phis)
                  ky=k0*indicem*dsin(thetas)*dsin(phis)
                  do i=-imaxk0,imaxk0
                     if (kx.ge.i*deltakx) i2=i
                  enddo
                  do j=-imaxk0,imaxk0
                     if (ky.ge.j*deltakx) j2=j
                  enddo
                  ii=imaxk0+i2+1
                  jj=imaxk0+j2+1
                
                  test=0
                  if (thetas.le.pi/2.d0) then
                     Emod11=cdabs(Ediffkzpos(ii,jj,1))**2
     $                    +cdabs(Ediffkzpos(ii,jj,2))**2
     $                    +cdabs(Ediffkzpos(ii,jj,3))**2
                     if (Emod11.eq.0.d0) test=test+1
                     Emod22=cdabs(Ediffkzpos(ii+1,jj+1,1))**2
     $                    +cdabs(Ediffkzpos(ii+1,jj+1,2))**2
     $                    +cdabs(Ediffkzpos(ii+1,jj+1,3))**2
                     if (Emod22.eq.0.d0)  test=test+1
                     Emod21=cdabs(Ediffkzpos(ii+1,jj,1))**2
     $                    +cdabs(Ediffkzpos(ii+1,jj,2))**2
     $                    +cdabs(Ediffkzpos(ii+1,jj,3))**2
                     if (Emod21.eq.0.d0) test=test+1
                     Emod12=cdabs(Ediffkzpos(ii,jj+1,1))**2
     $                    +cdabs(Ediffkzpos(ii,jj+1,2))**2
     $                    +cdabs(Ediffkzpos(ii,jj+1,3))**2
                     if (Emod12.eq.0.d0) test=test+1
                     if (test.eq.0) then
                        Emod=(Emod21-Emod11)*(kx-deltakx*i2)/deltakx
     $                       +(Emod12-Emod11)*(ky-deltakx*j2)/deltakx
     $                       +(Emod11+Emod22-Emod12-Emod21)*(kx-deltakx
     $                       *i2) /deltakx*(ky-deltakx*j2)/deltakx
     $                       +Emod11
                     else
                        Emod=(Emod11+Emod12+Emod21+Emod22)/dble(4-test)
                     endif
                     
                  else
                     Emod11=cdabs(Ediffkzneg(ii,jj,1))**2
     $                    +cdabs(Ediffkzneg(ii,jj,2))**2
     $                    +cdabs(Ediffkzneg(ii,jj,3))**2
                     if (Emod11.eq.0.d0) test=test+1
                     Emod22=cdabs(Ediffkzneg(ii+1,jj+1,1))**2
     $                    +cdabs(Ediffkzneg(ii+1,jj+1,2))**2
     $                    +cdabs(Ediffkzneg(ii+1,jj+1,3))**2
                     if (Emod22.eq.0.d0) test=test+1
                     Emod21=cdabs(Ediffkzneg(ii+1,jj,1))**2
     $                    +cdabs(Ediffkzneg(ii+1,jj,2))**2
     $                    +cdabs(Ediffkzneg(ii+1,jj,3))**2
                     if (Emod21.eq.0.d0)  test=test+1
                     Emod12=cdabs(Ediffkzneg(ii,jj+1,1))**2
     $                    +cdabs(Ediffkzneg(ii,jj+1,2))**2
     $                    +cdabs(Ediffkzneg(ii,jj+1,3))**2
                     if (Emod12.eq.0.d0) test=test+1
                     if (test.eq.0) then
                        Emod=(Emod21-Emod11)*(kx-deltakx*i2)/deltakx
     $                       +(Emod12-Emod11)*(ky-deltakx*j2)/deltakx
     $                       +(Emod11+Emod22-Emod12-Emod21)*(kx-deltakx
     $                       *i2) /deltakx*(ky-deltakx*j2)/deltakx
     $                       +Emod11
                     else
                        Emod=(Emod11+Emod12+Emod21+Emod22)/dble(4-test)
                     endif
                  endif      
                  poyntingfield(cnt) = Emod*k0*k0*k0*k0*c/8.d0/pi
     $                 *quatpieps0
                  if (nmatf.eq.0) write(50,*) Emod*k0*k0*k0*k0*c/8.d0
     $                 /pi*quatpieps0
                  
               enddo
            enddo
            close(50)
            close(51)
            close(52)
            if (nmatf.eq.2) then
               dim(1)=cnt
               dim(2)=max((ntheta+1)*nphi,nfft2d*nfft2d)
               datasetname='Poynting'
               call hdf5write1d(group_idff,datasetname,poyntingfield
     $              ,dim)
               datasetname='theta for Poynting'
               call hdf5write1d(group_idff,datasetname,thetafield,dim)
               datasetname='phi for Poynting'
               call hdf5write1d(group_idff,datasetname,phifield,dim)
            endif
c            write(*,*) 'imaxk000',imaxk0,ncote,deltakx,indicen,indice0
c     mise en memoire de kx et ky
            if (nmatf.eq.0) then
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0                  
                     kx=dble(i)*deltakx
                     ky=dble(j)*deltaky             
                     if (j.eq.0) write(520,*) kx
                     if (i.eq.0) write(521,*) ky
                  enddo
               enddo
c     ecriture en kx,ky du champ diffracte
               if (ncote.eq.0.or.ncote.eq.1) then
                  do i=-imaxk0,imaxk0
                     do j=-imaxk0,imaxk0              
                        kx=dble(i)*deltakx
                        ky=dble(j)*deltaky
                        ii=imaxk0+i+1
                        jj=imaxk0+j+1
                        if (k0*k0*indicen*indicen-kx*kx-ky*ky.gt.0.d0)
     $                       then 

c     write(*,*) 'dessus diff',ncote
                           write(500,*) dreal(Ediffkzpos(ii,jj,1))
     $                          ,dimag(Ediffkzpos(ii,jj,1))
                           write(501,*) dreal(Ediffkzpos(ii,jj,2))
     $                          ,dimag(Ediffkzpos(ii,jj,2))
                           write(502,*) dreal(Ediffkzpos(ii,jj,3))
     $                          ,dimag(Ediffkzpos(ii,jj,3))
                           Emod=cdabs(Ediffkzpos(ii,jj,1))**2
     $                          +cdabs(Ediffkzpos(ii,jj,2))**2
     $                          +cdabs(Ediffkzpos(ii,jj,3))**2
                           write(53,*) Emod*c/8/pi*quatpieps0
                           write(9000,*) Emod*c/8/pi*quatpieps0,i,j
     $                          ,dsqrt((kx*kx+ky*ky)/K0/k0)
                        else
                           write(500,*) 0.d0,0.d0
                           write(501,*) 0.d0,0.d0
                           write(502,*) 0.d0,0.d0
                           Emod=0.d0
                           write(53,*) Emod
                        endif
                     enddo
                  enddo
               endif
               if (ncote.eq.0.or.ncote.eq.-1) then
                  do i=-imaxk0,imaxk0
                     do j=-imaxk0,imaxk0
                        kx=dble(i)*deltakx
                        ky=dble(j)*deltaky
                        ii=imaxk0+i+1
                        jj=imaxk0+j+1
                        if (k0*k0*indice0*indice0-kx*kx-ky*ky.gt.0.d0)
     $                       then

                           write(506,*) dreal(Ediffkzneg(ii,jj,1))
     $                          ,dimag(Ediffkzneg(ii,jj,1))
                           write(507,*) dreal(Ediffkzneg(ii,jj,2))
     $                          ,dimag(Ediffkzneg(ii,jj,2))
                           write(508,*) dreal(Ediffkzneg(ii,jj,3))
     $                          ,dimag(Ediffkzneg(ii,jj,3))
                           Emod=cdabs(Ediffkzneg(ii,jj,1))**2
     $                          +cdabs(Ediffkzneg(ii,jj,2))**2
     $                          +cdabs(Ediffkzneg(ii,jj,3))**2
                           write(54,*) Emod*c/8/pi*quatpieps0
                        else
                           
                           write(506,*) 0.d0,0.d0
                           write(507,*) 0.d0,0.d0
                           write(508,*) 0.d0,0.d0
                           Emod=0.d0
                           write(54,*) Emod
                        endif
                     enddo
                  enddo
               endif
               close(500)
               close(501)
               close(502)
               close(506)
               close(507)
               close(508)
               close(53)
               close(54)
               close(520)
               close(521)
            elseif (nmatf.eq.2) then
               call writehdf5farfield(Ediffkzpos,Ediffkzneg,Eimagexpos
     $              ,Eimageypos,Eimagezpos,Eimageincxpos,kxy,deltakx
     $              ,deltaky,k0,imaxk0,nfft2d,indice0,indicen ,ncote
     $              ,name,group_idff)
            endif
         endif
         write(*,*) '************ END Csca g AND POYNTING ************'
         write(*,*) ' '
      endif

c*************************************************************
c     compute the far field with slow method for energy of lens
c*************************************************************
      if (nquickdiffracte.eq.0.and.nenergie+nlentille.ge.1) then

         write(*,*) '*************************************************'      
         write(*,*) '*** COMPUTE FAR FIELD WITH SLOW METHOD **********'
         write(*,*) '*************************************************'

         
c     initialization table
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)   
!$OMP DO SCHEDULE(STATIC)  COLLAPSE(2)
         do i=1,nfft2d
            do j=1,nfft2d
               Ediffkzpos(i,j,1)=0.d0
               Ediffkzpos(i,j,2)=0.d0
               Ediffkzpos(i,j,3)=0.d0
               Ediffkzneg(i,j,1)=0.d0
               Ediffkzneg(i,j,2)=0.d0
               Ediffkzneg(i,j,3)=0.d0
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  

         write(*,*) 'compute far field'
c     compute the diffracted field
         do i=-imaxk0,imaxk0               
            if (i.ge.0) then
               indicex=i+1
            else
               indicex=nfft2d+i+1
            endif
            kx=deltakx*dble(i)
            kxy(i+nfft2d2+1)=kx/k0
            xy(i+nfft2d2+1)=deltax*dble(i)*gross
            do j=-imaxk0,imaxk0
               
               if (j.ge.0) then
                  indicey=j+1
               else
                  indicey=nfft2d+j+1
               endif
               ky=deltaky*dble(j)
               ii=imaxk0+i+1
               jj=imaxk0+j+1

               if (ncote.eq.1.or.ncote.eq.0) then
c     calcul champ dessus
                  
                  if (k0*k0*indicen*indicen-kx*kx-ky*ky.gt.0.d0) then  
                     z=1.d0
                     kz=dsqrt(k0*k0*indicen*indicen-kx*kx-ky*ky)
                     call  tenseurmulticoucheloinfft(kx,ky,kz,z,zs(1)
     $                    ,k0,nepsmax,neps,dcouche,zcouche ,epscouche
     $                    ,Stenseur)
                     ctmp=cdexp(-icomp*(kx*xs(1)+ky*ys(1)))
                     Emx=(Stenseur(1,1)*FF(1)+Stenseur(1,2)*FF(2)
     $                    +Stenseur(1,3)*FF(3))*ctmp
                     Emy=(Stenseur(2,1)*FF(1)+Stenseur(2,2)*FF(2)
     $                    +Stenseur(2,3)*FF(3))*ctmp
                     Emz=(Stenseur(3,1)*FF(1)+Stenseur(3,2)*FF(2)
     $                    +Stenseur(3,3)*FF(3))*ctmp
c!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,kk,Stenseur,ctmp)   
c!$OMP DO SCHEDULE(STATIC) REDUCTION(+:Emx,Emy,Emz)      
                     do k=2,nbsphere
                        kk=3*(k-1)
                        if (zs(k).ne.zs(k-1)) then
                           call  tenseurmulticoucheloinfft(kx,ky,kz,z
     $                          ,zs(k),k0,nepsmax,neps,dcouche,zcouche
     $                          ,epscouche,Stenseur)                         
                        endif
                        ctmp=cdexp(-icomp*(kx*xs(k)+ky*ys(k)))
                        Emx=Emx+(Stenseur(1,1)*FF(kk+1)+Stenseur(1,2)
     $                       *FF(kk+2)+Stenseur(1,3)*FF(kk+3))*ctmp
                        Emy=Emy+(Stenseur(2,1)*FF(kk+1)+Stenseur(2,2)
     $                       *FF(kk+2)+Stenseur(2,3)*FF(kk+3))*ctmp
                        Emz=Emz+(Stenseur(3,1)*FF(kk+1)+Stenseur(3,2)
     $                       *FF(kk+2)+Stenseur(3,3)*FF(kk+3))*ctmp
                     enddo        
c!$OMP ENDDO 
c!$OMP END PARALLEL  
          
                     Ediffkzpos(ii,jj,1)=Emx
                     Ediffkzpos(ii,jj,2)=Emy
                     Ediffkzpos(ii,jj,3)=Emz
                  endif
               endif

c     calcul champ dessous
               if (ncote.eq.-1.or.ncote.eq.0) then
                  if (k0*k0*indice0*indice0-kx*kx-ky*ky.gt.0.d0) then    
                     z=-1.d0
                     kz=-dsqrt(k0*k0*indice0*indice0-kx*kx-ky*ky)
                     call  tenseurmulticoucheloinfft(kx,ky,kz,z,zs(1)
     $                    ,k0,nepsmax,neps,dcouche,zcouche ,epscouche
     $                    ,Stenseur)
                     ctmp=cdexp(-icomp*(kx*xs(1)+ky*ys(1)))
                     Emx=(Stenseur(1,1)*FF(1)+Stenseur(1,2)*FF(2)
     $                    +Stenseur(1,3)*FF(3))*ctmp
                     Emy=(Stenseur(2,1)*FF(1)+Stenseur(2,2)*FF(2)
     $                    +Stenseur(2,3)*FF(3))*ctmp
                     Emz=(Stenseur(3,1)*FF(1)+Stenseur(3,2)*FF(2)
     $                    +Stenseur(3,3)*FF(3))*ctmp
c!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,kk,ctmp,Stenseur)   
c!$OMP DO SCHEDULE(STATIC) REDUCTION(+:Emx,Emy,Emz)                         
                     do k=2,nbsphere
                        kk=3*(k-1)
                        if (zs(k).ne.zs(k-1)) then
                           call  tenseurmulticoucheloinfft(kx,ky,kz,z
     $                          ,zs(k),k0,nepsmax,neps,dcouche,zcouche
     $                          ,epscouche,Stenseur)
                        endif
                        ctmp=cdexp(-icomp*(kx*xs(k)+ky*ys(k)))
                        Emx=Emx+(Stenseur(1,1)*FF(kk+1)+Stenseur(1,2)
     $                       *FF(kk+2)+Stenseur(1,3)*FF(kk+3))*ctmp
                        Emy=Emy+(Stenseur(2,1)*FF(kk+1)+Stenseur(2,2)
     $                       *FF(kk+2)+Stenseur(2,3)*FF(kk+3))*ctmp
                        Emz=Emz+(Stenseur(3,1)*FF(kk+1)+Stenseur(3,2)
     $                       *FF(kk+2)+Stenseur(3,3)*FF(kk+3))*ctmp
                     enddo
c!$OMP ENDDO 
c!$OMP END PARALLEL            
                     Ediffkzneg(ii,jj,1)=Emx
                     Ediffkzneg(ii,jj,2)=Emy
                     Ediffkzneg(ii,jj,3)=Emz
                  endif
               endif

            enddo
         enddo


         if (nmatf.eq.0) then
            if (ncote.eq.0.or.ncote.eq.1) then
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0              
                     kx=dble(i)*deltakx
                     ky=dble(j)*deltaky
                     ii=imaxk0+i+1
                     jj=imaxk0+j+1
                     if (k0*k0*indicen*indicen-kx*kx-ky*ky.gt.0.d0) then
                        write(500,*) dreal(Ediffkzpos(ii,jj,1))
     $                       ,dimag(Ediffkzpos(ii,jj,1))
                        write(501,*) dreal(Ediffkzpos(ii,jj,2))
     $                       ,dimag(Ediffkzpos(ii,jj,2))
                        write(502,*) dreal(Ediffkzpos(ii,jj,3))
     $                       ,dimag(Ediffkzpos(ii,jj,3))
                        Emod=cdabs(Ediffkzpos(ii,jj,1))**2
     $                       +cdabs(Ediffkzpos(ii,jj,2))**2
     $                       +cdabs(Ediffkzpos(ii,jj,3))**2               
                     else
                        write(500,*) 0.d0,0.d0
                        write(501,*) 0.d0,0.d0
                        write(502,*) 0.d0,0.d0                     
                     endif
                  enddo
               enddo
            endif
            if (ncote.eq.0.or.ncote.eq.-1) then
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0              
                     kx=dble(i)*deltakx
                     ky=dble(j)*deltaky
                     ii=imaxk0+i+1
                     jj=imaxk0+j+1
                     if (k0*k0*indice0*indice0-kx*kx-ky*ky.gt.0.d0) then
                        write(506,*) dreal(Ediffkzneg(ii,jj,1))
     $                       ,dimag(Ediffkzneg(ii,jj,1))
                        write(507,*) dreal(Ediffkzneg(ii,jj,2))
     $                       ,dimag(Ediffkzneg(ii,jj,2))
                        write(508,*) dreal(Ediffkzneg(ii,jj,3))
     $                       ,dimag(Ediffkzneg(ii,jj,3))
                        Emod=cdabs(Ediffkzneg(ii,jj,1))**2
     $                       +cdabs(Ediffkzneg(ii,jj,2))**2
     $                       +cdabs(Ediffkzneg(ii,jj,3))**2
                     else
                        write(506,*) 0.d0,0.d0
                        write(507,*) 0.d0,0.d0
                        write(508,*) 0.d0,0.d0
                     endif
                  enddo
               enddo
            endif
         elseif (nmatf.eq.2) then
            call writehdf5farfield(Ediffkzpos,Ediffkzneg,Eimagexpos
     $           ,Eimageypos,Eimagezpos,Eimageincxpos,kxy,deltakx
     $           ,deltaky ,imaxk0,nfft2d,indice0,indicen ,ncote,name
     $           ,group_idff)
         endif
         write(*,*) '******* END FAR FIELD WITH SLOW METHOD **********'
         write(*,*) ' '

      endif

      
    
c     **********************************************************
c     calcul du champ total diffracte ainsi que du champ image
c     **********************************************************            
      if (nenergie.eq.1) then
         if (ncote.eq.0) then
            write(*,*) '**********************************************'      
            write(*,*) '******* COMPUTE ENERGY CONSERVATION **********'
            write(*,*) '**********************************************'
         else
            write(*,*) '***********************************************'      
            write(*,*) '**** COMPUTE DIFFRACTED FIELD WITH INCIDENT ***'
            write(*,*) '***********************************************'
         endif

         
c     on calcule e(k||) tel que E(r||,z)=Int ( e(k||) *exp(ik||.r||
c     +signe*i* gamma(k||)* z )dk||
c     signe=+1 e(k||) correspond au champ transmis total mesuré dans le
c     plan de c fourier du lentille de NA=1 (plan focal mis en 0)
c     signe=-1 correspond au champ reflechi
c     On a e(k||)=ediffkzpos(k||)/(-2i*pi*gamma)+etrans_reference(k||) signe=+1
c     On a e(k||)=ediffkzneg(k||)/(-2i*pi*gamma)+eref_reference(k||) signe=-1
c     ici la difficulté est que le champ diffracte est calculé comme :
c     Ediff(r)=epx(ik0r)/r * ediffkzpos(k||) alors que le champ de reference est :
c     Eref(r)=Int(etrans_reference(k||)exp(ik||.r||+igamma*z)
c     Dans le cas d'une onde plane, etrans_reference=dirac(k||-k||ref)*Tref
c     Numeriquement le dirac est egal à 0 partout et 1/dk||² sur le pixel kref.
         
         tmp=1.d0
         if (nquickdiffracte.eq.1) then
            deltax=aretecube      
         endif
         
c         write(*,*) 'energy',imaxk0,ncote,nfft2dtmp,nfft2d ,nepsmax,neps
c     $        ,k0,deltax,tmp,E0
c     Toutes les ondes Ondes gaussiennes:
         if (beam(1:11).eq.'gwavelinear' .or. beam(1:13).eq
     $        .'gwavecircular' .or. beam(1:7).eq.'speckle' .or.
     $        beam(1:8).eq.'gwaveiso') then

            write(*,*) 'wave based on FFT (Gauss, speckle)'
            
            call champtotalgauss(imaxk0,ncote,nfft2dtmp,nfft2d
     $           ,nepsmax,neps,epscouche,k0,deltax,tmp
     $           ,Ediffkzpos,Ediffkzneg,Efourierincxneg
     $           ,Efourierincyneg ,Efourierinczneg,Efourierincxpos
     $           ,Efourierincypos,Efourierinczpos,fluxreftot
     $           ,fluxtratot)

            efficacite=(fluxreftot+fluxtratot)/fluxinc
            efficaciteref=fluxreftot/fluxinc
            efficacitetrans=fluxtratot/fluxinc
            
         elseif (beam(1:11).eq.'pwavelinear') then
            call champtotallineaire(imaxk0,ncote,nfft2dtmp,nfft2d
     $           ,nepsmax,neps,epscouche,k0,deltax,tmp ,Ediffkzpos
     $           ,Ediffkzneg,zcouche,E0,ss,pp,theta,phi ,Efourierincxneg
     $           ,Efourierincyneg ,Efourierinczneg ,Efourierincxpos
     $           ,Efourierincypos,Efourierinczpos ,fluxreftot
     $           ,fluxtratot ,fluxinc,nstop ,infostr)
            
            efficacite=(fluxreftot+fluxtratot)/fluxinc
            efficaciteref=fluxreftot/fluxinc
            efficacitetrans=fluxtratot/fluxinc

         elseif (beam(1:13).eq.'pwavecircular') then
            write(*,*) 'circulaire'
            call champtotalcirc(imaxk0,ncote,nfft2dtmp,nfft2d ,nepsmax
     $           ,neps,epscouche,k0,deltax,tmp ,Ediffkzpos,Ediffkzneg
     $           ,zcouche,E0,ss,theta,phi ,Efourierincxneg
     $           ,Efourierincyneg ,Efourierinczneg ,Efourierincxpos
     $           ,Efourierincypos,Efourierinczpos ,fluxreftot
     $           ,fluxtratot ,fluxinc,nstop ,infostr)

            efficacite=(fluxreftot+fluxtratot)/fluxinc
            efficaciteref=fluxreftot/fluxinc
            efficacitetrans=fluxtratot/fluxinc

         elseif (beam(1:15).eq.'wavelinearmulti') then
            call champtotallineairemulti(imaxk0,ncote,nfft2dtmp ,nfft2d
     $           ,nepsmax,neps,epscouche,k0,deltax,tmp ,Ediffkzpos
     $           ,Ediffkzneg,zcouche,E0m,ssm ,ppm,thetam ,phim,nbinc
     $           ,Efourierincxneg ,Efourierincyneg ,Efourierinczneg
     $           ,Efourierincxpos ,Efourierincypos ,Efourierinczpos
     $           ,fluxreftot ,fluxtratot,fluxinc ,nstop ,infostr)

            
            efficacite=(fluxreftot+fluxtratot)/fluxinc
            efficaciteref=fluxreftot/fluxinc
            efficacitetrans=fluxtratot/fluxinc
            
c     rajoutter le speckle et iso avec Ar/(delta k)^2+Ediff/(-2 i pi
c     gamma)
         else
            write(*,*) 'onde pas encore faite'
            nstop=1
            infostr='Wave not done'
            return
         endif
         if (ncote.eq.0) then
            write(*,*) 'Incident flux   :',fluxinc
            write(*,*) 'Reflected flux  :',fluxreftot
            write(*,*) 'Transmitted flux:',fluxtratot
            write(*,*) 'Total flux      :',fluxreftot +fluxtratot
            
            write(*,*) 'Conservation of energy:',efficacite
            write(*,*) 'Absorptivity          :',1.d0-efficacite
            write(*,*) 'Reflextivity          :',efficaciteref
            write(*,*) 'Transmittivity        :',efficacitetrans
      
            write(*,*) '************ END ENERGY CONSERVATION ********'
            write(*,*) ' '
         else
            write(*,*) '**********************************************'      
            write(*,*) '***** END DIFFRACTED FIELD WITH INCIDENT *****'
            write(*,*) '**********************************************'
         endif
      endif
      
      if (nlentille.eq.1) then

         write(*,*) '*************************************************'      
         write(*,*) '************* COMPUTE MICROSCOPY ****************'
         write(*,*) '*************************************************'
         if (ncote.eq.0.or.ncote.eq.-1) then
            write(*,*) 'Microscopy half angle in reflexion   :'
     $           ,dasin(numaperref)*180.d0/pi
         endif
         if (ncote.eq.0.or.ncote.eq.1) then
            write(*,*) 'Microscopy half angle in transmission:'
     $           ,dasin(numapertra)*180.d0/pi
         endif

         if (ntypemic.ne.0.and.(nprochefft.ge.1.or.nlecture1.eq.1)) then
            write(*,*)
     $           'recompute Green function if large far field used'
            if (ninterp.eq.0) then
               call fonctiongreensurfcomp(hcc,tolinit,epsabs,xswf,yswf
     $              ,zswf,aretecube,k0,neps,nepsmax,dcouche,zcouche
     $              ,epscouche,ndipole ,nmax,n1m,nzm,nz,nbs,nmat,nmatim
     $              ,nplanm,Tabzn,a ,matind ,matindplan,matindice
     $              ,matrange,nt)
               call fonctiongreensurffft(nx,ny,nz,nx2,ny2,nxm,nym,n1m
     $              ,nzm,nplanm,nmatim,nbs,ntotalm,aretecube,a,matind
     $              ,matindplan,matindice,matrange,b11,b12,b13,b22 ,b23
     $              ,b31,b32,b33,a11,a12,a13,a22,a23,a31,a32,a33,planb)
            else
               call  fonctiongreensurfcompinterp(hcc,tolinit,epsabs,nx
     $              ,ny,nz,zswf,aretecube,k0,neps,nepsmax,dcouche
     $              ,zcouche,epscouche,ndipole,nmax,n1m,nzm,nz,nbs,nmat
     $              ,nmatim,nplanm,Tabzn,a,matind,matindplan,matindice
     $              ,matrange,ninter ,ninterp,nt)

               call fonctiongreensurfinterpfft(nx,ny,nz,nx2,ny2,nxm,nym
     $              ,n1m,nzm,nplanm,nmatim,nbs,ntotalm,ninter,ninterp
     $              ,aretecube,a ,matind ,matindplan ,matindice
     $              ,matrange,b11,b12,b13,b22,b23 ,b31 ,b32,b33,a11 ,a12
     $              ,a13 ,a22,a23,a31 ,a32,a33,planb)
            endif

         endif
            
         if (ntypemic.eq.1) then

            npolainc=0
            call microssurfbf(xi,xr,nbsphere,ndipole,nx,ny,nz,nx2,ny2
     $           ,nxm,nym,nzm,nplanm,nmatim,ntotalm,nmax,matindplan
     $           ,matindice ,Tabdip,b31 ,b32 ,b33,FF,FF0,FFloc,b11,b12
     $           ,b13,a11,a12,a13 ,a22,a23 ,a31 ,a32 ,a33 ,WRK,epscouche
     $           ,zcouche,neps,nepsmax ,xs,ys,zs,nlar,ldabi,polarisa
     $           ,methodeit ,nrig,ncote,tolinit ,aretecube ,npolainc
     $           ,nquicklens ,eps0,k0 ,P0 ,w0 ,nfft2d ,nproche
     $           ,Eimagexpos ,Eimageypos ,Eimagezpos, Eimageincxpos
     $           ,Eimageincypos ,Eimageinczpos, Efourierxpos,
     $           Efourierypos ,Efourierzpos, Efourierincxpos
     $           ,Efourierincypos, Efourierinczpos, Eimagexneg
     $           ,Eimageyneg ,Eimagezneg, Eimageincxneg,Eimageincyneg
     $           ,Eimageinczneg, Efourierxneg ,Efourieryneg,Efourierzneg
     $           , Efourierincxneg ,Efourierincyneg ,Efourierinczneg
     $           ,Ediffkzpos,Ediffkzneg, kxy ,xy,numaperref,numapertra
     $           ,numaperinc ,gross,zlensr,zlenst ,ntypemic , planf
     $           ,planb ,plan2f ,plan2b ,nmatf,file_id ,group_idmic
     $           ,nstop ,infostr)
            
         elseif (ntypemic.eq.2) then
            npolainc=0
            call microssurfdf(xi,xr,nbsphere,ndipole,nx,ny,nz,nx2,ny2
     $           ,nxm,nym,nzm,nplanm,nmatim,ntotalm,nmax,matindplan
     $           ,matindice ,Tabdip,b31 ,b32 ,b33,FF,FF0,FFloc,b11,b12
     $           ,b13,a11,a12,a13 ,a22,a23 ,a31 ,a32 ,a33 ,WRK,epscouche
     $           ,zcouche,neps,nepsmax ,xs,ys,zs,nlar,ldabi,polarisa
     $           ,methodeit ,nrig,ncote,tolinit ,aretecube ,npolainc
     $           ,nquicklens ,eps0,k0 ,P0 ,w0 ,nfft2d ,nproche
     $           ,Eimagexpos ,Eimageypos ,Eimagezpos, Eimageincxpos
     $           ,Eimageincypos ,Eimageinczpos, Efourierxpos,
     $           Efourierypos ,Efourierzpos, Efourierincxpos
     $           ,Efourierincypos, Efourierinczpos, Eimagexneg
     $           ,Eimageyneg ,Eimagezneg, Eimageincxneg,Eimageincyneg
     $           ,Eimageinczneg, Efourierxneg ,Efourieryneg,Efourierzneg
     $           , Efourierincxneg ,Efourierincyneg ,Efourierinczneg
     $           ,Ediffkzpos,Ediffkzneg, kxy ,xy,numaperref,numapertra
     $           ,numaperinc ,gross,zlensr,zlenst,ntypemic,planf,planb
     $           ,plan2f,plan2b ,nmatf,file_id ,group_idmic  ,nstop
     $           ,infostr)
         else
            write(*,*) 'holographic microscopy'
            call cpu_time(t1)
            call date_and_time(date,time,zone,values)
         if (nquickdiffracte.eq.1) deltax=aretecube
c     reduit l'ouverture numérique a numaper
         if (ncote.eq.0.or.ncote.eq.-1) then  
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,indice)   
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)             
            do i=-imaxk0,imaxk0
               do j=-imaxk0,imaxk0
                  ii=imaxk0+i+1
                  jj=imaxk0+j+1
                  kx=dble(i)*deltakx
                  ky=dble(j)*deltaky
                  if (indice0*indice0*k0*k0*numaperref*numaperref-kx*kx
     $                 -ky*ky.le.0.d0) then
                     indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                     Efourierincxneg(indice)=0.d0
                     Efourierincyneg(indice)=0.d0
                     Efourierinczneg(indice)=0.d0
                  endif                     
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         endif
         if (ncote.eq.0.or.ncote.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,indice)   
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            do i=-imaxk0,imaxk0
               do j=-imaxk0,imaxk0
                  ii=imaxk0+i+1
                  jj=imaxk0+j+1
                  kx=dble(i)*deltakx
                  ky=dble(j)*deltaky
                  if (indicen*indicen*k0*k0*numapertra*numapertra-kx*kx
     $                 -ky*ky.le.0.d0) then
                     indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                     Efourierincxpos(indice)=0.d0
                     Efourierincypos(indice)=0.d0
                     Efourierinczpos(indice)=0.d0
                  endif                     
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         endif


         
c     passe le champ diffracte lointain en amplitude e(k||)

         call diffractefft2dtoepos(Ediffkzpos,Ediffkzneg ,Efourierxpos
     $        ,Efourierypos,Efourierzpos,Efourierxneg ,Efourieryneg
     $        ,Efourierzneg,epscouche,nepsmax,neps ,numaperref
     $        ,numapertra,k0,deltax,imaxk0,nfft2dtmp,nfft2d,ncote ,nstop
     $        ,infostr)
         deltakx=2.d0*pi/(dble(nfft2d)*deltax)
         if (nstop.eq.1) return

         write(*,*) 'NA                       :',numaperref*indice0
     $        ,numapertra*indicen
         write(*,*) 'Focal point reflexion    :',zlensr,'m'
         write(*,*) 'Focal point transmission :',zlenst,'m'
         write(*,*) 'Number of point in NA    :',2*imaxk0+1
         write(*,*) 'Delta k                  :',deltakx,'m-1'

         if (zlenst.ne.0.d0.and.ncote.ge.0) then
            write(*,*) 'Change focal point in transmission', zlenst
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kz,ctmp,indice)   
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)          
            do i=-imaxk0,imaxk0
               do j=-imaxk0,imaxk0
                  ii=imaxk0+i+1
                  jj=imaxk0+j+1
                  kx=dble(i)*deltakx
                  ky=dble(j)*deltaky
                  if (indicen*indicen*k0*k0*numapertra*numapertra
     $                 *0.9999d0-kx*kx-ky*ky.gt.0.d0) then
                     kz=dsqrt(indicen*indicen*k0*k0-kx*kx-ky*ky) 
                     ctmp=cdexp(icomp*kz*zlenst)
                     indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                     Efourierxpos(indice)=Efourierxpos(indice)*ctmp
                     Efourierypos(indice)=Efourierypos(indice)*ctmp
                     Efourierzpos(indice)=Efourierzpos(indice)*ctmp
                     Efourierincxpos(indice)=Efourierincxpos(indice)
     $                    *ctmp
                     Efourierincypos(indice)=Efourierincypos(indice)
     $                    *ctmp
                     Efourierinczpos(indice)=Efourierinczpos(indice)
     $                    *ctmp
                  endif
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
         endif
         if (zlensr.ne.0.d0.and.ncote.le.0) then
            write(*,*) 'Change focal point in reflexion',zlensr
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kz,ctmp,indice)   
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)          
            do i=-imaxk0,imaxk0
               do j=-imaxk0,imaxk0
                  ii=imaxk0+i+1
                  jj=imaxk0+j+1
                  kx=dble(i)*deltakx
                  ky=dble(j)*deltaky
                  if (indice0*indice0*k0*k0*numaperref*numaperref
     $                 *0.9999d0-kx*kx-ky*ky.gt.0.d0) then
                     kz=dsqrt(indice0*indice0*k0*k0-kx*kx-ky*ky) 
                     ctmp=cdexp(-icomp*kz*zlensr)
                     indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                     Efourierxneg(indice)=Efourierxneg(indice)*ctmp
                     Efourieryneg(indice)=Efourieryneg(indice)*ctmp
                     Efourierzneg(indice)=Efourierzneg(indice)*ctmp
                     Efourierincxneg(indice)=Efourierincxneg(indice)
     $                    *ctmp
                     Efourierincyneg(indice)=Efourierincyneg(indice)
     $                    *ctmp
                     Efourierinczneg(indice)=Efourierinczneg(indice)
     $                    *ctmp
                  endif
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
         endif
         

c     sauve le champ diffracte total en kx,ky
c     sauve le champ dans l'espace de fourier total en x,y
         if (nmatf.eq.0) then
            do j=-imaxk0,imaxk0             
               kx=dble(j)*deltakx
               ky=dble(j)*deltaky
               write(650,*) kx/k0
               write(651,*) ky/k0
            enddo
            close(650)
            close(651)
         elseif (nmatf.eq.2) then
            do j=-imaxk0,imaxk0             
               kxy(j+imaxk0+1)=dble(j)*deltakx/k0
            enddo
            dim(1)=2*imaxk0+1
            dim(2)=nfft2d
            datasetname='kx Fourier'
            call hdf5write1d(group_idmic,datasetname,kxy,dim)
            datasetname='ky Fourier'
            call hdf5write1d(group_idmic,datasetname,kxy,dim)
         endif
c     do i=1,nfft2d*nfft2d
c     write(*,*) 'ptoui',i,Efourierincxneg(i),Efourierincyneg(i)
c     $           ,Efourierinczneg(i)
c     enddo
         if (ncote.eq.0.or.ncote.eq.-1) then

            if (gross.eq.-1.d0) then
               write(*,*) 'magnification -1',gross
               
               call passagefourierimage(Efourierincxneg ,Efourierincyneg
     $              ,Efourierinczneg,Eimageincxneg ,Eimageincyneg
     $              ,Eimageinczneg,nfft2d,nfft2d,imaxk0,deltakx ,deltax
     $              ,plan2b)
               call passagefourierimage(Efourierxneg,Efourieryneg
     $              ,Efourierzneg,Eimagexneg,Eimageyneg,Eimagezneg
     $              ,nfft2d,nfft2d,imaxk0,deltakx,deltax,plan2b)
            else
               write(*,*) 'Magnification',gross
               sidemic=-1.d0
               call passagefourierimagegross(Efourierincxneg
     $              ,Efourierincyneg,Efourierinczneg,Eimageincxneg
     $              ,Eimageincyneg,Eimageinczneg,nfft2d,nfft2d ,imaxk0
     $              ,deltakx,deltax,gross,k0,indice0,sidemic,plan2f
     $              ,plan2b)
               call passagefourierimagegross(Efourierxneg ,Efourieryneg
     $              ,Efourierzneg,Eimagexneg,Eimageyneg ,Eimagezneg
     $              ,nfft2d,nfft2d,imaxk0,deltakx,deltax ,gross,k0
     $              ,indice0,sidemic,plan2f,plan2b)
            endif

c            write(*,*) 'Number of point in NA',imaxk0*2+1
            if (nmatf.eq.0) then
               do j=-imaxk0,imaxk0
                  do i=-imaxk0,imaxk0
                     k=(i+nfft2d2+1)+nfft2d*(j+nfft2d2)
c     sauve le champ diffracté total (reflechi)
                     write(604,*) dsqrt(cdabs(Efourierincxneg(k))**2
     $                    +cdabs(Efourierincyneg(k))**2
     $                    +cdabs(Efourierinczneg(k))**2)
                     write(605,*) dreal(Efourierincxneg(k))
     $                    ,dimag(Efourierincxneg(k))
                     write(606,*) dreal(Efourierincyneg(k))
     $                    ,dimag(Efourierincyneg(k))
                     write(607,*) dreal(Efourierinczneg(k))
     $                    ,dimag(Efourierinczneg(k))
c     sauve le champ diffracté seul (reflechi)
                     write(612,*) dsqrt(cdabs(Efourierxneg(k))**2
     $                    +cdabs(Efourieryneg(k))**2
     $                    +cdabs(Efourierzneg(k))**2)
                     write(613,*) dreal(Efourierxneg(k))
     $                    ,dimag(Efourierxneg(k))
                     write(614,*) dreal(Efourieryneg(k))
     $                    ,dimag(Efourieryneg(k))
                     write(615,*) dreal(Efourierzneg(k))
     $                    ,dimag(Efourierzneg(k))
                  enddo
               enddo
               do indice=1,nfft2d*nfft2d
c     sauve le champ total dans le plan image.
                  write(620,*) dsqrt(cdabs(Eimageincxneg(indice))**2
     $                 +cdabs(Eimageincyneg(indice))**2
     $                 +cdabs(Eimageinczneg(indice))**2)
                  write(621,*) dreal(Eimageincxneg(indice))
     $                 ,dimag(Eimageincxneg(indice))
                  write(622,*) dreal(Eimageincyneg(indice))
     $                 ,dimag(Eimageincyneg(indice))
                  write(623,*) dreal(Eimageinczneg(indice))
     $                 ,dimag(Eimageinczneg(indice))
c     sauve le champ dans le plan de image.
                  write(628,*) dsqrt(cdabs(Eimagexneg(indice))**2
     $                 +cdabs(Eimageyneg(indice))**2
     $                 +cdabs(Eimagezneg(indice))**2)
                  write(629,*) dreal(Eimagexneg(indice))
     $                 ,dimag(Eimagexneg(indice))
                  write(630,*) dreal(Eimageyneg(indice))
     $                 ,dimag(Eimageyneg(indice))
                  write(631,*) dreal(Eimagezneg(indice))
     $                 ,dimag(Eimagezneg(indice))
               enddo
               close(604)
               close(605)
               close(606)
               close(607)
               close(612)
               close(613)
               close(614)
               close(615)
               close(620)
               close(621)
               close(622)
               close(623)
               close(628)
               close(629)
               close(630)
               close(631)
            elseif (nmatf.eq.2) then
               k=1
               name='Fourier kz<0'
               call writehdf5mic(Efourierxneg,Efourieryneg,Efourierzneg
     $              ,nfft2d,imaxk0,Ediffkzneg,k,name,group_idmic)
               name='Fourier+incident kz<0'
               call writehdf5mic(Efourierincxneg,Efourierincyneg
     $              ,Efourierinczneg,nfft2d,imaxk0,Ediffkzneg,k,name
     $              ,group_idmic)
               k=0
               name='Image kz<0'
               call writehdf5mic(Eimagexneg,Eimageyneg,Eimagezneg,nfft2d
     $              ,imaxk0,Ediffkzneg,k,name,group_idmic)
               name='Image+incident kz<0'
               call writehdf5mic(Eimageincxneg,Eimageincyneg
     $              ,Eimageinczneg,nfft2d,imaxk0,Ediffkzneg,k,name
     $              ,group_idmic)
            endif
         endif
         
         if (ncote.eq.0.or.ncote.eq.1) then

            if (gross.eq.-1.d0) then
               call passagefourierimage(Efourierincxpos
     $              ,Efourierincypos,Efourierinczpos,Eimageincxpos
     $              ,Eimageincypos,Eimageinczpos,nfft2d,nfft2d,imaxk0
     $              ,deltakx,deltax,plan2b)
               call passagefourierimage(Efourierxpos,Efourierypos
     $              ,Efourierzpos,Eimagexpos,Eimageypos,Eimagezpos
     $              ,nfft2d,nfft2d,imaxk0,deltakx,deltax,plan2b)
            else
               sidemic=1.d0
               call passagefourierimagegross(Efourierincxpos
     $              ,Efourierincypos,Efourierinczpos,Eimageincxpos
     $              ,Eimageincypos,Eimageinczpos,nfft2d,nfft2d,imaxk0
     $              ,deltakx,deltax,gross,k0,indicen,sidemic,plan2f
     $              ,plan2b)
               call passagefourierimagegross(Efourierxpos ,Efourierypos
     $              ,Efourierzpos,Eimagexpos,Eimageypos ,Eimagezpos
     $              ,nfft2d,nfft2d,imaxk0,deltakx ,deltax ,gross,k0
     $              ,indicen,sidemic,plan2f,plan2b)
            endif
            if (nmatf.eq.0) then
               do j=-imaxk0,imaxk0
                  do i=-imaxk0,imaxk0
                     k=(i+nfft2d2+1)+nfft2d*(j+nfft2d2)
c$$$  if (i.eq.0.and.j.eq.0) then
c$$$  write(*,*) 'xyzinc',Efourierincxpos(k)
c$$$  $                    ,Efourierincypos(k),Efourierinczpos(k)
c$$$  write(*,*) 'xyz',Efourierxpos(k)
c$$$  $                    ,Efourierypos(k),Efourierzpos(k)
c$$$  write(*,*) 'xyzdiff',Efourierincxpos(k)
c$$$  $                    -Efourierxpos(k),Efourierincypos(k)
c$$$  $                    -Efourierypos(k),Efourierinczpos(k)
c$$$  $                    -Efourierzpos(k)
c$$$  endif
c     sauve le champ diffracté total (transmis)
                     write(600,*) dsqrt(cdabs(Efourierincxpos(k))**2
     $                    +cdabs(Efourierincypos(k))**2
     $                    +cdabs(Efourierinczpos(k))**2)
                     write(601,*) dreal(Efourierincxpos(k))
     $                    ,dimag(Efourierincxpos(k))
                     write(602,*) dreal(Efourierincypos(k))
     $                    ,dimag(Efourierincypos(k))
                     write(603,*) dreal(Efourierinczpos(k))
     $                    ,dimag(Efourierinczpos(k))
c     
                     
c     sauve le champ diffracté seul (transmis)
                     write(608,*) dsqrt(cdabs(Efourierxpos(k))**2
     $                    +cdabs(Efourierypos(k))**2
     $                    +cdabs(Efourierzpos(k))**2)
                     write(609,*) dreal(Efourierxpos(k))
     $                    ,dimag(Efourierxpos(k))
                     write(610,*) dreal(Efourierypos(k))
     $                    ,dimag(Efourierypos(k))
                     write(611,*) dreal(Efourierzpos(k))
     $                    ,dimag(Efourierzpos(k))
                  enddo
               enddo
               do indice=1,nfft2d*nfft2d
c     sauve le champ total dans le plan image.
                  write(616,*) dsqrt(cdabs(Eimageincxpos(indice))**2
     $                 +cdabs(Eimageincypos(indice))**2
     $                 +cdabs(Eimageinczpos(indice))**2)
                  write(617,*) dreal(Eimageincxpos(indice))
     $                 ,dimag(Eimageincxpos(indice))
                  write(618,*) dreal(Eimageincypos(indice))
     $                 ,dimag(Eimageincypos(indice))
                  write(619,*) dreal(Eimageinczpos(indice))
     $                 ,dimag(Eimageinczpos(indice))
c     sauve le champ dans le plan image.
                  write(624,*) dsqrt(cdabs(Eimagexpos(indice))**2
     $                 +cdabs(Eimageypos(indice))**2
     $                 +cdabs(Eimagezpos(indice))**2)
                  write(625,*) dreal(Eimagexpos(indice))
     $                 ,dimag(Eimagexpos(indice))
                  write(626,*) dreal(Eimageypos(indice))
     $                 ,dimag(Eimageypos(indice))
                  write(627,*) dreal(Eimagezpos(indice))
     $                 ,dimag(Eimagezpos(indice))

               enddo
               close(600)
               close(601)
               close(602)
               close(603)
               close(608)
               close(609)
               close(610)
               close(611)
               
               close(616)
               close(617)
               close(618)
               close(619)
               close(624)
               close(625)
               close(626)
               close(627)
            elseif (nmatf.eq.2) then
               k=1
               name='Fourier kz>0'
               call writehdf5mic(Efourierxpos,Efourierypos,Efourierzpos
     $              ,nfft2d,imaxk0,Ediffkzpos,k,name,group_idmic)
               name='Fourier+incident kz>0'
               call writehdf5mic(Efourierincxpos,Efourierincypos
     $              ,Efourierinczpos,nfft2d,imaxk0,Ediffkzpos,k,name
     $              ,group_idmic)
               k=0
               name='Image kz>0'
               call writehdf5mic(Eimagexpos,Eimageypos,Eimagezpos,nfft2d
     $              ,imaxk0,Ediffkzpos,k,name,group_idmic)
               name='Image+incident kz>0'
               call writehdf5mic(Eimageincxpos,Eimageincypos
     $              ,Eimageinczpos,nfft2d,imaxk0,Ediffkzpos,k,name
     $              ,group_idmic)
            endif
         endif
         
         call cpu_time(t2)
         call date_and_time(date,time,zone,values2)
         message='to compute holography'
         call calculatedate(values2,values,t2,t1,message)

c     sauve les coordonnees cartesiennes:
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)  
         do i=1,nfft2d
            xy(i)=deltax*dble(i-nfft2d2-1)*gross
            kxy(i)=deltakx*dble(i-nfft2d2-1)/k0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         if (nmatf.eq.0) then
            open(309,file='ximage.mat')
            open(310,file='yimage.mat')
            open(311,file='kximage.mat')
            open(312,file='kyimage.mat')
            do i=1,nfft2d
               write(309,*) xy(i)
               write(310,*) xy(i)
               write(311,*) kxy(i)
               write(312,*) kxy(i)              
            enddo
            close(309)
            close(310)
            close(311)
            close(312)
         elseif (nmatf.eq.2) then
            dim(1)=nfft2d
            dim(2)=nfft2d
c            datasetname='kx Fourier'
c            call hdf5write1d(group_idmic,datasetname,kxy,dim)
            datasetname='x Image'
            call hdf5write1d(group_idmic,datasetname,xy,dim)
c             datasetname='ky Fourier'
c            call hdf5write1d(group_idmic,datasetname,kxy,dim)
            datasetname='y Image'
            call hdf5write1d(group_idmic,datasetname,xy,dim)
         endif
      endif
      write(*,*) '*************** END MICROSCOPY ******************'
      write(*,*)
      endif
      

c$$$  if (nforce.eq.1) then
c$$$  write(99,*) '****** Compute the optical force *********'
c$$$  c     Begin the computation of the force
c$$$  c     Compute the FFT of the dipole
c$$$  do k=1,nz
c$$$  do j=1,ny
c$$$  do i=1,nx
c$$$  ii=i+nx*(j-1)+nx*ny*(k-1)
c$$$  indice=i+nx2*(j-1)+nxy2*(k-1)
c$$$  if (Tabdip(ii).eq.0) then
c$$$  vectx(indice)=0.d0
c$$$  vecty(indice)=0.d0
c$$$  vectz(indice)=0.d0
c$$$  else
c$$$  jj=3*Tabdip(ii)
c$$$  vectx(indice)=FF(jj-2)
c$$$  vecty(indice)=FF(jj-1)
c$$$  vectz(indice)=FF(jj)
c$$$  endif
c$$$  vectx(i+nx+nx2*(j-1)+nxy2*(k-1))=0.d0
c$$$  vectx(i+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
c$$$  vectx(i+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
c$$$  vectx(i+nx+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
c$$$  vectx(i+nx+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
c$$$  vectx(i+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
c$$$  vectx(i+nx+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
c$$$  
c$$$  vecty(i+nx+nx2*(j-1)+nxy2*(k-1))=0.d0
c$$$  vecty(i+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
c$$$  vecty(i+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
c$$$  vecty(i+nx+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
c$$$  vecty(i+nx+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
c$$$  vecty(i+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
c$$$  vecty(i+nx+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
c$$$  
c$$$  vectz(i+nx+nx2*(j-1)+nxy2*(k-1))=0.d0
c$$$  vectz(i+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
c$$$  vectz(i+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
c$$$  vectz(i+nx+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
c$$$  vectz(i+nx+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
c$$$  vectz(i+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
c$$$  vectz(i+nx+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  
c$$$  CALL ZFFT3D(vectx,NX2,NY2,NZ2,1)
c$$$  CALL ZFFT3D(vecty,NX2,NY2,NZ2,1)
c$$$  CALL ZFFT3D(vectz,NX2,NY2,NZ2,1)
c$$$  c     do i=1,nx2*ny2*nz2
c$$$  c     write(99,*) 'vect',vectx(i),vecty(i),vectz(i)
c$$$  c     enddo
c$$$  c     Compute the x derivative of the local field
c$$$  test=1 
c$$$  
c$$$  if(beam(1:11).eq.'gwavelinear'.or.beam(1:13).eq
c$$$  $        .'gwavecircular'.or.beam(1:15).eq.'gparawavelinear'
c$$$  $        .or.beam(1:17).eq.'gparawavecircular') then
c$$$  open(150,file='dergaussy.mat',form='unformatted')
c$$$  open(151,file='dergaussz.mat',form='unformatted')
c$$$  elseif (beam(1:9).eq.'arbitrary') then
c$$$  call incidentarbitrarydercreate(namefileinc)
c$$$  endif
c$$$  call derivativefield(aretecube,k0,nx,ny,nz,nx2,ny2,nz2,nxy2
c$$$  $        ,ntotal,ntotalm,nmax,FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
c$$$  $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
c$$$  $        ,vectbx,vectby,vectbz,test,FF0)
c$$$  do i=1,ndipole
c$$$  ii=3*i
c$$$  if (Tabdip(i).ne.0) then
c$$$  j=Tabdip(i)  
c$$$  if (beam(1:11).eq.'pwavelinear') then
c$$$  call  ondeplaned(xs(j),ys(j),zs(j),k0,E0,ss,pp,theta
c$$$  $                 ,phi,test,Eder)  
c$$$  elseif(beam(1:13).eq.'pwavecircular') then
c$$$  call  ondecirced(xs(j),ys(j),zs(j),k0,E0,ss,theta,phi
c$$$  $                 ,test,Eder) 
c$$$  elseif(beam(1:11).eq.'gwavelinear') then
c$$$  call  gaussianchampd(xs(j),ys(j),zs(j),xgaus,ygaus
c$$$  $                 ,zgaus,theta,phi,w0,k0,ss,pp,E0,Eder,tol,nloin)
c$$$  write(150) Eder(1,2)
c$$$  write(150) Eder(2,2)
c$$$  write(150) Eder(3,2)
c$$$  write(151) Eder(1,3)
c$$$  write(151) Eder(2,3)
c$$$  write(151) Eder(3,3)
c$$$  elseif(beam(1:13).eq.'gwavecircular') then
c$$$  call  gaussianchampdcirc(xs(j),ys(j),zs(j),xgaus,ygaus
c$$$  $                 ,zgaus,theta,phi,w0,k0,ss,E0,Eder,tol,nloin)
c$$$  write(150) Eder(1,2)
c$$$  write(150) Eder(2,2)
c$$$  write(150) Eder(3,2)
c$$$  write(151) Eder(1,3)
c$$$  write(151) Eder(2,3)
c$$$  write(151) Eder(3,3)
c$$$  elseif(beam(1:15).eq.'gparawavelinear') then
c$$$  call  gaussianparalineard(xs(j),ys(j),zs(j),xgaus
c$$$  $                 ,ygaus,zgaus,theta,phi,w0,k0,ss,pp,E0,Eder,nstop
c$$$  $                 ,infostr)
c$$$  write(150) Eder(1,2)
c$$$  write(150) Eder(2,2)
c$$$  write(150) Eder(3,2)
c$$$  write(151) Eder(1,3)
c$$$  write(151) Eder(2,3)
c$$$  write(151) Eder(3,3)
c$$$  elseif(beam(1:17).eq.'gparawavecircular') then
c$$$  call  gaussianparacircd(xs(j),ys(j),zs(j),xgaus,ygaus
c$$$  $                 ,zgaus,theta,phi,w0,k0,ss,E0,Eder,nstop,infostr)
c$$$  write(150) Eder(1,2)
c$$$  write(150) Eder(2,2)
c$$$  write(150) Eder(3,2)
c$$$  write(151) Eder(1,3)
c$$$  write(151) Eder(2,3)
c$$$  write(151) Eder(3,3)
c$$$  elseif (beam(1:9).eq.'arbitrary') then
c$$$  c     call  ondeplaned(xs(j),ys(j),zs(j),k0,E0,ss,pp,theta
c$$$  c     $                 ,phi,test,Eder)  
c$$$  c     write(*,*) 'Eder1',Eder
c$$$  call incidentarbitraryder1(test,xs(j),ys(j),zs(j)
c$$$  $                 ,Eder)
c$$$  c     write(*,*) 'Eder2',Eder
c$$$  endif 
c$$$  Ederivex(j,1)=Eder(1,1)+FF0(ii-2)
c$$$  Ederivex(j,2)=Eder(2,1)+FF0(ii-1)
c$$$  Ederivex(j,3)=Eder(3,1)+FF0(ii)            
c$$$  c     write(*,*) 'i derivee',j
c$$$  c     write(*,*) 'Ederx',Eder(1,1),Eder(2,1),Eder(3,1)
c$$$  c     write(*,*) 'Edery',Eder(1,2),Eder(2,2),Eder(3,2)
c$$$  c     write(*,*) 'Ederz',Eder(1,3),Eder(2,3),Eder(3,3)
c$$$  
c$$$  endif
c$$$  enddo
c$$$  rewind(150)
c$$$  rewind(151)
c$$$  c     Compute the y derivative of the local field
c$$$  test=2
c$$$  call derivativefield(aretecube,k0,nx,ny,nz,nx2,ny2,nz2,nxy2
c$$$  $        ,ntotal,ntotalm,nmax,FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
c$$$  $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
c$$$  $        ,vectbx,vectby,vectbz,test,FF0)
c$$$  do i=1,ndipole
c$$$  ii=3*i
c$$$  if (Tabdip(i).ne.0) then
c$$$  j=Tabdip(i)    
c$$$  if (beam(1:11).eq.'pwavelinear') then
c$$$  call  ondeplaned(xs(j),ys(j),zs(j),k0,E0,ss,pp,theta
c$$$  $                 ,phi,test,Eder) 
c$$$  elseif(beam(1:13).eq.'pwavecircular') then
c$$$  call  ondecirced(xs(j),ys(j),zs(j),k0,E0,ss,theta,phi
c$$$  $                 ,test,Eder)   
c$$$  elseif(beam(1:11).eq.'gwavelinear') then
c$$$  read(150) Eder(1,2)
c$$$  read(150) Eder(2,2)
c$$$  read(150) Eder(3,2)                 
c$$$  elseif(beam(1:13).eq.'gwavecircular') then
c$$$  read(150) Eder(1,2)
c$$$  read(150) Eder(2,2)
c$$$  read(150) Eder(3,2)
c$$$  elseif(beam(1:15).eq.'gparawavelinear') then
c$$$  read(150) Eder(1,2)
c$$$  read(150) Eder(2,2)
c$$$  read(150) Eder(3,2)                 
c$$$  elseif(beam(1:17).eq.'gparawavecircular') then
c$$$  read(150) Eder(1,2)
c$$$  read(150) Eder(2,2)
c$$$  read(150) Eder(3,2)
c$$$  elseif (beam(1:9).eq.'arbitrary') then
c$$$  call incidentarbitraryder1(test,xs(j),ys(j),zs(j)
c$$$  $                 ,Eder)
c$$$  endif 
c$$$  Ederivey(j,1)=Eder(1,2)+FF0(ii-2)
c$$$  Ederivey(j,2)=Eder(2,2)+FF0(ii-1)
c$$$  Ederivey(j,3)=Eder(3,2)+FF0(ii)
c$$$  endif
c$$$  enddo
c$$$  c     Compute the z derivative of the local field
c$$$  test=3
c$$$  call derivativefield(aretecube,k0,nx,ny,nz,nx2,ny2,nz2,nxy2
c$$$  $        ,ntotal,ntotalm,nmax,FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
c$$$  $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
c$$$  $        ,vectbx,vectby,vectbz,test,FF0)
c$$$  do i=1,ndipole
c$$$  ii=3*i
c$$$  if (Tabdip(i).ne.0) then
c$$$  j=Tabdip(i) 
c$$$  if (beam(1:11).eq.'pwavelinear') then
c$$$  call  ondeplaned(xs(j),ys(j),zs(j),k0,E0,ss,pp,theta
c$$$  $                 ,phi,test,Eder)   
c$$$  write(99,*) 'prout',xs(j),ys(j),zs(j),k0,E0,ss,pp
c$$$  $                 ,theta,phi,test,Eder,j,i,ndipole
c$$$  elseif(beam(1:13).eq.'pwavecircular') then
c$$$  call  ondecirced(xs(j),ys(j),zs(j),k0,E0,ss,theta,phi
c$$$  $                 ,test,Eder)   
c$$$  elseif(beam(1:11).eq.'gwavelinear') then
c$$$  read(151) Eder(1,3)
c$$$  read(151) Eder(2,3)
c$$$  read(151) Eder(3,3)                 
c$$$  elseif(beam(1:13).eq.'gwavecircular') then
c$$$  read(151) Eder(1,3)
c$$$  read(151) Eder(2,3)
c$$$  read(151) Eder(3,3)
c$$$  elseif(beam(1:15).eq.'gparawavelinear') then
c$$$  read(151) Eder(1,3)
c$$$  read(151) Eder(2,3)
c$$$  read(151) Eder(3,3)                 
c$$$  elseif(beam(1:17).eq.'gparawavecircular') then
c$$$  read(151) Eder(1,3)
c$$$  read(151) Eder(2,3)
c$$$  read(151) Eder(3,3)
c$$$  elseif (beam(1:9).eq.'arbitrary') then
c$$$  c     call  ondeplaned(xs(j),ys(j),zs(j),k0,uncomp,0.d0,1.d0
c$$$  c     $                 ,0.d0,0.d0,test,Eder) 
c$$$  c     write(*,*) 'Ederz1',Eder
c$$$  call incidentarbitraryder1(test,xs(j),ys(j),zs(j)
c$$$  $                 ,Eder)
c$$$  c     write(*,*) 'Ederz2',Eder
c$$$  endif 
c$$$  Ederivez(j,1)=Eder(1,3)+FF0(ii-2)
c$$$  Ederivez(j,2)=Eder(2,3)+FF0(ii-1)
c$$$  Ederivez(j,3)=Eder(3,3)+FF0(ii)
c$$$  endif
c$$$  enddo
c$$$  close(150)
c$$$  close(151)
c$$$  c     comptute the optical force from the derivative of the local field
c$$$  c     and the dipole
c$$$  do i=1,3*nbsphere
c$$$  FF0(i)=dconjg(FF(i))
c$$$  enddo
c$$$  forcet(1)=0.d0
c$$$  forcet(2)=0.d0
c$$$  forcet(3)=0.d0
c$$$  
c$$$  do i=1,nbsphere 
c$$$  c     write(*,*) 'ooooooooooooo',i,nbsphere
c$$$  k=3*(i-1)
c$$$  force(i,1)=0.5d0*dreal(FF0(k+1)*Ederivex(i,1)+FF0(k+2)
c$$$  $           *Ederivex(i,2)+FF0(k+3)*Ederivex(i,3))*quatpieps0       
c$$$  force(i,2)=0.5d0*dreal(FF0(k+1)*Ederivey(i,1)+FF0(k+2)
c$$$  $           *Ederivey(i,2)+FF0(k+3)*Ederivey(i,3))*quatpieps0       
c$$$  force(i,3)=0.5d0*dreal(FF0(k+1)*Ederivez(i,1)+FF0(k+2)
c$$$  $           *Ederivez(i,2)+FF0(k+3)*Ederivez(i,3))*quatpieps0 
c$$$  c     write(*,*) 'champ',FF0(k+1),FF0(k+2),FF0(k+3)
c$$$  c     write(*,*) 'champdx',Ederivex(i,1),Ederivex(i,2),Ederivex(i
c$$$  c     $           ,3)
c$$$  c     write(*,*) 'champdy',Ederivey(i,1),Ederivey(i,2),Ederivey(i
c$$$  c     $           ,3)
c$$$  c     write(*,*) 'champdz',Ederivez(i,1),Ederivez(i,2),Ederivez(i
c$$$  c     $           ,3)
c$$$  c     write(*,*) 'force',force(i,1),force(i,2),force(i,3)
c$$$  forcet(1)=forcet(1)+force(i,1)
c$$$  forcet(2)=forcet(2)+force(i,2)
c$$$  forcet(3)=forcet(3)+force(i,3)
c$$$  enddo
c$$$  c     write(*,*) 'coucou'
c$$$  if (nstop == -1) then
c$$$  infostr = 'Calculation cancelled at the optical force'
c$$$  return
c$$$  endif
c$$$  c     save the density of force
c$$$  subunit=0
c$$$  if (nforced.eq.1) then
c$$$  do i=1,nbsphere
c$$$  subunit= subunit+1
c$$$  forcex(subunit) = force(i,1)
c$$$  forcey(subunit) = force(i,2)
c$$$  forcez(subunit) = force(i,3)
c$$$  write(60,*) force(i,1)
c$$$  write(61,*) force(i,2)
c$$$  write(62,*) force(i,3)
c$$$  enddo
c$$$  close(60)
c$$$  close(61)
c$$$  close(62)
c$$$  endif
c$$$  
c$$$  forcem=dsqrt(forcet(1)*forcet(1)+forcet(2)*forcet(2)+forcet(3)
c$$$  $        *forcet(3))
c$$$  write(99,*) 'optical force x',forcet(1)
c$$$  write(99,*) 'optical force y',forcet(2)
c$$$  write(99,*) 'optical force z',forcet(3)
c$$$  write(*,*) 'optical force x',forcet(1)
c$$$  write(*,*) 'optical force y',forcet(2)
c$$$  write(*,*) 'optical force z',forcet(3)
c$$$  forcemie=(MIECext-GSCA*MIECsca)/8.d0/pi
c$$$  write(99,*) 'modulus of the force',forcem,forcemie
c$$$  write(99,*) 'modulus of the force',forcem,forcemie
c$$$  endif 
c$$$  
c$$$  if (ntorque*nforce.eq.1) then
c$$$  write(99,*) '********* Compute the optical torque *********'
c$$$  c     Compute the optical torque on the particle 
c$$$  couplet(1)=0.d0
c$$$  couplet(2)=0.d0
c$$$  couplet(3)=0.d0
c$$$  c     compute the center of gravity of the object
c$$$  xg=0.d0
c$$$  zg=0.d0
c$$$  zg=0.d0
c$$$  do i=1,nbsphere
c$$$  xg=xg+xs(i)
c$$$  yg=yg+ys(i)
c$$$  zg=zg+zs(i)
c$$$  enddo
c$$$  xg=xg/dble(nbsphere)
c$$$  yg=yg/dble(nbsphere)
c$$$  zg=zg/dble(nbsphere)
c$$$  c     gamma=r x F+1/2 Re(p^* x p / alpha_0)
c$$$  do i=1,nbsphere
c$$$  kk=3*(i-1)
c$$$  
c$$$  do ii=1,3
c$$$  Em(ii)=2.d0/3.d0*icomp*k03*FF(kk+ii)+FFloc(kk+ii)
c$$$  enddo
c$$$  
c$$$  couple(i,1)=0.5d0*dreal(FF0(kk+2)*Em(3)-FF0(kk+3)*Em(2))
c$$$  $           *quatpieps0+((ys(i)-yg)*force(i,3)-(zs(i)-zg)*force(i
c$$$  $           ,2))
c$$$  couple(i,2)=0.5d0*dreal(-FF0(kk+1)*Em(3)+FF0(kk+3)*Em(1))
c$$$  $           *quatpieps0+((zs(i)-zg)*force(i,1)-(xs(i)-xg)*force(i
c$$$  $           ,3))
c$$$  couple(i,3)=0.5d0*dreal(FF0(kk+1)*Em(2)-FF0(kk+2)*Em(1))
c$$$  $           *quatpieps0+((xs(i)-xg)*force(i,2)-(ys(i)-yg)*force(i
c$$$  $           ,1))
c$$$  couplet(1)=couplet(1)+couple(i,1)
c$$$  couplet(2)=couplet(2)+couple(i,2)
c$$$  couplet(3)=couplet(3)+couple(i,3)               
c$$$  enddo
c$$$  c     save the density of torque
c$$$  if (nstop == -1) then
c$$$  infostr = 'Calculation cancelled at the optical torque'
c$$$  return
c$$$  endif
c$$$  subunit=0
c$$$  if (ntorqued.eq.1) then
c$$$  do i=1,nbsphere
c$$$  subunit= subunit+1
c$$$  torquex(subunit) = couple(i,1)
c$$$  torquey(subunit) = couple(i,2)
c$$$  torquez(subunit) = couple(i,3)
c$$$  write(63,*) couple(i,1)
c$$$  write(64,*) couple(i,2)
c$$$  write(65,*) couple(i,3)
c$$$  enddo
c$$$  endif
c$$$  
c$$$  couplem=dsqrt(couplet(1)*couplet(1)+couplet(2)*couplet(2)
c$$$  $        +couplet(3)*couplet(3))
c$$$  write(99,*) 'optical torque x',couplet(1)
c$$$  write(99,*) 'optical torque y',couplet(2)
c$$$  write(99,*) 'optical torque z',couplet(3)
c$$$  write(99,*) 'modulus of the optical torque',couplem
c$$$  write(*,*) 'couple',Cabs/8.d0/k0/pi*I0
c$$$  endif
c$$$  
c$$$  c     calcul sur des forces et couple sur differents objets si presents
c$$$  if (numberobjet.ne.1) then
c$$$  forcexmulti=0.d0
c$$$  forceymulti=0.d0
c$$$  forcezmulti=0.d0
c$$$  torquexmulti=0.d0
c$$$  torqueymulti=0.d0
c$$$  torquezmulti=0.d0    
c$$$  write(*,*) 'nbsphere=',nbsphere  
c$$$  do i=1,nbsphere
c$$$  is=tabmulti(i)
c$$$  forcexmulti(is)=forcexmulti(is)+force(i,1)
c$$$  forceymulti(is)=forceymulti(is)+force(i,2)
c$$$  forcezmulti(is)=forcezmulti(is)+force(i,3)
c$$$  torquexmulti(is)=torquexmulti(is)+couple(i,1)
c$$$  torqueymulti(is)=torqueymulti(is)+couple(i,2)
c$$$  torquezmulti(is)=torquezmulti(is)+couple(i,3)            
c$$$  enddo
c$$$  endif
      
c     close filesc
c     save the Poynting vecteur
      
c     save the density of the optical force
c     close(60)
c     close(61)
c     close(62)
c     save the density of optical torque
c     close(63)
c     close(64)
c     close(65)

c     output file
      close(99)

c     9999 format(201(d22.15,1x))

      if (beam(1:11).eq.'gwavelinear' .or .beam(1:13).eq
     $     .'gwavecircular' .or. beam(1:7).eq.'speckle' .or.
     $     beam(1:8).eq.'gwaveiso') then
      endif

 999  if (nmatf.eq.2) then
!     fermeture du fichier hdf5
         CALL h5gclose_f(group_idopt,error) 
         CALL h5gclose_f(group_iddip,error)
         CALL h5gclose_f(group_idff,error)
         CALL h5gclose_f(group_idmic,error)
c     CALL h5gclose_f(group_idof,error)
         CALL h5gclose_f(group_idnf,error)
         call hdf5close(file_id)
         write(*,*) 'close h5file'
      endif
      
      if (nstop == -1) then
         infostr = 'Calculation cancelled at the end!'
         return
      endif

      call cpu_time(tf)
      call date_and_time(date,time,zone,valuesf)

      write(*,*) 'COMPLETED'
      infostr='COMPLETED'
      write(*,*) 'end'
      end
