#ifdef USE_HDF5
      use HDF5
#endif
      implicit none
c     integer
      integer ii,jj,kk,i,j,k,nstop
      integer  nlocal,nmacro,nsection,nforce ,nforced
     $     ,ntorque,ntorqued,nsens,nproche,nlecture,nquickdiffracte
     $     ,nrig,ncote,ndiffracte,ninterp,ir,nenergie,nmatf,ntypemic

c     variables for the object
      integer nbsphere3,nbsphere,ndipole,IP(3),test,numberobjetmax
     $     ,numberobjet,ng
      parameter (numberobjetmax=20)
      integer nx,ny,nz,nx2,ny2,nxy2,nz2,nxm,nym,nzm,ntotal,nxmp,nymp
     $     ,nzmp
      integer subunit,nphi,ntheta
c     parameter (nxm=3,nym=3,nzm=3,nphi=72,ntheta=35)
      parameter (nxm=40,nym=40,nzm=40,nphi=144,ntheta=72)
c     definition of the size for the code
      INTEGER nmax, ntotalm, nmaxpp
      character(2) polarizability
c     variables for the positions
      double precision rayon,side,sidex,sidey ,sidez,hauteur
     $     ,xgmulti(numberobjetmax) ,ygmulti(numberobjetmax)
     $     ,zgmulti(numberobjetmax) ,rayonmulti(numberobjetmax),demiaxea
     $     ,demiaxeb,demiaxec ,thetaobj,phiobj,psiobj,lc,hc,density,psi
      double precision aretecube
      DOUBLE PRECISION, DIMENSION(nxm*nym*nzm) :: xs,ys,zs,xswf,yswf
     $     ,zswf
      DOUBLE PRECISION,DIMENSION((ntheta+1)*nphi)::thetafield,phifield
     $     ,poyntingfield
      double precision pi,lambda,lambda10n,k0,k03,epi,epr,c

c     variables for the material
      integer ierror
      double precision eps0
      double complex, dimension(nxm*nym*nzm,3,3) :: polarisa,epsilon
      double complex epsmulti(numberobjetmax)
     $     ,epsanimulti(numberobjetmax,3,3)
      character (64), DIMENSION(numberobjetmax) :: materiaumulti
      character(64) materiau,object,beam,namefileobj,namefileinc
     $     ,filereread
      character(3) trope

c     variables for the incident field and local field
      DOUBLE PRECISION, DIMENSION(nxm*nym*nzm) :: incidentfield,
     $     localfield,macroscopicfield,forcex,forcey,forcez, torquex
     $     ,torquey,torquez
      double precision forcexmulti(numberobjetmax)
     $     ,forceymulti(numberobjetmax),forcezmulti(numberobjetmax)
     $     ,torquexmulti(numberobjetmax),torqueymulti(numberobjetmax)
     $     ,torquezmulti(numberobjetmax)
      integer nbinc
      double precision ss,pp,theta,phi,I0,phim(10)
     $     ,thetam(10),ssm(10),ppm(10)
      double complex Eloc(3),Em(3),E0,uncomp,icomp,zzero,E0m(10)
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldx
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldy
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldz
      double complex, dimension(nxm*nym*nzm) :: localfieldx
      double complex, dimension(nxm*nym*nzm) :: localfieldy
      double complex, dimension(nxm*nym*nzm) :: localfieldz
      double complex, dimension(nxm*nym*nzm) :: incidentfieldx
      double complex, dimension(nxm*nym*nzm) :: incidentfieldy
      double complex, dimension(nxm*nym*nzm) :: incidentfieldz
      double complex propaesplibre(3,3)
      double complex, dimension(3*nxm*nym*nzm) :: FF,FF0,FFloc

c     Green function
      integer, dimension(nxm*nym*nzm) :: Tabdip,Tabmulti,Tabzn
      integer indice,nmat,nbs,n1m,nmatim,nplanm
      parameter(n1m=max(nxm,nym),nmatim =min(nxm,nym) *(2*max(nxm,nym)
     $     -min(nxm,nym)+1)/2)
      parameter (nplanm=nzm*(nzm+1)/2)
      parameter (ntotalm=4*nxm*nym,nbs=nzm*(nzm+1)*nmatim/2)
      double complex matrange(nbs,5)
      integer matindice(nplanm,nmatim),matind(0:2*n1m*n1m)
     $     ,matindplan(nzm,nzm)
      double complex a11(2*nxm,2*nym,nplanm),a12(2*nxm,2*nym,nplanm),
     $     a13(2*nxm,2*nym,nplanm),a22(2*nxm,2*nym,nplanm),a23(2*nxm,2
     $     *nym,nplanm),a31(2*nxm,2*nym,nplanm),a32(2*nxm,2*nym,nplanm)
     $     , a33(2*nxm,2*nym,nplanm),b11(4*nxm*nym),b12(4*nxm*nym),b13(4
     $     *nxm*nym),b22(4*nxm*nym),b23(4*nxm*nym),b31(4 *nxm*nym),b32(4
     $     *nxm*nym),b33(4*nxm*nym)

      integer neps,nepsmax
      parameter (nepsmax=8)
      double precision zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1)
      character (64), DIMENSION(0:nepsmax+1) :: materiaucouche

      double precision forcet(3),forcem,forcemie
      double precision couplet(3),couplem
      double complex Eder(3,3)
      
c     computation of the cross section
      integer iphi,itheta
      double precision MIECEXT,MIECABS,MIECSCA ,GSCA,Cext,normal(3)
     $     ,deltatheta,deltaphi,Csca,Cscai,Cabs,gasym,thetas,phis
      double complex ctmp
      
c     variables for the iterative method
      INTEGER ldabi, nlar
      integer nnnr,ncompte
      integer NLIM,ndim,nou,maxit,nstat,nloop,STEPERR
      DOUBLE PRECISION  NORM,TOL,DZNRM2,norm1,norm2,tolinit,tol1
      double complex ALPHA,BETA,GPETA,DZETA,R0RN

c     COMMON /ONTHEHEAP/ b,xr,xi,wrk
      double complex, dimension(3*nxm*nym*nzm) :: b,xr,xi
      double complex, dimension(3*nxm*nym*nzm,12) :: wrk
      
c     double complex wrk(*), xi(*), xr(*), b(*)
c     POINTER ( xr_p, xr ), ( b_p, b )
c     POINTER ( wrk_p, wrk ), ( xi_p, xi)

c     Poynting vector
      integer nr,nrmax,nw,nwmax
      double precision Poyntinginc

c     Info string
      character(64) infostr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     nouvelle variable a passer en argument d'entree
c     power et diametre
      double precision P0,w0,xgaus,ygaus,zgaus,quatpieps0
      character(12) methodeit

c     nouvelle variable de sortie Irra     
      double precision irra,efficacite,efficaciteref,efficacitetrans
      

c     nouvelle variable
      integer nloin

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Creation des nouvelles variables
      integer na

c     variable pour avoir l'image a travers la lentille
      integer nquicklens,nlentille,nobjet,nfft2d,nfft2d2
      parameter (nfft2d=1024)
      double precision kx,ky,kz,deltakx,deltaky,numaperref,numapertra
     $     ,deltax,gross,numaperinc,zlensr,zlenst
      double precision kxy(nfft2d),xy(nfft2d)
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
     $     *nfft2d),masque(nfft2d,nfft2d) ,Ediffkzpos(nfft2d,nfft2d,3)
     $     ,Ediffkzneg(nfft2d,nfft2d,3)

      character(LEN=100) :: h5file

c     constant
      c=299792458.d0
      quatpieps0=1.d0/(c*c*1.d-7)
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      
c     DATA INPUT
      lambda=632.8d0       !wavelength
      P0=1.d0              !power    
      w0=lambda*10.d0      !waist
      
  
c********************************************************
c     Define the polarizability
c********************************************************

c     polarizability='CM'  !Clausius Mossotti
      polarizability='RR'  !Clausius Mossotti with radiative reaction
c     polarizability='LA'  !Polarizability defines by Lakthakia
c     polarizability='LR'  !Polarizability defines by Draine
c     polarizability='GB'  !Polarizability with first Mie coefficient
c     polarizability='PS'  !Polarizability for a sphere with local correction
c********************************************************
c     End polarizability
c********************************************************
      
c*******************************************************
c     Define the incident wave
c*******************************************************
c     beam='gwavelinear'      !Linear Gaussian wave
c     beam='gwavecircular'    !Circular Gaussian wave
      beam='pwavelinear'      !Linear plane wave
c     beam='pwavecircular'    !Circular plane wave
c     beam='wavelinearmulti'  !Multiple linear plane wave
c     beam='gwaveiso'         !Iso focus wave
c     beam='speckle'          !Specke wave
c     beam='arbitrary'        !Arbitrary wave defines by the user
      
      if (beam(1:11).eq.'pwavelinear') then
         theta=30.d0
         phi=0.d0
         pp=0.d0
         ss=1.d0
      elseif (beam(1:13).eq.'pwavecircular') then
         theta=0.d0
         phi=0.d0
         ss=1.d0 
      elseif (beam(1:11).eq.'gwavelinear') then
         theta=0.d0
         phi=0.d0
         psi=0.d0      
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
         pp=psi
      elseif (beam(1:11).eq.'gwavecircular') then
         theta=0.d0
         phi=0.d0        
         ss=0.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
      elseif (beam(1:7).eq.'speckle') then        
         pp=1.d0
         ss=0.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
         ir=0
      elseif (beam(1:16).eq.'wavelinearmulti') then        
         nbinc=2
         ppm(1)=0.5d0
         ppm(2)=0.5d0
         ssm(1)=1.d0
         ssm(2)=1.d0
         thetam(1)=65.287098590488455d0
         thetam(2)=-65.287098590488455d0
         phim(1)=0.d0
         phim(2)=0.d0
         E0m(1)=(1.d0,0.d0)*cdexp(-icomp*4.d0*pi/3.d0)
         E0m(2)=(1.d0,0.d0)*cdexp(icomp*4.d0*pi/3.d0)
      elseif (beam(1:9).eq.'arbitrary') then
         
         namefileinc='incarbitrary.in'
         open(15,file=namefileinc,status='old',iostat=ierror)
         if (ierror.ne.0) then
            write(*,*) 'bad namefile for arbitrary'
            stop
         endif
         read(15,*) nx,ny,nz
         read(15,*) aretecube
         rewind(15)
         close(15)

         if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
            write(*,*) 'Size of the table too small'
            stop
         endif
         
      endif
c*******************************************************
c     End incident wave
c*******************************************************

      
c*******************************************************
c     Define the object
c*******************************************************
      object='sphere'               !Sphere
c     object='inhomosphere'         !Inhomogeneous sphere
c     object='cube'                 !Cube
c     object='cuboid1'              !Cuboid with side given
c     object='cuboid2'              !Cuboid with nx,ny,nz,aretecube given
c     object='inhomocuboid1'        !Inhomogeneous cuboid with side given
c     object='inhomocuboid2'        !Inhomogeneous cuboid with nx,ny,nz,aretecube given
c     object='ellipsoid'            !Ellipsoid
c     object='nspheres'             !Multiple spheres with same radius
c     object='cylinder'             !Cylinder
c     object='concentricsphere'     !Concentric sphere    
c     object='arbitrary'            !Arbitrary object



      if (object(1:6).eq.'sphere') then
         numberobjet=1
         rayonmulti(1)=200.d0         
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=200.d0
         materiaumulti(1)='xx'
         nnnr=20
         epsmulti(1)=(2.25d0,0.d0)
      elseif (object(1:4).eq.'cube') then
         numberobjet=1
         side=100.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         nnnr=20
         epsmulti(1)=(2.000d0,0.5d0)
      elseif (object(1:7).eq.'cuboid1') then
         numberobjet=1
         sidex=20.d0
         sidey=40.d0
         sidez=10.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         nnnr=10
         phiobj=0.d0
         thetaobj=0.d0
         psiobj=0.d0
         epsmulti(1)=(2.000d0,0.5d0)
      elseif (object(1:7).eq.'cuboid2') then
         numberobjet=1
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         aretecube=50.d0
         epsmulti(1)=(2.000d0,0.5d0)
      elseif (object(1:12).eq.'inhomosphere') then
         numberobjet=1
         hc=0.01d0
         lc=300.d0
         rayonmulti(1)=1000.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         ng=0
         nnnr=20
         epsmulti(1)=(1.01d0,0.0d0)
         materiaumulti(1)='xx'
      elseif (object(1:13).eq.'inhomocuboid1') then
         numberobjet=1
         hc=0.01d0
         lc=300.d0
         sidex=500.d0
         sidey=500.d0
         sidez=500.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         ng=0
         nnnr=20
         epsmulti(1)=(1.01d0,0.0d0)
         materiaumulti(1)='xx'
      elseif (object(1:13).eq.'inhomocuboid2') then
         numberobjet=1
         hc=0.01d0
         lc=300.d0
         aretecube=50.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         ng=0
         nnnr=20
         epsmulti(1)=(1.01d0,0.0d0)
         materiaumulti(1)='xx' 
      elseif (object(1:8).eq.'nspheres') then
         numberobjet=3
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=300.d0
         rayonmulti(1)=1500.d0
         epsmulti(1)=(1.3d0,0.d0)
         materiaumulti(1)='xx'
         xgmulti(2)=2500.d0
         ygmulti(2)=2500.d0
         zgmulti(2)=300.d0
         rayonmulti(2)=400.d0
         epsmulti(2)=(1.1d0,0.d0)
         materiaumulti(2)='xx'
         xgmulti(3)=-2500.d0
         ygmulti(3)=-2500.d0
         zgmulti(3)=300.d0
         rayonmulti(3)=400.d0
         epsmulti(3)=(1.1d0,0.d0)
         materiaumulti(3)='xx'
         nnnr=50
      elseif (object(1:9).eq.'ellipsoid') then
         demiaxea=10.d0
         demiaxeb=20.d0
         demiaxec=40.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         thetaobj=45.d0
         phiobj=45.d0
         psiobj=0.d0
         nnnr=10
         epsmulti(1)=(2.000d0,0.5d0)
      elseif (object(1:8).eq.'cylinder') then
         rayonmulti(1)=10.d0
         hauteur=20.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         thetaobj=0.d0
         phiobj=0.d0
         nnnr=10
         epsmulti(1)=(2.000d0,0.5d0)
      elseif (object(1:16).eq.'concentricsphere') then
         numberobjet=2      
         rayonmulti(1)=10.d0
         epsmulti(1)=(2.25d0,0.d0)
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         rayonmulti(2)=20.d0
         epsmulti(2)=(1.25d0,0.d0)        
         thetaobj=0.d0
         phiobj=0.d0
         nnnr=10
         materiaumulti(1)='xx'
         materiaumulti(2)='xx'
      elseif (object(1:9).eq.'arbitrary') then
         namefileobj='objet.in'
         numberobjet=1
         materiaumulti(1)='xx'
      endif
c*******************************************************
c     End object
c*******************************************************
      
c*******************************************************
c     Define if object is isotropic or not
c*******************************************************
      trope='iso'
c     trope='ani'


c*******************************************************
c     Define the iterative method used
c*******************************************************
      methodeit='GPBICG1'      
c     methodeit='GPBICG2'
c     methodeit='GPBICGplus'
c     methodeit='GPBICGsafe'          
c     methodeit='GPBICGAR'    
c     methodeit='GPBICGAR2'
c     methodeit='BICGSTARPLUS'
c     methodeit='GPBICOR'
c     methodeit='CORS'
c     methodeit='QMRCLA'          
c     methodeit='TFQMR'       
c     methodeit='CG'           
c     methodeit='BICGSTAB'    
c     methodeit='QMRBICGSTAB1' 
c     methodeit='QMRBICGSTAB2' 


c     Define the tolerance of the iterative method
      tolinit=1.d-4
c*******************************************************
c     End iterative method used
c*******************************************************

c*******************************************************
c     Define the multilayer
c*******************************************************
c     example for neps=1
c     epscouche(2)
c     ---------------------------- ZCOUCHE(1)
c     epscouche(1)
c     -----------------------------ZCOUCHE(0)
c     epscouche(0)
c*******************************************************
      neps=0
      epscouche(0)=2.25d0
      epscouche(1)=1.d0
      zcouche(0)=0.d0
      materiaucouche(0)='xx'
      materiaucouche(1)='xx'
c*******************************************************
c     End multilayer
c*******************************************************
      
c*******************************************************
c     define all the options
c*******************************************************
      nobjet=0                  ! 1 compute only the position of the dipole, all the other options are disabled.
      
c     nproche adjust the size of near field domain. Near field (0)
c     inside the object, (1) inside a cuboid which contains the object,
c     (2) inside the boxnx+2*nxmp,ny+2*nymp,nz+2*nzmp
      nproche=2
      
      nxmp=2                    ! if nproche=2 used then the addsize along x : nx+2*nxmp
      nymp=3                    ! if nproche=2 used then the addsize along y : ny+2*nymp
      nzmp=4                    ! if nproche=2 used then the addsize along z : nz+2*nzmp
      nlocal=1                  ! 0 do not compute the local field, 1 compute the local field
      nmacro=1                  ! 0 do not compute the macroscopic field, 1 compute the macroscopic field

c     1 reread or create a file which contains the local field at each
c     position. Avoid to compute again the local field if the
c     configuration is the same, i.e. keep the same local field.
      nlecture=0               
      filereread='toto'         ! name fo the file if reread the local field.

c     nrig adjust the ways used to compute the near field. (0) compute
c     rigorously the near field in solving the near field equation, (1)
c     use renormalized Born approximation, (2) use Born approximation,
c     (3) use Born series at order 1.
      nrig=0

c     ninterp =0 compute rigourosly the Green tensor, ninterp=n compute
c     the Green tensor at the position d/n and interpolate for get the
c     Green tensor at any position.
      ninterp=2      

c     nforce=1 ! calcul la force exercee par la lumiere sur l'objet.
c     nforced=1 ! calcul la densite de force dans l'objet.
c     ntorque=1 ! calcul le couple exerce par la lumiere sur l'objet.
c     ntorqued=1 ! calcul la desnite de couple dans l'objet l'objet.
      
      nsection=0                ! 0 do not compute the cross section, 1 compute the cross section. Possible only in homogeneous space.
      ncote=1                   ! 1 compute both side for the diffracted field and lens, 0 only below (z<0), 2 only above (z>0).
      ndiffracte=1              ! 0 do not compute the diffracted far field, 1 compute the diffracted far field.
      nquickdiffracte=1         ! 0 compute far field classically, 1 compute far field with FFT.      
      nlentille=1               ! Compute microscopy.
      nquicklens=1              ! Compute microscopy with FFT (1) or wihtout FFT (0).      
      numaperref= 1.3d0         ! Numerical aperture for the microscope in reflexion.
      numapertra= 0.9d0         ! Numerical aperture for the microscope in transmission.      
      zlensr=0.d0               ! Position of lens in reflexion.
      zlenst=0.d0               ! Position of lens in transmission.
      ntypemic=0                ! Type of microscope: O Holographic, 1 Bright field, 2 Dark field
      gross=100.d0              ! Manyfing factor for the microscope
      numaperinc=0.9d0          ! Numerical aperture for the condenser lens.

      nenergie=1                ! 0 Do not compute energy, 1 compute energy conservation.

      nmatf=2                   ! 1 Do not save the data, 0 save the data in mat file, 2 save the data in one hdf5 file.
      h5file='ifdda.h5'         ! name of the hdf5 file
     

c*******************************************************
c     End options
c*******************************************************
      

c     compute size when meshsize is given
      if (object(1:13).eq.'inhomocuboid2'.or.object(1:7).eq.'cuboid2')
     $     then
         nx=nxm-2*nxmp
         ny=nym-2*nymp
         nz=nzm-2*nzmp
      else
         if (nx+2*nxmp.gt.nxm) then
            write(*,*) 'pb with size: increase nxm'
            stop
         endif
         if (ny+2*nymp.gt.nym) then
            write(*,*) 'pb with size: increase nym'
            stop
         endif
         if (nz+2*nzmp.gt.nzm) then
            write(*,*) 'pb with size: increase nzm'
            stop
         endif
      endif
      
c*******************************************************
c     Call the routine cdmlibsurf
c*******************************************************
      call cdmlibsurf(
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
     $     kxy,xy,numaperref,numapertra,numaperinc,gross,zlensr,zlenst
     $     ,ntypemic,
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

      if (nstop.eq.0) then
         write(*,*) '***********************************************'
         write(*,*) 'Computation finished without problem:'
         write(*,*) '***********************************************'
      else
         write(*,*) '***********************************************'
         write(*,*) 'Programm finished with problem:'
         write(*,*) infostr
         write(*,*) '***********************************************'
      endif
      


      end
