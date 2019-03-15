      use HDF5
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
      parameter (nxm=20,nym=20,nzm=20,nphi=72,ntheta=36)
c     definition of the size for the code
      INTEGER nmax, ntotalm
      character(2) polarizability
c     variables for the positions
      double precision rayon,side,sidex,sidey ,sidez,hauteur
     $     ,xgmulti(numberobjetmax) ,ygmulti(numberobjetmax)
     $     ,zgmulti(numberobjetmax) ,rayonmulti(numberobjetmax),demiaxea
     $     ,demiaxeb,demiaxec ,thetaobj,phiobj,psiobj,lc,hc,density
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
      double precision ss,pp,theta,phi,I0,phim(20)
     $     ,thetam(20),ssm(20),ppm(20)
      double complex Eloc(3),Em(3),E0,uncomp,icomp,zzero,E0m(20)
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
      parameter (nfft2d=128)
      double precision kx,ky,kz,deltakx,deltaky,numaper,deltax,gross
     $     ,numaperinc,zlensr,zlenst
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
      
      
c     input
      lambda=632.8d0
      P0=1.d0            
      w0=lambda*10.d0    
      c=299792458.d0
      quatpieps0=1.d0/(c*c*1.d-7)
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
c********************************************************
c     Defini la polarisabilite choisie
c********************************************************

c      polarizability='CM'
      polarizability='RR'
c      methode='GB' pas fait
c      methode='LA' pas fait
c      methode='LR' pas fait


      
c*******************************************************
c     Defini l'onde incidente
c*******************************************************
c      beam='gwavelinear'
c      beam='gwavelinear'
c     beam='gwavecircular'
      beam='pwavelinear'
c      beam='pwavecircular'
c      beam='wavelinearmulti'
c      beam='gwaveiso'
c     beam='speckle'
c     beam='arbitrary' !pas fait
      if (beam(1:11).eq.'pwavelinear') then
c         theta=65.287098590488455d0
         theta=0.d0
         phi=0.d0
         pp=0.d0
         ss=1.d0
         write(*,*) 'totot'
      elseif (beam(1:13).eq.'pwavecircular') then
         theta=0.d0
         phi=0.d0
         ss=1.d0 
      elseif (beam(1:11).eq.'gwavelinear') then
         theta=0.d0
         phi=0.d0
         pp=1.d0
         ss=0.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
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
c     Defini l'objet
c*******************************************************
      object='sphere'
c      object='cube'
c      object='cuboid'
c      object='nspheres'     
c      object='ellipsoid'       
c      object='cylinder'           
c      object='concentricsphere'         
c      object='arbitrary'
c      namefileobj='arbitrary.in'
c      object='inhomo'

      if (object(1:6).eq.'sphere') then
         numberobjet=1
         rayonmulti(1)=500.d0         
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         nnnr=10
         epsmulti(1)=(1.1d0,0.d0)
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
         write(*,*) 'cuboid1'
      elseif (object(1:6).eq.'inhomo') then
         numberobjet=1
         hc=0.0d0
         lc=300.d0
         rayonmulti(1)=1000.d0
         ng=0
         nnnr=20
         epsmulti(1)=(1.01d0,0.0d0)
         materiaumulti(1)='xx'
      elseif (object(1:8).eq.'nspheres') then
         numberobjet=2
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         rayonmulti(1)=10.d0
         epsmulti(1)=(2.25d0,0.d0)
         materiaumulti(1)='xx'
         xgmulti(2)=30.d0
         ygmulti(2)=0.d0
         zgmulti(2)=0.d0
         rayonmulti(2)=15.d0
         epsmulti(2)=(1.25d0,0.d0)
         materiaumulti(2)='xx'
         nnnr=20
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
c     objet isotrope ou non
c*******************************************************
      trope='iso'
c     trope='ani'


c*******************************************************
c     defini la methode iterative utilisee
c*******************************************************
      methodeit='GPBICG1'      
c      methodeit='GPBICG2'      
c      methodeit='GPBICGAR1'    
c      methodeit='GPBICGAR2'
c      methodeit='GPBICGsafe'    
c      methodeit='QMRCLA'          
c      methodeit='TFQMR'       
c      methodeit='CG'           
c      methodeit='BICGSTAB'    
c      methodeit='QMRBICGSTAB1' 
c      methodeit='QMRBICGSTAB2' 
c      methodeit='GPBICOR' 

      tolinit=1.d-4

c*******************************************************
c     defini ce que l'on veut calculer
c*******************************************************
      nproche=1 !0 calcul le champ dans l'objet, 1 dans le cube
                !contenant l'objet 2 dans la boite nxm,nym,nzm
      nxmp=0 ! si boite plus large que l'objet suivant x
      nymp=0 ! si boite plus large que l'objet suivant y
      nzmp=0 ! si boite plus large que l'objet suivant z
      nlocal=1                  ! 0 ne calcul pas le champ local, 1 calcul le champ local
      nmacro=1 ! 0 ne calcul pas le champ macro, 1 calcul le champ macro
      nsection=1
      ndiffracte=1 !0  ne calcul pas le champ diffracte, 1 le calcul
      nquickdiffracte=1 ! 0 calcul le champ diffracte classique, 1 par FFT
      nrig=0 ! 0 calcul rigoureusement le champ, 1 Born renormalise
      ninterp=0 ! niveau d'interpolation si tenseur non rigoureux
      ncote=1 ! 1 calcul des deux cotes le champ diffracté, 0 dessous, 2 dessus
      numaper= 1.d0 ! ouverture numerique, 1 si on calcul des deux cotes
c     nforce=1 ! calcul la force exercee par la lumiere sur l'objet
c     nforced=1 ! calcul la densite de force dans l'objet
c     ntorque=1 ! calcul le couple exerce par la lumiere sur l'objet
c     ntorqued=1 ! calcul la desnite de couple dans l'objet l'objet
      nlentille=1 !calcul l'image a travers une lentille: 0 pas de lentille, 1 au dessus, -1 au dessous (on élcaire par en dessous)
      nquicklens=1!image a travers la lentille sans FFT (0) ou avec FFT(1)
      nenergie=0 !0: calcul du champ diffracte. 1 calcul du champ total
      nlecture=0 ! 0 ne relit pas, 1 relit ou cree le fichier FF
! (dipole) au premier coup
      nobjet=0                  ! si 1 ne fait que l'objet
      nmatf=0 ! 1 ne cre pas les fichiers .mat 0 sinon
      h5file='ifdda.h5'
      gross=100.d0              !grossissement
      numaperinc=0.9d0          !ouverture numerique de l'incident
      zlensr=0.d0               ! position lentille reflexion
      zlenst=0.d0               ! position lentille transmission
      ntypemic=0                !type microscope
c*******************************************************
c     Definition du multicouhe
c
c      ---------------------------

c     ----------------------------
c     epscouche(1)
c     -----------------------------ZCOUCHE(0)
c     epscouche(0)
   
          
c*******************************************************
      neps=0
c      epscouche(0)=1.d0
c      epscouche(1)=1.d0
      epscouche(0)=1.d0
      epscouche(1)=1.d0
c      epscouche(2)=1.d0
c      epscouche(2)=1.d0 
      zcouche(0)=-rayonmulti(1)*1.5d0
      zcouche(1)=rayonmulti(1)*1.5d0
     
      materiaucouche(0)='xx'
      materiaucouche(1)='xx'
      materiaucouche(2)='xx'
      write(*,*) 'data',lambda,beam,object,trope, materiaumulti,nnnr
     $     ,tolinit,methodeit,zcouche,epscouche
      write(*,*) 'numberobjet',numberobjet
c*******************************************************
c     fin de la definiton des parametres d'entres
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
     $     nxm, nym, nzm, nxmp, nymp, nzmp,
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
     $     kxy,xy,numaper,numaperinc,gross,zlensr,zlenst,ntypemic,
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
c     output
      write(*,*) 'nstop',nstop,infostr

      write(*,*) 'numberobjet',numberobjet

      stop
      
      if (nstop.eq.1) then
         write(*,*) infostr
         stop
      endif

      pi=dacos(-1.d0)
      if (materiau.ne.'xx') then     
         write(*,*) 'Relative permittivity',epsmulti(1)
      else
         if (trope.eq.'iso') then          
            write(*,*) 'Relative permittivity',epsmulti(1)
         else 
            do i=1,3
               do j=1,3
                  write(*,*) 'Relative permittivity',epsanimulti(1,i,j)
     $                 ,i,j
               enddo
            enddo
         endif
      endif

      write(*,*) 'efficacite',efficacite,efficaciteref,efficacitetrans

      write(*,*) 'Object under study ',object
      write(*,*) 'Nombre Object',numberobjet
      write(*,*) 'number of subunit for the object',nbsphere
      write(*,*) 'number of subunit for the mesh ',ndipole
      write(*,*) 'mesh size',aretecube
      write(*,*) 'lambda/(10n)',lambda/10.d0/cdabs(cdsqrt(epsmulti(1)))

      write(*,*) '******* Compute the incident field *******'
      write(*,*) 'Beam used',beam
      write(*,*) 'k0=',k0     
      write(*,*) 'theta=',theta
      write(*,*) 'phi=',phi
      write(*,*) 'Irradiance',Irra
      write(*,*) 'Field',E0
      I0=cdabs(E0)**2

      write(*,*) '***** Solve the linear system *****'
      write(*,*) 'Tolerance asked for the iterative method   ',tolinit
      write(*,*) 'Tolerance obtained for the iterative method',tol1
      write(*,*) 'Number of product Ax for the iterative method'
     $     ,ncompte,nloop
      nsection=1
      if (nsection.eq.1) then      
         rayon=rayonmulti(1)*1.d-9
         write(*,*) 'mie',epsmulti(1),rayon,lambda
       
         CALL CALLBHMIE(dreal(cdsqrt(epscouche(0))),epsmulti(1),rayon
     $        ,lambda,MIECEXT,MIECABS,MIECSCA,GSCA)
         write(*,*) 'MIECEXT',MIECEXT,MIECABS,MIECSCA,GSCA
c         I0=1.d0
c         epsmulti(1)=2.d0
c         lambda=lambda/dsqrt(2.d0)
c         CALL CALLBHMIE(I0,epsmulti(1),rayon,lambda,MIECEXT,MIECABS
c     $        ,MIECSCA,GSCA)
c         write(*,*) 'MIECEXT',MIECEXT,MIECABS,MIECSCA,GSCA

         write(*,*) 'force',(MIECEXT-GSCA*MIECSCA)/8.d0/pi*quatpieps0

         write(*,*) 'extinction cross section',Cext
         write(*,*) 'absorbing cross section ',Cabs
         write(*,*) 'scattering cross section',Csca
         write(*,*) 'cos',gasym
      endif
      if (ndiffracte.eq.1) then
         write(*,*) 'scattering cross section with integration',Cscai
         write(*,*) 'scattering asymetric parameter',gasym
      endif

      if (nforce.eq.1) then
         write(*,*) '****** Compute the optical force *********'
         write(*,*) 'optical force x',forcet(1)
         write(*,*) 'optical force y',forcet(2)
         write(*,*) 'optical force z',forcet(3)
         
         forcemie=(MIECext-GSCA*MIECsca)/8.d0/pi*I0*quatpieps0
         write(*,*) 'modulus of the force',forcem,'Mie',forcemie
      endif
      if (ntorque*nforce.eq.1) then
         write(*,*) '********* Compute the optical torque *********'
         write(*,*) 'optical torque x',couplet(1)
         write(*,*) 'optical torque y',couplet(2)
         write(*,*) 'optical torque z',couplet(3)
         write(*,*) 'modulus of the optical torque',couplem
         write(*,*) 'couple Mie',MIECABS/8.d0/k0/pi*I0*quatpieps0
      endif
c      if (numberobjet.ne.1) then
c         do i=1,numberobjet
c            write(*,*) 'numero',i
c            write(*,*) 'forcex',forcexmulti(i)
c            write(*,*) 'forcey',forceymulti(i)
c            write(*,*) 'forcez',forcezmulti(i)
c         enddo
c      endif


      end
