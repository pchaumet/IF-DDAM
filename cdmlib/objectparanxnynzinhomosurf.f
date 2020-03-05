      subroutine objetparanxnynzinhomosurf(eps,xs,ys,zs,xswf,yswf,zswf
     $     ,k0,aretecube ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm
     $     ,nym,nzm,nxmp,nymp,nzmp,methode ,na ,epsilon,polarisa,sidex
     $     ,sidey,sidez ,xg ,yg,zg,lc ,hc,ng,epsb,neps ,nepsmax ,dcouche
     $     ,zcouche ,epscouche,tabzn ,nmatf,file_id,group_iddip,infostr
     $     ,nstop)

#ifdef USE_HDF5
      use HDF5
#endif

      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,nxm ,nym ,nzm
     $     ,nxmp,nymp,nzmp,i,j,k,test ,IP(3),nnnr,dddis,inv,na,nstop
     $     ,nmatf
      double precision xs(nmax),ys(nmax),zs(nmax),xswf(nmax),yswf(nmax)
     $     ,zswf(nmax),k0,lc,hc,x,y,z,aretecube,pi
      double complex eps,polarisa(nmax,3,3) ,epsilon(nmax,3 ,3),ctmp
     $     ,epsb(nmax),icomp,eps0

      integer nk,ngraine,ng
      double precision kx,ky,kz,dkx,dky,dkz,coeff1,coeff2,coeff3 ,phase
     $     ,spec,moyenne,ecartype,lx,ly,lz,sunif,sidex,sidey,sidez
     $     ,xg ,yg,zg

      integer neps,nepsmax,nminc,nmaxc,numerocouche,tabzn(nmax)
      double precision dcouche(nepsmax),zcouche(0:nepsmax),zmin,zmax
      double complex epscouche(0:nepsmax+1)

      
      character(2) methode
      character(64) infostr
      integer *8 planb
      integer FFTW_BACKWARD,FFTW_ESTIMATE,FFTW_FORWARD
#ifndef USE_HDF5
      integer,parameter:: hid_t=4
#endif
      character(LEN=100) :: datasetname
      integer(hid_t) :: file_id
      integer(hid_t) :: group_iddip
      integer :: dim(4)
      integer error
      
      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64
      
c     Initialization
      nbsphere=0
      ndipole=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
      do i=1,nmax
         Tabdip(i)=0
         polarisa(i,1,1)=0.d0
         polarisa(i,1,2)=0.d0
         polarisa(i,1,3)=0.d0
         polarisa(i,2,1)=0.d0
         polarisa(i,2,2)=0.d0
         polarisa(i,2,3)=0.d0
         polarisa(i,3,1)=0.d0
         polarisa(i,3,2)=0.d0
         polarisa(i,3,3)=0.d0
         epsilon(i,1,1)=0.d0
         epsilon(i,1,2)=0.d0
         epsilon(i,1,3)=0.d0
         epsilon(i,2,1)=0.d0
         epsilon(i,2,2)=0.d0
         epsilon(i,2,3)=0.d0
         epsilon(i,3,1)=0.d0
         epsilon(i,3,2)=0.d0
         epsilon(i,3,3)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL   

      dddis=1
      inv=1
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      ngraine=4*ng+1

      nx=nxm-2*nxmp
      ny=nym-2*nymp
      nz=nzm-2*nzmp
      aretecube=aretecube*1.d-9
      sidex=aretecube*dble(nx)
      sidey=aretecube*dble(ny)
      sidez=aretecube*dble(nz)
      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      lc=lc*1.d-9
      
c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

      if (sidex.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidex=0!'
         return
      elseif (sidey.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidey=0!'
         return
      elseif (sidez.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidez=0!'
         return
      endif
      if (lc.le.0.d0) then
         infostr='coherence length equal to zero'
         nstop=1
         return
      endif      
      if (hc.le.0.d0) then
         infostr='standard deviation equal to zero'
         nstop=1
         return
      endif
      
c      write (*,*) 'na = ',na
c      write (*,*) 'nnnr = ',nnnr
c      write (*,*) 'pola = ',methode
      
c     verfie si on est bien multiple de 2 3 5 pour la discretisation,
c     car ma FFT est basee sur une decomposition en nombre premier de 2
c     3 5. Si on utilisait FFTW ceci disparaitrait.

      if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
         nstop=1
         infostr='Dimension Problem of the Box : Box too small!'
         write(99,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         write(*,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         return
      endif

      
c     test si l'objet est sur une couche ou plusieurs
      zmax=zg+sidez/2.d0
      zmin=zg-sidez/2.d0
      nminc=numerocouche(zmin,neps,nepsmax,zcouche)
      nmaxc=numerocouche(zmax,neps,nepsmax,zcouche)
c      write(*,*) 'obj',zmin,zmax,nminc,nmaxc,zcouche
      if (nmaxc-nminc.ge.1) then
c     shift the layers
         do k=nminc,nmaxc-1           
            do i=1,nz
               z=-sidez/2.d0+aretecube*(dble(i)-0.5d0)+zg
               if (zcouche(k).ge.z.and.zcouche(k).lt.z+aretecube) then
                  zcouche(k)=z+aretecube/2.d0
               endif              
            enddo
         enddo
      endif

c      write(*,*) 'obj2',nminc,nmaxc,zg
      lx=dble(nx)*aretecube
      ly=dble(ny)*aretecube
      lz=dble(nz)*aretecube
      dkx=2.d0*pi/lx
      dky=2.d0*pi/ly
      dkz=2.d0*pi/lz
      nk=0
      
      coeff1=0.125d0*hc*hc*lc*lc*lc/(pi*dsqrt(pi))
      coeff2=lc*lc
      coeff3=4.d0*dkx*dky*dkz
      do i=1,nz
         do j=1,ny
            do k=1,nx
               nk=nk+1
               if (k.ge.1 .and. k.le.(nx/2+1)) then
                  kx=dble(k-1)*dkx
               else
                  kx=dble(k-1-nx)*dkx
               endif
               if (j.ge.1 .and. j.le.(ny/2+1)) then
                  ky=dble(j-1)*dky
               else
                  ky=dble(j-1-ny)*dky
               endif
               if (i.ge.1 .and. i.le.(nz/2+1)) then
                  kz=dble(i-1)*dkz
               else
                  kz=dble(i-1-nz)*dkz
               endif
c     calcul du spec
               spec=coeff1*dexp(-0.25d0*coeff2*(kx*kx+ky*ky+kz*kz))
c     phase aleatoire
               phase=2.d0*pi*SUNIF(ngraine)
c     Amplitude complexe
               epsb(nk)=dsqrt(spec*coeff3)*cdexp(icomp *phase)
c     Spectre symetrique pour diffraction harmonique
               if (nk.eq.1) epsb(nk)=0.d0
               if (nk.ge.nx/2+2.and.nk.le.nx) epsb(nk)=0.d0
               if (nk.ge.(ny/2+1)*nx+1.and.nk.le.nx*ny) epsb(nk)=0.d0
               if (nk.ge.(nz/2+1)*nx*ny+1.and.nk.le.nx*ny*nz) epsb(nk)
     $              =0.d0

            enddo
         enddo
      enddo     
      
c     Profil des hauteurs
#ifdef USE_FFTW
      call dfftw_plan_dft_3d(planb, nx,ny,nz,epsb,epsb,FFTW_BACKWARD
     $     ,FFTW_ESTIMATE)
      call dfftw_execute_dft(planb, epsb, epsb)
#else  
      call fftsingletonz3d(epsb,NX,NY,NZ,FFTW_BACKWARD)
#endif
      moyenne=0.d0
      ecartype=0.d0
      do i=1,nk
         epsb(i)=dreal(epsb(i))+eps
         moyenne=moyenne+dreal(epsb(i))
         ecartype=ecartype+cdabs(epsb(i))**2.d0      
      enddo
        
      
      moyenne=moyenne/dble(nk)
      ecartype=ecartype/dble(nk)
    
      do i=1,nz
         do j=1,ny
            do k=1,nx
               
               x=-sidex/2.d0+aretecube*(dble(k)-0.5d0)
               y=-sidey/2.d0+aretecube*(dble(j)-0.5d0)
               z=-sidez/2.d0+aretecube*(dble(i)-0.5d0)                 
               
               if (j.eq.1.and.k.eq.1.and.nmatf.eq.0) write(22,*) z+zg
               if (i.eq.1.and.k.eq.1.and.nmatf.eq.0) write(21,*) y+yg
               if (j.eq.1.and.i.eq.1.and.nmatf.eq.0) write(20,*) x+xg
               
               ndipole=ndipole+1
               xswf(ndipole)=x+xg
               yswf(ndipole)=y+yg
               zswf(ndipole)=z+zg
               Tabzn(ndipole)=i
               
               nbsphere=nbsphere+1
               Tabdip(ndipole)=nbsphere
               xs(nbsphere)=x+xg
               ys(nbsphere)=y+yg
               zs(nbsphere)=z+zg
c     write(*,*) 'ggg',xs(nbsphere),ys(nbsphere)
c     $                    ,zs(nbsphere),xswf(ndipole),yswf(ndipole)
c     $                    ,zswf(ndipole),ndipole,nbsphere,zg
               eps=epsb(ndipole)
               eps0=epscouche(numerocouche(zs(nbsphere),neps ,nepsmax
     $              ,zcouche))
               call poladiffcomp(aretecube,eps,eps0,k0,dddis
     $              ,methode,ctmp)  
               polarisa(nbsphere,1,1)=ctmp
               polarisa(nbsphere,2,2)=ctmp
               polarisa(nbsphere,3,3)=ctmp
               epsilon(nbsphere,1,1)=eps
               epsilon(nbsphere,2,2)=eps
               epsilon(nbsphere,3,3)=eps
            enddo
         enddo
      enddo

      if (ndipole.gt.nmax) then
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
      endif
      if (nmatf.eq.0) then
         do i=1,nbsphere
            write(10,*) xs(i)
            write(11,*) ys(i)
            write(12,*) zs(i)
         enddo
      elseif (nmatf.eq.2) then
         
         dim(1)=nbsphere
         dim(2)=nmax
         datasetname="Dipole position x"
         call hdf5write1d(group_iddip,datasetname,xs,dim)
         datasetname="Dipole position y"
         call hdf5write1d(group_iddip,datasetname,ys,dim)
         datasetname="Dipole position z"
         call hdf5write1d(group_iddip,datasetname,zs,dim)
      endif

      
      close(10)
      close(11)
      close(12)
      close(20)
      close(21)
      close(22)
      write (*,*) ' OBJECT SPHERE INHOMOGENEOUS FINISHED'
      write (*,*) ' SEED',ng,'Average',moyenne,'Standard deviation'
     $     ,dsqrt(ecartype-moyenne*moyenne)

      write(99,*) ' SEED',ng
      write(99,*) 'Average',moyenne
      write(99,*) 'Standard deviation',dsqrt(ecartype-moyenne*moyenne)

      end
