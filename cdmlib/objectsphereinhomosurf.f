      subroutine objetsphereinhomosurf(eps,xs,ys,zs,xswf,yswf,zswf ,k0
     $     ,aretecube ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz
     $     ,methode ,na ,epsilon,polarisa,rayon,lc,hc,ng,epsb,neps
     $     ,nepsmax ,dcouche,zcouche,epscouche,tabzn,nùatf,infostr
     $     ,nstop)

      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,i,j,k,test
     $     ,IP(3),nnnr,dddis,inv,na,nstop,nmatf
      double precision xs(nmax),ys(nmax),zs(nmax),xswf(nmax),yswf(nmax)
     $     ,zswf(nmax),k0,lc,hc,x,y,z,zg,aretecube,pi
      double complex eps,polarisa(nmax,3,3) ,epsilon(nmax,3 ,3),ctmp
     $     ,epsb(nmax),icomp,eps0

      integer nk,ngraine,ng
      double precision rayon,kx,ky,kz,dkx,dky,dkz,coeff1,coeff2,coeff3
     $     ,phase,spec,moyenne,ecartype,lx,ly,lz,sunif

      integer neps,nepsmax,nminc,nmaxc,numerocouche,tabzn(nmax)
      double precision dcouche(nepsmax),zcouche(0:nepsmax),zmin,zmax
      double complex epscouche(0:nepsmax+1)

      
      character(2) methode
      character(64) infostr
      integer *8 planb
      integer FFTW_BACKWARD,FFTW_ESTIMATE,FFTW_FORWARD

      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64
      
      write(*,*) 'inhomogenous sphere',eps
c     Initialization
      nbsphere=0
      ndipole=0 
      Tabdip=0
      polarisa=0.d0
      epsilon=0.d0
      dddis=1
      inv=1
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      ngraine=4*ng+1

      lc=lc*1.d-9
      
c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

      rayon=rayon*1.d-9

      if (rayon.le.1.d-12) then
         nstop=1
         infostr='object sphere: radius=0!'
         return
      endif
      write (*,*) 'na = ',na
      write (*,*) 'nnnr = ',nnnr
      write (*,*) 'pola = ',methode
      
c     verfie si on est bien multiple de 2 3 5 pour la discretisation,
c     car ma FFT est basee sur une decomposition en nombre premier de 2
c     3 5. Si on utilisait FFTW ceci disparaitrait.

      nx=nnnr
      ny=nnnr
      nz=nnnr
c     size of the subunit
      aretecube=2.d0*rayon/dble(nnnr)

c     test si l'objet est sur une couche ou plusieurs
      zmax=rayon
      zmin=-rayon
      nminc=numerocouche(zmin,neps,nepsmax,zcouche)
      nmaxc=numerocouche(zmax,neps,nepsmax,zcouche)
      write(*,*) 'obj',zmin,zmax,nminc,nmaxc,zcouche
      if (nmaxc-nminc.ge.2) then
c     shift the layers
         zg=0.d0
         do k=nminc,nmaxc-1           
            do i=1,nnnr
               z=-rayon+aretecube*(dble(i)-0.5d0)
               if (zcouche(k).ge.z.and.zcouche(k).lt.z+aretecube) then
                  zcouche(k)=z+aretecube/2.d0
               endif              
            enddo
         enddo
      elseif (nmaxc-nminc.eq.1) then
c     object inside two layers
c     compute the position of down subunit: shift the object
         z=-rayon+aretecube*0.5d0
         zg=aretecube*nint(z/aretecube+0.5d0)-z-aretecube/2.d0
      endif

      write(*,*) 'obj2',nminc,nmaxc,zg
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
      call dfftw_plan_dft_3d(planb, nx,ny,nz,epsb,epsb,FFTW_BACKWARD
     $     ,FFTW_ESTIMATE)

      moyenne=0.d0
      ecartype=0.d0
      do i=1,nk
         epsb(i)=dreal(epsb(i))+eps
         moyenne=moyenne+dreal(epsb(i))
         ecartype=ecartype+cdabs(epsb(i))**2.d0      
      enddo
        
      
      moyenne=moyenne/dble(nk)
      ecartype=ecartype/dble(nk)
    
      write(*,*) 'moyenne',moyenne,ecartype
      if (na.eq.-1) then
         do i=1,nnnr         
            do j=1,nnnr              
               do k=1,nnnr                
                  x=-rayon+aretecube*(dble(k)-0.5d0)
                  y=-rayon+aretecube*(dble(j)-0.5d0)
                  z=-rayon+aretecube*(dble(i)-0.5d0)

                  if (j.eq.1.and.k.eq.1.and.nmatf.eq.0) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1.and.nmatf.eq.0) write(21,*) y
                  if (j.eq.1.and.i.eq.1.and.nmatf.eq.0) write(20,*) x

                  ndipole=ndipole+1
                  xswf(ndipole)=x
                  yswf(ndipole)=y
                  zswf(ndipole)=z+zg
                  Tabzn(ndipole)=i
                  
                  if (dsqrt(x*x+y*y+z*z).lt.rayon) then
                     nbsphere=nbsphere+1
                     Tabdip(ndipole)=nbsphere
                     xs(nbsphere)=x
                     ys(nbsphere)=y
                     zs(nbsphere)=z+zg
c                     write(*,*) 'ggg',xs(nbsphere),ys(nbsphere)
c     $                    ,zs(nbsphere),xswf(ndipole),yswf(ndipole)
c     $                    ,zswf(ndipole),ndipole,nbsphere,zg
                     eps=epsb(ndipole)
                     eps0=epscouche(numerocouche(zs(nbsphere),neps
     $                    ,nepsmax,zcouche))
                     call poladiffcomp(aretecube,eps,eps0,k0,dddis
     $                    ,methode,ctmp)  
                     polarisa(nbsphere,1,1)=ctmp
                     polarisa(nbsphere,2,2)=ctmp
                     polarisa(nbsphere,3,3)=ctmp
                     epsilon(nbsphere,1,1)=eps
                     epsilon(nbsphere,2,2)=eps
                     epsilon(nbsphere,3,3)=eps
                     
                     if (nmatf.eq.0) write(10,*) xs(nbsphere)
                     if (nmatf.eq.0) write(11,*) ys(nbsphere)
                     if (nmatf.eq.0) write(12,*) zs(nbsphere)
                  endif
               enddo
            enddo
         enddo
      elseif (na.eq.0) then
         do i=1,nnnr         
            do j=1,nnnr              
               do k=1,nnnr                
                  x=-rayon+aretecube*(dble(k)-0.5d0)
                  y=-rayon+aretecube*(dble(j)-0.5d0)
                  z=-rayon+aretecube*(dble(i)-0.5d0)

                  if (j.eq.1.and.k.eq.1.and.nmatf.eq.0) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1.and.nmatf.eq.0) write(21,*) y
                  if (j.eq.1.and.i.eq.1.and.nmatf.eq.0) write(20,*) x

                  ndipole=ndipole+1
                  nbsphere=nbsphere+1
                  Tabdip(ndipole)=nbsphere
                  Tabzn(ndipole)=i
                  xswf(ndipole)=x
                  yswf(ndipole)=y
                  zswf(ndipole)=z+zg
                  
                  if (dsqrt(x*x+y*y+z*z).lt.rayon) then                    
                     xs(nbsphere)=x
                     ys(nbsphere)=y
                     zs(nbsphere)=z+zg
c                     write(*,*) 'fff',xs(nbsphere),ys(nbsphere)
c     $                    ,zs(nbsphere),xswf(ndipole),yswf(ndipole)
c     $                    ,zswf(ndipole),ndipole,nbsphere,zg
                     eps=epsb(ndipole)
                     eps0=epscouche(numerocouche(zs(nbsphere),neps
     $                    ,nepsmax,zcouche))
                     call poladiffcomp(aretecube,eps,eps0,k0,dddis
     $                    ,methode,ctmp)
c                     write(*,*) 'polaref',aretecube,eps,eps0,k0,dddis
c     $                    ,methode,ctmp
c     write (*,*) 'nbsphere = ', nbsphere
                     polarisa(nbsphere,1,1)=ctmp
                     polarisa(nbsphere,2,2)=ctmp
                     polarisa(nbsphere,3,3)=ctmp
                     epsilon(nbsphere,1,1)=eps
                     epsilon(nbsphere,2,2)=eps
                     epsilon(nbsphere,3,3)=eps
                     
                     if (nmatf.eq.0) write(10,*) xs(nbsphere)
                     if (nmatf.eq.0) write(11,*) ys(nbsphere)
                     if (nmatf.eq.0) write(12,*) zs(nbsphere)
                  else
                     xs(nbsphere)=x
                     ys(nbsphere)=y
                     zs(nbsphere)=z
                     eps0=epscouche(numerocouche(zs(nbsphere),neps
     $                    ,nepsmax,zcouche))
                     epsilon(nbsphere,1,1)=eps0
                     epsilon(nbsphere,2,2)=eps0
                     epsilon(nbsphere,3,3)=eps0
                     if (nmatf.eq.0) write(10,*) xs(nbsphere)
                     if (nmatf.eq.0) write(11,*) ys(nbsphere)
                     if (nmatf.eq.0) write(12,*) zs(nbsphere) 
                  endif
               enddo
            enddo
         enddo
      else
         infostr='na should be equal to -1 or 0'
         nstop=1
         return
      endif
      if (ndipole.gt.nmax) then
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
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
