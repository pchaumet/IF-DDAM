      subroutine objetcubesurf(trope,eps,epsani,eps0,xs,ys,zs,xswf,yswf
     $     ,zswf,k0 ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny
     $     ,nz,methode ,epsilon,polarisa,side,xg,yg,zg,neps,nepsmax
     $     ,dcouche,zcouche ,epscouche,tabzn,nmatf,infostr,nstop)
      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,ii,jj,i,j,k
     $     ,test,IP(3),nnnr,dddis,inv,nnmax,nstop,nmatf
      double precision xs(nmax),ys(nmax),zs(nmax),xswf(nmax),yswf(nmax)
     $     ,zswf(nmax),k0,xg,yg,zg,x,y,z,aretecube,side,xc,yc,zc
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),ctmp,eps0


      integer neps,nepsmax,nminc,nmaxc,numerocouche,tabzn(nmax)
      double precision dcouche(nepsmax),zcouche(0:nepsmax),zmin,zmax
      double complex epscouche(0:nepsmax+1)


      character(2) methode
      character(3) trope
      character(64) infostr

c     Initialization
      nbsphere=0
      ndipole=0 
      Tabdip=0
      Tabzn=0
      polarisa=0.d0
      epsilon=0.d0
      dddis=1
      inv=1
      
c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

      side=side*1.d-9
      write(*,*) 'objetcube::side=',side
      if (side.eq.0) then
         write(*,*) 'objetcube::side is equal to zero: Aborting',side
         nstop=1
         infostr='object cube: side=0!'
         return
      endif

      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      
c     verfie si on est bien multiple de 2 3 5 pour la discretisation,
c     car ma FFT est basee sur une decomposition en nombre premier de 2
c     3 5. Si on utilisait FFTW ceci disparaitrait.

      nx=nnnr
      ny=nnnr
      nz=nnnr
c     size of the subunit
      aretecube=side/dble(nnnr)

c     test si l'objet est sur une couche ou plusieurs
      zmax=zg+side/2.d0
      zmin=zg-side/2.d0
      nminc=numerocouche(zmin,neps,nepsmax,zcouche)
      nmaxc=numerocouche(zmax,neps,nepsmax,zcouche)

      if (nmaxc-nminc.ge.2) then
c     shift the layers
         do k=nminc,nmaxc-1           
            do i=1,nnnr
               z=dble(i)*aretecube-(side+aretecube)/2.d0+zg
               if (zcouche(k).ge.z.and.zcouche(k).lt.z+aretecube) then
                  zcouche(k)=z+aretecube/2.d0
               endif              
            enddo
         enddo
      elseif (nmaxc-nminc.eq.1) then
c     object inside two layers
c     compute the position of down subunit
         z=-side/2.d0+aretecube*0.5d0
         zg=aretecube*nint(z/aretecube+zg/aretecube+0.5d0)-z-aretecube
     $        /2.d0
      endif
      xc=(side+aretecube)/2.d0
      yc=(side+aretecube)/2.d0
      zc=(side+aretecube)/2.d0

      do i=1,nnnr         
         do j=1,nnnr              
            do k=1,nnnr                
               x=dble(k)*aretecube-xc
               y=dble(j)*aretecube-yc
               z=dble(i)*aretecube-zc

               if (j.eq.1.and.k.eq.1.and.nmatf.eq.0) write(22,*) z+zg
               if (i.eq.1.and.k.eq.1.and.nmatf.eq.0) write(21,*) y+yg
               if (j.eq.1.and.i.eq.1.and.nmatf.eq.0) write(20,*) x+xg

               ndipole=ndipole+1               
               nbsphere=nbsphere+1
               Tabdip(ndipole)=nbsphere
               Tabzn(nbsphere)=i
               xs(nbsphere)=x+xg
               ys(nbsphere)=y+yg
               zs(nbsphere)=z+zg
               xswf(ndipole)=x+xg
               yswf(ndipole)=y+yg
               zswf(ndipole)=z+zg
               eps0=epscouche(numerocouche(zs(nbsphere),neps ,nepsmax
     $              ,zcouche))
               if (trope.eq.'iso') then
                  call poladiffcomp(aretecube,eps,eps0,k0,dddis,methode
     $                 ,ctmp)
                
                  polarisa(nbsphere,1,1)=ctmp
                  polarisa(nbsphere,2,2)=ctmp
                  polarisa(nbsphere,3,3)=ctmp
                  epsilon(nbsphere,1,1)=eps
                  epsilon(nbsphere,2,2)=eps
                  epsilon(nbsphere,3,3)=eps
               else
                  call polaepstenscomp(aretecube,epsani,eps0,k0,dddis
     $                 ,methode,inv,polaeps)
                  do ii=1,3
                     do jj=1,3
                        epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                        polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                     enddo
                  enddo
               endif
               if (nmatf.eq.0) write(10,*) xs(nbsphere)
               if (nmatf.eq.0) write(11,*) ys(nbsphere)
               if (nmatf.eq.0) write(12,*) zs(nbsphere)
               
            enddo
         enddo
      enddo
      if (ndipole.gt.nmax) then
         write(*,*) 'Size of the parameter too small',ndipole
     $        ,'should be smaller that',nmax
         write(*,*) 'Increase nxm nym and nzm'
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
      endif

      close(10)
      close(11)
      close(12)
      close(15)
      close(20)
      close(21)
      close(22)

      end
