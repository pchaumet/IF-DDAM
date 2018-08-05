      subroutine objetspheresurf(trope,eps,epsani,xs,ys,zs,xswf,yswf
     $     ,zswf,k0 ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny
     $     ,nz,methode ,na,epsilon,polarisa,rayon,xg,yg,zg,neps,nepsmax
     $     ,dcouche,zcouche,epscouche,tabzn,infostr,nstop)

      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,ii,jj,i,j,k
     $     ,test,IP(3),nnnr,dddis,inv,na,nstop
      double precision xs(nmax),ys(nmax),zs(nmax),xswf(nmax),yswf(nmax)
     $     ,zswf(nmax),k0,xg,yg,zg,x,y,z,aretecube
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),ctmp,eps0

      integer neps,nepsmax,nminc,nmaxc,numerocouche,tabzn(nmax)
      double precision dcouche(nepsmax),zcouche(0:nepsmax),zmin,zmax
      double complex epscouche(0:nepsmax+1)

      double precision rayon
      character(2) methode
      character(3) trope
      character(64) infostr
      write (*,*) 'pola0 = ',methode
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

      rayon=rayon*1.d-9
      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      if (rayon.le.1.d-12) then
         nstop=1
         infostr='object sphere: radius=0!'
         return
      endif
      write (*,*) 'na = ',na
      write (*,*) 'nnnr = ',nnnr
      write (*,*) 'pola = ',methode
      write(*,*) 'trop',trope
c     verfie si on est bien multiple de 2 3 5 pour la discretisation,
c     car ma FFT est basee sur une decomposition en nombre premier de 2
c     3 5. Si on utilisait FFTW ceci disparaitrait.

      nx=nnnr
      ny=nnnr
      nz=nnnr
c     size of the subunit
      aretecube=2.d0*rayon/dble(nnnr)

c     test si l'objet est sur une couche ou plusieurs
      zmax=zg+rayon
      zmin=zg-rayon
      nminc=numerocouche(zmin,neps,nepsmax,zcouche)
      nmaxc=numerocouche(zmax,neps,nepsmax,zcouche)
      write(*,*) 'obj',zmin,zmax,nminc,nmaxc,xg,yg,zg,zcouche
      if (nmaxc-nminc.ge.2) then
c     shift the layers
         do k=nminc,nmaxc-1           
            do i=1,nnnr
               z=-rayon+aretecube*(dble(i)-0.5d0)+zg
               if (zcouche(k).ge.z.and.zcouche(k).lt.z+aretecube) then
                  zcouche(k)=z+aretecube/2.d0
               endif              
            enddo
         enddo
      elseif (nmaxc-nminc.eq.1) then
c     object inside two layers
c     compute the position of down subunit: shift the object
         z=-rayon+aretecube*0.5d0
c         write(*,*) 'z',z,z/aretecube,zg,aretecube,z/aretecube+zg
c     $        /aretecube+0.5d0,nint(z/aretecube+zg /aretecube+0.5d0)
c     $        ,aint(z/aretecube+zg /aretecube+0.5d0)
         zg=aretecube*nint(z/aretecube+zg/aretecube+0.5d0)-z-aretecube
     $        /2.d0
c         write(*,*) z/aretecube+zg/aretecube+0.5d0,idint(z/aretecube+zg
c     $        /aretecube+0.5d0)
      endif

      write(*,*) 'obj2',nminc,nmaxc,zg
      if (na.eq.-1) then
         do i=1,nnnr 
            do j=1,nnnr              
               do k=1,nnnr                
                  x=-rayon+aretecube*(dble(k)-0.5d0)
                  y=-rayon+aretecube*(dble(j)-0.5d0)
                  z=-rayon+aretecube*(dble(i)-0.5d0)

                  if (j.eq.1.and.k.eq.1) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1) write(21,*) y+yg
                  if (j.eq.1.and.i.eq.1) write(20,*) x+xg

                  ndipole=ndipole+1
                  xswf(ndipole)=x+xg
                  yswf(ndipole)=y+yg
                  zswf(ndipole)=z+zg
                  Tabzn(ndipole)=i
                  if (dsqrt(x*x+y*y+z*z).lt.rayon) then
                     nbsphere=nbsphere+1
                     Tabdip(ndipole)=nbsphere                    
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
                     eps0=epscouche(numerocouche(zs(nbsphere),neps
     $                    ,nepsmax,zcouche))
                     if (trope.eq.'iso') then
                        call poladiffcomp(aretecube,eps,eps0,k0,dddis
     $                       ,methode,ctmp)  
                        polarisa(nbsphere,1,1)=ctmp
                        polarisa(nbsphere,2,2)=ctmp
                        polarisa(nbsphere,3,3)=ctmp
                        epsilon(nbsphere,1,1)=eps
                        epsilon(nbsphere,2,2)=eps
                        epsilon(nbsphere,3,3)=eps
                     else                        
                        call polaepstenscomp(aretecube,epsani,eps0,k0
     $                       ,dddis,methode,inv,polaeps)
                        do ii=1,3
                           do jj=1,3
                              epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                              polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           enddo
                        enddo
                     endif
                     write(10,*) xs(nbsphere)
                     write(11,*) ys(nbsphere)
                     write(12,*) zs(nbsphere)
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

                  if (j.eq.1.and.k.eq.1) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1) write(21,*) y+yg
                  if (j.eq.1.and.i.eq.1) write(20,*) x+xg

                  ndipole=ndipole+1
                  nbsphere=nbsphere+1
                  Tabdip(ndipole)=nbsphere
                  Tabzn(ndipole)=i
                  xswf(ndipole)=x+xg
                  yswf(ndipole)=y+yg
                  zswf(ndipole)=z+zg
                  if (dsqrt(x*x+y*y+z*z).lt.rayon) then                    
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
                     eps0=epscouche(numerocouche(zs(nbsphere),neps
     $                    ,nepsmax,zcouche))
                     if (trope.eq.'iso') then                       
                        call poladiffcomp(aretecube,eps,eps0,k0,dddis
     $                       ,methode,ctmp)
c                        write(*,*) 'polarefsphere',aretecube,eps,eps0,k0
c     $                       ,dddis,methode,ctmp
c                        write (*,*) 'nbsphere = ', nbsphere
                        polarisa(nbsphere,1,1)=ctmp
                        polarisa(nbsphere,2,2)=ctmp
                        polarisa(nbsphere,3,3)=ctmp
                        epsilon(nbsphere,1,1)=eps
                        epsilon(nbsphere,2,2)=eps
                        epsilon(nbsphere,3,3)=eps
                     else
                        call polaepstenscomp(aretecube,epsani,eps0,k0
     $                       ,dddis,methode,inv,polaeps)
                        do ii=1,3
                           do jj=1,3
                              epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                              polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           enddo
                        enddo
                     endif
                     write(10,*) xs(nbsphere)
                     write(11,*) ys(nbsphere)
                     write(12,*) zs(nbsphere)
                  else
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg   
                     eps0=epscouche(numerocouche(zs(nbsphere),neps
     $                    ,nepsmax,zcouche))
                     epsilon(nbsphere,1,1)=eps0
                     epsilon(nbsphere,2,2)=eps0
                     epsilon(nbsphere,3,3)=eps0
                     write(10,*) xs(nbsphere)
                     write(11,*) ys(nbsphere)
                     write(12,*) zs(nbsphere) 
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
      write (*,*) ' OBJECT SPHERE FINISHED'

      end
