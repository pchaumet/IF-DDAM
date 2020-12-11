      subroutine dipoleinctotal(imax,ncote,nfft2d,nfft2dmax,nepsmax
     $     ,neps,epscouche,k0,aretecube,NA,Ediffkzpos,Ediffkzneg
     $     ,zcouche,dcouche,E0,thetat,phit,xdip,ydip,zdip,Egausxref
     $     ,Egausyref ,Egauszref ,Egausxtra,Egausytra ,Egausztra,fluxref
     $     ,fluxtrans ,fluxinc ,nstop ,infostr)
      implicit none
      integer i,j,ii,jj,nfft2d,nfft2d2,nfft2dmax,ncote ,imax,indice
     $     ,nstop
      double precision kx,ky,kz,kkz,kp2,fluxref,fluxtrans,tmp,deltakx
     $     ,deltaky,aretecube,k0,NA,kxinc,kyinc,fluxinc,x,y,z,xdip,ydip
     $     ,zdip
      double complex Ediffkzpos(nfft2dmax,nfft2dmax,3)
     $     ,Ediffkzneg(nfft2dmax,nfft2dmax,3),Egausxref(nfft2dmax
     $     *nfft2dmax),Egausyref(nfft2dmax *nfft2dmax)
     $     ,Egauszref(nfft2dmax*nfft2dmax),Egausxtra(nfft2dmax
     $     *nfft2dmax),Egausytra(nfft2dmax*nfft2dmax)
     $     ,Egausztra(nfft2dmax*nfft2dmax)
      integer nepsmax,neps,imini,jmini
      double precision indice0,indicen,pi,ss,pp,theta,phi,thetat,phit,zz
     $     ,dcouche(nepsmax),zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),E0,icomp ,const,p(3),Ex,Ey
     $     ,Ez,Stenseur(3,3),ctmp
      character(64) infostr

c     initialisation
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      indice0=dsqrt(dreal(epscouche(0)))
      indicen=dsqrt(dreal(epscouche(neps+1)))
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      nfft2d2=nfft2d/2  
      
      x=0.d0
      y=0.d0
      write(*,*) 'dipole incident'
      write(*,*) 'delta k',deltakx,'m-1'
      write(*,*) 'Number of point in the numerical aperture',imax*2+1
      write(*,*) 'Size of FFT',nfft2d
      write(*,*) 'Size of the pixel',aretecube

      fluxtrans=0.d0
      fluxref=0.d0

      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0
      p(1)=E0*dsin(theta)*dcos(phi)
      p(2)=E0*dsin(theta)*dsin(phi)
      p(3)=E0*dcos(theta)
      write(*,*) 'dipole',p
      zz=1.d0
      
c     calcul du champ total
      if (ncote.eq.0.or.ncote.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
         do i=1,nfft2d*nfft2d
            Egausxtra(i)=0.d0
            Egausytra(i)=0.d0
            Egausztra(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kp2,kz) 
!$OMP& PRIVATE(indice,const,Stenseur,ctmp,Ex,Ey,Ez)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:fluxtrans) COLLAPSE(2)
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               kp2=kx*kx+ky*ky
               if (indicen*indicen*k0*k0*NA*NA-kp2.gt.0.d0)
     $              then
                  kz=dsqrt(indicen*indicen*k0*k0-kp2)
                  indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)

                  call tenseurmulticoucheloinfft(kx,ky,kz,zz,zdip,k0
     $                 ,nepsmax,neps,dcouche,zcouche,epscouche
     $                 ,Stenseur)

                  Ex=Stenseur(1,1)*p(1)+Stenseur(1,2)*p(2)+Stenseur(1,3)
     $                 *p(3)
                  Ey=Stenseur(2,1)*p(1)+Stenseur(2,2)*p(2)+Stenseur(2,3)
     $                 *p(3)
                  Ez=Stenseur(3,1)*p(1)+Stenseur(3,2)*p(2)+Stenseur(3,3)
     $                 *p(3)
                  ctmp=cdexp(-icomp*(kx*xdip+ky*ydip))
                  const=(-2.d0*icomp*pi*kz)
                  
                  Egausxtra(indice)=(Ediffkzpos(ii,jj,1)+Ex*ctmp)/const
                  Egausytra(indice)=(Ediffkzpos(ii,jj,2)+Ey*ctmp)/const
                  Egausztra(indice)=(Ediffkzpos(ii,jj,3)+Ez*ctmp)/const

                  fluxtrans=fluxtrans+(cdabs(Egausxtra(indice))**2
     $                 +cdabs(Egausytra(indice))**2
     $                 +cdabs(Egausztra(indice))**2)*kz
c                  write(*,*) 'rrrrr',kzt,fluxtrans
c                  write(*,*) 'champ1',Egausxtra(indice)
c     $                 ,Egausytra(indice),Egausztra(indice),i,j,indice
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
      endif
      if (ncote.eq.0.or.ncote.eq.-1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
         do i=1,nfft2d*nfft2d
            Egausxref(i)=0.d0
            Egausyref(i)=0.d0
            Egauszref(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kp2,kz,kkz) 
!$OMP& PRIVATE(indice,const,Stenseur,ctmp,Ex,Ey,Ez)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:fluxref) COLLAPSE(2)
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               kp2=kx*kx+ky*ky
               if (indice0*indice0*k0*k0*NA*NA-kp2.gt.0.d0)
     $              then
                  kz=dsqrt(indice0*indice0*k0*k0-kp2)
                  kkz=-kz
                  call tenseurmulticoucheloinfft(kx,ky,kkz,-zz,zdip,k0
     $                 ,nepsmax,neps,dcouche,zcouche,epscouche
     $                 ,Stenseur)

                  Ex=Stenseur(1,1)*p(1)+Stenseur(1,2)*p(2)+Stenseur(1,3)
     $                 *p(3)
                  Ey=Stenseur(2,1)*p(1)+Stenseur(2,2)*p(2)+Stenseur(2,3)
     $                 *p(3)
                  Ez=Stenseur(3,1)*p(1)+Stenseur(3,2)*p(2)+Stenseur(3,3)
     $                 *p(3)
                  ctmp=cdexp(-icomp*(kx*xdip+ky*ydip))

                  
                  indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                  const=(-2.d0*icomp*pi*kz)

                  Egausxref(indice)=(Ediffkzneg(ii,jj,1)+Ex*ctmp)/const
                  Egausyref(indice)=(Ediffkzneg(ii,jj,2)+Ey*ctmp)/const
                  Egauszref(indice)=(Ediffkzneg(ii,jj,3)+Ez*ctmp)/const


                  fluxref=fluxref+(cdabs(Egausxref(indice))**2
     $                 +cdabs(Egausyref(indice))**2
     $                 +cdabs(Egauszref(indice))**2)*kz
c                  write(*,*) 'fffrrr' ,fluxref,kz
c                  write(*,*) 'champ2',Egausxref(indice)
c     $                 ,Egausyref(indice),Egauszref(indice),i,j
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      endif

c      write(*,*) 'E0',E0

      tmp=4.d0*pi*pi*deltaky*deltakx/(k0*8.d0*pi*1.d-7*299792458.d0)
      fluxref=fluxref*tmp
      fluxtrans=fluxtrans*tmp
      return
      end

      
      
      end
