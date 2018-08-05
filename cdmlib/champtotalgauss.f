c     Author : l'equipe triso
c     
c     routine qui calcule le champ lointain total reflechi (Ediffkzneg)
c     et transmis (Ediffkzpos) qui correspond à E(r).

c     Attention ici on calcule E(r)=Int[e(k||) exp(i k||.r|| +ikz z) ]
c     rappel : phase stationnaire E(r)=-2*icomp*pi*kz*e(k||)*exp(ik0r)/r
c     Ici on travaille avec e(k||)

      subroutine champtotalgauss(imax,ncote,nfft2d,nfft2dmax,nepsmax
     $     ,neps,epscouche,k0,aretecube,NA,Ediffkzpos,Ediffkzneg
     $     ,Egausxref ,Egausyref ,Egauszref ,Egausxtra,Egausytra
     $     ,Egausztra,fluxref,fluxtrans)
      implicit none
      integer i,j,ii,jj,nfft2d,nfft2d2,nfft2dmax,ncote ,imax,indice
      double precision kx,ky,kz,kp2,fluxref,fluxtrans,tmp,deltakx
     $     ,deltaky,aretecube,k0,NA
      double complex Ediffkzpos(nfft2dmax,nfft2dmax,3)
     $     ,Ediffkzneg(nfft2dmax,nfft2dmax,3),Egausxref(nfft2dmax
     $     *nfft2dmax),Egausyref(nfft2dmax *nfft2dmax)
     $     ,Egauszref(nfft2dmax*nfft2dmax),Egausxtra(nfft2dmax
     $     *nfft2dmax),Egausytra(nfft2dmax*nfft2dmax)
     $     ,Egausztra(nfft2dmax*nfft2dmax)
      integer nepsmax,neps
      double precision indice0,indicen,pi
      double complex epscouche(0:nepsmax+1),icomp,const


c     initialisation
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      indice0=dsqrt(dreal(epscouche(0)))
      indicen=dsqrt(dreal(epscouche(neps+1)))
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      nfft2d2=nfft2d/2  
      fluxtrans=0.d0
      fluxref=0.d0

c      Ediffkzpos=0.d0
c      Ediffkzneg=0.d0
c      Egausxtra=0.d0
c      Egausytra=0.d0
c      Egausztra=0.d0
c     calcul du champ total
      if (ncote.eq.0.or.ncote.eq.1) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kp2,kz) 
!$OMP& PRIVATE(indice,const)
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
                  const=(-2.d0*icomp*pi*kz)
                  
c                  write(*,*) 'ttt',Ediffkzpos(ii,jj,1),Egausxtra(indice)
c     $                 ,Ediffkzpos(ii,jj,2),Egausytra(indice)
c     $                 ,Ediffkzpos(ii,jj,3),Egausztra(indice),indice
                  
                  Egausxtra(indice)=Egausxtra(indice)+Ediffkzpos(ii,jj
     $                 ,1)/const
                  Egausytra(indice)=Egausytra(indice)+Ediffkzpos(ii,jj
     $                 ,2)/const
                  Egausztra(indice)=Egausztra(indice)+Ediffkzpos(ii,jj
     $                 ,3)/const

c                  write(*,*) 'rrr',Ediffkzpos(ii,jj,1),Egausxtra(indice)
c     $                 ,Ediffkzpos(ii,jj,2),Egausytra(indice)
c     $                 ,Ediffkzpos(ii,jj,3),Egausztra(indice),indice

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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kp2,kz) 
!$OMP& PRIVATE(indice,const)
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

                  indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                  const=(-2.d0*icomp*pi*kz)

c                  write(*,*) 'rrr',Ediffkzneg(ii,jj,1),Egausxref(indice)
c     $                 ,Ediffkzneg(ii,jj,2),Egausyref(indice)
c     $                 ,Ediffkzneg(ii,jj,3),Egauszref(indice),indice
                  
                  Egausxref(indice)=Egausxref(indice)+Ediffkzneg(ii,jj
     $                 ,1)/const
                  Egausyref(indice)=Egausyref(indice)+Ediffkzneg(ii,jj
     $                 ,2)/const
                  Egauszref(indice)=Egauszref(indice)+Ediffkzneg(ii,jj
     $                 ,3)/const
                
            

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
  
   
      tmp=4.d0*pi*pi*deltaky*deltakx/(k0*8.d0*pi*1.d-7*299792458.d0)
     
      fluxref=fluxref*tmp
      fluxtrans=fluxtrans*tmp
      write(*,*) 'nnn',fluxref,fluxtrans,'irra',tmp
      return
      end
