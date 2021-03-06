c     Cette routine calcul le champ image à partir du cgamp de Fourier
c     pour un grossessment de -1. Attention le champ de Fourier est
c     rentré dans E(x,y,z) et le champ image est sorti dans le même
c     tableau.

      subroutine passagefourierimage2(Ex,Ey,Ez,nfft2d,nfftmax,imaxk0
     $     ,indiceopt,deltakx,deltax,planb)
      implicit none
      integer nfft2d,nfft2d2,nfftmax,i,j,kk,indicex,indicey,indeicez
     $     ,indice,imaxk0
      double precision fac,deltakx,deltaky,deltax,pi,indiceopt
     $     ,indiceoptrac
      double complex Ex(nfftmax*nfftmax),Ey(nfftmax*nfftmax),Ez(nfftmax
     $     *nfftmax),tmpx,tmpy ,tmpz
      integer FFTW_BACKWARD
      integer*8 planb
      FFTW_BACKWARD=+1
      pi=dacos(-1.d0)
      
      indiceoptrac=dsqrt(indiceopt)
      deltaky=deltakx
      deltax=2.d0*pi/dble(nfft2d)/deltakx
      fac=deltakx*deltaky

      nfft2d2=nfft2d/2

#ifdef USE_FFTW
      call dfftw_execute_dft(planb,Ex,Ex)
      call dfftw_execute_dft(planb,Ey,Ey)
      call dfftw_execute_dft(planb,Ez,Ez)
#else   
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
         call fftsingletonz2d(Ex,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
         call fftsingletonz2d(Ey,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION  
         call fftsingletonz2d(Ez,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif


      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indicex,indicey,indice,kk)   
!$OMP& PRIVATE(tmpx,tmpy,tmpz)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)     
      do i=-nfft2d2,nfft2d2-1
         do j=-nfft2d2,-1

            if (i.ge.0) then
               indicex=i+1
            else
               indicex=nfft2d+i+1
            endif
            if (j.ge.0) then
               indicey=j+1
            else
               indicey=nfft2d+j+1
            endif

            indice=indicex+nfft2d*(indicey-1)
            kk=i+nfft2d2+1+nfft2d*(j+nfft2d2)
            
            tmpx=Ex(kk)*fac
            tmpy=Ey(kk)*fac
            tmpz=Ez(kk)*fac
            Ex(kk)=Ex(indice)*fac
            Ey(kk)=Ey(indice)*fac
            Ez(kk)=Ez(indice)*fac
            Ex(indice)=tmpx*indiceoptrac
            Ey(indice)=tmpy*indiceoptrac
            Ez(indice)=tmpz*indiceoptrac
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      end
