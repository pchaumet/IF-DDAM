      subroutine passagefourierimagegross(Ex,Ey,Ez,Eimx,Eimy,Eimz,nfft2d
     $     ,nfftmax,imaxk0,deltakx,deltax,gross,k0,indiceopt,planb)
      implicit none
      integer nfft2d,nfft2d2,nfftmax,i,j,kk,indicex,indicey
     $     ,indice,imaxk0
      double precision tmp,deltakx,deltaky,deltax,pi,normal(3) ,gross,u1
     $     ,u2,kx,ky,k0,k0n,indiceopt,fac
      double complex Ex(nfftmax*nfftmax),Ey(nfftmax*nfftmax),Ez(nfftmax
     $     *nfftmax),Eimx(nfftmax*nfftmax),Eimy(nfftmax*nfftmax)
     $     ,Eimz(nfftmax*nfftmax),tmpx,tmpy ,tmpz
      integer FFTW_BACKWARD
      integer*8 planb
      FFTW_BACKWARD=+1


      pi=dacos(-1.d0)
         
      deltaky=deltakx
      deltax=2.d0*pi/dble(nfft2d)/deltakx
      fac=deltakx*deltaky/gross
      nfft2d2=nfft2d/2
      k0n=k0*indiceopt
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)        
      do i=1,nfftmax*nfftmax
         Eimx(i)=0.d0
         Eimy(i)=0.d0
         Eimz(i)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indicex,indicey,indice,kk)   
!$OMP& PRIVATE(kx,ky,normal,u1,u2,tmp)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)          
      do i=-imaxk0,imaxk0 
         do j=-imaxk0,imaxk0
            kx=deltakx*dble(i)
            if (i.ge.0) then
               indicex=i+1
            else
               indicex=nfft2d+i+1
            endif
            ky=deltaky*dble(j)
            if (j.ge.0) then
               indicey=j+1
            else
               indicey=nfft2d+j+1
            endif
            if (kx*kx+ky*ky.le.k0n*k0n) then
               normal(1)=kx/k0n
               normal(2)=ky/k0n
               normal(3)=dsqrt(1.d0-normal(1)*normal(1)-normal(2)
     $              *normal(2))
               
               u1=-normal(2)
               u2=normal(1)
               tmp=dsqrt(u1*u1+u2*u2)

               indice=indicex+nfft2d*(indicey-1)
               kk=i+nfft2d2+1+nfft2d*(j+nfft2d2)

               if (tmp.eq.0.d0) then
                  Eimx(indice)=Ex(kk)
                  Eimy(indice)=Ey(kk)
                  Eimz(indice)=Ez(kk)
               else
                  u1=u1/tmp
                  u2=u2/tmp
                  tmp=dasin(dsin(-dacos(normal(3)))/gross)
     $                 -dacos(normal(3))
                  Eimx(indice)=(u1*u1+(1.d0-u1*u1)*dcos(tmp))*Ex(kk)+u1
     $                 *u2*(1.d0-dcos(tmp))*Ey(kk)+u2*dsin(tmp)*Ez(kk)
                  Eimy(indice)=u1*u2*(1.d0-dcos(tmp))*Ex(kk)+(u2*u2
     $                 +(1.d0-u2*u2)*dcos(tmp))*Ey(kk)-u1*dsin(tmp)
     $                 *Ez(kk)
                  Eimz(indice)=-u2*dsin(tmp)*Ex(kk)+u1*dsin(tmp)*Ey(kk)
     $                 +dcos(tmp)*Ez(kk)
               endif
            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      call dfftw_execute_dft(planb,Eimx,Eimx)
      call dfftw_execute_dft(planb,Eimy,Eimy)
      call dfftw_execute_dft(planb,Eimz,Eimz)
      
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
            
            tmpx=Eimx(kk)*fac
            tmpy=Eimy(kk)*fac
            tmpz=Eimz(kk)*fac
            Eimx(kk)=Eimx(indice)*fac
            Eimy(kk)=Eimy(indice)*fac
            Eimz(kk)=Eimz(indice)*fac
            Eimx(indice)=tmpx
            Eimy(indice)=tmpy
            Eimz(indice)=tmpz
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      end
