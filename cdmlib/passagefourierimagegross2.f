      subroutine passagefourierimagegross2(Ex,Ey,Ez,nfft2d,nfftmax
     $     ,imaxk0,deltakx,deltax,gross,k0,indiceopt,side,planb)
      implicit none
      integer nfft2d,nfft2d2,nfftmax,i,j,kk,indicex,indicey
     $     ,indice,imaxk0
      double precision tmp,deltakx,deltaky,deltax,pi,gross,u(3)
     $     ,v(3),kx,ky,k0,k0n,indiceopt,fac,side,u1,u2,costmp,sintmp
      double complex Ex(nfftmax*nfftmax),Ey(nfftmax*nfftmax),Ez(nfftmax
     $     *nfftmax),tmpx,tmpy ,tmpz
      integer FFTW_BACKWARD
      integer*8 planb
      FFTW_BACKWARD=+1


      pi=dacos(-1.d0)
         
      deltaky=deltakx
      deltax=2.d0*pi/dble(nfft2d)/deltakx
      fac=deltakx*deltaky/gross
      nfft2d2=nfft2d/2
      k0n=k0*indiceopt


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indicex,indicey,indice,kk)   
!$OMP& PRIVATE(kx,ky,u,v,u1,u2,tmp,costmp,sintmp,tmpx,tmpy,tmpz)
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
             
               u(1)=kx/k0n
               u(2)=ky/k0n
               u(3)=dsqrt(1.d0-u(1)*u(1)-u(2)*u(2))*side

               v(1)=-kx/k0n/gross
               v(2)=-ky/k0n/gross
               v(3)=dsqrt(1.d0-v(1)*v(1)-v(2)*v(2))*side
               indice=indicex+nfft2d*(indicey-1)
               u1=u(2)*v(3)-u(3)*v(2)
               u2=-u(1)*v(3)+u(3)*v(1)
               costmp=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
               sintmp=dsqrt(u1*u1+u2*u2)
               u1=u1/sintmp
               u2=u2/sintmp

               if (sintmp.ne.0.d0) then

                  tmpx=Ex(indice)
                  tmpy=Ey(indice)
                  tmpz=Ez(indice)
                              
                  Ex(indice)=(u1*u1+(1.d0-u1*u1)*costmp)*tmpx+u1*u2 *(1
     $                 .d0-costmp)*tmpy+u2*sintmp*tmpz
                  Ey(indice)=u1*u2*(1.d0-costmp)*tmpx+(u2*u2+(1.d0 -u2
     $                 *u2)*costmp)*tmpy-u1*sintmp*tmpz
                  Ez(indice)=-u2*sintmp*tmpx+u1*sintmp*tmpy+costmp*tmpz
                  
              
               endif
            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL


      call dfftw_execute_dft(planb,Ex,Ex)
      call dfftw_execute_dft(planb,Ey,Ey)
      call dfftw_execute_dft(planb,Ez,Ez)
      
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
            Ex(indice)=tmpx
            Ey(indice)=tmpy
            Ez(indice)=tmpz
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      end
