      subroutine passagefourierimagegross2(Ex,Ey,Ez,nfft2d,nfftmax
     $     ,imaxk0,deltakx,deltax,gross,k0,indiceopt ,planb)
      implicit none
      integer nfft2d,nfft2d2,nfftmax,i,j,kk,indicex,indicey
     $     ,indice,imaxk0
      double precision tmp,deltakx,deltaky,deltax,pi,normal(3) ,gross,u1
     $     ,u2,kx,ky,k0,k0n,indiceopt,fac
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
!$OMP& PRIVATE(kx,ky,normal,u1,u2,tmp,tmpx,tmpy,tmpz)
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

               if (tmp.ne.0.d0) then
                  u1=u1/tmp
                  u2=u2/tmp
                  tmp=dasin(dsin(-dacos(normal(3)))/gross)
     $                 -dacos(normal(3))

                  tmpx=Ex(indice)
                  tmpy=Ey(indice)
                  tmpz=Ez(indice)
                              
                  Ex(indice)=(u1*u1+(1.d0-u1*u1)*dcos(tmp))*tmpx+u1*u2
     $                 *(1.d0-dcos(tmp))*tmpy+u2*dsin(tmp)*tmpz
                  Ey(indice)=u1*u2*(1.d0-dcos(tmp))*tmpx+(u2*u2+(1.d0
     $                 -u2*u2)*dcos(tmp))*tmpy-u1*dsin(tmp)*tmpz
                  Ez(indice)=-u2*dsin(tmp)*tmpx+u1*dsin(tmp)*tmpy
     $                 +dcos(tmp)*tmpz
                  
              
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
