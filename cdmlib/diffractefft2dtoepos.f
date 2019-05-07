      subroutine diffractefft2dtoepos(Ediffkzpos,Ediffkzneg,Efourierxpos
     $     ,Efourierypos,Efourierzpos,Efourierxneg,Efourieryneg
     $     ,Efourierzneg,epscouche,nepsmax,neps,numaperref ,numapertra
     $     ,k0,aretecube,imax ,nfft2d,nfft2dmax,ncote,nstop,infostr)
      implicit none
      integer i,j,ii,jj,nfft2dmax,nfft2d,nfft2d2,imax,ncote,indice
     $     ,nepsmax,neps,nstop
      double complex epscouche(0:nepsmax+1)
      double precision numaperref ,numapertra,k0,deltakx,deltaky,pi,kx
     $     ,ky,kz,indice0,indicen,aretecube
      double complex Ediffkzpos(nfft2dmax,nfft2dmax,3)
     $     ,Ediffkzneg(nfft2dmax,nfft2dmax,3),Efourierxpos(nfft2dmax
     $     *nfft2dmax),Efourierypos(nfft2dmax*nfft2dmax)
     $     ,Efourierzpos(nfft2dmax*nfft2dmax),Efourierxneg(nfft2dmax
     $     *nfft2dmax),Efourieryneg(nfft2dmax*nfft2dmax)
     $     ,Efourierzneg(nfft2dmax*nfft2dmax),ctmp,icomp

c     Info string
      character(64) infostr
      
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      nfft2d2=nfft2d/2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
      do i=1,nfft2d*nfft2d
         Efourierxpos(i)=0.d0
         Efourierypos(i)=0.d0
         Efourierzpos(i)=0.d0
         Efourierxneg(i)=0.d0
         Efourieryneg(i)=0.d0
         Efourierzneg(i)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
       if (nfft2d.gt.4096) then
         nstop=1
         infostr='nfft2d too large'
         return
      endif
      
      
      if (deltakx.ge.k0) then
         nstop=1
         infostr='In FFT  diffracted to epos nfft2d too small'
         return
      endif

      
      if (ncote.eq.0.or.ncote.eq.1) then        
         indicen=dsqrt(dreal(epscouche(neps+1)))
      endif

      if (ncote.eq.0.or.ncote.eq.-1) then            
         indice0=dsqrt(dreal(epscouche(0)))
      endif
      
  
      if (ncote.eq.0.or.ncote.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kz,indice,ctmp)   
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)          
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               if (indicen*indicen*k0*k0*numapertra*numapertra*0.9999d0
     $              -kx*kx-ky*ky.gt.0.d0) then
                  
                  kz=dsqrt(indicen*indicen*k0*k0-kx*kx-ky*ky) 
                 
                  indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                  ctmp=-2.d0*pi*icomp*kz
                  Efourierxpos(indice)=Ediffkzpos(ii,jj,1)/ctmp
                  Efourierypos(indice)=Ediffkzpos(ii,jj,2)/ctmp
                  Efourierzpos(indice)=Ediffkzpos(ii,jj,3)/ctmp   
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL         
      endif


      if (ncote.eq.0.or.ncote.eq.-1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kz,indice,ctmp)   
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)       
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               if (indice0*indice0*k0*k0*numaperref*numaperref*0.9999d0
     $              -kx*kx-ky*ky.gt.0.d0) then
                  kz=dsqrt(indice0*indice0*k0*k0-kx*kx-ky*ky) 
                  indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                  ctmp=-2.d0*pi*icomp*kz
                  Efourierxneg(indice)=Ediffkzneg(ii,jj,1)/ctmp
                  Efourieryneg(indice)=Ediffkzneg(ii,jj,2)/ctmp
                  Efourierzneg(indice)=Ediffkzneg(ii,jj,3)/ctmp   

               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL           
      endif

      
      end
