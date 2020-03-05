      subroutine fouriertoimage(deltakx,deltaky,gross,Eimagex,Eimagey
     $     ,Eimagez,Eimageincx,Eimageincy,Eimageincz,nfft2D,nfft2d2
     $     ,plan2b,plan2f)
      
      implicit none
      integer nfft2D,nfft2d2
      double precision deltakx,deltaky,gross
      double complex Eimagex(nfft2d*nfft2d),Eimagey(nfft2d*nfft2d)
     $     ,Eimagez(nfft2d*nfft2d),Eimageincx(nfft2d*nfft2d)
     $     ,Eimageincy(nfft2d *nfft2d),Eimageincz(nfft2d*nfft2d)

      integer i,j,indicex,indicey,indice,kk,FFTW_BACKWARD
      double precision tmp
      double complex ctmp
      integer*8 plan2f,plan2b
      FFTW_BACKWARD=+1
      tmp=deltakx*deltaky/gross

#ifdef USE_FFTW
      call dfftw_execute_dft(plan2b,Eimagex,Eimagex)
      call dfftw_execute_dft(plan2b,Eimagey,Eimagey)
      call dfftw_execute_dft(plan2b,Eimagez,Eimagez)
      call dfftw_execute_dft(plan2b,Eimageincx,Eimageincx)
      call dfftw_execute_dft(plan2b,Eimageincy,Eimageincy)
      call dfftw_execute_dft(plan2b,Eimageincz,Eimageincz)
#else 
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      call fftsingletonz2d(Eimageincx,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimageincy,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimageincz,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimagex,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimagey,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimagez,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif


!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,indicex,indicey,indice,kk,ctmp)
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
            
            ctmp=Eimagex(kk)
            Eimagex(kk)=Eimagex(indice)*tmp
            Eimagex(indice)=ctmp*tmp
            ctmp=Eimagey(kk)
            Eimagey(kk)=Eimagey(indice)*tmp
            Eimagey(indice)=ctmp*tmp
            ctmp=Eimagez(kk)
            Eimagez(kk)=Eimagez(indice)*tmp
            Eimagez(indice)=ctmp*tmp

            ctmp=Eimageincx(kk)
            Eimageincx(kk)=Eimageincx(indice)*tmp
            Eimageincx(indice)=ctmp*tmp
            ctmp=Eimageincy(kk)
            Eimageincy(kk)=Eimageincy(indice)*tmp
            Eimageincy(indice)=ctmp*tmp
            ctmp=Eimageincz(kk)
            Eimageincz(kk)=Eimageincz(indice)*tmp
            Eimageincz(indice)=ctmp*tmp
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL


      end
