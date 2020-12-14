      subroutine diffractefft2dsurf2lens(nbsphere,nx,ny,nz,nxm,nym,nzm
     $     ,nfft2d,nfft2dmax,tabfft2,k0,xs,ys,zs,aretecube,Eloinx,Eloiny
     $     ,Eloinz ,FF,imax,deltakx,deltaky,Ediffkzpos,NA,nepsmax ,neps
     $     ,dcouche ,zcouche ,epscouche,signe,nstop ,infostr,planf)
      implicit none
      integer nx,ny,nz,nxm,nym,nzm,nfft2dmax,nfft2d,nstop
      double precision xs(nxm*nym*nzm),ys(nxm*nym*nzm),zs(nxm*nym*nzm)
     $     ,aretecube,k0,zz,za,NA,signe
      double complex FF(3*nxm*nym*nzm),Ediffkzpos(nfft2dmax,nfft2dmax,3)

      integer nfft2d2,imax,i,j,k,tabfft2(nfft2dmax),indice,kk,ii,jj
      double precision deltakx,deltaky,var1,var2,kx,ky,kz,fac,pi
      double complex ctmp,ctmp1,icomp,Eloinx(nfft2dmax*nfft2dmax)
     $     ,Eloiny(nfft2dmax*nfft2dmax),Eloinz(nfft2dmax*nfft2dmax)
     $     ,Stenseur(3,3)

      integer nepsmax,neps
      double precision dcouche(nepsmax),zcouche(0:nepsmax),indice0
     $     ,indicen,indicem
      double complex epscouche(0:nepsmax+1)

      integer nbsphere,iii,iiii,jjjj
      double precision x,y,z
      double complex Em(3)

c     Info string
      character(64) infostr
      integer FFTW_FORWARD
      integer*8 planf

      FFTW_FORWARD=-1
      zz=signe
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      nfft2d2=nfft2d/2  
      var1=(xs(1)+dble(nfft2d2)*aretecube)*deltakx
      var2=(ys(1)+dble(nfft2d2)*aretecube)*deltaky

      if (signe.ge.0.d0) then
         if (dreal(epscouche(neps+1)).le.0.d0) then
            infostr='epsilon <0 higher layer computation lens'
            nstop=-1
            return
         endif
          if (dimag(epscouche(neps+1)).gt.0.d0) then
            infostr='absorbing layer to computate lens'
            nstop=-1
            return
         endif
         indicen=dsqrt(dreal(epscouche(neps+1)))
         indicem=indicen
      endif

      if (signe.le.0.d0) then     
         if (dreal(epscouche(0)).le.0.d0) then
            infostr='epsilon <0 higher layer computation lens'
            nstop=-1
            return
         endif
          if (dimag(epscouche(0)).gt.0.d0) then
            infostr='absorbing layer to computate lens'
            nstop=-1
            return
         endif
         indice0=dsqrt(dreal(epscouche(0)))
         indicem=indice0
      endif
      
      if (nfft2d.gt.65536) then
         nstop=-1
         infostr='FFT for the diffracted field too large'
         return
      endif

      do i=1,nfft2d
         if (i-nfft2d2-1.ge.0) then
            tabfft2(i)=i-nfft2d2
         else
            tabfft2(i)=nfft2d2+i
         endif
      enddo

      write(*,*) NA  
      write(*,*) k0
      write(*,*) deltakx
      write(*,*) indicen
      write(*,*) NA*k0*indicem/deltakx

c     imax=nint(NA*k0*indicem/deltakx)+1
      imax=nint(NA*k0/deltakx*max(dsqrt(dreal(epscouche(0)))
     $     ,dsqrt(dreal(epscouche(neps+1)))))+1
c      imax=nint(NA*k0*indice0/deltakx)+1
      write(*,*) 'imax',imax,nint(NA*k0*indice0/deltakx)+1,indicem
      write(*,*) 'pp',NA*k0*indice0/deltakx,NA,k0,indice0,deltakx
      if (2*imax+1.gt.nfft2d) then
         infostr='Size of FFT too small to compute the diffracted field'
         nstop=-1
         return
      endif
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)   
!$OMP DO SCHEDULE(STATIC)  COLLAPSE(2) 
         do i=1,nfft2d
            do j=1,nfft2d
               Ediffkzpos(i,j,1)=0.d0
               Ediffkzpos(i,j,2)=0.d0
               Ediffkzpos(i,j,3)=0.d0
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  

      do k=1,nz
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nfft2d*nfft2d
            Eloinx(i)=0.d0
            Eloiny(i)=0.d0
            Eloinz(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kk,indice)   
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2) 
         do i=1,nx
            do j=1,ny
               kk=i+nx*(j-1)+nx*ny*(k-1)
               indice=tabfft2(i)+nfft2d*(tabfft2(j)-1)
               Eloinx(indice)=FF(3*(kk-1)+1)                
               Eloiny(indice)=FF(3*(kk-1)+2)    
               Eloinz(indice)=FF(3*(kk-1)+3)  
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  

#ifdef USE_FFTW
         call dfftw_execute_dft(planf,Eloinx,Eloinx)
         call dfftw_execute_dft(planf,Eloiny,Eloiny)
         call dfftw_execute_dft(planf,Eloinz,Eloinz)  
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
         call fftsingletonz2d(Eloinx,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
         call fftsingletonz2d(Eloiny,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION  
         call fftsingletonz2d(Eloinz,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif         

         kk=1+nx*ny*(k-1)

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP& PRIVATE(i,j,ii,jj,kx,ky,kz,indice,ctmp1,Stenseur)   
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2) 
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               if (NA*NA*indicem*indicem*k0*k0*0.9999d0-kx*kx-ky
     $              *ky.gt.0.d0) then
c     write(*,*) 'ij',i,j
                  kz=signe*dsqrt(indicem*indicem*k0*k0-kx*kx-ky*ky)
                  indice=tabfft2(i+nfft2d2+1)+nfft2d*(tabfft2(j
     $                 +nfft2d2+1)-1)
c     rappel : E(r)=Int [ Ediffkzpos(k||)exp(ik||.r||+ikz z)] dk||
c     en champ lointain : E(r)=-2*pi*icomp*kz*e(k||) *exp(ik0r)/r
c
c 
                  ctmp1=cdexp(-icomp*(var1*dble(i)+var2*dble(j)))/(-2.d0
     $                 *pi*icomp*dabs(kz))
c                  ctmp1=cdexp(-icomp*(var1*dble(i)+var2*dble(j)))


                  call tenseurmulticoucheloinfft(kx,ky,kz,zz,zs(kk)
     $                 ,k0,nepsmax,neps,dcouche,zcouche,epscouche
     $                 ,Stenseur)
c     E(r)=A(k||)exp(ik0r)/r= Int (Stenseur(k||,r')*P(r') dr' )*exp(ik0r)/r
c     
c    donc  :  Ediffkzpos(k||) = A(k||)/(-2*icomp*pi*kz)
c

                  Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)+(Stenseur(1
     $                 ,1)*Eloinx(indice)+Stenseur(1,2)
     $                 *Eloiny(indice)+Stenseur(1,3)*Eloinz(indice))
     $                 *ctmp1

                  Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)+(Stenseur(2
     $                 ,1)*Eloinx(indice)+Stenseur(2,2)
     $                 *Eloiny(indice)+Stenseur(2,3)*Eloinz(indice))
     $                 *ctmp1
                  
                  Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)+(Stenseur(3
     $                 ,1)*Eloinx(indice)+Stenseur(3,2)
     $                 *Eloiny(indice)+Stenseur(3,3)*Eloinz(indice))
     $                 *ctmp1
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
      enddo

      end
