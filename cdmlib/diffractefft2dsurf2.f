c     Cette routine calcul le champ diffracté à partir des dipoles.  On
c     rentre les dipôles en FF et les champs diffractés par FFT sont
c     sortis dans Ediffkzpos pour kz>0 et Ediffkzneg pour kz<0. Eloin
c     tableaux temporaires.

      subroutine diffractefft2dsurf2(nbsphere,nx,ny,nz,nxm,nym,nzm
     $     ,nfft2d,nfft2dmax,tabfft2,k0,xs,ys,zs,aretecube,Eloinx,Eloiny
     $     ,Eloinz ,FF,imax,deltakx,deltaky,Ediffkzpos,Ediffkzneg,r,zz
     $     ,numaperref,numapertra,nepsmax ,neps,dcouche,zcouche
     $     ,epscouche ,ncote,nstop ,infostr,planf)
c On calcule le champ diffracté E(r) avec |r|=1
      implicit none
      integer nx,ny,nz,nxm,nym,nzm,nfft2dmax,nfft2d,ncote,nstop
      double precision xs(nxm*nym*nzm),ys(nxm*nym*nzm),zs(nxm*nym*nzm)
     $     ,aretecube,k0,r,zz,za,numaperref,numapertra
      double complex FF(3*nxm*nym*nzm),Ediffkzpos(nfft2dmax,nfft2dmax,3)
     $     ,Ediffkzneg(nfft2dmax,nfft2dmax,3)

      integer nfft2d2,imax,i,j,k,tabfft2(nfft2d),indice,kk,ii,jj
      double precision deltakx,deltaky,var1,var2,kx,ky,kz,fac,pi
      double complex ctmp,ctmp1,icomp,Eloinx(nfft2dmax*nfft2dmax)
     $     ,Eloiny(nfft2dmax*nfft2dmax),Eloinz(nfft2dmax*nfft2dmax)
     $     ,Stenseur(3,3)

      integer nepsmax,neps
      double precision dcouche(nepsmax),zcouche(0:nepsmax),indice0
     $     ,indicen
      double complex epscouche(0:nepsmax+1)

      integer nbsphere,iii,iiii,jjjj
      double precision x,y,z
      double complex Em(3)

c     Info string
      character(64) infostr
      integer FFTW_FORWARD
      integer*8 planf

      FFTW_FORWARD=-1
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      nfft2d2=nfft2d/2  
      var1=(xs(1)+dble(nfft2d2)*aretecube)*deltakx
      var2=(ys(1)+dble(nfft2d2)*aretecube)*deltaky

      if (nfft2d.gt.16384) then
         nstop=1
         infostr='nfft2d too large'
         return
      endif
      
      
      if (deltakx.ge.k0) then
         nstop=1
         infostr='In FFT Poynting nfft2d too small'
         return
      endif

      if (ncote.eq.0.or.ncote.eq.1) then
         if (dreal(epscouche(neps+1)).le.0.d0) then
            infostr='epsilon <0 higher layer computation far field'
            nstop=-1
            return
         endif
          if (dimag(epscouche(neps+1)).gt.0.d0) then
            infostr='absorbing layer to computate far field'
            nstop=-1
            return
         endif
         indicen=dsqrt(dreal(epscouche(neps+1)))
      endif

      if (ncote.eq.0.or.ncote.eq.-1) then     
         if (dreal(epscouche(0)).le.0.d0) then
            infostr='epsilon <0 higher layer computation far field'
            nstop=-1
            return
         endif
          if (dimag(epscouche(0)).gt.0.d0) then
            infostr='absorbing layer to compute far field'
            nstop=-1
            return
         endif
         indice0=dsqrt(dreal(epscouche(0)))
      endif
      


      do i=1,nfft2d
         if (i-nfft2d2-1.ge.0) then
            tabfft2(i)=i-nfft2d2
         else
            tabfft2(i)=nfft2d2+i
         endif
      enddo

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)   
!$OMP DO SCHEDULE(STATIC)  COLLAPSE(2) 
         do i=1,nfft2d
            do j=1,nfft2d
               Ediffkzpos(i,j,1)=0.d0
               Ediffkzpos(i,j,2)=0.d0
               Ediffkzpos(i,j,3)=0.d0
               Ediffkzneg(i,j,1)=0.d0
               Ediffkzneg(i,j,2)=0.d0
               Ediffkzneg(i,j,3)=0.d0
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

         if (ncote.eq.0.or.ncote.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  
!$OMP& PRIVATE(i,j,ii,jj,kx,ky,kz,indice,ctmp1,Stenseur) 
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)
            do i=-imax,imax
               do j=-imax,imax
                  ii=imax+i+1
                  jj=imax+j+1
                  kx=dble(i)*deltakx
                  ky=dble(j)*deltaky
                  if (indicen*indicen*k0*k0*numapertra*numapertra
     $                 *0.9999d0-kx*kx-ky*ky.gt.0.d0) then
                     
                     kz=dsqrt(indicen*indicen*k0*k0-kx*kx-ky*ky) 
                     indice=tabfft2(i+nfft2d2+1)+nfft2d*(tabfft2(j
     $                    +nfft2d2+1)-1)
                     ctmp1=cdexp(-icomp*(var1*dble(i)+var2 *dble(j)))

                     call tenseurmulticoucheloinfft(kx,ky,kz,zz,zs(kk)
     $                    ,k0,nepsmax,neps,dcouche,zcouche,epscouche
     $                    ,Stenseur)

                     Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)+(Stenseur(1
     $                    ,1)*Eloinx(indice)+Stenseur(1,2)
     $                    *Eloiny(indice)+Stenseur(1,3)*Eloinz(indice))
     $                    *ctmp1
                     Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)+(Stenseur(2
     $                    ,1)*Eloinx(indice)+Stenseur(2,2)
     $                    *Eloiny(indice)+Stenseur(2,3)*Eloinz(indice))
     $                    *ctmp1                     
                     Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)+(Stenseur(3
     $                    ,1)*Eloinx(indice)+Stenseur(3,2)
     $                    *Eloiny(indice)+Stenseur(3,3)*Eloinz(indice))
     $                    *ctmp1
                  endif
               enddo               
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
         endif


         if (ncote.eq.0.or.ncote.eq.-1) then
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,ii,jj,kx,ky,kz,indice,ctmp1,Stenseur) 
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)
            do i=-imax,imax
               do j=-imax,imax
                  ii=imax+i+1
                  jj=imax+j+1
                  kx=dble(i)*deltakx
                  ky=dble(j)*deltaky
                  if (indice0*indice0*k0*k0*numaperref*numaperref
     $                 *0.9999d0-kx*kx-ky*ky.gt.0.d0) then
                     
                     kz=-dsqrt(indice0*indice0*k0*k0-kx*kx-ky*ky) 
                     indice=tabfft2(i+nfft2d2+1)+nfft2d*(tabfft2(j
     $                    +nfft2d2+1)-1)
                     ctmp1=cdexp(-icomp*(var1*dble(i)+var2*dble(j)))
c                     write(*,*) 'dessous2',ctmp1,indice,kx,ky,kz,-zz
c     $                    ,zs(kk),k0
                     call tenseurmulticoucheloinfft(kx,ky,kz,-zz,zs(kk)
     $                    ,k0,nepsmax,neps,dcouche,zcouche,epscouche
     $                    ,Stenseur)
c                     write(*,*) 'Stenseur',Stenseur

                     Ediffkzneg(ii,jj,1)=Ediffkzneg(ii,jj,1)+(Stenseur(1
     $                    ,1)*Eloinx(indice)+Stenseur(1,2)
     $                    *Eloiny(indice)+Stenseur(1,3)*Eloinz(indice))
     $                    *ctmp1

                     Ediffkzneg(ii,jj,2)=Ediffkzneg(ii,jj,2)+(Stenseur(2
     $                    ,1)*Eloinx(indice)+Stenseur(2,2)
     $                    *Eloiny(indice)+Stenseur(2,3)*Eloinz(indice))
     $                    *ctmp1
                     
                     Ediffkzneg(ii,jj,3)=Ediffkzneg(ii,jj,3)+(Stenseur(3
     $                    ,1)*Eloinx(indice)+Stenseur(3,2)
     $                    *Eloiny(indice)+Stenseur(3,3)*Eloinz(indice))
     $                    *ctmp1
                  endif
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
         endif
      enddo
c     supression de exp(i k0 r)/r car on travaille avec les amplitudes
c     des ondes planes


      end
