      subroutine computegcfft2dsurf(imax,deltakx,deltaky,k0,P0,theta,phi
     $     ,nfft2d,Ediffkzpos,Ediffkzneg,gasym,Cscai,w0,nepsmax ,neps
     $     ,dcouche,zcouche ,epscouche,ncote)
c     Calcule la section efficace de diffusion avec le champ lointain
c     Ediffkzpos(k||).  Attention Ediffkzpos represente le champ
c     lointain en r avec |r|=1 dans la direction k=(k||,kz).  Le flux du
c     champ diffracte à travers un plan infini s'ecrit
c     Int(Pointing.z dr||)=1/(2mu0c)Int[ |E(k||)|² / (kz k0) dk||]
      implicit none
      integer i,j,imax,ii,jj,nfft2d,ncote
      double precision Cscai,gasym,kx,ky,kz,deltakx,deltaky,k0,Emod,I0
     $     ,k03,theta,phi,pi
      double complex Ediffkzpos(nfft2d,nfft2d,3),Ediffkzneg(nfft2d
     $     ,nfft2d,3)
      integer nepsmax,neps
      double precision dcouche(nepsmax),zcouche(0:nepsmax),indice0
     $     ,indicen,w0,P0
      double complex epscouche(0:nepsmax+1)

      Cscai=0.d0
      gasym=0.d0
      k03=k0*k0*k0
      pi=dacos(-1.d0)
      indice0=dsqrt(dreal(epscouche(0)))
      indicen=dsqrt(dreal(epscouche(neps+1)))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kx,ky,kz,ii,jj,Emod)   
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2) REDUCTION(+:Cscai,gasym)
      do i=-imax,imax
         do j=-imax,imax

            kx=dble(i)*deltakx
            ky=dble(j)*deltaky
            if (indice0*indice0*k0*k0-kx*kx-ky*ky.gt.0.d0) then 
                           
               kz=dsqrt(indice0*indice0*k0*k0-kx*kx-ky*ky)
               ii=imax+i+1
               jj=imax+j+1

               
               Emod=cdabs(Ediffkzpos(ii,jj,1))**2+cdabs(Ediffkzpos(ii,jj
     $              ,2))**2+cdabs(Ediffkzpos(ii,jj,3))**2
               Cscai=Cscai+deltakx*deltaky*Emod/kz
               gasym=gasym+deltakx*deltaky*Emod/kz*(kx *dsin(theta
     $              *pi/180.d0)*dcos(phi*pi/180.d0) +ky
     $              *dsin(theta*pi/180.d0)* dsin(phi *pi /180.d0)
     $              +kz*dcos(theta*pi/180.d0))/k0
               
               Emod=cdabs(Ediffkzneg(ii,jj,1))**2+cdabs(Ediffkzneg(ii,jj
     $              ,2))**2+cdabs(Ediffkzneg(ii,jj,3))**2
               
               Cscai=Cscai+deltakx*deltaky*Emod/kz
               gasym=gasym+deltakx*deltaky*Emod/kz*(kx *dsin(theta
     $              *pi/180.d0)*dcos(phi*pi /180.d0) +ky
     $              *dsin(theta*pi/180.d0)*dsin(phi*pi /180.d0)-kz
     $              *dcos(theta*pi/180.d0))/k0

            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL           
      gasym=gasym/cscai
      Cscai=Cscai/P0/k0*indice0*pi*w0*w0/(8.d0*pi*1.d-7*299792458.d0)

      end
