c     calcul pour une surface le produit: xr=xr+A xi avec FFT

c     Calcul le champ proche dans une boite (pleine) de dimension
c     nx,ny,ntp qui contient l'objet. E=E0+A*p fait avec des FFT.  En
c     entree xi est le dipole et en sortie xr le champ autour de l'objet
c     et le champ macro dans l'objet: champ macro dans la boite.
  
      subroutine produitfftmatvectsurplusboite(xi,xr,nbsphere,ndipole
     $     ,nx,ny,ntp ,nx2,ny2,nxm,nym,nzm,nplan,nplanm,ntotalm,nmax
     $     ,matindplan ,x1,x2,x3 ,b11,b21,b31,a11,a12,a13,a22,a23,a31
     $     ,a32 ,a33,planb,planf)
      implicit none
      integer i,j,m,ii,jj,nn,indice,ll,kk
      double complex x1(ntotalm),x2(ntotalm),x3(ntotalm),b11(ntotalm)
     $     ,b21(ntotalm),b31(ntotalm)
      integer nbsphere,ndipole,nx,ny,ntp,nx2,ny2,nxm,nym,nzm,nplanm
     $     ,nplan,ntotalm,nmax,matindplan(nplan,nplan)
      double complex xi(3*nmax),xr(3*nmax),a11(2*nxm,2*nym ,nplanm)
     $     ,a12(2*nxm,2*nym,nplanm),a13(2*nxm,2*nym,nplanm) ,a22(2*nxm,2
     $     *nym,nplanm),a23(2*nxm,2*nym,nplanm),a31(2*nxm,2 *nym,nplanm)
     $     ,a32(2*nxm,2*nym,nplanm),a33(2*nxm,2*nym,nplanm)
      integer*8 planf,planb
      double precision dntotal
      integer FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE,nxy
  
      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64
      nxy=nx2*ny2
      dntotal=dble(nxy)
  
c     calcul FFT du vecteur B
      do nn=1,ntp
 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)         
!$OMP DO SCHEDULE(STATIC)
         do i=1,ntotalm
            x1(i)=0.d0
            x2(i)=0.d0
            x3(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(i,j,jj,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)        
         do j=1,ny
            do i=1,nx
               indice=i+nx2*(j-1)                  
               jj=3*(i+nx*((j-1)+ny*(nn-1)))
               x1(indice)=xi(jj-2)
               x2(indice)=xi(jj-1)
               x3(indice)=xi(jj)
     
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
#ifdef USE_FFTW
         call dfftw_execute_dft(planb,x1,x1)
         call dfftw_execute_dft(planb,x2,x2)
         call dfftw_execute_dft(planb,x3,x3)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
         call fftsingletonz2d(x1,NX2,NY2,FFTW_BACKWARD)
!$OMP SECTION   
         call fftsingletonz2d(x2,NX2,NY2,FFTW_BACKWARD)
!$OMP SECTION  
         call fftsingletonz2d(x3,NX2,NY2,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif
         
         do ll=1,ntp
            kk=matindplan(ll,nn)
            if (ll.ge.nn) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2) 
               do j=1,ny2
                  do i=1,nx2
                     indice=i+nx2*(j-1)       
                   
                     b11(indice)=a11(i,j,kk)*x1(indice)+a12(i,j
     $                    ,kk)*x2(indice)+a13(i,j,kk)*x3(indice)
                     b21(indice)=a12(i,j,kk)*x1(indice)+a22(i,j
     $                    ,kk)*x2(indice)+a23(i,j,kk)*x3(indice)
                     b31(indice)=a31(i,j,kk)*x1(indice)+a32(i,j
     $                    ,kk)*x2(indice)+a33(i,j,kk)*x3(indice)
                  
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
            else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)          
               do j=1,ny2
                  do i=1,nx2
                     indice=i+nx2*(j-1)      
                     b11(indice)=a11(i,j,kk)*x1(indice)+a12(i,j
     $                    ,kk)*x2(indice)-a31(i,j,kk)*x3(indice)
                     b21(indice)=a12(i,j,kk)*x1(indice)+a22(i,j
     $                    ,kk)*x2(indice)-a32(i,j,kk)*x3(indice)
                     b31(indice)=-a13(i,j,kk)*x1(indice)-a23(i,j
     $                    ,kk)*x2(indice)+a33(i,j,kk)*x3(indice)
                   
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
            endif

#ifdef USE_FFTW
            call dfftw_execute_dft(planf,b11,b11)
            call dfftw_execute_dft(planf,b21,b21)
            call dfftw_execute_dft(planf,b31,b31)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
            call fftsingletonz2d(b11,NX2,NY2,FFTW_FORWARD)
!$OMP SECTION   
            call fftsingletonz2d(b21,NX2,NY2,FFTW_FORWARD)
!$OMP SECTION  
            call fftsingletonz2d(b31,NX2,NY2,FFTW_FORWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)     
            do j=1,ny
               do i=1,nx
                  indice=i+nx2*(j-1)           
                  ii=(i+nx*((j-1)+ny*(ll-1)))*3                
                  xr(ii-2)=xr(ii-2)+b11(indice)/dntotal
                  xr(ii-1)=xr(ii-1)+b21(indice)/dntotal
                  xr(ii)=xr(ii)+b31(indice)/dntotal
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
         enddo
      enddo

      
      return
      end
