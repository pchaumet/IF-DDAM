c     calcul pour une surface le produit: (I-alpha^t A^t)xi=xr avec FFT,
c     sachant que A^t=A
      subroutine produitfftmatvectsurmopttrans(xi,xr,nbsphere,ndipole,nx
     $     ,ny,ntp,nx2,ny2,nxm,nym,nzm,nplan,nplanm,ntotalm,nmax
     $     ,matindplan,x1,x2,x3,FF,b11,b21,b31,a11,a12,a13,a22
     $     ,a23,a31,a32,a33,polarisa,planb,planf)
      implicit none
      integer i,j,ii,jj,nn,indice,ll,kk
      double complex x1(ntotalm),x2(ntotalm),x3(ntotalm),b11(ntotalm)
     $     ,b21(ntotalm),b31(ntotalm)
      integer nbsphere,ndipole,nx,ny,ntp,nx2,ny2,nxm,nym,nzm,nplanm
     $     ,nplan,ntotalm,nmax,matindplan(nplan,nplan)
      double complex xi(3*nmax),xr(3*nmax),FF(3*nmax),a11(2*nxm,2*nym
     $     ,nplanm),a12(2*nxm,2*nym,nplanm),a13(2*nxm,2*nym,nplanm)
     $     ,a22(2*nxm,2*nym,nplanm),a23(2*nxm,2*nym,nplanm),a31(2*nxm,2
     $     *nym,nplanm),a32(2*nxm,2*nym,nplanm),a33(2*nxm,2*nym,nplanm)
     $     ,polarisa(nmax,3,3)
      integer*8 planf,planb
      double precision dntotal
      integer FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE,nxy
      
      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64
      nxy=nx2*ny2
      dntotal=dble(nxy)
      
c     calcul l'identite I xi
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
      do i=1,3*nbsphere
         FF(i)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
        
c     calcul FFT du vecteur B
      do nn=1,ntp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
      do i=1,nx2*ny2
         x1(i)=0.d0
         x2(i)=0.d0
         x3(i)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice,ii,jj)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)    
         do j=1,ny
            do i=1,nx
               indice=i+nx2*(j-1)                  
               ii=i+nx*((j-1)+ny*(nn-1))
               jj=3*ii
               x1(indice)=xi(jj-2)
               x2(indice)=xi(jj-1)
               x3(indice)=xi(jj)
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         call dfftw_execute_dft(planb,x1,x1)
         call dfftw_execute_dft(planb,x2,x2)
         call dfftw_execute_dft(planb,x3,x3)
         
         do ll=1,ntp
            kk=matindplan(ll,nn)
            if (ll.ge.nn) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
               do j=1,ny2
                  do i=1,nx2
                     indice=i+nx2*(j-1)       
c                     write(*,*) 'a',a11(i,j,kk),x1(indice),i,j ,kk   
                     b11(indice)=a11(i,j,kk)*x1(indice)+a12(i,j
     $                    ,kk)*x2(indice)+a13(i,j,kk)*x3(indice)
                     b21(indice)=a12(i,j,kk)*x1(indice)+a22(i,j
     $                    ,kk)*x2(indice)+a23(i,j,kk)*x3(indice)
                     b31(indice)=a31(i,j,kk)*x1(indice)+a32(i,j
     $                    ,kk)*x2(indice)+a33(i,j,kk)*x3(indice)
c                     write(*,*) 'bb1',b11(indice),i,j
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
            else
!$OMP PARALLEL  PRIVATE(i,j,indice)
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
c                     write(*,*) 'bbt1',b11(indice),i,j
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
            endif

            call dfftw_execute_dft(planf,b11,b11)
            call dfftw_execute_dft(planf,b21,b21)
            call dfftw_execute_dft(planf,b31,b31)
            

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice,ii)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            do j=1,ny
               do i=1,nx
                  indice=i+nx2*(j-1)           
                  ii=(i+nx*((j-1)+ny*(ll-1)))*3
                  FF(ii-2)=FF(ii-2)+b11(indice)/dntotal
                  FF(ii-1)=FF(ii-1)+b21(indice)/dntotal
                  FF(ii)=FF(ii)+b31(indice)/dntotal
c                  write(*,*) 'FF',FF(ii),ii,i,j,ll
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
         enddo
      enddo               

c     calcul de FF=A^t xi fait
c     calcul I-alpha^t FF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii)
!$OMP DO SCHEDULE(STATIC)
      do i=1,ndipole
         ii=3*i
         xr(ii-2)=xi(ii-2)-polarisa(i,1,1)*FF(ii-2)-polarisa(i,2,1)
     $        *FF(ii-1)-polarisa(i,3,1)*FF(ii)
         xr(ii-1)=xi(ii-1)-polarisa(i,1,2)*FF(ii-2)-polarisa(i,2,2)
     $        *FF(ii-1)-polarisa(i,3,2)*FF(ii)
         xr(ii)=xi(ii)-polarisa(i,1,3)*FF(ii-2)-polarisa(i,2,3)*FF(ii
     $        -1)-polarisa(i,3,3)*FF(ii)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      return
      end
