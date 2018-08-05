c     calcul pour une surface le produit: xr=xr+A xi avec FFT
      subroutine produitfftmatvectsurplus(xi,xr,nbsphere,ndipole,nx,ny
     $     ,ntp,nx2,ny2,nxm,nym,nzm,nplan,nplanm,ntotalm,nmax,matindplan
     $     ,Tabdip,x1,x2,x3,FF,b11,b21,b31,a11,a12,a13,a22,a23,a31,a32
     $     ,a33,planb,planf)
      implicit none
      integer i,j,m,ii,jj,nn,indice,ll,kk
      double complex x1(ntotalm),x2(ntotalm),x3(ntotalm),b11(ntotalm)
     $     ,b21(ntotalm),b31(ntotalm)
      integer nbsphere,ndipole,nx,ny,ntp,nx2,ny2,nxm,nym,nzm,nplanm
     $     ,nplan,ntotalm,nmax,matindplan(nplan,nplan),Tabdip(nxm*nym
     $     *nzm)
      double complex xi(3*nmax),xr(3*nmax),FF(3*nmax),a11(2*nxm,2*nym
     $     ,nplanm),a12(2*nxm,2*nym,nplanm),a13(2*nxm,2*nym,nplanm)
     $     ,a22(2*nxm,2*nym,nplanm),a23(2*nxm,2*nym,nplanm),a31(2*nxm,2
     $     *nym,nplanm),a32(2*nxm,2*nym,nplanm),a33(2*nxm,2*nym,nplanm)
      integer*8 planf,planb
      double precision dntotal
      integer FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE,nxy
      
      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64
      nxy=nx2*ny2
      dntotal=dble(nxy)

      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)         
!$OMP DO SCHEDULE(STATIC)     
      do i=1,3*nmax
         FF(i)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      
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

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(i,j,ii,jj,m,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)                
         do j=1,ny
            do i=1,nx
               indice=i+nx2*(j-1)                  
               ii=i+nx*((j-1)+ny*(nn-1))
               if (Tabdip(ii).ne.0) then
                  jj=3*Tabdip(ii)
                  m=(jj+2)/3
                  x1(indice)=xi(jj-2)
                  x2(indice)=xi(jj-1)
                  x3(indice)=xi(jj)
               endif                       
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
                     b31(indice)=-a13(i,j,kk)*x1(indice)-a23(i ,j
     $                    ,kk)*x2(indice)+a33(i,j,kk)*x3(indice)
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL              
            endif

            call dfftw_execute_dft(planf,b11,b11)
            call dfftw_execute_dft(planf,b21,b21)
            call dfftw_execute_dft(planf,b31,b31)
            
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            do j=1,ny
               do i=1,nx
                  indice=i+nx2*(j-1)           
                  ii=(i+nx*((j-1)+ny*(ll-1)))*3
                  FF(ii-2)=FF(ii-2)+b11(indice)
                  FF(ii-1)=FF(ii-1)+b21(indice)
                  FF(ii)=FF(ii)+b31(indice)
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL             
         enddo
      enddo

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii,jj,m)
!$OMP DO SCHEDULE(DYNAMIC)
      do i=1,ndipole
         ii=3*i
         if (Tabdip(i).ne.0) then
            jj=3*Tabdip(i)
            m=(jj+2)/3
            xr(jj-2)=xr(jj-2)+FF(ii-2)/dntotal
            xr(jj-1)=xr(jj-1)+FF(ii-1)/dntotal
            xr(jj)=xr(jj)+FF(ii)/dntotal
         endif
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL       
      end
