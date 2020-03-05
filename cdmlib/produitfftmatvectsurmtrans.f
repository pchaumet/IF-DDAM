c     calcul pour une surface le produit: (I-alpha^t A^t)xi=xr avec FFT,
c     sachant que A^t=A
      subroutine produitfftmatvectsurmtrans(xi,xr,nbsphere,ndipole,nx,ny
     $     ,ntp,nx2,ny2,nxm,nym,nzm,nplan,nplanm,ntotalm,nmax,matindplan
     $     ,Tabdip,x1,x2,x3,FF,b11,b21,b31,a11,a12,a13,a22,a23,a31,a32
     $     ,a33,polarisa,planb,planf)
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice,ii,jj,m)
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii,jj,m)
!$OMP DO SCHEDULE(STATIC)
      do i=1,ndipole
         ii=3*i
         if (Tabdip(i).ne.0) then
            jj=3*Tabdip(i)
            m=(jj+2)/3
            xr(jj-2)=xi(jj-2)-polarisa(m,1,1)*FF(jj-2)-polarisa(m,2,1)
     $           *FF(jj-1)-polarisa(m,3,1)*FF(jj)
            xr(jj-1)=xi(jj-1)-polarisa(m,1,2)*FF(jj-2)-polarisa(m,2,2)
     $           *FF(jj-1)-polarisa(m,3,2)*FF(jj)
            xr(jj)=xi(jj)-polarisa(m,1,3)*FF(jj-2)-polarisa(m,2,3)
     $           *FF(jj-1)-polarisa(m,3,3)*FF(jj)
         endif
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      return
      end
