      subroutine fonctiongreensurffft(nx,ny,ntp,nx2,ny2,nxm,nym,n1m
     $     ,nplan,nplanm,nmatim,nbs,ntotalm,aretecube,a,matind
     $     ,matindplan ,matindice ,matrange,b11,b12,b13,b22,b23,b31
     $     ,b32,b33,a11 ,a12,a13 ,a22,a23,a31,a32,a33,planb)
      implicit none
      integer i,j,nn,ll,jj,ii,kk,indice,nx,ny,ntp,nx2,ny2,nxm,nym,n1m
     $     ,nplan,nplanm,nmatim ,nbs,ntotalm,n1, FFTW_BACKWARD
      integer matindice(nplanm,nmatim) ,matind(0:2*n1m*n1m)
     $     ,matindplan(nplan,nplan)
      double complex matrange(nbs,5),Ixx,Ixy,Ixz,Izx,Izz,Sxx,Sxy,Sxz,Syy
     $     ,Szx,Syz,Szy
      double complex a11(2*nxm,2*nym ,nplanm),a12(2*nxm,2*nym,nplanm)
     $     ,a13(2*nxm,2*nym,nplanm) ,a22(2*nxm,2*nym,nplanm),a23(2*nxm,2
     $     *nym,nplanm),a31(2*nxm,2 *nym,nplanm),a32(2*nxm,2*nym,nplanm)
     $     ,a33(2*nxm,2*nym,nplanm),b11(ntotalm),b12(ntotalm)
     $     ,b13(ntotalm),b22(ntotalm),b23(ntotalm)
     $     ,b31(ntotalm) ,b32(ntotalm),b33(ntotalm) 
      double precision sphi,cphi,s2phi,c2phi,aretecube,a(0:2*n1m*n1m)
      integer*8 planb
      FFTW_BACKWARD=+1

      
c     initialise la matrice FFT
      do nn=1,ntp
         do ll=nn,ntp

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,indice,i,j,n1,kk)
!$OMP& PRIVATE(Ixx,Ixy,Ixz,Izz,Izx,sphi,cphi,s2phi,c2phi)
!$OMP& PRIVATE(Sxx,Sxy,Sxz,Syy,Syz,Szx,Szy)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)            
            do jj=1,ny2
               do ii=1,nx2
                  indice=ii+nx2*(jj-1)
                  if (ii.eq.nx+1.or.jj.eq.ny+1) then
                     b11(indice)=0.d0
                     b12(indice)=0.d0
                     b13(indice)=0.d0
                     b22(indice)=0.d0
                     b23(indice)=0.d0
                     b31(indice)=0.d0
                     b32(indice)=0.d0
                     b33(indice)=0.d0
                  else
                     if (ii.gt.nx) then
                        i=(ii-1)-nx2
                     else
                        i=ii-1
                     endif                    
                     if (jj.gt.ny) then
                        j=(jj-1)-ny2
                     else
                        j=jj-1
                     endif
                     n1=i*i+j*j
                     kk=matindice(matindplan(ll,nn),matind(n1))
                     Ixx=matrange(kk,1)
                     Ixy=matrange(kk,2)
                     Ixz=matrange(kk,3)
                     Izz=matrange(kk,4)
                     Izx=matrange(kk,5)
                     sphi=aretecube*dble(i)/a(n1)
                     cphi=aretecube*dble(j)/a(n1)
                     s2phi=2.d0*sphi*cphi
                     c2phi=cphi*cphi-sphi*sphi
                     Sxx=Ixx+c2phi*Ixy
                     Sxy=-s2phi*Ixy
                     Sxz=sphi*Ixz
                     Syy=Ixx-c2phi*Ixy
                     Syz=cphi*Ixz
                     Szx=sphi*Izx
                     Szy=cphi*Izx
c*******************************************************
                     b11(indice)=Sxx
                     b12(indice)=Sxy
                     b13(indice)=Sxz
                     b22(indice)=Syy
                     b23(indice)=Syz
                     b31(indice)=Szx
                     b32(indice)=Szy
                     b33(indice)=Izz
                  endif
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL              

            call dfftw_execute_dft(planb,b11,b11)
            call dfftw_execute_dft(planb,b12,b12)
            call dfftw_execute_dft(planb,b13,b13)
            call dfftw_execute_dft(planb,b22,b22)
            call dfftw_execute_dft(planb,b23,b23)
            call dfftw_execute_dft(planb,b31,b31)
            call dfftw_execute_dft(planb,b32,b32)
            call dfftw_execute_dft(planb,b33,b33)

            kk=matindplan(ll,nn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)             
            do jj=1,ny2
               do ii=1,nx2
                  indice=ii+nx2*(jj-1)
                  a11(ii,jj,kk)=b11(indice)
                  a12(ii,jj,kk)=b12(indice)
                  a13(ii,jj,kk)=b13(indice)
                  a22(ii,jj,kk)=b22(indice)
                  a23(ii,jj,kk)=b23(indice)
                  a31(ii,jj,kk)=b31(indice)
                  a32(ii,jj,kk)=b32(indice)
                  a33(ii,jj,kk)=b33(indice)
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL            
         enddo
      enddo


      end

