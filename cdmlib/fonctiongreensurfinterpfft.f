      subroutine fonctiongreensurfinterpfft(nx,ny,ntp,nx2,ny2,nxm,nym
     $     ,n1m,nplan,nplanm,nmatim,nbs,ntotalm,ninter,ninterp,aretecube
     $     ,a ,matind ,matindplan ,matindice ,matrange,b11,b12,b13 ,b22
     $     ,b23,b31 ,b32,b33,a11 ,a12,a13 ,a22,a23,a31,a32,a33,planb)
      implicit none
      integer i,j,nn,ll,jj,ii,kk,indice,nx,ny,ntp,nx2,ny2,nxm,nym,n1m
     $     ,nplan,nplanm,nmatim ,nbs,ntotalm,n1,FFTW_BACKWARD
      integer matindice(nplanm,nmatim) ,matind(0:2*n1m*n1m)
     $     ,matindplan(nplan,nplan),ninterp
      double complex matrange(nbs,5),Ixx,Ixy,Ixz,Izx,Izz,Sxx,Sxy,Sxz,Syy
     $     ,Szx,Syz,Szy
      double complex a11(2*nxm,2*nym ,nplanm),a12(2*nxm,2*nym,nplanm)
     $     ,a13(2*nxm,2*nym,nplanm) ,a22(2*nxm,2*nym,nplanm),a23(2*nxm,2
     $     *nym,nplanm),a31(2*nxm,2 *nym,nplanm),a32(2*nxm,2*nym,nplanm)
     $     ,a33(2*nxm,2*nym,nplanm),b11(ntotalm),b12(ntotalm)
     $     ,b13(ntotalm),b22(ntotalm),b23(ntotalm)
     $     ,b31(ntotalm) ,b32(ntotalm),b33(ntotalm) 
      double precision sphi,cphi,s2phi,c2phi,aretecube,a(0:2*n1m*n1m)

      integer nav1,nav2,nap1,nap2,n,nlong,kkav1,kkav2,kkap1,kkap2
      double precision long,If1,If2,If3,If4
      
      integer ninter
   
      double precision xa(10),dinterp
      double complex ya(10),y,dy
      integer*8 planb
      FFTW_BACKWARD=+1
      n=4
      ya=0.d0
      xa=0.d0
      dy=0.d0
      dinterp=dble(ninterp)
c     initialise la matrice FFT
      do nn=1,ntp
         do ll=nn,ntp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,indice,i,j,n1,kk)
!$OMP& PRIVATE(Ixx,Ixy,Ixz,Izz,Izx,sphi,cphi,s2phi,c2phi)
!$OMP& PRIVATE(Sxx,Sxy,Sxz,Syy,Syz,Szx,Szy)
!$OMP& FIRSTPRIVATE(xa,ya,dy)
!$OMP& PRIVATE(long,nlong,nav1,nav2,nap1,nap2,kkav1,kkav2,kkap1,kkap2)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)  
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
c     reconstruction des fonctions de Green par interpolation
                     n1=i*i+j*j
                     if (n1.eq.0) then
                        kk=matindice(matindplan(ll,nn),1)
                        Ixx=matrange(kk,1)
                        Ixy=matrange(kk,2)
                        Ixz=matrange(kk,3)
                        Izz=matrange(kk,4)
                        Izx=matrange(kk,5)
                        long=1.d300
                     else
                        long=dsqrt(dble(n1*ninterp*ninterp))
                        nlong=nint(long)
                        if (dabs(dble(nlong)-long).le.1.d-10) then
c     test pour savoir si long est un entier
                           kk=matindice(matindplan(ll,nn),nlong+1)
                           Ixx=matrange(kk,1)
                           Ixy=matrange(kk,2)
                           Ixz=matrange(kk,3)
                           Izz=matrange(kk,4)
                           Izx=matrange(kk,5)
                        else
                           nav1=floor(long)
                           nav2=nav1-1
                           nap1=ceiling(long)
                           nap2=nap1+1
                           
                           kkav1=matindice(matindplan(ll,nn),nav1+1)
                           kkav2=matindice(matindplan(ll,nn),nav2+1)
                           kkap1=matindice(matindplan(ll,nn),nap1+1)
                           kkap2=matindice(matindplan(ll,nn),nap2+1)

c     interpolation par une fonction rationelle a quatre abscisses

                           xa(1)=dble(nav2)
                           xa(2)=dble(nav1)
                           xa(3)=dble(nap1)
                           xa(4)=dble(nap2)


                           ya(1)=matrange(kkav2,1)
                           ya(2)=matrange(kkav1,1)
                           ya(3)=matrange(kkap1,1)
                           ya(4)=matrange(kkap2,1)
                           call ratint(xa,ya,n,long,Ixx,dy)

                           ya(1)=matrange(kkav2,2)
                           ya(2)=matrange(kkav1,2)
                           ya(3)=matrange(kkap1,2)
                           ya(4)=matrange(kkap2,2)
                           call ratint(xa,ya,n,long,Ixy,dy)

                           ya(1)=matrange(kkav2,3)
                           ya(2)=matrange(kkav1,3)
                           ya(3)=matrange(kkap1,3)
                           ya(4)=matrange(kkap2,3)
                           call ratint(xa,ya,n,long,Ixz,dy)

                           ya(1)=matrange(kkav2,4)
                           ya(2)=matrange(kkav1,4)
                           ya(3)=matrange(kkap1,4)
                           ya(4)=matrange(kkap2,4)
                           call ratint(xa,ya,n,long,Izz,dy)

                           ya(1)=matrange(kkav2,5)
                           ya(2)=matrange(kkav1,5)
                           ya(3)=matrange(kkap1,5)
                           ya(4)=matrange(kkap2,5)
                           call ratint(xa,ya,n,long,Izx,dy)


                        endif 
                     endif

                     sphi=dble(i)/long*dinterp
                     cphi=dble(j)/long*dinterp
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

#ifdef USE_FFTW
            call dfftw_execute_dft(planb,b11,b11)
            call dfftw_execute_dft(planb,b12,b12)
            call dfftw_execute_dft(planb,b13,b13)
            call dfftw_execute_dft(planb,b22,b22)
            call dfftw_execute_dft(planb,b23,b23)
            call dfftw_execute_dft(planb,b31,b31)
            call dfftw_execute_dft(planb,b32,b32)
            call dfftw_execute_dft(planb,b33,b33)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      call fftsingletonz3d(b11,0,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(b12,0,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(b13,0,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(b22,0,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(b23,0,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(b31,0,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(b32,0,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(b33,0,0,0,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL  
#endif

             
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
