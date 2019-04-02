c     routine qui calcul un speckle polarise centre en x0,y0,z0 dans une
c     boite parallélépipédique dans un multicouche de normal z.  Le
c     faisceau se propageant suivant z et le faisceau passe par un
c     polariseur oriente suivant x ou y.  psi =0 champ polarise suivant
c     x et psi =90 champ polarise selon y.

c     FF0 donne le champ total du faisceau speckle dans le multicouche
c     en xs,ys,zs.

c     Egausxref ,Egausyref ,Egauszref: champ speckle reflechi dans les
c     directions kx, ky.

c     Egausxtra, Egausytra, Egausztra : champ speckle transmi dans les
c     directions kx, ky.
      
c     Le faisceau speckle s'ecrit Einc(r)=Int[ A(k||)exp(ik||r||+ikz z)
c     dk||]. On calcule le flux incident à travers un plan infini avec :
c     fluxinc=1/(2mu0c) *4*pi²/k0 *Int[ |A(k||)|² kz dk|| ]

      subroutine specklesurf(eps,zcouche,neps,nepsmax,k0,Amp,numaper,IR
     $     ,xs,ys,zs,x0,y0,z0,psi,FF0,aretecube,nx,ny,nz,nxm,nym,nzm
     $     ,ndipole,nmax,nfft2d ,nfft2dmax ,Egausxref,Egausyref
     $     ,Egauszref,Egausxtra,Egausytra ,Egausztra,fluxinc ,fluxref
     $     ,fluxtrans,irra,nstop ,infostr,planb)
      implicit none
      integer neps,nepsmax,ndipole,nmax,nfft2d,nfft2dmax,nstop,nx,ny,nz
     $     ,IR,iref,nmx,nmy,nxm,nym,nzm,i
      double precision zcouche(0:nepsmax),waist,xs(nmax),ys(nmax)
     $     ,zs(nmax),Amp,k0,aretecube,numaper
      double complex eps(0:nepsmax+1),FF0(3*nmax)
      integer nfft2d2,nkx,nky,k,ii,jj,kk,indice,kkk,nnn,test
      double precision kx,ky,kz,kxx,kyy,kzz,deltak,euler(3,3),kzt
     $     ,eulerinv(3,3),kmax,kp2,pi,const,const1,var1,var2,fac
     $     ,fluxinc,fluxref,fluxtrans,irra,ran,sunif,x0,y0,z0,psi
      double complex Egausxref(nfft2dmax *nfft2dmax),Egausyref(nfft2dmax
     $     *nfft2dmax),Egauszref(nfft2dmax*nfft2dmax)
     $     ,Egausxtra(nfft2dmax *nfft2dmax),Egausytra(nfft2dmax
     $     *nfft2dmax),Egausztra(nfft2dmax*nfft2dmax),Ex,Ey,Ez,Axp,Ayp
     $     ,Azp,Axs,Ays,Azs,Ax,Ay,Az,const2,icomp,Arx,Ary,Arz,Atx,Aty
     $     ,Atz
      character(64) infostr
      integer FFTW_BACKWARD
      integer*8 planb
      FFTW_BACKWARD=+1
      
      write(*,*) 'Speckle illumination',x0,y0,z0,IR

      pi=dacos(-1.d0)    
      icomp=(0.d0,1.d0)
c     calcul du pas de discretisationen delta k
      nfft2d2=nfft2d/2
      kmax=k0*dsqrt(dreal(eps(0)))
      deltak=2.d0*pi/(dble(nfft2d)*aretecube)
      fluxinc=0.d0
      fluxref=0.d0
      fluxtrans=0.d0
      nmx=(nxm-nx)/2
      nmy=(nym-ny)/2

      iref=4*IR+1
      ran=SUNIF(iref)

      if (nx.ge.nfft2d) then
         nstop=1
         infostr='object larger than the FFT box along x'
      endif
      if (ny.eq.nfft2d) then
         nstop=1
         infostr='object larger than the FFT box along y'
      endif
      if (nstop.eq.1) return
     
      fac=deltak*deltak

c     initialisation de l'amplitude des ondes planes.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)   
         do i=1,nfft2dmax*nfft2dmax
            Egausxtra(i)=0.d0
            Egausytra(i)=0.d0
            Egausztra(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL   

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nkx,nky,kx,ky,kp2,kz)   
!$OMP& PRIVATE(ran,const,const1,Axp,Azp,Ays,Azs,indice)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)         
      do nkx=-nfft2d2,nfft2d2-1
         do nky=-nfft2d2,nfft2d2-1           
            kx=deltak*dble(nkx)
            ky=deltak*dble(nky)
            kp2=kx*kx+ky*ky
            if (kp2.le.kmax*kmax*numaper*0.999d0) then
c     le faisceau se deplace vers les z positifs
               kz=dsqrt(kmax*kmax-kp2)               
               ran=SUNIF(-1)
               const=cdexp(icomp*2.d0*pi*ran)
               
c     polarisation p
               const1=1.d0/dsqrt(kz*kz+kx*kx)
               Axp=kz*const1*const
               Azp=-kx*const1*const
c     polarisation s
               const1=1.d0/dsqrt(kz*kz+ky*ky)       
               Ays=kz*const1*const
               Azs=-ky*const1*const
c     composition des deux polarisations
               indice=nkx+nfft2d2+1+nfft2d*(nky+nfft2d2)
               Egausxtra(indice)=amp*dcos(psi)*Axp
               Egausytra(indice)=amp*dsin(psi)*Ays
               Egausztra(indice)=amp*(dcos(psi)*Azp+dsin(psi)*Azs)
            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

c     Commence boucle sur en z 
      do k=1,nz
c     kk indice de la cote suivant z
c         write(*,*) 'kk',kk
         kk=1+nx*ny*(k-1)

c     initialise
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)     
         do i=1,nfft2dmax*nfft2dmax
            Egausxref(i)=0.d0
            Egausyref(i)=0.d0
            Egauszref(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
c     Commence la boucle sur les delta k dans le repere x,y,z de la
c     surface

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nkx,nky,ii,jj,kx,ky,kp2,kz)   
!$OMP& PRIVATE(indice,Ax,Ay,Az,Ex,Ey,Ez,Arx,Ary,Arz,Atx,Aty,Atz,const2)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)         
         
         do nkx=-nfft2d2,nfft2d2-1
            do nky=-nfft2d2,nfft2d2-1

c     range pour la FFT suivant x
               if (nkx.ge.0) then
                  ii=nkx+1
               else
                  ii=nfft2d+nkx+1
               endif               
c     range pour la FFT suivant y
               if (nky.ge.0) then
                  jj=nky+1
               else
                  jj=nfft2d+nky+1
               endif

               kx=deltak*dble(nkx)
               ky=deltak*dble(nky)
               kp2=kx*kx+ky*ky
               if (kp2.le.kmax*kmax*numaper*0.999d0) then
                  indice=nkx+nfft2d2+1+nfft2d*(nky+nfft2d2)
                  Ax=Egausxtra(indice)
                  Ay=Egausytra(indice)
                  Az=Egausztra(indice)
                  kz=dsqrt(kmax*kmax-kp2)
c     appele du multicouche pour avoir le champ a l'origine (x=y=0)
                  call champmultifft(eps,zcouche,neps,nepsmax,zs(kk) ,k0
     $                 ,Ax,Ay,Az,kx,ky,infostr,nstop,Ex,Ey,Ez,Arx ,Ary
     $                 ,Arz,Atx,Aty,Atz)                               
c     const2=cdexp(icomp*(var1*dble(nkx)+var2*dble(nky)))
                  const2=cdexp(-icomp*(kx*x0+ky*y0+kz*z0))
                  indice=ii+nfft2d*(jj-1)
                  Egausxref(indice)=Ex*const2
                  Egausyref(indice)=Ey*const2
                  Egauszref(indice)=Ez*const2

               endif

            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

c     fin boucle en delta k
c     calcul de la FFT

         call dfftw_execute_dft(planb,Egausxref,Egausxref)
         call dfftw_execute_dft(planb,Egausyref,Egausyref)
         call dfftw_execute_dft(planb,Egauszref,Egauszref)

c     shift+remet dans FF0, le champ incident
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nkx,nky,ii,jj,indice,nnn,kkk)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)            
         do nkx=-nfft2d2+nmx,-nfft2d2+nmx+nx-1
            do nky=-nfft2d2+nmy,-nfft2d2+nmy+ny-1
c     range pour la FFT suivant x
            if (nkx.ge.0) then
               ii=nkx+1
            else
               ii=nfft2d+nkx+1
            endif
c     range pour la FFT suivant y

               if (nky.ge.0) then
                  jj=nky+1
               else
                  jj=nfft2d+nky+1
               endif  

               indice=ii+nfft2d*(jj-1)
c     nnn c'est le numero du dipole dans la boite object qui correspond
c     a une maille de la FFT nkx,nky=indice
               nnn=(nkx+nfft2d2+1-nmx)+nx*(nky+nfft2d2-nmy)+nx*ny*(k-1)
               kkk=3*(nnn-1)
               FF0(kkk+1)=Egausxref(indice)*fac
               FF0(kkk+2)=Egausyref(indice)*fac
               FF0(kkk+3)=Egauszref(indice)*fac
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      enddo
c*************************************
c     calcul des flux et faisceau reflechi et transmis
c*************************************
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
      do i=1,nfft2dmax*nfft2dmax
         Egausxref(i)=0.d0
         Egausyref(i)=0.d0
         Egauszref(i)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nkx,nky,kx,ky,kp2,kz,indice)   
!$OMP& PRIVATE(Ax,Ay,Az,Ex,Ey,Ez,Arx,Ary,Arz,Atx,Aty,Atz,const2,kzt)
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:fluxinc,fluxref,fluxtrans)
!$OMP&  COLLAPSE(2)           
      do nkx=-nfft2d2,nfft2d2-1
         do nky=-nfft2d2,nfft2d2-1

            kx=deltak*dble(nkx)
            ky=deltak*dble(nky)
            kp2=kx*kx+ky*ky
            if (kp2.le.kmax*kmax*numaper*0.999d0) then

               indice=nkx+nfft2d2+1+nfft2d*(nky+nfft2d2)
               Ax=Egausxtra(indice)
               Ay=Egausytra(indice)
               Az=Egausztra(indice)
c     appele du multicouche pour avoir le champ a l'origine (x=y=0)
               call champmultifft(eps,zcouche,neps,nepsmax,zs(1),k0
     $              ,Ax,Ay,Az,kx,ky,infostr,nstop,Ex,Ey,Ez,Arx ,Ary
     $              ,Arz,Atx,Aty,Atz)
               kz=dsqrt(kmax*kmax-kp2)
               const2=cdexp(-icomp*(kx*x0+ky*y0+kz*z0))
               Egausxref(indice)=Arx*const2
               Egausyref(indice)=Ary*const2
               Egauszref(indice)=Arz*const2
               Egausxtra(indice)=Atx*const2
               Egausytra(indice)=Aty*const2
               Egausztra(indice)=Atz*const2
c               write(*,*) 'champ',Egausxref(indice),Egausyref(indice)
c     $              ,Egauszref(indice),Egausxtra(indice)
c     $              ,Egausytra(indice),Egausztra(indice),nkx,nky,indice
c     calcul du fluxinc: somme des fluxinc de chacunes des ondes planes
             
               fluxinc=fluxinc+(cdabs(Ax)**2+cdabs(Ay)**2 +cdabs(Az)
     $              **2)*kz 
               fluxref=fluxref+(cdabs(Arx)**2+cdabs(Ary)**2 +cdabs(Arz)
     $              **2)*kz
               kzt=dreal(cdsqrt(eps(neps+1)*k0*k0-kp2))
               fluxtrans=fluxtrans+(cdabs(Atx)**2+cdabs(Aty)**2
     $              +cdabs(Atz)**2)*kzt
            
            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      irra=fac/(8.d0*pi*1.d-7*k0*299792458.d0)/k0*4.d0*pi*pi
      fluxinc=fluxinc*irra
      fluxref=fluxref*irra
      fluxtrans=fluxtrans*irra

      write(*,*) 'fluxxx',fluxinc,fluxref,fluxtrans,(fluxref+fluxtrans)
     $     /fluxinc

      irra=0.d0
      

      end
