c     routine qui calcul un faisceau gaussien incident de waist waist
c     centre en x0,y0,z0 dans une boite parallélépipédique dans un
c     multicouche de normal z.  On se place dans la base XX,YY,ZZ du
c     faisceau. Le faisceau se propageant suivant ZZ et le faisceau
c     passe par un polariseur oriente suivant XX ou YY, puis le faisceau
c     peut etre devie par un miroir de telle sorte que sa direction de
c     propagation soit decrite par theta, phi defnit par : ZZ=sin(theta)
c     cos(phi) x+sin(theta) sin(phi) y+ cos(theta) z.YY=ZZ ^ z/|ZZ ^ z|;
c     XX= YY ^ ZZ

c     psi =0 champ polarise suivant XX et psi =90 champ polarise selon
c     YY. YY reste toujours parallele a la surface.

c     FF0 donne le champ total du faisceau gaussien dans le multicouche
c     en xs,ys,zs.

c     Egausxref ,Egausyref ,Egauszref: champ gaussian reflechi dans les
c     directions kx, ky.

c     Egausxtra, Egausytra, Egausztra : champ gaussian transmi dans les
c     directions kx, ky.

c     Dans ce code qui utilise des FFT, l'objet est toujours placé
c     indice=1 : X=-N/2dx et Y=-N/2dy. Or, dans le programme principal
c     l'objet est placé en xs(1) ys(1). Donc nous translatons le
c     faisceau gaussien de d=-N/2*dx -xs(1) et pareil en y.
c
c     Le faisceau gaussien s'ecrit Einc(r)=Int[ A(k||)exp(ik||r||+ikz z) dk||]
c     On calcule le flux incident à travers un plan infini avec :
c     fluxinc=1/(2mu0c) *4*pi²/k0 *Int[ |A(k||)|² kz dk|| ]

      subroutine gaussianiso(eps,zcouche,neps,nepsmax,k0,waist,xx0,yy0
     $     ,zz0,psit,E0,xs,ys,zs,FF0,aretecube,nx,ny,nz ,ndipole,nmax
     $     ,nfft2d,nfft2dmax,Egausxref,Egausyref,Egauszref ,Egausxtra
     $     ,Egausytra,Egausztra,fluxinc ,fluxref ,fluxtrans ,irra ,na
     $     ,masque, nmasque,nstop ,infostr,planb)
      implicit none
      integer neps,nepsmax,ndipole,nmax,nfft2d,nfft2dmax,nstop,nx,ny,nz
     $     ,i,nmasque
      double precision zcouche(0:nepsmax),waist,xs(nmax),ys(nmax)
     $     ,zs(nmax),psit,k0,aretecube,na
      double complex eps(0:nepsmax+1),FF0(3*nmax),E0
      integer nfft2d2,nkx,nky,k,ii,jj,kk,indice,kkk,nnn,test
      double precision kx,ky,kz,kxx,kyy,kzz,deltak,kzt,kmax,kp2,psi,pi
     $     ,const,const1,var1 ,var2,xx0 ,yy0,zz0,x0,y0,z0,fac,fluxinc
     $     ,fluxref,fluxtrans,irra,Axp ,Azp,Ays,Azs
      double complex Egausxref(nfft2dmax*nfft2dmax),Egausyref(nfft2dmax
     $     *nfft2dmax) ,Egauszref(nfft2dmax*nfft2dmax)
     $     ,Egausxtra(nfft2dmax *nfft2dmax),Egausytra(nfft2dmax
     $     *nfft2dmax) ,Egausztra(nfft2dmax*nfft2dmax),E0x,E0y,E0z,Ex,Ey
     $     ,Ez ,Ax,Ay,Az,const2,icomp,Arx,Ary ,Arz ,Atx ,Aty ,Atz
     $     ,masque(nfft2dmax,nfft2dmax)

      character(64) infostr
      integer FFTW_BACKWARD
      integer*8 planb
      
      write(*,*) 'info',infostr
      infostr='tt'
      write(*,*) 'info',infostr
      na=1.d0
      write(*,*) 'nstop',nstop
c     changement unite angle, psi =0 defini pol p (suivant x)
      pi=dacos(-1.d0)
    
      psi=psit*pi/180.d0
      icomp=(0.d0,1.d0)
c     calcul du pas de discretisationen delta k
      nfft2d2=nfft2d/2
      kmax=k0*dsqrt(dreal(eps(0)))
      deltak=2.d0*pi/(dble(nfft2d)*aretecube)
      fluxinc=0.d0
      fluxref=0.d0
      fluxtrans=0.d0
      write(*,*) 'ddd',deltak,kmax,k0,waist,xx0,yy0 ,zz0,psit,E0
c      var1=(xs(1)+dble(nfft2d2)*aretecube)*deltak
c      var2=(ys(1)+dble(nfft2d2)*aretecube)*deltak
      
c     translation du faisceau gaussien pour que l'objet soit bien placé
c     pour le premier indice en -N/2dx et -N/2dy.
      x0=xx0-(xs(1)+dble(nfft2d2)*aretecube)
      y0=yy0-(ys(1)+dble(nfft2d2)*aretecube)
      z0=zz0
      write(*,*) 'xxxx',xx0,yy0,zz0,waist
      if (nx.ge.nfft2d) then
         nstop=1
         infostr='object larger than the FFT box along x'
      endif
      if (ny.eq.nfft2d) then
         nstop=1
         infostr='object larger than the FFT box along y'
      endif
      write(*,*) 'infi',infostr,deltak,waist,deltak*deltak*waist*waist
     $     /2.d0
      const=dexp(-deltak*deltak*waist*waist/2.d0)
      write(*,*) 'const',const,2.d0*pi/dsqrt(-dlog(0.5d0)/waist/waist
     $     *2.d0)/aretecube
      write(*,*) infostr
      if (const.le.0.5d0) then
         write(*,*) 'coucou',infostr
         infostr='Size of FFT too small for the Iso beam'
         write(*,*) infostr
         write(*,*) 'nfft nminimum',2.d0*pi/dsqrt(-dlog(0.5d0)/waist
     $        /waist*2.d0)/aretecube
         nstop=1
         return
      endif
      write(*,*) 'nmasque',nmasque

c     creation du masque
      if (nmasque.eq.1) then
         do nkx=-nfft2d2,nfft2d2-1
            do nky=-nfft2d2,nfft2d2-1
               kx=deltak*dble(nkx)
               ky=deltak*dble(nky)
               kp2=kx*kx+ky*ky
               if (kp2.le.na*na*kmax*kmax*0.999d0) then
                  kz=dsqrt(kmax*kmax-kp2)
                  if (sin(kz*zz0).le.0.d0) then
                     masque(nfft2d2+nkx+1,nfft2d2+nky+1)=-(1.d0,0.d0)
                  else
                     masque(nfft2d2+nkx+1,nfft2d2+nky+1)=(1.d0,0.d0)
                  endif               
               else
                  masque(nfft2d2+nkx+1,nfft2d2+nky+1)=(0.d0,0.d0)
               endif
            enddo
         enddo
      else
         nstop=1
         infostr='bad mask' 
      endif

      write(*,*) 'coucou',nstop,infostr
      
      if (nstop.eq.1) return


c     calcul de la matrice de rotation et de son inverse definition de
c     la base du faisceai: ZZ=sin(theta) cos(phi) x+sin(theta) sin(phi)
c     y+ cos(theta) z: YY=ZZ ^ z/|ZZ ^ z|; XX= YY ^ ZZ
c     psi defini la polarisation dans la base XX, YY et ZZ.

c     psi =0 champ polarise suivant XX et psi =90 champ polarise selon
c     YY. YY reste toujours parallele a la surface.

      
      fac=deltak*deltak



c     Boucle sur en z 
      do k=1,nz
c     kk indice de la cote suivant z
c         write(*,*) 'kk',kk
         kk=1+nx*ny*(k-1)

c     initialise
         Egausxref=0.d0
         Egausyref=0.d0
         Egauszref=0.d0
       
c     Commence la boucle sur les delta k dans le repere x,y,z de la
c     surface
         do nkx=-nfft2d2,nfft2d2-1

c     range pour la FFT suivant x
            if (nkx.ge.0) then
               ii=nkx+1
            else
               ii=nfft2d+nkx+1
            endif

            do nky=-nfft2d2,nfft2d2-1
c     range pour la FFT suivant y
               if (nky.ge.0) then
                  jj=nky+1
               else
                  jj=nfft2d+nky+1
               endif

               kx=deltak*dble(nkx)
               ky=deltak*dble(nky)
               kp2=kx*kx+ky*ky
               if (kp2.le.na*na*kmax*kmax*0.999d0) then
c     le faisceau se deplace vers les z positifs
                  kz=dsqrt(kmax*kmax-kp2)

c     Passe a la nouvelle base: zz est l'axe de propagtion du faisceau
                  kxx=kx
                  kyy=ky
                  kzz=kz
c     calcul amplitude de l'onde dans la nouvelle base: waist place en
c     zz0 c'est a dire dans la nouvelle base (celle du faisceau)

c                  const=dexp(-(kxx*kxx+kyy*kyy)*waist*waist/4.d0)*waist
c     $                 /2.d0/pi
                  const=dexp(-(kxx*kxx+kyy*kyy)*waist*waist/2.d0)*waist

c     const1 permet de normaliser le vecteur d'onde base XX, YY et ZZ.

c     polarisation p
                  const1=1.d0/dsqrt(kzz*kzz+kxx*kxx)
                  Axp=kzz*const*const1              
                  Azp=-kxx*const*const1
c     polarisation s
                  const1=1.d0/dsqrt(kzz*kzz+kyy*kyy)                 
                  Ays=kzz*const*const1
                  Azs=-kyy*const*const1 
c     composition des deux polarisations
                  Ax=E0*dcos(psi)*Axp
                  Ay=E0*dsin(psi)*Ays
                  Az=E0*(dcos(psi)*Azp+dsin(psi)*Azs)
c     changement de base pour retourner dans la base de la surface
                  E0x=Ax
                  E0y=Ay
                  E0z=Az
                  write(*,*) 'E000',cdabs(E0x)**2+cdabs(E0y)**2
     $                 +cdabs(E0z)**2,nkx,nky
c     appele du multicouche pour avoir le champ a l'origine (x=y=0)
                  call champmultifft(eps,zcouche,neps,nepsmax,zs(kk) ,k0
     $                 ,E0x,E0y,E0z,kx,ky,infostr,nstop,Ex,Ey,Ez,Arx
     $                 ,Ary,Arz,Atx,Aty,Atz)
               
                  indice=ii+nfft2d*(jj-1)
c     const2=cdexp(icomp*(var1*dble(nkx)+var2*dble(nky)))
                  const2=cdexp(-icomp*(kx*x0+ky*y0))*masque(nfft2d2+nkx
     $                 +1,nfft2d2+nky+1)
                  Egausxref(indice)=Ex*const2
                  Egausyref(indice)=Ey*const2
                  Egauszref(indice)=Ez*const2

               endif

            enddo
         enddo

c     fin boucle en delta k
c     calcul de la FFT
#ifdef USE_FFTW
         call dfftw_execute_dft(planb,Egausxref,Egausxref)
         call dfftw_execute_dft(planb,Egausyref,Egausyref)
         call dfftw_execute_dft(planb,Egauszref,Egauszref)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      call fftsingletonz2d(Egausxref,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Egausyref,0,0,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Egauszref,0,0,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif


c     shift+remet dans FF0, le champ incident
         do nkx=-nfft2d2,-nfft2d2+nx-1
c     range pour la FFT suivant x
            if (nkx.ge.0) then
               ii=nkx+1
            else
               ii=nfft2d+nkx+1
            endif

            do nky=-nfft2d2,-nfft2d2+ny-1
c     range pour la FFT suivant y

               if (nky.ge.0) then
                  jj=nky+1
               else
                  jj=nfft2d+nky+1
               endif  

               indice=ii+nfft2d*(jj-1)
c     nnn c'est le numero du dipole dans la boite object qui correspond
c     a une maille de la FFT nkx,nky=indice
               nnn=(nkx+nfft2d2+1)+nx*(nky+nfft2d2)+nx*ny*(k-1)
c     write(*,*) 'nnn',nnn,nx,ny,ii,jj,k
c     write(*,*) 'nnn2',nnn,nx,ny,ii,jj,k
               kkk=3*(nnn-1)
               FF0(kkk+1)=Egausxref(indice)*fac
               FF0(kkk+2)=Egausyref(indice)*fac
               FF0(kkk+3)=Egauszref(indice)*fac

               write(*,*) 'FF0',(cdabs(FF0(kkk+1))**2+cdabs(FF0(kkk+2))
     $              **2+cdabs(FF0(kkk+3))**2)/deltak/deltak/deltak
     $              /deltak ,zs(kk)
            enddo
         enddo

      enddo
c*************************************
c     calcul des flux et faisceau reflechi et transmis
c*************************************
      Egausxref=0.d0
      Egausyref=0.d0
      Egauszref=0.d0
      Egausxtra=0.d0
      Egausytra=0.d0
      Egausztra=0.d0
      i=0
      do nkx=-nfft2d2,nfft2d2-1
c     range pour la FFT suivant x
c         if (nkx.ge.0) then
c            ii=nkx+1
c         else
c            ii=nfft2d+nkx+1
c         endif

         do nky=-nfft2d2,nfft2d2-1
c     write(*,*) 'nkxy',nkx,nky
c     range pour la FFT suivant y
c            if (nky.ge.0) then
c               jj=nky+1
c            else
c               jj=nfft2d+nky+1
c            endif
            
            kx=deltak*dble(nkx)
            ky=deltak*dble(nky)
            kp2=kx*kx+ky*ky
            if (kp2.le.kmax*kmax*0.999d0) then
c     le faisceau se deplace vers les z positifs
               kz=dsqrt(kmax*kmax-kp2) 
c     Passe a la nouvelle base: zz est l'axe de propagtion du faisceau
               kxx=kx
               kyy=ky
               kzz=kz

c     calcul amplitude de l'onde dans la nouvelle base: waist place en
c     zz0 c'est a dire dans la nouvelle base (celle du faisceau)
c               const=dexp(-(kxx*kxx+kyy*kyy)*waist*waist/4.d0)*waist
c     $              /2.d0/pi

                  const=dexp(-(kxx*kxx+kyy*kyy)*waist*waist/2.d0)*waist

c     polarisation p
               const1=1.d0/dsqrt(kzz*kzz+kxx*kxx)
               Axp=kzz*const*const1              
               Azp=-kxx*const*const1
c     polarisation s
               const1=1.d0/dsqrt(kzz*kzz+kyy*kyy)             
               Ays=kzz*const*const1
               Azs=-kyy*const*const1 
c     composition des deux polarisations
               Ax=E0*dcos(psi)*Axp
               Ay=E0*dsin(psi)*Ays
               Az=E0*(dcos(psi)*Azp+dsin(psi)*Azs)
c     changement de base pour retourner dans la base de la surface
               E0x=Ax
               E0y=Ay
               E0z=Az
c     appele du multicouche pour avoir le champ a l'origine (x=y=0)
               call champmultifft(eps,zcouche,neps,nepsmax,zs(1),k0
     $              ,E0x,E0y,E0z,kx,ky,infostr,nstop,Ex,Ey,Ez,Arx ,Ary
     $              ,Arz,Atx,Aty,Atz)

c               indice=ii+nfft2d*(jj-1)
               indice=nkx+nfft2d2+1+nfft2d*(nky+nfft2d2)
               const2=cdexp(-icomp*(kx*xx0+ky*yy0))*masque(nfft2d2+nkx
     $              +1,nfft2d2+nky+1)
c               write(*,*) 'const2',const2,xx0,yy0,zz0
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
               fluxinc=fluxinc+(cdabs(E0x)**2+cdabs(E0y)**2 +cdabs(E0z)
     $              **2)*kz 
               fluxref=fluxref+(cdabs(Arx)**2+cdabs(Ary)**2 +cdabs(Arz)
     $              **2)*kz
               kzt=dreal(cdsqrt(eps(neps+1)*k0*k0-kp2))
               fluxtrans=fluxtrans+(cdabs(Atx)**2+cdabs(Aty)**2
     $              +cdabs(Atz)**2)*kzt
               i=i+1
            endif
         enddo
      enddo
      write(*,*) 'fluxinc',fluxinc,i,fac
      irra=fac/(8.d0*pi*1.d-7*299792458.d0)/k0*4.d0*pi*pi
c      irra=fac/(4.d0*pi*1.d-7*299792458.d0)/k0*4.d0*pi*pi
      fluxinc=fluxinc*irra
      fluxref=fluxref*irra
      fluxtrans=fluxtrans*irra

      write(*,*) 'fluxxx',fluxinc,fluxref,fluxtrans,(fluxref+fluxtrans)
     $     /fluxinc

      irra=fluxinc/pi/waist/waist
      

      end
