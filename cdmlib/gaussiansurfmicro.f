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

      subroutine gaussiansurfmicro(eps,zcouche,neps,nepsmax,k0,waist,xx0
     $     ,yy0,zz0,thetat,phit,psit,Amp,aretecube,nx,ny,nz ,nfft2d
     $     ,nfft2dmax,Egausxref,Egausyref,Egauszref,Egausxtra,Egausytra
     $     ,Egausztra ,nstop ,infostr)
      implicit none
      integer neps,nepsmax,nfft2d,nfft2dmax,nstop,nx,ny,nz ,i,ntest
      double precision zcouche(0:nepsmax),waist,thetat,phit,psit,Amp,k0
     $     ,aretecube,kdirx,kdiry,kdirz
      double complex eps(0:nepsmax+1)
      integer nfft2d2,nkx,nky,k,ii,jj,kk,indice,kkk,nnn,test
      double precision kx,ky,kz,kxx,kyy,kzz,deltak,euler(3,3),kzt
     $     ,eulerinv(3,3),kmax,kp2,theta,phi,psi,pi,const,const1,var1
     $     ,var2,xx0,yy0,zz0,x0,y0,z0,fac
      double complex Egausxref(nfft2dmax *nfft2dmax),Egausyref(nfft2dmax
     $     *nfft2dmax) ,Egauszref(nfft2dmax*nfft2dmax)
     $     ,Egausxtra(nfft2dmax *nfft2dmax),Egausytra(nfft2dmax
     $     *nfft2dmax) ,Egausztra(nfft2dmax*nfft2dmax),E0x,E0y,E0z,Ex
     $     ,Ey ,Ez ,Axp,Ayp ,Azp,Axs,Ays,Azs,Ax,Ay,Az,const2,icomp,Arx
     $     ,Ary ,Arz ,Atx,Aty ,Atz
      character(64) infostr
      write(*,*) 'nstop',nstop
c     changement unite angle, psi =0 defini pol p
      pi=dacos(-1.d0)
      theta=thetat*pi/180.d0
      phi=phit*pi/180.d0
      psi=psit*pi/180.d0
      icomp=(0.d0,1.d0)
c     calcul du pas de discretisationen delta k
      nfft2d2=nfft2d/2
      kmax=k0*dsqrt(dreal(eps(0)))
      deltak=2.d0*pi/(dble(nfft2d)*aretecube)
c     var1=(xs(1)+dble(nfft2d2)*aretecube)*deltak
c     var2=(ys(1)+dble(nfft2d2)*aretecube)*deltak
      
c     translation du faisceau gaussien pour que l'objet soit bien placé
c     pour le premier indice en -N/2dx et -N/2dy.

      x0=xx0
      y0=yy0
      z0=zz0
c     x0=xx0
c     y0=yy0
      ntest=0
c     tet si milieu homogene
      do i=1,neps+1
         if (eps(i).ne.eps(0)) ntest=1
      enddo

      if (nx.ge.nfft2d) then
         nstop=1
         infostr='object larger than the FFT box along x'
      endif
      if (ny.eq.nfft2d) then
         nstop=1
         infostr='object larger than the FFT box along y'
      endif
      
      const=dexp(-deltak*deltak*waist*waist/2.d0)
      write(*,*) 'const',const,2.d0*pi/dsqrt(-dlog(0.5d0)/waist/waist
     $     *2.d0)/aretecube
      if (const.le.0.5d0) then
         infostr='size of FFT too small of the Gaussian beam'
         write(*,*) 'nfft nminimum',2.d0*pi/dsqrt(-dlog(0.5d0)/waist
     $        /waist*2.d0)/aretecube
         nstop=1
      endif
      write(*,*) 'coucou',nstop,infostr
      if (nstop.eq.1) return




c     calcul de la matrice de rotation et de son inverse definition de
c     la base du faisceai: ZZ=sin(theta) cos(phi) x+sin(theta) sin(phi)
c     y+ cos(theta) z: YY=ZZ ^ z/|ZZ ^ z|; XX= YY ^ ZZ
c     psi defini la polarisation dans la base XX, YY et ZZ.

c     psi =0 champ polarise suivant XX et psi =90 champ polarise selon
c     YY. YY reste toujours parallele a la surface.

      euler(1,1)=-dcos(phi)*dcos(theta)
      euler(1,2)=-dsin(phi)*dcos(theta)
      euler(1,3)=dsin(theta)
      euler(2,1)=dsin(phi)
      euler(2,2)=-dcos(phi)
      euler(2,3)=0.d0
      euler(3,1)=dcos(phi)*dsin(theta)
      euler(3,2)=dsin(phi)*dsin(theta)
      euler(3,3)=dcos(theta)


      eulerinv(1,1)=euler(1,1)
      eulerinv(1,2)=euler(2,1)
      eulerinv(1,3)=euler(3,1)
      eulerinv(2,1)=euler(1,2)
      eulerinv(2,2)=euler(2,2)
      eulerinv(2,3)=euler(3,2)
      eulerinv(3,1)=euler(1,3)
      eulerinv(3,2)=euler(2,3)
      eulerinv(3,3)=euler(3,3)
      
      fac=deltak*deltak

      kdirx=dcos(phi)*dsin(theta)*k0
      kdiry=dsin(phi)*dsin(theta)*k0
      kdirz=dcos(theta)*k0



c     initialise
      Egausxref=0.d0
      Egausyref=0.d0
      Egauszref=0.d0
      Egausxtra=0.d0
      Egausytra=0.d0
      Egausztra=0.d0

      do nkx=-nfft2d2,nfft2d2-1
c     range pour la FFT suivant x
c         if (nkx.ge.0) then
c            ii=nkx+1
c         else
c            ii=nfft2d+nkx+1
c         endif

         do nky=-nfft2d2,nfft2d2-1

            kx=deltak*dble(nkx)
            ky=deltak*dble(nky)
            kp2=kx*kx+ky*ky
            if (kp2.le.kmax*kmax*0.999d0) then
c     le faisceau se deplace vers les z positifs
               kz=dsqrt(kmax*kmax-kp2) 
c     Passe a la nouvelle base: zz est l'axe de propagtion du faisceau
               kxx=euler(1,1)*kx+euler(1,2)*ky+euler(1,3)*kz
               kyy=euler(2,1)*kx+euler(2,2)*ky+euler(2,3)*kz
               kzz=euler(3,1)*kx+euler(3,2)*ky+euler(3,3)*kz

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
               Ax=amp*dcos(psi)*Axp
               Ay=amp*dsin(psi)*Ays
               Az=amp*(dcos(psi)*Azp+dsin(psi)*Azs)
c     changement de base pour retourner dans la base de la surface
               E0x=eulerinv(1,1)*Ax+eulerinv(1,2)*Ay+eulerinv(1,3)*Az
               E0y=eulerinv(2,1)*Ax+eulerinv(2,2)*Ay+eulerinv(2,3)*Az
               E0z=eulerinv(3,1)*Ax+eulerinv(3,2)*Ay+eulerinv(3,3)*Az
c     appele du multicouche pour avoir le champ a l'origine (x=y=0)
               call champmultifft(eps,zcouche,neps,nepsmax,0.d0,k0 ,E0x
     $              ,E0y,E0z,kx,ky,infostr,nstop,Ex,Ey,Ez,Arx ,Ary ,Arz
     $              ,Atx,Aty,Atz)

c               indice=ii+nfft2d*(jj-1)
               indice=nkx+nfft2d2+1+nfft2d*(nky+nfft2d2)
               const2=cdexp(-icomp*(kx*x0+ky*y0+kz*z0))
c               write(*,*) 'const2',const2,xx0,yy0,zz0
               Egausxref(indice)=Arx*const2
               Egausyref(indice)=Ary*const2
               Egauszref(indice)=Arz*const2
               Egausxtra(indice)=Atx*const2
               Egausytra(indice)=Aty*const2
               Egausztra(indice)=Atz*const2
            endif
         enddo
      enddo
 
      if (ntest.eq.0) then
         Egausxref=0.d0
         Egausyref=0.d0
         Egauszref=0.d0
      endif
      
      end
