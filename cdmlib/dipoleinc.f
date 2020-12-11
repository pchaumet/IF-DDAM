      subroutine dipoleinc(xdip,ydip,zdip,thetat,phit,x,y,z,aretecube,k0
     $     ,E0,Ex,Ey,Ez,neps,nepsmax,dcouche,zcouche ,epscouche,tolinit
     $     ,nstop ,infostr)
      implicit none
      integer nstop,neps,nepsmax
      double precision x,y,z,k0,kb,pi,xdip,ydip,zdip,theta,phi,thetat
     $     ,phit,aretecube,dist,dcouche(nepsmax),zcouche(0:nepsmax),vol
     $     ,tolinit
      double complex E0,Ex,Ey,Ez,icomp,uncomp,propaesplibre(3,3),p(3)
     $     ,const,epscouche(0:nepsmax+1)

      integer l,nt,numerocouche
      double complex Sxx,Sxy,Sxz,Syy,Syz ,Szx,Szy,Ixx,Ixy,Ixz ,Izx ,Izz
      double precision eps0,aa,sphi,cphi,s2phi,c2phi,hc,epsabs

      
      character(64) infostr
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      
      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0

c     calcul orientation et amplitude du dipole
      
      p(1)=E0*dsin(theta)*dcos(phi)
      p(2)=E0*dsin(theta)*dsin(phi)
      p(3)=E0*dcos(theta)

      hc=0.3d0
      epsabs=0.d0
      nt=1
      
      dist=dsqrt((xdip-x)*(xdip-x)+(ydip-y)*(ydip-y)+(zdip-z)*(zdip-z))

c     calcul ou est le dipôle dans le couche et le epsilon
      l=numerocouche(zdip,neps,nepsmax ,zcouche)
      eps0=dreal(epscouche(l))
      vol=1.d0/(4.d0*pi)
      
      if (dist.le.aretecube/100.d0) then
         kb=k0*dsqrt(eps0)
         const=-1.d0/3.d0/(aretecube*aretecube*aretecube)/epscouche(l)
     $        +icomp*2.d0/3.d0*kb*kb*kb/vol

         Ex=const*p(1)
         Ey=const*p(2)
         Ez=const*p(3)
      else

         call tenseurmulticouchecomp(hc,tolinit,epsabs,x,y,z,xdip,ydip
     $        ,zdip,k0,neps,dcouche,zcouche,epscouche,Ixx,Ixy,Ixz ,Izx
     $        ,Izz)
         aa=dsqrt((xdip-x)*(xdip-x)+(ydip-y)*(ydip-y))
         if (aa.eq.0) aa=1.d300
         sphi=aretecube*(x-xdip)/aa
         cphi=aretecube*(y-ydip)/aa
         s2phi=2.d0*sphi*cphi
         c2phi=cphi*cphi-sphi*sphi
         Sxx=Ixx+c2phi*Ixy
         Sxy=-s2phi*Ixy
         Sxz=sphi*Ixz
         Syy=Ixx-c2phi*Ixy
         Syz=cphi*Ixz
         Szx=sphi*Izx
         Szy=cphi*Izx
         
         Ex=(Sxx*p(1)+Sxy*p(2)+Sxz*p(3))*vol
         Ey=(Sxy*p(1)+Syy*p(2)+Syz*p(3))*vol
         Ez=(Szx*p(1)+Szy*p(2)+Izz*p(3))*vol
c         write(*,*) '***',Sxz,Syz,Izz,aa,x,y,z,xdip,ydip,zdip
      endif
      
      end
c****************************************************
c****************************************************
c****************************************************
      subroutine dipoleinside(xdip,ydip,zdip,xs,ys,zs,aretecube,nmax
     $     ,nbsphere,zcouche,neps,nepsmax,nstop ,infostr)
      implicit none
      integer i,nmax,nbsphere,ic,jc,kc,neps,nepsmax,nstop

      double precision xdip,ydip,zdip,aretecube,dist
      double precision xs(nmax),ys(nmax),zs(nmax),zcouche(0:nepsmax)
      character(64) infostr

c     test si le dipole eclairant est au moins a une demi maille de tous
c     les dipoles (dans ce cas exterieur a l'objet), si proche alors mis
c     a la position du dipole le plus proche.
c     +test si dipole eclairant proche d'une couche
      

      do i=1,nbsphere
         ic=nint((xdip-xs(i))/aretecube)
         jc=nint((ydip-ys(i))/aretecube)
         kc=nint((zdip-zs(i))/aretecube)
        
         
         if (ic.eq.0.and.jc.eq.0.and.kc.eq.0) then
            xdip=xs(i)
            ydip=ys(i)
            zdip=zs(i)
            write(*,*) 'Antenna inside the near field domain'
            write(*,*) 'location fits a subunit of the object'
            write(*,*) 'x=',xdip,'m'
            write(*,*) 'y=',ydip,'m'
            write(*,*) 'z=',zdip,'m'
            return
         endif
      enddo

      write(*,*) 'antenna outside the near field domain'

c     test si dipole proche d'une couche
      do i=0,neps
         dist=dabs(zdip-zcouche(i))
         if (dist.le.aretecube/4.d0) then
            infostr='Antenna too close to an interface' 
            nstop=1
         endif
      enddo
      
      
      end
