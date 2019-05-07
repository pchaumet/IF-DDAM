      subroutine anglecalculmic(theta,phi,indice0,k0,deltakx,deltaky)
      implicit none
      integer imini,jmini
      double precision theta,phi,k0,indice0,deltakx,deltaky
      double precision kxinc,kyinc,pi,thetat,phit
      pi=dacos(-1.d0)
      
c     calcul le kx et ky correspondant au theta et phi
      kxinc=k0*dsin(theta*pi/180.d0)*dcos(phi*pi/180.d0)*indice0
      kyinc=k0*dsin(theta*pi/180.d0)*dsin(phi*pi/180.d0)*indice0
c     calcul du kx et ky le plus proche
      imini=nint(kxinc/deltakx)
      jmini=nint(kyinc/deltaky)

c     recalcule l'incident pour fitter avec le maillage
      if (kxinc.eq.0.d0.and.kyinc.eq.0.d0) then
         phit=phi
         thetat=theta
      else
         kxinc=dble(imini)*deltakx
         kyinc=dble(jmini)*deltaky
         phit=datan2(kyinc,kxinc)*180.d0/pi
         if (theta.lt.0.d0) phit=phit-180.d0
         if (phit.gt.180.d0) phit=phit-360.d0
         if (phit.lt.-180.d0) phit=phit+360.d0
         thetat=dsign(dasin(dsqrt(kxinc*kxinc+kyinc*kyinc)/k0/indice0)
     $        *180/pi,theta)
      endif



      theta=thetat
      phi=phit
      return
      end
