c     routine pour le multicouche pour acceder a l'espace libre
c     propagateur d espace libre en utilisant la symetrie
      subroutine propesplibdiagI(z,z0,k0,eps,Txx,Tzz)
      implicit none
      double precision z,z0,k0,rab,rab2
      double complex eps,const1,const2,const3,Txx,Txy,Tzz,Txz,kc

      kc=k0*cdsqrt(eps)
      Rab=dabs(z-z0)

      if (Rab.eq.0.d0) then
         Txx=0.d0
         Tzz=0.d0
         Txy=0.d0
         Txz=0.d0
         goto 10
      else
         rab2=1.d0/Rab/Rab
         const1=(Rab*(1.d0,0.d0))**(-3.d0)-(0.d0,1.d0)*kc*rab2
         const2=kc*kc/Rab*(1.d0,0.d0)
         const3=cdexp((0.d0,1.d0)*kc*Rab)
         Txx=(const2-const1)*const3
         Tzz=2.d0*const1*const3
      endif
 10   return
      end
c*********************************************************************
      subroutine propesplibI(a,z,z0,k0,eps,Txx,Txy,Txz,Tzz)
      implicit none
      double precision a,z,z0,k0,rab,a2,rab2
      double complex eps,const1,const2,const3,const4,Txx,Txy,Tzz,Txz,kc

      kc=k0*cdsqrt(eps)
      a2=a*a
      Rab=dsqrt(a2+(z-z0)*(z-z0))

      if (Rab.eq.0.d0) then
         Txx=0.d0
         Tzz=0.d0
         Txy=0.d0
         Txz=0.d0
         goto 10
      else
         rab2=1.d0/Rab/Rab
         const1=(Rab*(1.d0,0.d0))**(-3.d0)-(0.d0,1.d0)*kc*rab2
         const2=kc*kc/Rab*(1.d0,0.d0)
         const3=cdexp((0.d0,1.d0)*kc*Rab)
         const4=3.d0*const1-const2
         Txx=((const2-const1)+a2*const4*rab2/2.d0)*const3
         Tzz=((const2-const1)+(z-z0)*(z-z0)*const4*rab2)*const3
         Txy=-a2*const4*rab2/2.d0*const3
         Txz=a*(z-z0)*const4*rab2*const3
      endif
 10   return
      end
