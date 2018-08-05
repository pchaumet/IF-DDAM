      subroutine poladiffcomp(rayon1,eps,eps0,k0,dddis, methode
     $     ,polarisabilite)
      implicit none
      integer dddis
      double precision rayon,volume,pi,k0,SLDR,rayon1
      double complex eps,rapport,pola,icomp,pola1
      double complex polarisabilite,k03,eps0,indice
      character(2) methode

c     dis: 0 si c est une sphere et 1 c est un element de discretisation
c     rayon=rayonon si dis=0 sinon facteur pour le cube -->sphere
c     methode: LA:lakthakia; CM: claussius-Mossotti;
c     CR: CM avec reaction de rayononnement;DB: dungey borhen;
c     DR: draine;DD: Draine 2

c      write(*,*) rayon,eps,eps0,k0,methode,dis

      pi=dacos(-1.d0)
      rapport=(eps-eps0)/(eps+2.d0*eps0)
      icomp=(0.d0,1.d0)
      indice=cdsqrt(eps0)
      k03=k0*k0*k0*indice*indice*indice
      rayon=rayon1
      
      
      if (dddis.eq.0) then
         volume=rayon*rayon*rayon
      elseif (dddis.eq.1) then
         rayon=rayon*((0.75d0/pi)**(1.d0/3.d0))
         volume=rayon*rayon*rayon
      else
         write(*,*) 'mauvaise valeur de dddis',dddis
         stop
      endif

c     polarisabilite de CM
      pola=rapport*volume

      if (methode.eq.'CM') then
         polarisabilite=eps0*pola
      elseif (methode.eq.'RR') then
         polarisabilite=eps0*pola/(1.d0-2.d0/3.d0*icomp*k03*pola)
c         polarisabilite=pola
      elseif (methode.eq.'GB') then
         polarisabilite=eps0*pola/(1.d0-2.d0/3.d0*icomp*k03*pola-
     *        k0*k0*indice*indice*pola/rayon)
      elseif (methode.eq.'LA') then         
         polarisabilite=eps0*pola/(1.d0-2.d0*rapport*((1.d0-icomp*k0
     $        *indice*rayon)*cdexp(icomp*k0*rayon*indice)-1.d0))        
      elseif (methode.eq.'LR') then
         SLDR=1.d0/5.d0
         pola1=pola/(1.d0+pola*(-1.8915316d0+eps/eps0* 0.1648469d0-eps
     $        /eps0*1.7700004d0*SLDR)*k0*indice*indice*k0/rayon)
         polarisabilite=eps0*pola1/(1.d0-2.d0/3.d0*icomp*k03*pola1)     
      else
         write(*,*) 'methode incorrecte'
         write(*,*) 'methode',methode
         stop
      endif

      end
