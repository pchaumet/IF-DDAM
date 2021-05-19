      subroutine shiftzcouche(neps,nepsmax,dcouche,zcouche,epscouche
     $     ,zshift)
      implicit none
      integer neps,nepsmax

      double precision dcouche(nepsmax),zcouche(0:nepsmax),zshift
      double complex epscouche(0:nepsmax+1)

      integer i,im
      double precision epsmax

      if (neps.eq.0) then
         zshift=zcouche(0)
         return
      endif

c     repere le plus fort epsilon et le numero de la couche
      epsmax=0.d0
      do i=0,neps+1
         epsmax=max(epsmax,cdabs(epscouche(i)))
      enddo
      do i=0,neps+1
         if (cdabs(epscouche(i)).ge.epsmax*0.999999999d0) im=i
      enddo
     
      zshift=zcouche(im)
      return
      end
