      integer function numerocouche(z,neps,nepsmax,zcouche)
      implicit none
      integer i,neps,nepsmax
      double precision z,zcouche(0:nepsmax)

      do i=0,neps
         if (z.le.zcouche(i)) then
            numerocouche=i
            return
         endif
      enddo
      numerocouche=neps+1
      return
      end
