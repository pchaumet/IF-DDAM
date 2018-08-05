      double complex function sommearg(a,b)
      implicit none
      double complex a,b
      double precision cr
      cr=dreal(b-a)
c      write(*,*) '****************************'
c      write(*,*) 'rout',a,b,ar,br,cr
      if (cr.le.-35.d0) then
         sommearg=a
      elseif (cr.ge.35.d0) then
         sommearg=b
      elseif (cr.ge.0.d0) then
         sommearg=b+cdlog(1.d0+cdexp(a-b))
      else
         sommearg=a+cdlog(1.d0+cdexp(b-a))
      endif
c      write(*,*) 'rout',sommearg
      end
