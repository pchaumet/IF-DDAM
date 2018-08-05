      SUBROUTINE ratint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      double precision x,xa(n),TINY
      double complex dy,y,ya(n)
      PARAMETER (NMAX=10,TINY=1.d-25)
      INTEGER i,m,ns
      double precision h,hh
      double complex c(NMAX),dd,d(NMAX),w,t,icomp
      ns=1
      icomp=(0.d0,1.d0)
      hh=dabs(x-xa(1))
      do 11 i=1,n
        h=dabs(x-xa(i))
        if (h.eq.0.d0)then
          y=ya(i)
          dy=0.d0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+(TINY,TINY)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.d0) stop
          dd=dreal(w)/dreal(dd)+icomp*dimag(w)/dimag(dd)
          d(i)=dreal(c(i+1))*dreal(dd)+icomp*dimag(c(i+1))*dimag(dd)
          c(i)=dreal(t)*dreal(dd)+icomp*dimag(t)*dimag(dd)
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
