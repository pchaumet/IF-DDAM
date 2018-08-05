c$$$
c$$$!====================================================================
c$$$! Spline interpolation
c$$$! Comments: values of function f(x) are calculated in n base points
c$$$! then: spline coefficients are computed
c$$$!       spline interpolation is computed in 2n-1 points, 
c$$$!       a difference sum|f(u)-ispline(u)| 
c$$$!====================================================================
c$$$implicit none
c$$$integer, parameter :: n=11      ! base points for interpolation
c$$$integer, parameter :: nint=21   ! compute interpolation in nint points
c$$$double precision xmin, xmax     ! given interval of x()
c$$$double precision, dimension (n) :: xi(n), yi(n), b(n), c(n), d(n)
c$$$double precision x, y, step, ys, error, errav
c$$$integer i
c$$$double precision f, ispline
c$$$
c$$$! open files
c$$$!open (unit=1, file='tablex1.dat')
c$$$!open (unit=2, file='tablex2.dat')
c$$$
c$$$xmin = 0.0
c$$$xmax = 2.0
c$$$
c$$$! step 1: generate xi and yi from f(x), xmin, xmax, n
c$$$step = (xmax-xmin)/(n-1)
c$$$do i=1,n
c$$$  xi(i) = xmin + step*float(i-1) 
c$$$  yi(i) = f(xi(i)) 
c$$$!  write (*,200) xi(i), yi(i)
c$$$end do
c$$$
c$$$!  step 2: call spline to calculate spline coeficients
c$$$call spline (xi, yi, b, c, d,n) 
c$$$
c$$$!  step 3: interpolation at nint points
c$$$errav = 0.0
c$$$step = (xmax-xmin)/(nint-1)
c$$$write(*,201) 
c$$$do i=1, nint
c$$$  x = xmin + step*float(i-1) 
c$$$  y = f(x) 
c$$$  ys = ispline(x, xi, yi, b, c, d, n)
c$$$  error = ys-y
c$$$  write (*,200) x, ys, error
c$$$! step 4: calculate quality of interpolation               
c$$$  errav = errav + abs(y-ys)/nint 
c$$$end do
c$$$write (*,202) errav
c$$$200 format (3f12.5)
c$$$201 format ('        x        spline      error')    
c$$$202 format ('           Average error',f12.5)
c$$$
c$$$end program main
c$$$
c$$$
c$$$!
c$$$!  Function f(x)
c$$$!
c$$$  function f(x)
c$$$  double precision f, x
c$$$  f = sin(x) 
c$$$  end function f

      subroutine spline (x, y, b, c, d, n)
!======================================================================
!     Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!     for cubic spline interpolation
!     s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!     for  x(i) <= x <= x(i+1)
!     Alex G: January 2010
!----------------------------------------------------------------------
!     input..
!     x = the arrays of data abscissas (in strictly increasing order)
!     y = the arrays of data ordinates
!     n = size of the arrays xi() and yi() (n>=2)
!     output..
!     b, c, d  = arrays of spline coefficients
!     comments ...
!     spline.f90 program is based on fortran version of program spline.f
!     the accompanying function fspline can be used for interpolation
!======================================================================
      implicit none
      integer n
      double precision x(n)
      double complex y(n), b(n), c(n), d(n),h,uncomp,icomp
      integer i, j, gap
      
      uncomp=(1.d0,0.d0)
      icomp=(0.d0,1.d0)
      gap = n-1
!     check input
      if ( n.lt. 2 ) return
      if ( n.lt.3 ) then
         b(1) = (y(2)-y(1))/(x(2)-x(1)) ! linear interpolation
         c(1) = 0.d0
         d(1) = 0.d0
         b(2) = b(1)
         c(2) = 0.d0
         d(2) = 0.d0
         return
      end if
!     
!     step 1: preparation
!     
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, gap
         d(i) = x(i+1) - x(i)
         b(i) = 2.d0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
      end do
!     
!     step 2: end conditions 
!     
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.d0
      c(n) = 0.d0
      if (n.ne. 3) then
         c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
         c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
         c(1) = c(1)*d(1)**2/(x(4)-x(1))
         c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      end if
!     
!     step 3: forward elimination 
!     
      do i = 2, n
         h = d(i-1)/b(i-1)
         b(i) = b(i) - h*d(i-1)
         c(i) = c(i) - h*c(i-1)
      end do
!     
!     step 4: back substitution
!     
      c(n) = c(n)/b(n)
      do j = 1, gap
         i = n-j
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
!     
!     step 5: compute spline coefficients
!     
      b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.d0*c(n))
      do i = 1, gap
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.d0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.d0*c(i)
      end do
      c(n) = 3.d0*c(n)
      d(n) = d(n-1)
      end subroutine spline

      double complex function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
      implicit none
      integer n
      double precision  u, x(n)
      double complex y(n), b(n), c(n), d(n)
      integer i, j, k
      double precision dx

!     if u is ouside the x() interval take a boundary value (left or right)
      if(u.le.x(1)) then
         ispline = y(1)
         return
      end if
      if(u .ge. x(n)) then
         ispline = y(n)
         return
      end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
      i = 1
      j = n+1
      do while (j .gt. i+1)
         k = (i+j)/2
         if(u .le. x(k)) then
            j=k
         else
            i=k
         end if
      end do
!     *
!     evaluate spline interpolation
!     *
      dx = u - x(i)
      ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      end function ispline
