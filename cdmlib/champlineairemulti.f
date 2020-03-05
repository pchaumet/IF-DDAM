      subroutine champlineairemulti(eps,disz,neps,nepsmax,x,y,z,k0,E0m
     $     ,ssm,ppm,thetam,phim,nba,infostr,nstop,Ex,Ey,Ez)

      implicit none
      integer nstop,neps,nepsmax,i,nba
      
      double precision disz(0:nepsmax),x,y,z,phit,thetat,phim(10)
     $     ,thetam(10),ssm(10),ppm(10),ss,pp,k0
      double complex eps(0:nepsmax+1),E0,E0m(10)
      
      double  complex Ex,Ey,Ez,Exs,Eys,Ezs
      character(64) infostr


      if (nba.gt.10) then
         nstop=1
         infostr='too many angle of incidence'
         return
      endif

      Ex=0.d0
      Ey=0.d0
      Ez=0.d0
      
      do i=1,nba
         pp=ppm(i)
c     ss=ssm(i)
         ss=1.d0
         thetat=thetam(i)
         phit=phim(i)
         E0=E0m(i)

         if (cdabs(E0).eq.0.d0) then
            nstop=1
            infostr='One of the plane wave has a null magnitude'
            return
         endif
         
         call champlineaire(eps,disz,neps,nepsmax,x,y,z,k0,E0,ss,pp
     $        ,thetat,phit,infostr,nstop,Exs,Eys,Ezs)

         Ex=Ex+Exs
         Ey=Ey+Eys
         Ez=Ez+Ezs
      enddo

      end
