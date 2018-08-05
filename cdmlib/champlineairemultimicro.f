      subroutine champlineairemultimicro(eps,disz,neps,nepsmax,x,y,k0
     $     ,E0m,ssm,ppm,thetam,phim,nba,infostr,nstop,Arx,Ary,Arz,Atx
     $     ,Aty,Atz)

      implicit none
      integer nstop,neps,nepsmax,i,nba
      
      double precision disz(0:nepsmax),x,y,phit,thetat,phim(20)
     $     ,thetam(20),ssm(20),ppm(20),ss,pp,k0
      double complex eps(0:nepsmax+1),E0,E0m(20),Arx,Ary,Arz,Atx,Aty,Atz
     $     ,Arxm,Arym,Arzm,Atxm,Atym,Atzm
      

      character(64) infostr


      if (nba.gt.20) then
         nstop=1
         infostr='too many angle of incidence'
         return
      endif

      Arx=0.d0
      Ary=0.d0
      Arz=0.d0
      Atx=0.d0
      Aty=0.d0
      Atz=0.d0
   
      
      do i=1,nba
         pp=ppm(i)
         ss=ssm(i)
         thetat=thetam(i)
         phit=phim(i)
         E0=E0m(i)
         call champlineairemicro(eps,disz,neps,nepsmax,x,y,k0,E0,ss,pp
     $        ,thetat,phit,infostr,nstop,Arxm,Arym,Arzm,Atxm,Atym,Atzm)

         Arx=Arx+Arxm
         Ary=Ary+Arym
         Arz=Arz+Arzm
         Atx=Atx+Atxm
         Aty=Aty+Atym
         Atz=Atz+Atzm

      enddo

      end
