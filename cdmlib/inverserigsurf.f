      subroutine inverserigsurf(xi,xr,nbsphere,ndipole,nx,ny,nz,nx2,ny2
     $     ,nxm,nym,nzm,nplanm,ntotalm,nmax,matindplan,Tabdip,b31 ,b32
     $     ,b33,FF,FF0,FFloc,b11,b12,b13,a11,a12,a13,a22,a23,a31,a32
     $     ,a33 ,WRK,nlar,ldabi,polarisa,methodeit,tol,tol1,nloop
     $     ,ncompte,planf,planb,nstop ,infostr)

      implicit none
      integer nbsphere,ndipole,nx,ny,nz,nx2,ny2 ,nxm,nym,nzm,nplanm
     $     ,ntotalm,nmax,nlar,ldabi,nloop,ncompte,nstop
      integer, dimension(nxm*nym*nzm) :: Tabdip
      integer matindplan(nzm,nzm)
      
      double precision tol,tol1

      double complex, dimension(3*nxm*nym*nzm) :: xr,xi
      double complex, dimension(3*nxm*nym*nzm,12) :: wrk
      double complex, dimension(3*nxm*nym*nzm) :: FF,FF0,FFloc
      double complex a11(2*nxm,2*nym,nplanm),a12(2*nxm,2*nym,nplanm),
     $     a13(2*nxm,2*nym,nplanm),a22(2*nxm,2*nym,nplanm),a23(2*nxm,2
     $     *nym,nplanm),a31(2*nxm,2*nym,nplanm),a32(2*nxm,2*nym,nplanm)
     $     , a33(2*nxm,2*nym,nplanm),b11(4*nxm*nym),b12(4*nxm*nym),b13(4
     $     *nxm*nym),b22(4*nxm*nym),b23(4*nxm*nym),b31(4 *nxm*nym),b32(4
     $     *nxm*nym),b33(4*nxm*nym)
      double complex, dimension(nxm*nym*nzm,3,3) :: polarisa

      character(64) infostr,message
      character(12) methodeit
c     ********************
      integer i,nlim,nbsphere3,ndim,nou,nstat,steperr,nt
      double precision NORM,t0,t1,t2,tole
      double complex ALPHA,BETA,GPETA,DZETA,R0RN,QMR1,QMR2,QMR3,QMR4
     $     ,QMR5,QMR6,QMR7,QMR8,QMR9
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer values(8),values2(8)
      double complex DOTS(4)
      integer *8 planf,planb
      
      ncompte=0
      nou=0
      nbsphere3=3*nbsphere
      ndim=nbsphere3
      nloop=0
      nlim=10000

      call cpu_time(t1)
      call date_and_time(date,time,zone,values)

      if (methodeit(1:7).eq.'GPBICG1') then
 2002    call GPBICG(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim,TOL
     $        ,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR)
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL           
         endif

         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32,a33,polarisa,planb,planf)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.ne.1) goto  2002
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
         write(*,*) 'tol1',tol1


      elseif (methodeit(1:7).eq.'GPBICG2') then
         write(*,*) 'methodeit',methodeit
 2009    call GPBICG2(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1           
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         endif
         
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         
         if (nstat.ne.1) goto  2009
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
      elseif (methodeit(1:10).eq.'GPBICGplus') then
         write(*,*) 'methodeit',methodeit
 2016    call GPBICGplus(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 

         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1           
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         endif
         
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         
         if (nstat.ne.1) goto  2016
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM   
         nloop=nloop+1
         
      elseif (methodeit(1:10).eq.'GPBICGsafe') then
         write(*,*) 'methodeit',methodeit
 2019    call GPBICGsafe(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 

         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1           
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         endif
         
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         
         if (nstat.ne.1) goto  2019
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM   
         nloop=nloop+1
         
      elseif (methodeit(1:8).eq.'GPBICGAR') then
         write(*,*) 'methodeit',methodeit
 2010    call GPBICGAR(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
         endif
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         
         
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         
         if (nstat.ne.1) goto  2010
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
      elseif (methodeit(1:9).eq.'GPBICGAR2') then
         write(*,*) 'methodeit',methodeit
 2011    call GPBICGAR2(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
         endif
         
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         
         if (nstat.ne.1) goto  2011
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
      elseif (methodeit(1:12).eq.'BICGSTARPLUS') then
         write(*,*) 'methodeit',methodeit
 2015    call GPBICGSTARPLUS(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP
     $        ,Nlim,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT
     $        ,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
         endif
         
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         
         if (nstat.ne.1) goto  2015
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         
      elseif (methodeit(1:7).eq.'GPBICOR') then
         write(*,*) 'methodeit',methodeit   
 2020    call GPBICOR(XI,XR,FF0,FFloc,ldabi,ndim,nlar,nou,WRK ,NLOOP
     $        ,Nlim,TOL,NORM,ALPHA,BETA,GPETA,DZETA,QMR1 ,QMR2,NSTAT
     $        ,STEPERR)

         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         ncompte=ncompte+1
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
         endif
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.ne.1) goto  2020
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
         
      elseif (methodeit(1:4).eq.'CORS') then
         write(*,*) 'methodeit',methodeit   
 2021    call CORS(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK ,NLOOP ,Nlim,TOL
     $        ,NORM,ALPHA,GPETA,DZETA,BETA,NSTAT ,STEPERR)

         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         ncompte=ncompte+1
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
         endif
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.ne.1) goto  2021
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM            
         
      elseif (methodeit(1:6).eq.'QMRCLA') then
         write(*,*) 'methodeit',methodeit            
 2003    call PIMZQMR(FFloc,XI,XR,FF0,WRK,NORM,LDABI,NDIM,NLAR,QMR1
     $        ,QMR2,QMR3,QMR4,QMR5,QMR6,QMR7,QMR8,QMR9,DOTS,NOU,NT
     $        ,nloop,NLIM,TOLE ,TOL ,NSTAT ,STEPERR)
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         
         ncompte=ncompte+1
         if (nstat.eq.1) then
            do i=1,nbsphere3
               xi(i)=FFloc(i)
            enddo
            nt=1
         endif
c     write(*,*) 'ncompte',ncompte,nt,tole
         if (nt.eq.1) then
            call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny
     $           ,nz,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $           ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11
     $           ,a12 ,a13,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
            
         elseif (nt.eq.2) then
c     calcul avec le transpose                   
            call produitfftmatvectsurmtrans(xi,xr,nbsphere,ndipole,nx
     $           ,ny,nz,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $           ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11
     $           ,a12 ,a13,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         endif
         
         if (nstat.ne.1) goto  2003
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM  

         
      elseif (methodeit(1:5).eq.'TFQMR') then
         write(*,*) 'methodeit',methodeit
 2004    call TFQMR(FFloc,Xi,XR,FF0,ldabi,ndim,nlar,nou,WRK,nloop
     $        ,NLIM,TOL,NORM,QMR1,QMR2,QMR3,QMR4,QMR5,QMR6,NSTAT
     $        ,STEPERR)
         
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               xi(i)=FFloc(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         endif
         
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)

         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2004

         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL          
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
         
      elseif (methodeit(1:5).eq.'CG') then
         write(*,*) 'methodeit',methodeit
 2005    call ZCG(XI,XR,FF0,NORM,WRK,QMR1,QMR2,QMR3,LDABI,NDIM,NLAR
     $        ,NOU,NSTAT,NLOOP,NLIM,TOLE,TOL)
         
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
         endif
         
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         
         if (nstat.ne.1) goto  2005
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM   
         
      elseif (methodeit(1:8).eq.'BICGSTAB') then
         write(*,*) 'methodeit',methodeit
 2006    call PIMZBICGSTAB(FFLOC,Xi,XR,FF0,ldabi,nlar,ndim,nou,WRK
     $        ,QMR1,QMR2,QMR3,NORM,TOL,nloop,nlim,NSTAT,STEPERR)
         
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               xi(i)=FFloc(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
         endif

         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)

         if (nstat.ne.1) goto  2006
         
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            FF(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(FF(i)*dconjg(FF(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM   
         
      elseif (methodeit(1:12).eq.'QMRBICGSTAB1') then
         write(*,*) 'methodeit',methodeit
         nt=1
 2007    call QMRBICGSTAB(FFloc,Xi,XR,FF0,ldabi,ndim,nlar,nou,WRK
     $        ,nloop,nloop,TOL,TOLE,NORM,QMR1,QMR2,QMR3,QMR4,QMR5
     $        ,QMR6,QMR7,QMR8,QMR9,NT,NSTAT ,STEPERR)
         
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
         endif
         
         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)
         
         if (nstat.ne.1) goto  2007
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM   
      elseif (methodeit(1:12).eq.'QMRBICGSTAB2') then
         write(*,*) 'methodeit',methodeit
         nt=2
 2008    call QMRBICGSTAB(FFloc,Xi,XR,FF0,ldabi,ndim,nlar,nou,WRK
     $        ,nloop,nloop,TOL,TOLE,NORM,QMR1,QMR2,QMR3,QMR4,QMR5
     $        ,QMR6,QMR7,QMR8,QMR9,NT,NSTAT ,STEPERR)
         
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO           
            do i=1,nbsphere3
               FFloc(i)=xi(i)       
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         endif

         call produitfftmatvectsurm(xi,xr,nbsphere,ndipole,nx,ny,nz
     $        ,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $        ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11,a12
     $        ,a13 ,a22,a23,a31,a32 ,a33,polarisa,planb,planf)

         if (nstat.ne.1) goto  2008
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
!$OMP DO  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)            
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM   
      else
         write(*,*) 'Iterative method not correct'
         nstop=1
         infostr='Iterative method not correct'
         return         
      endif

      call cpu_time(t2)
      call date_and_time(date,time,zone,values2)
      message='to solve Ax=b'
      call calculatedate(values2,values,t2,t1,message)

      write(*,*) 'Method iterative used             : ',methodeit
      write(*,*) 'Tolerance obtained                : ',tol1
      write(*,*) 'Tolerance asked                   : ',tol
      write(*,*) 'Number of product Ax needs        : ',ncompte
    
      end
