      subroutine bessel(z,norder,nseq,ierr,nseqlda,Jbessel)
      implicit none

      integer norder,nseq,ierr,nseqlda
      double complex z,Jbessel(nseqlda)

      integer nref,i,nmax,j
      double precision zm
      double complex Jn1,Jn,jbes,Ji
      
      ierr=0
      if (nseq.gt.nseqlda) ierr=1
      if (norder.lt.0) ierr=2
      if (nseq.le.0) ierr=3
      if (ierr.ne.0) return

      nmax=norder+nseq-1
      zm=cdabs(z)
      if (zm.eq.0) then
         do i=1,nseq
            Jbessel(i)=0.d0
         enddo
         if (norder.eq.0) Jbessel(1)=1.d0
         return
      endif

      if (zm.ge.20.d0) then
c     Calcul Bessel pour de grands arguments
         if (nseq.eq.1) then
            call besselpl(nmax,z,jbes,ierr)
            Jbessel(1)=jbes
         elseif (nseq.eq.2) then
            call besselpl(nmax,z,jbes,ierr)
            Jbessel(2)=jbes
            call besselpl(nmax-1,z,jbes,ierr)
            Jbessel(1)=jbes
         else
            call besselpl(nmax,z,jbes,ierr)
            Jbessel(nseq)=jbes
            call besselpl(nmax-1,z,jbes,ierr)
            Jbessel(nseq-1)=jbes
            do i=nmax-2,norder,-1
               j=i-nmax+nseq
               Jbessel(j)=dble(2*(i+1))/z*Jbessel(j+1)-Jbessel(j+2)
            enddo
         endif
c     Fin calcul Bessel pour de grands arguments
      else
c     Calcul Bessel pour  arguments  petits
         nref=idnint(zm**2/8.d0)

c     calcul z petit
         if (nref.eq.0) then
            if (nseq.eq.1) then
               call besselp(nmax,z,jbes,ierr)
               Jbessel(1)=jbes
            elseif (nseq.eq.2) then
               call besselp(nmax,z,jbes,ierr)
               Jbessel(2)=jbes
               call besselp(nmax-1,z,jbes,ierr)
               Jbessel(1)=jbes
            else
               call besselp(nmax,z,jbes,ierr)
               Jbessel(nseq)=jbes
               call besselp(nmax-1,z,jbes,ierr)
               Jbessel(nseq-1)=jbes
               do i=nmax-2,norder,-1
                  j=i-nmax+nseq
                  Jbessel(j)=dble(2*(i+1))/z*Jbessel(j+1)-Jbessel(j+2)
               enddo
            endif
            return
c fin calcul z petits
         else
            nref=max(nmax,nref)
            if (nref.eq.nmax) then
               call besselp(nmax,z,jbes,ierr)
               Jbessel(nseq)=jbes
               if (nseq.eq.1) return
               call besselp(nmax-1,z,jbes,ierr)
               Jbessel(nseq-1)=jbes
               if (nseq.eq.2) return
               do i=nmax-2,norder,-1
                  j=i-nmax+nseq
                  Jbessel(j)=dble(2*(i+1))/z*Jbessel(j+1)-Jbessel(j+2)
               enddo
            elseif (nref.eq.nmax+1) then
               call besselp(nref,z,Jn1,ierr)         
               call besselp(nref-1,z,Jn,ierr)
               Jbessel(nseq)=Jn
               do i=nref-2,norder,-1
                  Ji=dble(2*(i+1))/z*Jn-Jn1
                  Jn1=Jn
                  Jn=Ji
                  if (i.le.nmax) Jbessel(i-norder+1)=Ji
               enddo
            else
               call besselp(nref,z,Jn1,ierr)         
               call besselp(nref-1,z,Jn,ierr)            
               do i=nref-2,norder,-1
                  Ji=dble(2*(i+1))/z*Jn-Jn1
                  Jn1=Jn
                  Jn=Ji
                  if (i.le.nmax) Jbessel(i-norder+1)=Ji
               enddo
            endif
         endif
      endif

      end
c     ***********************************************************
      subroutine besselp(n,z,jbes,ierr)
      implicit none
      integer i,n,p,ierr
      double complex z,jbes,arg,arga,facz
      double precision facnp,expx,facp

      facnp=1.d0
      facz=-z*z/4.d0
      do i=2,n
         facnp=facnp*dble(i)
      enddo
      arg=(z/2.d0)**n/facnp
      jbes=arg
      do p=1,1000         
         arg=arg/dble(p*(n+p))*facz
         jbes=jbes+arg
         if (cdabs(arg/jbes).le.1.d-15) goto 10
      enddo
      ierr=5
 10   continue
      return
      
      end
c     ***********************************************
      subroutine besselpl(n,z,jbes,ierr)
      implicit none
      integer i,n,p,nr,nnr,ierr,nn4
      double complex z,jbes,arg,arga,facz,cosz,sinz,argexp
      double complex pn,qn,puis,icomp,facpuis,facpuis2,argp,argq,prodp
      double precision fac,pi
      jbes=0.d0
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      facz=cdsqrt(2.d0/pi/z)
      argexp=z-dble(n)*pi/2.d0-pi/4.d0
      cosz=(cdexp(icomp*argexp)+cdexp(-icomp*argexp))/2.d0
      sinz=(cdexp(icomp*argexp)-cdexp(-icomp*argexp))/2.d0/icomp
      facpuis=8.d0*z
      facpuis2=facpuis*facpuis
      if (facpuis.eq.0) then
         jbes=1.d0 
         return
      endif
      puis=(1.d0,0.d0)
      pn=1.d0
      qn=dble(4*n*n-1)/8.d0/z
      prodp=(1.d0,0.d0)
      nr=-1
      nn4=4*n*n
      do p=1,100
         
         nr=nr+2
         prodp=prodp*dble((nn4-(2*nr-1)*(2*nr-1))*(4*n*n-(2*nr+1)*(2 *nr
     $        +1)))/dble(nr*(nr+1))/facpuis2
         argp=(-1.d0)**dble(p)*prodp
         nnr=2*p+1
         argq=argp/dble(nnr)/facpuis*dble(nn4-(2*nnr-1)*(2*nnr-1))
         
         pn=pn+argp
         qn=qn+argq
         if (cdabs(argp/pn).le.1.d-15.and. cdabs(argq/qn).le.1.d-15)
     $        goto 10
      enddo
      ierr=6
 10   continue
      jbes=(pn*cosz-qn*sinz)*facz
      return
      end
