      subroutine Jhankel(z,norder,nseq,nspec,ierr,nseqlda,Hankel)
      implicit none

      integer norder,nseq,ierr,nseqlda,nspec
      double complex z,Hankel(nseqlda)

      integer nref,i,nmax,j
      double precision zm
      double complex Jn1,Jn,hank,Ji

      DOUBLE PRECISION ZR, ZI
      double complex JB0,icomp

      icomp=(0.d0,1.d0)
      ierr=0
      zm=cdabs(z)
      zr=dreal(z)
      zi=dimag(z)
      nmax=norder+nseq-1
      if (nseq.gt.nseqlda) ierr=1
      if (norder.lt.0) ierr=2
      if (nseq.le.0) ierr=3
      if (nspec.ne.1.and.nspec.ne.2) ierr=4
      if (zm.eq.0) ierr=5
      if (ierr.ne.0) return

           
c      if (zr*zr/70.d0+zi*zi/40.d0.ge.1.d0) then
         if (zr*zr/64.d0+zi*zi/36.d0.ge.1.d0) then
c     Calcul Hankel pour de grands arguments
         if (nseq.eq.1) then
            call hankelpl(nmax,z,nspec,hank,ierr)
            Hankel(1)=hank
         elseif (nseq.eq.2) then
            call hankelpl(nmax,z,nspec,hank,ierr)
            Hankel(2)=hank
            call hankelpl(nmax-1,z,nspec,hank,ierr)
            Hankel(1)=hank
         else
            call hankelpl(nmax,z,nspec,hank,ierr)
            Hankel(nseq)=hank
            call hankelpl(nmax-1,z,nspec,hank,ierr)
            Hankel(nseq-1)=hank
            do i=nmax-2,norder,-1
               j=i-nmax+nseq
               Hankel(j)=dble(2*(i+1))/z*Hankel(j+1)-Hankel(j+2)
            enddo
         endif
c     Fin calcul Bessel pour de grands arguments
      else
c     Calcul Bessel pour  arguments  petits
         if (nseq.eq.1) then
            call hankelp(nmax,z,nspec,hank,ierr)
            Hankel(1)=hank
         elseif (nseq.eq.2) then
            call hankelp(nmax,z,nspec,hank,ierr)
            Hankel(2)=hank
            call hankelp(nmax-1,z,nspec,hank,ierr)
            Hankel(1)=hank
         else
            call hankelp(nmax,z,nspec,hank,ierr)
            Hankel(nseq)=hank
            call hankelp(nmax-1,z,nspec,hank,ierr)
            Hankel(nseq-1)=hank
            do i=nmax-2,norder,-1
               j=i-nmax+nseq
               Hankel(j)=dble(2*(i+1))/z*Hankel(j+1)-Hankel(j+2)
            enddo
         endif
         return
      endif

      end
c     ***********************************************************
      subroutine hankelp(n,z,nspec,hank,ierr)
      implicit none
      integer i,n,p,nspec,ierr
      double complex z,hank,arg,arga,facz,icomp,const1,sommet,sommef
     $     ,csom1,csom2
      double precision facnp,expx,facp,som1,som2,euler,pi,sig

      euler=0.5772156649015328606065120d0
      pi=3.14159265358979323846d0
      icomp=(0.d0,1.d0)
      if (nspec.eq.1) sig=1.d0
      if (nspec.eq.2) sig=-1.d0
      
      
      facnp=1.d0
      som1=0.d0
      facz=-z*z/4.d0
      const1=1.d0+sig*icomp/pi*2.d0*(cdlog(z/2.d0)+euler)
      do i=2,n
         facnp=facnp*dble(i)
      enddo
      sommef=0.d0
      if (n.ge.2) then

         csom2=0.d0
         arga=facnp/dble(n)
         csom1=arga
         do i=1,n-2
            arga=arga*(z*z/4.d0)/dble(i*(n-i))
            csom1=csom1+arga
         enddo
         csom1=csom1+arga*(z*z/4.d0)/dble(n-1)
         csom1=csom1*(z/2.d0)**dble(-n)
         do i=1,n
            csom2=csom2+1.d0/dble(i)
         enddo
         csom2=(z/2.d0)**dble(n)*csom2/facnp
         sommef=-sig*icomp/pi*(csom1+csom2)
      elseif (n.eq.1) then
         sommef=-sig*icomp/pi*(2.d0/z+z/2.d0)
      endif

      som1=0.d0
      som2=0.d0
      do i=1,n
         som1=som1+1.d0/dble(i)
      enddo
      if (n.ge.1) then
         sommef=sommef+sig*icomp/pi*(z/2.d0)**dble(n)/facnp*som1
      endif

      arga=(z/2.d0)**dble(n)/facnp
      hank=arga*(const1-sig*icomp/pi*som1)
      do p=1,1000
         som1=som1+1.d0/dble(n+p)
         som2=som2+1.d0/dble(p)         
         arga=arga/dble(p*(n+p))*facz
         sommet=arga*(const1-sig*icomp/pi*(som1+som2))
         if (cdabs(sommet).le.cdabs(hank)*1.d-15) goto 10
         hank=hank+sommet
      enddo
      ierr=5
 10   continue

      hank=hank+sommef
      return
      end
c     ***********************************************
      subroutine hankelpl(n,z,nspec,hank,ierr)
      implicit none
      integer i,n,p,nr,nnr,nspec,ierr,k,kmax
      double complex z,hank,facz,argexp,exph,argp1,zz,c0,c1,c2,c3
      double complex pn,qn,icomp,facpuis,argp,argq,prodp,ctmp,sszz,ss2
     $     ,ss3,ss4,zz2,zz3,zz4
      double precision fac,pi,sig,zm,ss,dnr,dn

      hank=0.d0
      pi=3.14159265358979323846d0
      icomp=(0.d0,1.d0)
      if (nspec.eq.1) sig=1.d0
      if (nspec.eq.2) sig=-1.d0
      facz=cdsqrt(2.d0/pi/z)
      dn=dble(n)
      argexp=z-dn*pi/2.d0-pi/4.d0
      exph=cdexp(icomp*sig*argexp)

      facpuis=8.d0*z
      pn=(1.d0,0.d0)
      prodp=(1.d0,0.d0)
      nr=0
      dnr=0.d0
      zm=cdabs(z)
      kmax=dnint((zm+dsqrt(zm*zm+dble(n*n)))/2.d0)+1
      ctmp=icomp*sig/facpuis
      do k=1,kmax
         nr=nr+1
         dnr=dnr+1.d0
         prodp=prodp*dble((4*n*n-(2*nr-1)*(2*nr-1)))/dnr/facpuis
         argp=(sig*icomp)**dnr*prodp
         nr=nr+1
         dnr=dnr+1.d0
         prodp=prodp*dble((4*n*n-(2*nr-1)*(2*nr-1)))/dnr/facpuis
         argp1=(sig*icomp)**dnr*prodp
         pn=pn+argp+argp1
         if (cdabs(argp+argp1).le.cdabs(pn)*1.d-15) goto 10
      enddo

    
      nr=nr+1
      prodp=prodp*dble((4*n*n-(2*nr-1)*(2*nr-1)))/dble(nr)/facpuis
      argp1=(sig*icomp)**dble(nr)*prodp
      ss=dble(nr-n)-0.5d0
      zz=-2.d0*icomp*z*sig
      sszz=zz+ss
      ss2=ss*ss
      zz2=zz*zz
      c0=zz/sszz*(1.d0-zz/sszz**2-zz*(ss-2.d0*zz)/sszz**4-zz *(ss2-8.d0
     $     *ss*zz+6.d0*zz2)/sszz**6)-zz*(ss2*ss-22.d0 *ss2*zz+58.d0 *ss
     $     *zz2-24.d0*zz2*zz)/sszz**8
      c1=((sszz+1.d0)*c0-zz)/zz
      c2=(sszz*c1+c0-1.d0)/(2.d0*zz)
      c3=((sszz-1.d0)*c2+c1)/(3.d0*zz)
      dnr=dble(nr+n)
      c1=-(dn-0.5d0)*zz/(dnr-0.5d0)*c1
      c2=((dn-0.5d0)*(dn-1.5d0)*zz2/(dnr-0.5d0) /(dnr-1.5d0)) *c2
      c3=-(dn-0.5d0)*(dn-1.5d0)*(dn-2.5d0)*zz2*zz /(dnr-0.5d0) /(dnr
     $     -1.5d0)/(dnr-2.5d0)*c3
      argp1=argp1*(c0+c1+c2+c3)
      pn=pn+argp1

 10   continue
      hank=exph*pn*facz
      return
      end

