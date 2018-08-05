      subroutine champlineairemicro(eps,disz,neps,nepsmax,x,y,k0,E0,ss
     $     ,pp,thetat,phit,infostr,nstop,Arx,Ary,Arz,Atx,Aty,Atz)
c     corrige bug dans le calcul de la composante y
c     cf agarwall
c     on eclaire par en dessous.
      implicit none
      integer i,j,neps,nepsmax,max,nmax,nstop,ntest
      parameter (nmax=40+2)
      double precision disz(0:nepsmax),x,y,phi,theta,phit,thetat ,k0
     $     ,k0d,pi,s,p,ss,pp
      double precision kx,ky,kp2,k02,ka
      double complex eps(0:nepsmax+1),E0,E0s,E0p,E0x,E0y,Exx,Eyy,E0z
     $     ,Arx,Ary,Arz,Atx,Aty,Atz,icomp,uncomp,Esca0,Eprod0,Esca(nmax)
     $     ,Eprod(nmax),matsca(nmax,nmax),matprod(nmax,nmax)
     $     ,bprod(nmax),bsca(nmax)
      double complex w0,wa,wb,wi(nmax)

      integer lda,info,ipvt(nmax),job
      double complex work(nmax),deterz(2)

      double  complex EE(nmax,3)
      character(64) infostr


      
      if (2*(nepsmax+1).gt.nmax) then
         nstop=1
         infostr='too many layers'
         return
      endif
      if (pp.lt.0.d0.or.pp.gt.1.d0) then
         nstop=1
         infostr='problem in plane incident wave'
         return
      endif
      p=dsqrt(pp)
      s=dsqrt(1.d0-pp)

      ntest=0
c     tet si milieu homogene
      do i=1,neps+1
         if (eps(i).ne.eps(0)) ntest=1
      enddo

      k0d=k0
      if (k0.eq.0.d0) k0=0.1d0

      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0

      if (theta.eq.0.d0) then
         ka=k0*cdsqrt(eps(0))
         k02=k0*k0
         wa=cdsqrt(k02*eps(0)*uncomp)
         wb=cdsqrt(k02*eps(neps+1)*uncomp)
         w0=cdsqrt(k02*uncomp)
         do i=1,neps
            wi(i)=cdsqrt(eps(i)*k02*uncomp)
         enddo
         if (neps.eq.0) wi(1)=wb

         bsca(1)=E0*cdexp(icomp*wa*disz(0))
         bsca(2)=-eps(0)*E0*wi(1)*cdexp(icomp*wa*disz(0))
         max=neps*2+2
         do i=3,max
            bsca(i)=0.d0
            bprod(i)=0.d0
         enddo         
         do i=1,max
            do j=1,max
               matsca(i,j)=0.d0
            enddo
         enddo
       
         matsca(1,1)=-cdexp(-icomp*wa*disz(0))
         matsca(1,2)=cdexp(icomp*wi(1)*disz(0))
         matsca(1,3)=cdexp(-icomp*wi(1)*disz(0))
         matsca(2,1)=-wi(1)*eps(0)*cdexp(-icomp*wa*disz(0))
         matsca(2,2)=-wa*eps(1)*cdexp(icomp*wi(1)*disz(0))
         matsca(2,3)=wa*eps(1)*cdexp(-icomp*wi(1)*disz(0))
         if (neps.eq.0) goto 900
         matsca(max,max-2)=-eps(neps)*wb*cdexp(icomp*wi(neps)*disz(neps)
     $        )
         matsca(max,max-1)=eps(neps)*wb*cdexp(-icomp*wi(neps)*disz(neps)
     $        )
         matsca(max,max)=eps(neps+1)*wi(neps)*cdexp(icomp*wb*disz(neps))
         matsca(max-1,max-2)=cdexp(icomp*wi(neps)*disz(neps))
         matsca(max-1,max-1)=cdexp(-icomp*wi(neps)*disz(neps))
         matsca(max-1,max)=-cdexp(icomp*wb*disz(neps))

         do i=3,max-2,2
            j=(i-1)/2
            matsca(i,i-1)=cdexp(icomp*wi(j)*disz(j))
            matsca(i,i)=cdexp(-icomp*wi(j)*disz(j))
            matsca(i,i+1)=-cdexp(icomp*wi(j+1)*disz(j))
            matsca(i,i+2)=-cdexp(-icomp*wi(j+1)*disz(j))
            matsca(i+1,i-1)=-eps(j)*wi(j+1)*cdexp(icomp*wi(j)*disz(j))
            matsca(i+1,i)=eps(j)*wi(j+1)*cdexp(-icomp*wi(j)*disz(j))
            matsca(i+1,i+1)=eps(j+1)*wi(j)*cdexp(icomp*wi(j+1)*disz(j))
            matsca(i+1,i+2)=-eps(j+1)*wi(j)*cdexp(-icomp*wi(j+1)*disz(j)
     $           )
         enddo

 900     job=11
         LDA=nmax
         call zgefa(matsca,lda,max,ipvt,info)
         if (INFO.ne.0) then
            nstop=1
            infostr='problem  in the incident multilayer 1'
            return
         endif
         call zgedi(matsca,lda,max,ipvt,deterz,work,job)      
         if (INFO.ne.0) then
            nstop=1
            infostr='problem  in the incident multilayer 2'
            return
         endif
         
         do i=1,max
            Esca(i)=0.d0
            do j=1,max
               Esca(i)=Esca(i)+matsca(i,j)*bsca(j)
            enddo
         enddo

         Exx=p*Esca(1)
         Eyy=s*Esca(1)
         Arx=Exx*dcos(phi)-Eyy*dsin(phi)
         Ary=Eyy*dcos(phi)+Exx*dsin(phi)
         Arz=0.d0

         Exx=p*Esca(max)
         Eyy=s*Esca(max)
         Atx=Exx*dcos(phi)-Eyy*dsin(phi)
         Aty=Eyy*dcos(phi)+Exx*dsin(phi)
         Atz=0.d0

         if (ntest.eq.0) then
            Arx=0.d0
            Ary=0.d0
            Arz=0.d0
         endif
         
         goto 999
      endif

      E0s=E0*s
      E0p=E0*p
      
      Exx=E0p*dcos(theta)
      Eyy=E0s
      E0x=Exx*dcos(phi)-Eyy*dsin(phi)
      E0y=Eyy*dcos(phi)+Exx*dsin(phi)
      E0z=-E0p*dsin(theta)

      ka=k0*cdsqrt(eps(0))
      kx=ka*dsin(theta)*dcos(phi)
      ky=ka*dsin(theta)*dsin(phi)
      kp2=kx*kx+ky*ky
      k02=k0*k0
      wa=cdsqrt((k02*eps(0)-kp2)*uncomp)
      wb=cdsqrt((k02*eps(neps+1)-kp2)*uncomp)
      w0=cdsqrt((k02-kp2)*uncomp)
      if (dimag(w0).lt.0.d0) w0=-w0
      if (dimag(wa).lt.0.d0) wa=-wa    
      if (dimag(wb).lt.0.d0) wb=-wb
      
      Esca0=kx*E0x+ky*E0y
      Eprod0=kx*E0y-ky*E0x      
      
      do i=1,neps
         wi(i)=cdsqrt(eps(i)*k02-kp2*uncomp)
         if (dimag(wi(i)).lt.0.d0) wi(i)=-wi(i)
      enddo
      if (neps.eq.0) wi(1)=wb

      bsca(1)=Esca0*cdexp(icomp*wa*disz(0))
      bsca(2)=-eps(0)*Esca0*wi(1)*cdexp(icomp*wa*disz(0))
      bprod(1)=Eprod0*cdexp(icomp*wa*disz(0))
      bprod(2)=Eprod0*cdexp(icomp*wa*disz(0))*(kp2+w0*wa)

      max=neps*2+2
      do i=3,max
         bsca(i)=0.d0
         bprod(i)=0.d0
      enddo

      do i=1,max
         do j=1,max
            matsca(i,j)=0.d0
            matprod(i,j)=0.d0
         enddo
      enddo

      matsca(1,1)=-cdexp(-icomp*wa*disz(0))
      matsca(1,2)=cdexp(icomp*wi(1)*disz(0))
      matsca(1,3)=cdexp(-icomp*wi(1)*disz(0))
      matsca(2,1)=-wi(1)*eps(0)*cdexp(-icomp*wa*disz(0))
      matsca(2,2)=-wa*eps(1)*cdexp(icomp*wi(1)*disz(0))
      matsca(2,3)=wa*eps(1)*cdexp(-icomp*wi(1)*disz(0))

      if (neps.eq.0) goto 901
      matsca(max,max-2)=-eps(neps)*wb*cdexp(icomp*wi(neps)*disz(neps))
      matsca(max,max-1)=eps(neps)*wb*cdexp(-icomp*wi(neps)*disz(neps))
      matsca(max,max)=eps(neps+1)*wi(neps)*cdexp(icomp*wb*disz(neps))
      matsca(max-1,max-2)=cdexp(icomp*wi(neps)*disz(neps))
      matsca(max-1,max-1)=cdexp(-icomp*wi(neps)*disz(neps))
      matsca(max-1,max)=-cdexp(icomp*wb*disz(neps))

      do i=3,max-2,2
         j=(i-1)/2
         matsca(i,i-1)=cdexp(icomp*wi(j)*disz(j))
         matsca(i,i)=cdexp(-icomp*wi(j)*disz(j))
         matsca(i,i+1)=-cdexp(icomp*wi(j+1)*disz(j))
         matsca(i,i+2)=-cdexp(-icomp*wi(j+1)*disz(j))
         matsca(i+1,i-1)=-eps(j)*wi(j+1)*cdexp(icomp*wi(j)*disz(j))
         matsca(i+1,i)=eps(j)*wi(j+1)*cdexp(-icomp*wi(j)*disz(j))
         matsca(i+1,i+1)=eps(j+1)*wi(j)*cdexp(icomp*wi(j+1)*disz(j))
         matsca(i+1,i+2)=-eps(j+1)*wi(j)*cdexp(-icomp*wi(j+1)*disz(j))

      enddo

 901  job=11
      LDA=nmax
      call zgefa(matsca,lda,max,ipvt,info)
      if (INFO.ne.0) then
         nstop=1
         infostr='problem  in the incident multilayer 3'
         return
      endif
      call zgedi(matsca,lda,max,ipvt,deterz,work,job)      
      if (INFO.ne.0) then
         nstop=1
         infostr='problem  in the incident multilayer 4'
         return
      endif

      do i=1,max
         Esca(i)=0.d0
         do j=1,max
            Esca(i)=Esca(i)+matsca(i,j)*bsca(j)
         enddo
      enddo

      matprod(1,1)=-cdexp(-icomp*wa*disz(0))
      matprod(1,2)=cdexp(icomp*wi(1)*disz(0))
      matprod(1,3)=cdexp(-icomp*wi(1)*disz(0))
      matprod(2,1)=-(kp2-w0*wa)*cdexp(-icomp*wa*disz(0))
      matprod(2,2)=(kp2+w0*wi(1))*cdexp(icomp*wi(1)*disz(0))
      matprod(2,3)=(kp2-w0*wi(1))*cdexp(-icomp*wi(1)*disz(0))

      if (neps.eq.0) goto 902
      matprod(max,max-2)=(kp2+w0*wi(neps))*cdexp(icomp*wi(neps)
     $     *disz(neps))
      matprod(max,max-1)=(kp2-w0*wi(neps))*cdexp(-icomp*wi(neps)
     $     *disz(neps))
      matprod(max,max)=-(kp2+w0*wb)*cdexp(icomp*wb*disz(neps))
      matprod(max-1,max-2)=cdexp(icomp*wi(neps)*disz(neps))
      matprod(max-1,max-1)=cdexp(-icomp*wi(neps)*disz(neps))
      matprod(max-1,max)=-cdexp(icomp*wb*disz(neps))

      do i=3,max-2,2
         j=(i-1)/2
         matprod(i,i-1)=cdexp(icomp*wi(j)*disz(j))
         matprod(i,i)=cdexp(-icomp*wi(j)*disz(j))
         matprod(i,i+1)=-cdexp(icomp*wi(j+1)*disz(j))
         matprod(i,i+2)=-cdexp(-icomp*wi(j+1)*disz(j))

         matprod(i+1,i-1)=(kp2+w0*wi(j))*cdexp(icomp*wi(j)*disz(j))
         matprod(i+1,i)=(kp2-w0*wi(j))*cdexp(-icomp*wi(j)*disz(j))
         matprod(i+1,i+1)=-(kp2+w0*wi(j+1))*cdexp(icomp*wi(j+1)*disz(j))
         matprod(i+1,i+2)=-(kp2-w0*wi(j+1))*cdexp(-icomp*wi(j+1)*disz(j)
     $        )

      enddo

 902  job=11
      LDA=nmax
      call zgefa(matprod,lda,max,ipvt,info)
      if (INFO.ne.0) then
         nstop=1
         infostr='problem  in the incident multilayer 5'
         return
      endif
      call zgedi(matprod,lda,max,ipvt,deterz,work,job)      
      if (INFO.ne.0) then
         nstop=1
         infostr='problem  in the incident multilayer 6'
         return
      endif

      do i=1,max
         Eprod(i)=0.d0
         do j=1,max
            Eprod(i)=Eprod(i)+matprod(i,j)*bprod(j)
         enddo
      enddo

      EE(1,1)=(kx*Esca(1)-ky*Eprod(1))/kp2
      EE(1,2)=(kx*Eprod(1)+ky*Esca(1))/kp2
      EE(1,3)=-(-kx*EE(1,1)-ky*EE(1,2))/wa
      
      EE(max,1)=(kx*Esca(max)-ky*Eprod(max))/kp2
      EE(max,2)=(kx*Eprod(max)+ky*Esca(max))/kp2
      EE(max,3)=(-kx*EE(max,1)-ky*EE(max,2))/wb
      
      do i=2,max-1,2
         j=i/2
         EE(i,1)=(kx*Esca(i)-ky*Eprod(i))/kp2
         EE(i,2)=(kx*Eprod(i)+ky*Esca(i))/kp2
         EE(i,3)=(-kx*EE(i,1)-ky*EE(i,2))/wi(j)
         EE(i+1,1)=(kx*Esca(i+1)-ky*Eprod(i+1))/kp2
         EE(i+1,2)=(kx*Eprod(i+1)+ky*Esca(i+1))/kp2
         EE(i+1,3)=-(-kx*EE(i+1,1)-ky*EE(i+1,2))/wi(j)
      enddo


      if (dimag(wa).ne.0.d0) then
         Arx=0.d0
         Ary=0.d0
         Arz=0.d0
      else
         Arx=EE(1,1)*cdexp(icomp*(kx*x+ky*y))
         Ary=EE(1,2)*cdexp(icomp*(kx*x+ky*y))
         Arz=EE(1,3)*cdexp(icomp*(kx*x+ky*y))
      endif

      if (ntest.eq.0) then
         Arx=0.d0
         Ary=0.d0
         Arz=0.d0
      endif


      
      if (dimag(wb).ne.0.d0) then
         Atx=0.d0
         Aty=0.d0
         Atz=0.d0
      else
         Atx=EE(max,1)*cdexp(icomp*(kx*x+ky*y))
         Aty=EE(max,2)*cdexp(icomp*(kx*x+ky*y))
         Atz=EE(max,3)*cdexp(icomp*(kx*x+ky*y))
      endif
      
 999  k0=k0d
      return
      end





