      subroutine champmultifft(eps,disz,neps,nepsmax,z,k0,E0x,E0y,E0z
     $     ,kx,ky,infostr,nstop,Ex,Ey,Ez,Arx,Ary,Arz,Atx,Aty,Atz)

c     Calcule le champ dans un multicouche pour une onde plane incidente
c     de vecteur d'onde k polarisee linéairement (ss ou pp).

c     definition la base p,s,k est un triede direct
c     s=z vectoriel k
c     p=s vectoriel k
c     k= sin(theta)cos(phi)x+sin(theta)sin(phi)y+cos(theta)z
c     on eclaire par en dessous dans le sens des z positifs.

c     E0 est l'amplitude de l'onde incidente, E est le champ total en
c     (0,0,z), Ar est l'amplitude de l'onde reflechie et At et
c     l'amplitude de l'onde transmise.

      implicit none
      integer i,j,neps,nepsmax,max,nmax,nstop
      parameter (nmax=40+2)
      double precision disz(0:nepsmax),z,k0,k0d ,pi
      double precision kx,ky,kp2,k02
      double complex eps(0:nepsmax+1),E0x,E0y,E0z,E0 ,icomp,uncomp,Esca0
     $     ,Eprod0,Esca(nmax),Eprod(nmax),matsca(nmax ,nmax)
     $     ,matprod(nmax,nmax) ,bprod(nmax),bsca(nmax)
      double complex w0,wa,wb,wi(nmax)

      integer lda,info,ipvt(nmax),job
      double complex work(nmax),deterz(2)

      double  complex EE(nmax,3),Ex,Ey,Ez,Arx,Ary ,Arz,Atx,Aty,Atz
      character(64) infostr
      
      if (2*(nepsmax+1).gt.nmax) then
         nstop=1
         infostr='too many layers'
         return
      endif
     
      k0d=k0
      if (k0.eq.0.d0) k0=0.1d0

c      write(*,*) 'ppp',eps,disz,neps,nepsmax,z,k0,E0x,E0y,E0z ,kx,ky

      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      kp2=kx*kx+ky*ky

      if (kp2.eq.0.d0) then
         k02=k0*k0
         wa=cdsqrt(k02*eps(0)*uncomp)
         wb=cdsqrt(k02*eps(neps+1)*uncomp)
         w0=cdsqrt(k02*uncomp)
         do i=1,neps
            wi(i)=cdsqrt(eps(i)*k02*uncomp)
         enddo
         if (neps.eq.0) wi(1)=wb

         bsca(1)=cdexp(icomp*wa*disz(0))
         bsca(2)=-eps(0)*wi(1)*cdexp(icomp*wa*disz(0))
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

         
         if (z.le.disz(0)) then
            E0=cdexp(icomp*wa*z)+Esca(1)*cdexp(-icomp*wa*z)
         endif

         if (z.gt.disz(neps)) then
            E0=Esca(max)*cdexp(icomp*wb*z)
         endif

         do i=0,neps-1
            if (z.gt.disz(i).and.z.le.disz(i+1)) then
               j=2*i+2
               E0=Esca(j)*cdexp(icomp*wi(i+1)*z)+Esca(j+1)*cdexp(-icomp
     $              *wi(i+1)*z)
            endif
         enddo    
        
         Ex=E0x*E0
         Ey=E0y*E0
         Ez=0.d0

         Arx=E0x*Esca(1)
         Ary=E0y*Esca(1)
         Arz=0.d0
         Atx=E0x*Esca(max)
         Aty=E0y*Esca(max)
         Atz=0.d0

         goto 999
      endif

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
      if (z.le.disz(0)) then
         Ex=E0x*cdexp(icomp*wa*z)+EE(1,1)*cdexp(-icomp*wa*z)
         Ey=E0y*cdexp(icomp*wa*z)+EE(1,2)*cdexp(-icomp*wa*z)
         Ez=E0z*cdexp(icomp*wa*z)+EE(1,3)*cdexp(-icomp*wa*z)
      endif

      if (z.gt.disz(neps)) then
         Ex=EE(max,1)*cdexp(icomp*wb*z)
         Ey=EE(max,2)*cdexp(icomp*wb*z)
         Ez=EE(max,3)*cdexp(icomp*wb*z)
      endif

      do i=0,neps-1
         if (z.gt.disz(i).and.z.le.disz(i+1)) then
            j=2*i+2
            Ex=EE(j,1)*cdexp(icomp*wi(i+1)*z)+EE(j+1,1)
     $           *cdexp(-icomp*wi(i+1)*z)
            Ey=EE(j,2)*cdexp(icomp*wi(i+1)*z)+EE(j+1,2)
     $           *cdexp(-icomp*wi(i+1)*z)
            Ez=EE(j,3)*cdexp(icomp*wi(i+1)*z)+EE(j+1,3)
     $           *cdexp(-icomp*wi(i+1)*z)
         endif
      enddo

      Arx=EE(1,1)
      Ary=EE(1,2)
      Arz=EE(1,3)

      Atx=EE(max,1)
      Aty=EE(max,2)
      Atz=EE(max,3)


 999  k0=k0d
      return
      end





