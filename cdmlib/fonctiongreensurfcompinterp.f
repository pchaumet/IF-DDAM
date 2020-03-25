      subroutine fonctiongreensurfcompinterp(hc,epsrel,epsabs,nx,ny,nz
     $     ,zs,aretecube,k0,neps,nepsmax,dcouche,zcouche,epscouche
     $     ,nbsphere,nmax,n1m,nplan,ntp,nbs,nmat,nmatim,nplanm,Tabzn,a
     $     ,matind ,matindplan,matindice,matrange,ninter,ninterp,nt)
      implicit none
      integer i,j,k,ii,kav
      integer neps,nepsmax,nbsphere,nmax,nmat,n1m,nplan,nplanm,nbs
     $     ,nmatim,nmati,n1max,n2max,np1,np2,matindice(nplanm,nmatim)
     $     ,matind(0:2*n1m*n1m),matindplan(nplan,nplan) ,Tabzn(nmax),n11
     $     ,n22,ntp,nt,nx,ny,nz,ninter,ninterp
      double precision aretecube,k0,xxp,yyp,a(0:2*n1m*n1m),hc,epsrel
     $     ,epsabs,x0,y0,aretecubeint,epsmax,zs(nmax),pi,zz(ntp)
      double complex matrange(nbs,5),Ixx,Ixy,Ixz,Izx,Izz
      double precision dcouche(nepsmax),zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),stenseur(3,3)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)  
      do i=0,2*n1m*n1m
         matind(i)=0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)       
      do i=1,nbs
         matrange(i,1)=0.0d0
         matrange(i,2)=0.0d0
         matrange(i,3)=0.0d0
         matrange(i,4)=0.0d0
         matrange(i,5)=0.0d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)      
      do i=1,nplanm
         do k=1,nmatim
            matindice(i,k)=0
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
c      write(*,*) 'couc3'
      k=0
      do i=1,nplan
         do j=i,nplan
            k=k+1
            matindplan(i,j)=k
            matindplan(j,i)=k
         enddo
      enddo
      zz(1)=zs(1)
      kav=1
      do i=2,nbsphere
         k=tabzn(i)
         if (k.eq.kav) goto 123
         zz(k)=zs(i)
         kav=k
 123  enddo


      nmat=0
      nmati=0
      epsabs=0.d0

      pi=dacos(-1.d0)
      x0=0.d0
      y0=0.d0
      aretecubeint=aretecube/dble(ninterp)
      ninter=ceiling(dsqrt(dble(nx*nx+ny*ny)))*ninterp
      
      do k=0,ninter
         nmati=nmati+1
         matind(k)=nmati
         a(k)=aretecubeint*dble(k)
      enddo

      do np1=1,nz
         do np2=np1,nz

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,xxp,nmat,Ixx,Ixy,Ixz,Izx,Izz)   
!$OMP DO SCHEDULE(STATIC)             
            do k=0,ninter
               xxp=a(k)             
               nmat=matind(k)+(matindplan(np2,np1)-1)*nmati
               matindice(matindplan(np2,np1),matind(k))=nmat

               call tenseurmulticouchecomp(hc,epsrel,epsabs,x0,y0
     $              ,zz(np1),xxp,y0,zz(np2),k0,neps,dcouche,zcouche
     $              ,epscouche,Ixx,Ixy,Ixz,Izx,Izz)

c               write(*,*) k,matind(k),matindplan(np2,np1),nmati,nmat,Ixx
c     $              ,np1,np2                  
               matrange(nmat,1)=Ixx
               matrange(nmat,2)=Ixy
               matrange(nmat,3)=-Izx
               matrange(nmat,4)=Izz
               matrange(nmat,5)=-Ixz
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                     
         enddo
      enddo
      a(0)=1.d300  
c**********************************************************
c     verifie dimension des tableaux
c**********************************************************
      nmat=matindplan(nplan,nplan)*nmati
      write(*,*) 'nbsphere',nbsphere,' nmat',nmat,' nmati',nmati,'ntp'
     $     ,ntp
      write(*,*) 'nmax',nmax,' nbs',nbs,' nmatim',nmatim,'ntp',ntp
      if (nbsphere.gt.nmax) then
         write(*,*) 'decoupe trop la sphere',nbsphere,nmax
         stop
      endif

      if (nmat.gt.nbs) then
         write(*,*) 'pas assez de place pour propa',nmat,nbs
         stop
      endif
      
      if (ntp.gt.nplan) then
         write(*,*) 'trop de plan',ntp,nplan
         stop
      endif
      
      if (nmati.gt.nmatim) then
         write(*,*) 'grandir matindice',nmati,nmatim
         stop
      endif
      
      end
