      subroutine fonctiongreensurfcomp(hc,epsrel,epsabs,xs,ys,zs
     $     ,aretecube,k0,neps,nepsmax,dcouche,zcouche,epscouche,nbsphere
     $     ,nmax,n1m,nplan,ntp,nbs,nmat,nmatim,nplanm ,Tabzn,a,matind
     $     ,matindplan,matindice,matrange,nt)
      implicit none
      integer i,j,k,ii,nplanm,nplan
      integer neps,nepsmax,nbsphere,nmax,nmat,n1m,nbs,nmatim,nmati
     $     ,n1max,n2max,n1,n2,np1,np2,matindice(nplanm,nmatim)
     $     ,matind(0:2*n1m*n1m),matindplan(nplan,nplan),Tabzn(nmax),n11
     $     ,n22,ntp,nt
      double precision xs(nmax),ys(nmax),zs(nmax),aretecube,k0,xxp,yyp
     $     ,a(0:2*n1m*n1m),hc,epsrel,epsabs      
      double complex matrange(nbs,5),Ixx,Ixy,Ixz,Izx,Izz
      double precision dcouche(nepsmax),zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),stenseur(3,3)
c     write(*,*) 'couccou',n1m,nplan,nmatim,nepsmax,nmax,nbs
      double precision t1,t2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)     
      do i=0,2*n1m*n1m
         matind(i)=0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
c      write(*,*) 'couc nesp',neps,nbs
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

c     write(*,*) 'couc',nplan,nmatim
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      do i=1,nplanm
         do k=1,nmatim
            matindice(i,k)=0
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
c      write(*,*) 'couc'
      k=0
      do i=1,nplan
         do j=i,nplan
            k=k+1
            matindplan(i,j)=k
            matindplan(j,i)=k
         enddo
      enddo
c      write(*,*) 'couc'
      nmat=0
      nmati=0
      n1max=0
      n2max=0
      n2=0
      epsabs=0.d0
      call cpu_time(t1)
      do i=1,nbsphere
c         write(*,*) 'iter',i,nbsphere
c     write(*,*) n2,nmat,nmati,i
         np1=Tabzn(i) 
         do j=i,nbsphere
c            write(*,*) 'iter',j,nbsphere,xs(i),ys(j)
            np2=Tabzn(j)
            xxp=xs(i)-xs(j)
            yyp=ys(i)-ys(j)
            n11=dnint(dabs(xxp)/aretecube)
            n22=dnint(dabs(yyp)/aretecube)
c     write(*,*) 'nn',n11,n22
            n1max=max(n1max,n11)
            n2max=max(n2max,n22)
            n1=n11*n11+n22*n22
            n2=max(n2,n1)
            if (matind(n1).eq.0) then
               nmati=nmati+1
               matind(n1)=nmati
               a(n1)=dsqrt(xxp*xxp+yyp*yyp)
            endif
c            write(*,*) 'indeic',nmati,np1,np2,matind(n1),n1,nmat
c     $           ,matindice(matindplan(np2,np1),matind(n1))
c     a t on deja calculer le tenseur?
            if (matindice(matindplan(np2,np1),matind(n1)).ne.0) goto
     $           5003
c     non donc on calcule cet element
            nmat=nmat+1
            matindice(matindplan(np2,np1),matind(n1))=nmat
c            matindice(np1,np2,matind(n1))=nmat
c            write(*,*) 'matindice',matindice(matindplan(np2,np1)
c     $           ,matind(n1))
c            write(*,*) 'np',np1,np2,'nmat',nmat,'ij',i,j
c********************************************appelle tenseur**********  
c            write(*,*) hc,epsrel,epsabs,xs(i),ys(i) ,zs(i),xs(j),ys(j)
c     $           ,zs(j),k0,neps,dcouche,zcouche ,epscouche
            call tenseurmulticouchecomp(hc,epsrel,epsabs,xs(i),ys(i)
     $           ,zs(i),xs(j),ys(j),zs(j),k0,neps,dcouche,zcouche
     $           ,epscouche,Ixx,Ixy,Ixz,Izx,Izz)
c            write(*,*) 'ij',i,j,n1,matind(n1),nmat
c            write(*,*) 'green',Ixx,Ixy,Ixz,Izx,Izz
c            write(*,*) 'np1 np2',np1,np2

               matrange(nmat,1)=Ixx
               matrange(nmat,2)=Ixy
               matrange(nmat,3)=-Izx
               matrange(nmat,4)=Izz
               matrange(nmat,5)=-Ixz
            
c            endif
c            write(*,*) 'nmat',nmat
c     fin remplissage diago                              
 5003    enddo
      enddo
      call cpu_time(t2)
      write(*,*) 'temps',t2-t1
      a(0)=1.d300  
c**********************************************************
c     verifie dimension des tableaux
c**********************************************************
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
      
      if (n1max.gt.n1m) then
         write(*,*) 'trop d element lateral',n1max,n2max,n1m         
         stop
      endif         

      end
