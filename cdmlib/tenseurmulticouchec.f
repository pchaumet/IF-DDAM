      subroutine tenseurmulticouchecomp(hc,epsrel,epsabs,x,y,z,xa,ya,za
     $     ,k0,neps,dcouche,zcouche,epscouche,Ixx,Ixy,Ixz,Izx,Izz)
      implicit none
      integer i,nepsmax,neps,ntest
      parameter (nepsmax=20)
      double precision dcouche(nepsmax),zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),epsref
      double complex Ixx,Ixy,Ixz,Izx,Izz
      double precision x,y,z,xa,ya,za,a,k0,pi,xxp,yyp,hc,epsrel,epsabs
      double complex icomp,icomppi,uncomp,zero
      common/datalog/icomp,icomppi,uncomp,zero
!$OMP THREADPRIVATE(/datalog/)
      icomp=(0.d0,1.d0)
      pi=dacos(-1.d0)
      icomppi=icomp*pi
      uncomp=(1.d0,0.d0)
      zero=(-1.d20,0.d0)

c      write(*,*) 'multi',x,y,z,xa,ya,za,k0,neps
c      write(*,*) 'couche',dcouche,zcouche,epscouche



      xxp=x-xa
      yyp=y-ya
      a=dsqrt(xxp*xxp+yyp*yyp)

      epsref=epscouche(0)
      ntest=0
c     tet si milieu homogene
      do i=1,neps+1
         if (epscouche(i).ne.epsref) ntest=1
      enddo
      if (ntest.eq.0) then
         call propesplibI(a,z,za,k0,epsref,Ixx,Ixy,Ixz,Izz)
         Ixx=Ixx/epsref
         Ixy=Ixy/epsref
         Ixz=Ixz/epsref
         Izz=Izz/epsref
         Izx=Ixz
         return
      endif

      if (a.eq.0.d0) then
         Ixy=0.d0
         Izx=0.d0
         Ixz=0.d0
         call Idiagmulticouchecomp(hc,epsrel,epsabs,k0,z,za,neps
     $        ,dcouche,zcouche,epscouche,Ixx,Izz)
      else
         call Imulticouchecomp(hc,epsrel,epsabs,k0,a,z,za,neps,dcouche
     $        ,zcouche,epscouche,Ixx,Ixy,Ixz,Izx,Izz)
c         write(*,*) 'rr',Ixx,Ixy,Ixz,Izx,Izz
      endif
      end
c******************************************************************
c******************************************************************
c******************************************************************
      subroutine Idiagmulticouchecomp(hc,epsrel,epsabs,k0,z,za,neps
     $     ,dcouche,zcouche,epscouche,Ixx,Izz)
      implicit none
      integer nepsmax,neps,testloin,neval
      parameter (nepsmax=20)
      double precision dcouche(nepsmax),zcouche(0:nepsmax),epsmax,k0,z
     $     ,za,a,alpham,hc
      double complex epscouche(0:nepsmax+1)
      double complex Ixx,Izz,Txx,Tzz

c     declaration lie a l'integration
      integer nnn
      double precision pi,aa,bb,kmax,hkmax     

c     declaration des variables du programme
      integer i,j,nbint

      integer nlda,ishanks,lordref
      parameter (nlda=5)
      double complex Ielliptique(nlda),Iinfini(nlda)
      double precision Imin,epsabs,epsrel,erreur(nlda)

      double complex icomp
      integer nnepsmax,nneps,nc,no
      double precision k02,zz,zza,zzcouche(0:nepsmax),zmax,zmin

      double complex eepscouche(0:nepsmax+1)
      
      common/donneemultin/nneps,nc,no
      common/donneemulti/k02,kmax,hkmax,a,zz,zza
      common/zcoucheroutine/zzcouche
      common/epscouroutine/eepscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      external Idiagmultiecomp,Idiagmultiicomp,Idiagmultiedessouscomp
     $     ,Idiagmultiidessouscomp,Idiagmultiedessuscomp
     $     ,Idiagmultiidessuscomp ,Idiagmultiilogcomp

      
c     initialisation des donnees
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      nbint=2
      ishanks=0

      do i=1,nbint
         Ielliptique(i)=0.d0
         Iinfini(i)=0.d0
      enddo

c     repere dans quelle couche est le dipole
      nc=0
      do i=0,neps
         if (za.ge.zcouche(i)) then
            nc=nc+1
         else
            goto 10
         endif
      enddo
 10   no=0
 
c     repere dans quelle couche est l'observation
      do i=0,neps
         if (z.ge.zcouche(i)) then
            no=no+1
         else
            goto 20
         endif
      enddo

c      write(*,*) 'numero couche dipole',nc,no,zcouche,'zz',z,za

 20   if (no.eq.0.and.nc.eq.0) then
         zmin=dabs(z+za-2.d0*zcouche(0))
         goto 30
      elseif (no.eq.neps+1.and.nc.eq.neps+1) then
         zmin=dabs(z+za-2.d0*zcouche(neps))
         goto 30
      endif
c      write(*,*) 'numero couche observation',no,z,neps
c     calcul de zmax
c$$$      if(no.le.nc) then
c$$$         if (no.eq.0) then
c$$$            zmin=min(dabs(z-za),-z-za+2.d0*zcouche(nc))
c$$$            zmax=zmin
c$$$            write(*,*) 'coucou1',-z+za,-z-za+2.d0*zcouche(nc)
c$$$         else
c$$$            if (no.eq.nc) then
c$$$               zmin=-z-za+2.d0*zcouche(nc)
c$$$               write(*,*) 'zzz',zmin
c$$$            else
c$$$               zmin=min(dabs(z-za),-z-za+2.d0*zcouche(nc))
c$$$            endif
c$$$            zmax=zmin
c$$$
c$$$            zmin=min(zmin,z+za-2.d0*zcouche(no-1))
c$$$            zmax=max(zmax,z+za-2.d0*zcouche(no-1))
c$$$
c$$$            zmin=min(zmin,dabs(z-za)-2.d0*zcouche(no-1)+2.d0*zcouche(nc)
c$$$     $           )
c$$$            zmax=max(zmax,dabs(z-za)-2.d0*zcouche(no-1)+2.d0*zcouche(nc)
c$$$     $           )
c$$$c     write(*,*) 'coucou4',dabs(z-za)-2.d0*zcouche(no-1)+2.d0
c$$$c     $           *zcouche(nc),zmin,zmax
c$$$         endif
c$$$      elseif(no.gt.nc) then
c$$$         if (no.eq.neps+1) then
c$$$            zmin=min(dabs(z-za),z+za-2.d0*zcouche(nc-1))
c$$$            zmax=zmin
c$$$            write(*,*) 'coucou10',-za+z,z+za-2.d0*zcouche(nc-1)
c$$$         else
c$$$            zmin=min(dabs(z-za),z+za-2.d0*zcouche(nc-1))
c$$$            zmax=zmin
c$$$            write(*,*) 'coucou2',-za+z,z+za-2.d0*zcouche(nc)
c$$$            zmin=min(zmin,-z-za+2.d0*zcouche(no))
c$$$            zmax=max(zmax,-z-za+2.d0*zcouche(no))
c$$$            write(*,*) 'coucou3',-z-za+2.d0*zcouche(no)
c$$$            zmin=min(zmin,-z+za+2.d0*zcouche(no)-2.d0*zcouche(nc-1))
c$$$            zmax=max(zmax,-z+za+2.d0*zcouche(no)-2.d0*zcouche(nc-1))
c$$$            write(*,*) 'coucou4',-z+za+2.d0*zcouche(no)-2.d0*zcouche(nc
c$$$     $           -1),zmin,zmax
c$$$         endif
c$$$      endif
      if (nc.eq.no) then
         zmin=min(dabs(z+za-2.d0*zcouche(nc-1)),dabs(-z-za+2.d0
     $        *zcouche(nc)))
         goto 30
      elseif (nc.lt.no) then
         zmin=dabs(za-zcouche(nc))
         goto 30
      elseif (nc.gt.no) then
         zmin=dabs(za-zcouche(nc-1))
         goto 30
      endif
      
c     repette les donnees a cause du common de merde
 30   k02=k0*k0
      a=0.d0
      zz=z
      zza=za
      nnepsmax=nepsmax
      nneps=neps
      zzcouche(0)=zcouche(0)
      eepscouche(0)=epscouche(0)
      do i=1,neps
         zzcouche(i)=zcouche(i)
         eepscouche(i)=epscouche(i)
      enddo
      eepscouche(neps+1)=epscouche(neps+1)

c     commence l'integration
      epsmax=0.d0
      do i=0,neps+1
         epsmax=max(epsmax,dsqrt(dabs(dreal(epscouche(i)))))
      enddo
c     kmax pour le chemin elliptique +k0 pour plus de surete
      kmax=k0*epsmax+k0/2.d0
c     kmax=1.5d0*k0
c     signe moins parce qu'on passe par la partie imaginaire negative
c     pour eviter les poles
      hkmax=-hc*kmax

      if (no.eq.nc) then
         call propesplibdiagI(z,za,k0,epscouche(nc),Txx,Tzz)
         EPSABS=max(epsabs,min(cdabs(Txx),cdabs(Tzz))*EPSREL)
      endif
c      write(*,*) 'EPSABS',no,nc,epsrel,EPSABS


c     PREMIERE INTEGRALE DE 0 A KMAX AVEC CHEMIN ELLIPTIQUE
      aa=-pi/2.d0
      bb=pi/2.d0
      nnn=2
      
      if (nc.eq.0) then        
c         write(*,*) 'coucou2 e'

         call intgausskronrodpattersonmulti(aa,bb,ishanks,Ielliptique
     $        ,erreur,epsrel,EPSABS,lordref,Idiagmultiedessouscomp,nlda
     $        ,nnn)
c         write(*,*) 'dessous',Ielliptique
      elseif (nc.eq.neps+1) then

         call intgausskronrodpattersonmulti(aa,bb,ishanks,Ielliptique
     $        ,erreur,epsrel,EPSABS,lordref,Idiagmultiedessuscomp,nlda
     $        ,nnn)
c         write(*,*) 'dessus',Ielliptique
       
      else

         call intgausskronrodpattersonmulti(aa,bb,ishanks,Ielliptique
     $        ,erreur,epsrel,EPSABS,lordref,Idiagmultiecomp,nlda,nnn)
c         write(*,*) 'milieu',Ielliptique
      endif
      NEVAL=2**lordref-1
c      write(*,*) 'Ielliptique',Ielliptique,NEVAL
c      write(99,*) NEVAL
c     DEUXIEME INTEGRALE DE KMAX A L'INFINI AVEC CHEMIN EN DROITE LIGNE
c     defini la borne superieur de l'integrale par raaport au mini de
c     l'exponentiel.

      Imin=min(cdabs(Ielliptique(1)),cdabs(Ielliptique(2)))
      epsabs=max(epsabs,Imin*epsrel)
      
c     zmin plus grand que 2 lambda: plus d'onde evanescente
      if (zmin.ge.12.d0/k0) then
         Iinfini=0.d0
         NEVAL=0
c         write(*,*) 'evite infini',zmin
         goto 40
      endif
c      write(*,*) 'comp',zmin,kmax,-dlog(min(epsrel,epsabs/Imin))/zmin
      alpham=-dlog(min(epsrel,epsabs/Imin))/zmin+kmax
      if (zcouche(neps).ge.350.d0/alpham) then
c         write(*,*) 'depasse les capacites machine'
c         write(*,*) 'zcouche',zcouche(neps),alpham,350.d0/alpham
c         write(*,*) 'zmin alpham',zmin,alpham
         alpham=350.d0/zcouche(neps)
c         write(*,*) 'alpham possible',alpham
      endif
      
c      write(*,*) 'alpham',kmax,alpham,kmax/k0,alpham/k0      
      aa=kmax
      bb=alpham
      testloin=0
      if (nc.eq.0) then
         
         
         call Idiagmultiidessouscomp(aa,nnn,nlda,Iinfini)                  
         do i=1,nnn
            if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin =testloin
     $           +1
            Iinfini(i)=0.d0
         enddo
         NEVAL=1
c         write(*,*) 'testloin',testloin,nnn
         if (testloin.eq.nnn) goto 40

         call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $        ,erreur,epsrel,EPSABS,lordref,Idiagmultiidessouscomp,nlda
     $        ,nnn)
         NEVAL=2**lordref-1
      elseif (nc.eq.neps+1) then

         call Idiagmultiidessuscomp(aa,nnn,nlda,Iinfini)         
         do i=1,nnn
            if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin =testloin
     $           +1
            Iinfini(i)=0.d0
         enddo
         NEVAL=1
c         write(*,*) 'testloin',Iinfini,testloin,nnn,aa,bb
         if (testloin.eq.nnn) goto 40
         call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $        ,erreur,epsrel,EPSABS,lordref,Idiagmultiidessuscomp,nlda
     $        ,nnn)
         NEVAL=2**lordref-1
      else
c     write(*,*) 'coucou'

         call Idiagmultiicomp(aa,nnn,nlda,Iinfini)         
         do i=1,nnn
            if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin =testloin
     $           +1
            Iinfini(i)=0.d0
         enddo
         NEVAL=1

        

         if (testloin.eq.nnn) goto 40 

         if (no.eq.nc) then

            call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $           ,erreur,epsrel,EPSABS,lordref,Idiagmultiilogcomp,nlda
     $           ,nnn)
            NEVAL=2**lordref-1
         else

            call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $           ,erreur,epsrel,EPSABS,lordref,Idiagmultiicomp,nlda
     $           ,nnn)
            NEVAL=2**lordref-1
         endif
      endif
c      write(*,*) 'Iinfini',Iinfini,NEVAL
      
 40   Ixx=(Ielliptique(1)+Iinfini(1))*0.5d0*icomp
      Izz=icomp*(Ielliptique(2)+Iinfini(2))
c     write(*,*) Ixx,Izz
      if (no.eq.nc) then
c         write(*,*) 'propa',Txx,Tzz,'surf',Ixx,Izz
         Ixx=Ixx+Txx
         Izz=Izz+Tzz
      endif
      write(98,*) NEVAL
c     divise par epsilon de la couche ou il y a le dipole
      Ixx=Ixx/epscouche(nc)
      Izz=Izz/epscouche(nc)


      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Idiagmultiecomp(theta,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision theta
      double complex S11,P11,S21,P21,S12,S22,icomp
      double complex k,k2,wronskien,integrant(nlda)

      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)

c      write(*,*) 'routine ancienne',k02,kmax,hkmax,a,z,za,neps,nc,no     

c     calcul de k complex avec theta
      k=kmax/2.d0*(dsin(theta)+1.d0)+icomp*hkmax*dcos(theta)
      k2=k*k
c      write(*,*) 'kkk',k,k2
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
c         write(*,*) 'eplus',e_plus,icomp*(w(i)+w(i+1))*zcouche(i),w(i)
c     $        ,w(i+1),zcouche(i)
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp
c         write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c         write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire         
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca1',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod1',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no.and.no.ne.nc) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo
c      write(*,*) 'obs1',mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
c     $     ,mat_sca_obs22,mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
c     $     ,mat_prod_obs22
c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-mat_sca_sous22
c      write(*,*) 'rr1',mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      mat_prod_11=mat_prod_sur11
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      det=mat_sca_11*mat_sca_22-mat_sca_12*mat_sca_21
c      write(*,*) 'det',det,mat_sca_11,mat_sca_22,mat_sca_12
c     $     ,mat_sca_21
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22/det
      mat_sca_22=ctmp/det
      mat_sca_12=-mat_sca_12/det
      mat_sca_21=-mat_sca_21/det

      ctmp=mat_sca_11*B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=mat_sca_11*B_sca_12+mat_sca_12*B_sca_22
      mat_sca_11=ctmp      
      ctmp=mat_sca_21*B_sca_11+mat_sca_22*B_sca_21
      mat_sca_22=mat_sca_21*B_sca_12+mat_sca_22*B_sca_22
      mat_sca_21=ctmp
c      write(*,*) 'tt1',mat_sca_11,mat_sca_22, mat_sca_12,mat_sca_21
      det=mat_prod_11*mat_prod_22-mat_prod_12*mat_prod_21
c      write(*,*) 'det',det
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22/det
      mat_prod_22=ctmp/det
      mat_prod_12=-mat_prod_12/det
      mat_prod_21=-mat_prod_21/det

      mat_prod_21=mat_prod_21*B_prod_11+mat_prod_22*B_prod_21
      mat_prod_11=mat_prod_11*B_prod_11+mat_prod_12*B_prod_21
      mat_prod_22=0.d0
      mat_prod_12=0.d0
c      write(*,*) 'toto1',mat_prod_11,mat_prod_22,mat_prod_21,mat_prod_21
c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp

      elseif (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c         write(*,*) 'dessus',S11,S12,S21,S22,P11,P21,ctmp,mat_sca_11
c     $        ,mat_sca_12
      else
c     le point d'observation est dans le multicouche
         if (no.lt.nc) then
c     en dessous du dipole
c            write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22

            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs22*mat_sca_21*ctmp
            S22=mat_sca_obs22*mat_sca_22*ctmp
            P21=mat_prod_obs22*mat_prod_21*ctmp
         elseif (no.gt.nc) then
c     au dessus du dipole
c            write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs11*mat_sca_11*ctmp
            P11=mat_prod_obs11*mat_prod_11*ctmp
            S12=mat_sca_obs11*mat_sca_12*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs21*mat_sca_11*ctmp
            S22=mat_sca_obs21*mat_sca_12*ctmp
            P21=mat_prod_obs21*mat_prod_11*ctmp
c            write(*,*) 'S2222',S22
         elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c            write(*,*) 'identique'
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
            S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c            write(*,*) 'S22',mat_sca_obs22*mat_sca_21,B_sca_21
c     $           ,mat_sca_obs22*mat_sca_22,B_sca_22,mat_prod_obs22
c     $           *mat_prod_21,B_prod_21
            P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c     avec l'espace libre
c            S21=mat_sca_obs22*mat_sca_21*ctmp
c            S22=mat_sca_obs22*mat_sca_22*ctmp
c            P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c            ctmp=cdexp(icomp*w(no)*z)
c            S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c            S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c            P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
c            ctmp=cdexp(-icomp*w(no)*z)
c            S21=(mat_sca_obs21*mat_sca_11)*ctmp
c            S22=(mat_sca_obs21*mat_sca_12)*ctmp
c            P21=(mat_prod_obs21*mat_prod_11)*ctmp
         endif
      endif


c     calcul element differentiel
      wronskien=kmax/2.d0*dcos(theta)-icomp*hkmax*dsin(theta)
c      write(*,*) 'wr',wronskien
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'SSS',S11,P11,S21,P21,S12,S22
      integrant(1)=k*(S11+P11+S21+P21)*wronskien
      integrant(2)=(-S12+S22)*wronskien*k/w(no)

c      write(*,*) 'theta',theta
c      write(*,*)  integrant(1),integrant(2),integrant(3),integrant(4)
c      write(*,*) '*********************************'
c      write(920,*) (theta),dreal(integrant(1))
c      write(921,*) theta,dimag(integrant(1))
c      write(922,*) (theta),dreal(integrant(2))
c      write(923,*) (theta),dimag(integrant(2))
     
      end


c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Idiagmultiicomp(k,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision k,k2
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      k2=k*k
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
c         w(i)=cdsqrt(-k2*(1.d0,0.d0))
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp

c         write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c         write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c         write(*,*) 'prod',mat_prod_11,mat_prod_12 ,mat_prod_21
c     $        ,mat_prod_22
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22
c     evite le cas no=nc fait au dessus
         if (i.eq.no.and.no.ne.nc) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c      write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      det=mat_sca_11*mat_sca_22-mat_sca_12*mat_sca_21
c      write(*,*) 'det',det,mat_sca_11,mat_sca_22,mat_sca_12,mat_sca_21
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22/det
      mat_sca_22=ctmp/det
      mat_sca_12=-mat_sca_12/det
      mat_sca_21=-mat_sca_21/det
      
c      write(*,*) 'inverse s',mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
c      write(*,*) 'chiant1',mat_sca_12

      ctmp=mat_sca_11*B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=mat_sca_11*B_sca_12+mat_sca_12*B_sca_22
      mat_sca_11=ctmp      
      ctmp=mat_sca_21*B_sca_11+mat_sca_22*B_sca_21
      mat_sca_22=mat_sca_21*B_sca_12+mat_sca_22*B_sca_22
      mat_sca_21=ctmp

c      write(*,*) 'chiant2',mat_sca_12,mat_sca_11,B_sca_12,mat_sca_12
c     $     ,B_sca_22

      det=mat_prod_11*mat_prod_22-mat_prod_12*mat_prod_21
c      write(*,*) 'det',det
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22/det
      mat_prod_22=ctmp/det
      mat_prod_12=-mat_prod_12/det
      mat_prod_21=-mat_prod_21/det

c      write(*,*) 'inverse p',mat_prod_11,mat_prod_12,mat_prod_21
c     $     ,mat_prod_22
c      write(*,*) 'rrr',(mat_prod_obs22*mat_prod_21*B_prod_11 +B_prod_21
c     $     *(1.d0+mat_prod_22*mat_prod_obs22))*k*cdexp(-icomp*w(no)*z)
      mat_prod_21=mat_prod_21*B_prod_11+mat_prod_22*B_prod_21
      mat_prod_11=mat_prod_11*B_prod_11+mat_prod_12*B_prod_21
      mat_prod_22=0.d0
      mat_prod_12=0.d0


c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c         write(*,*) 'dessus',S11,S12,S21,S22,P11,P21
c         write(*,*) 'ctmp',ctmp,mat_sca_11,mat_sca_12,mat_prod_11
c         write(*,*)'exp',cdexp(-icomp*w(nc)*za)*cdexp(icomp*w(neps+1)*z)
c         write(*,*) 'matsca12',nu_plus*w(nc)*k2
      else
c     le point d'observation est dans le multicouche
         if (no.lt.nc) then
c     en dessous du dipole
c            write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs22*mat_sca_21*ctmp
            S22=mat_sca_obs22*mat_sca_22*ctmp
            P21=mat_prod_obs22*mat_prod_21*ctmp
         elseif (no.gt.nc) then
c     au dessus du dipole
c            write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs11*mat_sca_11*ctmp
            P11=mat_prod_obs11*mat_prod_11*ctmp
            S12=mat_sca_obs11*mat_sca_12*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs21*mat_sca_11*ctmp
            S22=mat_sca_obs21*mat_sca_12*ctmp
            P21=mat_prod_obs21*mat_prod_11*ctmp
         elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c            write(*,*) 'identique'
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)

            S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
            S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
            P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp
c            if (cdabs(S21)/cdabs(B_sca_21*ctmp).le.1.d-14) S21=0.d0
c            if (cdabs(S22)/cdabs(B_sca_22*ctmp).le.1.d-14) S22=0.d0
c            if (cdabs(P21)/cdabs(B_prod_21*ctmp).le.1.d-14) P21=0.d0

c            write(*,*) 'toto1',mat_sca_obs22*mat_sca_21,B_sca_21
c            write(*,*) 'toto1',mat_sca_obs22*mat_sca_22,B_sca_22
c            write(*,*) 'toto1',mat_prod_obs22*mat_prod_21,B_prod_21


c            write(*,*) 'ee',mat_prod_obs22,mat_prod_21,mat_prod_22
c     $           ,B_prod_11,B_prod_21,(mat_prod_obs22*mat_prod_21
c     $           +B_prod_21)*ctmp*k
c            write(*,*) 'rrr2',P21*k
c            S21=(mat_sca_obs22*mat_sca_21)*ctmp-w(nc)*cdexp(icomp*w(nc)
c     $           *(za-z)) 
c            S22=(mat_sca_obs22*mat_sca_22)*ctmp-k2*cdexp(icomp*w(nc)
c     $           *(za-z)) 
c            P21=(mat_prod_obs22*mat_prod_21)*ctmp-epscouche(nc)*k02
c     $           /w(nc)*cdexp(icomp*w(nc) *(za-z)) 


c            write(*,*) 'cc',k,S11,P11,S12,S21,S22,P21
c            write(*,*) 'tt',mat_prod_obs22 *mat_prod_21*ctmp,B_prod_21
c     $           *ctmp,dimag(k*P21),dimag(k*(S11+P11+S21 +P21))
c     avec l'espace libre
c            S21=mat_sca_obs22*mat_sca_21*ctmp
c            S22=mat_sca_obs22*mat_sca_22*ctmp
c            P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c            ctmp=cdexp(icomp*w(no)*z)
c            S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c            S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c            P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
c            ctmp=cdexp(-icomp*w(no)*z)
c            S21=(mat_sca_obs21*mat_sca_11)*ctmp
c            S22=(mat_sca_obs21*mat_sca_12)*ctmp
c            P21=(mat_prod_obs21*mat_prod_11)*ctmp
         endif
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22
      integrant(1)=k*(S11+P11+S21+P21)
      integrant(2)=(-S12+S22)*k/w(no)
 

c      write(30,*) k,dreal(integrant(1))
c      write(31,*) k,dimag(integrant(1))
c      write(32,*) k,dreal(integrant(2))
c      write(33,*) k,dimag(integrant(2))

      end
c******************************************************************
c******************************************************************
c******************************************************************
      subroutine Imulticouchecomp(hc,epsrel,epsabs,k0,a,z,za,neps
     $     ,dcouche,zcouche,epscouche,Ixx,Ixy,Ixz,Izx,Izz)
      implicit none
      integer nepsmax,neps,testloin,ndecoupe,ishanks,neval
      parameter (nepsmax=20)
      double precision dcouche(nepsmax),zcouche(0:nepsmax),epsmax,k0,z
     $     ,za,a,alpham,apara,hc
      double complex epscouche(0:nepsmax+1)
      double complex Ixx,Ixy,Ixz,Izx,Izz,Txx,Txy,Txz,Tzz

c     declaration lie a l'integration
      integer nnn
      double precision pi,aa,bb,kmax,hkmax     

c     declaration des variables du programme
      integer i,j,nbint,nlda,n,lordref
      parameter(nlda=5)
      double complex Ielliptique(nlda),Iinfini(nlda),Iinfinihp(nlda)
     $     ,Iinfinihm(nlda)
      double precision Imin,epsabs,epsrel,erreur(nlda)


      double complex icomp
      integer nnepsmax,nneps,nc,no
      double precision k02,zz,zza,zzcouche(0:nepsmax),zmax,zmin

      double complex eepscouche(0:nepsmax+1)

      
      common/donneemultin/nneps,nc,no
      common/donneemulti/k02,kmax,hkmax,apara,zz,zza
      common/zcoucheroutine/zzcouche
      common/epscouroutine/eepscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      external Imultiecomp,Imultiicomp,ImultiiH_pluscomp
     $     ,ImultiiH_moinscomp,Imultiilogcomp,ImultiiH_pluslogcomp
     $     ,ImultiiH_moinslogcomp
      external Imultiedessouscomp,Imultiidessouscomp
     $     ,ImultiiH_plusdessouscomp,ImultiiH_moinsdessouscomp
      external Imultiedessuscomp,Imultiidessuscomp
     $     ,ImultiiH_plusdessuscomp,ImultiiH_moinsdessuscomp
      
c     initialisation des donnees
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      nnn=5
      ishanks=0
      do i=1,nnn
         Ielliptique(i)=0.d0
         Iinfini(i)=0.d0
      enddo

c     repere dans quelle couche est le dipole
      nc=0
      do i=0,neps
         if (za.ge.zcouche(i)) then
            nc=nc+1
         else
            goto 10
         endif
      enddo
 10   no=0
c     write(*,*) 'numero couche dipole',nc,za
c     repere dans quelle couche est l'observation
      do i=0,neps
         if (z.ge.zcouche(i)) then
            no=no+1
         else
            goto 20
         endif
      enddo

c     write(*,*) 'numero couche observation',no,z,neps
c     calcul de zmax
 20   if (no.eq.nc) then
         call propesplibI(a,z,za,k0,epscouche(nc),Txx,Txy,Txz,Tzz)
         EPSABS=max(epsabs,min(cdabs(Txx),cdabs(Tzz),cdabs(Txy)
     $        ,cdabs(Txz))*EPSREL)
      endif

      if (no.eq.0.and.nc.eq.0) then
         zmin=dabs(z+za-2.d0*zcouche(0))
         goto 30
      elseif (no.eq.neps+1.and.nc.eq.neps+1) then
         zmin=dabs(z+za-2.d0*zcouche(neps))
         goto 30
      endif

      if (no.le.nc) then
         if (no.eq.0) then
            zmin=min(dabs(z-za),dabs(-z-za+2.d0*zcouche(nc)))
            zmax=zmin
c     write(*,*) 'coucou1',-z+za,-z-za+2.d0*zcouche(nc)
         else
            if (no.eq.nc) then
               zmin=min(dabs(-z-za+2.d0*zcouche(nc)),dabs(-z-za+2.d0
     $              *zcouche(nc+1)))
c               write(*,*) 'zmin',zmin
            else
               zmin=min(dabs(z-za),-z-za+2.d0*zcouche(nc))
            endif
            zmax=zmin
            zmin=min(zmin,dabs(z+za-2.d0*zcouche(no-1)))
            zmax=max(zmax,dabs(z+za-2.d0*zcouche(no-1)))
            zmin=min(zmin,dabs(z-za)-2.d0*zcouche(no-1)+2.d0*zcouche(nc)
     $           )
            zmax=max(zmax,dabs(z-za)-2.d0*zcouche(no-1)+2.d0*zcouche(nc)
     $           )
            
         endif
      elseif(no.gt.nc) then
         if (no.eq.neps+1) then
            zmin=min(dabs(z-za),z+za-2.d0*zcouche(nc-1))
            zmax=zmin
c     write(*,*) 'coucou1',-za+z,z+za-2.d0*zcouche(nc-1)
         else
            zmin=min(dabs(z-za),z+za-2.d0*zcouche(nc-1))
            zmax=zmin
c     write(*,*) 'coucou2',-za+z,z+za-2.d0*zcouche(nc)
            zmin=min(zmin,-z-za+2.d0*zcouche(no))
            zmax=max(zmax,-z-za+2.d0*zcouche(no))
c     write(*,*) 'coucou3',-z-za+2.d0*zcouche(no)
            zmin=min(zmin,-z+za+2.d0*zcouche(no)-2.d0*zcouche(nc-1))
            zmax=max(zmax,-z+za+2.d0*zcouche(no)-2.d0*zcouche(nc-1))
c     write(*,*) 'coucou4',-z+za+2.d0*zcouche(no)-2.d0*zcouche(nc
c     $           -1),zmin,zmax
         endif
      endif

c     write(*,*) 'zzzz',zmin,dabs(z-za),a

c     repette les donnees a cause du common de merde
 30   k02=k0*k0
      apara=a
      zz=z
      zza=za
      nnepsmax=nepsmax
      nneps=neps
      zzcouche(0)=zcouche(0)
      eepscouche(0)=epscouche(0)
      do i=1,neps
         zzcouche(i)=zcouche(i)
         eepscouche(i)=epscouche(i)
      enddo
      eepscouche(neps+1)=epscouche(neps+1)

c     commence l'integration
      epsmax=0.d0
      do i=0,neps+1
         epsmax=max(epsmax,dsqrt(dabs(dreal(epscouche(i)))))
      enddo
c     kmax pour le chemin elliptique +k0 pour plus de surete
      kmax=k0*epsmax+k0/2.d0
c      write(*,*) 'donnees',k0,kmax,epsmax,kmax*a
c     signe moins parce qu'on passe par la partie imaginaire negative
c     pour eviter les poles
      hkmax=-min(hc*kmax,1.d0/a)
c      write(*,*) 'ggg',a*kmax/6.28d0/2.d0,dint(a*kmax/6.28d0/2.d0)+1
c     hkmax=-hc*kmax

c     PREMIERE INTEGRALE DE 0 A KMAX AVEC CHEMIN ELLIPTIQUE
c     regarde si il va y avoir beaucoup d'osillation a grand ou pas
      if (a*kmax/6.28d0/2.d0.le.1.d0) then
         
         aa=-pi/2.d0
         bb=pi/2.d0
         
c         write(*,*) 'nc',nc,neps
         if (nc.eq.0) then
            call intgausskronrodpattersonmulti(aa,bb,ishanks,Ielliptique
     $           ,erreur,epsrel,EPSABS,lordref,Imultiedessouscomp,nlda
     $           ,nnn)
            NEVAL=2**lordref-1
         elseif (nc.eq.neps+1) then

            call intgausskronrodpattersonmulti(aa,bb,ishanks,Ielliptique
     $           ,erreur,epsrel,EPSABS,lordref,Imultiedessuscomp,nlda
     $           ,nnn)
            NEVAL=2**lordref-1
         else
            call intgausskronrodpattersonmulti(aa,bb,ishanks,Ielliptique
     $           ,erreur,epsrel,EPSABS,lordref,Imultiecomp,nlda,nnn)
            NEVAL=2**lordref-1
c
c            write(*,*) 'Ielliptique',Ielliptique,NEVAL
c            call intgausskronrodpattersonmultisub(aa,bb,ishanks
c     $           ,Ielliptique,erreur,epsrel,EPSABS,neval,Imultiecomp
c     $           ,nlda,nnn)
c            write(*,*) 'Ielliptique',Ielliptique,neval
c            stop


         endif
      else        
         ndecoupe=dint(a*kmax/6.28d0/2.d0)+1
         Ielliptique=0.d0
c         write(*,*) 'ndecoupe',ndecoupe,nc
         NEVAL=0
         if (nc.eq.0) then
            
            do i=1,ndecoupe              
               aa=-pi/2.d0+dble(i-1)/dble(ndecoupe)*pi
               bb=-pi/2.d0+dble(i)/dble(ndecoupe)*pi               
               call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $              ,erreur,epsrel,EPSABS,lordref,Imultiedessouscomp
     $              ,nlda ,nnn)
c               write(*,*) 'ttt',aa,bb,Iinfini
               do j=1,nnn
                  Ielliptique(j)=Ielliptique(j)+Iinfini(j)
               enddo
            enddo
            NEVAL=NEVAL+2**lordref-1
         elseif (nc.eq.neps+1) then
            do i=1,ndecoupe
               aa=-pi/2.d0+dble(i-1)/dble(ndecoupe)*pi
               bb=-pi/2.d0+dble(i)/dble(ndecoupe)*pi
               call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $              ,erreur,epsrel,EPSABS,lordref,Imultiedessuscomp,nlda
     $              ,nnn)
               NEVAL=NEVAL+2**lordref-1
c               write(*,*) 'ttt',aa,bb,nnn,Iinfini
               do j=1,nnn
                  Ielliptique(j)=Ielliptique(j)+Iinfini(j)
               enddo
            enddo
         else
            do i=1,ndecoupe
               aa=-pi/2.d0+dble(i-1)/dble(ndecoupe)*pi
               bb=-pi/2.d0+dble(i)/dble(ndecoupe)*pi
               call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $              ,erreur,epsrel,EPSABS,lordref,Imultiecomp,nlda
     $              ,nnn)
               NEVAL=NEVAL+2**lordref-1
c               write(*,*) 'nnn',2**lordref-1,i,ndecoupe
               do j=1,nnn
                  Ielliptique(j)=Ielliptique(j)+Iinfini(j)
               enddo
            enddo
         endif
         
      endif
c     do i=1,nnn
c        write(*,*) 'Ielliptique1',Ielliptique(i),NEVAL
c     enddo
c     stop
c      write(99,*) NEVAL
c      write(*,*) 'zmin',zmin
c     zmin plus grand que 2 lambda: plus d'onde evanescente
      if (zmin.ge.12.d0/k0) then
         Ixx=Ielliptique(1)*0.5d0*icomp
         Ixy=Ielliptique(2)*0.5d0*icomp
         Ixz=-Ielliptique(3)
         Izx=-Ielliptique(4)
         Izz=Ielliptique(5)*icomp
         NEVAL=0
         goto 40
      endif
c     DEUXIEME INTEGRALE DE KMAX A L'INFINI AVEC CHEMIN EN DROITE LIGNE
c     defini la borne superieur de l'integrale par raaport au maxi de
c     l'exponentiel.
      Imin=min(cdabs(Ielliptique(1)), cdabs(Ielliptique(2)),
     $     cdabs(Ielliptique(3)), cdabs(Ielliptique(4))
     $     ,cdabs(Ielliptique(5)))
      epsabs=max(epsabs,Imin*epsrel)

c     compare par rapport a zmin car c'est cela qui defini la
c     decroissance de l'esponentiel..
      if (zmin.ge.a) then
         alpham=-dlog(min(epsrel,epsabs/Imin))/zmin+kmax
c         write(*,*) 'alpham',alpham,alpham/k0
         if (zcouche(neps).ge.350.d0/alpham) then
c            write(*,*) 'depasse les capacites machine'
c            write(*,*) 'zcouche',zcouche(neps),350.d0/alpham
c            write(*,*) 'zmin alpham',zmin,alpham
            alpham=350.d0/zcouche(neps)
c            write(*,*) 'alpham possible',alpham
c     stop
         endif
         aa=kmax
         bb=alpham
         testloin=0

         if (nc.eq.0) then
            
            call Imultiedessouscomp(aa,nnn,nlda,Iinfini)         
            do i=1,nnn
               if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin
     $              =testloin+1
               Iinfini(i)=0.d0
            enddo
            NEVAL=1
            if (testloin.eq.nnn) goto 50

            call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $           ,erreur,epsrel,EPSABS,lordref,Imultiidessouscomp,nlda
     $           ,nnn)
            NEVAL=2**lordref-1
         elseif  (nc.eq.neps+1) then

            
            call Imultiidessuscomp(aa,nnn,nlda,Iinfini)         
            do i=1,nnn
               if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin
     $              =testloin+1
               Iinfini(i)=0.d0
            enddo
            NEVAL=1
            if (testloin.eq.nnn) goto 50

            call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $           ,erreur,epsrel,EPSABS,lordref,Imultiidessuscomp,nlda
     $           ,nnn)
            NEVAL=2**lordref-1
         else
            call Imultiicomp(aa,nnn,nlda,Iinfini)         
            do i=1,nnn
               if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin
     $              =testloin+1
               Iinfini(i)=0.d0
            enddo
            NEVAL=1
            if (testloin.eq.nnn) goto 50

            if (no.eq.nc) then

               call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $              ,erreur,epsrel,EPSABS,lordref,Imultiilogcomp,nlda
     $              ,nnn)
               NEVAL=2**lordref-1
            else

               call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfini
     $              ,erreur,epsrel,EPSABS,lordref,Imultiicomp,nlda,nnn)
               NEVAL=2**lordref-1
            endif
         endif
c         do i=1,nnn
c            write(*,*) 'Iinfini',Iinfini(i),i,NEVAL
c         enddo
 50      Ixx=(Ielliptique(1)+Iinfini(1))*0.5d0*icomp
         Ixy=(Ielliptique(2)+Iinfini(2))*0.5d0*icomp
         Ixz=-(Ielliptique(3)+Iinfini(3))
         Izx=-(Ielliptique(4)+Iinfini(4))
         Izz=(Ielliptique(5)+Iinfini(5))*icomp
      else
         alpham=-dlog(epsrel)/a
c         write(*,*) 'alphamH',alpham,alpham/k0
         if (zcouche(neps).ge.350.d0/alpham) then
c            write(*,*) 'depasse les capacites machine'
c            write(*,*) 'zcouche',zcouche(neps),350.d0/alpham
c            write(*,*) 'zmin alpham',zmin,alpham
            alpham=350.d0/zcouche(neps)
c            write(*,*) 'alpham possible',alpham
c     stop
         endif
         aa=0.d0
         bb=alpham
         if (nc.eq.0) then
            call ImultiiH_plusdessouscomp(aa,nnn,nlda,Iinfini)         
            do i=1,nnn
               if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin
     $              =testloin+1
               Iinfini(i)=0.d0
            enddo
            NEVAL=1
            if (testloin.eq.nnn) goto 60

            call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfinihp
     $           ,erreur,epsrel,EPSABS,lordref,ImultiiH_plusdessouscomp
     $           ,nlda ,nnn)
            NEVAL=2**lordref-1
c            write(*,*) 'oucou1'
         elseif (nc.eq.neps+1) then
            call ImultiiH_plusdessuscomp(aa,nnn,nlda,Iinfini)         
            do i=1,nnn
               if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin
     $              =testloin+1
               Iinfini(i)=0.d0
            enddo
            NEVAL=1
            if (testloin.eq.nnn) goto 60

            call intgausskronrodpattersonmulti(aa,bb,ishanks,Iinfinihp
     $           ,erreur,epsrel,EPSABS,lordref,ImultiiH_plusdessuscomp
     $           ,nlda ,nnn)
            NEVAL=2**lordref-1
c            write(*,*) 'oucou2'
         else
            call ImultiiH_pluscomp(aa,nnn,nlda,Iinfini)    
            do i=1,nnn
c               write(*,*) 'ttt',cdabs(Iinfini(i))*(bb-aa),epsabs,i,nnn
               if (cdabs(Iinfini(i))*(bb-aa).le.epsabs) testloin
     $              =testloin+1
               Iinfini(i)=0.d0
            enddo
            NEVAL=1
            if (testloin.eq.nnn) goto 60

            if (no.eq.nc) then
c               write(*,*) 'oucou3'
               call intgausskronrodpattersonmulti(aa,bb,ishanks
     $              ,Iinfinihp,erreur,epsrel,EPSABS,lordref
     $              ,ImultiiH_pluslogcomp,nlda,nnn)
               NEVAL=2**lordref-1
            else
c               write(*,*) 'oucou4'
               call intgausskronrodpattersonmulti(aa,bb,ishanks
     $              ,Iinfinihp,erreur,epsrel,EPSABS,lordref
     $              ,ImultiiH_pluscomp,nlda,nnn)
               NEVAL=2**lordref-1
            endif
         endif
c         do i=1,nnn
c            write(*,*) 'Iinfinihp',Iinfinihp(i),i,NEVAL
c         enddo


         aa=0.d0
         bb=alpham
         if (nc.eq.0) then
            call intgausskronrodpattersonmulti(aa,bb,ishanks ,Iinfinihm
     $           ,erreur,epsrel,EPSABS,lordref,ImultiiH_moinsdessouscomp
     $           ,nlda,nnn)
            NEVAL=2**lordref-1

         elseif (nc.eq.neps+1) then
            call intgausskronrodpattersonmulti(aa,bb,ishanks ,Iinfinihm
     $           ,erreur,epsrel,EPSABS,lordref ,ImultiiH_moinsdessuscomp
     $           ,nlda,nnn)
            NEVAL=2**lordref-1
         else
            if (no.eq.nc) then
               call intgausskronrodpattersonmulti(aa,bb,ishanks
     $              ,Iinfinihm,erreur,epsrel,EPSABS,lordref
     $              ,ImultiiH_moinslogcomp,nlda,nnn)
               NEVAL=2**lordref-1
            else
               call intgausskronrodpattersonmulti(aa,bb,ishanks
     $              ,Iinfinihm,erreur,epsrel,EPSABS,lordref
     $              ,ImultiiH_moinscomp,nlda,nnn)
               NEVAL=2**lordref-1
            endif
         endif
c         do i=1,nnn
c            write(*,*) 'Iinfinihm',Iinfinihm(i),i,NEVAL
c         enddo
         do i=1,nnn
            Iinfini(i)=(Iinfinihp(i)-Iinfinihm(i))/2.d0
c            write(*,*) 'Iinfinit',Iinfini(i),i,NEVAL
         enddo
         
 60      Ixx=(Ielliptique(1)+icomp*Iinfini(1))*0.5d0*icomp
         Ixy=(Ielliptique(2)+icomp*Iinfini(2))*0.5d0*icomp
         Ixz=-(Ielliptique(3)+icomp*Iinfini(3))
         Izx=-(Ielliptique(4)+icomp*Iinfini(4))
         Izz=(Ielliptique(5)+icomp*Iinfini(5))*icomp
      endif

c 40   write(*,*) 'multi',Ixx,Ixy,Ixz,Izx,Izz
 40   if (no.eq.nc) then
c     write(*,*) 'libre',Txx,Txy,Txz,Tzz
         Ixx=Ixx+Txx
         Ixy=Ixy+Txy
         Ixz=Ixz+Txz
         Izx=Izx+Txz
         Izz=Izz+Tzz
      endif   
      write(98,*) NEVAL
c     divise par epsilon de la couche ou il y a le dipole
      Ixx=Ixx/epscouche(nc)
      Izz=Izz/epscouche(nc)
      Ixz=Ixz/epscouche(nc)
      Izx=Izx/epscouche(nc)
      Ixy=Ixy/epscouche(nc)


      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Imultiecomp(theta,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision theta
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex k,k2,wronskien

      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex JB0,JB1,JB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0

c      write(*,*) 'routine',k02,kmax,hkmax,a,z,za,neps,nc,no     

c     calcul de k complex avec theta
      zr=kmax/2.d0*(dsin(theta)+1.d0)
      zi=hkmax*dcos(theta)
      k=zr+icomp*zi
c     calcul des fonctions de Bessel a l'ordre 0,1,2
      zr=zr*a
      zi=zi*a
      call ZBESJ(ZR, ZI, FNU, 1, N, CYR, CYI, NZ, IERR)     
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de bessel',NZ,IERR
         stop
      endif
      JB0=cyr(1)+icomp*cyi(1)
      JB1=cyr(2)+icomp*cyi(2)
      JB2=cyr(3)+icomp*cyi(3)
      k2=k*k
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)
      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp
c         write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c         write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no.and.no.ne.nc) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      det=mat_sca_11*mat_sca_22-mat_sca_12*mat_sca_21
c      write(*,*) 'det',det,mat_sca_11,mat_sca_22,mat_sca_12,mat_sca_21
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22/det
      mat_sca_22=ctmp/det
      mat_sca_12=-mat_sca_12/det
      mat_sca_21=-mat_sca_21/det

      ctmp=mat_sca_11*B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=mat_sca_11*B_sca_12+mat_sca_12*B_sca_22
      mat_sca_11=ctmp      
      ctmp=mat_sca_21*B_sca_11+mat_sca_22*B_sca_21
      mat_sca_22=mat_sca_21*B_sca_12+mat_sca_22*B_sca_22
      mat_sca_21=ctmp

      det=mat_prod_11*mat_prod_22-mat_prod_12*mat_prod_21
c      write(*,*) 'det',det
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22/det
      mat_prod_22=ctmp/det
      mat_prod_12=-mat_prod_12/det
      mat_prod_21=-mat_prod_21/det

      mat_prod_21=mat_prod_21*B_prod_11+mat_prod_22*B_prod_21
      mat_prod_11=mat_prod_11*B_prod_11+mat_prod_12*B_prod_21
      mat_prod_22=0.d0
      mat_prod_12=0.d0

c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c         write(*,*) 'dessus',S11,S12,S21,S22,P11,P21,ctmp,mat_sca_11
c     $        ,mat_sca_12
      else
c     le point d'observation est dans le multicouche
         if (no.lt.nc) then
c     en dessous du dipole
c            write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22

            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs22*mat_sca_21*ctmp
            S22=mat_sca_obs22*mat_sca_22*ctmp
            P21=mat_prod_obs22*mat_prod_21*ctmp
         elseif (no.gt.nc) then
c     au dessus du dipole
c            write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs11*mat_sca_11*ctmp
            P11=mat_prod_obs11*mat_prod_11*ctmp
            S12=mat_sca_obs11*mat_sca_12*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs21*mat_sca_11*ctmp
            S22=mat_sca_obs21*mat_sca_12*ctmp
            P21=mat_prod_obs21*mat_prod_11*ctmp
         elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c            write(*,*) 'identique'
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
            S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
            P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c     avec l'espace libre
c            S21=mat_sca_obs22*mat_sca_21*ctmp
c            S22=mat_sca_obs22*mat_sca_22*ctmp
c            P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c            ctmp=cdexp(icomp*w(no)*z)
c            S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c            S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c            P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
c            ctmp=cdexp(-icomp*w(no)*z)
c            S21=(mat_sca_obs21*mat_sca_11)*ctmp
c            S22=(mat_sca_obs21*mat_sca_12)*ctmp
c            P21=(mat_prod_obs21*mat_prod_11)*ctmp
         endif
      endif


c     calcul element differentiel
      wronskien=kmax/2.d0*dcos(theta)-icomp*hkmax*dsin(theta)
c      write(*,*) 'wr',wronskien
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'SP',S11,P11,S21,P21,S12,S22


      integrant(1)=k*(S11+P11+S21+P21)*wronskien*JB0
      integrant(2)=k*(S11-P11+S21-P21)*wronskien*JB2
      integrant(3)=(S12+S22)*wronskien*JB1
      integrant(4)=k2/w(no)*(-S11+S21)*wronskien*JB1
      integrant(5)=(-S12+S22)*k/w(no)*wronskien*JB0

c      write(10,*) (theta),dreal(integrant(1))
c      write(11,*) (theta),dimag(integrant(1))
c      write(12,*) (theta),dreal(integrant(2))
c      write(13,*) (theta),dimag(integrant(2))

      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Imultiicomp(k,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision k,k2
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

c     declarations liees a netlib
      integer n,nz
      parameter (n=3)
      double precision y,alp,xx
      dimension y(n)


      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      k2=k*k
      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c      write(*,*) '**************',k

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp

c         write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c         write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c         write(*,*) 'prod',mat_prod_11,mat_prod_12 ,mat_prod_21
c     $        ,mat_prod_22
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no.and.no.ne.nc) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c      write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      det=mat_sca_11*mat_sca_22-mat_sca_12*mat_sca_21
c      write(*,*) 'det',det,mat_sca_11,mat_sca_22,mat_sca_12,mat_sca_21
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22/det
      mat_sca_22=ctmp/det
      mat_sca_12=-mat_sca_12/det
      mat_sca_21=-mat_sca_21/det
      
c      write(*,*) 'inverse s',mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
c      write(*,*) 'chiant1',mat_sca_12

      ctmp=mat_sca_11*B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=mat_sca_11*B_sca_12+mat_sca_12*B_sca_22
      mat_sca_11=ctmp      
      ctmp=mat_sca_21*B_sca_11+mat_sca_22*B_sca_21
      mat_sca_22=mat_sca_21*B_sca_12+mat_sca_22*B_sca_22
      mat_sca_21=ctmp

c      write(*,*) 'chiant2',mat_sca_12,mat_sca_11,B_sca_12,mat_sca_12
c     $     ,B_sca_22

      det=mat_prod_11*mat_prod_22-mat_prod_12*mat_prod_21
c      write(*,*) 'det',det
c      write(99,*) k,dreal(det),dimag(det)
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22/det
      mat_prod_22=ctmp/det
      mat_prod_12=-mat_prod_12/det
      mat_prod_21=-mat_prod_21/det

c      write(*,*) 'inverse p',mat_prod_11,mat_prod_12,mat_prod_21
c     $     ,mat_prod_22

      mat_prod_21=mat_prod_21*B_prod_11+mat_prod_22*B_prod_21
      mat_prod_11=mat_prod_11*B_prod_11+mat_prod_12*B_prod_21
      mat_prod_22=0.d0
      mat_prod_12=0.d0

c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c         write(*,*) 'dessus',S11,S12,S21,S22,P11,P21
c         write(*,*) 'ctmp',ctmp,mat_sca_11,mat_sca_12,mat_prod_11
c         write(*,*)'exp',cdexp(-icomp*w(nc)*za)*cdexp(icomp*w(neps+1)*z)
c         write(*,*) 'matsca12',nu_plus*w(nc)*k2
      else
c     le point d'observation est dans le multicouche
         if (no.lt.nc) then
c     en dessous du dipole
c            write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs22*mat_sca_21*ctmp
            S22=mat_sca_obs22*mat_sca_22*ctmp
            P21=mat_prod_obs22*mat_prod_21*ctmp
         elseif (no.gt.nc) then
c     au dessus du dipole
c            write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs11*mat_sca_11*ctmp
            P11=mat_prod_obs11*mat_prod_11*ctmp
            S12=mat_sca_obs11*mat_sca_12*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs21*mat_sca_11*ctmp
            S22=mat_sca_obs21*mat_sca_12*ctmp
            P21=mat_prod_obs21*mat_prod_11*ctmp
         elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c            write(*,*) 'identique'
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
c     les oscillations viennent de la. Je somme deux choses qui sont
c     quasi identiques au signe pres. la somme est proche de zero et
c     entache d'erreur numerique:mat_sca_obs22*mat_sca_21+B_sca_21
            S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
            S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
            P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c            write(99,9999) k,dreal(mat_sca_obs22*mat_sca_21+B_sca_21)

c     avec l'espace libre
c            S21=mat_sca_obs22*mat_sca_21*ctmp
c            S22=mat_sca_obs22*mat_sca_22*ctmp
c            P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c            ctmp=cdexp(icomp*w(no)*z)
c            S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c            S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c            P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
c            ctmp=cdexp(-icomp*w(no)*z)
c            S21=(mat_sca_obs21*mat_sca_11)*ctmp
c            S22=(mat_sca_obs21*mat_sca_12)*ctmp
c            P21=(mat_prod_obs21*mat_prod_11)*ctmp
         endif
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22

      alp=0.d0
      xx=a*k
      call DBESJ (XX,ALP,N,Y,NZ)
      if (NZ.ne.0) then
         write(*,*) 'PB calcul de BESSEL'
         stop
      endif

c      write(99,9999) k,dreal(Y(1)),dreal(S11),dreal(P11),dreal(S21)
c     $     ,dreal(P21)
      integrant(1)=k*(S11+P11+S21+P21)*Y(1)
      integrant(2)=k*(S11-P11+S21-P21)*Y(3)
      integrant(3)=(S12+S22)*Y(2)
      integrant(4)=k2/w(no)*(-S11+S21)*Y(2)
      integrant(5)=(-S12+S22)*k/w(no)*Y(1)

c      write(20,*) k,integrant(1)
c      write(21,*) k,integrant(2)
c      write(22,*) k,integrant(3)
c      write(23,*) k,integrant(4)
 9999 format(201(d22.15,1x))
      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine ImultiiH_pluscomp(kint,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision kint
      double complex S11,P11,S21,P21,S12,S22,icomp,k,k2,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex HB0,HB1,HB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0
      k=kmax+icomp*kint

c     premiere espece pour la fonction de Hankel
      zr=kmax*a
      zi=kint*a
      call ZBESH(ZR, ZI, FNU, 1, 1, N, CYR, CYI, NZ, IERR)
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de Hankel1',NZ,IERR
         write(*,*) 'argument',zr,zi
         stop
      endif
      HB0=cyr(1)+icomp*cyi(1)
      HB1=cyr(2)+icomp*cyi(2)
      HB2=cyr(3)+icomp*cyi(3)

      k2=k*k      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp

c         write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c         write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c         write(*,*) 'prod',mat_prod_11,mat_prod_12 ,mat_prod_21
c     $        ,mat_prod_22
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no.and.no.ne.nc) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c      write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      det=mat_sca_11*mat_sca_22-mat_sca_12*mat_sca_21
c      write(*,*) 'det',det,mat_sca_11,mat_sca_22,mat_sca_12,mat_sca_21
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22/det
      mat_sca_22=ctmp/det
      mat_sca_12=-mat_sca_12/det
      mat_sca_21=-mat_sca_21/det
      
c      write(*,*) 'inverse s',mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
c      write(*,*) 'chiant1',mat_sca_12

      ctmp=mat_sca_11*B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=mat_sca_11*B_sca_12+mat_sca_12*B_sca_22
      mat_sca_11=ctmp      
      ctmp=mat_sca_21*B_sca_11+mat_sca_22*B_sca_21
      mat_sca_22=mat_sca_21*B_sca_12+mat_sca_22*B_sca_22
      mat_sca_21=ctmp

c      write(*,*) 'chiant2',mat_sca_12,mat_sca_11,B_sca_12,mat_sca_12
c     $     ,B_sca_22

      det=mat_prod_11*mat_prod_22-mat_prod_12*mat_prod_21
c      write(*,*) 'det',det
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22/det
      mat_prod_22=ctmp/det
      mat_prod_12=-mat_prod_12/det
      mat_prod_21=-mat_prod_21/det

c      write(*,*) 'inverse p',mat_prod_11,mat_prod_12,mat_prod_21
c     $     ,mat_prod_22

      mat_prod_21=mat_prod_21*B_prod_11+mat_prod_22*B_prod_21
      mat_prod_11=mat_prod_11*B_prod_11+mat_prod_12*B_prod_21
      mat_prod_22=0.d0
      mat_prod_12=0.d0

c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c         write(*,*) 'dessus',S11,S12,S21,S22,P11,P21
c         write(*,*) 'ctmp',ctmp,mat_sca_11,mat_sca_12,mat_prod_11
c         write(*,*)'exp',cdexp(-icomp*w(nc)*za)*cdexp(icomp*w(neps+1)*z)
c         write(*,*) 'matsca12',nu_plus*w(nc)*k2
      else
c     le point d'observation est dans le multicouche
         if (no.lt.nc) then
c     en dessous du dipole
c            write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs22*mat_sca_21*ctmp
            S22=mat_sca_obs22*mat_sca_22*ctmp
            P21=mat_prod_obs22*mat_prod_21*ctmp
         elseif (no.gt.nc) then
c     au dessus du dipole
c            write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs11*mat_sca_11*ctmp
            P11=mat_prod_obs11*mat_prod_11*ctmp
            S12=mat_sca_obs11*mat_sca_12*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs21*mat_sca_11*ctmp
            S22=mat_sca_obs21*mat_sca_12*ctmp
            P21=mat_prod_obs21*mat_prod_11*ctmp
         elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c            write(*,*) 'identique'
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)

            S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
            S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
            P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c     avec l'espace libre
c            S21=mat_sca_obs22*mat_sca_21*ctmp
c            S22=mat_sca_obs22*mat_sca_22*ctmp
c            P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c            ctmp=cdexp(icomp*w(no)*z)
c            S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c            S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c            P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
c            ctmp=cdexp(-icomp*w(no)*z)
c            S21=(mat_sca_obs21*mat_sca_11)*ctmp
c            S22=(mat_sca_obs21*mat_sca_12)*ctmp
c            P21=(mat_prod_obs21*mat_prod_11)*ctmp
         endif
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22
 
      integrant(1)=k*(S11+P11+S21+P21)*HB0
      integrant(2)=k*(S11-P11+S21-P21)*HB2
      integrant(3)=(S12+S22)*HB1
      integrant(4)=k2/w(no)*(-S11+S21)*HB1
      integrant(5)=(-S12+S22)*k/w(no)*HB0


c      write(30,*) kint,integrant(1)
c      write(31,*) kint,integrant(2)
c      write(32,*) kint,integrant(3)
c      write(33,*) kint,integrant(4)
      
      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine ImultiiH_moinscomp(kint,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision kint
      double complex S11,P11,S21,P21,S12,S22,icomp,k,k2,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex HB0,HB1,HB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0
      k=kmax-icomp*kint

c     deuxieme espece pour la fonction de Hankel
      zr=kmax*a
      zi=-kint*a
      call ZBESH(ZR, ZI, FNU, 1, 2, N, CYR, CYI, NZ, IERR)
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de Hankel2',NZ,IERR
         write(*,*) 'argument',zr,zi
         stop
      endif
      HB0=cyr(1)+icomp*cyi(1)
      HB1=cyr(2)+icomp*cyi(2)
      HB2=cyr(3)+icomp*cyi(3)

      k2=k*k      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp

c         write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c         write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c         write(*,*) 'prod',mat_prod_11,mat_prod_12 ,mat_prod_21
c     $        ,mat_prod_22
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp


         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no.and.no.ne.nc) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c      write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      det=mat_sca_11*mat_sca_22-mat_sca_12*mat_sca_21
c      write(*,*) 'det',det,mat_sca_11,mat_sca_22,mat_sca_12,mat_sca_21
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22/det
      mat_sca_22=ctmp/det
      mat_sca_12=-mat_sca_12/det
      mat_sca_21=-mat_sca_21/det
      
c      write(*,*) 'inverse s',mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
c      write(*,*) 'chiant1',mat_sca_12

      ctmp=mat_sca_11*B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=mat_sca_11*B_sca_12+mat_sca_12*B_sca_22
      mat_sca_11=ctmp      
      ctmp=mat_sca_21*B_sca_11+mat_sca_22*B_sca_21
      mat_sca_22=mat_sca_21*B_sca_12+mat_sca_22*B_sca_22
      mat_sca_21=ctmp

c      write(*,*) 'chiant2',mat_sca_12,mat_sca_11,B_sca_12,mat_sca_12
c     $     ,B_sca_22

      det=mat_prod_11*mat_prod_22-mat_prod_12*mat_prod_21
c      write(*,*) 'det',det
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22/det
      mat_prod_22=ctmp/det
      mat_prod_12=-mat_prod_12/det
      mat_prod_21=-mat_prod_21/det

c      write(*,*) 'inverse p',mat_prod_11,mat_prod_12,mat_prod_21
c     $     ,mat_prod_22

      mat_prod_21=mat_prod_21*B_prod_11+mat_prod_22*B_prod_21
      mat_prod_11=mat_prod_11*B_prod_11+mat_prod_12*B_prod_21
      mat_prod_22=0.d0
      mat_prod_12=0.d0

c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c         write(*,*) 'dessus',S11,S12,S21,S22,P11,P21
c         write(*,*) 'ctmp',ctmp,mat_sca_11,mat_sca_12,mat_prod_11
c         write(*,*)'exp',cdexp(-icomp*w(nc)*za)*cdexp(icomp*w(neps+1)*z)
c         write(*,*) 'matsca12',nu_plus*w(nc)*k2
      else
c     le point d'observation est dans le multicouche
         if (no.lt.nc) then
c     en dessous du dipole
c            write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs22*mat_sca_21*ctmp
            S22=mat_sca_obs22*mat_sca_22*ctmp
            P21=mat_prod_obs22*mat_prod_21*ctmp
         elseif (no.gt.nc) then
c     au dessus du dipole
c            write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs11*mat_sca_11*ctmp
            P11=mat_prod_obs11*mat_prod_11*ctmp
            S12=mat_sca_obs11*mat_sca_12*ctmp
            ctmp=cdexp(-icomp*w(no)*z)
            S21=mat_sca_obs21*mat_sca_11*ctmp
            S22=mat_sca_obs21*mat_sca_12*ctmp
            P21=mat_prod_obs21*mat_prod_11*ctmp
         elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c            write(*,*) 'identique'
            ctmp=cdexp(icomp*w(no)*z)
            S11=mat_sca_obs12*mat_sca_21*ctmp
            P11=mat_prod_obs12*mat_prod_21*ctmp
            S12=mat_sca_obs12*mat_sca_22*ctmp
            ctmp=cdexp(-icomp*w(no)*z)

            S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
            S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
            P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c     avec l'espace libre
c            S21=mat_sca_obs22*mat_sca_21*ctmp
c            S22=mat_sca_obs22*mat_sca_22*ctmp
c            P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c            ctmp=cdexp(icomp*w(no)*z)
c            S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c            S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c            P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
c            ctmp=cdexp(-icomp*w(no)*z)
c            S21=(mat_sca_obs21*mat_sca_11)*ctmp
c            S22=(mat_sca_obs21*mat_sca_12)*ctmp
c            P21=(mat_prod_obs21*mat_prod_11)*ctmp
         endif
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22
 
      integrant(1)=k*(S11+P11+S21+P21)*HB0
      integrant(2)=k*(S11-P11+S21-P21)*HB2
      integrant(3)=(S12+S22)*HB1
      integrant(4)=k2/w(no)*(-S11+S21)*HB1
      integrant(5)=(-S12+S22)*k/w(no)*HB0


c      write(40,*) kint,integrant(1)
c      write(41,*) kint,integrant(2)
c      write(42,*) kint,integrant(3)
c      write(43,*) kint,integrant(4)

      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Idiagmultiedessouscomp(theta,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision theta
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex k,k2,wronskien

      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)

c      write(*,*) 'routine',k02,kmax,hkmax,a,z,za,neps,nc,no     

c     calcul de k complex avec theta
      k=kmax/2.d0*(dsin(theta)+1.d0)+icomp*hkmax*dcos(theta)
      k2=k*k
c      write(*,*) k,k2
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)
      enddo

c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca2',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod2',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22
c         write(*,*) 'i2',i,no,nc,neps
         if (i.eq.no) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo


c      write(*,*) 'obs2',mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
c     $     ,mat_sca_obs22,mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
c     $     ,mat_prod_obs22

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=0.d0
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-(1.d0,0.d0)
c      write(*,*) 'rr2',mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      mat_prod_11=mat_prod_sur11
      mat_prod_12=0.d0
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-(1.d0,0.d0)


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_21=mat_sca_21/mat_sca_11
      mat_sca_12=B_sca_12/mat_sca_11
      mat_sca_11=B_sca_11/mat_sca_11
      mat_sca_22=mat_sca_21*B_sca_12-B_sca_22
      mat_sca_21=mat_sca_21*B_sca_11-B_sca_21
      mat_prod_21=mat_prod_21/mat_prod_11
      mat_prod_21=mat_prod_21*B_prod_11-B_prod_21
      mat_prod_11=B_prod_11/mat_prod_11

c     remonte au point d'observation
      if (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c         write(*,*) 'dessus',S11,S12,S21,S22,P11,P21,ctmp,mat_sca_11
c     $        ,mat_sca_12
     
      elseif (no.gt.nc) then
c     au dessus du dipole
c            write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs11*mat_sca_11*ctmp
         P11=mat_prod_obs11*mat_prod_11*ctmp
         S12=mat_sca_obs11*mat_sca_12*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs21*mat_sca_11*ctmp
         S22=mat_sca_obs21*mat_sca_12*ctmp
         P21=mat_prod_obs21*mat_prod_11*ctmp
c     write(*,*) 'S2222',S22
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c            write(*,*) 'identique'
c         ctmp=cdexp(icomp*w(no)*z)
c         S11=mat_sca_obs12*mat_sca_21*ctmp
c         P11=mat_prod_obs12*mat_prod_21*ctmp
c         S12=mat_sca_obs12*mat_sca_22*ctmp
c         ctmp=cdexp(-icomp*w(no)*z)
c         S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c         S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c         write(*,*) 'S22',S22,mat_sca_obs22*mat_sca_22,B_sca_22,ctmp
c         P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp
c         ctmp=cdexp(-icomp*w(0)*z)

c     avec l'espace libre
c            S21=mat_sca_obs22*mat_sca_21*ctmp
c            S22=mat_sca_obs22*mat_sca_22*ctmp
c            P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c            ctmp=cdexp(icomp*w(no)*z)
c            S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c            S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c            P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
            S11=0.d0
            S12=0.d0
            P11=0.d0
            ctmp=cdexp(-icomp*w(no)*z)
            S21=(mat_sca_obs21*mat_sca_11)*ctmp
            S22=(mat_sca_obs21*mat_sca_12)*ctmp
            P21=(mat_prod_obs21*mat_prod_11)*ctmp
         

      endif


c     calcul element differentiel
      wronskien=kmax/2.d0*dcos(theta)-icomp*hkmax*dsin(theta)
c      write(*,*) 'wr',wronskien
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'SSS',S11,P11,S21,P21,S12,S22
      integrant(1)=k*(S11+P11+S21+P21)*wronskien
      integrant(2)=(-S12+S22)*wronskien*k/w(no)

c      write(*,*) 'theta',theta
c      write(*,*)  integrant(1),integrant(2)
c      write(*,*) '*********************************'
c      write(10,*) (theta),dreal(integrant(1))
c      write(11,*) (theta),dimag(integrant(1))
c      write(12,*) (theta),dreal(integrant(2))
c      write(13,*) (theta),dimag(integrant(2))

      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Idiagmultiidessouscomp(k,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision k,k2
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      k2=k*k
      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
      enddo

c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp


         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22
c     evite le cas no=nc fait au dessus
         if (i.eq.no) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=0.d0
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-(1.d0,0.d0)

      mat_prod_11=mat_prod_sur11
      mat_prod_12=0.d0
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-(1.d0,0.d0)


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c      write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_21=mat_sca_21/mat_sca_11
      mat_sca_12=B_sca_12/mat_sca_11
      mat_sca_11=B_sca_11/mat_sca_11
      mat_sca_22=mat_sca_21*B_sca_12-B_sca_22
      mat_sca_21=mat_sca_21*B_sca_11-B_sca_21
      mat_prod_21=mat_prod_21/mat_prod_11
      mat_prod_21=mat_prod_21*B_prod_11-B_prod_21
      mat_prod_11=B_prod_11/mat_prod_11

      if (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c         write(*,*) 'dessus',S11,S12,S21,S22,P11,P21,ctmp,mat_sca_11
c     $        ,mat_sca_12
         
      elseif (no.gt.nc) then
c     au dessus du dipole
c         write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $        ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs11*mat_sca_11*ctmp
         P11=mat_prod_obs11*mat_prod_11*ctmp
         S12=mat_sca_obs11*mat_sca_12*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs21*mat_sca_11*ctmp
         S22=mat_sca_obs21*mat_sca_12*ctmp
         P21=mat_prod_obs21*mat_prod_11*ctmp
c     write(*,*) 'S2222',S22
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c         write(*,*) 'identique'
c     ctmp=cdexp(icomp*w(no)*z)
c     S11=mat_sca_obs12*mat_sca_21*ctmp
c     P11=mat_prod_obs12*mat_prod_21*ctmp
c     S12=mat_sca_obs12*mat_sca_22*ctmp
c     ctmp=cdexp(-icomp*w(no)*z)
c     S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c     S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c     write(*,*) 'S22',S22,mat_sca_obs22*mat_sca_22,B_sca_22,ctmp
c     P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp
c     ctmp=cdexp(-icomp*w(0)*z)

c     avec l'espace libre
c     S21=mat_sca_obs22*mat_sca_21*ctmp
c     S22=mat_sca_obs22*mat_sca_22*ctmp
c     P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c         ctmp=cdexp(icomp*w(no)*z)
c         S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c         S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c         P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
         S11=0.d0
         S12=0.d0
         P11=0.d0
         ctmp=cdexp(-icomp*w(no)*z)
         S21=(mat_sca_obs21*mat_sca_11)*ctmp
         S22=(mat_sca_obs21*mat_sca_12)*ctmp
         P21=(mat_prod_obs21*mat_prod_11)*ctmp
c         write(*,*) 'S22',S22,mat_sca_obs22*mat_sca_22,B_sca_22,ctmp

      endif

c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) 'fin',S11,P11,S21,P21,S12,S22
      integrant(1)=k*(S11+P11+S21+P21)
      integrant(2)=(-S12+S22)*k/w(no)

c      write(60,*) k,integrant(1)
c      write(61,*) k,integrant(2)
c      write(62,*) k,integrant(3)
c      write(63,*) k,integrant(4)

      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Idiagmultiedessuscomp(theta,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision theta
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex k,k2,wronskien

      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)

c      write(*,*) 'routine ancienne',k02,kmax,hkmax,a,z,za,neps,nc,no     

c     calcul de k complex avec theta
      k=kmax/2.d0*(dsin(theta)+1.d0)+icomp*hkmax*dcos(theta)
      k2=k*k
c      write(*,*) 'kkk****************',k,theta
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
c     write(*,*) 'eplus',e_plus,icomp*(w(i)+w(i+1))*zcouche(i),w(i)
c     $        ,w(i+1),zcouche(i)
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp
c     write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c     write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire         
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c     write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c     write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c     stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c      write(*,*) 'obs1',mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
c     $     ,mat_sca_obs22,mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
c     $     ,mat_prod_obs22
c     calcul de la matrice coupee
      mat_sca_11=(1.d0,0.d0)
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=0.d0
      mat_sca_22=-mat_sca_sous22
c      write(*,*) 'rr1',mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      mat_prod_11=(1.d0,0.d0)
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=0.d0
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)
c      write(*,*) 't',mat_sca_obs22/mat_sca_22
c     inverse de la matrice coupee et produit par B
      mat_sca_12=-mat_sca_12/mat_sca_22
      mat_sca_11=B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=B_sca_12+mat_sca_12*B_sca_22
      mat_sca_21=B_sca_21/mat_sca_22
      mat_sca_22=B_sca_22/mat_sca_22

      mat_prod_12=-mat_prod_12/mat_prod_22
      mat_prod_22=1.d0/mat_prod_22
      mat_prod_21=mat_prod_22*B_prod_21
      mat_prod_11=B_prod_11+mat_prod_12*B_prod_21
    

c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp

      elseif (no.lt.nc) then
c     en dessous du dipole
c     write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22

         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs22*mat_sca_21*ctmp
         S22=mat_sca_obs22*mat_sca_22*ctmp
         P21=mat_prod_obs22*mat_prod_21*ctmp
         
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c         write(*,*) 'identique'
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
c         ctmp=cdexp(-icomp*w(no)*z)
c         S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c         S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c         write(*,*) 'S21',mat_sca_obs22*mat_sca_21,B_sca_21
c         write(*,*) 'S22',mat_sca_obs22*mat_sca_22,B_sca_22
c         write(*,*) 'P21',mat_prod_obs22*mat_prod_21,B_prod_21
c         write(*,*) mat_sca_obs22,mat_sca_21,mat_sca_22,mat_prod_21
c         P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp
         S21=0.d0
         S22=0.d0
         P21=0.d0
c     avec l'espace libre
c     S21=mat_sca_obs22*mat_sca_21*ctmp
c     S22=mat_sca_obs22*mat_sca_22*ctmp
c     P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c     ctmp=cdexp(icomp*w(no)*z)
c     S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c     S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c     P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
c     ctmp=cdexp(-icomp*w(no)*z)
c     S21=(mat_sca_obs21*mat_sca_11)*ctmp
c     S22=(mat_sca_obs21*mat_sca_12)*ctmp
c     P21=(mat_prod_obs21*mat_prod_11)*ctmp
         
      endif


c     calcul element differentiel
      wronskien=kmax/2.d0*dcos(theta)-icomp*hkmax*dsin(theta)
c      write(*,*) 'wr',wronskien
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'SSS',S11,P11,S12,S22,S21,P21
      integrant(1)=k*(S11+P11+S21+P21)*wronskien
      integrant(2)=(-S12+S22)*wronskien*k/w(no)

c     write(*,*) 'theta',theta
c     write(*,*)  integrant(1),integrant(2),integrant(3),integrant(4)
c     write(*,*) '*********************************'
c      write(110,*) (theta),dreal(integrant(1))
c      write(111,*) (theta),dimag(integrant(1))
c      write(112,*) (theta),dreal(integrant(2))
c      write(113,*) (theta),dimag(integrant(2))

      end

c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Idiagmultiidessuscomp(k,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision k,k2
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      k2=k*k
      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp

c     write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c     write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     write(*,*) 'prod',mat_prod_11,mat_prod_12 ,mat_prod_21
c     $        ,mat_prod_22
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c     write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c     write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c     stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo



c     calcul de la matrice coupee
      mat_sca_11=(1.d0,0.d0)
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=0.d0
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=(1.d0,0.d0)
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=0.d0
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c      write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_12=-mat_sca_12/mat_sca_22
      mat_sca_11=B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=B_sca_12+mat_sca_12*B_sca_22
      mat_sca_21=B_sca_21/mat_sca_22
      mat_sca_22=B_sca_22/mat_sca_22

      mat_prod_12=-mat_prod_12/mat_prod_22
      mat_prod_22=1.d0/mat_prod_22
      mat_prod_21=mat_prod_22*B_prod_21
      mat_prod_11=B_prod_11+mat_prod_12*B_prod_21
    
c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.lt.nc) then
c     en dessous du dipole
c     write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs22*mat_sca_21*ctmp
         S22=mat_sca_obs22*mat_sca_22*ctmp
         P21=mat_prod_obs22*mat_prod_21*ctmp
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c         write(*,*) 'identique'
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
c         ctmp=cdexp(-icomp*w(no)*z)
c         S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c         S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c         P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp
         S21=0.d0
         S22=0.d0
         P21=0.d0
c     avec l'espace libre
c     S21=mat_sca_obs22*mat_sca_21*ctmp
c     S22=mat_sca_obs22*mat_sca_22*ctmp
c     P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c     ctmp=cdexp(icomp*w(no)*z)
c     S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c     S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c     P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
c     ctmp=cdexp(-icomp*w(no)*z)
c     S21=(mat_sca_obs21*mat_sca_11)*ctmp
c     S22=(mat_sca_obs21*mat_sca_12)*ctmp
c     P21=(mat_prod_obs21*mat_prod_11)*ctmp
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22,k*(S11+P11+S21+P21),(-S12
c     $     +S22)*k/w(no)
      integrant(1)=k*(S11+P11+S21+P21)
      integrant(2)=(-S12+S22)*k/w(no)

c      write(110,*) (k),dreal(integrant(1))
c      write(111,*) (k),dimag(integrant(1))
c      write(112,*) (k),dreal(integrant(2))
c      write(113,*) (k),dimag(integrant(2))


      end


c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Imultiedessouscomp(theta,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision theta
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex k,k2,wronskien

      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex JB0,JB1,JB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0

c     write(*,*) 'routine',k02,kmax,hkmax,a,z,za,neps,nc,no     

c     calcul de k complex avec theta
      zr=kmax/2.d0*(dsin(theta)+1.d0)
      zi=hkmax*dcos(theta)
      k=zr+icomp*zi
c     calcul des fonctions de Bessel a l'ordre 0,1,2
      zr=zr*a
      zi=zi*a
      call ZBESJ(ZR, ZI, FNU, 1, N, CYR, CYI, NZ, IERR)     
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de bessel',NZ,IERR
         stop
      endif
      JB0=cyr(1)+icomp*cyi(1)
      JB1=cyr(2)+icomp*cyi(2)
      JB2=cyr(3)+icomp*cyi(3)
      k2=k*k
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)
         
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c     write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c     write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=0.d0
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-(1.d0,0.d0)

      mat_prod_11=mat_prod_sur11
      mat_prod_12=0.d0
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-(1.d0,0.d0)


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_21=mat_sca_21/mat_sca_11
      mat_sca_12=B_sca_12/mat_sca_11
      mat_sca_11=B_sca_11/mat_sca_11
      mat_sca_22=mat_sca_21*B_sca_12-B_sca_22
      mat_sca_21=mat_sca_21*B_sca_11-B_sca_21
      mat_prod_21=mat_prod_21/mat_prod_11
      mat_prod_21=mat_prod_21*B_prod_11-B_prod_21
      mat_prod_11=B_prod_11/mat_prod_11

c     remonte au point d'observation
      if (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c     write(*,*) 'dessus',S11,S12,S21,S22,P11,P21,ctmp,mat_sca_11
c     $        ,mat_sca_12
      elseif (no.gt.nc) then
c     au dessus du dipole
c     write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs11*mat_sca_11*ctmp
         P11=mat_prod_obs11*mat_prod_11*ctmp
         S12=mat_sca_obs11*mat_sca_12*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs21*mat_sca_11*ctmp
         S22=mat_sca_obs21*mat_sca_12*ctmp
         P21=mat_prod_obs21*mat_prod_11*ctmp
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c     write(*,*) 'identique'
c     ctmp=cdexp(icomp*w(no)*z)
c     S11=mat_sca_obs12*mat_sca_21*ctmp
c     P11=mat_prod_obs12*mat_prod_21*ctmp
c     S12=mat_sca_obs12*mat_sca_22*ctmp
c     ctmp=cdexp(-icomp*w(no)*z)
c     S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c     S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c     P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c     avec l'espace libre
c     S21=mat_sca_obs22*mat_sca_21*ctmp
c     S22=mat_sca_obs22*mat_sca_22*ctmp
c     P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c         ctmp=cdexp(icomp*w(no)*z)
c         S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c         S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c         P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
         S11=0.d0
         S12=0.d0
         P11=0.d0
         ctmp=cdexp(-icomp*w(no)*z)
         S21=(mat_sca_obs21*mat_sca_11)*ctmp
         S22=(mat_sca_obs21*mat_sca_12)*ctmp
         P21=(mat_prod_obs21*mat_prod_11)*ctmp
      endif


c     calcul element differentiel
      wronskien=kmax/2.d0*dcos(theta)-icomp*hkmax*dsin(theta)
c     write(*,*) 'wr',wronskien
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) S11,P11,S21,P21,S12,S22


      integrant(1)=k*(S11+P11+S21+P21)*wronskien*JB0
      integrant(2)=k*(S11-P11+S21-P21)*wronskien*JB2
      integrant(3)=(S12+S22)*wronskien*JB1
      integrant(4)=k2/w(no)*(-S11+S21)*wronskien*JB1
      integrant(5)=(-S12+S22)*k/w(no)*wronskien*JB0


c     write(10,*) dreal(theta),integrant(1)
c     write(11,*) dreal(theta),integrant(2)
c     write(12,*) dreal(theta),integrant(3)
c     write(13,*) dreal(theta),integrant(4)

      end

c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Imultiidessouscomp(k,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision k,k2
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

c     declarations liees a netlib
      integer n,nz
      parameter (n=3)
      double precision y,alp,xx
      dimension y(n)


      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      k2=k*k
      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c     write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c     write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=0.d0
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-(1.d0,0.d0)

      mat_prod_11=mat_prod_sur11
      mat_prod_12=0.d0
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-(1.d0,0.d0)


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c     write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_21=mat_sca_21/mat_sca_11
      mat_sca_12=B_sca_12/mat_sca_11
      mat_sca_11=B_sca_11/mat_sca_11
      mat_sca_22=mat_sca_21*B_sca_12-B_sca_22
      mat_sca_21=mat_sca_21*B_sca_11-B_sca_21
      mat_prod_21=mat_prod_21/mat_prod_11
      mat_prod_21=mat_prod_21*B_prod_11-B_prod_21
      mat_prod_11=B_prod_11/mat_prod_11

c     remonte au point d'observation
      if (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
c     write(*,*) 'dessus',S11,S12,S21,S22,P11,P21
c     write(*,*) 'ctmp',ctmp,mat_sca_11,mat_sca_12,mat_prod_11
c     write(*,*)'exp',cdexp(-icomp*w(nc)*za)*cdexp(icomp*w(neps+1)*z)
c     write(*,*) 'matsca12',nu_plus*w(nc)*k2
      elseif (no.gt.nc) then
c     au dessus du dipole
c     write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs11*mat_sca_11*ctmp
         P11=mat_prod_obs11*mat_prod_11*ctmp
         S12=mat_sca_obs11*mat_sca_12*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs21*mat_sca_11*ctmp
         S22=mat_sca_obs21*mat_sca_12*ctmp
         P21=mat_prod_obs21*mat_prod_11*ctmp
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c     write(*,*) 'identique'
c     ctmp=cdexp(icomp*w(no)*z)
c     S11=mat_sca_obs12*mat_sca_21*ctmp
c     P11=mat_prod_obs12*mat_prod_21*ctmp
c     S12=mat_sca_obs12*mat_sca_22*ctmp
c     ctmp=cdexp(-icomp*w(no)*z)
c     les oscillations viennent de la. Je somme deux choses qui sont
c     quasi identiques au signe pres. la somme est proche de zero et
c     entache d'erreur numerique:mat_sca_obs22*mat_sca_21+B_sca_21
c     S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c     S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c     P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c     write(99,9999) k,dreal(mat_sca_obs22*mat_sca_21+B_sca_21)

c     avec l'espace libre
c     S21=mat_sca_obs22*mat_sca_21*ctmp
c     S22=mat_sca_obs22*mat_sca_22*ctmp
c     P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c         ctmp=cdexp(icomp*w(no)*z)
c         S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c         S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c         P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
         S11=0.d0
         S12=0.d0
         P11=0.d0
         ctmp=cdexp(-icomp*w(no)*z)
         S21=(mat_sca_obs21*mat_sca_11)*ctmp
         S22=(mat_sca_obs21*mat_sca_12)*ctmp
         P21=(mat_prod_obs21*mat_prod_11)*ctmp
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) 'fin',S11,P11,S21,P21,S12,S22

      alp=0.d0
      xx=a*k
      call DBESJ (XX,ALP,N,Y,NZ)
      if (NZ.ne.0) then
         write(*,*) 'PB calcul de BESSEL'
         stop
      endif

c     write(99,9999) k,dreal(Y(1)),dreal(S11),dreal(P11),dreal(S21)
c     $     ,dreal(P21)
      integrant(1)=k*(S11+P11+S21+P21)*Y(1)
      integrant(2)=k*(S11-P11+S21-P21)*Y(3)
      integrant(3)=(S12+S22)*Y(2)
      integrant(4)=k2/w(no)*(-S11+S21)*Y(2)
      integrant(5)=(-S12+S22)*k/w(no)*Y(1)


c      write(70,*) k,integrant(1)
c      write(71,*) k,integrant(2)
c      write(72,*) k,integrant(3)
c      write(73,*) k,integrant(4)
 9999 format(201(d22.15,1x))
      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine ImultiiH_plusdessouscomp(kint,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision kint
      double complex S11,P11,S21,P21,S12,S22,icomp,k,k2,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex HB0,HB1,HB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0
      k=kmax+icomp*kint

c     premiere espece pour la fonction de Hankel
      zr=kmax*a
      zi=kint*a
      call ZBESH(ZR, ZI, FNU, 1, 1, N, CYR, CYI, NZ, IERR)
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de Hankel1',NZ,IERR
         write(*,*) 'argument',zr,zi
         stop
      endif
      HB0=cyr(1)+icomp*cyi(1)
      HB1=cyr(2)+icomp*cyi(2)
      HB2=cyr(3)+icomp*cyi(3)

      k2=k*k      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp


         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c     write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c     write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=0.d0
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-(1.d0,0.d0)

      mat_prod_11=mat_prod_sur11
      mat_prod_12=0.d0
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-(1.d0,0.d0)


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c     write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_21=mat_sca_21/mat_sca_11
      mat_sca_12=B_sca_12/mat_sca_11
      mat_sca_11=B_sca_11/mat_sca_11
      mat_sca_22=mat_sca_21*B_sca_12-B_sca_22
      mat_sca_21=mat_sca_21*B_sca_11-B_sca_21
      mat_prod_21=mat_prod_21/mat_prod_11
      mat_prod_21=mat_prod_21*B_prod_11-B_prod_21
      mat_prod_11=B_prod_11/mat_prod_11

c     remonte au point d'observation
      if (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
      elseif (no.gt.nc) then
c     au dessus du dipole
c     write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs11*mat_sca_11*ctmp
         P11=mat_prod_obs11*mat_prod_11*ctmp
         S12=mat_sca_obs11*mat_sca_12*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs21*mat_sca_11*ctmp
         S22=mat_sca_obs21*mat_sca_12*ctmp
         P21=mat_prod_obs21*mat_prod_11*ctmp
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c     write(*,*) 'identique'
c     ctmp=cdexp(icomp*w(no)*z)
c     S11=mat_sca_obs12*mat_sca_21*ctmp
c     P11=mat_prod_obs12*mat_prod_21*ctmp
c     S12=mat_sca_obs12*mat_sca_22*ctmp
c     ctmp=cdexp(-icomp*w(no)*z)

c     S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c     S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c     P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c     avec l'espace libre
c     S21=mat_sca_obs22*mat_sca_21*ctmp
c     S22=mat_sca_obs22*mat_sca_22*ctmp
c     P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c         ctmp=cdexp(icomp*w(no)*z)
c         S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c         S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c         P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
         S11=0.d0
         S12=0.d0
         P11=0.d0
         ctmp=cdexp(-icomp*w(no)*z)
         S21=(mat_sca_obs21*mat_sca_11)*ctmp
         S22=(mat_sca_obs21*mat_sca_12)*ctmp
         P21=(mat_prod_obs21*mat_prod_11)*ctmp
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) 'fin',S11,P11,S21,P21,S12,S22
      
      integrant(1)=k*(S11+P11+S21+P21)*HB0
      integrant(2)=k*(S11-P11+S21-P21)*HB2
      integrant(3)=(S12+S22)*HB1
      integrant(4)=k2/w(no)*(-S11+S21)*HB1
      integrant(5)=(-S12+S22)*k/w(no)*HB0
    
c      write(80,*) kint,integrant(1)
c      write(81,*) kint,integrant(2)
c      write(82,*) kint,integrant(3)
c      write(83,*) kint,integrant(4)

      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine ImultiiH_moinsdessouscomp(kint,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision kint
      double complex S11,P11,S21,P21,S12,S22,icomp,k,k2,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex HB0,HB1,HB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0
      k=kmax-icomp*kint

c     deuxieme espece pour la fonction de Hankel
      zr=kmax*a
      zi=-kint*a
      call ZBESH(ZR, ZI, FNU, 1, 2, N, CYR, CYI, NZ, IERR)
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de Hankel2',NZ,IERR
         write(*,*) 'argument',zr,zi
         stop
      endif
      HB0=cyr(1)+icomp*cyi(1)
      HB1=cyr(2)+icomp*cyi(2)
      HB2=cyr(3)+icomp*cyi(3)

      k2=k*k      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour interface dessus
      mat_sca_sur11=(1.d0,0.d0)
      mat_sca_sur12=0.d0
      mat_sca_sur21=0.d0
      mat_sca_sur22=(1.d0,0.d0)

      mat_prod_sur11=(1.d0,0.d0)
      mat_prod_sur12=0.d0
      mat_prod_sur21=0.d0
      mat_prod_sur22=(1.d0,0.d0)

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i+1)-w(i))*zcouche(i))
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))/ctmp
         ctmp=2.d0*w(i)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i+1)-w(i))/ctmp

         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins

         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sur11+mat_sca_12*mat_sca_sur21
         mat_sca_sur21=mat_sca_21*mat_sca_sur11+mat_sca_22
     $        *mat_sca_sur21
         mat_sca_sur11=ctmp
         ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
         mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $        *mat_sca_sur22
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12*mat_prod_sur21
         mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $        *mat_prod_sur21
         mat_prod_sur11=ctmp
         ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12*mat_prod_sur22
         mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $        *mat_prod_sur22
         mat_prod_sur12=ctmp

c     write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c     write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

c     evite le cas no=nc fait au dessus
         if (i.eq.no) then
            mat_sca_obs11=mat_sca_sur11
            mat_sca_obs12=mat_sca_sur12
            mat_sca_obs21=mat_sca_sur21
            mat_sca_obs22=mat_sca_sur22
            mat_prod_obs11=mat_prod_sur11
            mat_prod_obs12=mat_prod_sur12
            mat_prod_obs21=mat_prod_sur21
            mat_prod_obs22=mat_prod_sur22
         endif

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=0.d0
      mat_sca_21=mat_sca_sur21
      mat_sca_22=-(1.d0,0.d0)

      mat_prod_11=mat_prod_sur11
      mat_prod_12=0.d0
      mat_prod_21=mat_prod_sur21
      mat_prod_22=-(1.d0,0.d0)


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c     write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_21=mat_sca_21/mat_sca_11
      mat_sca_12=B_sca_12/mat_sca_11
      mat_sca_11=B_sca_11/mat_sca_11
      mat_sca_22=mat_sca_21*B_sca_12-B_sca_22
      mat_sca_21=mat_sca_21*B_sca_11-B_sca_21
      mat_prod_21=mat_prod_21/mat_prod_11
      mat_prod_21=mat_prod_21*B_prod_11-B_prod_21
      mat_prod_11=B_prod_11/mat_prod_11

c     remonte au point d'observation
      if (no.eq.neps+1) then
c     le point d'observation est dessus
         ctmp=cdexp(icomp*w(neps+1)*z)
         S11=mat_sca_11*ctmp
         S12=mat_sca_12*ctmp
         S21=0.d0
         S22=0.d0

         P11=mat_prod_11*ctmp
         P21=0.d0
      elseif (no.gt.nc) then
c     au dessus du dipole
c     write(*,*) 'sur',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs11*mat_sca_11*ctmp
         P11=mat_prod_obs11*mat_prod_11*ctmp
         S12=mat_sca_obs11*mat_sca_12*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs21*mat_sca_11*ctmp
         S22=mat_sca_obs21*mat_sca_12*ctmp
         P21=mat_prod_obs21*mat_prod_11*ctmp
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c     write(*,*) 'identique'
c     ctmp=cdexp(icomp*w(no)*z)
c     S11=mat_sca_obs12*mat_sca_21*ctmp
c     P11=mat_prod_obs12*mat_prod_21*ctmp
c     S12=mat_sca_obs12*mat_sca_22*ctmp
c     ctmp=cdexp(-icomp*w(no)*z)

c     S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c     S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c     P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

c     avec l'espace libre
c     S21=mat_sca_obs22*mat_sca_21*ctmp
c     S22=mat_sca_obs22*mat_sca_22*ctmp
c     P21=mat_prod_obs22*mat_prod_21*ctmp

c     identique normalement
c         ctmp=cdexp(icomp*w(no)*z)
c         S11=(mat_sca_obs11*mat_sca_11-B_sca_11)*ctmp
c         S12=(mat_sca_obs11*mat_sca_12-B_sca_12)*ctmp
c         P11=(mat_prod_obs11*mat_prod_11-B_prod_11)*ctmp
         S11=0.d0
         S12=0.d0
         P11=0.d0
         ctmp=cdexp(-icomp*w(no)*z)
         S21=(mat_sca_obs21*mat_sca_11)*ctmp
         S22=(mat_sca_obs21*mat_sca_12)*ctmp
         P21=(mat_prod_obs21*mat_prod_11)*ctmp
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) 'fin',S11,P11,S21,P21,S12,S22
      
      integrant(1)=k*(S11+P11+S21+P21)*HB0
      integrant(2)=k*(S11-P11+S21-P21)*HB2
      integrant(3)=(S12+S22)*HB1
      integrant(4)=k2/w(no)*(-S11+S21)*HB1
      integrant(5)=(-S12+S22)*k/w(no)*HB0
    

c     write(40,*) kint,integrant(1)
c     write(41,*) kint,integrant(2)
c     write(42,*) kint,integrant(3)
c     write(43,*) kint,integrant(4)

      end

c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Imultiedessuscomp(theta,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision theta
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex k,k2,wronskien

      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex JB0,JB1,JB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0

c     write(*,*) 'routine',k02,kmax,hkmax,a,z,za,neps,nc,no     

c     calcul de k complex avec theta
      zr=kmax/2.d0*(dsin(theta)+1.d0)
      zi=hkmax*dcos(theta)
      k=zr+icomp*zi
c     calcul des fonctions de Bessel a l'ordre 0,1,2
      zr=zr*a
      zi=zi*a
      call ZBESJ(ZR, ZI, FNU, 1, N, CYR, CYI, NZ, IERR)     
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de bessel',NZ,IERR
         stop
      endif
      JB0=cyr(1)+icomp*cyi(1)
      JB1=cyr(2)+icomp*cyi(2)
      JB2=cyr(3)+icomp*cyi(3)
      k2=k*k
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)
         
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp
c     write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c     write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c     write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c     write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c     stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     calcul de la matrice coupee
      mat_sca_11=(1.d0,0.d0)
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=0.d0
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=(1.d0,0.d0)
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=0.d0
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_12=-mat_sca_12/mat_sca_22
      mat_sca_11=B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=B_sca_12+mat_sca_12*B_sca_22
      mat_sca_21=B_sca_21/mat_sca_22
      mat_sca_22=B_sca_22/mat_sca_22

      mat_prod_12=-mat_prod_12/mat_prod_22
      mat_prod_22=1.d0/mat_prod_22
      mat_prod_21=mat_prod_22*B_prod_21
      mat_prod_11=B_prod_11+mat_prod_12*B_prod_21

c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.lt.nc) then
c     en dessous du dipole
c     write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22

         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs22*mat_sca_21*ctmp
         S22=mat_sca_obs22*mat_sca_22*ctmp
         P21=mat_prod_obs22*mat_prod_21*ctmp

      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c     write(*,*) 'identique'
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
c         ctmp=cdexp(-icomp*w(no)*z)
c         S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c         S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c         P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp
         S21=0.d0
         S22=0.d0
         P21=0.d0
      endif


c     calcul element differentiel
      wronskien=kmax/2.d0*dcos(theta)-icomp*hkmax*dsin(theta)
c     write(*,*) 'wr',wronskien
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) S11,P11,S21,P21,S12,S22


      integrant(1)=k*(S11+P11+S21+P21)*wronskien*JB0
      integrant(2)=k*(S11-P11+S21-P21)*wronskien*JB2
      integrant(3)=(S12+S22)*wronskien*JB1
      integrant(4)=k2/w(no)*(-S11+S21)*wronskien*JB1
      integrant(5)=(-S12+S22)*k/w(no)*wronskien*JB0


c     write(10,*) dreal(theta),integrant(1)
c     write(11,*) dreal(theta),integrant(2)
c     write(12,*) dreal(theta),integrant(3)
c     write(13,*) dreal(theta),integrant(4)

      end

c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Imultiidessuscomp(k,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision k,k2
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

c     declarations liees a netlib
      integer n,nz
      parameter (n=3)
      double precision y,alp,xx
      dimension y(n)


      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      k2=k*k
      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     write(*,*) '**************',k

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp

c     write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c     write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     write(*,*) 'prod',mat_prod_11,mat_prod_12 ,mat_prod_21
c     $        ,mat_prod_22
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c     write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c     write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c     stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo

      mat_sca_11=(1.d0,0.d0)
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=0.d0
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=(1.d0,0.d0)
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=0.d0
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c     write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B

      mat_sca_12=-mat_sca_12/mat_sca_22
      mat_sca_11=B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=B_sca_12+mat_sca_12*B_sca_22
      mat_sca_21=B_sca_21/mat_sca_22
      mat_sca_22=B_sca_22/mat_sca_22

      mat_prod_12=-mat_prod_12/mat_prod_22
      mat_prod_22=1.d0/mat_prod_22
      mat_prod_21=mat_prod_22*B_prod_21
      mat_prod_11=B_prod_11+mat_prod_12*B_prod_21


c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.lt.nc) then
c     en dessous du dipole
c     write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs22*mat_sca_21*ctmp
         S22=mat_sca_obs22*mat_sca_22*ctmp
         P21=mat_prod_obs22*mat_prod_21*ctmp
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c     write(*,*) 'identique'
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
c         ctmp=cdexp(-icomp*w(no)*z)
c     les oscillations viennent de la. Je somme deux choses qui sont
c     quasi identiques au signe pres. la somme est proche de zero et
c     entache d'erreur numerique:mat_sca_obs22*mat_sca_21+B_sca_21
c         S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c         S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c         P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp
         S21=0.d0
         S22=0.d0
         P21=0.d0
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) 'fin',S11,P11,S21,P21,S12,S22

      alp=0.d0
      xx=a*k
      call DBESJ (XX,ALP,N,Y,NZ)
      if (NZ.ne.0) then
         write(*,*) 'PB calcul de BESSEL'
         stop
      endif

c     write(99,9999) k,dreal(Y(1)),dreal(S11),dreal(P11),dreal(S21)
c     $     ,dreal(P21)
      integrant(1)=k*(S11+P11+S21+P21)*Y(1)
      integrant(2)=k*(S11-P11+S21-P21)*Y(3)
      integrant(3)=(S12+S22)*Y(2)
      integrant(4)=k2/w(no)*(-S11+S21)*Y(2)
      integrant(5)=(-S12+S22)*k/w(no)*Y(1)


c     write(20,*) k,integrant(1)
c     write(21,*) k,integrant(2)
c     write(22,*) k,integrant(3)
c     write(23,*) k,integrant(4)
 9999 format(201(d22.15,1x))
      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine ImultiiH_plusdessuscomp(kint,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision kint
      double complex S11,P11,S21,P21,S12,S22,icomp,k,k2,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex HB0,HB1,HB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0
      k=kmax+icomp*kint

c     premiere espece pour la fonction de Hankel
      zr=kmax*a
      zi=kint*a
      call ZBESH(ZR, ZI, FNU, 1, 1, N, CYR, CYI, NZ, IERR)
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de Hankel1',NZ,IERR
         write(*,*) 'argument',zr,zi
         stop
      endif
      HB0=cyr(1)+icomp*cyi(1)
      HB1=cyr(2)+icomp*cyi(2)
      HB2=cyr(3)+icomp*cyi(3)

      k2=k*k      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp

c     write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c     write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     write(*,*) 'prod',mat_prod_11,mat_prod_12 ,mat_prod_21
c     $        ,mat_prod_22
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c     write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c     write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c     stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     calcul de la matrice coupee
      mat_sca_11=(1.d0,0.d0)
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=0.d0
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=(1.d0,0.d0)
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=0.d0
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c     write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_12=-mat_sca_12/mat_sca_22
      mat_sca_11=B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=B_sca_12+mat_sca_12*B_sca_22
      mat_sca_21=B_sca_21/mat_sca_22
      mat_sca_22=B_sca_22/mat_sca_22

      mat_prod_12=-mat_prod_12/mat_prod_22
      mat_prod_22=1.d0/mat_prod_22
      mat_prod_21=mat_prod_22*B_prod_21
      mat_prod_11=B_prod_11+mat_prod_12*B_prod_21

c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.lt.nc) then
c     en dessous du dipole
c     write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs22*mat_sca_21*ctmp
         S22=mat_sca_obs22*mat_sca_22*ctmp
         P21=mat_prod_obs22*mat_prod_21*ctmp
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c     write(*,*) 'identique'
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
c         ctmp=cdexp(-icomp*w(no)*z)
c         S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c         S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c         P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp

         S21=0.d0
         S22=0.d0
         P21=0.d0
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) 'fin',S11,P11,S21,P21,S12,S22
      
      integrant(1)=k*(S11+P11+S21+P21)*HB0
      integrant(2)=k*(S11-P11+S21-P21)*HB2
      integrant(3)=(S12+S22)*HB1
      integrant(4)=k2/w(no)*(-S11+S21)*HB1
      integrant(5)=(-S12+S22)*k/w(no)*HB0


c     write(30,*) kint,integrant(1)
c     write(31,*) kint,integrant(2)
c     write(32,*) kint,integrant(3)
c     write(33,*) kint,integrant(4)

      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine ImultiiH_moinsdessuscomp(kint,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision kint
      double complex S11,P11,S21,P21,S12,S22,icomp,k,k2,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex HB0,HB1,HB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
      icomp=(0.d0,1.d0)
      FNU=0.d0
      k=kmax-icomp*kint

c     deuxieme espece pour la fonction de Hankel
      zr=kmax*a
      zi=-kint*a
      call ZBESH(ZR, ZI, FNU, 1, 2, N, CYR, CYI, NZ, IERR)
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de Hankel2',NZ,IERR
         write(*,*) 'argument',zr,zi
         stop
      endif
      HB0=cyr(1)+icomp*cyi(1)
      HB1=cyr(2)+icomp*cyi(2)
      HB2=cyr(3)+icomp*cyi(3)

      k2=k*k      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c     write(*,*) 'www',w(i),k,kmax
c     write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=(1.d0,0.d0)
      mat_sca_sous12=0.d0
      mat_sca_sous21=0.d0
      mat_sca_sous22=(1.d0,0.d0)

      mat_prod_sous11=(1.d0,0.d0)
      mat_prod_sous12=0.d0
      mat_prod_sous21=0.d0
      mat_prod_sous22=(1.d0,0.d0)

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=cdexp(-icomp*(w(i)+w(i+1))*zcouche(i))
         e_moins=cdexp(icomp*(w(i)-w(i+1))*zcouche(i))
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))/ctmp
         deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))/ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=(w(i+1)+w(i))/ctmp
         deltas_moins=(w(i)-w(i+1))/ctmp

c     write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $        ,deltas_plus,deltas_moins
         
         mat_sca_11=e_moins*deltap_plus
         mat_sca_12=e_plus*deltap_moins
         mat_sca_21=deltap_moins/e_plus
         mat_sca_22=deltap_plus/e_moins
c     write(*,*) 'sca',mat_sca_11,mat_sca_12 ,mat_sca_21,mat_sca_22
         mat_prod_11=e_moins*deltas_plus
         mat_prod_12=-e_plus*deltas_moins
         mat_prod_21=-deltas_moins/e_plus
         mat_prod_22=deltas_plus/e_moins
c     write(*,*) 'prod',mat_prod_11,mat_prod_12 ,mat_prod_21
c     $        ,mat_prod_22
c     produit des deux matrices; scalaire
         ctmp=mat_sca_11*mat_sca_sous11+mat_sca_12*mat_sca_sous21
         mat_sca_sous21=mat_sca_21*mat_sca_sous11+mat_sca_22
     $        *mat_sca_sous21
         mat_sca_sous11=ctmp
         ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
         mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $        *mat_sca_sous22
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12*mat_prod_sous21
         mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $        *mat_prod_sous21
         mat_prod_sous11=ctmp
         ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12*mat_prod_sous22
         mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $        *mat_prod_sous22
         mat_prod_sous12=ctmp

c     write(*,*) 'soussca',mat_sca_sous11,mat_sca_sous12
c     $        ,mat_sca_sous21,mat_sca_sous22
c     write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c     stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo

c     calcul de la matrice coupee
      mat_sca_11=(1.d0,0.d0)
      mat_sca_12=-mat_sca_sous12
      mat_sca_21=0.d0
      mat_sca_22=-mat_sca_sous22

      mat_prod_11=(1.d0,0.d0)
      mat_prod_12=-mat_prod_sous12
      mat_prod_21=0.d0
      mat_prod_22=-mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-cdexp(-icomp*w(nc)*za)
      nu_moins=-cdexp(icomp*w(nc)*za)

c     write(*,*) 'nu',nu_plus,nu_moins

      B_sca_11=-nu_plus*w(nc)
      B_sca_12=nu_plus*k2
      B_sca_21=nu_moins*w(nc)
      B_sca_22=nu_moins*k2

      B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
      B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
      mat_sca_12=-mat_sca_12/mat_sca_22
      mat_sca_11=B_sca_11+mat_sca_12*B_sca_21
      mat_sca_12=B_sca_12+mat_sca_12*B_sca_22
      mat_sca_21=B_sca_21/mat_sca_22
      mat_sca_22=B_sca_22/mat_sca_22

      mat_prod_12=-mat_prod_12/mat_prod_22
      mat_prod_22=1.d0/mat_prod_22
      mat_prod_21=mat_prod_22*B_prod_21
      mat_prod_11=B_prod_11+mat_prod_12*B_prod_21


c     remonte au point d'observation
      if (no.eq.0) then
c     le point d'observation est dessous
         ctmp=cdexp(-icomp*w(0)*z)
         S11=0.d0
         S12=0.d0
         S21=mat_sca_21*ctmp
         S22=mat_sca_22*ctmp

         P11=0.d0
         P21=mat_prod_21*ctmp
      elseif (no.lt.nc) then
c     en dessous du dipole
c     write(*,*) 'sous',mat_sca_obs12,mat_sca_obs22,mat_sca_21
c     $           ,mat_sca_22
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
         ctmp=cdexp(-icomp*w(no)*z)
         S21=mat_sca_obs22*mat_sca_21*ctmp
         S22=mat_sca_obs22*mat_sca_22*ctmp
         P21=mat_prod_obs22*mat_prod_21*ctmp
      elseif (no.eq.nc) then
c     dans la meme couche que le dipole
c     write(*,*) 'identique'
         ctmp=cdexp(icomp*w(no)*z)
         S11=mat_sca_obs12*mat_sca_21*ctmp
         P11=mat_prod_obs12*mat_prod_21*ctmp
         S12=mat_sca_obs12*mat_sca_22*ctmp
c         ctmp=cdexp(-icomp*w(no)*z)
c         S21=(mat_sca_obs22*mat_sca_21+B_sca_21)*ctmp
c         S22=(mat_sca_obs22*mat_sca_22+B_sca_22)*ctmp
c         P21=(mat_prod_obs22*mat_prod_21+B_prod_21)*ctmp
         S21=0.d0
         S22=0.d0
         P21=0.d0
         
      endif


c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c     write(*,*) 'fin',S11,P11,S21,P21,S12,S22
      
      integrant(1)=k*(S11+P11+S21+P21)*HB0
      integrant(2)=k*(S11-P11+S21-P21)*HB2
      integrant(3)=(S12+S22)*HB1
      integrant(4)=k2/w(no)*(-S11+S21)*HB1
      integrant(5)=(-S12+S22)*k/w(no)*HB0


c     write(40,*) kint,integrant(1)
c     write(41,*) kint,integrant(2)
c     write(42,*) kint,integrant(3)
c     write(43,*) kint,integrant(4)

      end
c     partie log
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Idiagmultiilogcomp(k,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision k,k2,pi
      double complex S11,P11,S21,P21,S12,S22,icomp,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp,ctmp1,ctmp2
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det,zero

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1),sommearg
     $     ,uncomp,icomppi

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
      common/datalog/icomp,icomppi,uncomp,zero
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
!$OMP THREADPRIVATE(/datalog/)
      k2=k*k
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=0.d0
      mat_sca_sous12=zero
      mat_sca_sous21=zero
      mat_sca_sous22=0.d0

      mat_prod_sous11=0.d0
      mat_prod_sous12=zero
      mat_prod_sous21=zero
      mat_prod_sous22=0.d0

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=-icomp*(w(i)+w(i+1))*zcouche(i)
c         write(*,*) 'eplus',e_plus,icomp*(w(i)+w(i+1))*zcouche(i),w(i)
c     $        ,w(i+1),zcouche(i)
         e_moins=icomp*(w(i)-w(i+1))*zcouche(i)
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=cdlog((epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $        /ctmp)
         deltap_moins=cdlog((epscouche(i+1)*w(i)-epscouche(i)*w(i+1))
     $        /ctmp)
c         write(*,*) 'delta',deltap_moins,epscouche(i+1)*w(i),epscouche(i
c     $        )*w(i+1),ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=cdlog((w(i+1)+w(i))/ctmp)
         deltas_moins=cdlog((w(i)-w(i+1))/ctmp)
c         write(*,*) 'ddd',cdexp(e_plus),cdexp(e_moins),cdexp(deltap_plus
c     $        ),cdexp(deltap_moins),cdexp(deltas_plus)
c     $        ,cdexp(deltas_moins)
         
         mat_sca_11=e_moins+deltap_plus
         mat_sca_12=e_plus+deltap_moins
         mat_sca_21=deltap_moins-e_plus
         mat_sca_22=deltap_plus-e_moins
c         write(*,*) 'sca',(mat_sca_11),(mat_sca_12)
c     $        ,(mat_sca_21),(mat_sca_22)
         mat_prod_11=e_moins+deltas_plus
         mat_prod_12=icomppi+e_plus+deltas_moins
         mat_prod_21=icomppi+deltas_moins-e_plus
         mat_prod_22=deltas_plus-e_moins
c     produit des deux matrices; scalaire
         ctmp1=mat_sca_11+mat_sca_sous11
         ctmp2=mat_sca_12+mat_sca_sous21
         ctmp=sommearg(ctmp1,ctmp2)
c         write(*,*) 'bbbb',ctmp,ctmp1,ctmp2
         ctmp1=mat_sca_21+mat_sca_sous11
         ctmp2=mat_sca_22+mat_sca_sous21
         mat_sca_sous21=sommearg(ctmp1,ctmp2)
         mat_sca_sous11=ctmp

         ctmp1=mat_sca_11+mat_sca_sous12
         ctmp2=mat_sca_12+mat_sca_sous22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sous12
         ctmp2=mat_sca_22+mat_sca_sous22
         mat_sca_sous22=sommearg(ctmp1,ctmp2)
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp1=mat_prod_11+mat_prod_sous11
         ctmp2=mat_prod_12+mat_prod_sous21
         ctmp=sommearg(ctmp1,ctmp2)
         
         ctmp1=mat_prod_21+mat_prod_sous11
         ctmp2=mat_prod_22+mat_prod_sous21
         mat_prod_sous21=sommearg(ctmp1,ctmp2)
         mat_prod_sous11=ctmp

         ctmp1=mat_prod_11+mat_prod_sous12
         ctmp2=mat_prod_12+mat_prod_sous22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sous12
         ctmp2=mat_prod_22+mat_prod_sous22
         mat_prod_sous22=sommearg(ctmp1,ctmp2)
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',cdexp(mat_sca_sous11),cdexp(mat_sca_sous12
c     $        ),cdexp(mat_sca_sous21),cdexp(mat_sca_sous22)
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
c            write(*,*) 'obs', mat_sca_obs22,cdexp(mat_sca_obs22)
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=0.d0
      mat_sca_sur12=zero
      mat_sca_sur21=zero
      mat_sca_sur22=0.d0

      mat_prod_sur11=0.d0
      mat_prod_sur12=zero
      mat_prod_sur21=zero
      mat_prod_sur22=0.d0

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=-icomp*(w(i)+w(i+1))*zcouche(i)
         e_moins=icomp*(w(i+1)-w(i))*zcouche(i)
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=cdlog((epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $        /ctmp)
         deltap_moins=cdlog((epscouche(i)*w(i+1)-epscouche(i+1)*w(i))
     $        /ctmp)
         ctmp=2.d0*w(i)
         deltas_plus=cdlog((w(i+1)+w(i))/ctmp)
         deltas_moins=cdlog((w(i+1)-w(i))/ctmp)

         mat_sca_11=e_moins+deltap_plus
         mat_sca_12=e_plus+deltap_moins
         mat_sca_21=deltap_moins-e_plus
         mat_sca_22=deltap_plus-e_moins

         mat_prod_11=e_moins+deltas_plus
         mat_prod_12=icomppi+e_plus+deltas_moins
         mat_prod_21=icomppi+deltas_moins-e_plus
         mat_prod_22=deltas_plus-e_moins
c     produit des deux matrices; scalaire
         ctmp1=mat_sca_11+mat_sca_sur11
         ctmp2=mat_sca_12+mat_sca_sur21
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sur11
         ctmp2=mat_sca_22+mat_sca_sur21
         mat_sca_sur21=sommearg(ctmp1,ctmp2)
         mat_sca_sur11=ctmp

         ctmp1=mat_sca_11+mat_sca_sur12
         ctmp2=mat_sca_12+mat_sca_sur22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sur12
         ctmp2=mat_sca_22+mat_sca_sur22
         mat_sca_sur22=sommearg(ctmp1,ctmp2)
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp1=mat_prod_11+mat_prod_sur11
         ctmp2=mat_prod_12+mat_prod_sur21
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sur11
         ctmp2=mat_prod_22+mat_prod_sur21
         mat_prod_sur21=sommearg(ctmp1,ctmp2)
         mat_prod_sur11=ctmp

         ctmp1=mat_prod_11+mat_prod_sur12
         ctmp2=mat_prod_12+mat_prod_sur22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sur12
         ctmp2=mat_prod_22+mat_prod_sur22
         mat_prod_sur22=sommearg(ctmp1,ctmp2)
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22


      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=icomppi+mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=icomppi+mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=icomppi+mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=icomppi+mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-icomp*w(nc)*za+icomppi
      nu_moins=icomp*w(nc)*za+icomppi


      ctmp1=cdlog(w(nc))
      ctmp2=dlog(k2)
      B_sca_11=icomppi+nu_plus+ctmp1
      B_sca_12=nu_plus+ctmp2
      B_sca_21=nu_moins+ctmp1
      B_sca_22=nu_moins+ctmp2
      ctmp=cdlog(epscouche(nc)*k02/w(nc))
c      write(*,*) 'B_sca_22',B_sca_22,cdexp(B_sca_22)
      B_prod_11=icomppi+nu_plus+ctmp
      B_prod_21=nu_moins+ctmp

c     inverse de la matrice coupee et produit par B
      ctmp1=mat_sca_11+mat_sca_22
      ctmp2=mat_sca_12+mat_sca_21+icomppi
      det=sommearg(ctmp1,ctmp2)
c      write(*,*) 'det',(det)
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22-det
      mat_sca_22=ctmp-det
      mat_sca_12=icomppi+mat_sca_12-det
      mat_sca_21=icomppi+mat_sca_21-det

      ctmp1=mat_sca_11+B_sca_11
      ctmp2=mat_sca_12+B_sca_21
      ctmp=sommearg(ctmp1,ctmp2)
      ctmp1=mat_sca_11+B_sca_12
      ctmp2=mat_sca_12+B_sca_22
      mat_sca_12=sommearg(ctmp1,ctmp2)
      mat_sca_11=ctmp

      ctmp1=mat_sca_21+B_sca_11
      ctmp2=mat_sca_22+B_sca_21
      ctmp=sommearg(ctmp1,ctmp2)
c      write(*,*) 'aaa',ctmp,ctmp1,ctmp2
      ctmp1=mat_sca_21+B_sca_12
      ctmp2=mat_sca_22+B_sca_22
      mat_sca_22=sommearg(ctmp1,ctmp2)
      mat_sca_21=ctmp


      ctmp1=mat_prod_11+mat_prod_22
      ctmp2=icomppi+mat_prod_12+mat_prod_21
      det=sommearg(ctmp1,ctmp2)
c      write(*,*) 'det',(det)
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22-det
      mat_prod_22=ctmp-det
      mat_prod_12=icomppi+mat_prod_12-det
      mat_prod_21=icomppi+mat_prod_21-det

      ctmp1=mat_prod_21+B_prod_11
      ctmp2=mat_prod_22+B_prod_21
      mat_prod_21=sommearg(ctmp1,ctmp2)
      ctmp1=mat_prod_11+B_prod_11
      ctmp2=mat_prod_12+B_prod_21
      mat_prod_11=sommearg(ctmp1,ctmp2)

      ctmp=icomp*w(no)*z

      S11=mat_sca_obs12+mat_sca_21+ctmp
      P11=mat_prod_obs12+mat_prod_21+ctmp
      S12=mat_sca_obs12+mat_sca_22+ctmp
      ctmp=-icomp*w(no)*z
      ctmp1=mat_sca_obs22+mat_sca_21+ctmp
      ctmp2=B_sca_21+ctmp
      S21=sommearg(ctmp1,ctmp2)
      ctmp1=mat_sca_obs22+mat_sca_22+ctmp
c     write(*,*) 'rr',cdexp(mat_sca_obs22),cdexp(mat_sca_22)
      ctmp2=B_sca_22+ctmp
      S22=sommearg(ctmp1,ctmp2)
c     write(*,*) 'S22',cdexp(S22),cdexp(ctmp1),cdexp(ctmp2)
c     $           ,cdexp(ctmp)
      ctmp1=mat_prod_obs22+mat_prod_21+ctmp
      ctmp2=B_prod_21+ctmp
      P21=sommearg(ctmp1,ctmp2)

c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22
c      write(*,*) 'SS',S11,P11,S21,P21,S12,S22
      S11=cdexp(S11)
      S12=cdexp(S12)
      S21=cdexp(S21)
      S22=cdexp(S22)
      P11=cdexp(P11)
      P21=cdexp(P21)
c      write(*,*) 'SSS',S11,P11,S21,P21,S12,S22

      integrant(1)=k*(S11+P11+S21+P21)
      integrant(2)=(-S12+S22)*k/w(no)

c      write(70,*) k,integrant(1)
c      write(71,*) k,integrant(2)
c      write(72,*) k,integrant(3)
c      write(73,*) k,integrant(4)
      end

c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine Imultiilogcomp(k,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision k,k2
      double complex S11,P11,S21,P21,S12,S22,icomp,icomppi,uncomp,zero
     $     ,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp,ctmp1,ctmp2,sommearg
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

c     declarations liees a netlib
      integer n,nz
      parameter (n=3)
      double precision y,alp,xx
      dimension y(n)


      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
      common/datalog/icomp,icomppi,uncomp,zero
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
!$OMP THREADPRIVATE(/datalog/)

      k2=k*k
      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo
c     initialise pour l'interface 0    
      mat_sca_sous11=0.d0
      mat_sca_sous12=zero
      mat_sca_sous21=zero
      mat_sca_sous22=0.d0

      mat_prod_sous11=0.d0
      mat_prod_sous12=zero
      mat_prod_sous21=zero
      mat_prod_sous22=0.d0

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=-icomp*(w(i)+w(i+1))*zcouche(i)
c         write(*,*) 'eplus',e_plus,icomp*(w(i)+w(i+1))*zcouche(i),w(i)
c     $        ,w(i+1),zcouche(i)
         e_moins=icomp*(w(i)-w(i+1))*zcouche(i)
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=cdlog((epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $        /ctmp)
         deltap_moins=cdlog((epscouche(i+1)*w(i)-epscouche(i)*w(i+1))
     $        /ctmp)
c         write(*,*) 'delta',deltap_moins,epscouche(i+1)*w(i),epscouche(i
c     $        )*w(i+1),ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=cdlog((w(i+1)+w(i))/ctmp)
         deltas_moins=cdlog((w(i)-w(i+1))/ctmp)
c         write(*,*) 'ddd',cdexp(e_plus),cdexp(e_moins),cdexp(deltap_plus
c     $        ),cdexp(deltap_moins),cdexp(deltas_plus)
c     $        ,cdexp(deltas_moins)
         
         mat_sca_11=e_moins+deltap_plus
         mat_sca_12=e_plus+deltap_moins
         mat_sca_21=deltap_moins-e_plus
         mat_sca_22=deltap_plus-e_moins
c         write(*,*) 'sca',(mat_sca_11),(mat_sca_12)
c     $        ,(mat_sca_21),(mat_sca_22)
         mat_prod_11=e_moins+deltas_plus
         mat_prod_12=icomppi+e_plus+deltas_moins
         mat_prod_21=icomppi+deltas_moins-e_plus
         mat_prod_22=deltas_plus-e_moins
c     produit des deux matrices; scalaire
         ctmp1=mat_sca_11+mat_sca_sous11
         ctmp2=mat_sca_12+mat_sca_sous21
         ctmp=sommearg(ctmp1,ctmp2)
c         write(*,*) 'bbbb',ctmp,ctmp1,ctmp2
         ctmp1=mat_sca_21+mat_sca_sous11
         ctmp2=mat_sca_22+mat_sca_sous21
         mat_sca_sous21=sommearg(ctmp1,ctmp2)
         mat_sca_sous11=ctmp

         ctmp1=mat_sca_11+mat_sca_sous12
         ctmp2=mat_sca_12+mat_sca_sous22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sous12
         ctmp2=mat_sca_22+mat_sca_sous22
         mat_sca_sous22=sommearg(ctmp1,ctmp2)
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp1=mat_prod_11+mat_prod_sous11
         ctmp2=mat_prod_12+mat_prod_sous21
         ctmp=sommearg(ctmp1,ctmp2)
         
         ctmp1=mat_prod_21+mat_prod_sous11
         ctmp2=mat_prod_22+mat_prod_sous21
         mat_prod_sous21=sommearg(ctmp1,ctmp2)
         mat_prod_sous11=ctmp

         ctmp1=mat_prod_11+mat_prod_sous12
         ctmp2=mat_prod_12+mat_prod_sous22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sous12
         ctmp2=mat_prod_22+mat_prod_sous22
         mat_prod_sous22=sommearg(ctmp1,ctmp2)
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',cdexp(mat_sca_sous11),cdexp(mat_sca_sous12
c     $        ),cdexp(mat_sca_sous21),cdexp(mat_sca_sous22)
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
c            write(*,*) 'obs', mat_sca_obs22,cdexp(mat_sca_obs22)
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=0.d0
      mat_sca_sur12=zero
      mat_sca_sur21=zero
      mat_sca_sur22=0.d0

      mat_prod_sur11=0.d0
      mat_prod_sur12=zero
      mat_prod_sur21=zero
      mat_prod_sur22=0.d0

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=-icomp*(w(i)+w(i+1))*zcouche(i)
         e_moins=icomp*(w(i+1)-w(i))*zcouche(i)
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=cdlog((epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $        /ctmp)
         deltap_moins=cdlog((epscouche(i)*w(i+1)-epscouche(i+1)*w(i))
     $        /ctmp)
         ctmp=2.d0*w(i)
         deltas_plus=cdlog((w(i+1)+w(i))/ctmp)
         deltas_moins=cdlog((w(i+1)-w(i))/ctmp)

         mat_sca_11=e_moins+deltap_plus
         mat_sca_12=e_plus+deltap_moins
         mat_sca_21=deltap_moins-e_plus
         mat_sca_22=deltap_plus-e_moins

         mat_prod_11=e_moins+deltas_plus
         mat_prod_12=icomppi+e_plus+deltas_moins
         mat_prod_21=icomppi+deltas_moins-e_plus
         mat_prod_22=deltas_plus-e_moins
c     produit des deux matrices; scalaire
         ctmp1=mat_sca_11+mat_sca_sur11
         ctmp2=mat_sca_12+mat_sca_sur21
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sur11
         ctmp2=mat_sca_22+mat_sca_sur21
         mat_sca_sur21=sommearg(ctmp1,ctmp2)
         mat_sca_sur11=ctmp

         ctmp1=mat_sca_11+mat_sca_sur12
         ctmp2=mat_sca_12+mat_sca_sur22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sur12
         ctmp2=mat_sca_22+mat_sca_sur22
         mat_sca_sur22=sommearg(ctmp1,ctmp2)
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp1=mat_prod_11+mat_prod_sur11
         ctmp2=mat_prod_12+mat_prod_sur21
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sur11
         ctmp2=mat_prod_22+mat_prod_sur21
         mat_prod_sur21=sommearg(ctmp1,ctmp2)
         mat_prod_sur11=ctmp

         ctmp1=mat_prod_11+mat_prod_sur12
         ctmp2=mat_prod_12+mat_prod_sur22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sur12
         ctmp2=mat_prod_22+mat_prod_sur22
         mat_prod_sur22=sommearg(ctmp1,ctmp2)
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22


      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=icomppi+mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=icomppi+mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=icomppi+mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=icomppi+mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-icomp*w(nc)*za+icomppi
      nu_moins=icomp*w(nc)*za+icomppi


      ctmp1=cdlog(w(nc))
      ctmp2=dlog(k2)
      B_sca_11=icomppi+nu_plus+ctmp1
      B_sca_12=nu_plus+ctmp2
      B_sca_21=nu_moins+ctmp1
      B_sca_22=nu_moins+ctmp2
      ctmp=cdlog(epscouche(nc)*k02/w(nc))
c      write(*,*) 'B_sca_22',B_sca_22,cdexp(B_sca_22)
      B_prod_11=icomppi+nu_plus+ctmp
      B_prod_21=nu_moins+ctmp

c     inverse de la matrice coupee et produit par B
      ctmp1=mat_sca_11+mat_sca_22
      ctmp2=mat_sca_12+mat_sca_21+icomppi
      det=sommearg(ctmp1,ctmp2)
c      write(*,*) 'det',(det)
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22-det
      mat_sca_22=ctmp-det
      mat_sca_12=icomppi+mat_sca_12-det
      mat_sca_21=icomppi+mat_sca_21-det

      ctmp1=mat_sca_11+B_sca_11
      ctmp2=mat_sca_12+B_sca_21
      ctmp=sommearg(ctmp1,ctmp2)
      ctmp1=mat_sca_11+B_sca_12
      ctmp2=mat_sca_12+B_sca_22
      mat_sca_12=sommearg(ctmp1,ctmp2)
      mat_sca_11=ctmp

      ctmp1=mat_sca_21+B_sca_11
      ctmp2=mat_sca_22+B_sca_21
      ctmp=sommearg(ctmp1,ctmp2)
c      write(*,*) 'aaa',ctmp,ctmp1,ctmp2
      ctmp1=mat_sca_21+B_sca_12
      ctmp2=mat_sca_22+B_sca_22
      mat_sca_22=sommearg(ctmp1,ctmp2)
      mat_sca_21=ctmp


      ctmp1=mat_prod_11+mat_prod_22
      ctmp2=icomppi+mat_prod_12+mat_prod_21
      det=sommearg(ctmp1,ctmp2)
c      write(*,*) 'det',(det)
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22-det
      mat_prod_22=ctmp-det
      mat_prod_12=icomppi+mat_prod_12-det
      mat_prod_21=icomppi+mat_prod_21-det

      ctmp1=mat_prod_21+B_prod_11
      ctmp2=mat_prod_22+B_prod_21
      mat_prod_21=sommearg(ctmp1,ctmp2)
      ctmp1=mat_prod_11+B_prod_11
      ctmp2=mat_prod_12+B_prod_21
      mat_prod_11=sommearg(ctmp1,ctmp2)

      ctmp=icomp*w(no)*z

      S11=mat_sca_obs12+mat_sca_21+ctmp
      P11=mat_prod_obs12+mat_prod_21+ctmp
      S12=mat_sca_obs12+mat_sca_22+ctmp
      ctmp=-icomp*w(no)*z
      ctmp1=mat_sca_obs22+mat_sca_21+ctmp
      ctmp2=B_sca_21+ctmp
      S21=sommearg(ctmp1,ctmp2)
      ctmp1=mat_sca_obs22+mat_sca_22+ctmp
c     write(*,*) 'rr',cdexp(mat_sca_obs22),cdexp(mat_sca_22)
      ctmp2=B_sca_22+ctmp
      S22=sommearg(ctmp1,ctmp2)
c     write(*,*) 'S22',cdexp(S22),cdexp(ctmp1),cdexp(ctmp2)
c     $           ,cdexp(ctmp)
      ctmp1=mat_prod_obs22+mat_prod_21+ctmp
      ctmp2=B_prod_21+ctmp
      P21=sommearg(ctmp1,ctmp2)

c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22
c      write(*,*) 'SS',S11,P11,S21,P21,S12,S22
      S11=cdexp(S11)
      S12=cdexp(S12)
      S21=cdexp(S21)
      S22=cdexp(S22)
      P11=cdexp(P11)
      P21=cdexp(P21)

      alp=0.d0
      xx=a*k
      call DBESJ (XX,ALP,N,Y,NZ)
      if (NZ.ne.0) then
         write(*,*) 'PB calcul de BESSEL'
         stop
      endif

c      write(99,9999) k,dreal(Y(1)),dreal(S11),dreal(P11),dreal(S21)
c     $     ,dreal(P21)
      integrant(1)=k*(S11+P11+S21+P21)*Y(1)
      integrant(2)=k*(S11-P11+S21-P21)*Y(3)
      integrant(3)=(S12+S22)*Y(2)
      integrant(4)=k2/w(no)*(-S11+S21)*Y(2)
      integrant(5)=(-S12+S22)*k/w(no)*Y(1)

c      write(70,*) k,integrant(1)
c      write(71,*) k,integrant(2)
c      write(72,*) k,integrant(3)
c      write(73,*) k,integrant(4)
 9999 format(201(d22.15,1x))
      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine ImultiiH_pluslogcomp(kint,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision kint
      double complex S11,P11,S21,P21,S12,S22,icomp,k,k2,icomppi,uncomp
     $     ,zero,sommearg,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp,ctmp1,ctmp2
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex HB0,HB1,HB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
      common/datalog/icomp,icomppi,uncomp,zero
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
!$OMP THREADPRIVATE(/datalog/)

      FNU=0.d0
      k=kmax+icomp*kint

c     premiere espece pour la fonction de Hankel
      zr=kmax*a
      zi=kint*a
      call ZBESH(ZR, ZI, FNU, 1, 1, N, CYR, CYI, NZ, IERR)
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de Hankel1',NZ,IERR
         write(*,*) 'argument',zr,zi
         stop
      endif
      HB0=cyr(1)+icomp*cyi(1)
      HB1=cyr(2)+icomp*cyi(2)
      HB2=cyr(3)+icomp*cyi(3)

      k2=k*k      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=0.d0
      mat_sca_sous12=zero
      mat_sca_sous21=zero
      mat_sca_sous22=0.d0

      mat_prod_sous11=0.d0
      mat_prod_sous12=zero
      mat_prod_sous21=zero
      mat_prod_sous22=0.d0

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=-icomp*(w(i)+w(i+1))*zcouche(i)
c         write(*,*) 'eplus',e_plus,icomp*(w(i)+w(i+1))*zcouche(i),w(i)
c     $        ,w(i+1),zcouche(i)
         e_moins=icomp*(w(i)-w(i+1))*zcouche(i)
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=cdlog((epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $        /ctmp)
         deltap_moins=cdlog((epscouche(i+1)*w(i)-epscouche(i)*w(i+1))
     $        /ctmp)
c         write(*,*) 'delta',deltap_moins,epscouche(i+1)*w(i),epscouche(i
c     $        )*w(i+1),ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=cdlog((w(i+1)+w(i))/ctmp)
         deltas_moins=cdlog((w(i)-w(i+1))/ctmp)
c         write(*,*) 'ddd',cdexp(e_plus),cdexp(e_moins),cdexp(deltap_plus
c     $        ),cdexp(deltap_moins),cdexp(deltas_plus)
c     $        ,cdexp(deltas_moins)
         
         mat_sca_11=e_moins+deltap_plus
         mat_sca_12=e_plus+deltap_moins
         mat_sca_21=deltap_moins-e_plus
         mat_sca_22=deltap_plus-e_moins
c         write(*,*) 'sca',(mat_sca_11),(mat_sca_12)
c     $        ,(mat_sca_21),(mat_sca_22)
         mat_prod_11=e_moins+deltas_plus
         mat_prod_12=icomppi+e_plus+deltas_moins
         mat_prod_21=icomppi+deltas_moins-e_plus
         mat_prod_22=deltas_plus-e_moins
c     produit des deux matrices; scalaire
         ctmp1=mat_sca_11+mat_sca_sous11
         ctmp2=mat_sca_12+mat_sca_sous21
         ctmp=sommearg(ctmp1,ctmp2)
c         write(*,*) 'bbbb',ctmp,ctmp1,ctmp2
         ctmp1=mat_sca_21+mat_sca_sous11
         ctmp2=mat_sca_22+mat_sca_sous21
         mat_sca_sous21=sommearg(ctmp1,ctmp2)
         mat_sca_sous11=ctmp

         ctmp1=mat_sca_11+mat_sca_sous12
         ctmp2=mat_sca_12+mat_sca_sous22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sous12
         ctmp2=mat_sca_22+mat_sca_sous22
         mat_sca_sous22=sommearg(ctmp1,ctmp2)
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp1=mat_prod_11+mat_prod_sous11
         ctmp2=mat_prod_12+mat_prod_sous21
         ctmp=sommearg(ctmp1,ctmp2)
         
         ctmp1=mat_prod_21+mat_prod_sous11
         ctmp2=mat_prod_22+mat_prod_sous21
         mat_prod_sous21=sommearg(ctmp1,ctmp2)
         mat_prod_sous11=ctmp

         ctmp1=mat_prod_11+mat_prod_sous12
         ctmp2=mat_prod_12+mat_prod_sous22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sous12
         ctmp2=mat_prod_22+mat_prod_sous22
         mat_prod_sous22=sommearg(ctmp1,ctmp2)
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',cdexp(mat_sca_sous11),cdexp(mat_sca_sous12
c     $        ),cdexp(mat_sca_sous21),cdexp(mat_sca_sous22)
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
c            write(*,*) 'obs', mat_sca_obs22,cdexp(mat_sca_obs22)
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=0.d0
      mat_sca_sur12=zero
      mat_sca_sur21=zero
      mat_sca_sur22=0.d0

      mat_prod_sur11=0.d0
      mat_prod_sur12=zero
      mat_prod_sur21=zero
      mat_prod_sur22=0.d0

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=-icomp*(w(i)+w(i+1))*zcouche(i)
         e_moins=icomp*(w(i+1)-w(i))*zcouche(i)
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=cdlog((epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $        /ctmp)
         deltap_moins=cdlog((epscouche(i)*w(i+1)-epscouche(i+1)*w(i))
     $        /ctmp)
         ctmp=2.d0*w(i)
         deltas_plus=cdlog((w(i+1)+w(i))/ctmp)
         deltas_moins=cdlog((w(i+1)-w(i))/ctmp)

         mat_sca_11=e_moins+deltap_plus
         mat_sca_12=e_plus+deltap_moins
         mat_sca_21=deltap_moins-e_plus
         mat_sca_22=deltap_plus-e_moins

         mat_prod_11=e_moins+deltas_plus
         mat_prod_12=icomppi+e_plus+deltas_moins
         mat_prod_21=icomppi+deltas_moins-e_plus
         mat_prod_22=deltas_plus-e_moins
c     produit des deux matrices; scalaire
         ctmp1=mat_sca_11+mat_sca_sur11
         ctmp2=mat_sca_12+mat_sca_sur21
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sur11
         ctmp2=mat_sca_22+mat_sca_sur21
         mat_sca_sur21=sommearg(ctmp1,ctmp2)
         mat_sca_sur11=ctmp

         ctmp1=mat_sca_11+mat_sca_sur12
         ctmp2=mat_sca_12+mat_sca_sur22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sur12
         ctmp2=mat_sca_22+mat_sca_sur22
         mat_sca_sur22=sommearg(ctmp1,ctmp2)
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp1=mat_prod_11+mat_prod_sur11
         ctmp2=mat_prod_12+mat_prod_sur21
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sur11
         ctmp2=mat_prod_22+mat_prod_sur21
         mat_prod_sur21=sommearg(ctmp1,ctmp2)
         mat_prod_sur11=ctmp

         ctmp1=mat_prod_11+mat_prod_sur12
         ctmp2=mat_prod_12+mat_prod_sur22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sur12
         ctmp2=mat_prod_22+mat_prod_sur22
         mat_prod_sur22=sommearg(ctmp1,ctmp2)
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=icomppi+mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=icomppi+mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=icomppi+mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=icomppi+mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-icomp*w(nc)*za+icomppi
      nu_moins=icomp*w(nc)*za+icomppi


      ctmp1=cdlog(w(nc))
      ctmp2=cdlog(k2)
      B_sca_11=icomppi+nu_plus+ctmp1
      B_sca_12=nu_plus+ctmp2
      B_sca_21=nu_moins+ctmp1
      B_sca_22=nu_moins+ctmp2
      ctmp=cdlog(epscouche(nc)*k02/w(nc))
c      write(*,*) 'B_sca_22',B_sca_22,cdexp(B_sca_22)
      B_prod_11=icomppi+nu_plus+ctmp
      B_prod_21=nu_moins+ctmp

c     inverse de la matrice coupee et produit par B
      ctmp1=mat_sca_11+mat_sca_22
      ctmp2=mat_sca_12+mat_sca_21+icomppi
      det=sommearg(ctmp1,ctmp2)
c      write(*,*) 'det',(det)
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22-det
      mat_sca_22=ctmp-det
      mat_sca_12=icomppi+mat_sca_12-det
      mat_sca_21=icomppi+mat_sca_21-det

      ctmp1=mat_sca_11+B_sca_11
      ctmp2=mat_sca_12+B_sca_21
      ctmp=sommearg(ctmp1,ctmp2)
      ctmp1=mat_sca_11+B_sca_12
      ctmp2=mat_sca_12+B_sca_22
      mat_sca_12=sommearg(ctmp1,ctmp2)
      mat_sca_11=ctmp

      ctmp1=mat_sca_21+B_sca_11
      ctmp2=mat_sca_22+B_sca_21
      ctmp=sommearg(ctmp1,ctmp2)
c      write(*,*) 'aaa',ctmp,ctmp1,ctmp2
      ctmp1=mat_sca_21+B_sca_12
      ctmp2=mat_sca_22+B_sca_22
      mat_sca_22=sommearg(ctmp1,ctmp2)
      mat_sca_21=ctmp


      ctmp1=mat_prod_11+mat_prod_22
      ctmp2=icomppi+mat_prod_12+mat_prod_21
      det=sommearg(ctmp1,ctmp2)
c      write(*,*) 'det',(det)
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22-det
      mat_prod_22=ctmp-det
      mat_prod_12=icomppi+mat_prod_12-det
      mat_prod_21=icomppi+mat_prod_21-det

      ctmp1=mat_prod_21+B_prod_11
      ctmp2=mat_prod_22+B_prod_21
      mat_prod_21=sommearg(ctmp1,ctmp2)
      ctmp1=mat_prod_11+B_prod_11
      ctmp2=mat_prod_12+B_prod_21
      mat_prod_11=sommearg(ctmp1,ctmp2)

      ctmp=icomp*w(no)*z

      S11=mat_sca_obs12+mat_sca_21+ctmp
      P11=mat_prod_obs12+mat_prod_21+ctmp
      S12=mat_sca_obs12+mat_sca_22+ctmp
      ctmp=-icomp*w(no)*z
      ctmp1=mat_sca_obs22+mat_sca_21+ctmp
      ctmp2=B_sca_21+ctmp
      S21=sommearg(ctmp1,ctmp2)
      ctmp1=mat_sca_obs22+mat_sca_22+ctmp
c     write(*,*) 'rr',cdexp(mat_sca_obs22),cdexp(mat_sca_22)
      ctmp2=B_sca_22+ctmp
      S22=sommearg(ctmp1,ctmp2)
c     write(*,*) 'S22',cdexp(S22),cdexp(ctmp1),cdexp(ctmp2)
c     $           ,cdexp(ctmp)
      ctmp1=mat_prod_obs22+mat_prod_21+ctmp
      ctmp2=B_prod_21+ctmp
      P21=sommearg(ctmp1,ctmp2)

c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22
c      write(*,*) 'SS',S11,P11,S21,P21,S12,S22
      S11=cdexp(S11)
      S12=cdexp(S12)
      S21=cdexp(S21)
      S22=cdexp(S22)
      P11=cdexp(P11)
      P21=cdexp(P21)

 
      integrant(1)=k*(S11+P11+S21+P21)*HB0
      integrant(2)=k*(S11-P11+S21-P21)*HB2
      integrant(3)=(S12+S22)*HB1
      integrant(4)=k2/w(no)*(-S11+S21)*HB1
      integrant(5)=(-S12+S22)*k/w(no)*HB0

c      write(230,*) kint,integrant(1)
c      write(231,*) kint,integrant(2)
c      write(232,*) kint,integrant(3)
c      write(233,*) kint,integrant(4)

      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
      subroutine ImultiiH_moinslogcomp(kint,nnn,nlda,integrant)
      implicit none
      integer i,nnn,nlda
      double precision kint
      double complex S11,P11,S21,P21,S12,S22,icomp,k,k2,icomppi,uncomp
     $     ,zero,sommearg,integrant(nlda)
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp,ctmp1,ctmp2
      double complex mat_sca_sous11,mat_sca_sous12,mat_sca_sous21
     $     ,mat_sca_sous22
      double complex mat_sca_11,mat_sca_12,mat_sca_21,mat_sca_22
      double complex mat_prod_sous11,mat_prod_sous12,mat_prod_sous21
     $     ,mat_prod_sous22
      double complex mat_prod_11,mat_prod_12,mat_prod_21,mat_prod_22
      double complex mat_sca_sur11,mat_sca_sur12,mat_sca_sur21
     $     ,mat_sca_sur22
      double complex mat_prod_sur11,mat_prod_sur12,mat_prod_sur21
     $     ,mat_prod_sur22
      double complex B_sca_11,B_sca_12,B_sca_21,B_sca_22
      double complex B_prod_11,B_prod_21

      double complex mat_sca_obs11,mat_sca_obs12,mat_sca_obs21
     $     ,mat_sca_obs22
      double complex mat_prod_obs11,mat_prod_obs12,mat_prod_obs21
     $     ,mat_prod_obs22

      double complex nu_plus,nu_moins,det

      integer nepsmax,neps,nc,no
      parameter(nepsmax=20)

      double precision k02,kmax,hkmax,a,z,za,zcouche(0:nepsmax)
      double complex epscouche(0:nepsmax+1),w(0:nepsmax+1)

      integer N,IERR,NZ
      parameter (N=3)
      DOUBLE PRECISION ZR, ZI, CYR(n), CYI(n), FNU
      double complex HB0,HB1,HB2

      common/donneemultin/neps,nc,no      
      common/donneemulti/k02,kmax,hkmax,a,z,za
      common/zcoucheroutine/zcouche
      common/epscouroutine/epscouche
      common/datalog/icomp,icomppi,uncomp,zero
!$OMP THREADPRIVATE(/donneemultin/)
!$OMP THREADPRIVATE(/donneemulti/)
!$OMP THREADPRIVATE(/zcoucheroutine/)
!$OMP THREADPRIVATE(/epscouroutine/)
!$OMP THREADPRIVATE(/datalog/)

      FNU=0.d0
      k=kmax-icomp*kint

c     deuxieme espece pour la fonction de Hankel
      zr=kmax*a
      zi=-kint*a
      call ZBESH(ZR, ZI, FNU, 1, 2, N, CYR, CYI, NZ, IERR)
      if (IERR.ne.0.or.NZ.ne.0) then
         write(*,*) 'probleme dans la fonction de Hankel2',NZ,IERR
         write(*,*) 'argument',zr,zi
         stop
      endif
      HB0=cyr(1)+icomp*cyi(1)
      HB1=cyr(2)+icomp*cyi(2)
      HB2=cyr(3)+icomp*cyi(3)

      k2=k*k      
c     calcul des wz pour toutes les couches
      do i=0,neps+1
         w(i)=cdsqrt(epscouche(i)*k02-k2)
         if (dimag(w(i)).lt.0.d0) w(i)=-w(i)      
c         write(*,*) 'www',w(i),k,kmax
c         write(*,*) 'eps',epscouche(i),k02,nepsmax,neps
      enddo

c     initialise pour l'interface 0    
      mat_sca_sous11=0.d0
      mat_sca_sous12=zero
      mat_sca_sous21=zero
      mat_sca_sous22=0.d0

      mat_prod_sous11=0.d0
      mat_prod_sous12=zero
      mat_prod_sous21=zero
      mat_prod_sous22=0.d0

c     calcul matrice de l'interface 0 a nc-1
      do i=0,nc-1
         e_plus=-icomp*(w(i)+w(i+1))*zcouche(i)
c         write(*,*) 'eplus',e_plus,icomp*(w(i)+w(i+1))*zcouche(i),w(i)
c     $        ,w(i+1),zcouche(i)
         e_moins=icomp*(w(i)-w(i+1))*zcouche(i)
         ctmp=2.d0*epscouche(i+1)*w(i)
         deltap_plus=cdlog((epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $        /ctmp)
         deltap_moins=cdlog((epscouche(i+1)*w(i)-epscouche(i)*w(i+1))
     $        /ctmp)
c         write(*,*) 'delta',deltap_moins,epscouche(i+1)*w(i),epscouche(i
c     $        )*w(i+1),ctmp
         ctmp=2.d0*w(i+1)
         deltas_plus=cdlog((w(i+1)+w(i))/ctmp)
         deltas_moins=cdlog((w(i)-w(i+1))/ctmp)
c         write(*,*) 'ddd',cdexp(e_plus),cdexp(e_moins),cdexp(deltap_plus
c     $        ),cdexp(deltap_moins),cdexp(deltas_plus)
c     $        ,cdexp(deltas_moins)
         
         mat_sca_11=e_moins+deltap_plus
         mat_sca_12=e_plus+deltap_moins
         mat_sca_21=deltap_moins-e_plus
         mat_sca_22=deltap_plus-e_moins
c         write(*,*) 'sca',(mat_sca_11),(mat_sca_12)
c     $        ,(mat_sca_21),(mat_sca_22)
         mat_prod_11=e_moins+deltas_plus
         mat_prod_12=icomppi+e_plus+deltas_moins
         mat_prod_21=icomppi+deltas_moins-e_plus
         mat_prod_22=deltas_plus-e_moins
c     produit des deux matrices; scalaire
         ctmp1=mat_sca_11+mat_sca_sous11
         ctmp2=mat_sca_12+mat_sca_sous21
         ctmp=sommearg(ctmp1,ctmp2)
c         write(*,*) 'bbbb',ctmp,ctmp1,ctmp2
         ctmp1=mat_sca_21+mat_sca_sous11
         ctmp2=mat_sca_22+mat_sca_sous21
         mat_sca_sous21=sommearg(ctmp1,ctmp2)
         mat_sca_sous11=ctmp

         ctmp1=mat_sca_11+mat_sca_sous12
         ctmp2=mat_sca_12+mat_sca_sous22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sous12
         ctmp2=mat_sca_22+mat_sca_sous22
         mat_sca_sous22=sommearg(ctmp1,ctmp2)
         mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
         ctmp1=mat_prod_11+mat_prod_sous11
         ctmp2=mat_prod_12+mat_prod_sous21
         ctmp=sommearg(ctmp1,ctmp2)
         
         ctmp1=mat_prod_21+mat_prod_sous11
         ctmp2=mat_prod_22+mat_prod_sous21
         mat_prod_sous21=sommearg(ctmp1,ctmp2)
         mat_prod_sous11=ctmp

         ctmp1=mat_prod_11+mat_prod_sous12
         ctmp2=mat_prod_12+mat_prod_sous22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sous12
         ctmp2=mat_prod_22+mat_prod_sous22
         mat_prod_sous22=sommearg(ctmp1,ctmp2)
         mat_prod_sous12=ctmp

c         write(*,*) 'soussca',cdexp(mat_sca_sous11),cdexp(mat_sca_sous12
c     $        ),cdexp(mat_sca_sous21),cdexp(mat_sca_sous22)
c         write(*,*) 'sousprod',mat_prod_sous11,mat_prod_sous12
c     $        ,mat_prod_sous21,mat_prod_sous22
c         stop
         if (i+1.eq.no) then
            mat_sca_obs11=mat_sca_sous11
            mat_sca_obs12=mat_sca_sous12
            mat_sca_obs21=mat_sca_sous21
            mat_sca_obs22=mat_sca_sous22
c            write(*,*) 'obs', mat_sca_obs22,cdexp(mat_sca_obs22)
            mat_prod_obs11=mat_prod_sous11
            mat_prod_obs12=mat_prod_sous12
            mat_prod_obs21=mat_prod_sous21
            mat_prod_obs22=mat_prod_sous22
         endif
         
      enddo


c     initialise pour interface dessus
      mat_sca_sur11=0.d0
      mat_sca_sur12=zero
      mat_sca_sur21=zero
      mat_sca_sur22=0.d0

      mat_prod_sur11=0.d0
      mat_prod_sur12=zero
      mat_prod_sur21=zero
      mat_prod_sur22=0.d0

c     calcul matrice de l'interface neps-1 a nc
      do i=neps,nc,-1
         e_plus=-icomp*(w(i)+w(i+1))*zcouche(i)
         e_moins=icomp*(w(i+1)-w(i))*zcouche(i)
         ctmp=2.d0*epscouche(i)*w(i+1)
         deltap_plus=cdlog((epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $        /ctmp)
         deltap_moins=cdlog((epscouche(i)*w(i+1)-epscouche(i+1)*w(i))
     $        /ctmp)
         ctmp=2.d0*w(i)
         deltas_plus=cdlog((w(i+1)+w(i))/ctmp)
         deltas_moins=cdlog((w(i+1)-w(i))/ctmp)

         mat_sca_11=e_moins+deltap_plus
         mat_sca_12=e_plus+deltap_moins
         mat_sca_21=deltap_moins-e_plus
         mat_sca_22=deltap_plus-e_moins

         mat_prod_11=e_moins+deltas_plus
         mat_prod_12=icomppi+e_plus+deltas_moins
         mat_prod_21=icomppi+deltas_moins-e_plus
         mat_prod_22=deltas_plus-e_moins
c     produit des deux matrices; scalaire
         ctmp1=mat_sca_11+mat_sca_sur11
         ctmp2=mat_sca_12+mat_sca_sur21
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sur11
         ctmp2=mat_sca_22+mat_sca_sur21
         mat_sca_sur21=sommearg(ctmp1,ctmp2)
         mat_sca_sur11=ctmp

         ctmp1=mat_sca_11+mat_sca_sur12
         ctmp2=mat_sca_12+mat_sca_sur22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_sca_21+mat_sca_sur12
         ctmp2=mat_sca_22+mat_sca_sur22
         mat_sca_sur22=sommearg(ctmp1,ctmp2)
         mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
         ctmp1=mat_prod_11+mat_prod_sur11
         ctmp2=mat_prod_12+mat_prod_sur21
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sur11
         ctmp2=mat_prod_22+mat_prod_sur21
         mat_prod_sur21=sommearg(ctmp1,ctmp2)
         mat_prod_sur11=ctmp

         ctmp1=mat_prod_11+mat_prod_sur12
         ctmp2=mat_prod_12+mat_prod_sur22
         ctmp=sommearg(ctmp1,ctmp2)
         ctmp1=mat_prod_21+mat_prod_sur12
         ctmp2=mat_prod_22+mat_prod_sur22
         mat_prod_sur22=sommearg(ctmp1,ctmp2)
         mat_prod_sur12=ctmp

c         write(*,*) 'sursca',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c         write(*,*) 'surprod',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22

      enddo

c     calcul de la matrice coupee
      mat_sca_11=mat_sca_sur11
      mat_sca_12=icomppi+mat_sca_sous12
      mat_sca_21=mat_sca_sur21
      mat_sca_22=icomppi+mat_sca_sous22

      mat_prod_11=mat_prod_sur11
      mat_prod_12=icomppi+mat_prod_sous12
      mat_prod_21=mat_prod_sur21
      mat_prod_22=icomppi+mat_prod_sous22


c     calcul de la matrice B
      nu_plus=-icomp*w(nc)*za+icomppi
      nu_moins=icomp*w(nc)*za+icomppi


      ctmp1=cdlog(w(nc))
      ctmp2=cdlog(k2)
      B_sca_11=icomppi+nu_plus+ctmp1
      B_sca_12=nu_plus+ctmp2
      B_sca_21=nu_moins+ctmp1
      B_sca_22=nu_moins+ctmp2
      ctmp=cdlog(epscouche(nc)*k02/w(nc))
c      write(*,*) 'B_sca_22',B_sca_22,cdexp(B_sca_22)
      B_prod_11=icomppi+nu_plus+ctmp
      B_prod_21=nu_moins+ctmp

c     inverse de la matrice coupee et produit par B
      ctmp1=mat_sca_11+mat_sca_22
      ctmp2=mat_sca_12+mat_sca_21+icomppi
      det=sommearg(ctmp1,ctmp2)
c      write(*,*) 'det',(det)
      ctmp=mat_sca_11
      mat_sca_11=mat_sca_22-det
      mat_sca_22=ctmp-det
      mat_sca_12=icomppi+mat_sca_12-det
      mat_sca_21=icomppi+mat_sca_21-det

      ctmp1=mat_sca_11+B_sca_11
      ctmp2=mat_sca_12+B_sca_21
      ctmp=sommearg(ctmp1,ctmp2)
      ctmp1=mat_sca_11+B_sca_12
      ctmp2=mat_sca_12+B_sca_22
      mat_sca_12=sommearg(ctmp1,ctmp2)
      mat_sca_11=ctmp

      ctmp1=mat_sca_21+B_sca_11
      ctmp2=mat_sca_22+B_sca_21
      ctmp=sommearg(ctmp1,ctmp2)
c      write(*,*) 'aaa',ctmp,ctmp1,ctmp2
      ctmp1=mat_sca_21+B_sca_12
      ctmp2=mat_sca_22+B_sca_22
      mat_sca_22=sommearg(ctmp1,ctmp2)
      mat_sca_21=ctmp


      ctmp1=mat_prod_11+mat_prod_22
      ctmp2=icomppi+mat_prod_12+mat_prod_21
      det=sommearg(ctmp1,ctmp2)
c      write(*,*) 'det',(det)
      ctmp=mat_prod_11
      mat_prod_11=mat_prod_22-det
      mat_prod_22=ctmp-det
      mat_prod_12=icomppi+mat_prod_12-det
      mat_prod_21=icomppi+mat_prod_21-det

      ctmp1=mat_prod_21+B_prod_11
      ctmp2=mat_prod_22+B_prod_21
      mat_prod_21=sommearg(ctmp1,ctmp2)
      ctmp1=mat_prod_11+B_prod_11
      ctmp2=mat_prod_12+B_prod_21
      mat_prod_11=sommearg(ctmp1,ctmp2)

      ctmp=icomp*w(no)*z

      S11=mat_sca_obs12+mat_sca_21+ctmp
      P11=mat_prod_obs12+mat_prod_21+ctmp
      S12=mat_sca_obs12+mat_sca_22+ctmp
      ctmp=-icomp*w(no)*z
      ctmp1=mat_sca_obs22+mat_sca_21+ctmp
      ctmp2=B_sca_21+ctmp
      S21=sommearg(ctmp1,ctmp2)
      ctmp1=mat_sca_obs22+mat_sca_22+ctmp
c     write(*,*) 'rr',cdexp(mat_sca_obs22),cdexp(mat_sca_22)
      ctmp2=B_sca_22+ctmp
      S22=sommearg(ctmp1,ctmp2)
c     write(*,*) 'S22',cdexp(S22),cdexp(ctmp1),cdexp(ctmp2)
c     $           ,cdexp(ctmp)
      ctmp1=mat_prod_obs22+mat_prod_21+ctmp
      ctmp2=B_prod_21+ctmp
      P21=sommearg(ctmp1,ctmp2)

c     calcul element differentiel
c     calcul des integrants partie reelle et virtuelle
c      write(*,*) 'fin',S11,P11,S21,P21,S12,S22
c      write(*,*) 'SS',S11,P11,S21,P21,S12,S22
      S11=cdexp(S11)
      S12=cdexp(S12)
      S21=cdexp(S21)
      S22=cdexp(S22)
      P11=cdexp(P11)
      P21=cdexp(P21)

 
      integrant(1)=k*(S11+P11+S21+P21)*HB0
      integrant(2)=k*(S11-P11+S21-P21)*HB2
      integrant(3)=(S12+S22)*HB1
      integrant(4)=k2/w(no)*(-S11+S21)*HB1
      integrant(5)=(-S12+S22)*k/w(no)*HB0

c      write(40,*) kint,integrant(1)
c      write(41,*) kint,integrant(2)
c      write(42,*) kint,integrant(3)
c      write(43,*) kint,integrant(4)

      end
