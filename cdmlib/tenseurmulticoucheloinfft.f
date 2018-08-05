      subroutine tenseurmulticoucheloinfft(kx,ky,kz,z,za,k0,nepsmax,neps
     $     ,dcouche,zcouche,epscouche,Stenseur)
      implicit none
      integer nepsmax,neps,nc,no
      integer i,j
      double precision dcouche(nepsmax),zcouche(0:nepsmax),r,pi,k02
      double complex epscouche(0:nepsmax+1)
      double complex Stenseur(3,3),w(0:nepsmax+1),k2(0:nepsmax+1)
      double precision z,za,k0,kd,kd2,kx,ky,kz,kp2
      
      double complex e_moins,e_plus,deltap_plus,deltap_moins,deltas_plus
     $     ,deltas_moins,ctmp,ctmp1
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

      double complex icomp,nu_plus,nu_moins,det,S11,P11,S21,P21,S12,S22

      
      icomp=(0.d0,1.d0)
      pi=dacos(-1.d0)
      k02=k0*k0    

c     repere dans quelle couche est le dipole
      nc=0
      do i=0,neps
         if (za.ge.zcouche(i)) then
            nc=nc+1
         else
            goto 10
         endif
      enddo


c     repere dans quelle couche est l'observation
 10   if (z.gt.zcouche(neps)) then
         no=neps+1
      elseif (z.lt.zcouche(0)) then
         no=0
      else
         write(*,*) 'problem',z,zcouche(neps),zcouche(0)
         stop
      endif
         
c     write(*,*) 'numero couche observation',no,z,neps


c     *************** champ au desssus *********************
      if (no.eq.neps+1) then

         if (dimag(epscouche(neps+1)).ne.0) then
            write(*,*) 'milieu du champ lointain absorbant'
            stop
         endif

c     calcul module pour le champ diffracte         
         kd=k0*dsqrt(dreal(epscouche(neps+1)))
         kd2=kd*kd
         kp2=kx*kx+ky*ky

         do i=0,neps+1
            k2(i)=epscouche(i)*k0*k0
            w(i)=cdsqrt(k2(i)-kp2)
            if (dimag(w(i)).lt.0.d0) w(i)=-w(i)
            if (cdabs(w(i)).eq.0.d0) w(i)=k0/1.d10
c     write(*,*) 'w',w(i),i,kz
         enddo
         if (nc.eq.0) then
c     write(*,*) 'dipole dessous'
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
               deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $              /ctmp
               deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))
     $              /ctmp
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
     $              *mat_sca_sur21
               mat_sca_sur11=ctmp
               ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
               mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $              *mat_sca_sur22
               mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
               ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12
     $              *mat_prod_sur21
               mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $              *mat_prod_sur21
               mat_prod_sur11=ctmp
               ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12
     $              *mat_prod_sur22
               mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $              *mat_prod_sur22
               mat_prod_sur12=ctmp

c     write(*,*) 'sursca2',mat_sca_sur11,mat_sca_sur12
c     $        ,mat_sca_sur21,mat_sca_sur22
c     write(*,*) 'surprod2',mat_prod_sur11,mat_prod_sur12
c     $        ,mat_prod_sur21,mat_prod_sur22
c     write(*,*) 'i2',i,no,nc,neps


            enddo

            mat_sca_11=mat_sca_sur11
            mat_prod_11=mat_prod_sur11

c     calcul de la matrice B
            nu_plus=-cdexp(-icomp*w(nc)*za)
            B_sca_11=-nu_plus*w(nc)
            B_sca_12=nu_plus*kp2
            B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
            mat_sca_12=B_sca_12/mat_sca_11
            mat_sca_11=B_sca_11/mat_sca_11
            mat_prod_11=B_prod_11/mat_prod_11

            S11=mat_sca_11
            S12=mat_sca_12
            P11=mat_prod_11

            if (kp2.eq.0.d0) then
               Stenseur(1,1)=S11*w(neps+1)
               Stenseur(1,2)=0.d0
               Stenseur(1,3)=0.d0
               Stenseur(2,1)=0.d0
               Stenseur(2,2)=S11*w(neps+1)
               Stenseur(2,3)=0.d0
               Stenseur(3,1)=0.d0
               Stenseur(3,2)=0.d0
               Stenseur(3,3)=-S12
            else
               Stenseur(1,1)=(kx*kx*S11+ky*ky*P11)/kp2*w(neps+1)
               Stenseur(1,2)=kx*ky*(S11-P11)/kp2*w(neps+1)
               Stenseur(1,3)=kx*S12/kp2*w(neps+1)
               Stenseur(2,1)=Stenseur(1,2)
               Stenseur(2,2)=(ky*ky*S11+kx*kx*P11)/kp2*w(neps+1)
               Stenseur(2,3)=ky*S12/kp2*w(neps+1)
               Stenseur(3,1)=-kx*S11
               Stenseur(3,2)=-ky*S11
               Stenseur(3,3)=-S12
            endif

         elseif (nc.eq.neps+1) then
c     write(*,*) 'dipole dessus'
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
               deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $              /ctmp
               deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))
     $              /ctmp
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
     $              *mat_sca_sous21
               mat_sca_sous11=ctmp
               ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
               mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $              *mat_sca_sous22
               mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
               ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12
     $              *mat_prod_sous21
               mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $              *mat_prod_sous21
               mat_prod_sous11=ctmp
               ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12
     $              *mat_prod_sous22
               mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $              *mat_prod_sous22
               mat_prod_sous12=ctmp

            enddo


            mat_sca_22=-mat_sca_sous22
            mat_prod_22=-mat_prod_sous22
c     calcul de la matrice B
            nu_moins=-cdexp(icomp*w(nc)*za)
            B_sca_21=nu_moins*w(nc)
            B_sca_22=nu_moins*kp2
            B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)
            
            mat_sca_21=B_sca_21/mat_sca_22
            mat_sca_22=B_sca_22/mat_sca_22
            mat_prod_21=B_prod_21/mat_prod_22
            
            
            S11=mat_sca_sous12*mat_sca_21
            P11=mat_prod_sous12*mat_prod_21
            S12=mat_sca_sous12*mat_sca_22

c     write(*,*) 'aa kp2',a,kp2,kx,ky,kz,r
c     calcul le tenseur
            
            if (kp2.eq.0.d0) then
               ctmp1=cdexp(icomp*(-kz*za))*kd2
               Stenseur(1,1)=S11*w(neps+1)+ctmp1
               Stenseur(1,2)=0.d0
               Stenseur(1,3)=0.d0
               Stenseur(2,1)=0.d0
               Stenseur(2,2)=S11*w(neps+1)+ctmp1
               Stenseur(2,3)=0.d0
               Stenseur(3,1)=0.d0
               Stenseur(3,2)=0.d0
               Stenseur(3,3)=-S12
            else
               ctmp1=cdexp(icomp*(-kz*za))*kd2
               Stenseur(1,1)=(kx*kx*S11+ky*ky*P11)/kp2*w(neps+1)+ctmp1
     $              *(1.d0-kx*kx/kd2)
               Stenseur(1,2)=kx*ky*(S11-P11)/kp2*w(neps+1)+ctmp1*(-kx*ky
     $              /kd2)
               Stenseur(1,3)=kx*S12/kp2*w(neps+1)+ctmp1*(-kx*kz/kd2)
               Stenseur(2,1)=Stenseur(1,2)
               Stenseur(2,2)=(ky*ky*S11+kx*kx*P11)/kp2*w(neps+1)+ctmp1
     $              *(1.d0-ky*ky/kd2)
               Stenseur(2,3)=ky*S12/kp2*w(neps+1)+ctmp1*(-ky*kz/kd2)
               Stenseur(3,1)=-kx*S11+ctmp1*(-kz*kx/kd2)
               Stenseur(3,2)=-ky*S11+ctmp1*(-kz*ky/kd2)
               Stenseur(3,3)=-S12+ctmp1*(1.d0-kz*kz/kd2)
            endif

         else
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
               deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $              /ctmp
               deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))
     $              /ctmp
               ctmp=2.d0*w(i+1)
               deltas_plus=(w(i+1)+w(i))/ctmp
               deltas_moins=(w(i)-w(i+1))/ctmp
c     write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $           ,deltas_plus,deltas_moins,ctmp
               
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
     $              *mat_sca_sous21
               mat_sca_sous11=ctmp
               ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
               mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $              *mat_sca_sous22
               mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
               ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12
     $              *mat_prod_sous21
               mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $              *mat_prod_sous21
               mat_prod_sous11=ctmp
               ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12
     $              *mat_prod_sous22
               mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $              *mat_prod_sous22
               mat_prod_sous12=ctmp
               

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
               deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $              /ctmp
               deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))
     $              /ctmp
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
     $              *mat_sca_sur21
               mat_sca_sur11=ctmp
               ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
               mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $              *mat_sca_sur22
               mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
               ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12
     $              *mat_prod_sur21
               mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $              *mat_prod_sur21
               mat_prod_sur11=ctmp
               ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12
     $              *mat_prod_sur22
               mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $              *mat_prod_sur22
               mat_prod_sur12=ctmp
               
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
            B_sca_12=nu_plus*kp2
            B_sca_21=nu_moins*w(nc)
            B_sca_22=nu_moins*kp2
            
            B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
            B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)
            
c     write(*,*) 'Bprod',B_prod_11,B_prod_21,k02

c     inverse de la matrice coupee et produit par B
            det=mat_sca_11*mat_sca_22-mat_sca_12*mat_sca_21
c     write(*,*) 'det',det,mat_sca_11,mat_sca_22,mat_sca_12,mat_sca_21
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
c     write(*,*) 'det',det
            ctmp=mat_prod_11
            mat_prod_11=mat_prod_22/det
            mat_prod_22=ctmp/det
            mat_prod_12=-mat_prod_12/det
            mat_prod_21=-mat_prod_21/det
            
            mat_prod_21=mat_prod_21*B_prod_11+mat_prod_22*B_prod_21
            mat_prod_11=mat_prod_11*B_prod_11+mat_prod_12*B_prod_21
c     write(*,*) 'B11',mat_prod_11,B_prod_11,mat_prod_12,B_prod_21
            mat_prod_22=0.d0
            mat_prod_12=0.d0
            
c     le point d'observation est dessus
            S11=mat_sca_11
            S12=mat_sca_12
            P11=mat_prod_11

c     calcul le tenseur
            if (kp2.eq.0.d0) then
               Stenseur(1,1)=S11*w(neps+1)
               Stenseur(1,2)=0.d0
               Stenseur(1,3)=0.d0
               Stenseur(2,1)=0.d0
               Stenseur(2,2)=S11*w(neps+1)
               Stenseur(2,3)=0.d0
               Stenseur(3,1)=0.d0
               Stenseur(3,2)=0.d0
               Stenseur(3,3)=-S12
            else
               Stenseur(1,1)=(kx*kx*S11+ky*ky*P11)/kp2*w(neps+1)
               Stenseur(1,2)=kx*ky*(S11-P11)/kp2*w(neps+1)
               Stenseur(1,3)=kx*S12/kp2*w(neps+1)
               Stenseur(2,1)=Stenseur(1,2)
               Stenseur(2,2)=(ky*ky*S11+kx*kx*P11)/kp2*w(neps+1)
               Stenseur(2,3)=ky*S12/kp2*w(neps+1)
               Stenseur(3,1)=-kx*S11
               Stenseur(3,2)=-ky*S11
               Stenseur(3,3)=-S12
            endif

         endif

c************************************************************
c     ATTAQUE CAS OU ON REGARDE DESSOUS
c************************************************************
      elseif (no.eq.0) then
c     write(*,*) 'dessous lointain'
         if (dimag(epscouche(0)).ne.0) then
            write(*,*) 'milieu du champ lointain absorbant'
            stop
         endif
         
c     calcul module pour le champ diffracte
         kd=k0*dsqrt(dreal(epscouche(0)))
         kd2=kd*kd
         kp2=kx*kx+ky*ky
c*************************************************************
c     calcul des k2 et wz pour toutes les couches
         do i=0,neps+1
            k2(i)=epscouche(i)*k0*k0
            w(i)=cdsqrt(k2(i)-kp2)
            if (dimag(w(i)).lt.0.d0) w(i)=-w(i)
            if (cdabs(w(i)).eq.0.d0) w(i)=k0/1.d10
c     write(*,*) 'w',w(i),i,kz
         enddo

         if (nc.eq.0) then
c            write(*,*) 'dipole dessous'
            
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
               deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $              /ctmp
               deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))
     $              /ctmp
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
     $              *mat_sca_sur21
               mat_sca_sur11=ctmp
               ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
               mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $              *mat_sca_sur22
               mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
               ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12
     $              *mat_prod_sur21
               mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $              *mat_prod_sur21
               mat_prod_sur11=ctmp
               ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12
     $              *mat_prod_sur22
               mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $              *mat_prod_sur22
               mat_prod_sur12=ctmp

            enddo

            mat_sca_11=mat_sca_sur11          
            mat_prod_11=mat_prod_sur11
            
c     calcul de la matrice B
            nu_plus=-cdexp(-icomp*w(nc)*za)
            B_sca_11=-nu_plus*w(nc)
            B_sca_12=nu_plus*kp2
            B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)

c     inverse de la matrice coupee et produit par B
            mat_sca_12=B_sca_12/mat_sca_11
            mat_sca_11=B_sca_11/mat_sca_11
            mat_prod_11=B_prod_11/mat_prod_11


            S21=mat_sca_sur21*mat_sca_11
            S22=mat_sca_sur21*mat_sca_12
            P21=mat_prod_sur21*mat_prod_11
            if (kp2.eq.0.d0) then
               ctmp1=cdexp(icomp*(-kz*za))*kd2
               Stenseur(1,1)=S21*w(0)+ctmp1
               Stenseur(1,2)=0.d0
               Stenseur(1,3)=0.d0
               Stenseur(2,1)=0.d0
               Stenseur(2,2)=S21*w(0)+ctmp1
               Stenseur(2,3)=0.d0
               Stenseur(3,1)=0.d0
               Stenseur(3,2)=0.d0
               Stenseur(3,3)=S22
            else
               ctmp1=cdexp(icomp*(-kz*za))*kd2
               Stenseur(1,1)=(kx*kx*S21+ky*ky*P21)/kp2*w(0)+ctmp1*(1.d0
     $              -kx*kx/kd2)
               Stenseur(1,2)=kx*ky*(S21-P21)/kp2*w(0)+ctmp1*(-kx*ky/kd2)
               Stenseur(1,3)=kx*S22/kp2*w(0)+ctmp1*(-kx*kz/kd2)
               Stenseur(2,1)=Stenseur(1,2)
               Stenseur(2,2)=(ky*ky*S21+kx*kx*P21)/kp2*w(0)+ctmp1*(1.d0
     $              -ky*ky/kd2)
               Stenseur(2,3)=ky*S22/kp2*w(0)+ctmp1*(-ky*kz/kd2)
               Stenseur(3,1)=kx*S21+ctmp1*(-kz*kx/kd2)
               Stenseur(3,2)=ky*S21+ctmp1*(-kz*ky/kd2)
               Stenseur(3,3)=S22+ctmp1*(1.d0-kz*kz/kd2)
            endif

         elseif (nc.eq.neps+1) then
c            write(*,*) 'dipole dessus'
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
               deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $              /ctmp
               deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))
     $              /ctmp
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
     $              *mat_sca_sous21
               mat_sca_sous11=ctmp
               ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
               mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $              *mat_sca_sous22
               mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
               ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12
     $              *mat_prod_sous21
               mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $              *mat_prod_sous21
               mat_prod_sous11=ctmp
               ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12
     $              *mat_prod_sous22
               mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $              *mat_prod_sous22
               mat_prod_sous12=ctmp

            enddo

            mat_sca_12=-mat_sca_sous12
            mat_sca_22=-mat_sca_sous22
            mat_prod_22=-mat_prod_sous22

c     calcul de la matrice B
            
            nu_moins=-cdexp(icomp*w(nc)*za)
            B_sca_21=nu_moins*w(nc)
            B_sca_22=nu_moins*kp2
            B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)

            mat_sca_21=B_sca_21/mat_sca_22
            mat_sca_22=B_sca_22/mat_sca_22
            mat_prod_21=B_prod_21/mat_prod_22
            
            S21=mat_sca_21
            S22=mat_sca_22
            P21=mat_prod_21
            if (kp2.eq.0.d0) then
               Stenseur(1,1)=S21*w(0)
               Stenseur(1,2)=0.d0
               Stenseur(1,3)=0.d0
               Stenseur(2,1)=0.d0
               Stenseur(2,2)=S21*w(0)
               Stenseur(2,3)=0.d0
               Stenseur(3,1)=0.d0
               Stenseur(3,2)=0.d0
               Stenseur(3,3)=S22
            else
               Stenseur(1,1)=(kx*kx*S21+ky*ky*P21)/kp2*w(0)
               Stenseur(1,2)=kx*ky*(S21-P21)/kp2*w(0)
               Stenseur(1,3)=kx*S22/kp2*w(0)
               Stenseur(2,1)=Stenseur(1,2)
               Stenseur(2,2)=(ky*ky*S21+kx*kx*P21)/kp2*w(0)
               Stenseur(2,3)=ky*S22/kp2*w(0)
               Stenseur(3,1)=kx*S21
               Stenseur(3,2)=ky*S21
               Stenseur(3,3)=S22
            endif
         else
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
               deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $              /ctmp
               deltap_moins=(epscouche(i+1)*w(i)-epscouche(i)*w(i+1))
     $              /ctmp
               ctmp=2.d0*w(i+1)
               deltas_plus=(w(i+1)+w(i))/ctmp
               deltas_moins=(w(i)-w(i+1))/ctmp
c     write(*,*) 'ddd',e_plus,e_moins,deltap_plus,deltap_moins
c     $           ,deltas_plus,deltas_moins,ctmp
               
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
     $              *mat_sca_sous21
               mat_sca_sous11=ctmp
               ctmp=mat_sca_11*mat_sca_sous12+mat_sca_12*mat_sca_sous22
               mat_sca_sous22=mat_sca_21*mat_sca_sous12+mat_sca_22
     $              *mat_sca_sous22
               mat_sca_sous12=ctmp
c     produit des deux matrices; vectorielle
               ctmp=mat_prod_11*mat_prod_sous11+mat_prod_12
     $              *mat_prod_sous21
               mat_prod_sous21=mat_prod_21*mat_prod_sous11+mat_prod_22
     $              *mat_prod_sous21
               mat_prod_sous11=ctmp
               ctmp=mat_prod_11*mat_prod_sous12+mat_prod_12
     $              *mat_prod_sous22
               mat_prod_sous22=mat_prod_21*mat_prod_sous12+mat_prod_22
     $              *mat_prod_sous22
               mat_prod_sous12=ctmp
               
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
               deltap_plus=(epscouche(i+1)*w(i)+epscouche(i)*w(i+1))
     $              /ctmp
               deltap_moins=(epscouche(i)*w(i+1)-epscouche(i+1)*w(i))
     $              /ctmp
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
     $              *mat_sca_sur21
               mat_sca_sur11=ctmp
               ctmp=mat_sca_11*mat_sca_sur12+mat_sca_12*mat_sca_sur22
               mat_sca_sur22=mat_sca_21*mat_sca_sur12+mat_sca_22
     $              *mat_sca_sur22
               mat_sca_sur12=ctmp
c     produit des deux matrices; vectorielle
               ctmp=mat_prod_11*mat_prod_sur11+mat_prod_12
     $              *mat_prod_sur21
               mat_prod_sur21=mat_prod_21*mat_prod_sur11+mat_prod_22
     $              *mat_prod_sur21
               mat_prod_sur11=ctmp
               ctmp=mat_prod_11*mat_prod_sur12+mat_prod_12
     $              *mat_prod_sur22
               mat_prod_sur22=mat_prod_21*mat_prod_sur12+mat_prod_22
     $              *mat_prod_sur22
               mat_prod_sur12=ctmp
               
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
            B_sca_12=nu_plus*kp2
            B_sca_21=nu_moins*w(nc)
            B_sca_22=nu_moins*kp2
            
            B_prod_11=-nu_plus*epscouche(nc)*k02/w(nc)
            B_prod_21=nu_moins*epscouche(nc)*k02/w(nc)
            
c     write(*,*) 'Bprod',B_prod_11,B_prod_21,k02

c     inverse de la matrice coupee et produit par B
            det=mat_sca_11*mat_sca_22-mat_sca_12*mat_sca_21
c     write(*,*) 'det',det,mat_sca_11,mat_sca_22,mat_sca_12,mat_sca_21
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
c     write(*,*) 'det',det
            ctmp=mat_prod_11
            mat_prod_11=mat_prod_22/det
            mat_prod_22=ctmp/det
            mat_prod_12=-mat_prod_12/det
            mat_prod_21=-mat_prod_21/det
            
            mat_prod_21=mat_prod_21*B_prod_11+mat_prod_22*B_prod_21
            mat_prod_11=mat_prod_11*B_prod_11+mat_prod_12*B_prod_21
c     write(*,*) 'B11',mat_prod_11,B_prod_11,mat_prod_12,B_prod_21
            mat_prod_22=0.d0
            mat_prod_12=0.d0
            
c     le point d'observation est dessous
            S21=mat_sca_21
            S22=mat_sca_22
            P21=mat_prod_21
c            write(*,*) 'SSS4',S21,S22,P21,kx,ky,kp2,kz
c     calcul le tenseur

            if (kp2.eq.0.d0) then
               Stenseur(1,1)=S21*w(0)
               Stenseur(1,2)=0.d0
               Stenseur(1,3)=0.d0
               Stenseur(2,1)=0.d0
               Stenseur(2,2)=S21*w(0)
               Stenseur(2,3)=0.d0
               Stenseur(3,1)=0.d0
               Stenseur(3,2)=0.d0
               Stenseur(3,3)=S22
            else
               Stenseur(1,1)=(kx*kx*S21+ky*ky*P21)/kp2*w(0)
               Stenseur(1,2)=kx*ky*(S21-P21)/kp2*w(0)
               Stenseur(1,3)=kx*S22/kp2*w(0)
               Stenseur(2,1)=Stenseur(1,2)
               Stenseur(2,2)=(ky*ky*S21+kx*kx*P21)/kp2*w(0)
               Stenseur(2,3)=ky*S22/kp2*w(0)
               Stenseur(3,1)=kx*S21
               Stenseur(3,2)=ky*S21
               Stenseur(3,3)=S22
            endif
c            write(*,*) 'Stenseur',Stenseur,kx*S12/kp2*kz
         endif

      else
         write(*,*) 'observation dans le multicouche',z,za,k0,nepsmax
     $        ,neps,dcouche,zcouche,epscouche
         stop
      endif


c     divise par epsilon de la couche ou il y a le dipole
      do i=1,3
         do j=1,3
            Stenseur(i,j)=Stenseur(i,j)/epscouche(nc)
         enddo
      enddo


      end
