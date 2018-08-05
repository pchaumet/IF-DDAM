c     Author : l'equipe triso
c     
c     routine qui calcule le champ lointain total reflechi (Ediffkzneg)
c     et transmis (Ediffkzpos) qui correspond à E(r).

c     Attention ici on calcule E(r)=Int[e(k||) exp(i k||.r|| +ikz z) ]
c     rappel : phase stationnaire E(r)=-2*icomp*pi*kz*e(k||)*exp(ik0r)/r
c     Ici on travaille avec e(k||)

      subroutine champtotallineairemulti(imax,ncote,nfft2d,nfft2dmax
     $     ,nepsmax,neps,epscouche,k0,aretecube,NA,Ediffkzpos,Ediffkzneg
     $     ,zcouche,E0m,ssm,ppm,thetam,phim,nba,Egausxref ,Egausyref
     $     ,Egauszref ,Egausxtra,Egausytra ,Egausztra,fluxref ,fluxtrans
     $     ,fluxinc,nstop ,infostr)
      implicit none
      integer i,j,ii,jj,nfft2d,nfft2d2,nfft2dmax,ncote ,imax,indice
     $     ,nstop,im
      double precision kx,ky,kz,kp2,fluxref,fluxtrans,tmp,deltakx
     $     ,deltaky,aretecube,k0,NA,kxinc,kyinc,x,y,fluxinc
      double complex Ediffkzpos(nfft2dmax,nfft2dmax,3)
     $     ,Ediffkzneg(nfft2dmax,nfft2dmax,3),Egausxref(nfft2dmax
     $     *nfft2dmax),Egausyref(nfft2dmax *nfft2dmax)
     $     ,Egauszref(nfft2dmax*nfft2dmax),Egausxtra(nfft2dmax
     $     *nfft2dmax),Egausytra(nfft2dmax*nfft2dmax)
     $     ,Egausztra(nfft2dmax*nfft2dmax)
      integer nepsmax,neps,imini,jmini,nba
      double precision indice0,indicen,pi
      double precision phim(20) ,thetam(20),ssm(20),ppm(20)
      double complex epscouche(0:nepsmax+1),zcouche(0:nepsmax),icomp
     $     ,const,Arx,Ary,Arz,Atx,Aty,Atz,E0m(20)
      character(64) infostr

c     initialisation
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      indice0=dsqrt(dreal(epscouche(0)))
      indicen=dsqrt(dreal(epscouche(neps+1)))
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      nfft2d2=nfft2d/2  
      fluxtrans=0.d0
      fluxref=0.d0
      fluxinc=0.d0
      x=0.d0
      y=0.d0

c     Dans le domaine de fourier: je somme le flux de toutes les ondes
c     incidentes
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kx,ky,kp2,kz,im)   
!$OMP& PRIVATE(imini,jmini,kxinc,kyinc)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:fluxinc) COLLAPSE(2)
      do i=-imax,imax
         do j=-imax,imax
            kx=dble(i)*deltakx
            ky=dble(j)*deltaky
            kp2=kx*kx+ky*ky
            do im=1,nba
               kxinc=k0*dsin(thetam(im)*pi/180.d0)*dcos(phim(im)*pi
     $              /180.d0)*indice0
               kyinc=k0*dsin(thetam(im)*pi/180.d0)*dsin(phim(im)*pi
     $              /180.d0)*indice0
c     calcul du kx et ky le plus proche
               imini=nint(kxinc/deltakx)
               jmini=nint(kyinc/deltaky)
               if (indice0*indice0*k0*k0*NA*NA-kp2.gt.0.d0) then
                  kz=dsqrt(indice0*indice0*k0*k0-kp2)               
                  if (i.eq.imini.and.j.eq.jmini) then
                     fluxinc=fluxinc+(cdabs(E0m(im)/deltakx/deltaky)
     $                    **2.d0)*kz
                  endif
               endif
            enddo
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

c     calcul du champ total
      if (ncote.eq.0.or.ncote.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
         do i=1,nfft2d*nfft2d
            Egausxtra(i)=0.d0
            Egausytra(i)=0.d0
            Egausztra(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kp2,kz,im,jmini) 
!$OMP& PRIVATE(indice,const,Arx,Ary,Arz,Atx,Aty,Atz,kxinc,kyinc,imini)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:fluxtrans) COLLAPSE(2)
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               kp2=kx*kx+ky*ky
               if (indicen*indicen*k0*k0*NA*NA-kp2.gt.0.d0)
     $              then
                  kz=dsqrt(indicen*indicen*k0*k0-kp2)
                  
                  indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                  const=(-2.d0*icomp*pi*kz)
                  
                  Egausxtra(indice)=Ediffkzpos(ii,jj,1)/const
                  Egausytra(indice)=Ediffkzpos(ii,jj,2)/const
                  Egausztra(indice)=Ediffkzpos(ii,jj,3)/const
                  do im=1,nba
c     calcul le kx et ky correspondant au theta et phi
                     kxinc=k0*dsin(thetam(im)*pi/180.d0)*dcos(phim(im)
     $                    *pi/180.d0)*indice0
                     kyinc=k0*dsin(thetam(im)*pi/180.d0)*dsin(phim(im)
     $                    *pi/180.d0)*indice0
c     calcul du kx et ky le plus proche
                     imini=nint(kxinc/deltakx)
                     jmini=nint(kyinc/deltaky)
                     if (i.eq.imini.and.j.eq.jmini) then
                        write(*,*) 'iminiii'
                        call  champlineairemicro(epscouche,zcouche ,neps
     $                       ,nepsmax,x,y,k0,E0m(im),ssm(im),ppm(im)
     $                       ,thetam(im),phim(im),infostr,nstop,Arx,Ary
     $                       ,Arz,Atx,Aty ,Atz)
                        Egausxtra(indice)=Egausxtra(indice)+Atx/deltakx
     $                       /deltaky
                        Egausytra(indice)=Egausytra(indice)+Aty/deltakx
     $                       /deltaky
                        Egausztra(indice)=Egausztra(indice)+Atz/deltakx
     $                       /deltaky
                        write(*,*) 'centre',Egausxtra(indice),Atx
     $                       /deltakx/deltaky,Egausytra(indice),Aty
     $                       /deltakx/deltaky,Egausztra(indice),Atz
     $                       /deltakx/deltaky,Atx,Aty,Atz
                     endif
                  enddo
c     write(*,*) 'champ00',kz/k0,Egausxtra(indice)
c     $                 ,Egausytra(indice),Egausztra(indice),i,j,indice
c     $                 ,Ediffkzpos(ii,jj,1)/const,Ediffkzpos(ii,jj,2)
c     $                 /const ,Ediffkzpos(ii,jj,3)/const
c     write(*,*) 'rrr',Ediffkzpos(ii,jj,1),Egausxtra(indice)
c     $                 ,Ediffkzpos(ii,jj,2),Egausytra(indice)
c     $                 ,Ediffkzpos(ii,jj,3),Egausztra(indice),indice

                  fluxtrans=fluxtrans+(cdabs(Egausxtra(indice))**2
     $                 +cdabs(Egausytra(indice))**2
     $                 +cdabs(Egausztra(indice))**2)*kz
c     write(*,*) 'rrrrr',kzt,fluxtrans
c     write(*,*) 'champ1',Egausxtra(indice)
c     $                 ,Egausytra(indice),Egausztra(indice),i,j,indice
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL         
      endif

      if (ncote.eq.0.or.ncote.eq.-1) then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
         do i=1,nfft2d*nfft2d
            Egausxref(i)=0.d0
            Egausyref(i)=0.d0
            Egauszref(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,kx,ky,kp2,kz,im,jmini) 
!$OMP& PRIVATE(indice,const,Arx,Ary,Arz,Atx,Aty,Atz,kxinc,kyinc,imini)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:fluxref) COLLAPSE(2)
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               kp2=kx*kx+ky*ky
               if (indice0*indice0*k0*k0*NA*NA-kp2.gt.0.d0)
     $              then
                  kz=dsqrt(indice0*indice0*k0*k0-kp2)

                  indice=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                  const=(-2.d0*icomp*pi*kz)



                  Egausxref(indice)=Ediffkzneg(ii,jj,1)/const
                  Egausyref(indice)=Ediffkzneg(ii,jj,2)/const
                  Egauszref(indice)=Ediffkzneg(ii,jj,3)/const
                  do im=1,nba
c     calcul le kx et ky correspondant au theta et phi
                     kxinc=k0*dsin(thetam(im)*pi/180.d0)*dcos(phim(im)
     $                    *pi/180.d0)
                     kyinc=k0*dsin(thetam(im)*pi/180.d0)*dsin(phim(im)
     $                    *pi/180.d0)
c     calcul du kx et ky le plus proche
                     imini=nint(kxinc/deltakx)
                     jmini=nint(kyinc/deltaky)
c                     write(*,*) 'imini',imini,jmini,kxinc/deltakx
c     $                    ,dasin(dble(imini)*deltakx/k0)*180.d0/pi
                     if (i.eq.imini.and.j.eq.jmini) then
                        call  champlineairemicro(epscouche,zcouche ,neps
     $                       ,nepsmax,x,y,k0,E0m(im),ssm(im),ppm(im)
     $                       ,thetam(im),phim(im),infostr,nstop,Arx,Ary
     $                       ,Arz,Atx,Aty,Atz)
                        Egausxref(indice)=Egausxref(indice)+Arx/deltakx
     $                       /deltaky
                        Egausyref(indice)=Egausyref(indice)+Ary/deltakx
     $                       /deltaky
                        Egauszref(indice)=Egauszref(indice)+Arz/deltakx
     $                       /deltaky
                        
                     endif
                  enddo
                  fluxref=fluxref+(cdabs(Egausxref(indice))**2
     $                 +cdabs(Egausyref(indice))**2
     $                 +cdabs(Egauszref(indice))**2)*kz
c     write(*,*) 'fffrrr' ,fluxref,kz
c     write(*,*) 'champ2',Egausxref(indice)
c     $                 ,Egausyref(indice),Egauszref(indice),i,j
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL         
      endif
      
      
      tmp=4.d0*pi*pi*deltaky*deltakx/(k0*8.d0*pi*1.d-7*299792458.d0)

      fluxinc=fluxinc*tmp
      fluxref=fluxref*tmp
      fluxtrans=fluxtrans*tmp
      write(*,*) 'nnn',fluxref,fluxtrans,tmp
      return
      end
