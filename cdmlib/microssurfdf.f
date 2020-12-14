      subroutine microssurfdf(xi,xr,nbsphere,ndipole,nx,ny,nz,nx2,ny2
     $     ,nxm,nym,nzm,nplanm,nmatim,ntotalm,nmax,matindplan,matindice
     $     ,Tabdip,b31 ,b32 ,b33,FF,FF0,FFloc,b11,b12,b13,a11,a12,a13
     $     ,a22,a23 ,a31 ,a32 ,a33 ,WRK,epscouche,zcouche,neps,nepsmax
     $     ,xs,ys,zs,nlar,ldabi,polarisa,epsilon ,methodeit ,nrig,ncote
     $     ,tolinit ,aretecube ,npolainc ,nquicklens ,eps0,k0 ,P0 ,w0
     $     ,nfft2d,tabfft2,nproche ,Eimagexpos ,Eimageypos ,Eimagezpos,
     $     Eimageincxpos ,Eimageincypos ,Eimageinczpos, Efourierxpos,
     $     Efourierypos ,Efourierzpos, Efourierincxpos ,Efourierincypos,
     $     Efourierinczpos, Eimagexneg ,Eimageyneg ,Eimagezneg,
     $     Eimageincxneg,Eimageincyneg ,Eimageinczneg, Efourierxneg
     $     ,Efourieryneg,Efourierzneg, Efourierincxneg ,Efourierincyneg
     $     ,Efourierinczneg,Ediffkzpos,Ediffkzneg, kxy ,xy,numaperref
     $     ,numapertra,numaperinc,imaxk0,gross,zlensr,zlenst ,ntypemic ,
     $     planf ,planb ,plan2f ,plan2b ,nmatf,file_id ,group_idmic
     $     ,nstop ,infostr)

#ifdef USE_HDF5
      use HDF5
#endif

      implicit none
      integer nbsphere,ndipole,nx,ny,nz,nx2,ny2 ,nxm,nym,nzm,nplanm
     $     ,ntotalm,nmax,nmatim,nlar,ldabi,nfft2d,nfft2d2,nstop,nmatf
     $     ,ntypemic,ipol,nbsphere3,ncompte,nloop,nproche,nfft2dtmp
      integer matindice(nplanm,nmatim) ,matindplan(nzm,nzm)
     $     ,tabfft2(nfft2d)
      double complex a11(2*nxm,2*nym,nplanm),a12(2*nxm,2*nym,nplanm),
     $     a13(2*nxm,2*nym,nplanm),a22(2*nxm,2*nym,nplanm),a23(2*nxm,2
     $     *nym,nplanm),a31(2*nxm,2*nym,nplanm),a32(2*nxm,2*nym,nplanm)
     $     , a33(2*nxm,2*nym,nplanm),b11(4*nxm*nym),b12(4*nxm*nym),b13(4
     $     *nxm*nym),b22(4*nxm*nym),b23(4*nxm*nym),b31(4 *nxm*nym),b32(4
     $     *nxm*nym),b33(4*nxm*nym)
      double complex, dimension(3*nxm*nym*nzm) :: xr,xi
      double complex, dimension(3*nxm*nym*nzm,12) :: wrk
      double complex, dimension(nxm*nym*nzm,3,3) :: polarisa,epsilon
      double complex, dimension(3*nxm*nym*nzm) :: FF,FF0,FFloc
      integer, dimension(nxm*nym*nzm) :: Tabdip
      character(12) methodeit
      double precision tol,tol1,tolinit,aretecube,eps0,k0,k02,P0,irra,pi
     $     ,numaperref,numapertra,numaperinc,numaperk,kxinc,kyinc,I0
     $     ,deltax,deltakx ,deltaky,w0 ,deltak,pp,ss,tmp,rloin,sidemic
     $     ,deltakm,xmin,xmax,ymin,ymax
      integer i,j,ii,jj,k,kk,idelta,jdelta,ideltam,nrig,npolainc
     $     ,nquicklens,npol,imaxk0,niter,niterii,imul,imaxinc
      DOUBLE PRECISION,DIMENSION(nxm*nym*nzm)::xs,ys,zs
      double precision kxy(nfft2d),xy(nfft2d),gross,kx,ky,deltatheta,phi
     $     ,theta,zlensr,zlenst
      double complex Eimagexpos(nfft2d*nfft2d),Eimageypos(nfft2d*nfft2d)
     $     ,Eimagezpos(nfft2d*nfft2d),Eimageincxpos(nfft2d*nfft2d)
     $     ,Eimageincypos(nfft2d*nfft2d),Eimageinczpos(nfft2d*nfft2d)
     $     ,Efourierxpos(nfft2d*nfft2d),Efourierypos(nfft2d*nfft2d)
     $     ,Efourierzpos(nfft2d*nfft2d),Efourierincxpos(nfft2d*nfft2d)
     $     ,Efourierincypos(nfft2d*nfft2d),Efourierinczpos(nfft2d
     $     *nfft2d),Eimagexneg(nfft2d*nfft2d),Eimageyneg(nfft2d*nfft2d)
     $     ,Eimagezneg(nfft2d*nfft2d),Eimageincxneg(nfft2d*nfft2d)
     $     ,Eimageincyneg(nfft2d*nfft2d),Eimageinczneg(nfft2d*nfft2d)
     $     ,Efourierxneg(nfft2d*nfft2d),Efourieryneg(nfft2d*nfft2d)
     $     ,Efourierzneg(nfft2d*nfft2d),Efourierincxneg(nfft2d*nfft2d)
     $     ,Efourierincyneg(nfft2d*nfft2d),Efourierinczneg(nfft2d
     $     *nfft2d)
      double complex Ediffkzpos(nfft2d ,nfft2d,3),Ediffkzneg(nfft2d
     $     ,nfft2d,3)
      double complex E0,icomp,ctmp1
      integer  neps,nepsmax,ncote,ikxinc,jkyinc,indice,indicex,indicey
      double complex epscouche(0:nepsmax+1)
      double precision dcouche(nepsmax),zcouche(0:nepsmax),x,y,indice0
     $     ,indicen,z,kz
      double complex Arx,Ary,Arz,Atx,Aty,Atz,Emx,Emy,Emz,ctmp,Stenseur(3
     $     ,3),Em(3),Eloc(3),epsani(3,3),zfocus
      integer*8 planf,planb,plan2f,plan2b,planfn,planbn
      integer FFTW_FORWARD,FFTW_ESTIMATE,FFTW_BACKWARD,nsens
     $     ,numerocouche
      character(64) infostr,beam

      character(LEN=100) :: datasetname

#ifndef USE_HDF5
      integer,parameter:: hid_t=4
#endif
      integer(hid_t) :: file_id
      integer(hid_t) :: group_idmic
      integer :: dim(4)
      integer error
      
      write(*,*) 'Dark field microscope'

      
c     initialise
      pi=dacos(-1.d0)
      beam='pwavelinear'
      nfft2d2=nfft2d/2
      npolainc=0
      icomp=(0.d0,1.d0)
      nbsphere3=3*nbsphere
      numaperk=k0*numaperinc
      k02=k0*k0
      x=0.d0
      y=0.d0
      indice0=dsqrt(dreal(epscouche(0)))
      indicen=dsqrt(dreal(epscouche(neps+1)))

      if (numaperinc.ge.1.d0) then
         infostr='NA inc strictly between 0 and 1'
         nstop=1
         return
      endif

c     calcul de deltak      
      xmax=-1.d300
      xmin=1.d300
      ymax=-1.d300
      ymin=1.d300
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:xmax,ymax)
!$OMP& REDUCTION(min:xmin,ymin)      
      do i=1,nbsphere
         xmax=max(xmax,xs(i))
         xmin=min(xmin,xs(i))
         ymax=max(ymax,ys(i))
         ymin=min(ymin,ys(i))     
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      deltakm=pi/max(xmax-xmin,ymax-ymin)
      deltax=aretecube
      
      if (nquicklens.eq.1) then
         
       deltakx=2.d0*pi/(aretecube*dble(nfft2d))

         if (deltakx.ge.numaperk) then
            nstop=1
          infostr='In FFT lens nfft2d too small:Increase size of FFT 1'
            return
         endif        
         imul=idint(deltakm/deltakx)
         if (imul.eq.0) then
            nstop=1
          infostr='In FFT lens nfft2d too small:Increase size of FFT 2'
            return            
         endif
 223     deltak=deltakx*dble(imul)
         write(*,*) 'change delta k incident:',deltak,'m-1',imul
         imaxinc=nint(k0/deltak)+1
         if (imaxinc.le.3) then
            imul=imul-1
            if (imul.eq.0) then
               nstop=1
          infostr='In FFT lens nfft2d too small:Increase size of FFT 3'
               return
            endif
            goto 223
         endif
         write(*,*) 'Step size delta k diffracted: ',deltak,'m-1'
      else
         k=0
 224     deltakx=2.d0*pi/(dble(nfft2d)*aretecube)/dble(2**k)
         imaxk0=nint(k0/deltakx)+1            
         if (imaxk0.le.5) then
            k=k+1
            write(*,*) 'Change delta k diffracted:',deltakx,'m-1',k
            goto 224
         endif
         write(*,*) 'Final delta k diffracted',deltakx,'m-1'

         k=0
 222     deltak=deltakx*dnint(deltakm/deltakx)/dble(2**k)
         imaxinc=nint(k0/deltak)+1
         if (imaxinc.le.2) then
            k=k+1
            write(*,*) 'Change delta k incident:',deltak,'m-1',k
            goto 222
         endif
         write(*,*) 'Final delta k incident : ',deltak,'m-1'
      endif

      deltax=aretecube
      deltatheta=deltak/k0*numaperinc
      ideltam=max(int(2.d0*pi/deltatheta)+1,8)
      deltaky=deltakx

      if (nfft2d.gt.65536) then
         nstop=1
         infostr='nfft2d too large'
         return
      endif
      if (deltak.ge.numaperk) then
         nstop=1
         infostr='In FFT lens nfft2d too small'
         return
      endif

      write(*,*) 'Number of incidence',ideltam,deltakx,imaxinc

c     initalise
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
      do i=1,nfft2d*nfft2d
         Eimageincxpos(i)=0.d0
         Eimageincypos(i)=0.d0
         Eimageinczpos(i)=0.d0
         Eimagexpos(i)=0.d0
         Eimageypos(i)=0.d0
         Eimagezpos(i)=0.d0
         Eimageincxneg(i)=0.d0
         Eimageincyneg(i)=0.d0
         Eimageinczneg(i)=0.d0
         Eimagexneg(i)=0.d0
         Eimageyneg(i)=0.d0
         Eimagezneg(i)=0.d0   
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
      
      npol=1
      niter=ideltam*npol
      if (npolainc.eq.0) then
         npol=2
         niter=ideltam*npol
      endif

c     calcul puissance
      P0=P0/dble(npol*ideltam)
      call irradiancesurf(P0,w0,E0,irra,epscouche(0))
      I0=cdabs(E0)**2
      niterii=0
       write(*,*) 'Magnitude for each plane wave:',E0
      do ipol=1,npol
         if (npolainc.eq.1) then
            ss=1.d0
            pp=0.d0
         elseif (npolainc.eq.2) then
            ss=1.d0
            pp=1.d0
         else
            if (ipol.eq.1) then
               ss=1.d0
               pp=0.d0
            else
               ss=1.d0
               pp=1.d0
            endif
         endif

         
c     sommation 
         do idelta=0,ideltam-1
c         do idelta=0,1
            niterii=niterii+1
            write(*,*) '*** incidence',niterii,'/',niter,' *****'
            phi=dble(idelta)*2.d0*pi/dble(ideltam)
            theta=dasin(numaperinc)
            kxinc=dcos(phi)*k0*numaperinc
            kyinc=dsin(phi)*k0*numaperinc
            
c     calcul champ incident

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
            do i=1,nbsphere           
               call champlineairekxky(epscouche,zcouche,neps ,nepsmax
     $              ,xs(i),ys(i),zs(i),k0,E0,ss,pp,kxinc ,kyinc ,infostr
     $              ,nstop,FF0(3*i-2),FF0(3*i-1) ,FF0(3*i))
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL
            if (nstop.eq.1) return
c     calcul champ local
            if (nrig.eq.0) then
               tol=tolinit
               do i=1,nbsphere3
                  xi(i)=FF0(i)
               enddo
               
               if (nproche.eq.-1) then
                  call inverserigsurf(xi,xr,nbsphere,ndipole,nx,ny ,nz
     $                 ,nx2,ny2,nxm,nym,nzm,nplanm,ntotalm,nmax
     $                 ,matindplan,Tabdip,b31,b32 ,b33,FF,FF0 ,FFloc
     $                 ,b11,b12,b13,a11,a12,a13,a22,a23,a31 ,a32 ,a33
     $                 ,WRK,nlar,ldabi,polarisa ,methodeit,tol,tol1
     $                 ,nloop ,ncompte ,planf ,planb,nstop ,infostr)
                  if (nstop.eq.1) return
               else
                  call inverserigsurfopt(xi,xr,nbsphere,ndipole,nx ,ny
     $                 ,nz,nx2,ny2,nxm,nym,nzm,nplanm,ntotalm ,nmax
     $                 ,matindplan,b31,b32 ,b33,FF,FF0,FFloc ,b11,b12
     $                 ,b13,a11,a12,a13,a22,a23,a31 ,a32 ,a33 ,WRK,nlar
     $                 ,ldabi,polarisa,methodeit,tol ,tol1 ,nloop
     $                 ,ncompte ,planf,planb,nstop ,infostr)
                  if (nstop.eq.1) return
               endif
            elseif (nrig.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)          
               do i=1,nbsphere3
                  FFloc(i)=FF0(i)
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL


            elseif (nrig.eq.2) then 

               nsens=-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,kk,Em,Eloc,ii,jj,epsani,eps0,z)   
!$OMP DO SCHEDULE(STATIC) 
               do k=1,nbsphere
                  z=zs(k)
                  kk=3*(k-1)
                  Em(1)=FF0(kk+1)
                  Em(2)=FF0(kk+2)
                  Em(3)=FF0(kk+3)
                  do ii=1,3
                     do jj=1,3
                        epsani(ii,jj)=epsilon(k,ii,jj)
                     enddo
                  enddo
                  eps0=epscouche(numerocouche(z,neps,nepsmax,zcouche))
                  call local_macro_surf(Eloc,Em,epsani,eps0,aretecube,k0
     $                 ,nsens)
                  FFloc(kk+1)=Eloc(1)
                  FFloc(kk+2)=Eloc(2)
                  FFloc(kk+3)=Eloc(3)
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
         
            elseif (nrig.eq.3) then
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
               do i=1,nbsphere3
                  xr(i)=FF0(i)           
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)  
!$OMP DO SCHEDULE(STATIC)
               do i=1,nbsphere
                  k=3*(i-1)
                  xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k
     $                 +2)+polarisa(i,1,3)*xr(k+3)
                  xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k
     $                 +2)+polarisa(i,2,3)*xr(k+3)
                  xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k
     $                 +2)+polarisa(i,3,3)*xr(k+3)
               enddo        
!$OMP ENDDO 
!$OMP END PARALLEL

               call produitfftmatvectsurplus(xi,xr,nbsphere,ndipole,nx
     $              ,ny,nz,nx2,ny2,nxm,nym,nzm,nzm,nplanm,ntotalm,nmax
     $              ,matindplan,Tabdip,b31,b32,b33,FF,b11,b12,b13,a11
     $              ,a12,a13,a22,a23,a31,a32 ,a33,planb,planf)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)       
               do i=1,nbsphere
                  k=3*(i-1)
                  FFloc(k+1)=xr(k+1)
                  FFloc(k+2)=xr(k+2)
                  FFloc(k+3)=xr(k+3)       
               enddo            
!$OMP ENDDO 
!$OMP END PARALLEL      
               
            endif

c     dipole a partir champ local
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC)         
            do i=1,nbsphere
               k=3*(i-1)
               FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1,2)
     $              *FFloc(k+2)+polarisa(i,1,3)*FFloc(k+3)
               FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2,2)
     $              *FFloc(k+2)+polarisa(i,2,3)*FFloc(k+3)
               FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3,2)
     $              *FFloc(k+2)+polarisa(i,3,3)*FFloc(k+3) 
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL  

            nfft2dtmp=nfft2d
            if (nquicklens.eq.1) then
               rloin=1.d0
c               write(*,*) 'ff local',FF
               call diffractefft2dsurf2(nbsphere,nx,ny,nz,nxm,nym ,nzm
     $              ,nfft2dtmp,nfft2d,tabfft2,k0,xs,ys,zs,aretecube
     $              ,Efourierxpos,Efourierypos,Efourierzpos,FF ,imaxk0
     $              ,deltakx,deltaky,Ediffkzpos,Ediffkzneg ,rloin,rloin
     $              ,numaperref,numapertra,nepsmax ,neps,dcouche,zcouche
     $              ,epscouche,ncote ,nstop ,infostr,plan2f)
               if (nstop.eq.1) return
            else
c     compute the diffracted field
               do i=-imaxk0,imaxk0               
                  if (i.ge.0) then
                     indicex=i+1
                  else
                     indicex=nfft2d+i+1
                  endif
                  kx=deltakx*dble(i)
                  do j=-imaxk0,imaxk0
                     
                     if (j.ge.0) then
                        indicey=j+1
                     else
                        indicey=nfft2d+j+1
                     endif
                     ky=deltaky*dble(j)
                     ii=imaxk0+i+1
                     jj=imaxk0+j+1

                     if (ncote.eq.1.or.ncote.eq.0) then
c     calcul champ dessus
                        
                        if (k02*indicen*indicen*numapertra*numapertra
     $                       *0.9999d0-kx*kx-ky*ky.gt.0.d0) then  
                           z=1.d0
                           kz=dsqrt(k02*indicen*indicen-kx*kx-ky*ky)

                           Emx=0.d0
                           Emy=0.d0
                           Emz=0.d0
                           call  tenseurmulticoucheloinfft(kx,ky ,kz,z
     $                          ,zs(1),k0,nepsmax,neps ,dcouche,zcouche
     $                          ,epscouche ,Stenseur)
                           ctmp=cdexp(-icomp*(kx*xs(1)+ky*ys(1)))
                           Emx=(Stenseur(1,1)*FF(1)+Stenseur(1,2) *FF(2)
     $                          +Stenseur(1,3)*FF(3))*ctmp
                           Emy=(Stenseur(2,1)*FF(1)+Stenseur(2,2) *FF(2)
     $                          +Stenseur(2,3)*FF(3))*ctmp
                           Emz=(Stenseur(3,1)*FF(1)+Stenseur(3,2) *FF(2)
     $                          +Stenseur(3,3)*FF(3))*ctmp
                           do k=2,nbsphere
                              kk=3*(k-1)
                              if (zs(k).ne.zs(k-1)) then
                                 call tenseurmulticoucheloinfft(kx ,ky
     $                                ,kz,z,zs(k),k0,nepsmax ,neps
     $                                ,dcouche,zcouche ,epscouche
     $                                ,Stenseur) 
                              endif
                              ctmp=cdexp(-icomp*(kx*xs(k)+ky *ys(k)))
                              Emx=Emx+(Stenseur(1,1)*FF(kk+1)
     $                             +Stenseur(1,2)*FF(kk+2) +Stenseur(1
     $                             ,3)*FF(kk+3))*ctmp
                              Emy=Emy+(Stenseur(2,1)*FF(kk+1)
     $                             +Stenseur(2,2)*FF(kk+2) +Stenseur(2
     $                             ,3)*FF(kk+3))*ctmp
                              Emz=Emz+(Stenseur(3,1)*FF(kk+1)
     $                             +Stenseur(3,2)*FF(kk+2) +Stenseur(3
     $                             ,3)*FF(kk+3))*ctmp
                           enddo        
                           Ediffkzpos(ii,jj,1)=Emx
                           Ediffkzpos(ii,jj,2)=Emy
                           Ediffkzpos(ii,jj,3)=Emz
                           
                        endif
                     endif

c     calcul champ dessous
                     if (ncote.eq.-1.or.ncote.eq.0) then
                        if (k02*indice0*indice0-kx*kx-ky*ky.gt.0.d0)
     $                       then    
                           z=-1.d0
                           kz=-dsqrt(k02*indice0*indice0-kx*kx-ky*ky)
                           call  tenseurmulticoucheloinfft(kx,ky,kz,z
     $                          ,zs(1),k0,nepsmax,neps ,dcouche,zcouche
     $                          ,epscouche ,Stenseur)
                           ctmp=cdexp(-icomp*(kx*xs(1)+ky*ys(1)))
                           Emx=(Stenseur(1,1)*FF(1)+Stenseur(1,2) *FF(2)
     $                          +Stenseur(1,3)*FF(3))*ctmp
                           Emy=(Stenseur(2,1)*FF(1)+Stenseur(2,2) *FF(2)
     $                          +Stenseur(2,3)*FF(3))*ctmp
                           Emz=(Stenseur(3,1)*FF(1)+Stenseur(3,2) *FF(2)
     $                          +Stenseur(3,3)*FF(3))*ctmp
                           do k=2,nbsphere
                              kk=3*(k-1)
                              if (zs(k).ne.zs(k-1)) then
                                 call tenseurmulticoucheloinfft(kx ,ky
     $                                ,kz,z,zs(k),k0,nepsmax ,neps
     $                                ,dcouche,zcouche ,epscouche
     $                                ,Stenseur)
                              endif
                              ctmp=cdexp(-icomp*(kx*xs(k)+ky *ys(k)))
                              Emx=Emx+(Stenseur(1,1)*FF(kk+1)
     $                             +Stenseur(1,2)*FF(kk+2) +Stenseur(1
     $                             ,3)*FF(kk+3))*ctmp
                              Emy=Emy+(Stenseur(2,1)*FF(kk+1)
     $                             +Stenseur(2,2)*FF(kk+2) +Stenseur(2
     $                             ,3)*FF(kk+3))*ctmp
                              Emz=Emz+(Stenseur(3,1)*FF(kk+1)
     $                             +Stenseur(3,2)*FF(kk+2) +Stenseur(3
     $                             ,3)*FF(kk+3))*ctmp
                           enddo
                           Ediffkzneg(ii,jj,1)=Emx
                           Ediffkzneg(ii,jj,2)=Emy
                           Ediffkzneg(ii,jj,3)=Emz
                        endif
                     endif
                  enddo
               enddo
            endif

c     calcul image
c     passe le champ diffracte lointain en amplitude e(k||)
            call diffractefft2dtoeposfour(Ediffkzpos,Ediffkzneg
     $           ,Efourierxpos,Efourierypos,Efourierzpos ,Efourierxneg
     $           ,Efourieryneg,Efourierzneg ,epscouche,nepsmax,neps
     $           ,numaperref,numapertra,k0,deltax ,deltakx,imaxk0
     $           ,nfft2dtmp,nfft2d ,ncote ,nstop ,infostr)
            if (nstop.eq.1) return
            tmp=indice0*k0
            call deltakroutine(kxinc,kyinc,deltakx,deltaky,tmp,ikxinc
     $           ,jkyinc)
c            write(*,*) 'ff ij',kxinc,kyinc,deltakx,deltaky,k0,ikxinc
c     $           ,jkyinc
            if (ncote.eq.0.or.ncote.eq.-1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
               do i=1,nfft2d*nfft2d
                  Efourierincxneg(i)=0.d0
                  Efourierincyneg(i)=0.d0
                  Efourierinczneg(i)=0.d0
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indicex,indicey,indice)
!$OMP& PRIVATE(kx,ky,Arx,Ary,Arz,Atx,Aty,Atz,zfocus)       
!$OMP DO SCHEDULE(STATIC)  COLLAPSE(2)                 
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0
                     kx=dble(i)*deltakx
                     if (i.ge.0) then
                        indicex=i+1
                     else
                        indicex=nfft2d+i+1
                     endif
                     ky=deltaky*dble(j)
                     if (j.ge.0) then
                        indicey=j+1
                     else
                        indicey=nfft2d+j+1
                     endif
                     if (indice0*indice0*k02*numaperref*numaperref
     $                    *0.9999d0-kx*kx-ky*ky.gt.0.d0) then
                        kz=dsqrt(indice0*indice0*k02-kx*kx-ky*ky)
                        zfocus=cdexp(-icomp*kz*zlensr)
                        indice=indicex+nfft2d*(indicey-1)
                        Efourierxneg(indice)=Efourierxneg(indice)
     $                       *zfocus
                        Efourieryneg(indice)=Efourieryneg(indice)
     $                       *zfocus
                        Efourierzneg(indice)=Efourierzneg(indice)
     $                       *zfocus

                        if (i.eq.ikxinc.and.j.eq.jkyinc) then
                           call champlineairemicrokxky(epscouche,
     $                          zcouche,neps,nepsmax,x,y,k0,E0 ,ss ,pp
     $                          ,kx,ky,infostr ,nstop ,Arx,Ary,Arz ,Atx
     $                          ,Aty ,Atz)
                           Efourierincxneg(indice) =Efourierxneg(indice)
     $                          +Arx/deltakx/deltakx*icomp*zfocus
                           Efourierincyneg(indice) =Efourieryneg(indice)
     $                          +Ary/deltakx/deltakx*icomp*zfocus
                           Efourierinczneg(indice) =Efourierzneg(indice)
     $                          +Arz/deltakx/deltakx*icomp*zfocus
                        else
                           Efourierincxneg(indice) =Efourierxneg(indice)
                           Efourierincyneg(indice) =Efourieryneg(indice)
                           Efourierinczneg(indice) =Efourierzneg(indice)
                        endif                        
                     endif
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL
               if (nstop.eq.1) return

               if (gross.eq.-1.d0) then
                  
                  call passagefourierimage2(Efourierincxneg
     $                 ,Efourierincyneg,Efourierinczneg ,nfft2d ,nfft2d
     $                 ,imaxk0,indice0,deltakx ,deltax ,plan2b)
                  call passagefourierimage2(Efourierxneg ,Efourieryneg
     $                 ,Efourierzneg,nfft2d,nfft2d ,imaxk0,indice0
     $                 ,deltakx ,deltax,plan2b)
                  
               else
                  sidemic=-1.d0
                  call passagefourierimagegross2(Efourierincxneg
     $                 ,Efourierincyneg,Efourierinczneg,nfft2d ,nfft2d
     $                 ,imaxk0 ,deltakx,deltax,gross,k0,indice0
     $                 ,numaperref,sidemic ,plan2f,plan2b)
                  call passagefourierimagegross2(Efourierxneg
     $                 ,Efourieryneg,Efourierzneg,nfft2d,nfft2d ,imaxk0
     $                 ,deltakx,deltax ,gross,k0,indice0,numaperref
     $                 ,sidemic,plan2f ,plan2b)
               endif
            endif
            if (ncote.eq.0.or.ncote.eq.1) then
c     ajoute onde plane
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
               do i=1,nfft2d*nfft2d
                  Efourierincxpos(i)=0.d0
                  Efourierincypos(i)=0.d0
                  Efourierinczpos(i)=0.d0
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL

               
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indicex,indicey,indice)
!$OMP& PRIVATE(kx,ky,Arx,Ary,Arz,Atx,Aty,Atz,zfocus)
!$OMP DO SCHEDULE(STATIC)  COLLAPSE(2)                    
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0
                     kx=dble(i)*deltakx
                     if (i.ge.0) then
                        indicex=i+1
                     else
                        indicex=nfft2d+i+1
                     endif
                     ky=deltaky*dble(j)
                     if (j.ge.0) then
                        indicey=j+1
                     else
                        indicey=nfft2d+j+1
                     endif
                     if (indicen*indicen*k02*numapertra*numapertra
     $                    *0.9999d0-kx*kx-ky*ky.gt.0.d0) then
                        kz=dsqrt(indicen*indicen*k02-kx*kx-ky *ky)
                        zfocus=cdexp(icomp*kz*zlenst)
                        indice=indicex+nfft2d*(indicey-1)
                        Efourierxpos(indice)=Efourierxpos(indice)
     $                       *zfocus
                        Efourierypos(indice)=Efourierypos(indice)
     $                       *zfocus
                        Efourierzpos(indice)=Efourierzpos(indice)
     $                       *zfocus
                        
                        if (i.eq.ikxinc.and.j.eq.jkyinc) then
                           call champlineairemicrokxky(epscouche,
     $                          zcouche,neps,nepsmax,x,y,k0,E0 ,ss ,pp
     $                          ,kx,ky,infostr ,nstop ,Arx,Ary,Arz ,Atx
     $                          ,Aty ,Atz)
                           Efourierincxpos(indice) =Efourierxpos(indice)
     $                          +Atx/deltakx /deltakx*icomp*zfocus
                           Efourierincypos(indice) =Efourierypos(indice)
     $                          +Aty/deltakx /deltakx*icomp*zfocus
                           Efourierinczpos(indice) =Efourierzpos(indice)
     $                          +Atz/deltakx /deltakx*icomp*zfocus
                        else
                           Efourierincxpos(indice) =Efourierxpos(indice)
                           Efourierincypos(indice) =Efourierypos(indice)
                           Efourierinczpos(indice) =Efourierzpos(indice)
                        endif
                     endif
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL
               if (nstop.eq.1) return

               if (gross.eq.-1.d0) then
                  call passagefourierimage2(Efourierincxpos
     $                 ,Efourierincypos,Efourierinczpos,nfft2d ,nfft2d
     $                 ,imaxk0,indicen ,deltakx ,deltax ,plan2b)
                  call passagefourierimage2(Efourierxpos ,Efourierypos
     $                 ,Efourierzpos,nfft2d,nfft2d ,imaxk0,indicen
     $                 ,deltakx ,deltax,plan2b)
               else
                  sidemic=1.d0
                  call passagefourierimagegross2(Efourierincxpos
     $                 ,Efourierincypos,Efourierinczpos ,nfft2d ,nfft2d
     $                 ,imaxk0 ,deltakx,deltax ,gross,k0 ,indicen
     $                 ,numapertra ,sidemic,plan2f,plan2b)
                  call passagefourierimagegross2(Efourierxpos
     $                 ,Efourierypos,Efourierzpos,nfft2d,nfft2d ,imaxk0
     $                 ,deltakx,deltax ,gross,k0,indicen,numapertra
     $                 ,sidemic,plan2f ,plan2b)
               endif
            endif
c     sommation de toutes les images incohérentes.

            if (ncote.eq.0.or.ncote.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(indice)   
!$OMP DO SCHEDULE(STATIC)    
               do indice=1,nfft2d*nfft2d
                  Eimagexpos(indice)=Efourierxpos(indice)
     $                 *dconjg(Efourierxpos(indice)) +Eimagexpos(indice)
                  Eimageypos(indice)=Efourierypos(indice)
     $                 *dconjg(Efourierypos(indice)) +Eimageypos(indice)
                  Eimagezpos(indice)=Efourierzpos(indice)
     $                 *dconjg(Efourierzpos(indice)) +Eimagezpos(indice)
                  Eimageincxpos(indice)=Efourierincxpos(indice)
     $                 *dconjg(Efourierincxpos(indice))
     $                 +Eimageincxpos(indice)
                  Eimageincypos(indice)=Efourierincypos(indice)
     $                 *dconjg(Efourierincypos(indice))
     $                 +Eimageincypos(indice)
                  Eimageinczpos(indice)=Efourierinczpos(indice)
     $                 *dconjg(Efourierinczpos(indice))
     $                 +Eimageinczpos(indice)
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL
            endif
            if (ncote.eq.0.or.ncote.eq.-1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(indice)   
!$OMP DO SCHEDULE(STATIC)    
               do indice=1,nfft2d*nfft2d
                  Eimagexneg(indice)=Efourierxneg(indice)
     $                 *dconjg(Efourierxneg(indice)) +Eimagexneg(indice)
                  Eimageyneg(indice)=Efourieryneg(indice)
     $                 *dconjg(Efourieryneg(indice)) +Eimageyneg(indice)
                  Eimagezneg(indice)=Efourierzneg(indice)
     $                 *dconjg(Efourierzneg(indice)) +Eimagezneg(indice)
                  Eimageincxneg(indice)=Efourierincxneg(indice)
     $                 *dconjg(Efourierincxneg(indice))
     $                 +Eimageincxneg(indice)
                  Eimageincyneg(indice)=Efourierincyneg(indice)
     $                 *dconjg(Efourierincyneg(indice))
     $                 +Eimageincyneg(indice)
                  Eimageinczneg(indice)=Efourierinczneg(indice)
     $                 *dconjg(Efourierinczneg(indice))
     $                 +Eimageinczneg(indice)
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL
            endif
            if (nstop == -1) then
               infostr = 'Calculation cancelled during iterative method'
               return
            endif
         enddo         
      enddo


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(indice)   
!$OMP DO SCHEDULE(STATIC)        
      do indice=1,nfft2d*nfft2d
         Eimagexpos(indice)=cdsqrt(Eimagexpos(indice))
         Eimageypos(indice)=cdsqrt(Eimageypos(indice))
         Eimagezpos(indice)=cdsqrt(Eimagezpos(indice))
         Eimageincxpos(indice)=cdsqrt(Eimageincxpos(indice))
         Eimageincypos(indice)=cdsqrt(Eimageincypos(indice))
         Eimageinczpos(indice)=cdsqrt(Eimageinczpos(indice))
         Eimagexneg(indice)=cdsqrt(Eimagexneg(indice))
         Eimageyneg(indice)=cdsqrt(Eimageyneg(indice))
         Eimagezneg(indice)=cdsqrt(Eimagezneg(indice))
         Eimageincxneg(indice)=cdsqrt(Eimageincxneg(indice))
         Eimageincyneg(indice)=cdsqrt(Eimageincyneg(indice))
         Eimageinczneg(indice)=cdsqrt(Eimageinczneg(indice))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      deltax=2.d0*pi/dble(nfft2d)/deltakx

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)   
      do i=-nfft2d2,nfft2d2-1
         xy(i+nfft2d2+1)=deltax*dble(i)*gross
         kxy(i+nfft2d2+1)=deltakx*dble(i)/k0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      if (nmatf.eq.0) then
         open(401,file='ximage.mat')
         open(402,file='yimage.mat')

         do i=-nfft2d2,nfft2d2-1
            write(401,*) xy(i+nfft2d2+1)
            write(402,*) xy(i+nfft2d2+1)
         enddo
         close(401)
         close(402)
         do i=-imaxk0,imaxk0
            kx=deltakx*dble(i)
            write(650,*) kx/k0
            write(651,*) kx/k0
         enddo

         open(410,file='kxincidentdf.mat')
         open(411,file='kyincidentdf.mat')
         do idelta=0,ideltam-1
            phi=dble(idelta)*2.d0*pi/dble(ideltam)
            theta=dasin(numaperinc)
            kxinc=dcos(phi)*numaperinc
            kyinc=dsin(phi)*numaperinc
            write(410,*) kxinc
            write(411,*) kyinc
         enddo
         close(410)
         close(411)

         if (ncote.eq.0.or.ncote.eq.1) then
            
            open(301,file='imagedfpos.mat')
            open(302,file='imagedfxpos.mat')
            open(303,file='imagedfypos.mat')
            open(304,file='imagedfzpos.mat')
            open(305,file='imageincdfpos.mat')
            open(306,file='imageincdfxpos.mat')
            open(307,file='imageincdfypos.mat')
            open(308,file='imageincdfzpos.mat')

            do i=1,nfft2D*nfft2D
               write(301,*) dsqrt(dreal(Eimagexpos(i)**2.d0
     $              +Eimageypos(i)**2.d0+Eimagezpos(i)**2.d0))
               write(302,*) dreal(Eimagexpos(i))
               write(303,*) dreal(Eimageypos(i))
               write(304,*) dreal(Eimagezpos(i))
               write(305,*) dsqrt(dreal(Eimageincxpos(i)**2.d0
     $              +Eimageincypos(i)**2.d0+Eimageinczpos(i)**2.d0))
               write(306,*) dreal(Eimageincxpos(i))
               write(307,*) dreal(Eimageincypos(i))
               write(308,*) dreal(Eimageinczpos(i))
            enddo
            close(301)
            close(302)
            close(303)
            close(304)
            close(305)
            close(306)
            close(307)
            close(308)
         endif

         if (ncote.eq.0.or.ncote.eq.-1) then
            
            open(301,file='imagedfneg.mat')
            open(302,file='imagedfxneg.mat')
            open(303,file='imagedfyneg.mat')
            open(304,file='imagedfzneg.mat')
            open(305,file='imageincdfneg.mat')
            open(306,file='imageincdfxneg.mat')
            open(307,file='imageincdfyneg.mat')
            open(308,file='imageincdfzneg.mat')

            do i=1,nfft2D*nfft2D
               write(301,*) dsqrt(dreal(Eimagexneg(i)**2.d0
     $              +Eimageyneg(i)**2.d0+Eimagezneg(i)**2.d0))
               write(302,*) dreal(Eimagexneg(i))
               write(303,*) dreal(Eimageyneg(i))
               write(304,*) dreal(Eimagezneg(i))
               write(305,*) dsqrt(dreal(Eimageincxneg(i)**2.d0
     $              +Eimageincyneg(i)**2.d0+Eimageinczneg(i)**2.d0))
               write(306,*) dreal(Eimageincxneg(i))
               write(307,*) dreal(Eimageincyneg(i))
               write(308,*) dreal(Eimageinczneg(i))
            enddo
            close(301)
            close(302)
            close(303)
            close(304)
            close(305)
            close(306)
            close(307)
            close(308)
         endif
      elseif (nmatf.eq.2) then
         dim(1)=nfft2d
         dim(2)=nfft2d
         datasetname='x Image'
         call hdf5write1d(group_idmic,datasetname,xy,dim)
         k=0
         if (ncote.eq.0.or.ncote.eq.1) then
            datasetname='Image dark field kz>0'
            call writehdf5mic(Eimagexpos,Eimageypos,Eimagezpos,nfft2d
     $           ,imaxk0,Ediffkzpos,k,datasetname,group_idmic)
            datasetname='Image+incident dark field kz>0'
            call writehdf5mic(Eimageincxpos,Eimageincypos,Eimageinczpos
     $           ,nfft2d,imaxk0,Ediffkzpos,k,datasetname,group_idmic)
         endif
         if (ncote.eq.0.or.ncote.eq.-1) then
            datasetname='Image dark field kz<0'
            call writehdf5mic(Eimagexneg,Eimageyneg,Eimagezneg,nfft2d
     $           ,imaxk0,Ediffkzpos,k,datasetname,group_idmic)
            datasetname='Image+incident dark field kz<0'
            call writehdf5mic(Eimageincxneg,Eimageincyneg,Eimageinczneg
     $           ,nfft2d,imaxk0,Ediffkzpos,k,datasetname,group_idmic)
         endif

         do idelta=0,ideltam-1
            phi=dble(idelta)*2.d0*pi/dble(ideltam)
            theta=dasin(numaperinc)
            kxy(idelta+1)=dcos(phi)*numaperinc
         enddo
         dim(1)=ideltam
         dim(2)=nfft2d
         datasetname='kx incident df'
         call hdf5write1d(group_idmic,datasetname,kxy,dim)
         do idelta=0,ideltam-1
            phi=dble(idelta)*2.d0*pi/dble(ideltam)
            theta=dasin(numaperinc)
            kxy(idelta+1)=dsin(phi)*numaperinc
         enddo
         dim(1)=ideltam
         dim(2)=nfft2d
         datasetname='ky incident df'
         call hdf5write1d(group_idmic,datasetname,kxy,dim)
         
      endif
      

      
      end
      
