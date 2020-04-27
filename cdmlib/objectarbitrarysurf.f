      subroutine objetarbitrarysurf(trope,eps,epsani,eps0,xs,ys,zs,xswf
     $     ,yswf,zswf,k0,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx
     $     ,ny,nz,methode,namefile,na,epsilon,polarisa,neps,nepsmax
     $     ,dcouche,zcouche,epscouche,tabzn,nmatf,file_id,group_iddip
     $     ,infostr,nstop)
#ifdef USE_HDF5
      use HDF5
#endif
      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,ii,jj,i,j,k
     $     ,test,IP(3),nnnr,dddis,inv,nx1,ny1,nz1,na,nstop,ierror
      double precision xs(nmax),ys(nmax),zs(nmax),xswf(nmax),yswf(nmax)
     $     ,zswf(nmax),k0,x,y,z,aretecube,zg
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),ctmp,eps0
      integer neps,nepsmax,nminc,nmaxc,numerocouche,tabzn(nmax),nmatf
      double precision dcouche(nepsmax),zcouche(0:nepsmax),zmin,zmax
     $     ,xmin,ymin
      double complex epscouche(0:nepsmax+1)
      character(2) methode
      character(3) trope
      character(64) namefile
      character(64) infostr
#ifndef USE_HDF5
      integer,parameter:: hid_t=4
#endif
      character(LEN=100) :: datasetname
      integer(hid_t) :: file_id
      integer(hid_t) :: group_iddip
      integer :: dim(4)
      integer error
      
c     Initialization
      nbsphere=0
      ndipole=0 
      Tabdip=0
      Tabzn=0
      polarisa=0.d0
      epsilon=0.d0
      dddis=1
      inv=1
c     mesh 
      if (nmatf.eq.0) then
         open(20,file='x.mat')
         open(21,file='y.mat')
         open(22,file='z.mat')  
c     discretization of the object under study
         open(10,file='xc.mat')
         open(11,file='yc.mat')
         open(12,file='zc.mat')  
      endif
c     read the input file
      open(15,file=namefile,status='old',iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary object: name of file does not exist'
         nstop=1
         return
      endif
      read(15,*) nx,ny,nz
      read(15,*) aretecube
      aretecube=aretecube*1.d-9
      ndipole=nx*ny*nz
      if (ndipole.gt.nmax) then
         write(*,*) 'Size of the parameter too small',ndipole
     $        ,'should be smaller that',nmax
         write(*,*) 'Increase nxm nym and nzm'
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
      endif
      ii=0
      zmin=1.d300
      xmin=1.d300
      ymin=1.d300
      zmax=-1.d300
      do i=1,nz
         do j=1,ny
            do k=1,nx
               ii=ii+1
               read(15,*) xs(ii),ys(ii),zs(ii)            
               xs(ii)=xs(ii)*1.d-9
               ys(ii)=ys(ii)*1.d-9
               zs(ii)=zs(ii)*1.d-9
c               write(*,*) 'xyzn',xs(ii),ys(ii),zs(ii),ii
               xmin=min(xmin,xs(ii))
               ymin=min(ymin,ys(ii))
               zmin=min(zmin,zs(ii))
               zmax=max(zmax,zs(ii))
            enddo
         enddo
      enddo

      zmin=zmin-aretecube/2.d0
      zmax=zmax+aretecube/2.d0
      nminc=numerocouche(zmin,neps,nepsmax,zcouche)
      nmaxc=numerocouche(zmax,neps,nepsmax,zcouche)
      write(*,*) 'position of the z minimum:',zmin
      write(*,*) 'position of the z maximum:',zmax
      if (nmaxc-nminc.ge.2) then
         infostr='object inside three differents layers'
         nstop=1
         return
      elseif (nmaxc-nminc.eq.1) then
c     object inside two layers
c     compute the position of down subunit
         z=zs(1)
         zg=aretecube*nint(z/aretecube+0.5d0)-z-aretecube/2.d0
         write(*,*) z,zg
      endif
      write(*,*) '**********************************************'
      write(*,*) 'You choose an arbitrary object that you define'
      write(*,*) 'The code assumes that you define correctly' 
      write(*,*) 'a cubic lattice for your object'
      write(*,*) '**********************************************'
      nx1=nx
      ny1=ny
      nz1=nz

      ndipole=0
      do i=1,nz1
         do j=1,ny1
            do k=1,nx1
               x=xmin+dble(k-1)*aretecube
               y=ymin+dble(j-1)*aretecube
               z=zmin+dble(i-1)*aretecube+zg+aretecube/2.d0
               

               if (j.eq.1.and.k.eq.1.and.nmatf.ne.1) write(22,*) zs(1)
     $              +dble(i-1)*aretecube
               if (i.eq.1.and.k.eq.1.and.nmatf.ne.1) write(21,*) ys(1)
     $              +dble(j-1)*aretecube
               if (j.eq.1.and.i.eq.1.and.nmatf.ne.1) write(20,*) xs(1)
     $              +dble(k-1)*aretecube

               ndipole=ndipole+1
               nbsphere=nbsphere+1
               Tabdip(ndipole)=nbsphere
               Tabzn(nbsphere)=i
               xs(nbsphere)=x
               ys(nbsphere)=y
               zs(nbsphere)=z
               xswf(ndipole)=x
               yswf(ndipole)=y
               zswf(ndipole)=z
                  
               if (i.le.nz.and.j.le.ny.and.k.le.nx) then                
                  eps0=epscouche(numerocouche(zs(nbsphere),neps ,nepsmax
     $                 ,zcouche))
                  if (trope.eq.'iso') then
                     read(15,*) eps
c                     write(*,*) 'eps',aretecube,eps,eps0,k0,dddis
c     $                    ,methode,ctmp
                     call poladiffcomp(aretecube,eps,eps0,k0,dddis
     $                    ,methode,ctmp)  
                     polarisa(nbsphere,1,1)=ctmp
                     polarisa(nbsphere,2,2)=ctmp
                     polarisa(nbsphere,3,3)=ctmp
c                     write(*,*) 'pola',ctmp
                     epsilon(nbsphere,1,1)=eps
                     epsilon(nbsphere,2,2)=eps
                     epsilon(nbsphere,3,3)=eps
                  else
                     do ii=1,3
                        do jj=1,3
                           read(15,*) epsani(ii,jj)
                        enddo
                     enddo
                     call polaepstenscomp(aretecube,epsani,eps0,k0,dddis
     $                    ,methode,inv,polaeps)
                     do ii=1,3
                        do jj=1,3
                           polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                        enddo
                     enddo
                  endif
               endif
            enddo
         enddo
      enddo
      if (ndipole.gt.nmax) then
         write(*,*) 'Size of the parameter too small',ndipole
     $        ,'should be smaller that',nmax
         write(*,*) 'Increase nxm nym and nzm'
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
      endif
      nx=nx1
      ny=ny1
      nz=nz1

      if (nmatf.eq.0) then
         do i=1,nbsphere
            write(10,*) xs(i)
            write(11,*) ys(i)
            write(12,*) zs(i)
         enddo
      elseif (nmatf.eq.2) then
         
         dim(1)=nbsphere
         dim(2)=nmax
         datasetname="Dipole position x"
         call hdf5write1d(group_iddip,datasetname,xs,dim)
         datasetname="Dipole position y"
         call hdf5write1d(group_iddip,datasetname,ys,dim)
         datasetname="Dipole position z"
         call hdf5write1d(group_iddip,datasetname,zs,dim)
      endif
      close(10)
      close(11)
      close(12)
      close(15)
      close(20)
      close(21)
      close(22)

      end
