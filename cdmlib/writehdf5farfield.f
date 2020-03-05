      subroutine writehdf5farfield(Ediffkzpos,Ediffkzneg,Eimagexpos
     $     ,Eimageypos,Eimagezpos,Eimageincxpos,kxy,deltakx,deltaky,k0
     $     ,imaxk0,nfft2d,indice0,indicen,ncote,name,group_idff)
#ifdef USE_HDF5
      use HDF5

      implicit none
      integer nfft2d,imaxk0,ncote
      double precision  kxy(nfft2d),k0,indice0,indicen,deltakx,deltaky
      double complex Eimagexpos(nfft2d*nfft2d),Eimageypos(nfft2d*nfft2d)
     $     ,Eimagezpos(nfft2d*nfft2d),Ediffkzpos(nfft2d ,nfft2d,3)
     $     ,Ediffkzneg(nfft2d,nfft2d,3),Eimageincxpos(nfft2d*nfft2d)

      integer i,j,k,ii,jj
      double precision kx,ky,c,pi,quatpieps0

      
      character(40) :: name
      character(LEN=100) :: datasetname
      integer(hid_t) :: group_idff
      integer :: dim(4)

      c=299792458.d0
      pi=dacos(-1.d0)
      quatpieps0=1.d0/(c*c*1.d-7)

      do i=-imaxk0,imaxk0
         kxy(i+imaxk0+1)=dble(i)*deltakx
      enddo
      dim(1)=2*imaxk0+1
      dim(2)=nfft2d
      datasetname='kx Poynting'
      call hdf5write1d(group_idff,datasetname,kxy,dim)
      datasetname='ky Poynting'
      call hdf5write1d(group_idff,datasetname,kxy,dim)

      if (ncote.eq.0.or.ncote.eq.1) then
         k=0
         do i=-imaxk0,imaxk0
            do j=-imaxk0,imaxk0
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               k=k+1
               if (k0*k0-kx*kx-ky*ky.gt.0.d0) then                      
                  ii=imaxk0+i+1
                  jj=imaxk0+j+1
                  Eimagexpos(k)=Ediffkzpos(ii,jj,1)
                  Eimageypos(k)=Ediffkzpos(ii,jj,2)
                  Eimagezpos(k)=Ediffkzpos(ii,jj,3)

                  Eimageincxpos(k)=(cdabs(Ediffkzpos(ii,jj,1)) **2
     $                 +cdabs(Ediffkzpos(ii,jj,2))**2
     $                 +cdabs(Ediffkzpos(ii,jj,3))**2)*k0*k0*k0 *k0 *c/8
     $                 /pi*quatpieps0
               else
                  Eimagexpos(k)=0.d0
                  Eimageypos(k)=0.d0
                  Eimagezpos(k)=0.d0
               endif
            enddo
         enddo
         dim(1)=(2*imaxk0+1)**2
         dim(2)=nfft2d*nfft2d
         datasetname ='Diffracted field kz>0 x component real part'
         call hdf5write1d(group_idff,datasetname ,dreal(Eimagexpos),dim)
         datasetname ='Diffracted field kz>0 x component imaginary part'
         call hdf5write1d(group_idff,datasetname ,dimag(Eimagexpos),dim)
         datasetname ='Diffracted field kz>0 y component real part'
         call hdf5write1d(group_idff,datasetname ,dreal(Eimagexpos),dim)
         datasetname ='Diffracted field kz>0 y component imaginary part'
         call hdf5write1d(group_idff,datasetname ,dimag(Eimagexpos),dim)
         datasetname ='Diffracted field kz>0 z component real part'
         call hdf5write1d(group_idff,datasetname ,dreal(Eimagexpos),dim)
         datasetname ='Diffracted field kz>0 z component imaginary part'
         call hdf5write1d(group_idff,datasetname ,dimag(Eimagexpos),dim)
         
         datasetname='Poynting positive'
         call hdf5write1d(group_idff,datasetname ,dreal(Eimageincxpos)
     $        ,dim)
      endif
      if (ncote.eq.0.or.ncote.eq.-1) then
         k=0
         do i=-imaxk0,imaxk0
            do j=-imaxk0,imaxk0
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               k=k+1
               if (k0*k0-kx*kx-ky*ky.gt.0.d0) then                      
                  ii=imaxk0+i+1
                  jj=imaxk0+j+1
                  Eimagexpos(k)=Ediffkzneg(ii,jj,1)
                  Eimageypos(k)=Ediffkzneg(ii,jj,2)
                  Eimagezpos(k)=Ediffkzneg(ii,jj,3)

                  Eimageincxpos(k)=(cdabs(Ediffkzneg(ii,jj,1))
     $                 **2+cdabs(Ediffkzneg(ii,jj,2))**2
     $                 +cdabs(Ediffkzneg(ii,jj,3))**2)*k0*k0*k0
     $                 *k0 *c/8/pi*quatpieps0
               else
                  Eimagexpos(k)=0.d0
                  Eimageypos(k)=0.d0
                  Eimagezpos(k)=0.d0
               endif
            enddo
         enddo
         dim(1)=(2*imaxk0+1)**2
         dim(2)=nfft2d*nfft2d
         datasetname ='Diffracted field kz<0 x component real part'
         call hdf5write1d(group_idff,datasetname ,dreal(Eimagexpos),dim)
         datasetname ='Diffracted field kz<0 x component imaginary part'
         call hdf5write1d(group_idff,datasetname ,dimag(Eimagexpos),dim)
         datasetname ='Diffracted field kz<0 y component real part'
         call hdf5write1d(group_idff,datasetname ,dreal(Eimagexpos),dim)
         datasetname ='Diffracted field kz<0 y component imaginary part'
         call hdf5write1d(group_idff,datasetname ,dimag(Eimagexpos),dim)
         datasetname ='Diffracted field kz<0 z component real part'
         call hdf5write1d(group_idff,datasetname ,dreal(Eimagexpos),dim)
         datasetname ='Diffracted field kz<0 z component imaginary part'
         call hdf5write1d(group_idff,datasetname ,dimag(Eimagexpos),dim)
         
         datasetname='Poynting negative'
         call hdf5write1d(group_idff,datasetname,dreal(Eimageincxpos)
     $        ,dim)
      endif

#endif
      end
