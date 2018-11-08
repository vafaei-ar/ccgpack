program map_converter

  USE healpix_types
  USE head_fits,  ONLY : write_minimal_header
  USE fitstools, ONLY : write_bintab,getsize_fits,input_map,output_map
  USE pix_tools, ONLY : convert_inplace,ring2nest,nest2ring, npix2nside,convert_ring2nest,xy2pix_nest,pix2ang_nest,pix2ang_ring !  use pix_tools,  ONLY : convert_inplace, ring2nest, nest2ring, npix2nside,pix2ang_ring,convert_ring2nest
  Use paramfile_io, ONLY : paramfile_handle, parse_init, parse_string  
  USE extension,  ONLY : getArgument, nArguments  
  use misc_utils, ONLY : assert_alloc, fatal_error

  Implicit none

 INTEGER(I4B) nloop,ip,jp,b,nface,ops
  integer(i4b) :: nside,npix,npixr,nmin,nmax,npp,dow1,equator,min_peak,directional,biask,wac,nh,sh,thetai,phii,nthetai,nphii
  integer(i4b) :: i,j,k,face_num,nn,boi,np,nf,pow,aa,n,nth,ktsh,xp,yp,pa,nupo,lk,spo,lp

  real*8 thetae,phie,gam,theta,phi
  integer,allocatable :: a(:),p(:),q(:)
 
  real(sp), allocatable :: inmap(:,:)
  integer(i4b),allocatable :: peak(:)    !pointer 

  integer(i4b) :: ordering,order_type
  integer(i4b) :: nlheader
  integer(i4b) :: nmaps = 1

  character(len=80), dimension(1:120) :: header
  character(len=filenamelen) :: parafile,infile, &
                                outfile
  character(len=40) :: coordsys,sktsh,str,salpha,sthetad,sphid,snam
  character*25 fnam,cface,nam,scheme,comnam,opss
  integer(i4b) :: status

  type(paramfile_handle) :: handle

  CHARACTER(LEN=*), PARAMETER :: code    = 'map_converter'
  character(len=*), parameter :: VERSION = HEALPIX_VERSION

  Integer(i4b)  numpa,nump(2),pole,n_sh,min_k,max_k,l,kk,siz,xc,yc,icent,j1,j2,npatch

  real  var,dummy,ns,base,alpha,di,x,y,ds,eql1,eql2,cc2
  Integer(i4b) ,Allocatable ::  npk(:),peak_c(:),peak_list(:),north_list(:),south_list(:)
  Real ,Allocatable :: patch(:,:,:),pcri(:),face(:,:)
 
  !--------------------------------------------------------------------------
100 FORMAT(A)

!  PRINT *,'This program is writed for Calculatin Bias Factor of fits CMB maps.'

  if (nArguments() == 3) then  !Ring/Nested together
     call getArgument(1,parafile)
     infile=trim(adjustl(parafile))
     call getArgument(2,parafile)
     outfile=trim(adjustl(parafile))
     call getArgument(3,parafile)
     scheme=trim(adjustl(parafile))
  else if (nArguments() == 4) then  !Ring/Nested to Patche
     call getArgument(1,parafile)
     infile=trim(adjustl(parafile))
     call getArgument(2,parafile)
     outfile=trim(adjustl(parafile))
     call getArgument(3,parafile)
     scheme=trim(adjustl(parafile))
     call getArgument(4,parafile)
     read(parafile,*) npatch
  else if (nArguments() == 5 .or. nArguments() == 6) then
     call getArgument(1,parafile)
     infile=trim(adjustl(parafile))
     call getArgument(2,parafile)
     outfile=trim(adjustl(parafile))
     call getArgument(3,parafile)
     scheme=trim(adjustl(parafile))
     call getArgument(4,parafile)
     read(parafile,*) npatch
     call getArgument(5,parafile)
     comnam=trim(adjustl(parafile))
!     call getArgument(6,parafile)
!     opss=trim(adjustl(parafile))
  else
     Print*, 'Please, argument!'
     Print*, 'Input;Output;out scheme;n(# of maps=12n);common name;ops'
     STOP
  end if

!print *, 'Size of patch (Degree)?'
alpha=60
!print *, 'Enter 1 if you want separated patches sequences and 2 for dual?'
spo=1
!print *, 'Dismiss equator (10 degree)?'
equator=0

          cc2=-11     
          base=2.0

          di=npatch     !tedade morabbaA dar zelE
          pa=di**2        !tedade morabbaA dar face
          numpa=12*pa

!  handle=parse_init (parafile)

!     --- gets the file name for the map ---

!  infile=parse_string(handle,'infile',default='.fits', &
!       descr='Input file name:', &
!       filestatus='old')
!  PRINT *,''

  !     --- finds out the pixel number of the map and its ordering --- 

if (nArguments() == 3 .or. nArguments() == 4) then
  npix = getsize_fits(infile ,ordering=ordering, nside=nside, coordsys = coordsys)
  order_type=ordering

  if (nside.eq.0) then
     print*,'Keyword NSIDE not found in FITS header!'
  endif
  if (nside.ne.npix2nside(npix)) then
     print*,'FITS header keyword NSIDE does not correspond'
     print*,'to the size of the map!'
  endif

  !     --- check ordering scheme ---

  if ((order_type.ne.1).and.(order_type.ne.2)) then
     PRINT*,'The ordering scheme of the map must be RING or NESTED.'
     PRINT*,'No ordering specification is given in the FITS-header!'
  endif

  allocate(inmap(0:npix-1,1),stat=status)
  call assert_alloc(status,code,'map')

  !    --- Reading ---

  call input_map(infile,inmap,npix,nmaps,0.0,header)

  inmap(0:npix-1,1) = inmap(0:npix-1,1)

else if (nArguments() == 5 .or. nArguments() == 6) then

          nmaps=1

  do nface=1,12

      if (nface==1) then
        fnam=trim(adjustl(comnam))//'1'
        open (1,file=fnam)
        lp=0

        do
          read (1,*,end=101) dummy
          lp=lp+1
        end do
101     Rewind(1)  

        nside=di*lp
        npix=12*nside**2
        ns=nside
        nloop=log(ns)/log(base)
!        lp=ns/di        !tedad pixele zelE har morabba
        n=log(ns)/log(base)

        allocate(inmap(0:npix-1,1),stat=status)
        call assert_alloc(status,code,'map')

        allocate(face(nside,nside),pcri(numpa))
      end if



!      np=mod(nface-1,pa)+1
!floor((nface-0.05)/pa)+1
    do np=1,pa
      write(cface,'(i2)') pa*(nface-1)+np
!      fnam=trim(adjustl(nam))//trim(adjustl(cface))
      fnam=trim(adjustl(comnam))//trim(adjustl(cface))
      open (1,file=fnam)
      do j=1,lp
        read (1,*) (face(i+int(lp*floor((np-0.05)/di)),j+int(lp*(mod(np-1.,di)))),i=1,lp)
      end do !i
      close(1)
    end do

      do i=1,nside
        do j=1,nside
          ns=nside
          np=0
          ip=i
          jp=j
          do k=1,nloop
            np=(int((ip*2-0.01)/ns)+2*int((jp*2-0.01)/ns))*(ns/2)**2+np
            ns=ns/2
            ip=ip-int((ip-0.01)/ns)*ns
            jp=jp-int((jp-0.01)/ns)*ns
          end do !k
          np=np+(nface-1)*nside**2
          inmap(np,1)=face(i,j)
        end do !j
      end do !i
  end do !face

end if

     Print*, 'Do you want to setup mask? 1 if yes.'
     read*, equator

! Masking
  if(equator==1) then

     Print*, 'Theta1 (Mask more that it)?'
     read*, eql1
     Print*, 'Theta2 (Mask less that it)?'
     read*, eql2
     Print*, 'Mask value?'
     read*, cc2

    do i=0,npix-1

if (nArguments() == 3 .or. nArguments() == 4) then

  if (order_type==1) then
      call pix2ang_ring(nside,i,thetae,phie)
  else if (order_type==2) then
      call pix2ang_nest(nside,i,thetae,phie)
  end if

else if (nArguments() == 5 .or. nArguments() == 6) then
      call pix2ang_nest(nside,i,thetae,phie)
end if

      thetae=thetae*180.0/pi
      if(thetae>eql1 .and. thetae<eql2) then
        inmap(i,1)=cc2
        npixr=npixr-1
      end if
    end do
  end if


if (scheme=='nested') then

  if (order_type==2) then
    print *, 'Input map is already NESTED!'
    goto 102
!    STOP
  end if

  if (nArguments() < 5 ) then
     write(*,*) code//'> Converting to NESTED numbering scheme.'
     call convert_inplace(ring2nest,inmap(:,1))
  end if
	 order_type=2	   
         ordering=order_type

102  call write_minimal_header(header,'map',&
       nside=nside , order = ordering , &
       creator = code, version = version)
  nlheader=SIZE(header)
  call write_bintab(inmap, npix, nmaps, header, nlheader, outfile)

else if (scheme=='ring') then

  if (order_type==1) then
    print *, 'Input map is already RING!'
    goto 103
!    STOP
  end if

     write(*,*) code//'> Converting to RING numbering scheme.'
     call convert_inplace(nest2ring,inmap(:,1))
	 order_type=1	    
         ordering=order_type

103  call write_minimal_header(header,'map',&
       nside=nside , order = ordering , &
       creator = code, version = version)

  nlheader=SIZE(header)
  call write_bintab(inmap, npix, nmaps, header, nlheader, outfile)

else if (scheme=='patch') then

  if (order_type==1) then
     write(*,*) code//'> Converting to NESTED numbering scheme.'
     call convert_inplace(ring2nest,inmap(:,1))
	 order_type=2	   
         ordering=order_type
  end if

          ns=nside
          lp=ns/di        !tedad pixele zelE har morabba
          n=log(ns)/log(base)
          np=nside**2


allocate(a(n),p(n),q(n),patch(numpa,lp,lp),pcri(numpa)) 
  do i=0,npix-1

          boi=i
          x=0 ;y=0

          nf=int(boi/np)+1
          boi=boi-(nf-1)*np
          do j=1,n
            pow=np/(4**(j-1))
            a(j)=int((boi*4)/pow)
            boi=boi-(a(j)*pow/4)
            if (a(j)==0) then
              p(j)=1
              q(j)=1
            end if
            if (a(j)==1) then
              p(j)=1
              q(j)=2
            end if
            if (a(j)==2) then
              p(j)=2
              q(j)=1
            end if
            if (a(j)==3) then
              p(j)=2
              q(j)=2
            end if
            x=x+(nside/2**j)*(p(j)-1)
            y=y+(nside/2**j)*(q(j)-1)

          end do   ! j
!           call pix2xy_nest(nside,x,y,nf,i)
!		  x=x
!		  y=y
	  ! x & y are positions in face 'nf'	  

          npatch=int(x/real(lp))+di*int(y/real(lp))+1  !the number of patch in face
          nn=(nf-1)*pa+npatch  !the number of patch
          xp=MOD(x,real(lp))+1.0   !position in patch
          yp=MOD(y,real(lp))+1.0

     if (xp<=real(lp).and.yp<=real(lp))  patch(nn,xp,yp)=inmap(i,1)
     if (patch(nn,xp,yp)==cc2)  pcri(nn)=pcri(nn)+1 


  end do  !i

if (spo==1) then
  do nn=1,numpa

    write(sthetad,'(i4)') nn
    sthetad=trim(adjustl(outfile))//trim(adjustl(sthetad))
  open(21,file=sthetad)
    do i=1,lp
      write(21,*) (patch(nn,i,j),j=1,lp)  
    end do
  close(21)

  end do
end if

if (spo==2) then
  do nn=1,numpa/2

    write(sthetad,'(i4)') nn
    sthetad='f_'//trim(adjustl(sthetad))
  open(21,file=sthetad)
    do i=1,lp
      write(21,*) (patch(nn,i,j),j=1,lp)  
    end do
  close(21)

      nf=int((real(nn)-0.001)/real(pa))
      npatch=mod(nn-1,pa)+1
      ops=pa*(11-nf+(-1)**real(nf+1))+npatch

    write(sthetad,'(i4)') nn
    sthetad='s_'//trim(adjustl(sthetad))
  open(21,file=sthetad)
    do i=1,lp
      write(21,*) (patch(nn,i,j),j=1,lp)  
    end do
  close(21)

  end do
end if

else if (scheme=='angle') then

	  do i=0,npix-1
  if (order_type==1) then
	    call pix2ang_ring(nside,i,theta,phi)
  else if (order_type==2) then
	    call pix2ang_nest(nside,i,theta,phi)
  end if
	    write (1,*) theta,phi,inmap(i,1)
	  end do

end if

if (.false.) then

      do nn=1,npp
	j=peak_list(nn)

          nf=int((real(j)-0.001)/real(pa))
          npatch=mod(j-1,pa)+1
          xc=int((mod(real(npatch)-0.5,real(di)))*real(lp))
          yc=int((int(real(npatch)/real(di)-0.05)+0.5)*real(lp))
!          call xy2pix_nest(nside,xc,yc,nf,icent)
!          call pix2ang_nest(nside,icent,theta(j),phi(j))
!          theta(j)=int(theta(j)*180/pi)
!          phi(j)=int(phi(j)*180/pi) 
      
       end do   

end if
end program map_converter

