subroutine sky_to_patch(inmap,npix,npatch,patch,numpa,lp)
Implicit none

integer(kind=8), intent (inout) :: npatch,npix,lp,numpa
real, intent (in) :: inmap(npix)
Real ,intent (out) :: patch(numpa,lp,lp)
integer(kind=8) :: nside,i,j,k,nn,boi,np,nf,pow,n,xp,yp,pa
integer(kind=8),allocatable :: a(:),p(:),q(:)
real(kind=8)  ns,base,di,x,y

base=dfloat(2)
di=npatch     !tedade morabbaA dar zelE
pa=di**2        !tedade morabbaA dar face
!numpa=12*pa
!npix=12*nside**2
nside = int(sqrt(npix/dfloat(12)))
ns=nside
np=nside**2
!lp=ns/di        !tedad pixele zelE har morabba
n=log(ns)/log(base)

allocate(a(n),p(n),q(n)) 
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
	npatch=int(x/dfloat(lp))+di*int(y/dfloat(lp))+1  !the number of patch in face
	nn=(nf-1)*pa+npatch  !the number of patch
	xp=MOD(x,dfloat(lp))+1.0   !position in patch
	yp=MOD(y,dfloat(lp))+1.0

	patch(nn,xp,yp)=inmap(i+1)

end do  !i

end


subroutine patch_to_sky(patch,inmap,npix,npatch,numpa,lp)
Implicit none

integer(kind=8), intent (inout) :: npatch,npix,lp,numpa
Real ,intent (in) :: patch(numpa,lp,lp)
real, intent (out) :: inmap(npix)
integer(kind=8) :: nside,i,j,k,np,n,nloop,pa
integer(kind=8) :: ip,jp,nface
integer(kind=8),allocatable :: a(:),p(:),q(:)
Real ,allocatable :: face(:,:)
real(kind=8) :: ns,base,di
integer, parameter:: dp=kind(0.d0)                   ! double precision

base=2.0
di=npatch     !tedade morabbaA dar zelE
pa=di**2        !tedade morabbaA dar face
!numpa=12*pa
!npix=12*nside**2
nside = int(sqrt(npix/dfloat(12)))
ns=nside
np=nside**2
!lp=ns/di        !tedad pixele zelE har morabba
n=log(ns)/log(base)
nloop=log(ns)/log(base)

allocate(a(n),p(n),q(n)) 
allocate(face(nside,nside))

do nface=1,12
!      np=mod(nface-1,pa)+1
!floor((nface-0.05)/pa)+1
  do np=1,pa
    do j=1,lp
      do i=1,lp
        face(i+int(lp*floor((np-0.05_dp)/di)),j+int(lp*(mod(np-1.,di)))) = patch(pa*(nface-1)+np,i,j)
      end do
    end do 
  end do

  do i=1,nside
    do j=1,nside
      ns=nside
      np=0
      ip=i
      jp=j
      do k=1,nloop
        np=(int((ip*2-0.01_dp)/ns)+2*int((jp*2-0.01_dp)/ns))*(ns/2)**2+np
        ns=ns/2
        ip=ip-int((ip-0.01_dp)/ns)*ns
        jp=jp-int((jp-0.01_dp)/ns)*ns
      end do !k
      np=np+(nface-1)*nside**2
      inmap(np+1)=face(j,i)
    end do !j
  end do !i
end do !face

end
