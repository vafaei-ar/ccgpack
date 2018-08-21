subroutine sky_to_patch(inmap,npix,npatch,patch,numpa,lp)
Implicit none

integer, intent (inout) :: npatch,npix,lp,numpa
real, intent (in) :: inmap(npix)
Real ,intent (out) :: patch(numpa,lp,lp)
integer :: nside,i,j,k,face_num,nn,boi,np,nf,pow,n,xp,yp,pa
integer,allocatable :: a(:),p(:),q(:)
real  ns,base,di,x,y,ds,eql,cri,cc1,cc2

base=2.0
di=npatch     !tedade morabbaA dar zelE
pa=di**2        !tedade morabbaA dar face
!numpa=12*pa
!npix=12*nside**2
nside = int(sqrt(npix/12.))
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
	npatch=int(x/real(lp))+di*int(y/real(lp))+1  !the number of patch in face
	nn=(nf-1)*pa+npatch  !the number of patch
	xp=MOD(x,real(lp))+1.0   !position in patch
	yp=MOD(y,real(lp))+1.0

	patch(nn,xp,yp)=inmap(i)

end do  !i

end
