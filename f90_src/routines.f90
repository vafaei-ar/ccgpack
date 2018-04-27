
!##############   STANDARD  ##############
subroutine make_standard(nside,map,mskv,npixr,mean,var)

implicit none

integer, intent (in) :: nside
Real*8, intent (in) :: mskv
Real*8, intent (inout) :: map(nside,nside)
integer, intent (out) :: npixr
Real*8, intent (out) :: mean,var

integer  i,j

npixr=0
mean=0
do i=1,nside
	do j=1,nside
	  if (map(i,j)/=mskv) then
	    mean=mean+map(i,j)
	    npixr=npixr+1
	  end if
	end do
end do

mean=mean/npixr
map=map-mean

var=0
do i=1,nside
	do j=1,nside
	  if (map(i,j)/=mskv) then
	    var=var+map(i,j)**2
	  end if
	end do
end do
var=sqrt(var/npixr)
map=map/var

end

!##############   Threshold finder  ##############
subroutine thsh_f(rows,columns,d,lt,ut,nt,thsh)
implicit none

integer, intent (in) :: rows,columns
Real*8, intent (in) :: lt,ut
Real*8, intent (in) :: d(rows,columns)
integer, intent (out) :: nt
integer, intent (out) :: thsh(2,rows*columns)

integer  i,j

nt=0

do i=2,rows-1
  do j=2,columns-1
    if(d(i,j)>lt .and. d(i,j)<ut) then
      nt=nt+1
      thsh(1,nt)=i
      thsh(2,nt)=j
    end if
  end do 
end do

end

!##############   Peak finder  ##############
subroutine peak_f(rows,columns,d,lt,ut,np,peaks)
implicit none

integer, intent (in) :: rows,columns
Real*8, intent (in) :: lt,ut
Real*8, intent (in) :: d(rows,columns)
integer, intent (out) :: np
integer, intent (out) :: peaks(2,rows*columns)

integer  i,j
Real*8 temp

np=0

do i=2,rows-1
  do j=2,columns-1
    if(d(i,j)>lt .and. d(i,j)<ut) then

      temp=max(d(i+1,j),d(i,j+1),d(i-1,j),d(i,j-1), &
& d(i+1,j+1),d(i+1,j-1),d(i-1,j+1),d(i-1,j-1))
      if (d(i,j)>temp)  then
        np=np+1
        peaks(1,np)=i
        peaks(2,np)=j
      end if
    end if
  end do 
end do

end

!##############   Cross finder  ##############
subroutine cross_f(rows,columns,d,lt,ut,nc1,nc2,cr1,cr2)
implicit none

integer, intent (in) :: rows,columns
Real*8, intent (in) :: lt,ut
Real*8, intent (in) :: d(rows,columns)
integer, intent (out) :: nc1,nc2
integer, intent (out) :: cr1(2,rows*columns),cr2(2,rows*columns)

integer  i,j
Real*8 temp

nc1=0; nc2=0

do i=2,rows-1
  do j=2,columns-1

! First Direction
    if(d(i,j)>lt .and. d(i-1,j)<lt) then
	     nc1=nc1+1
	     cr1(1,nc1)=i
	     cr1(2,nc1)=j
		end if
! Second Direction
    if(d(i,j)>lt .and. d(i,j-1)<lt) then
		   nc2=nc2+1
		   cr2(1,nc2)=i
		   cr2(2,nc2)=j
		end if

  end do 
end do
end

!##############   Ridge finder  ##############
subroutine ridge_f(rows,columns,d,lt,ut,nr1,nr2,rdg1,rdg2)
implicit none

integer, intent (in) :: rows,columns
Real*8, intent (in) :: lt,ut
Real*8, intent (in) :: d(rows,columns)
integer, intent (out) :: nr1,nr2
integer, intent (out) :: rdg1(2,rows*columns),rdg2(2,rows*columns)

integer  i,j
Real*8 temp

nr1=0; nr2=0

do i=2,rows-1
  do j=2,columns-1

! First Direction
    if(d(i,j)>lt .and. d(i-1,j)<lt) then
			if (d(i,j+1)<d(i,j) .and. d(i,j-1)<d(i,j)) then
		    nr1=nr1+1
		    rdg1(1,nr1)=i
		    rdg1(2,nr1)=j
			end if
    end if

! Second Direction
    if(d(i,j)>lt .and. d(i,j-1)<lt) then
			if (d(i+1,j)<d(i,j) .and. d(i-1,j)<d(i,j)) then
		    nr2=nr2+1
		    rdg2(1,nr2)=i
		    rdg2(2,nr2)=j
			end if
    end if

  end do 
end do
end


!##############   P-P Correlation  ##############
subroutine   ppcf(ds,sze,x,y,n_pks,ksi,vksi,n_sh,mth,lng)
implicit none

integer, intent (in) :: sze,n_pks,n_sh,mth,lng
Real*8, intent (in) :: ds
integer, intent (in) :: x(n_pks),y(n_pks)
Real*8, intent (out) :: ksi(0:lng),vksi(0:lng)

Integer  i,j,k
Real*8  gam,ran1,ran2
integer DD(0:lng),DR(0:lng),RR(0:lng)
integer x_sh(n_sh),y_sh(n_sh)

 call init_random_seed()

ksi=0; vksi=0
dd=0; dr=0; rr=0

do i=1,n_sh
  call random_number(ran1) 
  call random_number(ran2)
  x_sh(i)=int(ran1*dble(sze)) 
  y_sh(i)=int(ran2*dble(sze))
end do

do i=1,n_pks-1
  do j=i+1,n_pks
!    if (max(abs(x(i)-x(j)),abs(y(i)-y(j)))>lng )cycle
    gam=sqrt( real( (x(i)-x(j))**2 + (y(i)-y(j))**2 ) )
    k=floor(gam/ds)
    if (k<=lng) DD(k)=DD(k)+1          
  enddo
enddo

do i=1,n_pks
  do j=1,n_sh 
!    if (max(abs(x(i)-x_sh(j)),abs(y(i)-y_sh(j)))>lng )cycle
    gam=sqrt( real( (x(i)-x_sh(j))**2 + (y(i)-y_sh(j))**2 ) )
    k=floor(gam/ds)
    if (k<=lng) DR(k)=DR(k)+1  
  end do
end do

if (mth==2 .or. mth==3) then
do i=1,n_sh 
  do j=i+1,n_sh 
    gam=sqrt( real( (x_sh(i)-x_sh(j))**2 + (y_sh(i)-y_sh(j))**2 ) )
    k=floor(gam/ds)
    if (k<=lng) RR(k)=RR(k)+1  
  end do
end do
end if

do k=1, lng 
  if(DR(k).ne.0) then
if (mth==1) ksi(k)=2.*(n_sh*1./n_pks)*DD(k)*1./DR(k)-1.
if (mth==2) ksi(k)=4.*(DD(k)*RR(k)*1./DR(k)**2)-1.
if (mth==3) ksi(k)=(DD(k)*1./RR(k))*(n_sh*1./n_pks)**2-2*(DR(k)*1./RR(k))*(n_sh*1./n_pks) + 3.

    vksi(k)=2.*(n_sh*1./n_pks)*(  ((1.-2*DD(k)/(n_pks*(n_pks-1)))/DR(k))**2 + &
	& (DD(k)*(1.-2*DR(k)/n_sh*(n_sh-1))/DR(k)**2)**2  )
  end if
end do

end

!##############   P-P Correlation  ##############
subroutine   ppcf_getrand(ds,sze,x,y,n_pks,ksi,vksi,n_sh,x_sh,y_sh,mth,lng)
implicit none

integer, intent (in) :: sze,n_pks,n_sh,mth,lng
Real*8, intent (in) :: ds
integer, intent (in) :: x(n_pks),y(n_pks)
integer, intent (in) :: x_sh(n_sh),y_sh(n_sh)
Real*8, intent (out) :: ksi(0:lng),vksi(0:lng)

Integer  i,j,k
Real*8  gam,ran1,ran2
integer DD(0:lng),DR(0:lng),RR(0:lng)
integer x_sh(n_sh),y_sh(n_sh)

 call init_random_seed()

ksi=0; vksi=0
dd=0; dr=0; rr=0

do i=1,n_pks-1
  do j=i+1,n_pks
!    if (max(abs(x(i)-x(j)),abs(y(i)-y(j)))>lng )cycle
    gam=sqrt( real( (x(i)-x(j))**2 + (y(i)-y(j))**2 ) )
    k=floor(gam/ds)
    if (k<=lng) DD(k)=DD(k)+1          
  enddo
enddo

do i=1,n_pks
  do j=1,n_sh 
!    if (max(abs(x(i)-x_sh(j)),abs(y(i)-y_sh(j)))>lng )cycle
    gam=sqrt( real( (x(i)-x_sh(j))**2 + (y(i)-y_sh(j))**2 ) )
    k=floor(gam/ds)
    if (k<=lng) DR(k)=DR(k)+1  
  end do
end do

if (mth==2 .or. mth==3) then
do i=1,n_sh 
  do j=i+1,n_sh 
    gam=sqrt( real( (x_sh(i)-x_sh(j))**2 + (y_sh(i)-y_sh(j))**2 ) )
    k=floor(gam/ds)
    if (k<=lng) RR(k)=RR(k)+1  
  end do
end do
end if

do k=1, lng 
  if(DR(k).ne.0) then
if (mth==1) ksi(k)=2.*(n_sh*1./n_pks)*DD(k)*1./DR(k)-1.
if (mth==2) ksi(k)=4.*(DD(k)*RR(k)*1./DR(k)**2)-1.
if (mth==3) ksi(k)=(DD(k)*1./RR(k))*(n_sh*1./n_pks)**2-2*(DR(k)*1./RR(k))*(n_sh*1./n_pks) + 3.

    vksi(k)=2.*(n_sh*1./n_pks)*(  ((1.-2*DD(k)/(n_pks*(n_pks-1)))/DR(k))**2 + &
	& (DD(k)*(1.-2*DR(k)/n_sh*(n_sh-1))/DR(k)**2)**2  )
  end if
end do

end

!##############   RANDOM SEED GERERATOR  ##############
subroutine init_random_seed()
  use iso_fortran_env, only: int64
implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
   + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
   + dt(3) * 24_int64 * 60 * 60 * 1000 &
   + dt(5) * 60 * 60 * 1000 &
   + dt(6) * 60 * 1000 + dt(7) * 1000 &
   + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for Real*8 work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed

!##############   TIME STAMP  ##############
subroutine timestamp ( )

!   call timestamp ( )

implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

!##############   CROSS  ##############
subroutine cross_old(p00,p0m,pm0,thsh,mskv,msrx,msry)

implicit none

real*8, intent(in) :: p00,p0m,pm0,thsh,mskv
 logical, intent(out) :: msrx,msry

msrx=(pm0<thsh) .and. (p00>thsh) .and. (p00/=mskv)
msry=(p0m<thsh) .and. (p00>thsh) .and. (p00/=mskv)

end

!##############   F-F Correlation  ##############
subroutine   ffcf(ds,sze,fl1,n_pks1,fl2,n_pks2,ksi,vksi,n_sh,mth,lng)

implicit none

integer, intent (in) :: sze,n_pks1,n_pks2,n_sh,mth,lng
Real*8, intent (in) :: ds
integer, intent (in) :: fl1(2,n_pks1)
integer, intent (in) :: fl2(2,n_pks2)
Real*8, intent (out) :: ksi(0:lng),vksi(0:lng)

Integer  i,j,k
Real*8  gam,ran1,ran2
Real*8 DD(0:lng),DR1(0:lng),DR2(0:lng),RR(0:lng)
Integer x_sh(n_sh),y_sh(n_sh)

 call init_random_seed()

ksi=0; vksi=0
dd=0; dr1=0; dr2=0; rr=0

do i=1,n_sh
  call random_number(ran1) 
  call random_number(ran2)
  x_sh(i)=int(ran1*dble(sze)) 
  y_sh(i)=int(ran2*dble(sze))
end do

do i=1,n_pks1
  do j=1,n_pks2
!    if (max(abs(x(i)-x(j)),abs(y(i)-y(j)))>lng )cycle
    gam=sqrt( real( (fl1(1,i)-fl2(1,j))**2 + (fl1(2,i)-fl2(2,j))**2 ) )
    k=floor(gam/ds)
    if (k<=lng) DD(k)=DD(k)+1          
  enddo
enddo
DD=1.*DD/(n_pks1*n_pks2)

do i=1,n_pks1
  do j=1,n_sh 
!    if (max(abs(x(i)-x_sh(j)),abs(y(i)-y_sh(j)))>lng )cycle
    gam=sqrt( real((fl1(1,i)-x_sh(j))**2 + (fl1(2,i)-y_sh(j))**2) )
    k=floor(gam/ds)
    if (k<=lng) DR1(k)=DR1(k)+1  
  end do
end do
DR1=1.*DR1/(n_pks1*n_sh)

do i=1,n_pks2
  do j=1,n_sh 
!    if (max(abs(x(i)-x_sh(j)),abs(y(i)-y_sh(j)))>lng )cycle
    gam=sqrt( real((fl2(1,i)-x_sh(j))**2 + (fl2(2,i)-y_sh(j))**2) )
    k=floor(gam/ds)
    if (k<=lng) DR2(k)=DR2(k)+1  
  end do
end do
DR2=1.*DR2/(n_pks2*n_sh)

if (mth==2 .or. mth==3) then
print*, "This part is out of service for now!"
!do i=1,n_sh 
!  do j=i+1,n_sh 
!    gam=sqrt( real( (x_sh(i)-x_sh(j))**2 + (y_sh(i)-y_sh(j))**2 ) )
!    k=floor(gam/ds)
!    if (k<=lng) RR(k)=RR(k)+1  
!  end do
!end do
end if

do k=1, lng 
  if(DR1(k)+DR2(k).ne.0) then
if (mth==1) ksi(k)=2*DD(k)/(DR1(k)+DR2(k))-1
!if (mth==2) ksi(k)=4.*(DD(k)*RR(k)*1./DR(k)**2)-1.
!if (mth==3) ksi(k)=(DD(k)*1./RR(k))*(n_sh*1./n_pks)**2-2*(DR(k)*1./RR(k))*(n_sh*1./n_pks) + 3.
vksi(k)=0
!    vksi(k)=2.*(n_sh*1./n_pks)*(  ((1.-2*DD(k)/(n_pks*(n_pks-1)))/DR(k))**2 + &
!	& (DD(k)*(1.-2*DR(k)/n_sh*(n_sh-1))/DR(k)**2)**2  )
  end if
end do

end

!##############  Correlation  ##############
subroutine cf(map,l,n_p,mskv,cor,vcor)
implicit none

integer, intent(in) :: l
real*8, intent(in) :: map(l,l)
real*8, intent(in) :: n_p
real*8, intent(in) :: mskv
real*8, intent(out) :: cor(0:l),vcor(0:l)

integer :: i,i1,j1,i2,j2,k,kp,limit
integer :: nu(0:l)
real*8 :: r1,r2,r3,gam
real*8 :: c(0:l,int(3*n_p/l))

limit=int(3*n_p/l)
cor=0
vcor=0
nu=0

i=0
Do while (i<n_p/100)
i=i+1

  call random_number(r1)
  call random_number(r2)
  call random_number(r3)
	r1=floor((l*l-1)*r1)
	r2=floor((l-1)/10.*r2)
	r3=floor((l-1)/10.*r3)
  i1=mod(int(r1),l)+1
  j1=floor(r1/l)+1
  i2=mod(int(i1+r2+1),l)+1
  j2=mod(int(j1+r3+1),l)+1

  if(map(i1,j1)/=mskv.and.map(i2,j2)/=mskv) then
    gam=sqrt( real( (i1-i2)**2 + (j1-j2)**2 ) )
    k=floor(gam)
		if (k<=l .and. nu(k)<limit) then
			nu(k)=nu(k)+1
			kp=nu(k)
			c(k,kp)=map(i1,j1)*map(i2,j2)
		end if
  end if
end do !i

Do while (i<n_p)
i=i+1

  call random_number(r1)
  call random_number(r2)
	r1=floor((l*l-1)*r1)
	r2=floor((l*l-1)*r2)
  i1=mod(int(r1),l)+1
  j1=floor(r1/l)+1
  i2=mod(int(r2),l)+1
  j2=floor(r2/l)+1

  if(map(i1,j1)/=mskv.and.map(i2,j2)/=mskv) then

    gam=sqrt( real( (i1-i2)**2 + (j1-j2)**2 ) )
    k=floor(gam)
		if (k<=l .and. nu(k)<limit) then
			nu(k)=nu(k)+1
			kp=nu(k)
			c(k,kp)=map(i1,j1)*map(i2,j2)
		end if
  end if
end do !i

! mean 
  do k=0,l
    kp=nu(k)
    do i=1,kp
      cor(k)=cor(k)+c(k,i)
    end do
    if(kp.ne.0.0) cor(k)=cor(k)/real(kp)
  end do

! error 
  do k=0,l
    kp=nu(k)
    do i=1,kp
      vcor(k)=vcor(k)+(c(k,i)-cor(k))**2
    end do
    if(kp.ne.0.0) vcor(k)=sqrt(vcor(k)/real(kp))
  end do

end 

!##############  Cross Correlation  ##############
subroutine cross_cf(map1,map2,l,n_p,mskv,cor,vcor)
implicit none

integer, intent(in) :: l
real*8, intent(in) :: map1(l,l),map2(l,l)
real*8, intent(in) :: n_p
real*8, intent(in) :: mskv
real*8, intent(out) :: cor(0:l),vcor(0:l)

integer :: i,i1,j1,i2,j2,k,kp,limit
integer :: nu(0:l)
real*8 :: r1,r2,r3,gam
real*8 :: c(0:l,int(3*n_p/l))

limit=int(3*n_p/l)
cor=0
vcor=0
nu=0

i=0
Do while (i<n_p/100)
i=i+1

11 call random_number(r1)
  call random_number(r2)
  call random_number(r3)
	r1=floor((l*l-1)*r1)
	r2=floor((l-1)/10.*r2)
	r3=floor((l-1)/10.*r3)
  i1=mod(int(r1),l)+1
  j1=floor(r1/l)+1
  i2=i1+r2+1
  j2=j1+r3+1

	if (i2>l .or. j2>l) goto 11

  if(map1(i1,j1)/=mskv.and.map2(i2,j2)/=mskv) then
    gam=sqrt( real( (i1-i2)**2 + (j1-j2)**2 ) )
    k=floor(gam)
		if (k<=l .and. nu(k)<limit) then
			nu(k)=nu(k)+1
			kp=nu(k)
			c(k,kp)=map1(i1,j1)*map2(i2,j2)
		end if
  end if
end do !i

Do while (i<n_p)
i=i+1

  call random_number(r1)
  call random_number(r2)
	r1=floor((l*l-1)*r1)
	r2=floor((l*l-1)*r2)
  i1=mod(int(r1),l)+1
  j1=floor(r1/l)+1
  i2=mod(int(r2),l)+1
  j2=floor(r2/l)+1

  if(map1(i1,j1)/=mskv.and.map2(i2,j2)/=mskv) then

    gam=sqrt( real( (i1-i2)**2 + (j1-j2)**2 ) )
    k=floor(gam)
		if (k<=l .and. nu(k)<limit) then

			nu(k)=nu(k)+1
			kp=nu(k)
			c(k,kp)=map1(i1,j1)*map2(i2,j2)
		end if

  end if
end do !i

! mean 
  do k=0,l
    kp=nu(k)
    do i=1,kp
      cor(k)=cor(k)+c(k,i)
    end do
    if(kp.ne.0.0) cor(k)=cor(k)/real(kp)
  end do

! error 
  do k=0,l
    kp=nu(k)
    do i=1,kp
      vcor(k)=vcor(k)+(c(k,i)-cor(k))**2
    end do
    if(kp.ne.0.0) vcor(k)=sqrt(vcor(k)/real(kp))
  end do

end 

!##############  VDD  ##############
subroutine vdd(map,l,w,dx,dy,vx,vy)

integer, intent(in) :: l,w
Real*8, intent(in) :: map(l,w)
Real*8, intent(out) :: dx(l,w),dy(l,w),vx,vy

integer i,j,npix,nind

do i=1,l
  nind=mod(i,l)+1
  dx(i,:)=map(i,:)-map(nind,:)
end do
do j=1,w
  nind=mod(j,w)+1
  dy(:,j)=map(:,j)-map(:,nind)
end do

mx=0
vx=0
my=0
vy=0

do i=1,l
  do j=1,w
    mx=dx(i,j)+mx
  end do
end do
mx=mx/((l)*w)

do i=1,l
  do j=1,w
    vx=(dx(i,j)-mx)**2+vx
  end do
end do
vx=sqrt(vx/((l)*w))

do i=1,l
  do j=1,w
    my=dy(i,j)+my
  end do
end do
my=my/((w)*l)

do i=1,l
  do j=1,w
    vy=(dy(i,j)-my)**2+vy
  end do
end do
vy=sqrt(vy/((w)*l))

end subroutine

!##############  Beamer  ##############
subroutine beamer(num,map,siz,fwhm,beamfun)

integer, intent (in) :: num
real*8, intent (inout) :: map(num,num)
real*8, intent (in) :: siz,fwhm
 character(256), intent (in) :: beamfun

integer lensav,ier,ldim,lenwrk,m_fft,l_fft
real*8 nu_band,diam,optele,FWHM_rad,pi,sigma_beam
real*8 delta_wavelength_x,l_k
real*8 ,allocatable ::WSAVE(:),WORK(:)
COMPLEX(8) ,allocatable :: cmap(:,:)

pi=4*atan(1.0)

nu_band=fwhm         !147e9    !Hz
diam=6    !m
optele=70     !Degree
ldim=num
m_fft=num
l_fft=num

lensav=2*(L_fft+M_fft) + INT(LOG(REAL(L_fft))/log(2.0))+ INT(LOG(REAL(M_fft))/log(2.0)) + 8
lenwrk=2*L_fft*M_fft+10

allocate(cmap(ldim,m_fft),WSAVE(lensav),WORK(lenwrk))

delta_wavelength_x=2*pi/(siz*pi/180)

if (beamfun=="gaussian") then
  FWHM_rad=pi*FWHM/(60.0*180.0)
  sigma_beam=FWHM_rad/sqrt(8*log(2.0))
else if (beamfun=="airy") then
  sigma_beam=360.*nu_band*diam/(optele*299792458.)
end if

 cmap = map
 call cfft2i ( l_fft,m_fft, wsave, lensav, ier )
 call cfft2f ( ldim, l_fft, m_fft, cmap, wsave, lensav, work, lenwrk, ier )

!********************** for me ********************
do i=1,num
do j=1,num

! this true when we use L as radian
l_k=sqrt((i*1.0)**2+(j*1.0)**2)*delta_wavelength_x

if (beamfun=="gaussian") then
  beam= exp(-(sigma_beam**2)*l_k*(l_k+1)/2.)
else if (beamfun=="airy") then
  if (l_k<sigma_beam) then
    beam=(2./pi)*(acos(l_k/sigma_beam)-(l_k/sigma_beam) &
& *sqrt(1-(l_k/sigma_beam)**2))
  else 
    beam=0
  end if
end if
 cmap(i,j)=cmap(i,j)*beam    !*(size/2*pi)**2)
enddo
enddo

 call cfft2b ( ldim, l_fft, m_fft, cmap, wsave, lensav, work, lenwrk, ier )

map = real(cmap)    

end

