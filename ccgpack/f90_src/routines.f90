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


