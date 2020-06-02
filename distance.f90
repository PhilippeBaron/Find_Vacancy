!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine SEmin1(At21,x21,y21,z21, id1,Rid)

COMMON / RAD / Radx, Rady, Radz, Radm, PRATIOx, PRATIOy, PRATIOz
COMMON / SPH_MESH / Xe, Ye, r, Ni, idref

!!===========================================
! The main purpose of this subroutine is
! calculate the minimum coordinate (Superellipse)
! of the particle (2D)
! Outputs:
! Rid
!!===========================================
double precision	At21(3,3)
double precision	x21, y21, z21

double precision	Radx, Rady, Radz, Radm
double precision	PRATIOx, PRATIOy, PRATIOz

double precision	r
integer				Ni, idref
double precision	Xe(600), Ye(600)

double precision, parameter :: pi = 3.141592653589793

integer				idx(8)
integer				j, k
double precision	Xe21, Ye21
double precision	Re21
double precision	Rid
integer				ind, id1
!!===========================================

!!-------------------------------------------
!! part. 2 respect to part. 1 coordinates
!!-------------------------------------------
Rid = 1000000.0

idx = (/Ni + 1, Ni + 1 + idref, 2*Ni + 1, 2*Ni + 1 + idref, &
	3*Ni + 1, 3*Ni + 1 + idref, 4*Ni + 1, 4*Ni + 1 + idref/)

ind = 0
do j = 1,8
	ind = ind + 1
	k = idx(j)
	!! part 2 relative to part 1
	Xe21 = x21 + At21(1,1)*Xe(k) + At21(2,1)*Ye(k) !! + At21(3,1)*Ze(k)
	Ye21 = y21 + At21(1,2)*Xe(k) + At21(2,2)*Ye(k) !! + At21(3,2)*Ze(k)
	!! Ze21 = z21 + At21(1,3)*Xe(k) + At21(2,3)*Ye(k) + At21(3,3)*Ze(k)

	Re21 = ( abs(Xe21/PRATIOx)**r + abs(Ye21/PRATIOy)**r )**(1.0/r)

	if(Re21 .LT. Rid )then
		!! position of the minimum (1st mesh)
		id1 = ind
		Rid = Re21
	endif

end do
!!-------------------------------------------

!!===========================================
end subroutine SEmin1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SEmin2(At21,x21,y21,z21,id1, id2,Rid)

COMMON / RAD / Radx, Rady, Radz, Radm, PRATIOx, PRATIOy, PRATIOz
COMMON / SPH_MESH / Xe, Ye, r, Ni, idref

!!===========================================
! The main purpose of this subroutine is
! calculate the minimum coordinate (Superellipse)
! of the particle (2D)
! Outputs:
! Rid
!!===========================================
double precision	At21(3,3)
double precision	x21, y21, z21

double precision	Radx, Rady, Radz, Radm
double precision	PRATIOx, PRATIOy, PRATIOz

double precision	r
integer				Ni, idref
double precision	Xe(600), Ye(600)

double precision, parameter :: pi = 3.141592653589793

integer				idx(10)
double precision	Xe21, Ye21
double precision	Re21
double precision	Rid
integer				del, ind, id1, id2
integer				i0, i1, i2
!!===========================================

!!-------------------------------------------
!! part. 2 respect to part. 1 coordinates
!!-------------------------------------------
Rid = 1000000.0

del = 10
idx = (/Ni + 1 - idref, Ni + 1, Ni + 1 + idref, 2*Ni + 1, 3*Ni + 1 - idref, &
	3*Ni + 1, 3*Ni + 1 + idref, 4*Ni + 1, 5*Ni + 1 - idref, 5*Ni + 1/)
i1 = idx(id1)
i2 = idx(id1 + 2)

ind = 0
do k = i1, i2, del
	ind = ind + 1
	!! part 2 relative to part 1
	Xe21 = x21 + At21(1,1)*Xe(k) + At21(2,1)*Ye(k) !! + At21(3,1)*Ze(k)
	Ye21 = y21 + At21(1,2)*Xe(k) + At21(2,2)*Ye(k) !! + At21(3,2)*Ze(k)
	!! Ze21 = z21 + At21(1,3)*Xe(k) + At21(2,3)*Ye(k) + At21(3,3)*Ze(k)

	Re21 = ( abs(Xe21/PRATIOx)**r + abs(Ye21/PRATIOy)**r )**(1.0/r)

	if(Re21 .LT. Rid )then
		!! position of the minimum 
		id2 = ind
		Rid = Re21
	endif

end do
!!-------------------------------------------

!!===========================================
end subroutine SEmin2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SEmin3(At21,x21,y21,z21,id1,id2, id3,Rid)

COMMON / RAD / Radx, Rady, Radz, Radm, PRATIOx, PRATIOy, PRATIOz
COMMON / SPH_MESH / Xe, Ye, r, Ni, idref

!!===========================================
! The main purpose of this subroutine is
! calculate the minimum coordinate (Superellipse)
! of the particle (2D)
! Outputs:
! Rid
!!===========================================
double precision	At21(3,3)
double precision	x21, y21, z21

double precision	Radx, Rady, Radz, Radm
double precision	PRATIOx, PRATIOy, PRATIOz

double precision	r
integer				Ni, idref
double precision	Xe(600), Ye(600)

double precision, parameter :: pi = 3.141592653589793

integer				idx(10)
double precision	Xe21, Ye21
double precision	Re21
double precision	Rid
integer				del, ind, id1, id2, id3
integer				i0, i1, i2
!!===========================================

!!-------------------------------------------
!! part. 2 respect to part. 1 coordinates
!!-------------------------------------------
Rid = 1000000.0

del = 10
idx = (/Ni + 1 - idref, Ni + 1, Ni + 1 + idref, 2*Ni + 1, 3*Ni + 1 - idref, &
	3*Ni + 1, 3*Ni + 1 + idref, 4*Ni + 1, 5*Ni + 1 - idref, 5*Ni + 1/)

i1 = idx(id1) + (id2 - 1)*del

ind = 0
do k = i1 - del, i1 + del
	ind = ind + 1
	!! part 2 relative to part 1
	Xe21 = x21 + At21(1,1)*Xe(k) + At21(2,1)*Ye(k) !! + At21(3,1)*Ze(k)
	Ye21 = y21 + At21(1,2)*Xe(k) + At21(2,2)*Ye(k) !! + At21(3,2)*Ze(k)
	!! Ze21 = z21 + At21(1,3)*Xe(k) + At21(2,3)*Ye(k) + At21(3,3)*Ze(k)

	Re21 = ( abs(Xe21/PRATIOx)**r + abs(Ye21/PRATIOy)**r )**(1.0/r)

	if(Re21 .LT. Rid )then
		!! position of the minimum
		id3 = ind
		Rid = Re21
	endif

end do
!!-------------------------------------------

!!===========================================
end subroutine SEmin3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CellList(Cid,idc,ind,cell,pcell)

COMMON / XYZ / RX, RY, RZ
COMMON / BOXLEN / xbox, ybox, rcut
COMMON / NUM_PAR / N
COMMON / CELLS / delx, dely, rcell, Nc, Nppc, Nref
!!===========================================
! The main purpose of this subroutine is
! generate the cell list
! Outputs:
! Cid,idc
!!===========================================
integer				N
integer,parameter  ::  Nmax = 1000
double precision	RX(Nmax), RY(Nmax), RZ(Nmax)
integer				Nc, Nppc, Nref
double precision	rcell
double precision	delx, dely
double precision	xbox, ybox, rcut
integer				idc(Nc,Nc)
integer				Cid(Nc*Nc,9)

integer				i, j, k1

double precision	Xc(Nc), Yc(Nc)
integer				ind(Nc*Nc)
integer				cell(Nc*Nc,Nppc) !cell(Nc*Nc,Nppc)
integer				pcell(N)
integer				cx, cy, cp

double precision, parameter :: pi = 3.141592653589793
!!===========================================

!!-------------------------------------------
!! Cell list
!!-------------------------------------------
do i = 1, Nc
	!! Cell position (center)
	Xc(i) = -xbox/2.0 + (i - 1)*delx + delx/2.0
	Yc(i) = -ybox/2.0 + (i - 1)*dely + dely/2.0
enddo

k1 = 0
do i = 1, Nc
	do j = 1, Nc
		k1 = k1 + 1
		!! Cell index
		idc(i,j) =  k1
	enddo
enddo

k1 = 0
do j = 1, Nc
	do i = 1, Nc
		k1 = k1 + 1

		if( i .EQ. 1 .AND. (j .GT. 1 .AND. j .LT. Nc) )then

			Cid(k1,:) = (/ idc(j,i), idc(j-1,i+1), idc(j,i+1), idc(j+1,i+1), &
				idc(j-1,i), idc(j+1,i), idc(j-1,Nc), idc(j,Nc), idc(j+1,Nc) /)

		elseif( j .EQ. 1 .AND. (i .GT. 1 .AND. i .LT. Nc) )then

			Cid(k1,:) = (/ idc(j,i), idc(Nc,i+1), idc(j,i+1), idc(j+1,i+1), &
				idc(Nc,i), idc(j+1,i), idc(Nc,i-1), idc(j,i-1), idc(j+1,i-1) /)

		elseif(i .EQ. 1 .AND. j .EQ. 1)then

			Cid(k1,:) = (/ idc(j,i), idc(Nc,i+1), idc(j,i+1), idc(j+1,i+1), &
				idc(Nc,i), idc(j+1,i), idc(Nc,Nc), idc(j,Nc), idc(j+1,Nc) /)

		elseif(i .EQ. 1 .AND. j .EQ. Nc)then

			Cid(k1,:) = (/ idc(j,i), idc(j-1,i+1), idc(j,i+1), idc(1,i+1), &
				idc(j-1,i), idc(1,i), idc(j-1,Nc), idc(j,Nc), idc(1,Nc) /)

		elseif(i .EQ. Nc .AND. j .EQ. 1)then

			Cid(k1,:) = (/ idc(j,i), idc(Nc,1), idc(j,1), idc(j+1,1), idc(Nc,i), &
				idc(j+1,i), idc(Nc,i-1), idc(j,i-1), idc(j+1,i-1) /)

		elseif( i .EQ. Nc .AND. (j .GT. 1 .AND. j .LT. Nc) ) then

			Cid(k1,:) = (/ idc(j,i), idc(j-1,1), idc(j,1), idc(j+1,1), idc(j-1,i), &
				idc(j+1,i), idc(j-1,i-1), idc(j,i-1), idc(j+1,i-1) /)

		elseif( j .EQ. Nc .AND. (i .GT. 1 .AND. i .LT. Nc) ) then

			Cid(k1,:) = (/ idc(j,i), idc(j-1,i+1), idc(j,i+1), idc(1,i+1), &
				idc(j-1,i), idc(1,i), idc(j-1,i-1), idc(j,i-1), idc(1,i-1) /)

		elseif(i .EQ. Nc .AND. j .EQ. Nc) then

			Cid(k1,:) = (/ idc(j,i), idc(j-1,1), idc(j,1), idc(1,1), idc(j-1,i), &
				idc(1,i), idc(j-1,i-1), idc(j,i-1), idc(1,i-1) /)

		else

			Cid(k1,:) = (/ idc(j,i), idc(j-1,i+1), idc(j,i+1), idc(j+1,i+1), &
				idc(j-1,i), idc(j+1,i), idc(j-1,i-1), idc(j,i-1), idc(j+1,i-1) /)

		endif

	enddo
enddo
!!-------------------------------------------
!! Cell position of every particle (Cell list)
ind(:) = 0
cell(:,:) = Nref
pcell(:) = 0

do i = 1, N
	cx = floor((xbox/2.0 + RX(i))/delx)		!! cell position along x
	cy = floor((ybox/2.0 + RY(i))/dely) + 1	!! cell position along y
	cp = cx*Nc + cy							!! cell position of the part. i
	ind(cp) = ind(cp) + 1	!! number of parts. in the cell
    WRITE(*,*) i, " ", cx, " ", cy, " ", cp, " ", ind(cp)
	cell(cp,ind(cp)) = i	!! list of particles in the cell
	pcell(i) =  cp			!! cell position of the parts
enddo
!!-------------------------------------------

!!===========================================
end subroutine CellList
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SortVect(N, X, Xs)
!!===========================================
! The main purpose of this subroutine is
! sort the elements of a vector X
! and report the sorted part. index
! Outputs:
! Xs
!!===========================================
integer				N
integer				X(N), Xs(N)

integer				i, j
integer				myPM
double precision	m

double precision, parameter :: pi = 3.141592653589793
!!===========================================
Xs = X

do i = 1, N
	do j = 1, N
		if(Xs(i) .LT. Xs(j)) then
			m = Xs(j)
			Xs(j) = Xs(i)
			Xs(i) = m
		endif
	enddo
enddo
!!===========================================

end subroutine SortVect
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





