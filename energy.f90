!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE Energy (rx1,ry1,rz1,chi1,eta1,xi1,zeta1,I, V)

COMMON / XYZ / RX, RY, RZ
COMMON / ORIENT / chi, zeta
COMMON / NUM_PAR / N
COMMON / BOXLEN / xbox, ybox, rcut
COMMON / RAD / Radx, Rady, Radz, Radm, PRATIOx, PRATIOy, PRATIOz
COMMON / SPH_MESH / Xe, Ye, r, Ni, idref
COMMON / ELECT / Bpw, kappa
COMMON / FIELD / lambdaV, dg, de, fcmr, fcmi, F0, Lx, Ly, Lz

!!*******************************************************************
!!** RETURNS THE POTENTIAL ENERGY OF ATOM I WITH ALL OTHER ATOMS.  **
!!**                                                               **
!!** PRINCIPAL VARIABLES:                                          **
!!**                                                               **
!!** INTEGER I					THE ATOM OF INTEREST               **
!!** INTEGER N					THE NUMBER OF ATOMS                **
!!** REAL    RX(N),RY(N),RZ(N)	THE ATOM POSITIONS                 **
!!** REAL    rx1,ry1,rz1		THE COORDINATES OF ATOM I          **
!!** REAL    ORX(N),ORY(N),ORZ(N) ORIENTATION                      **
!!** REAL    ORXI,ORYI,ORZI		ORIENTATION OF ATOM I              **
!!** REAL    V					THE POTENTIAL ENERGY OF ATOM I     **
!!** REAL    W					THE VIRIAL OF ATOM I               **
!!**                                                               **
!!** USAGE:                                                        **
!!**                                                               **
!!** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
!!** DURING A TRIAL MOVE OF ATOM I. IT IS CALLED BEFORE AND        **
!!** AFTER THE RANDOM DISPLACEMENT OF I.                           **
!!*******************************************************************

integer				N
integer,parameter::	Nmax = 1000 !, NBINSMAX = 1000
double precision	RX(Nmax), RY(Nmax), RZ(Nmax)
double precision	chi(Nmax), zeta(Nmax)
double precision	V
double precision	POTREF

double precision	xbox, ybox, rcut

double precision	Radx, Rady, Radz, Radm
double precision	PRATIOx, PRATIOy, PRATIOz

double precision	r
integer				Ni, idref
double precision	Xe(600), Ye(600)
double precision	Rid1, Rid2, Rid3
integer				id1, id2, id3
double precision	Rjd1, Rjd2, Rjd3
integer				jd1, jd2, jd3

integer				I, J, k

double precision	rx1, ry1, rz1
double precision	chi1, eta1, xi1, zeta1
double precision	At1(3,3)

double precision	rx2, ry2, rz2
double precision	chi2, eta2, xi2, zeta2
double precision	At2(3,3)

double precision	rx21, ry21, rz21, r21

double precision	chi21, eta21, xi21, zeta21
double precision	At21(3,3)
double precision	x21, y21, z21, re21
double precision	hx21, hy21, hz21

double precision	chi12, eta12, xi12, zeta12
double precision	At12(3,3)
double precision	x12, y12, z12
double precision	hx12, hy12, hz12

double precision	Bpw, kappa

double precision	lambdaV, dg, de, fcmr(3)
double precision	fcmi(3), F0, Lx, Ly, Lz
double precision	Vdep
double precision	Epx1, Epy1, Epz1
double precision	Vdd

double precision, parameter :: pi = 3.141592653589793
!!******************************************************************

V  = 0.d0
POTREF = 0.d0

!!===========================================
!! Particle (1) I - position and orientation
!! rx1, ry1, rz1, chi1, eta1, xi1, zeta1
!! Transformation matrix
call Atransf(chi1,eta1,xi1,zeta1,At1)
!!===========================================
!! Time averaged particle-field interaction
!! Udep of particle (1) I
call Udep(rx1,ry1,rz1,At1, Vdep,Epx1,Epy1,Epz1)
!write(*,*) Vdep
!!write(*,*) 'Ep1 = ', Epx1, Epy1, Epz1
!!===========================================

V = V + 0.42*Vdep - POTREF

!!******************************************************************
!!** LOOP OVER ALL MOLECULES EXCEPT I  **

DO 100 J = 1, N
	!!===========================================
	IF ( I .NE. J ) THEN
	!!===========================================
		!! Particle I -
		!! rx1, ry1, rz1
		!! chi1, eta1, xi1, zeta1
		!!===========================================
		!!===========================================
		!! J-Particle (2) position and orientation
		rx2 = RX(J)
		ry2 = RY(J)
		rz2 = RZ(J)
		chi2 = chi(J)
		eta2 = 0.d0
		xi2 = 0.d0
		zeta2 = zeta(J)
		!!-------------------------------------------
		!! Transformation matrix
		!! call Atransf(chi2,eta2,xi2,zeta2,At2)

		!!-------------------------------------------
		!! Periodic boundary condition (PBC)
		!!-------------------------------------------
		!! position of part.2 relative to part.1
		rx21 = rx2 - rx1
		ry21 = ry2 - ry1
		rz21 = 0.d0 !! rz2 - rz1
		!! x-direction (PBC)
		rx21  = rx21 - ANINT( rx21/xbox )*xbox
		!! y-direction
		ry21  = ry21 - ANINT( ry21/ybox )*ybox
		!! particle position (PBC)
		rx2 = rx1 + rx21
		ry2 = ry1 + ry21
		!! Distance between particle centers
		!! (Relative to part.1 - lab coord)
		r21 = sqrt( rx21**2.d0 + ry21**2.d0 ) !! + rz21**2.d0 )
		!!-------------------------------------------

		!!-------------------------------------------
		!! Transformation matrix
		!! call Atransf(chi1,eta1,xi1,zeta1,At1)
		call Atransf(chi2,eta2,xi2,zeta2,At2)
		!!-------------------------------------------
		!! part. 2 respect to part. 1 coordinates
		!!-------------------------------------------
		!! quaternions of part. 2 respect to part. 1
		call quat21(chi2,eta2,xi2,zeta2,chi1,eta1,xi1,zeta1,chi21,eta21,xi21,zeta21)
		!! Transformation matrix
		call Atransf(chi21,eta21,xi21,zeta21,At21)
		!! write(*,*) 'q21', chi21,eta21,xi21,zeta21
		!!-------------------------------------------
		!! Relative position of part. 2 relative to part. 1
		x21 = At1(1,1)*rx21 + At1(1,2)*ry21 ! + At1(1,3)*rz21
		y21 = At1(2,1)*rx21 + At1(2,2)*ry21 ! + At1(2,3)*rz21
		z21 = 0.d0 !! At1(3,1)*rx21 + At1(3,2)*ry21 + At1(3,3)*rz21
		!! write(*,*) 'X21 ->', x21, y21, z21
		!!-------------------------------------------


		!!-------------------------------------------
		!!-------------------------------------------
		!!-------------------------------------------
		!!-------------------------------------------
		!! Overlap
		!!-------------------------------------------
		!! First approximation
		if( r21 .GT. rcut )then !! rcut is for hard particle interaction
			V = V + 0.d0
			GOTO 400
		else
			!! Second approximation: Distance to a rectangular part.
			!! Transformation matrix
			!! call Atransf(chi1,eta1,xi1,zeta1,At1)
			!! call Atransf(chi2,eta2,xi2,zeta2,At2)
			!!:::::::::::::::::::::::::::::::::::::::::::
			!! part. 2 respect to part. 1 coordinates
			!!:::::::::::::::::::::::::::::::::::::::::::
			!! quaternions of part. 2 respect to part. 1
			!! call quat21(chi2,eta2,xi2,zeta2,chi1,eta1,xi1,zeta1,chi21,eta21,xi21,zeta21)
			!! Transformation matrix
			!! call Atransf(chi21,eta21,xi21,zeta21,At21)
			!! write(*,*) 'q21', chi21,eta21,xi21,zeta21
			!!:::::::::::::::::::::::::::::::::::::::::::
			!! Relative position of part. 2 relative to part. 1
			!! x21 = At1(1,1)*rx21 + At1(1,2)*ry21 ! + At1(1,3)*rz21
			!! y21 = At1(2,1)*rx21 + At1(2,2)*ry21 ! + At1(2,3)*rz21
			!! z21 = 0.d0 !! At1(3,1)*rx21 + At1(3,2)*ry21 + At1(3,3)*rz21
			!! write(*,*) 'X21 ->', x21, y21, z21
			!!:::::::::::::::::::::::::::::::::::::::::::
			!! part. 2 respect to part. 1 coordinates
			!!:::::::::::::::::::::::::::::::::::::::::::
			!! quaternions of part. 2 respect to part. 1
			call quat21(chi1,eta1,xi1,zeta1,chi2,eta2,xi2,zeta2,chi12,eta12,xi12,zeta12)
			!! Transformation matrix
			!! A12 is equal to transpose A21
			!!:::::::::::::::::::::::::::::::::::::::::::
			!! Relative position of part. 1 relative to part. 2
			x12 =-At2(1,1)*rx21 - At2(1,2)*ry21 !! - At2(1,3)*rz21
			y12 =-At2(2,1)*rx21 - At2(2,2)*ry21 !! - At2(2,3)*rz21
			z12 = 0.d0 !! -At2(3,1)*rx21 - At2(3,2)*ry21 - At2(3,3)*rz21
			!! write(*,*) 'X12 ->', x12, y12, z12
			!!:::::::::::::::::::::::::::::::::::::::::::

			!!:::::::::::::::::::::::::::::::::::::::::::
			!! Minimum distance between prticles
			!!:::::::::::::::::::::::::::::::::::::::::::
			hy21 = (abs(y21) - PRATIOy) - ( (PRATIOx*abs(At21(1,2)))**(r/(r - 1)) + &
(PRATIOy*abs(At21(2,2)))**(r/(r - 1)) + (PRATIOz*abs(At21(3,2)))**(r/(r - 1)) )**((r - 1)/r)
			hy12 = (abs(y12) - PRATIOy) - ( (PRATIOx*abs(At21(2,1)))**(r/(r - 1)) + &
(PRATIOy*abs(At21(2,2)))**(r/(r - 1)) + (PRATIOz*abs(At21(2,3)))**(r/(r - 1)) )**((r - 1)/r)

			IF( (hy21 .GT. 0.d0) .OR. (hy12 .GT. 0.d0) )THEN
				V = V + 0.d0
				GOTO 400
			ELSE
				hx21 = (abs(x21) - PRATIOx) - ( (PRATIOx*abs(At21(1,1)))**(r/(r - 1)) + &
(PRATIOy*abs(At21(2,1)))**(r/(r - 1)) + (PRATIOz*abs(At21(3,1)))**(r/(r - 1)) )**((r - 1)/r)
				hx12 = (abs(x12) - PRATIOx) - ( (PRATIOx*abs(At21(1,1)))**(r/(r - 1)) + &
(PRATIOy*abs(At21(1,2)))**(r/(r - 1)) + (PRATIOz*abs(At21(1,3)))**(r/(r - 1)) )**((r - 1)/r)

				if( (hx21 .GT. 0.d0) .OR. (hx12 .GT. 0.d0) )then
					V = V + 0.d0
					GOTO 400
				else
					!! Third approximation: particle mesh
					!!...........................................
					!! 1st mesh
					At12(1,1:3) = (/ At21(1,1), At21(2,1), At21(3,1) /)
					At12(2,1:3) = (/ At21(1,2), At21(2,2), At21(3,2) /)
					At12(3,1:3) = (/ At21(1,3), At21(2,3), At21(3,3) /)

					call SEmin1(At21,x21,y21,z21, id1,Rid1)
					call SEmin1(At12,x12,y12,z12, jd1,Rjd1)
					!!...........................................
					if( (Rid1 .LE. 1.0) .OR. (Rjd1 .LE. 1.0) ) then
						V = V + 1.0e6
						GOTO 200
					else
						!!...........................................
						!! 2nd mesh
						call SEmin2(At21,x21,y21,z21,id1, id2,Rid2)
						call SEmin2(At12,x12,y12,z12,jd1, jd2,Rjd2)
						!!...........................................
						if( (Rid2 .LE. 1.0) .OR. (Rjd2 .LE. 1.0) )then
							V = V + 1.0e6
							GOTO 200
						else
							!!...........................................
							!! 3rd mesh
							call SEmin3(At21,x21,y21,z21,id1,id2, id3,Rid3)
							call SEmin3(At12,x12,y12,z12,jd1,jd2, jd3,Rjd3)
							!!...........................................
							if( (Rid3 .LE. 1.0) .OR. (Rjd3 .LE. 1.0) )then
								V = V + 1.0e6
								GOTO 200
							else
								V = V + 0.d0
								GOTO 400
							endif
							!!...........................................
						endif
						!!...........................................
					endif
					!!...........................................
				endif
			ENDIF

		endif
		!!-------------------------------------------
		!!-------------------------------------------
		!!-------------------------------------------
		!!-------------------------------------------

		!!:::::::::::::::::::::::::::::::::::::::::::
		!! (induced) dipole-dipole interaction
		!!...........................................
400		re21 = sqrt( (x21/PRATIOx)**2.d0 + (y21/PRATIOy)**2.d0 ) !! + (z21/PRATIOz)**2.d0 )
		if( re21 .LT. rcut*2 )then
			call Uddp(rx2,ry2,rz2,At2,x21,y21,z21,At21,Epx1,Epy1,Epz1, Vdd)
			V = V + 0.42*Vdd
		endif
		!!:::::::::::::::::::::::::::::::::::::::::::

	ENDIF
	!!===========================================


100	CONTINUE
!!******************************************************************
!write(*,*) V
!write(*,*) '==========================================='

200 RETURN
end subroutine Energy
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

