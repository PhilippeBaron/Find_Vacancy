!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Udep(rx, ry, rz, At, Vdep, Epx, Epy, Epz)

COMMON / FIELD / lambdaV, dg, de, fcmr, fcmi, F0, Lx, Ly, Lz

!!===========================================
! The main purpose of this subroutine is
! calculate the dielectriphoretic potential
! and electric field in a particle
! Outputs:
! Vdep, Epx, Epy, Epz
!!===========================================

double precision	lambdaV, dg, de, fcmr(3)
double precision	fcmi(3), F0, Lx, Ly, Lz

double precision	rx, ry, rz
double precision	At(3,3)

double precision	Epx, Epy, Epz
double precision	Vdep
double precision, parameter :: pi = 3.141592653589793

!!=====================================================
!! Electric field (lab. coordinates) - Morgan (2001)
!! J. Phys. D: Appl. Phys. 34, 1553-1561 (2001)
!! Variables are normalized by part. radius

call Epart(rx,ry,rz,At, Epx,Epy,Epz)

!write(*,*) rx, rz, Epx, Epy, Epz
!!=====================================================

!!Time averaged dipole-field interaction potential

Vdep = -1.d0*lambdaV*( &
		  fcmr(1)*Epx**2.d0  &
		+ fcmr(2)*Epy**2.d0 &
		+ fcmr(3)*Epz**2.d0 )/F0**2.d0

!write(*,*) rz, Vdep

end subroutine Udep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Epart(rx,ry,rz,At, Epx,Epy,Epz)

COMMON / FIELD / lambdaV, dg, de, fcmr, fcmi, F0, Lx, Ly, Lz

!!===========================================
! The main purpose of this subroutine is
! calculate the electric field in a particle
! Outputs:
! Epx, Epy, Epz
!!===========================================

double precision	lambdaV, de, dg, fcmr(3)
double precision	fcmi(3), F0, Lx, Ly, Lz

double precision	rx, ry, rz
double precision	At(3,3)

double precision	Ex, Ey, Ez
double precision	Epx, Epy, Epz

double precision, parameter :: pi = 3.141592653589793
!!=====================================================
!! Electric field (lab. coordinates) - Morgan (2001)
!! J. Phys. D: Appl. Phys. 34, 1553-1561 (2001)
!! Variables are normalized by part. radius


Ex = ( atan(sin((pi*(rx+dg))/(2.d0*dg) + pi/4.d0)/sinh((pi*rz)/(2.d0*dg))) &
	- atan(cos((pi*(rx+dg))/(2.d0*dg) + pi/4.d0)/sinh((pi*rz)/(2.d0*dg))) )

Ey = 0.d0

Ez = (0.50)*log( ( (cosh((pi*rz)/(2.d0*dg)) + cos((pi*(rx+dg))/(2.d0*dg) + pi/4.d0))/ &
	(cosh((pi*rz)/(2.d0*dg)) - cos((pi*(rx+dg))/(2.d0*dg) + pi/4.d0)) )* &
	( (cosh((pi*rz)/(2.d0*dg)) + sin((pi*(rx+dg))/(2.d0*dg) + pi/4.d0))/ &
	(cosh((pi*rz)/(2.d0*dg)) - sin((pi*(rx+dg))/(2.d0*dg) + pi/4.d0)) ) )

!write(*,*) rx, ry, rz, Ex, Ez
!!=====================================================

!! Electric field in part. coordinates
Epx = Ex*At(1,1) + Ey*At(1,2) + Ez*At(1,3)
Epy = Ex*At(2,1) + Ey*At(2,2) + Ez*At(2,3)
Epz = Ex*At(3,1) + Ey*At(3,2) + Ez*At(3,3)

!write(*,*) rx, ry, rz, Ex, Ez

end subroutine Epart
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Xicoord3(rx,ry,rz,x,y,z,xi)

!!===========================================
! The main purpose of this subroutine is
! calculate the position xi of the particle
! Outputs:
! xi
!!===========================================

double precision	rx, ry, rz
double precision	x, y, z

double precision	a0, a1, a2, a3

double precision	xi

double complex		jj
double complex		xic

double precision, parameter :: pi = 3.141592653589793

jj = ( 0.d0, 1.d0 )
!!=====================================================
!! Polynomial constants to calculate xi
!! xi is always related to lagest axis of ellipsoid

a0 = ry**2.d0*rz**2.d0*x**2.d0 + rx**2.d0*rz**2.d0*y**2.d0 + &
	rx**2.d0*ry**2.d0*z**2.d0 - rx**2.d0*ry**2.d0*rz**2.d0
a1 = ry**2.d0*x**2.d0 + rz**2.d0*x**2.d0 + rx**2.d0*y**2.d0 + &
	rz**2.d0*y**2.d0 + rx**2.d0*z**2.d0 + ry**2.d0*z**2.d0 - &
	rx**2.d0*ry**2.d0 - rx**2.d0*rz**2.d0 - ry**2.d0*rz**2.d0
a2 = x**2.d0 + y**2.d0 + z**2.d0 - rx**2.d0 - ry**2.d0 - rz**2.d0
a3 = - 1.d0


xic =  - a2/(3.0*a3)  + (1.0 + jj*sqrt(3.0))*(- a2**2.0 + 3.0*a1*a3)/ &
	(3.0*2.0**(2.0/3.0)*a3*(- 2.0*a2**3.0 + 9.0*a1*a2*a3 - 27.0*a0*a3**2.0 + &
	cdsqrt( (0.d0,0.d0) + 4.0*(- a2**2.0 + 3.0*a1*a3)**3.0 + &
	(- 2.0*a2**3.0 + 9.0*a1*a2*a3 - 27.0*a0*a3**2.0)**2.0) )**(1.0/3.0)) - &
	((1.0 - jj*sqrt(3.0))*(- 2.0*a2**3.0 + 9.0*a1*a2*a3 - 27.0*a0*a3**2.0 + &
	cdsqrt( (0.d0,0.d0) + 4.0*(- a2**2.0 + 3.0*a1*a3)**3.0 + &
	(- 2.0*a2**3.0 + 9.0*a1*a2*a3 - 27.0*a0*a3**2.0)**2.0) )**(1.0/3.0)) / &
	(6.0*2.0**(1.0/3.0)*a3)

xi = real(xic)

!!write(*,*) xic, xi
!!=====================================================

end subroutine Xicoord3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine XiCoord(rx,ry,x,y,xi)
!!===========================================
! The main purpose of this subroutine is
! calculate the xi-coordinate (ellipsoidal coord)
! of the point x, y
! Outputs:
! xi
!!===========================================
double precision	rx, ry
double precision	x, y

double precision	a0, a1

double precision	xi
double precision, parameter :: pi = 3.141592653589793
!!===========================================

a0 = (ry*x)**2.d0 + (rx*y)**2.d0 - (rx*ry)**2.d0
a1 = x**2.d0 + y**2.d0 - rx**2.d0 - ry**2.d0

xi = 0.50*(a1 + sqrt(4.d0*a0 + a1**2.d0) )
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end subroutine XiCoord
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Dxi(rx,ry,rz,x,y,z,xi,DxiDx,DxiDy,DxiDz)

!!===========================================
! The main purpose of this subroutine is
! calculate the derivative of xi respect
! to x, y, and z
! Outputs:
! Epx, Epy, Epz
!!===========================================

double precision	rx, ry, rz
double precision	x, y, z
double precision	xi

double precision	DxiDx, DxiDy, DxiDz

double precision, parameter :: pi = 3.141592653589793

!!=====================================================
DxiDx = ((2.0*x)/(rx**2.0 + xi))*(1.0/(x**2.0/(rx**2.0 + xi)**2.0 + &
y**2.0/(ry**2.0 + xi)**2.0 + z**2.0/(rz**2.0 + xi)**2.0))
DxiDy = ((2.0*y)/(ry**2.0 + xi))*(1.0/(x**2.0/(rx**2.0 + xi)**2.0 + &
y**2.0/(ry**2.0 + xi)**2.0 + z**2.0/(rz**2.0 + xi)**2.0))
DxiDz = ((2.0*z)/(rz**2.0 + xi))*(1.0/(x**2.0/(rx**2.0 + xi)**2.0 + &
y**2.0/(ry**2.0 + xi)**2.0 + z**2.0/(rz**2.0 + xi)**2.0))
!!=====================================================

end subroutine Dxi
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine FunXi(rx,ry,rz,x,fun1,fun2,fun3)

!!===========================================
! The main purpose of this subroutine is
! calculate the functions
! Outputs:
! fun1, fun2, fun3
!!===========================================

double precision	rx, ry, rz
double precision	x

double precision	fun1, fun2, fun3

double precision, parameter :: pi = 3.141592653589793

!!=====================================================
fun1 = (rx*ry*rz/2.0)/((x+rx**2.0)*sqrt((x+rx**2.0)*(x+ry**2.0)*(x+rz**2.0)))

fun2 = (rx*ry*rz/2.0)/((x+ry**2.0)*sqrt((x+rx**2.0)*(x+ry**2.0)*(x+rz**2.0)))

fun3 = (rx*ry*rz/2.0)/((x+rz**2.0)*sqrt((x+rx**2.0)*(x+ry**2.0)*(x+rz**2.0)))
!!=====================================================

end subroutine FunXi
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LXi(Lx,Ly,Lz,xi,Lxi_x,Lxi_y,Lxi_z)

!!===========================================
! The main purpose of this subroutine is
! calculate the functions
! Outputs:
! L(xi)_x,L(xi)_y,L(xi)_z
! instead of the integral
!!===========================================

double precision	Lx, Ly, Lz
double precision	xi

double precision	k1x, k1y, k1z
double precision	g1x, g1y, g1z
double precision	x

double precision	Lxi_x, Lxi_y, Lxi_z

double precision, parameter :: pi = 3.141592653589793

!!=====================================================
x = log(xi)

k1x = 1.02
g1x = 0.75
Lxi_x = Lx/(1.0 + exp(-k1x*x + g1x))


!BY Rachel 2019-01-28 - ALLOWS FOR 3rd dimension
k1y = 1.07
g1y = -0.3
Lxi_y = Ly/(1.0 + exp(-k1y*x + g1y))

k1z = 1.07
g1z = -0.3
Lxi_z = Lz/(1.0 + exp(-k1z*x + g1z))
!Lxi_y = Lxi_z Deleted 2019-01-28 removed by Rachel

!!=====================================================

end subroutine LXi
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Uddp(rx2,ry2,rz2,At2,x21,y21,z21,At21,Epx1,Epy1,Epz1,Vdd)

COMMON / RAD / Radx, Rady, Radz, Radm, PRATIOx, PRATIOy, PRATIOz
COMMON / FIELD / lambdaV, dg, de, fcmr, fcmi, F0, Lx, Ly, Lz

!!===========================================
! The main purpose of this subroutine is
! calculate the dipole-dipole
! interaction energy
! Outputs:
! Vdd
!!===========================================
double precision	Radx, Rady, Radz, Radm
double precision	PRATIOx, PRATIOy, PRATIOz

double precision	lambdaV, de, dg, fcmr(3)
double precision	fcmi(3), F0, Lx, Ly, Lz

double precision	rx2, ry2, rz2
double precision	At2(3,3)
double precision	x21, y21, z21
double precision	At21(3,3)
double precision	Epx1,Epy1,Epz1

double precision	Epx2,Epy2,Epz2

double precision	Xip2

double precision	DxiDx, DxiDy, DxiDz

double precision	fun1, fun2, fun3

double precision	Lxi_x, Lxi_y, Lxi_z

double precision	Exx, Eyx, Ezx
double precision	Exy, Eyy, Ezy
double precision	Exz, Eyz, Ezz

double precision	Erx, Ery, Erz
double precision	Eix, Eiy, Eiz
double precision	Erx2, Ery2, Erz2
double precision	Eix2, Eiy2, Eiz2

double precision	prx, pry, prz
double precision	pix, piy, piz

double precision	Vdd

double precision, parameter :: pi = 3.141592653589793

!!=====================================================
!! Field in the 2nd particle
call Epart(rx2,ry2,rz2,At2, Epx2,Epy2,Epz2)
!!-------------------------------------------
!! Position xi of the second particle
!! respect to the first particle
!call Xicoord3(PRATIOx,PRATIOy,PRATIOz,x21,y21,z21, Xip2)
call Xicoord(PRATIOx,PRATIOy,x21,y21, Xip2)
!!-------------------------------------------
!! Derivative of xi respect to x, y, z
call Dxi(PRATIOx,PRATIOy,PRATIOz,x21,y21,z21,Xip2, DxiDx,DxiDy,DxiDz)
!!-------------------------------------------
!! Integrant of shape parameter
call FunXi(PRATIOx,PRATIOy,PRATIOz,Xip2,fun1,fun2,fun3)
!!-------------------------------------------
!! Shape parameter as function of xi
call LXi(Lx,Ly,Lz,Xip2, Lxi_x,Lxi_y,Lxi_z)
!!-------------------------------------------
!! Dipole field in 1st particle
Exx = 3.d0*Epx1*( (Lxi_x - Lx) + x21*fun1*DxiDx )
Eyx = 3.d0*Epy1*( y21*fun2*DxiDx )
Ezx = 3.d0*Epz1*( z21*fun3*DxiDx )

Exy = 3.d0*Epx1*( x21*fun1*DxiDy )
Eyy = 3.d0*Epy1*( (Lxi_y - Ly) + y21*fun2*DxiDy )
Ezy = 3.d0*Epz1*( z21*fun3*DxiDy )

Exz = 3.d0*Epx1*( x21*fun1*DxiDz )
Eyz = 3.d0*Epy1*( y21*fun2*DxiDz )
Ezz = 3.d0*Epz1*( (Lxi_z - Lz) + z21*fun3*DxiDz )
!!-------------------------------------------
!! Dipole field (particle 1 - part.1 coord)
!! Real part
Erx = fcmr(1)*Exx + fcmr(2)*Eyx + fcmr(3)*Ezx
Ery = fcmr(1)*Exy + fcmr(2)*Eyy + fcmr(3)*Ezy
Erz = fcmr(1)*Exz + fcmr(2)*Eyz + fcmr(3)*Ezz

!! Imaginary part.
Eix = fcmi(1)*Exx + fcmi(2)*Eyx + fcmi(3)*Ezx
Eiy = fcmi(1)*Exy + fcmi(2)*Eyy + fcmi(3)*Ezy
Eiz = fcmi(1)*Exz + fcmi(2)*Eyz + fcmi(3)*Ezz
!!-------------------------------------------
!!! Dipole field (particle 1 - part.2 coord)
!!! Real part
Erx2 = At21(1,1)*Erx + At21(1,2)*Ery + At21(1,3)*Erz
Ery2 = At21(2,1)*Erx + At21(2,2)*Ery + At21(2,3)*Erz
Erz2 = At21(3,1)*Erx + At21(3,2)*Ery + At21(3,3)*Erz
!
!!! Imaginary part.
Eix2 = At21(1,1)*Eix + At21(1,2)*Eiy + At21(1,3)*Eiz
Eiy2 = At21(2,1)*Eix + At21(2,2)*Eiy + At21(2,3)*Eiz
Eiz2 = At21(3,1)*Eix + At21(3,2)*Eiy + At21(3,3)*Eiz

!! Since At12 = transpose(At21) then

!! Dipole field (particle 1 - part.2 coord)
!! Real part
!Erx2 = At21(1,1)*Erx + At21(2,1)*Ery + At21(3,1)*Erz
!Ery2 = At21(1,2)*Erx + At21(2,2)*Ery + At21(3,2)*Erz
!Erz2 = At21(1,3)*Erx + At21(2,3)*Ery + At21(3,3)*Erz

!! Imaginary part.
!Eix2 = At21(1,1)*Eix + At21(2,1)*Eiy + At21(3,1)*Eiz
!Eiy2 = At21(1,2)*Eix + At21(2,2)*Eiy + At21(3,2)*Eiz
!Eiz2 = At21(2,3)*Eix + At21(2,3)*Eiy + At21(3,3)*Eiz
!!-------------------------------------------
!! polarization part. 2 ( in part. 2 coordinates)
!! real part
prx = fcmr(1)*Epx2
pry = fcmr(2)*Epy2
prz = fcmr(3)*Epz2

!! imaginary part
pix = fcmi(1)*Epx2
piy = fcmi(2)*Epy2
piz = fcmi(3)*Epz2
!!-------------------------------------------
!! Dipole-dipole interaction energy
Vdd = -lambdaV*((prx*Erx2 + pry*Ery2 + prz*Erz2) + &
(pix*Eix2 + piy*Eiy2 + piz*Eiz2))/F0**2.0
Vdd = Vdd/2
!!=====================================================

end subroutine Uddp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



