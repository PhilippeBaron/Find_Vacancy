!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! MONTE CARLO SIMULATION PROGRAM IN THE CONSTANT-NVT ENSEMBLE.
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM MCNVT

COMMON / XYZ / RX, RY, RZ
COMMON / ORIENT / chi, zeta
COMMON / NUM_PAR / N
COMMON / BOXLEN / xbox, ybox, rcut
COMMON / RAD / Radx, Rady, Radz, Radm, PRATIOx, PRATIOy, PRATIOz
COMMON / IST / Istart, st
COMMON / SPH_MESH / Xe, Ye, r, Ni, idref
COMMON / ELECT / Bpw, kappa
COMMON / FIELD / lambdaV, dg, de, fcmr, fcmi, F0, Lx, Ly, Lz

!*******************************************************************
!** THIS PROGRAM TAKES A CONFIGURATION OF COLLOIDAL PARTICLES     **
!** AND PERFORMS A CONVENTIONAL NVT MC SIMULATION. THE BOX IS OF  **
!** VARYING LENGTH THERE ARE NO LOOKUP TABLES.                    **
!**                                                               **
!** PRINCIPAL VARIABLES:                                          **
!**                                                               **
!** INTEGER N                   NUMBER OF MOLECULES               **
!** INTEGER Nstep               MAXIMUM NUMBER OF CYCLES          **
!** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
!** REAL    th(N),phi(N),psi(N)	ORIENTATION                       **
!** REAL    PRATIO              PART. ASPECT RATIO                **
!** REAL    Radx, Rady, Radz    PART. RADIUS                      **
!** REAL    rcut                REDUCED CUTOFF DISTANCE           **
!** REAL    dRmax               REDUCED MAXIMUM DISPLACEMENT      **
!** REAL    V                   THE POTENTIAL ENERGY              **
!**                                                               **
!*******************************************************************

integer				N, Nv

integer, parameter :: Nmax = 1000

double precision	RX(Nmax), RY(Nmax), RZ(Nmax)

double precision	dRmax

double precision	V, Vnew, Vold, Vend, VN, delV
double precision	rxold, ryold, rzold, rxnew, rynew, &
						rznew, POTREF

double precision	AVV, ACV, ACVSQ, FLV
double precision	ACM, ACATMA, RATIO
double precision	RANF

integer				st, Nstep, Iprint, Isave, Iratio, &
					Istart

integer				dummy, I, J, A, B, C

logical				ovrlap
character			CNfile*30, CNout*60, INFout*30, &
					HISTfile*30
character			CNoutHead*30, S_id_rho*30
integer				id_rho

double precision	dum1

double precision	chi(Nmax), zeta(Nmax)
double precision	chiold, chinew
double precision	etaold, etanew
double precision	xiold, xinew
double precision	zetaold, zetanew
double precision	norm

double precision	PRATIOx, PRATIOy, PRATIOz
double precision	Radx, Rady, Radz, Radm
double precision	dQmax

integer				i1, j1, i2, j2
integer				k1, k2, k3
double precision	xbox, ybox, rcut
double precision	rho

real				TIMER, T1, T2
integer				hour, minute, second

double precision	Bpw, kappa

double precision	lambdaV, lambda, dg, de, fcmr(3)
double precision	fcmi(3), F0, Lx, Ly, Lz

double precision, parameter :: pi = 3.141592653589793

!! particle mesh
!! remember to change int1 allso in adistance.f90 
double precision	r
integer				Ni, idref
double precision	Xe(600), Ye(600)
!! Box dimensions
double precision	Xb(200), Yb(200)
integer				Nvid(200)

dummy = -6
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!	** READ INPUT DATA **
OPEN(901,file = 'run.txt')
READ(901,*)		!! Parameters for the MC-NVT code
READ(901,*)		!! Number of volumen increases (Nv)
READ(901,*) Nv
READ(901,*)		!! Number of particles (N)
READ(901,*) N
READ(901,*)		!! Number of steps (Nstep)
READ(901,*) Nstep
READ(901,*)		!! Number of steps between output lines (Iprint)
READ(901,*) Iprint
READ(901,*)		!! Number of steps between data saves (Isave)
READ(901,*) Isave
READ(901,*)		!! Interval for update of Max. Displacetmet/Rotation (Iratio)
READ(901,*) Iratio
READ(901,*)		!! Start collecting data after (Istart) steps
READ(901,*) Istart
READ(901,*)		!! Configuration input file name (CNfile)
READ(901,*) CNfile
READ(901,*)		!! Info output file name (INFout)
READ(901,*) INFout
READ(901,*)		!! Configuration output file name header (CNout)
READ(901,*) CNoutHead
READ(901,*)		!! Simulated histogram data (hist)
READ(901,*) HISTfile
READ(901,*)		!! Superellipse parameter (r = n)
READ(901,*) r
READ(901,*)		!! x-semi axis (Radx)
READ(901,*) Radx
READ(901,*)		!! y-semi axis (Rady)
READ(901,*) Rady
READ(901,*)		!! z-semi axis (Radz)
READ(901,*) Radz
READ(901,*)		!! x-box size r/a (xbox)
READ(901,*) xbox
READ(901,*)		!! y-box size r/a (ybox)
READ(901,*) ybox
READ(901,*)		!! Hard particle cut-off distance r/a (rcut)
READ(901,*) rcut
READ(901,*)		!! Number of mesh part. points (Ni)
READ(901,*) Ni
READ(901,*)		!! Mesh part. point - 1st approach (idref)
READ(901,*) idref
READ(901,*)		!! particle step size r/a (dRmax),   max=0.2
READ(901,*) dRmax
READ(901,*)		!! angular step size r/a (dQmax), max=0.1
READ(901,*) dQmax
READ(901,*)		!! Electrostatic interaction constant
READ(901,*) Bpw
READ(901,*)		!! Debye length
READ(901,*) kappa
READ(901,*)		!! Electrode gap (r/a) dg =
READ(901,*) dg
READ(901,*)		!! Electrode width (r/a) dg =
READ(901,*) de
READ(901,*)		!! lambda
READ(901,*) lambda
READ(901,*)		!! Shape parameter x-AXIS --> Lx =
READ(901,*) Lx
READ(901,*)		!! Shape parameter y-AXIS --> Ly =
READ(901,*) Ly
READ(901,*)		!! Shape parameter z-AXIS --> Lz =
READ(901,*) Lz
READ(901,*)		!! CLAUSSIUS-MOSSOTTI FACTOR at 0 Hz
READ(901,*) F0
READ(901,*)		!! REAL CLAUSSIUS-MOSSOTTI FACTOR OF x'- AXIS --> fcm1 =
READ(901,*) fcmr(1)
READ(901,*)		!! REAL CLAUSSIUS-MOSSOTTI FACTOR OF y'- AXIS --> fcm2 =
READ(901,*) fcmr(2)
READ(901,*)		!! REAL CLAUSSIUS-MOSSOTTI FACTOR OF z'- AXIS --> fcm3 =
READ(901,*) fcmr(3)
READ(901,*)		!! Imaginary CLAUSSIUS-MOSSOTTI FACTOR OF x'- AXIS --> ficm1 =
READ(901,*) fcmi(1)
READ(901,*)		!! Imaginary CLAUSSIUS-MOSSOTTI FACTOR OF y'- AXIS --> ficm2 =
READ(901,*) fcmi(2)
READ(901,*)		!! Imaginary CLAUSSIUS-MOSSOTTI FACTOR OF z'- AXIS --> ficm3 =
READ(901,*) fcmi(3)
CLOSE (901)
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Radm = MIN(Radx,Rady,Radz)
PRATIOx = Radx/Radm
PRATIOy = Rady/Radm
PRATIOz = Radz/Radm

!!write(*,*) PRATIOx, PRATIOy, PRATIOz
!!******************************************************************
!! Particle Mesh
open(4, FILE = 'mesh.txt'  )
do I = 1, 6*Ni
	read(4, *) Xe(I), Ye(I)
enddo
close(4)
!!******************************************************************
!! Box length
open(5, FILE = 'box.txt'  )
do I = 1, Nv
	read(5, *) Nvid(I), Xb(I) , Yb(I)
enddo
close(5)
!!******************************************************************

!!******************************************************************
OPEN(2, FILE = INFout  )
OPEN(7000, FILE = 'SIMINFO.TXT' )

!!** BEGIN CLOCKING THE PROGRAM **
CALL CPU_TIME(T1)

!! WRITE INPUT DATA
write(2,'('' NUMBER OF PARTICLES       '',I10   )') N
write(2,'('' NUMBER OF CYCLES          '',I10   )') Nstep
write(2,'('' OUTPUT FREQUENCY          '',I10   )') Iprint
write(2,'('' SAVE FREQUENCY            '',I10   )') Isave
write(2,'('' RATIO UPDATE FREQUENCY    '',I10   )') Iratio
write(2,'('' START COLLECTING DATA     '',I10   )') Istart
write(2,'('' CONFIGURATION FILE  NAME  '',A     )') CNfile
write(2,'('' CONFIGURATION OUT FILE    '',A     )') CNout
write(2,'('' PARTICLE X-RADIUS, NANOMETER  '',F10.4 )') Radx
write(2,'('' PARTICLE Y-RADIUS, NANOMETER  '',F10.4 )') Rady
write(2,'('' PARTICLE Z-RADIUS, NANOMETER  '',F10.4 )') Radz
write(2,'('' PARTICLE ASPECT RATIO X     '',F10.4 )') (PRATIOx)
write(2,'('' PARTICLE ASPECT RATIO Y     '',F10.4 )') (PRATIOy)
write(2,'('' PARTICLE ASPECT RATIO Z     '',F10.4 )') (PRATIOz)
write(2,'('' POTENTIAL CUTOFF, RADIUS  '',F10.4 )') rcut

!!********************************************************************
!! Initial configuration
OPEN( 3, FILE = CNfile )
DO I = 1, N
	READ(3, *)dum1, RX(I), RY(I), RZ(I), chi(I), zeta(I)
END DO
CLOSE(3)
!!********************************************************************

!! This loop is done at the end of the program
!! Is to change the box size continuously
!! Irho = part. concentration index for rho
id_rho = 1 !! Modify this to continue another simulation
DO WHILE(id_rho .LE. Nv)
	xbox = Xb(id_rho)
	ybox = Yb(id_rho)
	!!**************** density ************************
	rho = real(N)/(xbox*ybox)
	!!**********************************************************
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPEN(53, FILE = 'box_size.txt' )

if(id_rho .LT. 10) then
	write(S_id_rho, '(I1)' ) id_rho !! concentration
elseif(id_rho .GE. 10 .AND. id_rho .LT. 100) then
	write(S_id_rho, '(I2)' ) id_rho !! concentration
elseif(id_rho .GE. 100 .AND. id_rho .LT. 1000) then
	write(S_id_rho, '(I3)' ) id_rho !! concentration
elseif(id_rho .GE. 1000 .AND. id_rho .LT. 10000) then
	write(S_id_rho, '(I4)' ) id_rho !! concentration
endif

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! cut-off must be lower than the simulation box
IF( ( (xbox .LT. ybox) .AND. (rcut .GT. xbox / 2.0) ) &
	.OR. ( (ybox .LT. xbox) .AND. (rcut .GT. ybox / 2.0) ) ) &
	STOP ' CUT-OFF TOO LARGE '
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! Zero accumulatos
CNout = trim(CNoutHead)//'_'//trim(S_id_rho)//'.txt'
OPEN ( 1000, FILE = CNout )
CNout = 'mc_order_'//trim(S_id_rho)//'.txt'
OPEN (2000, file = CNout )

POTREF =  0.d0

ACATMA = 0.0
ACM    = 0.0
ACV    = 0.0
ACVSQ  = 0.0
FLV    = 0.0

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(id_rho .EQ. 1)then
	!! CALCULATE INITIAL ENERGY AND CHECK FOR OVERLAPS
	CALL SUMUP ( ovrlap, V )
	IF ( ovrlap ) THEN
		STOP ' OVERLAP IN INITIAL CONFIGURATION '
	ENDIF
endif
!!write(*, *)'INITIAL V=', V

write(2, *)'INITIAL V=', V
write(2,'(//'' START OF MARKOV CHAIN               ''//)')
write(2,'(''  NMOVE		dRmax		RATIO		V/N      ''/)')
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!!*******************************************************************
!!** LOOPS OVER ALL CYCLES AND ALL MOLECULES                       **
!!*******************************************************************
k3 = 0 !! number of samples for g(r)

open(999, FILE = 'voltage.txt')
DO 100 st = 1, Nstep
!!  write(*,*) st
    !lambdaV = (124.d0)*((3.d0/1.d0)*1.d0)**2.d0 !!
    if (st .LE. (0.6*Nstep)) then
    	lambdaV = (124.5636071)*(3.d0+(3.5/(0.6*Nstep))*st)**2.d0
    else if (st .LE. (0.8*Nstep)) then
      	lambdaV = (124.5636071)*((9.d0/(0.8*Nstep))*st)**2.d0
    else
      	lambdaV = (124.5636071)*((11.25/(Nstep))*st)**2.d0
    endif
	
    write(999,*) SQRT(lambdaV/124.5636071) 

    
	DO 99 I = 1, N
		!!========================================
		!! LOOP OVER PARTICLES
		!!========================================

		rxold = RX(I)
		ryold = RY(I)
		rzold = RZ(I)
		chiold = chi(I)
		etaold = 0.d0
		xiold = 0.d0
		zetaold = zeta(I)

		!! Calculate the Energy of I-part in the old configuration
		call Energy(rxold,ryold,rzold,chiold,etaold,xiold,zetaold,I, Vold)
		!!Vold = 0.d0

		!! Move I-part
		rxnew = rxold + ( 2.0 * RANF( dummy ) - 1.0 ) * dRmax
		rynew = ryold + ( 2.0 * RANF( dummy ) - 1.0 ) * dRmax
		rznew = rzold !! + ( 2.0 * RANF( dummy ) - 1.0 ) * dRmax
		chinew = chiold + (2.d0 * RANF(dummy) - 1.d0 ) * dQmax
		etanew = 0.d0
		xinew = 0.d0
		zetanew = zetaold + (2.d0 * RANF(dummy) - 1.d0 ) * dQmax
		!! Normalization
		norm = dsqrt(chinew**2.d0 + etanew**2.d0 + xinew**2.d0 + zetanew**2.d0)
		chinew = chinew/norm
		etanew = etanew/norm
		xinew = xinew/norm
		zetanew = zetanew/norm

		!!--------------------------------------
		!! Periodic boundary condition
		!!--------------------------------------
		rxnew = rxnew - ANINT( rxnew/xbox )*xbox
		rynew = rynew - ANINT( rynew/ybox )*ybox
		!!--------------------------------------

		!! Calculate the Energy of I-part in the new configuration
		call Energy(rxnew,rynew,rznew,chinew,etanew,xinew,zetanew,I, Vnew)

		!!--------------------------------------
		!! Check for acceptance
		!!--------------------------------------
		delV  = Vnew - Vold

		IF ( delV .LT. 75.0 ) THEN

			IF ( delV .LE. 0.d0 ) THEN

				V = V + delV

				RX(I) = rxnew
				RY(I) = rynew
				RZ(I) = rznew
				chi(I) = chinew
				zeta(I) = zetanew
				ACATMA = ACATMA + 1.0

			ELSEIF ( EXP( - delV ) .GT. RANF(dummy) ) THEN
				V = V + delV

				RX(I)  = rxnew
				RY(I)  = rynew
				RZ(I)  = rznew
				chi(I) = chinew
				zeta(I) = zetanew
				ACATMA = ACATMA + 1.0

			ENDIF

		ENDIF
		!!--------------------------------------
		!!--------------------------------------
		!write(*,*) RX(I), RZ(I), phi(I), Vnew, delV, V

		ACM = ACM + 1.0

		!! Calculate instantaneous values
		VN = ( V ) / REAL ( N )

		!! Accumulate averages
		ACV = ACV + VN
		ACVSQ = ACVSQ + VN*VN
		!!========================================
		!! ENDS LOOP OVER PARTICLES
		!!========================================
99      CONTINUE

	!!========================================
	!! PERFORM PERIODIC OPERATIONS
	!!========================================
	IF ( MOD ( st, Iratio ) .EQ. 0 ) THEN
		!!--------------------------------------
		!! Adjust maximum displacement
		!!--------------------------------------
		RATIO = ACATMA / REAL ( N * Iratio )

		if( RATIO .GT. 0.5 ) then
			dRmax  = dRmax  * 1.05
			dQmax = dQmax * 1.05
			if(dRmax > 0.05)then
				dRmax = 0.05
			endif
			if(dQmax > 0.025)then
				dQmax = 0.025
			endif
		else
			dRmax  = dRmax  * 0.95
			dQmax = dQmax * 0.95
		endif
		!!--------------------------------------
			ACATMA = 0.0
	ENDIF

	IF ( MOD ( st, Iprint ) .EQ. 0 ) THEN
		!! Write out runtime information
		write(2,'(I8,2F20.7,E20.7)') INT(ACM), dRmax, RATIO, VN
	ENDIF

	IF ( (st.GT.Istart) .AND. (MOD(st, Isave) .EQ. 0) ) THEN
		!! Write out the configuration at intervals
		!!......................................
		!! To save data (step, energy, Nematic order, Global psi6)
		write ( 2000, '(I12, F20.7)' ) st, V !, Sx ! Gpsi8

		!! To save data (step, position, orientation, local psi6)
		do I = 1, N
			write ( 1000, '(I12, I5, 5F20.7)' ) st, I, RX(I), RY(I), RZ(I), &
				chi(I), zeta(I) !!, psi8(I) !!V, Sx, Vold, Vnew, delV
		end do
	ENDIF

100	CONTINUE
!!*******************************************************************
!! ENDS THE LOOP OVER CYCLES
!!*******************************************************************
close(999)
CLOSE (1000)
CLOSE (2000)

write(2,'(//'' END OF MARKOV CHAIN          ''//)')

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!Checks the final calue of the potential energy
CALL SUMUP ( ovrlap, Vend )
!!write(*,*) 'FINAL V=', int(V/1e6)

IF ( ABS ( Vend - V ) .GT. 1.0E-03 ) THEN
	write(2,'('' PROBLEM WITH ENERGY,'')')
	write(2,'('' Vend              = '', E20.8)') Vend
	write(2,'('' V                 = '', E20.8)') V
ENDIF
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!** CALCULATE AND WRITE OUT RUNNING AVERAGES **
AVV   = ACV / ACM
ACVSQ = ( ACVSQ / ACM ) - AVV ** 2.0

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! CALCULATE FLUCTUATIONS
IF ( ACVSQ .GT. 0.0 ) FLV = SQRT ( ACVSQ )

write(2,'(/'' AVERAGES ''/ )')
write(2,'('' <V/N>   = '',F20.8)') AVV
write(2,'(/'' FLUCTUATIONS ''/)')
write(2,'('' FLUCTUATION IN <V/N> = '',F20.8)') FLV
write(2,'(/'' END OF SIMULATION '')')
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


write(53, '(I12, 3F20.7)') id_rho, xbox, ybox, rho
id_rho = id_rho + 1
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENDDO
CLOSE (53)


!** FINISH CLOCKING THE PROGRAM **
 400	CALL CPU_TIME(T2)
        TIMER  = T2 - T1
        hour   = 0
        minute = 0
        second = 0
        IF ( TIMER .GE. 3600 ) hour = INT( TIMER / 3600 )
        IF ( TIMER .GE. 60 ) minute = INT( (TIMER - 3600 * hour) / 60 )
        second = INT( ( TIMER - 3600 * hour - 60 * minute ) )
        write(*,60) hour, ':', minute, ':', second,' ELAPSED (HH:MM:SS)'
60      FORMAT(I3,A,I2.2,A,I2.2,A)

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STOP
END
