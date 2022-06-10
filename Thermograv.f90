!--------------------------------------------------------------------------------------------
! Program package
! Thermo-gravitational effects: Titan's subsurface liquids
! Written by: Sugata Tan - 2021
! Funded by NASA SSW Grant 80NSSC19K0792
!
! Two subroutines that are equation of state (EOS) specific, to be used with EOS chosen by users, 
! are not provided here:
!	(1) DENSITY	-	to calculate the density of a phase at T, P, and composition
!	(2) FCOEF	-	to calculate the fugacity coefficient of a phase at T, P, and 
!                                composition
!
! A subroutine to calculate the inverse of an (N x N) matrix is not provided here
!	INVERSE(A,N,D) - to calculate the inverse of matrix A of dimension (N x N); the outputs
!                                are the inverse A and determinant D
!--------------------------------------------------------------------------------------------
! Defining variables universally used across the main program and its subroutines

  MODULE universal_vars

	IMPLICIT NONE
	DOUBLE PRECISION :: P0, dh, dT, grav, R, pi
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: molWt
	INTEGER :: nComp

  END MODULE universal_vars
!--------------------------------------------------------------------------------------------
! The main program

  PROGRAM Thermograv

	USE universal_vars

	PARAMETER (nC = 3)	! Titan's fluid is assumed to be a ternary mixture of
                             	! N2(1)/CH4(2)/C2H6(3)
	DOUBLE PRECISION :: MW(nC), Z0_0(nC), Z0_1(nC), Z0_2(nC), Z(nC), X(nC), Y(nC)
	DOUBLE PRECISION :: T0, q, rho
	DOUBLE PRECISION :: h, hmin, hmax, T, TL, TU, P
	INTEGER :: lit, i, n

	nComp = nC
	ALLOCATE(molWt(nComp))
	R = 83.1447d0		! Gas constant in [bar.cc/mol/K]

! Titan's fluid; the parameters of pure components depend on the equation of state in use 
! (not included here)

	DATA MW /28.018d0, 16.043d0, 30.07d0/		! The molecular weights of the components:
							! N2, CH4, C2H6
! The composition at surface in equilibrium with the atmosphere in low latitudes

	DATA Z0_0 /0.069676714d0, 0.367302904d0, 0.563020382d0/	

! The composition at surface in equilibrium with the atmosphere in high latitudes

	DATA Z0_1 /0.203892325d0, 0.700433111d0, 0.095674565d0/	

! Titan
	
	grav = 1.352d0		! [m/s2]	gravity
	P0 = 1.467d0		! [bar]		pressure at the surface
	molWt = MW		! molecular weight of the components

! Cases (the value in the last order is the active one for calculation, swap for the desired case)

	lit = 0			! subsurface liquid in water-ice crust in low latitudes
	lit = 1			! subsurface liquid in water-ice crust in high latitudes
	lit = 2			! subsurface liquid in methane-clathrate crust in low latitudes
	lit = 3			! subsurface liquid in methane-clathrate crust in high latitudes

! Layering

	hmin = 0.d0			! from the surface
	dh = -1.d0			! [meter]	depth increment set to be 1 meter

	SELECT CASE (lit)
			
		CASE (0)		! subsurface liquid in water-ice crust in low latitudes

			T0 = 94.d0	! [K]		temperature of surface in low latitudes
			Z = Z0_0
			hmax = -3.5d4	! [meter]	maximum depth of 35 km
			q = 13.8d-3	! [W/m2]	heat flux

		CASE (1)		! subsurface liquid in water-ice crust in high latitudes

			T0 = 90.d0	! [K]		temperature of surface in low latitudes
			Z = Z0_1
			hmax = -3.5d4	! [meter]	maximum depth of 35 km 
			q = 13.8d-3	! [W/m2]	heat flux

		CASE (2)		! subsurface liquid in CH4-clathrate crust in low latitudes
			
			T0 = 94.d0	! [K]		temperature of surface in low latitudes
			Z = Z0_0
			hmax = -6.d3	! [meter]	maximum depth of 6 km 
			q = 5.425d-3	! [W/m2]	heat flux

		CASE (3)		! subsurface liquid in CH4-clathrate in high latitudes
			
			T0 = 90.d0	! [K]		temperature of surface in high latitudes
			Z = Z0_1
			hmax = -6.d3	! [meter]	maximum depth of 6 km 
			q = 5.425d-3	! [W/m2]	heat flux

	END SELECT

	OPEN(2,file='results.txt')	! Opening a text file to record the results

! Calculating property profiles of subsurface liquids from hmin to hmax
	X = Z; Y = Z			! Initial guess for the composition at the layer boundaries 

	DO h = hmin, hmax, dh

		ht = h+dh/2.d0		! the middle point of a layer

		! Temperature profile obtained from Eq (A1)

		IF (lit > 1) THEN	! methane-clathrate crust
			T = T0-q*ht/0.5d0
			TL = T0-q*h/0.5d0
			TU = T0-q*(h+dh)/0.5d0
			n = 100
		ELSE			! water-ice crust
			T = (T0**0.0248d0-0.0248d0*ht*q/10.d0**2.7154d0)**(1.d0/0.0248d0)
			TL = (T0**0.0248d0-0.0248d0*h*q/10.d0**2.7154d0)**(1.d0/0.0248d0)
			TU = (T0**0.0248d0-0.0248d0*(h+dh)*q/10.d0**2.7154d0)**(1.d0/0.0248d0)
			n = 1000
		END IF

		! Calculate the pressure, density, and composition profiles

		dT = TU-TL	! Temperature increment between lower and upper boundaries of layer

		! EOS-specific subroutine, using EOS chosen by user, not provided here
		CALL DENSITY(1,P0,T,Z,rho)
		! Initial value to start the calculation of pressure gradient	
		P = P0-rho*sum(z*MW)*grav*dh*1.d-2	

		! Calculate the pressure [bar] and composition Y at the depth of h
		CALL P_COMP(T,P,Z,X,Y)
		! Calculate the density in [mol/cc] at the depth of h	 
		CALL DENSITY(1,P,T,Y,rho)	
		rho = rho*sum(Y*MW)*1.d3	! Convert to the unit of [kg/m3]

		! At this point, X = Z is the composition of the upper boundary and Y is of the 
		! lower boundary of the layer

		! Print out the results and store in the output file every n meters of depth 

		i = INT(h); i = ABS(i)
		IF(MOD(i,n) == (n+INT(dh))) THEN
  			 WRITE(2,"(X,F9.2,2(2X,F10.5),3X,3(F15.9,3X),F9.4)") DABS(h+dh), T, P, & 
					&Y(1), Y(2), Y(3), rho
			 PRINT*, DABS(h+dh)	! For display on screen
		END IF

		! Updating for the next iteration
		Z = Y
		Y = X
		P0 = P

	END DO
	STOP

END PROGRAM Thermograv

!--------------------------------------------------------------------------------------------
! Subroutine to calculate pressure and composition due to thermo-gravitational effects using the 
! algorithm for bubble points at T and Z (see the main paper)

SUBROUTINE P_COMP(T,P,Z,X,Y)

	USE universal_vars

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: T, Z(nComp)
	DOUBLE PRECISION, INTENT(INOUT) :: P, X(nComp),Y(nComp)
	DOUBLE PRECISION :: YY(nComp), K(nComp), KK(nComp), dK_dP(nComp), dK_dY(nComp,nComp)
	DOUBLE PRECISION :: diff, det, delta(nComp+1), func(nComp+1), jac(nComp+1,nComp+1)
	INTEGER :: i, j
! Newton-Raphson method to solve the phase-equilibrium equations

	diff = 1.d-4; det=1.d0
	X = Z
	DO WHILE (det.gt.1.0D-10)

! Calculate K factors   
 
		CALL KFACT(P,T,X,Y,K)
                                                     
! Calculate derivative of K factor with respect to composition (mol fraction)

		YY = Y
		DO i = 1,nComp
			DO j = 1,nComp
				YY(j)=Y(j)*(1.d0+diff)
				CALL KFACT(P,T,X,YY,KK)
				YY=Y
				dK_dY(i,j) = (KK(i)-K(i))/(diff*Y(j))
			END DO
		END DO

! Calculate derivative of K factor with respect to P or T

		CALL KFACT(P*(1.d0+diff),T,X,Y,KK)
		dK_dP = (KK-K)/(diff*P)

! The zero functions
 
        	DO i = 1,nComp
	       	 	func(i) = -Y(i)+K(i)*Z(i)	! energy balance (equifugacity)
		END DO
		func(nComp+1) = -1.d0+SUM(Y)		! consistency equation

! Build the Jacobian for bubble-point calculations
 
        	DO i = 1,nComp
	      		DO j=1,nComp
				jac(i,j) = Z(i)*dK_dY(i,j)
		  	END DO
		  	jac(i,i) = jac(i,i)-1.0D0
		  	jac(i,nComp+1)=Z(i)*dK_dP(i)
		  	jac(nComp+1,i)= 1.d0
		END DO
		jac(nComp+1,nComp+1) = 0.0D0

! Using inverse of Jacobian to solve the linear equations
 
		CALL INVERSE(jac,nComp+1,det)		! The INVERSE subroutine is not included

		delta = 0.d0
		DO i = 1,nComp+1
			DO j = 1,nComp+1
				delta(i) = delta(i)+jac(i,j)*func(j)
			END DO
		END DO

		det = 0.d0
		DO i = 1,nComp+1
			det = det + delta(i)**2
		END DO
		det = DSQRT(det/(nComp+1))

! Updating the variables

		P = P-delta(nComp+1)
		Y = Y-delta(1:nComp)

	END DO

END SUBROUTINE P_COMP
!--------------------------------------------------------------------------------------------
! Subroutine to calculate the K-factor at P, T, X, Y
! Originally, it is for vapor-liquid equilibria, but now modified for thermo-gravity purposes	
	
SUBROUTINE KFACT(P,T,X,Y,K)

	USE universal_vars

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: P, T, X(nComp), Y(nComp)
  	DOUBLE PRECISION, INTENT(OUT) :: K(nComp)
	DOUBLE PRECISION :: ln_fugacityCoef_v(nComp)
	DOUBLE PRECISION :: ln_fugacityCoef_l(nComp), Qnet(nComp)
	INTEGER :: phase

	phase = 1	! For common vapor-liquid equilibria, change to vapor phase (phase = 0)

! EOS-specific subroutine to calculate phase fugacity coefficient (not provided here)

	CALL FCOEF(phase,P,T,Y,ln_fugacityCoef_v)	
	CALL FCOEF(1,P0,T,X,ln_fugacityCoef_l)                    

! Calculate K-factor     
                             
	K = DEXP(ln_fugacityCoef_l-ln_fugacityCoef_v)

! Inclusion of thermo-gravitational effects
		
	CALL Qdiff(1,P,T,Y,Qnet)	! Thermal diffusion in subsurface liquid at P, T, Y
	K = K*dexp(-molWt*1.d-2*grav*(dh)/R/T-Qnet/T*dT)*P0/P	

END SUBROUTINE KFACT

!------------------------------------------------------------------------------
! Subroutine to calculate the net heat diffusion in a phase at P, T, X
! Based on Firoozabadi et al. (2000)

SUBROUTINE Qdiff(phase,P,T,X,Qnet)

	USE universal_vars
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::P, T, X(nComp)
	DOUBLE PRECISION, INTENT(OUT) :: Qnet(nComp)
	INTEGER, INTENT(IN) :: phase
	DOUBLE PRECISION :: V(nComp), H(nComp), rho, rho1, rho2, ln_fugacityCoef(nComp)
	DOUBLE PRECISION :: d_lnPhi_dP(nComp), d_lnPhi_dT(nComp), U(nComp)

! Calculate the partial molar volume

	CALL FCOEF(phase,P*1.0001d0,T,X,ln_fugacityCoef)
	CALL FCOEF(phase,P*0.9999d0,T,X,d_lnPhi_dP)
	d_lnPhi_dP = (ln_fugacityCoef-d_lnPhi_dP)/2.d-4/P
	V = (d_lnPhi_dP + 1.d0/P)*R*T	
	CALL DENSITY(phase,P,T,X,rho)

! Calculate the partial molar enthalpy

	CALL FCOEF(phase,P,T*1.0001d0,X,ln_fugacityCoef)
	CALL FCOEF(phase,P,T*0.9999d0,X,d_lnPhi_dT)
	d_lnPhi_dT = (ln_fugacityCoef-d_lnPhi_dT)/2.d-4/T
	H = -T*d_lnPhi_dT

! Calculate the partial internal energy

	U = H - P*d_lnPhi_dP

! Calculate the output

	Qnet = (-U+V*sum(X*U)/sum(X*V))/4.d0

END SUBROUTINE Qdiff
