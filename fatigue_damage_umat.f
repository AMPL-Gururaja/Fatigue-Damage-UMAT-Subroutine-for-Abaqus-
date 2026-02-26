!!!!!!!!!!! UMAT suboutine for mechanical analysis !!!!!!!!!
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

C      INCLUDE 'ABA_PARAM.INC'
      parameter (nprecd=2)
      implicit real*8(a-h,o-z)

C
      CHARACTER*80 CMNAME
C      
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

      DIMENSION STRAN_TOTAL(6),Y_THERMO(4),DAMAGE_RATE(4),DAMAGE(4),DAMAGE_OLD(4),STRESS_OLD(6)
      real E11_0, E22_0, G12_0, G23_0, NU12, NU23, GAMMA11, M11,LAMBDA11,DELTA11,TEMP_FACTOR,TEMP_DEPENDENCY
      real GAMMA22, M22,LAMBDA22,DELTA22,GAMMA12, M12,LAMBDA12,DELTA12,HYSTERESIS_INCREMENT,STRESS_AVG
      real GAMMA23, M23,LAMBDA23,DELTA23,FREQUENCY,NU13,G13_0, TERM11, TERM22
      real*8 E, nu, rho ,cp, DELTA_DENOM	  	  
	  DOUBLE PRECISION Q_friction,w_diss_visco, delta_T,w_diss_damage,w_diss_total,delta_T_total
      logical is_complete_cycle
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, HALF=0.5D0, FOUR=4.D0, NINE=9.D0, TOLER=1.D-12)
      REAL*8 DDSDDE_INV(6,6), STRAN_NEW(6), prony(3), D_prony(3), strain_ve_elem_old(3,6), strain_ve_elem(3,6)
      REAL*8 STRAN_OLD(6), strain_ve(6), strain_el(6), strain_ve_old(6), dstrain_ve_elem(3,6)
      INTEGER INFO_INV, I, J, idx
      REAL*8 DSTRAN_ACTUAL(6), sigma_ve_elem_old(3,6), sigma_ve_elem_new(3,6)
	  
      integer CYCLE_JUMP, N_STEPS, i_step, i_comp
      real*8 dN, D_old(4), D_new(4), D_temp(4), max_dD, dD
      real*8 MODULUS, NU21, NU31, NU32
      real*8 N_CYCLE, N_CYCLE_OLD, N_CYCLE_JUMP, N_CYCLE_STEP
      parameter (CYCLE_JUMP = 90)
	  
c *** get Properties, and Damage model parameters (matrix)
      E11_0 = PROPS(1)
      E22_0 = PROPS(2)
      G12_0 = PROPS(3)
      G23_0 = PROPS(4)
      NU12  = PROPS(5)
      NU23  = PROPS(6)
      rho   = PROPS(7)
      cp    = PROPS(8)
      prony(1) = PROPS(9)
      prony(2) = PROPS(10)
      prony(3) = PROPS(11)
      D_prony(1) = PROPS(12)
      D_prony(2) = PROPS(13)
      D_prony(3) = PROPS(14)
      GAMMA11 = PROPS(15)
      M11 = PROPS(16)
      LAMBDA11 = PROPS(17)
      DELTA11= PROPS(18)
      GAMMA22 = PROPS(19)
      M22 = PROPS(20)
      LAMBDA22 = PROPS(21)
      DELTA22 = PROPS(22)
      GAMMA12 = PROPS(23)
      M12 = PROPS(24)
      LAMBDA12 = PROPS(25)
      DELTA12 = PROPS(26)
      GAMMA23 = PROPS(27)
      M23 = PROPS(28)
      LAMBDA23 = PROPS(29)
      DELTA23 = PROPS(30)
      FREQUENCY = PROPS(31)		
	  
c ***    For transversely isotropic material: ! NU13 = NU12, G13 = G12
      NU13 = NU12
      G13_0 = G12_0
	  	  
c ***    Strain calculation (Step 1)
      STRAN_TOTAL(1) = STRAN(1) + DSTRAN(1)
      STRAN_TOTAL(2) = STRAN(2) + DSTRAN(2)
      STRAN_TOTAL(3) = STRAN(3) + DSTRAN(3)
      STRAN_TOTAL(4) = STRAN(4) + DSTRAN(4)
      STRAN_TOTAL(5) = STRAN(5) + DSTRAN(5)
      STRAN_TOTAL(6) = STRAN(6) + DSTRAN(6)

      DO I = 1, 6
         strain_ve(I) = STATEV(I)
         strain_el(I) = STRAN_TOTAL(I)-strain_ve(I)
		 strain_ve_old(I) = strain_ve(I)
      END DO

      w_diss_visco = 0.0D0
	  
c ***    Initialize stiffness matrix
      DO I = 1, 3
          DO J = 1, 6
              idx = 14 + (I-1)*6 + J
			  strain_ve_elem_old(I,J) = STATEV(idx)
          END DO
      END DO

c ***    Initialize stiffness matrix
      DO I = 1,3
        DO J = 1,6
          ! Assign correct modulus for each component
          SELECT CASE (J)
            CASE (1)
              MODULUS = E11_0
            CASE (2)
              MODULUS = E22_0
            CASE (3)
              MODULUS = E22_0
            CASE (4)
              MODULUS = G12_0
            CASE (5)
              MODULUS = G13_0
            CASE (6)
              MODULUS = G23_0
          END SELECT
          strain_ve_elem(I,J) = strain_ve_elem_old(I,J) + (DTIME*prony(I))*(strain_el(J) - strain_ve_elem_old(I,J))
          dstrain_ve_elem(I,J) = strain_ve_elem(I,J) - strain_ve_elem_old(I,J)
          sigma_ve_elem_old(I,J) = MODULUS * D_prony(I) * strain_ve_elem_old(I,J)
          sigma_ve_elem_new(I,J) = MODULUS * D_prony(I) * strain_ve_elem(I,J)
          w_diss_visco = w_diss_visco + 0.5*(sigma_ve_elem_old(I,J)+sigma_ve_elem_new(I,J))*dstrain_ve_elem(I,J)
          idx = 14 + (I-1)*6 + J
          STATEV(idx) = strain_ve_elem(I,J)
        END DO
      END DO

      DO I = 1, 6
          strain_ve(I) = 0.0	  
          DO J = 1, 3
              strain_ve(I) = strain_ve(I) + strain_ve_elem(J,I)
          END DO
		  STATEV(I) = strain_ve(I)
		  strain_el(I) = STRAN_TOTAL(I) - strain_ve(I)
		  STATEV(6+I) = strain_el(I)
      END DO

      STATEV(13) = STATEV(13) + w_diss_visco
	  delta_T = w_diss_visco / (rho * cp)*1000000
	  STATEV(14) = STATEV(14) + delta_T

      DO I = 1, NTENS
         DDSDDT(I) = 0.0D0
         DRPLDE(I) = 0.0D0
      END DO
      DRPLDT = 0.0D0
	  
c ***    Initialize state variable for Damage
      DAMAGE_OLD(1) = STATEV(33)
      DAMAGE_OLD(2) = STATEV(34)
      DAMAGE_OLD(3) = STATEV(35)
      DAMAGE_OLD(4) = STATEV(36)

	  
c ***    Calculate damaged elastic moduli
      E11 = E11_0* (ONE - DAMAGE_OLD(1))
      E22 = E22_0* (ONE - DAMAGE_OLD(2))
      E33 = E22  ! Transversely isotropic: E33 = E22
      G12 = G12_0* (ONE - DAMAGE_OLD(3))
      G13 = G12  ! Transversely isotropic: G13 = G12
      G23 = G23_0 * (ONE - DAMAGE_OLD(4))

c ***    Ensure non-negative moduli
      IF (E11 .LE. TOLER) E11 = E11_0 * TOLER
      IF (E22 .LE. TOLER) E22 = E22_0 * TOLER
      IF (G12 .LE. TOLER) G12 = G12_0 * TOLER
      IF (G23 .LE. TOLER) G23 = G23_0 * TOLER
      E33 = E22
      G13 = G12
	  
c ***    Calculate NU21 and NU31 from symmetry conditions
      NU21 = NU12 * E22 / E11
      NU31 = NU13 * E33 / E11
      NU32 = NU23

c ***    Build damaged stiffness matrix for transversely isotropic material
      DELTA_DENOM = ONE - NU12*NU21 - NU23*NU32 - NU13*NU31- TWO*NU21*NU32*NU13
	  
c ***    Initialize stiffness matrix
      DO I = 1, NTENS
          DO J = 1, NTENS
              DDSDDE(I,J) = ZERO
          END DO
      END DO
	  
!     CALCULATE THE TANGENT STIFFNESS (Step 2)
	  DDSDDE(1,1) = E11*(ONE - NU23*NU32) / DELTA_DENOM
	  DDSDDE(2,2) = E22*(ONE - NU13*NU31) / DELTA_DENOM
	  DDSDDE(3,3) = E33*(ONE - NU12*NU21) / DELTA_DENOM
	  DDSDDE(1,2) = E11*(NU21 + NU31*NU23) / DELTA_DENOM
	  DDSDDE(1,3) = E11*(NU31 + NU21*NU32) / DELTA_DENOM
	  DDSDDE(2,1) = DDSDDE(1,2)
	  DDSDDE(2,3) = E22*(NU32 + NU12*NU31) / DELTA_DENOM
	  DDSDDE(3,1) = DDSDDE(1,3)
	  DDSDDE(3,2) = DDSDDE(2,3)
	  DDSDDE(4,4) = G12
	  DDSDDE(5,5) = G13
	  DDSDDE(6,6) = G23
	  
c ***   UPDATE THE STRESS (Step 3)
       DO K1=1,NTENS
         STRESS(K1)=0.0D0
         DO K2=1,NTENS
            STRESS(K1)=STRESS(K1)+DDSDDE(K2,K1)*strain_el(K2)
        ENDDO
      ENDDO	  

c ***Calculate thermodynamic dual variables (Y_ij) based on strain energy derivatives (Step 5)
      
c ***Y11 calculation
      TERM11 = STRAN_TOTAL(1) + NU21 * STRAN_TOTAL(2) + NU31 * STRAN_TOTAL(3)
      IF (TERM11 .GT. ZERO) THEN
          Y_THERMO(1) = HALF * E11 * STRAN_TOTAL(1) * TERM11 / DELTA_DENOM
      ELSE
          Y_THERMO(1) = ZERO
      END IF
	  
c ***Y22 calculation  
      TERM22 = STRAN_TOTAL(2) + NU12 * STRAN_TOTAL(1) + NU32 * STRAN_TOTAL(3)
      IF (TERM22 .GT. ZERO) THEN
          Y_THERMO(2) = HALF * E22 * STRAN_TOTAL(2) * TERM22 / DELTA_DENOM
      ELSE
          Y_THERMO(2) = ZERO
      END IF
	  
c ***Shear components     
	  Y_THERMO(3) = HALF * G12 * STRAN_TOTAL(4)**2  ! γ12
	  Y_THERMO(4) = HALF * G23 * STRAN_TOTAL(6)**2  ! γ23

c ***Calculate number of cycles
      N_CYCLE_STEP = DTIME*FREQUENCY
	  STATEV(38)=STATEV(37)+N_CYCLE_STEP

c *** Calculate damage (Step 6)
c ***d11 evolution
      IF (Y_THERMO(1) .GT. TOLER .AND. M11 .GT. TOLER) THEN
          TERM1 = (GAMMA11 * DAMAGE_OLD(1)**M11)
          TERM2 = LAMBDA11 * Y_THERMO(1) * EXP(-DELTA11 * STATEV(6))
          DAMAGE_RATE(1) = (TERM1 + TERM2)
      ELSE
          DAMAGE_RATE(1) = ZERO
      END IF

c ***d22 evolution
      IF (Y_THERMO(2) .GT. TOLER .AND. M22 .GT. TOLER) THEN
          TERM1 = (GAMMA22 * DAMAGE_OLD(2)**M22)
          TERM2 = LAMBDA22 * Y_THERMO(2) * EXP(-DELTA22 * STATEV(6))
          DAMAGE_RATE(2) = (TERM1 + TERM2)
      ELSE
          DAMAGE_RATE(2) = ZERO
      END IF

c ***d12 evolution
      IF (Y_THERMO(3) .GT. TOLER .AND. M12 .GT. TOLER) THEN
          TERM1 = (GAMMA12 * DAMAGE_OLD(3)**M12)
          TERM2 = LAMBDA12 * Y_THERMO(3) * EXP(-DELTA12 * STATEV(6))
          DAMAGE_RATE(3) = (TERM1 + TERM2) 
      ELSE
          DAMAGE_RATE(3) = ZERO
      END IF

c ***d23 evolution
      IF (Y_THERMO(4) .GT. TOLER .AND. M23 .GT. TOLER) THEN
          TERM1 = (GAMMA23 * DAMAGE_OLD(4)**M23)
          TERM2 = LAMBDA23 * Y_THERMO(4) * EXP(-DELTA23 * STATEV(6))
          DAMAGE_RATE(4) = (TERM1 + TERM2)
      ELSE
          DAMAGE_RATE(4) = ZERO
      END IF

c ***Update damage variables using explicit integration (Step 7)
      STATEV(33)= DAMAGE_OLD(1) + DAMAGE_RATE(1) * (STATEV(38)-STATEV(37))
      STATEV(34)= DAMAGE_OLD(2) + DAMAGE_RATE(2) * (STATEV(38)-STATEV(37))	  
      STATEV(35)= DAMAGE_OLD(3) + DAMAGE_RATE(3) * (STATEV(38)-STATEV(37))
      STATEV(36)= DAMAGE_OLD(4) + DAMAGE_RATE(4) * (STATEV(38)-STATEV(37))  

!     Ensure Damage is Non-Decreasing
      IF (STATEV(33) .LT. DAMAGE_OLD(1)) then
         STATEV(33) = DAMAGE_OLD(1)
      END IF

      IF (STATEV(34) .LT. DAMAGE_OLD(2)) then
         STATEV(34) = DAMAGE_OLD(2)
      END IF

      IF (STATEV(35) .LT. DAMAGE_OLD(3)) then
         STATEV(35) = DAMAGE_OLD(3)
      END IF

      IF (STATEV(36) .LT. DAMAGE_OLD(4)) then
         STATEV(36) = DAMAGE_OLD(4)
      END IF
	  
	  IF (STATEV(33) .GT. 0.99D0) STATEV(33) = 0.99D0
	  IF (STATEV(34) .GT. 0.99D0) STATEV(34) = 0.99D0
	  IF (STATEV(35) .GT. 0.99D0) STATEV(35) = 0.99D0
	  IF (STATEV(36) .GT. 0.99D0) STATEV(36) = 0.99D0


      ! Check if we're currently in a cooldown period
      IF (STATEV(42) .EQ. 1.0D0) THEN
          ! In cooldown - check if cooldown period has ended
          IF (STATEV(38) .GE. STATEV(41)) THEN
              STATEV(42) = 0.0D0  ! End cooldown
          END IF
      END IF

      ! Only perform cycle jump check if:
      ! 1. We have enough cycles (>= 150)
      ! 2. We're not in a cooldown period
      ! 3. We haven't already jumped in this increment
      IF (STATEV(38) .GE. 150.0D0 .AND. STATEV(42) .EQ. 0.0D0) THEN
          
          ! Prepare for cycle jump check
          do i_comp = 1, 4
            D_old(i_comp) = DAMAGE_OLD(i_comp)
          end do
          N_CYCLE_OLD = STATEV(37)
          N_CYCLE = STATEV(38)

          ! Compute damage rates for cycle jump evaluation
          if (Y_THERMO(1) .GT. TOLER .AND. M11 .GT. TOLER) then
            TERM1 = (GAMMA11 * D_old(1)**M11)
            TERM2 = LAMBDA11 * Y_THERMO(1) * EXP(-DELTA11 * N_CYCLE)
            DAMAGE_RATE(1) = (TERM1 + TERM2)
          else
            DAMAGE_RATE(1) = ZERO
          end if
          if (Y_THERMO(2) .GT. TOLER .AND. M22 .GT. TOLER) then
            TERM1 = (GAMMA22 * D_old(2)**M22)
            TERM2 = LAMBDA22 * Y_THERMO(2) * EXP(-DELTA22 * N_CYCLE)
            DAMAGE_RATE(2) = (TERM1 + TERM2)
          else
            DAMAGE_RATE(2) = ZERO
          end if
          if (Y_THERMO(3) .GT. TOLER .AND. M12 .GT. TOLER) then
            TERM1 = (GAMMA12 * D_old(3)**M12)
            TERM2 = LAMBDA12 * Y_THERMO(3) * EXP(-DELTA12 * N_CYCLE)
            DAMAGE_RATE(3) = (TERM1 + TERM2)
          else
            DAMAGE_RATE(3) = ZERO
          end if
          if (Y_THERMO(4) .GT. TOLER .AND. M23 .GT. TOLER) then
            TERM1 = (GAMMA23 * D_old(4)**M23)
            TERM2 = LAMBDA23 * Y_THERMO(4) * EXP(-DELTA23 * N_CYCLE)
            DAMAGE_RATE(4) = (TERM1 + TERM2)
          else
            DAMAGE_RATE(4) = ZERO
          end if

          ! Calculate predicted damage after one cycle
          do i_comp = 1, 4
            D_new(i_comp) = D_old(i_comp) + DAMAGE_RATE(i_comp) * (N_CYCLE - N_CYCLE_OLD)
            if (D_new(i_comp) .GT. 0.99d0) D_new(i_comp) = 0.99d0
            if (D_new(i_comp) .LT. 0.0d0) D_new(i_comp) = 0.0d0
          end do

          ! Find maximum damage change
          max_dD = 0.0d0
          do i_comp = 1, 4
            dD = abs(D_new(i_comp) - D_old(i_comp))
            if (dD .GT. max_dD) max_dD = dD
          end do

          ! Check if cycle jump criteria is met
          if (max_dD .LT. 0.01d0) then
              ! Perform cycle jump
              N_CYCLE_JUMP = STATEV(37) + CYCLE_JUMP
              
              ! Update damage using the jumped cycle number
              do i_comp = 1, 4
                D_new(i_comp) = D_old(i_comp)
              end do
              
              ! Recalculate damage rates with jumped cycle number
              TERM1 = (GAMMA11 * D_new(1)**M11)
              TERM2 = LAMBDA11 * Y_THERMO(1) * EXP(-DELTA11 * N_CYCLE_JUMP)
              DAMAGE_RATE(1) = (TERM1 + TERM2)
              TERM1 = (GAMMA22 * D_new(2)**M22)
              TERM2 = LAMBDA22 * Y_THERMO(2) * EXP(-DELTA22 * N_CYCLE_JUMP)
              DAMAGE_RATE(2) = (TERM1 + TERM2)
              TERM1 = (GAMMA12 * D_new(3)**M12)
              TERM2 = LAMBDA12 * Y_THERMO(3) * EXP(-DELTA12 * N_CYCLE_JUMP)
              DAMAGE_RATE(3) = (TERM1 + TERM2)
              TERM1 = (GAMMA23 * D_new(4)**M23)
              TERM2 = LAMBDA23 * Y_THERMO(4) * EXP(-DELTA23 * N_CYCLE_JUMP)
              DAMAGE_RATE(4) = (TERM1 + TERM2)
              
              ! Update damage with jumped cycle
              do i_comp = 1, 4
                D_new(i_comp) = D_new(i_comp) + DAMAGE_RATE(i_comp)
                if (D_new(i_comp) .GT. 0.99d0) D_new(i_comp) = 0.99d0
              end do
              
              ! Update state variables
              do i_comp = 1, 4
                STATEV(i_comp+32) = D_new(i_comp)
                if (STATEV(i_comp+32) .LT. DAMAGE_OLD(i_comp)) STATEV(i_comp+32) = DAMAGE_OLD(i_comp)
              end do
              
              ! Set cycle number to jumped value
              STATEV(37) = N_CYCLE_JUMP
              STATEV(38) = N_CYCLE_JUMP
              
              ! Start cooldown period (600 cycles from jumped cycle)
              STATEV(41) = N_CYCLE_JUMP + CYCLE_JUMP  ! Cooldown ends after 600 cycles
              STATEV(42) = 1.0D0                  ! Start cooldown
              STATEV(43) = N_CYCLE_JUMP           ! Record the jump cycle
              
          else
              ! No jump - normal update
              do i_comp = 1, 4
                STATEV(i_comp+32) = D_new(i_comp)
                if (STATEV(i_comp+32) .LT. DAMAGE_OLD(i_comp)) STATEV(i_comp+32) = DAMAGE_OLD(i_comp)
              end do
              STATEV(37) = N_CYCLE
          end if
          
      else
          ! Normal damage update (no cycle jump check)
          do i_comp = 1, 4
            D_old(i_comp) = DAMAGE_OLD(i_comp)
          end do
          N_CYCLE_OLD = STATEV(37)
          N_CYCLE = STATEV(38)
          
          if (Y_THERMO(1) .GT. TOLER .AND. M11 .GT. TOLER) then
            TERM1 = (GAMMA11 * D_old(1)**M11)
            TERM2 = LAMBDA11 * Y_THERMO(1) * EXP(-DELTA11 * N_CYCLE)
            DAMAGE_RATE(1) = (TERM1 + TERM2)
          else
            DAMAGE_RATE(1) = ZERO
          end if
          if (Y_THERMO(2) .GT. TOLER .AND. M22 .GT. TOLER) then
            TERM1 = (GAMMA22 * D_old(2)**M22)
            TERM2 = LAMBDA22 * Y_THERMO(2) * EXP(-DELTA22 * N_CYCLE)
            DAMAGE_RATE(2) = (TERM1 + TERM2)
          else
            DAMAGE_RATE(2) = ZERO
          end if
          if (Y_THERMO(3) .GT. TOLER .AND. M12 .GT. TOLER) then
            TERM1 = (GAMMA12 * D_old(3)**M12)
            TERM2 = LAMBDA12 * Y_THERMO(3) * EXP(-DELTA12 * N_CYCLE)
            DAMAGE_RATE(3) = (TERM1 + TERM2)
          else
            DAMAGE_RATE(3) = ZERO
          end if
          if (Y_THERMO(4) .GT. TOLER .AND. M23 .GT. TOLER) then
            TERM1 = (GAMMA23 * D_old(4)**M23)
            TERM2 = LAMBDA23 * Y_THERMO(4) * EXP(-DELTA23 * N_CYCLE)
            DAMAGE_RATE(4) = (TERM1 + TERM2)
          else
            DAMAGE_RATE(4) = ZERO
          end if
          
          do i_comp = 1, 4
            STATEV(i_comp+32) = D_old(i_comp) + DAMAGE_RATE(i_comp) * (N_CYCLE - N_CYCLE_OLD)
            if (STATEV(i_comp+32) .GT. 0.99d0) STATEV(i_comp+32) = 0.99d0
            if (STATEV(i_comp+32) .LT. 0.0d0) STATEV(i_comp+32) = 0.0d0
            if (STATEV(i_comp+32) .LT. DAMAGE_OLD(i_comp)) STATEV(i_comp+32) = DAMAGE_OLD(i_comp)
          end do
          STATEV(37) = N_CYCLE		
      end if

      ! Final check to ensure damage is non-decreasing
      do i_comp = 1, 4
        if (STATEV(i_comp+32) .LT. DAMAGE_OLD(i_comp)) STATEV(i_comp+32) = DAMAGE_OLD(i_comp)
      end do

      ! Detect new step and reset cycle-related STATEV
      IF (STATEV(45) .NE. KSTEP) THEN
         STATEV(38) = 0.0D0
         STATEV(39) = 0.0D0
         STATEV(40) = 0.0D0
         STATEV(44) = 0.0D0
		 STATEV(13) = 0.0D0
		 STATEV(14) = 0.0D0
         STATEV(45) = KSTEP
      END IF
	  	  
C --- Damage dissipation calculation (map 4 unique Y_THERMO to 6 components)
      w_diss_damage = 0.0D0
      w_diss_damage = w_diss_damage + Y_THERMO(1) * (STATEV(33) - DAMAGE_OLD(1))
      w_diss_damage = w_diss_damage + Y_THERMO(2) * (STATEV(34) - DAMAGE_OLD(2)) * 2  ! 22 and 33
      w_diss_damage = w_diss_damage + Y_THERMO(3) * (STATEV(35) - DAMAGE_OLD(3)) * 2  ! 12 and 13
      w_diss_damage = w_diss_damage + Y_THERMO(4) * (STATEV(36) - DAMAGE_OLD(4))      ! 23

C --- Total dissipated energy and temperature rise
      w_diss_total = w_diss_visco + w_diss_damage
      STATEV(39) = STATEV(39) + w_diss_total
      delta_T_total = w_diss_total / (rho * cp)*1000000
      STATEV(40) = STATEV(40) + delta_T_total
	  STATEV(37) = STATEV(38)
	  RPL=w_diss_total
	  STATEV(44)=STATEV(40)+DTEMP
      RETURN
      END