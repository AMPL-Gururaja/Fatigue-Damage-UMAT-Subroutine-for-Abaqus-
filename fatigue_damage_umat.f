!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CF-ABS COUPLED FATIGUE DAMAGE + THERMAL SUBROUTINES
!!
!!  Subroutines:   UMAT  (mechanical + fatigue damage)
!!                 UMATHT (thermal constitutive — damage-dep. conductivity)
!!                 CALC_DAMAGE_RATE (three-phase CDM helper)
!!
!!  Element type:  C3D8T  (fully coupled temperature-displacement)
!!  Procedure:     *COUPLED TEMPERATURE-DISPLACEMENT, DELTMX=<value>
!!  Abaqus ver.:   2019 and later
!!
!! ── PROPS LAYOUT (NPROPS = 50) ──────────────────────────────────────────
!!  #   Name          Example value   Units            Notes
!!  1   E11_0         9605.0          MPa              fibre-dir. modulus
!!  2   E22_0         4500.0          MPa              transverse modulus
!!  3   G12_0         3500.0          MPa              in-plane shear
!!  4   G23_0         2800.0          MPa              transverse shear
!!  5   NU12          0.270           —
!!  6   NU23          0.350           —
!!  7   RHO           1.15e-9         tonne/mm³        density
!!  8   CP            1.40e6          mJ/(tonne·K)     specific heat
!!  9   PRONY_1       0.010           1/s              Prony relax rate 1
!!  10  PRONY_2       0.001           1/s              Prony relax rate 2
!!  11  PRONY_3       0.0001          1/s              Prony relax rate 3
!!  12  D_PRONY_1     0.08            —                Prony amplitude 1
!!  13  D_PRONY_2     0.15            —                Prony amplitude 2
!!  14  D_PRONY_3     0.05            —                Prony amplitude 3
!!  15  GAMMA_D11     8.0e-3          —                CDM Phase-II coeff D11
!!  16  Q_D11         4.5482          —                CDM Paris exponent D11
!!  17  LAMBDA_D11    3.9e-1          —                CDM Phase-I coeff D11
!!  18  DELTA_D11     9.2858          —                CDM decay rate D11
!!  19  Y_TH_D11      1.2498e-3       MPa              ERR threshold D11
!!  20  GAMMA_D22     5.1052e-3       —                CDM D22
!!  21  Q_D22         4.5482          —
!!  22  LAMBDA_D22    4.1086e-2       —
!!  23  DELTA_D22     1.1776e-1       —
!!  24  Y_TH_D22      1.2498e-3       MPa
!!  25  GAMMA_D12     7.0e-3          —                CDM D12
!!  26  Q_D12         4.5482          —
!!  27  LAMBDA_D12    5.0e-2          —
!!  28  DELTA_D12     1.1776e-1       —
!!  29  Y_TH_D12      1.2498e-3       MPa
!!  30  GAMMA_D23     4.0e-3          —                CDM D23
!!  31  Q_D23         4.5482          —
!!  32  LAMBDA_D23    3.5e-2          —
!!  33  DELTA_D23     1.1776e-1       —
!!  34  Y_TH_D23      1.2498e-3       MPa
!!  35  FREQUENCY     15.0            Hz               cyclic loading freq.
!!  36  ALPHA_E       0.002           1/K              thermal softening E
!!  37  ALPHA_G       0.003           1/K              thermal softening G
!!  38  T_REF         298.0           K                stress-free ref. temp
!!  39  NCYC_JUMP     500.0           cycles           cycle-jump size
!!  40  NCYC_MIN      1000.0          cycles           min cycles before jump
!!  41  JUMP_TOL      0.005           —                max dD per step for jump
!!  42  D_F           0.07153         —                failure damage threshold
!!  43  K_TERT        2.90            —                Phase-III exponent
!!  44  P_GAMMA       0.1513          —                Y-scaling exponent
!!  45  P_LAMBDA      2.6091          —
!!  46  P_DELTA       3.1844          —
!!  47  K_FIBER_AX    17.0            N/(s·K)          fibre axial conductivity
!!  48  K_FIBER_TR    7.0             N/(s·K)          fibre transverse conduct.
!!  49  K_MATRIX      0.17            N/(s·K)          ABS matrix conductivity
!!  50  V_FIBER       0.13            —                fibre volume fraction
!!  51  HV_VISCO      (stress-dep.)   MPa/cycle        calibrated viscous heat
!!  52  BETA_CAL      (stress-dep.)   /cycle           convection decay
!!  53  SLOPE_CAL     (stress-dep.)   K/cycle          late-phase drift
!!  STATEV(49) = T_CAL  step-local temperature rise [K]  <- PLOT THIS
!!
!!  Unit note: 1 N/(s·K) ≡ 1 W/(m·K) in the Abaqus mm-tonne-s-N system.
!!
!! ── STATEV MAP (NSTATV = 46) ────────────────────────────────────────────
!!  1-6   strain_ve(6)         cumulative VE strain vector
!!  7-12  strain_el(6)         elastic strain vector
!!  13    W_DISS_VISCO         cumul. viscous dissipation [MPa]
!!  14    DELTA_T_VISCO        cumul. viscous self-heating [K]
!!  15-32 strain_ve_elem(3,6)  per-Prony VE strains
!!  33    D11     ← SDV33 : plot vs SDV37 in Abaqus Viewer
!!  34    D22
!!  35    D12
!!  36    D23
!!  37    N_CYCLE ← SDV37 : cycle counter (x-axis)
!!  38    W_DISS_TOTAL         cumul. total dissipation [MPa]
!!  39    DELTA_T_TOTAL        cumul. total self-heating [K]
!!  40    TEMP_CURRENT         current nodal temperature [K]
!!  41    W_DISS_DAMAGE        cumul. damage dissipation [MPa]
!!  42-43 (reserved)
!!  44    DELTA_T_DAMAGE       3-phase temp: SDV44 vs SDV37  [K]
!!  45    KSTEP_OLD
!!  46    RATE_D11             dD11/dN at current increment
!!  47    SIGMA_PEAK_CYCLE     peak |sigma_11| in current cycle (resets each cycle)
!!  48    N_FLOOR_OLD          last integer cycle boundary (internal)
!!
!! ── ABAQUS INPUT FILE TEMPLATE ──────────────────────────────────────────
!!  *MATERIAL, NAME=CF_ABS_FATIGUE
!!  *USER MATERIAL, CONSTANTS=50, TYPE=THERMOMECHANICAL
!!  9605.0, 4500.0, 3500.0, 2800.0, 0.270, 0.350, 1.15E-9, 1.40E6,
!!  0.010,  0.001,  0.0001, 0.08,   0.15,  0.05,
!!  8.0E-3, 4.5482, 3.9E-1, 9.2858, 1.2498E-3,
!!  5.1052E-3, 4.5482, 4.1086E-2, 1.1776E-1, 1.2498E-3,
!!  7.0E-3, 4.5482, 5.0E-2, 1.1776E-1, 1.2498E-3,
!!  4.0E-3, 4.5482, 3.5E-2, 1.1776E-1, 1.2498E-3,
!!  15.0, 0.002, 0.003, 298.0, 500.0, 1000.0, 0.005,
!!  0.07153, 2.90, 0.1513, 2.6091, 3.1844,
!!  17.0, 7.0, 0.17, 0.13
!!  *DEPVAR
!!  46
!!  *INITIAL CONDITIONS, TYPE=TEMPERATURE
!!  ALL_NODES, 298.0
!!  *STEP, INC=500000, NLGEOM=NO
!!  *COUPLED TEMPERATURE-DISPLACEMENT, DELTMX=5.0
!!  <total_time>, <total_time>, 1.0E-6, <max_time_inc>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!============================================================================
!  UMAT — mechanical + fatigue damage constitutive subroutine
!  Modified from cf_abs_1elem_umat.f:
!    (1) Thermal softening f_E(T), f_G(T) applied to all moduli
!    (2) TEMP_NOW = TEMP + DTEMP used in modulus degradation
!    (3) PROPS extended to 50 entries (thermal props at 47-50)
!============================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C --- Working arrays
      REAL*8 STRAN_TOTAL(6), STRAIN_VE(6), STRAIN_EL(6)
      REAL*8 STRAIN_VE_ELEM(3,6), STRAIN_VE_ELEM_OLD(3,6)
      REAL*8 PRONY(3), D_PRONY(3)
      REAL*8 EXP_PRONY(3), ONE_MINUS_EXP(3)
      REAL*8 MODULUS_VEC(6), Y_THERMO(6)
      REAL*8 DAMAGE_OLD(4), DAMAGE_RATE(4), RATE_JUMP(4), D_NEW(4)
      REAL*8 GAMMA_D(4), Q_D(4), LAMBDA_D(4), DELTA_D(4), Y_TH_D(4)

C --- Scalars
      REAL*8 E11_0, E22_0, G12_0, G23_0, G13_0
      REAL*8 NU12, NU23, NU13, NU21, NU31, NU32
      REAL*8 RHO, CP, FREQUENCY
      REAL*8 ALPHA_E, ALPHA_G, T_REF, TEMP_NOW, F_E, F_G
      REAL*8 NCYC_JUMP, NCYC_MIN, JUMP_TOL, K_TERT, D_F
      REAL*8 P_GAMMA, P_LAMBDA, P_DELTA
      REAL*8 E11, E22, E33, G12, G13, G23
      REAL*8 DELTA_DENOM, INV_DENOM
      REAL*8 W_DISS_VISCO, W_DISS_DAMAGE, W_DISS_TOTAL
      REAL*8 HV_VISCO, BETA_CAL, SLOPE_CAL, N_LOCAL_CAF, T_CAL_CAF
      REAL*8 N_CYCLE, N_CYCLE_STEP, N_JUMP, N_CYC_TOTAL
      REAL*8 DD_MAX, DD, MODULUS, SIGMA_OLD, SIGMA_NEW, DSTRAIN_VE_IJ
      REAL*8 TERM_N
      REAL*8 SIGMA_PEAK, N_FLOOR_OLD, N_FLOOR_NEW, DN_CYC
      REAL*8 Y_PEAK(4)
      INTEGER I, J, IDX, I_COMP, K1, K2

C --- STATEV indices
      INTEGER ISV_VE_STRAIN, ISV_EL_STRAIN
      INTEGER ISV_WVISCO, ISV_DT_VISCO, ISV_VE_ELEM
      INTEGER ISV_D11, ISV_NCYC, ISV_WDISS, ISV_DT_TOTAL
      INTEGER ISV_TEMP_CURRENT, ISV_WDAMAGE, ISV_DT_DAMAGE
      INTEGER ISV_KSTEP, ISV_RATE_D11

      PARAMETER (ISV_VE_STRAIN    = 1 )
      PARAMETER (ISV_EL_STRAIN    = 7 )
      PARAMETER (ISV_WVISCO       = 13)
      PARAMETER (ISV_DT_VISCO     = 14)
      PARAMETER (ISV_VE_ELEM      = 15)
      PARAMETER (ISV_D11          = 33)
      PARAMETER (ISV_NCYC         = 37)
      PARAMETER (ISV_WDISS        = 38)
      PARAMETER (ISV_DT_TOTAL     = 39)
      PARAMETER (ISV_TEMP_CURRENT = 40)
      PARAMETER (ISV_WDAMAGE      = 41)
      PARAMETER (ISV_DT_DAMAGE    = 44)
      PARAMETER (ISV_KSTEP        = 45)
      PARAMETER (ISV_RATE_D11     = 46)
      PARAMETER (ISV_SIGMA_PEAK   = 47)
      PARAMETER (ISV_NFLOOR_D     = 48)

      REAL*8 ZERO, ONE, TWO, HALF, TOLER
      PARAMETER (ZERO  = 0.0D0)
      PARAMETER (ONE   = 1.0D0)
      PARAMETER (TWO   = 2.0D0)
      PARAMETER (HALF  = 0.5D0)
      PARAMETER (TOLER = 1.0D-12)

!======================================================================
!  STEP 1a — Read material properties
!======================================================================
      E11_0      = PROPS(1)
      E22_0      = PROPS(2)
      G12_0      = PROPS(3)
      G23_0      = PROPS(4)
      NU12       = PROPS(5)
      NU23       = PROPS(6)
      RHO        = PROPS(7)
      CP         = PROPS(8)
      PRONY(1)   = PROPS(9)
      PRONY(2)   = PROPS(10)
      PRONY(3)   = PROPS(11)
      D_PRONY(1) = PROPS(12)
      D_PRONY(2) = PROPS(13)
      D_PRONY(3) = PROPS(14)

      GAMMA_D(1)  = PROPS(15)
      Q_D(1)      = PROPS(16)
      LAMBDA_D(1) = PROPS(17)
      DELTA_D(1)  = PROPS(18)
      Y_TH_D(1)   = PROPS(19)
      GAMMA_D(2)  = PROPS(20)
      Q_D(2)      = PROPS(21)
      LAMBDA_D(2) = PROPS(22)
      DELTA_D(2)  = PROPS(23)
      Y_TH_D(2)   = PROPS(24)
      GAMMA_D(3)  = PROPS(25)
      Q_D(3)      = PROPS(26)
      LAMBDA_D(3) = PROPS(27)
      DELTA_D(3)  = PROPS(28)
      Y_TH_D(3)   = PROPS(29)
      GAMMA_D(4)  = PROPS(30)
      Q_D(4)      = PROPS(31)
      LAMBDA_D(4) = PROPS(32)
      DELTA_D(4)  = PROPS(33)
      Y_TH_D(4)   = PROPS(34)

      FREQUENCY = PROPS(35)
      ALPHA_E   = PROPS(36)
      ALPHA_G   = PROPS(37)
      T_REF     = PROPS(38)
      NCYC_JUMP = PROPS(39)
      NCYC_MIN  = PROPS(40)
      JUMP_TOL  = PROPS(41)
      D_F       = PROPS(42)
      K_TERT    = PROPS(43)
      P_GAMMA   = PROPS(44)
      P_LAMBDA  = PROPS(45)
      P_DELTA   = PROPS(46)
C     Props 47-50 are thermal (read by UMATHT only)

      NU13  = NU12
      G13_0 = G12_0

C --- Zero coupled Jacobian arrays
      DO I = 1, NTENS
         DDSDDT(I) = ZERO
      END DO
      DRPLDT = ZERO

!======================================================================
!  STEP 1b — Read state variables from previous increment
!======================================================================
      DO I_COMP = 1, 4
         DAMAGE_OLD(I_COMP) = STATEV(ISV_D11 + I_COMP - 1)
         IF (DAMAGE_OLD(I_COMP) .GT. D_F)
     1      DAMAGE_OLD(I_COMP) = D_F
      END DO

      DO I = 1, 3
         DO J = 1, 6
            IDX = ISV_VE_ELEM + (I-1)*6 + (J-1)
            STRAIN_VE_ELEM_OLD(I,J) = STATEV(IDX)
         END DO
      END DO

      N_CYCLE      = STATEV(ISV_NCYC)
      SIGMA_PEAK   = STATEV(ISV_SIGMA_PEAK)
      N_FLOOR_OLD  = STATEV(ISV_NFLOOR_D)

!======================================================================
!  STEP 2 — Thermal softening functions f_E(T) and f_G(T)
!  Activated by non-zero ALPHA_E / ALPHA_G in PROPS(36/37)
!  TEMP_NOW = temperature at end of increment (from Abaqus thermal pass)
!======================================================================
      TEMP_NOW = TEMP + DTEMP

      IF (ABS(ALPHA_E) .GT. TOLER) THEN
         F_E = EXP(-ALPHA_E * (TEMP_NOW - T_REF))
         IF (F_E .GT. ONE)  F_E = ONE
         IF (F_E .LT. 0.1D0) F_E = 0.1D0
      ELSE
         F_E = ONE
      END IF

      IF (ABS(ALPHA_G) .GT. TOLER) THEN
         F_G = EXP(-ALPHA_G * (TEMP_NOW - T_REF))
         IF (F_G .GT. ONE)  F_G = ONE
         IF (F_G .LT. 0.1D0) F_G = 0.1D0
      ELSE
         F_G = ONE
      END IF

!======================================================================
!  STEP 3 — Damage- and temperature-degraded elastic moduli
!======================================================================
      IF (DAMAGE_OLD(1) .GE. D_F) THEN
         E11 = 0.1D0 * E11_0 * F_E
      ELSE
         E11 = E11_0 * F_E * (ONE - DAMAGE_OLD(1))
      END IF
      IF (DAMAGE_OLD(2) .GE. D_F) THEN
         E22 = 0.1D0 * E22_0 * F_E
      ELSE
         E22 = E22_0 * F_E * (ONE - DAMAGE_OLD(2))
      END IF
      IF (DAMAGE_OLD(3) .GE. D_F) THEN
         G12 = 0.1D0 * G12_0 * F_G
      ELSE
         G12 = G12_0 * F_G * (ONE - DAMAGE_OLD(3))
      END IF
      IF (DAMAGE_OLD(4) .GE. D_F) THEN
         G23 = 0.1D0 * G23_0 * F_G
      ELSE
         G23 = G23_0 * F_G * (ONE - DAMAGE_OLD(4))
      END IF

      E33 = E22;  G13 = G12
      NU21 = NU12 * E22 / E11
      NU31 = NU13 * E33 / E11
      NU32 = NU23

      DELTA_DENOM = ONE - NU12*NU21 - NU23*NU32
     1            - NU13*NU31 - TWO*NU21*NU32*NU13
      IF (ABS(DELTA_DENOM) .LT. TOLER) DELTA_DENOM = TOLER
      INV_DENOM = ONE / DELTA_DENOM

!======================================================================
!  STEP 4 — Tangent stiffness DDSDDE
!======================================================================
      DO I = 1, NTENS
         DO J = 1, NTENS
            DDSDDE(I,J) = ZERO
         END DO
      END DO
      DDSDDE(1,1) = E11 * (ONE  - NU23*NU32) * INV_DENOM
      DDSDDE(2,2) = E22 * (ONE  - NU13*NU31) * INV_DENOM
      DDSDDE(3,3) = E33 * (ONE  - NU12*NU21) * INV_DENOM
      DDSDDE(1,2) = E11 * (NU21 + NU31*NU23) * INV_DENOM
      DDSDDE(1,3) = E11 * (NU31 + NU21*NU32) * INV_DENOM
      DDSDDE(2,3) = E22 * (NU32 + NU12*NU31) * INV_DENOM
      DDSDDE(2,1) = DDSDDE(1,2)
      DDSDDE(3,1) = DDSDDE(1,3)
      DDSDDE(3,2) = DDSDDE(2,3)
      DDSDDE(4,4) = G12
      DDSDDE(5,5) = G13
      DDSDDE(6,6) = G23

!======================================================================
!  STEP 5 — Strain decomposition and viscoelastic update
!======================================================================
      DO I = 1, 6
         STRAN_TOTAL(I) = STRAN(I) + DSTRAN(I)
      END DO

      MODULUS_VEC(1) = E11_0;  MODULUS_VEC(2) = E22_0
      MODULUS_VEC(3) = E22_0;  MODULUS_VEC(4) = G12_0
      MODULUS_VEC(5) = G13_0;  MODULUS_VEC(6) = G23_0

C --- VE Pass 1: initial elastic strain from old Prony strains
      DO I = 1, 6
         STRAIN_VE(I) = ZERO
         DO J = 1, 3
            STRAIN_VE(I) = STRAIN_VE(I) + STRAIN_VE_ELEM_OLD(J,I)
         END DO
         STRAIN_EL(I) = STRAN_TOTAL(I) - STRAIN_VE(I)
      END DO

      DO I = 1, 3
         EXP_PRONY(I)     = EXP(-PRONY(I) * DTIME)
         ONE_MINUS_EXP(I) = ONE - EXP_PRONY(I)
      END DO

C --- Update Prony element strains; accumulate viscous dissipation
      W_DISS_VISCO = ZERO
      DO J = 1, 6
         MODULUS = MODULUS_VEC(J)
         DO I = 1, 3
            STRAIN_VE_ELEM(I,J) =
     1           EXP_PRONY(I)     * STRAIN_VE_ELEM_OLD(I,J)
     2         + ONE_MINUS_EXP(I) * STRAIN_EL(J)
            DSTRAIN_VE_IJ = STRAIN_VE_ELEM(I,J)
     1                    - STRAIN_VE_ELEM_OLD(I,J)
            SIGMA_OLD = MODULUS * D_PRONY(I) * STRAIN_VE_ELEM_OLD(I,J)
            SIGMA_NEW = MODULUS * D_PRONY(I) * STRAIN_VE_ELEM(I,J)
            W_DISS_VISCO = W_DISS_VISCO
     1                   + HALF * (SIGMA_OLD + SIGMA_NEW) * DSTRAIN_VE_IJ
         END DO
      END DO

C --- VE Pass 2: updated elastic strain from new Prony strains
      DO J = 1, 6
         STRAIN_VE(J) = ZERO
         DO I = 1, 3
            STRAIN_VE(J) = STRAIN_VE(J) + STRAIN_VE_ELEM(I,J)
         END DO
         STRAIN_EL(J) = STRAN_TOTAL(J) - STRAIN_VE(J)
      END DO

!======================================================================
!  STEP 6 — Cauchy stress
!======================================================================
      DO K1 = 1, NTENS
         STRESS(K1) = ZERO
         DO K2 = 1, NTENS
            STRESS(K1) = STRESS(K1) + DDSDDE(K1,K2) * STRAIN_EL(K2)
         END DO
      END DO

C --- Track running peak axial stress within current cycle
      SIGMA_PEAK = MAX(SIGMA_PEAK, ABS(STRESS(1)))

!======================================================================
!  STEP 7 — Thermodynamic energy release rates  Y_k
!  Normal components (k=1,2,3): tension gate applied
!  Shear components (k=4,6):    always active
!======================================================================
      TERM_N = STRAIN_EL(1) + NU21*STRAIN_EL(2) + NU31*STRAIN_EL(3)
      IF (TERM_N .GT. ZERO .AND. STRAIN_EL(1) .GT. ZERO) THEN
         Y_THERMO(1) = HALF * E11 * STRAIN_EL(1) * TERM_N * INV_DENOM
      ELSE
         Y_THERMO(1) = ZERO
      END IF

      TERM_N = STRAIN_EL(2) + NU12*STRAIN_EL(1) + NU32*STRAIN_EL(3)
      IF (TERM_N .GT. ZERO .AND. STRAIN_EL(2) .GT. ZERO) THEN
         Y_THERMO(2) = HALF * E22 * STRAIN_EL(2) * TERM_N * INV_DENOM
      ELSE
         Y_THERMO(2) = ZERO
      END IF

      TERM_N = STRAIN_EL(3) + NU13*STRAIN_EL(1) + NU23*STRAIN_EL(2)
      IF (TERM_N .GT. ZERO .AND. STRAIN_EL(3) .GT. ZERO) THEN
         Y_THERMO(3) = HALF * E33 * STRAIN_EL(3) * TERM_N * INV_DENOM
      ELSE
         Y_THERMO(3) = ZERO
      END IF

      Y_THERMO(4) = HALF * G12 * STRAIN_EL(4)**2
      Y_THERMO(5) = HALF * G13 * STRAIN_EL(5)**2
      Y_THERMO(6) = HALF * G23 * STRAIN_EL(6)**2

!======================================================================
!  STEP 8 — Cycle counter advance
!======================================================================
      N_CYCLE_STEP = DTIME * FREQUENCY
      N_CYCLE = N_CYCLE + N_CYCLE_STEP

!======================================================================
!  STEP 8b — Cycle boundary and peak-Y for damage law
!  Damage fires ONCE per completed cycle using peak-cycle sigma, so
!  Y is consistent with the calibration (which used sigma_max not
!  time-averaged instantaneous sigma). Without this, Phase II is
!  ~8x too slow because instantaneous Y < Y_TH for 66% of each cycle.
!======================================================================
      N_FLOOR_NEW = AINT(N_CYCLE)
      DN_CYC      = N_FLOOR_NEW - N_FLOOR_OLD

      IF (DN_CYC .GT. TOLER) THEN

C     Compute peak-cycle Y for all four CDM channels
         Y_PEAK(1) = HALF * SIGMA_PEAK**2 / E11_0
         Y_PEAK(2) = Y_PEAK(1)
         Y_PEAK(3) = Y_PEAK(1)
         Y_PEAK(4) = Y_PEAK(1)

!======================================================================
!  STEP 9 — Three-phase CDM damage rates + forward-Euler update
!  Y_PEAK(k) = sigma_peak^2/(2*E0) consistent with calibration
!======================================================================
         CALL CALC_DAMAGE_RATE(Y_PEAK(1), DAMAGE_OLD(1), Q_D(1),
     1      GAMMA_D(1), LAMBDA_D(1), DELTA_D(1), Y_TH_D(1), K_TERT, D_F,
     2      P_GAMMA, P_LAMBDA, P_DELTA, N_CYCLE, TOLER, DAMAGE_RATE(1))

         CALL CALC_DAMAGE_RATE(Y_PEAK(2), DAMAGE_OLD(2), Q_D(2),
     1      GAMMA_D(2), LAMBDA_D(2), DELTA_D(2), Y_TH_D(2), K_TERT, D_F,
     2      P_GAMMA, P_LAMBDA, P_DELTA, N_CYCLE, TOLER, DAMAGE_RATE(2))

         CALL CALC_DAMAGE_RATE(Y_PEAK(3), DAMAGE_OLD(3), Q_D(3),
     1      GAMMA_D(3), LAMBDA_D(3), DELTA_D(3), Y_TH_D(3), K_TERT, D_F,
     2      P_GAMMA, P_LAMBDA, P_DELTA, N_CYCLE, TOLER, DAMAGE_RATE(3))

         CALL CALC_DAMAGE_RATE(Y_PEAK(4), DAMAGE_OLD(4), Q_D(4),
     1      GAMMA_D(4), LAMBDA_D(4), DELTA_D(4), Y_TH_D(4), K_TERT, D_F,
     2      P_GAMMA, P_LAMBDA, P_DELTA, N_CYCLE, TOLER, DAMAGE_RATE(4))

         DO I_COMP = 1, 4
            D_NEW(I_COMP) = DAMAGE_OLD(I_COMP)
     1                    + DAMAGE_RATE(I_COMP) * DN_CYC
         END DO

!======================================================================
!  STEP 10 — Monotonicity + failure caps
!======================================================================
         DO I_COMP = 1, 4
            IF (D_NEW(I_COMP) .GT. D_F) D_NEW(I_COMP) = D_F
            IF (D_NEW(I_COMP) .LT. DAMAGE_OLD(I_COMP))
     1         D_NEW(I_COMP) = DAMAGE_OLD(I_COMP)
         END DO

!======================================================================
!  STEP 11 — Cycle-jump acceleration
!  Fires when: DD_max < JUMP_TOL  AND  N > NCYC_MIN  AND  NCYC_JUMP > 0
!======================================================================
      N_JUMP = ZERO
      DD_MAX = ZERO
      DO I_COMP = 1, 4
         DD = D_NEW(I_COMP) - DAMAGE_OLD(I_COMP)
         IF (DD .GT. DD_MAX) DD_MAX = DD
      END DO

         IF (DD_MAX .LT. JUMP_TOL .AND. N_CYCLE .GT. NCYC_MIN .AND.
     1       NCYC_JUMP .GT. TOLER) THEN

            N_JUMP  = NCYC_JUMP
            N_CYCLE = N_CYCLE + N_JUMP

            CALL CALC_DAMAGE_RATE(Y_PEAK(1), D_NEW(1), Q_D(1),
     1         GAMMA_D(1), LAMBDA_D(1), DELTA_D(1), Y_TH_D(1), K_TERT,
     2         D_F, P_GAMMA, P_LAMBDA, P_DELTA, N_CYCLE, TOLER, RATE_JUMP(1))
            CALL CALC_DAMAGE_RATE(Y_PEAK(2), D_NEW(2), Q_D(2),
     1         GAMMA_D(2), LAMBDA_D(2), DELTA_D(2), Y_TH_D(2), K_TERT,
     2         D_F, P_GAMMA, P_LAMBDA, P_DELTA, N_CYCLE, TOLER, RATE_JUMP(2))
            CALL CALC_DAMAGE_RATE(Y_PEAK(3), D_NEW(3), Q_D(3),
     1         GAMMA_D(3), LAMBDA_D(3), DELTA_D(3), Y_TH_D(3), K_TERT,
     2         D_F, P_GAMMA, P_LAMBDA, P_DELTA, N_CYCLE, TOLER, RATE_JUMP(3))
            CALL CALC_DAMAGE_RATE(Y_PEAK(4), D_NEW(4), Q_D(4),
     1         GAMMA_D(4), LAMBDA_D(4), DELTA_D(4), Y_TH_D(4), K_TERT,
     2         D_F, P_GAMMA, P_LAMBDA, P_DELTA, N_CYCLE, TOLER, RATE_JUMP(4))

            DO I_COMP = 1, 4
               DD = D_NEW(I_COMP) + RATE_JUMP(I_COMP) * N_JUMP
               IF (DD .GT. D_F)           DD = D_F
               IF (DD .LT. D_NEW(I_COMP)) DD = D_NEW(I_COMP)
               D_NEW(I_COMP) = DD
            END DO
         END IF


C     --- Reset per-cycle peak accumulator after damage update ---
         SIGMA_PEAK = ZERO

      ELSE
C     Mid-cycle: no damage update; carry D_NEW = D_OLD
         DO I_COMP = 1, 4
            D_NEW(I_COMP) = DAMAGE_OLD(I_COMP)
         END DO
         N_JUMP = ZERO
      END IF

!======================================================================
!  STEP 12 — Dissipation energies + RPL heat source
!======================================================================
      N_CYC_TOTAL = N_CYCLE_STEP + N_JUMP

      W_DISS_DAMAGE =
     1      Y_THERMO(1) * (D_NEW(1) - DAMAGE_OLD(1))
     2    + Y_THERMO(2) * (D_NEW(2) - DAMAGE_OLD(2)) * TWO
     3    + Y_THERMO(4) * (D_NEW(3) - DAMAGE_OLD(3)) * TWO
     4    + Y_THERMO(6) * (D_NEW(4) - DAMAGE_OLD(4))

      W_DISS_TOTAL = W_DISS_VISCO + W_DISS_DAMAGE

C --- Read calibrated temperature parameters
      HV_VISCO  = PROPS(51)   ! viscous heat [MPa/cycle] from IRT
      BETA_CAL  = PROPS(52)   ! convection decay [/cycle]
      SLOPE_CAL = PROPS(53)   ! late-phase slope [K/cycle]

C --- T_CAL: calibrated per-step temperature using TIME(1) (step time)
C --- TIME(1) resets to 0 at each new *STEP -> T_CAL resets automatically
      N_LOCAL_CAF = TIME(1) * FREQUENCY
      T_CAL_CAF = MAX(ZERO,
     1    (PROPS(51)/(BETA_CAL*RHO*CP))
     2    *(ONE-EXP(-BETA_CAL*N_LOCAL_CAF))
     3    + SLOPE_CAL*N_LOCAL_CAF)

C --- RPL: replace near-zero Prony visco with calibrated Hv + damage heat
      IF (N_CYCLE_STEP .GT. TOLER) THEN
         RPL = (HV_VISCO
     1        + W_DISS_DAMAGE / MAX(N_CYC_TOTAL, TOLER)) * FREQUENCY
      ELSE
         RPL = ZERO
      END IF
      DO I = 1, NTENS
         DRPLDE(I) = STRESS(I) * FREQUENCY
      END DO

!======================================================================
!  STEP 13 — Write all STATEV for next increment
!======================================================================
      DO J = 1, 6
         STATEV(ISV_VE_STRAIN + J - 1) = STRAIN_VE(J)
         STATEV(ISV_EL_STRAIN + J - 1) = STRAIN_EL(J)
      END DO

      STATEV(ISV_WVISCO)   = STATEV(ISV_WVISCO)   + W_DISS_VISCO
      STATEV(ISV_DT_VISCO) = STATEV(ISV_DT_VISCO)
     1                     + W_DISS_VISCO / (RHO * CP)

      DO I = 1, 3
         DO J = 1, 6
            IDX = ISV_VE_ELEM + (I-1)*6 + (J-1)
            STATEV(IDX) = STRAIN_VE_ELEM(I,J)
         END DO
      END DO

      DO I_COMP = 1, 4
         STATEV(ISV_D11 + I_COMP - 1) = D_NEW(I_COMP)
      END DO

      STATEV(ISV_NCYC) = N_CYCLE

      STATEV(ISV_WDISS)    = STATEV(ISV_WDISS)    + W_DISS_TOTAL
      STATEV(ISV_DT_TOTAL) = STATEV(ISV_DT_TOTAL)
     1                     + W_DISS_TOTAL / (RHO * CP)

      STATEV(ISV_TEMP_CURRENT) = TEMP_NOW

      STATEV(ISV_WDAMAGE) = STATEV(ISV_WDAMAGE) + W_DISS_DAMAGE
      STATEV(ISV_DT_DAMAGE) = STATEV(ISV_WDAMAGE) / (RHO * CP)

      STATEV(ISV_KSTEP)      = DBLE(KSTEP)
      STATEV(ISV_RATE_D11)   = DAMAGE_RATE(1)
      STATEV(ISV_SIGMA_PEAK) = SIGMA_PEAK
      STATEV(ISV_NFLOOR_D)   = N_FLOOR_NEW
      STATEV(49)             = T_CAL_CAF  ! step-local T rise

      RETURN
      END


!============================================================================
!  UMATHT — thermal constitutive subroutine
!  Computes: damage-dependent anisotropic conductivity k(D)
!            heat flux  FLUX = -k(D) . grad(T)
!            specific internal energy U = cp * T
!            heat capacity DUDT = cp
!
!  PROPS read here (same array as UMAT):
!    PROPS(7)  = RHO     [tonne/mm³]
!    PROPS(8)  = CP      [mJ/(tonne·K)]
!    PROPS(47) = K_FIBER_AXIAL  [N/(s·K) = W/(m·K) in SI]
!    PROPS(48) = K_FIBER_TRANS  [N/(s·K)]
!    PROPS(49) = K_MATRIX       [N/(s·K)]
!    PROPS(50) = V_FIBER        [—]
!
!  STATEV read here (shared with UMAT, do NOT write in UMATHT):
!    STATEV(33) = D11  (axial damage)
!    STATEV(34) = D22  (transverse damage)
!
!  Physical model:
!    Axial direction (1 = fibre direction):
!      k11(D) = k_f_ax * V_f + k_m * (1-V_f) * (1 - D11)
!               Rule of mixtures; matrix contribution degrades with D11
!    Transverse (2 and 3 = transverse isotropy plane):
!      k_HS = 2*k_m*k_f_tr / (k_m + k_f_tr)   [Hashin harmonic mean]
!      k22 = k33 = k_HS * (1 - D22)
!
!  Note on sign convention (Abaqus):
!    FLUX(I) = -k(I,I) * DTEMDX(I)   (Fourier's law, isotropic per direction)
!    DFLUXDG(I,J) = d FLUX(I) / d DTEMDX(J)  -->  DFLUXDG(I,I) = -k(I,I)
!============================================================================
      SUBROUTINE UMATHT(U, DUDT, DUDG, FLUX, DFLUDT, DFLUXDG,
     1 STATEV, TEMP, DTEMP, DTEMDX, TIME, DTIME,
     2 PREDEF, DPRED, CMNAME, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, PNEWDT, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME

      DIMENSION FLUX(NTENS), DFLUDT(NTENS), DFLUXDG(NTENS,NTENS)
      DIMENSION DUDG(NTENS), DTEMDX(NTENS), STATEV(NSTATV)
      DIMENSION PREDEF(1), DPRED(1), PROPS(NPROPS)
      DIMENSION COORDS(3), TIME(2)

C --- Local variables
      REAL*8 RHO, CP, K_FAX, K_FTR, K_MAT, V_F
      REAL*8 D11, D22, D_AVG
      REAL*8 K11, K22, K33, K_HASHIN
      REAL*8 TEMP_NOW
      REAL*8 ZERO, ONE, TWO, TOLER, R_FLOOR
      INTEGER I, J

      PARAMETER (ZERO   = 0.0D0)
      PARAMETER (ONE    = 1.0D0)
      PARAMETER (TWO    = 2.0D0)
      PARAMETER (TOLER  = 1.0D-12)
      PARAMETER (R_FLOOR = 0.05D0)

!----------------------------------------------------------------------
!  Read material constants from PROPS
!----------------------------------------------------------------------
      RHO   = PROPS(7)
      CP    = PROPS(8)
      K_FAX = PROPS(47)
      K_FTR = PROPS(48)
      K_MAT = PROPS(49)
      V_F   = PROPS(50)

!----------------------------------------------------------------------
!  Read current damage from STATEV (written by UMAT, do not modify here)
!----------------------------------------------------------------------
      D11 = STATEV(33)
      D22 = STATEV(34)

C --- Safety floor: damage can be negative initially (first call = 0)
      IF (D11 .LT. ZERO) D11 = ZERO
      IF (D22 .LT. ZERO) D22 = ZERO

C --- Average damage for heat capacity reduction (optional, mild effect)
      D_AVG = 0.5D0 * (D11 + D22)
      IF (D_AVG .GT. 0.5D0) D_AVG = 0.5D0

!----------------------------------------------------------------------
!  Damage-dependent anisotropic thermal conductivity
!----------------------------------------------------------------------
C --- Axial conductivity (rule of mixtures, matrix degrades with D11)
      K11 = K_FAX * V_F
     1    + K_MAT * (ONE - V_F) * MAX(ONE - D11, R_FLOOR)

C --- Transverse conductivity (Hashin harmonic mean, degrades with D22)
      IF ((K_MAT + K_FTR) .GT. TOLER) THEN
         K_HASHIN = TWO * K_MAT * K_FTR / (K_MAT + K_FTR)
      ELSE
         K_HASHIN = K_MAT
      END IF
      K22 = K_HASHIN * MAX(ONE - D22, R_FLOOR)
      K33 = K22

!----------------------------------------------------------------------
!  Temperature at end of increment
!----------------------------------------------------------------------
      TEMP_NOW = TEMP + DTEMP

!----------------------------------------------------------------------
!  Specific internal energy per unit mass [mJ/tonne]
!  U = cp * T  (simple linear model; absolute T in Kelvin)
!----------------------------------------------------------------------
      U    = CP * TEMP_NOW
      DUDT = CP * (ONE - 0.1D0 * D_AVG)

C --- dU/d(grad T) = 0 (standard assumption)
      DO I = 1, NTENS
         DUDG(I) = ZERO
      END DO

!----------------------------------------------------------------------
!  Heat flux vector:  FLUX(I) = -K(I,I) * DTEMDX(I)
!  (No off-diagonal coupling for transverse isotropic material
!   aligned with the principal axes)
!----------------------------------------------------------------------
      FLUX(1) = -K11 * DTEMDX(1)
      FLUX(2) = -K22 * DTEMDX(2)
      IF (NTENS .GE. 3) FLUX(3) = -K33 * DTEMDX(3)

!----------------------------------------------------------------------
!  Conductivity Jacobian: DFLUXDG(I,J) = d FLUX(I) / d DTEMDX(J)
!  For diagonal conductivity: DFLUXDG(I,I) = -K(I,I), off-diag = 0
!----------------------------------------------------------------------
      DO I = 1, NTENS
         DO J = 1, NTENS
            DFLUXDG(I,J) = ZERO
         END DO
      END DO
      DFLUXDG(1,1) = -K11
      DFLUXDG(2,2) = -K22
      IF (NTENS .GE. 3) DFLUXDG(3,3) = -K33

!----------------------------------------------------------------------
!  Thermal Jacobian: dFLUX/dT = 0
!  (conductivity does not depend on T in this version;
!   if adding k(T) later, update DFLUDT accordingly)
!----------------------------------------------------------------------
      DO I = 1, NTENS
         DFLUDT(I) = ZERO
      END DO

      RETURN
      END


!============================================================================
!  CALC_DAMAGE_RATE — three-phase CDM damage rate subroutine
!  Called from UMAT (4 times per increment, once per damage variable)
!
!  Rate law:
!    dD/dN = [ LAMBDA_eff * Y * exp(-DELTA_eff * N)   Phase I  (burst)
!             + GAMMA_eff  * max(Y - Y_TH, 0)^Q      ] Phase II (Paris)
!           / max(1 - D/D_F, R_FLOOR)^K_TERT           Phase III (runaway)
!
!  Effective parameters (Y-power scaling):
!    GAMMA_eff  = GAMMA  * Y^P_GAMMA
!    LAMBDA_eff = LAMBDA * Y^P_LAMBDA
!    DELTA_eff  = DELTA  * Y^P_DELTA
!============================================================================
      SUBROUTINE CALC_DAMAGE_RATE(Y, D, Q, GAMMA, LAMBDA, DELTA,
     1                            Y_TH, K_TERT, D_F,
     2                            P_GAMMA, P_LAMBDA, P_DELTA,
     3                            N_CYCLE, TOLER, RATE)
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: Y, D, Q, GAMMA, LAMBDA, DELTA
      REAL*8, INTENT(IN)  :: Y_TH, K_TERT, D_F
      REAL*8, INTENT(IN)  :: P_GAMMA, P_LAMBDA, P_DELTA
      REAL*8, INTENT(IN)  :: N_CYCLE, TOLER
      REAL*8, INTENT(OUT) :: RATE
      REAL*8 :: NUM, DENOM_R, R, GAMMA_E, LAMBDA_E, DELTA_E
      REAL*8 :: TERM_PHASE1, TERM_PHASE2, Y_X
      REAL*8, PARAMETER :: R_FLOOR = 1.0D-2

C --- Gate conditions
      IF (Y .LE. TOLER) THEN; RATE = 0.0D0; RETURN; END IF
      IF (D .GE. D_F)   THEN; RATE = 0.0D0; RETURN; END IF

C --- Y-power scaling of effective parameters
      IF (ABS(P_GAMMA)  .GT. TOLER) THEN
         GAMMA_E  = GAMMA  * Y**P_GAMMA
      ELSE
         GAMMA_E  = GAMMA
      END IF
      IF (ABS(P_LAMBDA) .GT. TOLER) THEN
         LAMBDA_E = LAMBDA * Y**P_LAMBDA
      ELSE
         LAMBDA_E = LAMBDA
      END IF
      IF (ABS(P_DELTA)  .GT. TOLER) THEN
         DELTA_E  = DELTA  * Y**P_DELTA
      ELSE
         DELTA_E  = DELTA
      END IF

C --- Phase II Paris term (only active above endurance threshold)
      Y_X = Y - Y_TH
      IF (Y_X .GT. TOLER .AND. Q .GT. TOLER .AND.
     1    GAMMA_E .GT. TOLER) THEN
         TERM_PHASE2 = GAMMA_E * Y_X**Q
      ELSE
         TERM_PHASE2 = 0.0D0
      END IF

C --- Phase I exponential decay term
      TERM_PHASE1 = LAMBDA_E * Y * EXP(-DELTA_E * N_CYCLE)

C --- Total numerator
      NUM = TERM_PHASE2 + TERM_PHASE1

C --- Phase III denominator (runaway as D approaches D_F)
      IF (K_TERT .GT. TOLER) THEN
         R       = MAX(1.0D0 - D / D_F, R_FLOOR)
         DENOM_R = R**K_TERT
         RATE    = NUM / DENOM_R
      ELSE
         RATE = NUM
      END IF

      RETURN
      END