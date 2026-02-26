!!!!!!!!!!! UMAT suboutine for static mechanical analysis !!!!!!!!!
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

      DIMENSION FLOW(6), STRSS(6), STRANT(6), DSTRESS(6), CDTHREE(6,6)
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      DATA NEWTON,TOLER/10,1.D-6/
	  
c *** get Properties, and Damage model parameters (matrix)
      E11_0 = PROPS(1)
      E22_0 = PROPS(2)
      E33_0 = PROPS(3)
      NU12  = PROPS(4)
      NU23  = PROPS(5)	
	  G12_0 = PROPS(6)
	  G13_0 = PROPS(7)
	  G23_0 = PROPS(8)
	  
c ***    For transversely isotropic material: ! NU13 = NU12, G13 = G12
      NU13 = NU12
	  	  
c ***    Strain calculation (Step 1)
      STRANT(1) = STRAN(1) + DSTRAN(1)
      STRANT(2) = STRAN(2) + DSTRAN(2)
      STRANT(3) = STRAN(3) + DSTRAN(3)
      STRANT(4) = STRAN(4) + DSTRAN(4)
      STRANT(5) = STRAN(5) + DSTRAN(5)
      STRANT(6) = STRAN(6) + DSTRAN(6)

	  
c ***    Calculate damaged elastic moduli
      E11 = E11_0* (ONE - STATEV(9))
      E22 = E22_0* (ONE - STATEV(9))
      E33 = E33_0* (ONE - STATEV(9))
      G12 = G12_0* (ONE - STATEV(9))
      G13 = G13_0* (ONE - STATEV(9))
      G23 = G23_0 * (ONE - STATEV(9))

	  
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
	  DDSDDE(4,4) = G23
	  DDSDDE(5,5) = G13
	  DDSDDE(6,6) = G12
	  
	  TRVAL = DSTRAN(1)+DSTRAN(2)+DSTRAN(3)
	  
	  
	  
c ***   UPDATE THE STRESS (Step 3)
       DO K1=1,NTENS
         STRESS(K1)=0.0D0
         DO K2=1,NTENS
            STRESS(K1)=STRESS(K1)+DDSDDE(K2,K1)*STRANT(K2)
        ENDDO
      ENDDO	  


c ***  Total dissipated energy and temperature rise

        SMISES=(STRESS(1)-STRESS(2))*(STRESS(1)-STRESS(2)) +
     1          (STRESS(2)-STRESS(3))*(STRESS(2)-STRESS(3)) +
     1          (STRESS(3)-STRESS(1))*(STRESS(3)-STRESS(1))
        DO 90 K1=NDI+1,NTENS
              SMISES=SMISES+SIX*STRESS(K1)*STRESS(K1)
 90     CONTINUE
        SVM=SQRT(SMISES/TWO)
	  
	  I1 = STRESS(1) + STRESS(2) + STRESS(3)
	  Tmi = 40.96
	  Cmi = 64
	  RATIO = (Cmi/Tmi)
	  TERM1 = (I1*(RATIO-1)) 
	  TERM2 = (RATIO-1)*2 *I1*I1
	  TERM3 = 4*(RATIO*SVM*SVM)
	  TERM4 = TERM2+TERM3
	  TERM5 = SQRT(TERM4)
	  TERM6 = TERM1+TERM5
	  TERM7 = TERM6/(RATIO*2)
   
	  EQSTRESS = SVM 
	  STATEV(3) = EQSTRESS 
	  STATEV(4) = Tmi + STATEV(5)        
	
	  f = EQSTRESS - STATEV(4)	
	  
      IF(f.GE.0) THEN  
        DAMAGE = (STATEV(4))/Tmi
        RM=EQSTRESS/Tmi       
        dKm = EQSTRESS - Tmi
        STATEV(5) = dKm
        STATEV(1) = DAMAGE
        STATEV(2) = RM        
        Dm=1-exp(0.33*(1-((STATEV(3))/Tmi)))   
        STATEV(9) = Dm         
      END IF 	  
	  
      RETURN
      END