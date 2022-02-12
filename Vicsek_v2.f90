!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%***************************************************************************%%
!%**  PROGRAM         ACTIVE MOTION SIMULATION                             **%%
!%**  AUTHOR          ALPIXELS                                             **%%
!%**  DATE            JANUARY 2022                                         **%%
!%**  COPYRIGHT       GNU/LGPL-V3                                          **%%
!%**                                                                       **%% 
!%**  ALGHORITHM      VICSEK MODEL                                         **%%
!%**  OBS             FIRST VERSION                                        **%%
!%***************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE AMSVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: S       = KIND(1.0)   !VARIABLES PRECISION
 INTEGER, PARAMETER:: NPARTX  = 2000        !MAXIMUM NUMBER OF PARTICLES

 INTEGER, PARAMETER:: NPART   = 500         !NUMBER OF PARTICLES

 REAL(S), PARAMETER:: V0      = 0.5         !VELOCITY
 REAL(S), PARAMETER:: ETA     = 0.5         !ANGLE FLUCTUATION MAGNITUDE (RADIANS)
 REAL(S), PARAMETER:: BOXX    = 10.0        !SIMULATION BOX X LENGHT
 REAL(S), PARAMETER:: RINT    = 1.0         !INTERACTION RANGE BETWEEN PARTICLES
 REAL(S), PARAMETER:: DT      = 0.2         !TIME STEP
 REAL(S), PARAMETER:: PI      = ACOS(-1.0)  !PI NUMBER
 
 INTEGER:: IAM,AMSTEPS

 REAL(S):: RX(NPARTX),RY(NPARTX),THETA(NPARTX)
 REAL(S):: VX(NPARTX),VY(NPARTX)
 REAL(S):: BOXY
END MODULE AMSVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM ACTIVEMOTION
 USE AMSVAR
 IMPLICIT NONE

 CALL BEGINAMS                              !BEGIN SIMULATION

 STOP
END PROGRAM ACTIVEMOTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SUBROUTINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE BEGINAMS
 USE AMSVAR
 IMPLICIT NONE

 CALL INITCNF                               !BUILT INITIAL CONFIGURATION

 RETURN 
END SUBROUTINE BEGINAMS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INITCNF
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(S):: RN

 BOXY=BOXX

 DO I=1,NPART                               !SQUARE SIMULATION BOX
    CALL RANDOM_NUMBER(RN)
    RX(I)=RN*BOXX                           !X RANDOM POSITION 
    CALL RANDOM_NUMBER(RN)
    RY(I)=RN*BOXY                           !Y RANDOM POSITION

    CALL RANDOM_NUMBER(RN)
    THETA(I)=2.0*PI*RN                      !THETA IN [0,2PI]
    VX(I)=V0*COS(THETA(I))                  !VELOCITY IN X DIRECTION
    VY(I)=V0*SIN(THETA(I))                  !VELOCITY IN Y DIRECTION
 ENDDO

 OPEN(UNIT=10,FILE='AMSnap.xyz')

 WRITE(10,*)NPART
 WRITE(10,*)'FRAME',1
 DO I=1,NPART
    WRITE(10,*)'C',RX(I),RY(I),0.0
 ENDDO

 CLOSE(UNIT=10)
 RETURN
END SUBROUTINE INITCNF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
