!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%***************************************************************************%%
!%**  PROGRAM         ACTIVE MOTION SIMULATION                             **%%
!%**  AUTHOR          ALPIXELS                                             **%%
!%**  DATE            JANUARY 2022                                         **%%
!%**  COPYRIGHT       GNU/LGPL-V3                                          **%%
!%**                                                                       **%% 
!%**  ALGHORITHM      VICSEK MODEL                                         **%%
!%**  OBS             THIRD VERSION                                        **%%
!%***************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE AMSVAR
 INTEGER, PARAMETER:: S      = KIND(1.0)    !VARIABLES PRECISION
 INTEGER, PARAMETER:: NPARTX = 2500         !MAXIMUM NUMBER OF PARTICLES

 REAL(S), PARAMETER:: PI     = ACOS(-1.0)   !PI NUMBER
 
 INTEGER:: NPART,IAMS,AMSTEPS

 REAL(S):: BOXX,BOXY,DT,ETA,RINT,VC
 REAL(S):: RX(NPARTX),RY(NPARTX),THETA(NPARTX)
 REAL(S):: VX(NPARTX),VY(NPARTX)
END MODULE AMSVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM ACTIVEMOTIONSIM
 USE AMSVAR
 IMPLICIT NONE

 CALL BEGINAMS
 
 DO IAMS=1,AMSTEPS
 ENDDO

 STOP
END PROGRAM ACTIVEMOTIONSIM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SUBROUTINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE BEGINAMS
 USE AMSVAR
 IMPLICIT NONE

 CALL READINP                               !READ SIMULATION PARAMETERS
 CALL INITCNF                               !BUILT INITIAL CONFIGURATION
 RETURN
END SUBROUTINE BEGINAMS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE READINP
 USE AMSVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='AMS.inp',STATUS='OLD')

 READ(10,*)NPART                            !NUMBER OF PARTICLES 
 READ(10,*)BOXX                             !BOX SIMULATION SIZE
 READ(10,*)VC                               !CONSTANT VELOCITY
 READ(10,*)DT                               !TIME STEP
 READ(10,*)ETA                              !ANGLE MAXIMUM VARIATION
 READ(10,*)RINT                             !INTERACTION RANGE

 CLOSE(UNIT=10)

 RETURN
END SUBROUTINE READINP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INITCNF
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I

 DO I=1,NPART
 ENDDO

 RETURN
END SUBROUTINE INITCNF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
