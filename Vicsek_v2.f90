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
 INTEGER, PARAMETER:: AMSTEPS = 200         !NUMBER OF TIME STEPS

 REAL(S), PARAMETER:: V0      = 0.5         !VELOCITY
 REAL(S), PARAMETER:: ETA     = 0.5         !ANGLE FLUCTUATION MAGNITUDE (RADIANS)
 REAL(S), PARAMETER:: BOXX    = 10.0        !SIMULATION BOX X LENGHT
 REAL(S), PARAMETER:: RINT    = 1.0         !INTERACTION RANGE BETWEEN PARTICLES
 REAL(S), PARAMETER:: DT      = 0.2         !TIME STEP
 REAL(S), PARAMETER:: PI      = ACOS(-1.0)  !PI NUMBER
 
 INTEGER:: IAM,IFRAME

 REAL(S):: RX(NPARTX),RY(NPARTX),THETA(NPARTX)
 REAL(S):: VX(NPARTX),VY(NPARTX),MTHETA(NPARTX)
 REAL(S):: BOXY
END MODULE AMSVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM ACTIVEMOTION
 USE AMSVAR
 IMPLICIT NONE

 CALL BEGINAMS                              !BEGIN SIMULATION

 DO IAM=1,AMSTEPS
    CALL AMSMOTION                          !PARTICLE MOTION
    CALL MEANANGLE                          !MEAN ANGLE OF NEIGHBOURS WITHIN RINT
    CALL NEWVELOCI                          !UPDATE VELOCITIES 
    IF(MOD(IAM,2) .EQ. 0)CALL AMSMOVIE      !TAKE A FRAME FOR A MOVIE
 ENDDO

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
 IFRAME=0

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

 OPEN(UNIT=11,FILE='AMSFilm.xyz')

 RETURN
END SUBROUTINE INITCNF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE AMSMOTION
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(S):: RXN,RYN

 !LOOK IF BOUNDARY CONDITIONS CAN BE IMPROVED ALONG CENTERING BOX SIMULATION
 DO I=1,NPART
    RXN=RX(I) + VX(I)*DT                    !X POSITION @ t+ðt
    RYN=RY(I) + VY(I)*DT                    !Y POSITION @ t+ðt
    RXN=MOD(RXN,BOXX)                       !X PERIODIC BOUNDARY CONDITION
    RYN=MOD(RYN,BOXY)                       !Y PERIODIC BOUNDARY CONDITION
    RX(I)=RXN                               !UPDATE X POSITION @ t + ðt
    RY(I)=RYN                               !UPDATE Y POSITION @ t + ðt
 ENDDO

 RETURN
END SUBROUTINE AMSMOTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE MEANANGLE
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I,J
 REAL(S):: X,Y,DX,DY,RIJ
 REAL(S):: SX,SY,RN

 SX=0.0
 SY=0.0

 DO I=1,NPART
    X=RX(I)
    Y=RY(I)
    DO J=1,NPART
       DX=X-RX(J)
       DY=Y-RY(J)
       RIJ=SQRT(DX*DX + DY*DY)               !**MINIMUM IMAGE CONVENTION!?
       IF(RIJ .LT. RINT)THEN
         SX=SX + COS(THETA(J))
         SY=SY + SIN(THETA(J))
         MTHETA(J)=ATAN2(SY,SX)
       ENDIF
    ENDDO
 ENDDO

 DO I=1,NPART
    CALL RANDOM_NUMBER(RN)
    THETA(I)=MTHETA(I) + ETA*(RN - 0.5)
 ENDDO

 RETURN
END SUBROUTINE MEANANGLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE NEWVELOCI
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I
 
 DO I=1,NPART
    VX(I)=V0*COS(THETA(I))
    VY(I)=V0*SIN(THETA(I))
 ENDDO

 RETURN
END SUBROUTINE NEWVELOCI
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE AMSMOVIE
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I

 IFRAME=IFRAME + 1
 WRITE(11,*)NPART
 WRITE(11,*)'FRAME',IFRAME

 DO I=1,NPART
    WRITE(11,*)'C',RX(I),RY(I),0.0D0
 ENDDO

 RETURN
END SUBROUTINE AMSMOVIE
