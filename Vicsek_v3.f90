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
 INTEGER:: IFRAME

 REAL(S):: BOXX,BOXY,DT,ETA,RINT,VC
 REAL(S):: RX(NPARTX),RY(NPARTX),THETA(NPARTX)
 REAL(S):: VX(NPARTX),VY(NPARTX),DRIJ(NPARTX,NPARTX)
 REAL(S):: MNANG(NPARTX)

END MODULE AMSVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM ACTIVEMOTIONSIM
 USE AMSVAR
 IMPLICIT NONE

 CALL BEGINAMS                              !SET UP SIMULATION PARAMETERS
 
 DO IAMS=1,AMSTEPS                          !SIMULATION
    CALL IJDIST                             !SEPARATION DISTANCE BETWEEN I-J
    CALL MEANAG                             !MEAN ANGLE BETWEEN INTERACTION PARTICLES
    CALL AMMOVE                             !MOVE PARTICLES
    CALL AMFILM                             !SAVE A SNAPSHOT TO BUILD A MOVIE
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

 OPEN(UNIT=12,FILE='AMFilm.xyz')            !FILE WITH SNAPSHOTS FOR MOVIE

 RETURN
END SUBROUTINE BEGINAMS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE READINP
 USE AMSVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='AMS.inp',STATUS='OLD')

 READ(10,*)NPART                            !NUMBER OF PARTICLES 
 READ(10,*)AMSTEPS                          !NUMBER SIMULATION STEPS
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
 REAL(S):: RN

 BOXY=BOXX                                  !SQUARE SIMULATION BOX
 IFRAME=0                                   !NUMBER OF SNAPSHOT

 DO I=1,NPART
    CALL RANDOM_NUMBER(RN)                  !RANDOM NUMBER FROM PLANAR DISTRIBUTION
    RX(I)=RN*BOXX                           !X RANDOM POSITION IN BOX
    CALL RANDOM_NUMBER(RN)                  !RANDOM NUMBER FROM PLANAR DISTRIBUTION
    RY(I)=RN*BOXY                           !Y RANDOM POSITION IN BOX

    CALL RANDOM_NUMBER(RN)                  !RANDOM NUMBER FROM PLANAR DISTRIBUTION
    THETA(I)=(RN-0.5)*2.0*PI                !RANDOM VALUE IN [-PI,PI]

    CALL RANDOM_NUMBER(RN)                  !RANDOM NUMBER FROM PLANAR DISTRIBUTION
    VX(I)=VC*COS(THETA(I))                  !X VELOCITY ACCORDING TO RANDOM ORIENTATION
    CALL RANDOM_NUMBER(RN)                  !RANDOM NUMBER FROM PLANAR DISTRIBUTION
    VY(I)=VC*SIN(THETA(I))                  !Y VELOCITY ACCORDING TO RANDOM ORIENTATION
 ENDDO

 OPEN(UNIT=11,FILE='AMSnap.xyz')            !INITIAL CONFIGURATION SNAPSHOT 

 WRITE(11,*)NPART
 WRITE(11,*)'FRAME',1

 DO I=1,NPART
    WRITE(11,*)'C',RX(I),RY(I),0.0
 ENDDO 

 CLOSE(UNIT=11) 

 RETURN
END SUBROUTINE INITCNF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE IJDIST
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I,J
 REAL(S):: X,Y,DX,DY,RIJ

 DO I=1,NPART
    X=RX(I)
    Y=RY(I)
    DO J=1,NPART
       DX=ABS(X-RX(J))
       DY=ABS(Y-RY(J))

       IF(DX .GT. (BOXX-DX))DX=BOXX-DX
       IF(DY .GT. (BOXY-DY))DY=BOXY-DY

       RIJ=SQRT(DX*DX + DY*DY)
       DRIJ(I,J)=RIJ
    ENDDO
 ENDDO
 RETURN
END SUBROUTINE IJDIST
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE MEANAG
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I,J,NCOUNT 
 REAL(S):: LANG(NPART),MNSIN,MNCOS
!NOTE: THIS ROUTINE CAN BE DONE MORE SIMPLE, CHECKOUT IJDIST, COULD BE UPDATE 
!      DIRECTLY IN THAT STEP...

 DO I=1,NPART
    NCOUNT=0
    DO J=1,NPART
       IF(DRIJ(I,J) .LE. RINT)THEN          !IF PARTICLE I INTERACTS WITH J
         NCOUNT=NCOUNT + 1                  !ORIENTATION MUST CHANGE 
         LANG(NCOUNT)=THETA(J)              !ACCORDING TO AN AVERAGE
       ENDIF
    ENDDO

    MNSIN=0.0                               !SUM TO COMPUTE AVERAGE OF SIN
    MNCOS=0.0                               !SUM TO COMPUTE AVERAGE OF COS
    DO J=1,NCOUNT
       MNCOS=MNCOS + COS(LANG(J))
       MNSIN=MNSIN + SIN(LANG(J))
    ENDDO
    MNCOS=MNCOS/REAL(NCOUNT)                !AVERAGE COS
    MNSIN=MNSIN/REAL(NCOUNT)                !AVERAGE SIN

    MNANG(I)=ATAN2(MNSIN,MNCOS)             !NEW ORIENTATION FOR I PARTICLE
 ENDDO

 RETURN
END SUBROUTINE MEANAG
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE AMMOVE
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(S):: RN,RXN,RYN

 DO I=1,NPART
    CALL RANDOM_NUMBER(RN)                  !RANDOM NUMBER FROM PLANAR DIST
    RN=(RN-0.5)*PI                      !RANDOM VALUE IN [-PI,PI]
    THETA(I)=MNANG(I) + ETA*RN              !NEW ANGLE FOR I PARTICLE

    VX(I)=VC*COS(THETA(I))                  !NEW X VELOCITY
    VY(I)=VC*SIN(THETA(I))                  !NEW Y VELOCITY

    RXN=RX(I) + VX(I)*DT                    !X POSITION AT t + ðt
    RYN=RY(I) + VY(I)*DT                    !Y POSITION AT t + ðt

    RX(I)=MOD((RXN+BOXX),BOXX)              !X PERIODIC BOUNDARY CONDITION
    RY(I)=MOD((RYN+BOXY),BOXY)              !Y PERIODIC BOUNDARY CONDITION
 ENDDO

 RETURN
END SUBROUTINE AMMOVE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE AMFILM
 USE AMSVAR
 IMPLICIT NONE
 INTEGER:: I

 IFRAME=IFRAME + 1

 WRITE(12,*)NPART
 WRITE(12,*)'FRAME',IFRAME
 DO I=1,NPART
    WRITE(12,*)'C',RX(I),RY(I),0.0
 ENDDO

 RETURN
END SUBROUTINE AMFILM
