!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%**
!%**
!%**
!%**
!%**
!%**
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE VICVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: S      = KIND(1.0)    !VARIABLES PRECISION
 INTEGER, PARAMETER:: NPARTX = 2500         !MAXIMUM NUMBER OF PARTICLES
 INTEGER, PARAMETER:: NCELX  = 100           !MAXIMUM NUMBER OF CELLS IN EACH DIRECTION
 INTEGER, PARAMETER:: TCELX  = NCELX**2 + 4*(NCELX+1)
 
 REAL(S), PARAMETER:: SIGMA  = 1.0          !MOLECULAR DIAMETER
 REAL(S), PARAMETER:: PI     = ACOS(-1.0)   !PI NUMBER
 REAL(S), PARAMETER:: V      = 0.1          !CONSTANT VELOCITY

 INTEGER:: NPART,NCEL,TCEL,IAM,AMSTEPS
 INTEGER:: NEIGH(NCELX**2,8)                !NUMBER OF CELLS AND NEIGHBOURS OF EACH ONE
 INTEGER:: PCELL(NPARTX),CELL(TCELX,NPARTX)
 INTEGER:: INFCL(TCELX),IFRAME

 REAL(S):: BOXX,BOXY,AMP
 REAL(S):: RX(NPARTX),RY(NPARTX)
 REAL(S):: ANG(NPARTX),NANG(NPARTX)

END MODULE VICVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM VICSEK
 USE VICVAR
 IMPLICIT NONE
 
 CALL BEGINSIM

 DO IAM=1,AMSTEPS
    CALL UPDCELLS                           !IDENTIFY PARTICLE IN CELLS 
    CALL UPDGHOST                           !UPDATE GHOST CELLS
    CALL INTERACT                           !UPDATE ALIGMENT
    CALL MOVEPART                           !MOVE PARTICLES
    CALL ACMMOVIE                           !ACTIVE MOTION MOVIE
 ENDDO

 STOP
END PROGRAM VICSEK
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE BEGINSIM
 IMPLICIT NONE
 
 CALL READINP                               !READ SIMULATION PARAMETERS
 CALL ICONFIG                               !BUILD RANDOM INITIAL CONFIGURATION
 CALL SETPARM                               !SET SIMULATION PARAMETERS
 CALL NEIGHBR                               !SET A MAP OF NEIGHBOURS OF EACH CELL

 OPEN(UNIT=12,FILE='AMMovie.xyz')

 RETURN
END SUBROUTINE BEGINSIM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE READINP
 USE VICVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='AM.inp')
 REWIND(10)

 READ(10,*)NPART                            !NUMBER OF PARTICLES
 READ(10,*)BOXX                             !BOX LENGHT
 READ(10,*)AMSTEPS                          !SIMULATION TIME STEPS

 CLOSE(UNIT=10)

 RETURN
END SUBROUTINE READINP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE SETPARM
 USE VICVAR
 IMPLICIT NONE

 NCEL=INT(BOXX/SIGMA)                       !NUMBER OF CELL IN X DIRECTION
 TCEL=NCEL*NCEL + 4*(NCEL + 1)              !TOTAL NUMBER OF CELLS
 AMP=0.5                                    !NOISE AMPLITUDE
 IFRAME=0

 RETURN
END SUBROUTINE SETPARM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ICONFIG
 USE VICVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(S):: RN

 BOXY=BOXX                                  !SQUARE SIMULATION BOX

 DO I=1,NPART
    CALL RANDOM_NUMBER(RN) 
    RX(I)=RN*BOXX
    CALL RANDOM_NUMBER(RN) 
    RY(I)=RN*BOXY
    CALL RANDOM_NUMBER(RN)
    ANG(I)=2.0*PI*RN - PI                   !RANDOM ANGLE IN[-PI,PI]
 ENDDO


 OPEN(UNIT=11,FILE="IniCnf.xyz")
 WRITE(11,*)NPART
 WRITE(11,*)'FRAME',1

 DO I=1,NPART
    WRITE(11,*)'C',RX(I),RY(I),0.0
 ENDDO

 CLOSE(UNIT=11)

 RETURN
END SUBROUTINE ICONFIG
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE NEIGHBR
 USE VICVAR 
 IMPLICIT NONE
 INTEGER:: I,X,Y

 !CELLS ARE LABELLED FROM 0,0 TO NCEL,NCEL
 DO I=1,NCEL*NCEL
    X=MOD((I-1),NCEL)                       !INIT AT X DIRECTION
    Y=INT((I-1)/NCEL)                       !MOVE UPWARDS IN Y DIRECTION

    NEIGH(I,1)=X + 1 + Y*NCEL + 1           !RIGHT NEIGHBOUR 
    NEIGH(I,2)=X - 1 + Y*NCEL + 1           !LEFT NEIGHBOUR 
    NEIGH(I,3)=X + 1 + (Y+1)*NCEL + 1       !DIAGONAL RIGHT+UP NEIGHBOUR 
    NEIGH(I,4)=X - 1 + (Y-1)*NCEL + 1       !DIAGONAL LEFT+DOWN NEIGHBOUR 
    NEIGH(I,5)=X + 1 + (Y-1)*NCEL + 1       !DIAGONAL RIGHT+DOWN NEIGHBOUR 
    NEIGH(I,6)=X - 1 + (Y+1)*NCEL + 1       !DIAGONAL LEFT+UP NEIGHBOUR 
    NEIGH(I,7)=X + (Y-1)*NCEL + 1           !DOWN NEIGHBOUR
    NEIGH(I,8)=X + (Y+1)*NCEL + 1           !UP NEIGHBOUR

    IF(Y .EQ. (NCEL-1))THEN                 !GHOST CELLS ON ARISTA FOR PBC
      NEIGH(I,8)=X+1 + NCEL*NCEL
      NEIGH(I,3)=X+2 + NCEL*NCEL
      NEIGH(I,6)=X + NCEL*NCEL
    ENDIF

    IF(Y .EQ. 0)THEN
      NEIGH(I,7)=X+1 + NCEL + NCEL*NCEL
      NEIGH(I,5)=X+2 + NCEL + NCEL*NCEL
      NEIGH(I,4)=X + NCEL + NCEL*NCEL
    ENDIF

    IF(X .EQ. (NCEL-1))THEN
      NEIGH(I,1)=Y+1 + 2*NCEL + NCEL*NCEL
      NEIGH(I,3)=Y+2 + 2*NCEL + NCEL*NCEL
      NEIGH(I,5)=Y + 2*NCEL + NCEL*NCEL
    ENDIF

    IF(X .EQ. 0)THEN
      NEIGH(I,2)=Y+1 + 3*NCEL + NCEL*NCEL
      NEIGH(I,6)=Y+2 + 3*NCEL + NCEL*NCEL
      NEIGH(I,4)=Y + 3*NCEL + NCEL*NCEL
    ENDIF
 ENDDO

 !GHOST CELLS ON VERTEX FOR PBC
 NEIGH(NCEL**2,3)=1 + 4*NCEL + NCEL*NCEL
 NEIGH(NCEL,5)=2 + 4*NCEL + NCEL*NCEL
 NEIGH(1,4)=3 + 4*NCEL + NCEL*NCEL
 NEIGH(NCEL*(NCEL-1)+1,6)=4 + 4*NCEL + NCEL*NCEL

 RETURN
END SUBROUTINE NEIGHBR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE UPDCELLS
 USE VICVAR
 IMPLICIT NONE
 INTEGER:: I,XCELL,YCELL,CINDX
 
 DO I=1,TCEL
    INFCL(I)=0
 ENDDO

 DO I=1,NPART
    IF(RX(I) .GE. BOXX)THEN                 !BOUNDARY CONDITIONS (HARD WAY!?)
      RX(I)=RX(I) - BOXX
    ELSEIF(RX(I) .LT. 0.0)THEN
      RX(I)=RX(I) + BOXX
    ENDIF
    IF(RY(I) .GE. BOXY)THEN
      RY(I)=RY(I) - BOXY
    ELSEIF(RY(I) .LT. 0.0)THEN
      RY(I)=RY(I) + BOXY
    ENDIF

    XCELL=INT(RX(I)/SIGMA)
    YCELL=INT(RY(I)/SIGMA)
    CINDX=XCELL + YCELL*NCEL + 1            !CELL INDEX

    PCELL(I)=CINDX
    INFCL(CINDX)=INFCL(CINDX) + 1
    CELL(CINDX,INFCL(CINDX))=I
 ENDDO

 RETURN
END SUBROUTINE UPDCELLS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE UPDGHOST
 USE VICVAR
 IMPLICIT NONE
 INTEGER:: I,J,BOX,UP,DOWN
 INTEGER:: AUX,IPART
 !WARNING: POSSIBLE ERROR HERE
!*******************
! DO I=1,NPART
!    PCELL(I)=NCEL**2
!    UP(I)=PCELL
!    DOWN(I)=NCEL*(NCEL-1)
! ENDDO

 BOX=NCEL*NCEL
 UP=NCEL**2
 DOWN=NCEL*(NCEL-1)
!******************
 DO I=1,NCEL
    !GHOST CELLS ABOVE Y=L
    AUX=INFCL(I)
    INFCL(UP+I)=AUX
    DO J=1,AUX
       IPART=CELL(I,J)
       CELL(UP+I,J)=IPART + NPART
       RX(IPART+NPART)=RX(IPART)
       RY(IPART+NPART)=RX(IPART) + BOXY
       ANG(IPART+NPART)=ANG(IPART)
    ENDDO
    !GHOST CELLS BELOW Y=0
    AUX=INFCL(DOWN+I)
    INFCL(BOX+NCEL+I)=AUX
    DO J=1,AUX
       IPART=CELL(DOWN+1,J)
       CELL(BOX+NCEL+I,J)=IPART + NPART
       RX(IPART+NPART)=RX(IPART)
       RY(IPART+NPART)=RY(IPART) - BOXY
       ANG(IPART+NPART)=ANG(IPART)
    ENDDO
    !GHOST CELLS AT RIGHT X=NCEL-1
    AUX=INFCL((I-1)*NCEL+1)
    INFCL(BOX+2*NCEL+I)=AUX
    DO J=1,AUX
       IPART=CELL((I-1)*NCEL+1,J)
       CELL(BOX+2*NCEL+I,J)=IPART + NPART
       RX(IPART+NPART)=RX(IPART) + BOXX
       RY(IPART+NPART)=RY(IPART)
       ANG(IPART+NPART)=ANG(IPART)
    ENDDO
    !GHOST CELLS AT LEFT OF X=0
    AUX=INFCL(I*NCEL)
    INFCL(BOX+3*NCEL+1)=AUX
    DO J=1,AUX
       IPART=CELL(I*NCEL,J)
       CELL(BOX+3*NCEL+I,J)=IPART + NPART
       RX(IPART+NPART)=RX(IPART) - BOXX
       RY(IPART+NPART)=RY(IPART)
       ANG(IPART+NPART)=ANG(IPART)
    ENDDO
    !GHOST CELL IN X=NCEL Y=NCEL
    AUX=INFCL(1)
    INFCL(NCEL**2+4*NCEL+1)=AUX
    DO J=1,AUX
       IPART=CELL(1,J)
       CELL(NCEL**2+4*NCEL+1,J)=IPART + NPART
       RX(IPART+NPART)=RX(IPART) + BOXX
       RY(IPART+NPART)=RY(IPART) + BOXY
       ANG(IPART+NPART)=ANG(IPART)
    ENDDO
    !GHOST CELL IN X=NCELL Y-1
    AUX=INFCL(BOX-NCEL+1)
    INFCL(BOX+4*NCEL+2)=AUX
    DO J=1,AUX
       IPART=CELL(BOX-NCEL+1,J)
       CELL(BOX+4*NCEL+2,J)=IPART + NPART
       RX(IPART+NPART)=RX(IPART) + BOXX
       RY(IPART+NPART)=RY(IPART) - BOXY
       ANG(IPART+NPART)=ANG(IPART)
    ENDDO
    !GHOST CELL IN X=-1 Y=-1 
    AUX=INFCL(NCEL**2)
    INFCL(NCEL**2+4*NCEL+3)=AUX
    DO J=1,AUX
       IPART=CELL(NCEL**2,J)
       CELL(BOX+4*NCEL+3,J)=IPART + NPART
       RX(IPART+NPART)=RX(IPART) - BOXX
       RY(IPART+NPART)=RY(IPART) - BOXY
       ANG(IPART+NPART)=ANG(IPART)
    ENDDO
    !GHOST CELL IN X=-1 Y=NCEL
    AUX=INFCL(NCEL)
    INFCL(TCEL)=AUX
    DO J=1,AUX
       IPART=CELL(NCEL,J)
       CELL(TCEL,J)=IPART + NPART
       RX(IPART+NPART)=RX(IPART) - BOXX
       RY(IPART+NPART)=RY(IPART) + BOXY
       ANG(IPART+NPART)=ANG(IPART)
    ENDDO

 ENDDO

 RETURN
END SUBROUTINE UPDGHOST
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INTERACT
 USE VICVAR
 IMPLICIT NONE
 INTEGER:: I,J,K,IPART
 INTEGER:: COUNT,AUX
 REAL(S):: SM,CM,DX,DY,RIJ

 DO I=1,NPART
    COUNT=1
    SM=SIN(ANG(I))
    CM=COS(ANG(I))

    DO J=1,8
       AUX=NEIGH(PCELL(I),J)
       DO K=1,INFCL(AUX)
          IPART=CELL(AUX,K)
          DX=RX(I) - RX(IPART)
          DY=RY(I) - RY(IPART)
          RIJ=SQRT(DX*DX + DY*DY)
          !PARTICLES INTERACTS IF THEY ARE CLOSE
          IF(RIJ .LE. SIGMA)THEN
            COUNT=COUNT+1
            SM=SM + SIN(ANG(IPART))
            CM=CM + COS(ANG(IPART))
          ENDIF
       ENDDO
    ENDDO

    SM=SM/REAL(COUNT)
    CM=CM/REAL(COUNT)
    NANG(I)=ATAN2(SM,CM)                    !COMPUTES ANGLE
 ENDDO

 RETURN
END SUBROUTINE INTERACT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE MOVEPART
 USE VICVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(S):: RN,ANGNOISE

 DO I=1,NPART
    CALL RANDOM_NUMBER(RN)
    ANGNOISE=-0.5*AMP + RN*AMP
    ANG(I)=NANG(I) + ANGNOISE
    RX(I)=RX(I) + V*COS(ANG(I))
    RY(I)=RY(I) + V*SIN(ANG(I))
 ENDDO
 RETURN
END SUBROUTINE MOVEPART
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ACMMOVIE
 USE VICVAR
 IMPLICIT NONE
 INTEGER:: I 

 IFRAME=IFRAME + 1
 WRITE(12,*)NPART
 WRITE(12,*)'FRAM',IFRAME

 DO I=1,NPART
    WRITE(12,*)'C',RX(I),RY(I),0.0
 ENDDO

 RETURN
END SUBROUTINE ACMMOVIE
