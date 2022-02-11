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
 INTEGER, PARAMETER:: NPARTX = 2000         !MAXIMUM NUMBER OF PARTICLES
 INTEGER, PARAMETER:: NCELX  = 50           !MAXIMUM NUMBER OF CELLS IN EACH DIRECTION
 
 REAL(S), PARAMETER:: SIGMA  = 1.0          !MOLECULAR DIAMETER
 REAL(S), PARAMETER:: PI     = ACOS(-1.0)   !PI NUMBER

 INTEGER:: NPART,NCEL,TCEL
 INTEGER:: NEIGH(NCELX**2,8)                !NUMBER OF CELLS AND NEIGHBOURS OF EACH ONE

 REAL(S):: BOXX,BOXY
 REAL(S):: RX(NPARTX),RY(NPARTX)
 REAL(S):: ANG(NPARTX)

END MODULE VICVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM VICSEK
 USE VICVAR
 IMPLICIT NONE
 
 CALL BEGINSIM

 STOP
END PROGRAM VICSEK
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE BEGINSIM
 IMPLICIT NONE
 
 CALL READINP                               !READ SIMULATION PARAMETERS
 CALL ICONFIG                               !BUILD RANDON INITIAL CONFIGURATION
 CALL SETPARM                               !SET SIMULATION PARAMETERS
 CALL NEIGHBR                               !SET A MAP OF NEIGHBOURS OF EACH CELL

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
 
