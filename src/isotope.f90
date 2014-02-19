!-----------------------------------------------------------------------------
!    Copyright (c) 2013-2014 NDACC/IRWG
!    This file is part of sfit.
!
!    sfit is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    sfit is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with sfit.  If not, see <http://www.gnu.org/licenses/>
!-----------------------------------------------------------------------------

      MODULE ISOTOPE

      USE PARAMS
      USE MOLCPARAM
      USE RETVPARAM
      USE VIBFCN
      USE DATAFILES

      IMPLICIT NONE

      INTEGER                                :: NISOSEP = 0, NISOVMR
      INTEGER, DIMENSION(ISOMAX)             :: OLDID, OLDISO, NEWID, NEWISO
      REAL(DOUBLE), DIMENSION(LAYMAX,ISOMAX) :: NEWVMR
      REAL(DOUBLE)                           :: ISOSCALE(ISOMAX)
      LOGICAL                                :: USEISO = .FALSE.
      CHARACTER (LEN=7), DIMENSION(ISOMAX)   :: NEWNAME

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE RDISOFILE (LUN )

      ! --- LUN IS DETAIL OUTPUT FILE, WHEN CALLED FROM SFIT

      CHARACTER (LEN=7), DIMENSION(ISOMAX)      :: OLDNAME
      INTEGER, DIMENSION(ISOMAX)                :: NEWMODE
      INTEGER                                   :: I, J, K, LUN
      REAL(DOUBLE), DIMENSION(ISOMAX)           :: NEWMASS, NEWTDEP
      INTEGER, DIMENSION(2,MAXVIBVALS,ISOMAX)   :: NEWIVIB
      LOGICAL                                   :: ISFILE

      OLDID = -1
      OLDISO = -1
      NEWID = -1
      NEWISO = -1
      OLDNAME = 'OTHER'
      NEWNAME = 'OTHER'
      NEWVMR = 0.0D0
      NEWMASS = 0.0D0
      NEWMODE = 0

      IF(.NOT. USEISO)THEN
         WRITE(LUN,35)
         ISFILE= .FALSE.
         RETURN
      ENDIF

!  --- CHECK FOR FILE - ISOTOPE SEPARATION INPUT FILE
      INQUIRE(FILE=TFILE(9), EXIST=ISFILE)
      IF( .NOT. ISFILE )THEN
          WRITE(LUN,30) trim(TFILE(9))
          WRITE( 0, 30) trim(TFILE(9))
          CALL SHUTDOWN
          STOP 1
      ENDIF

      CALL FILEOPEN( 9, 3 )

      READ(9,*) NISOSEP, NISOVMR
      IF( NISOSEP .GT. ISOMAX )STOP 'TOO MANY ISOTOPE SEPARATIONS REQUESTED'
      DO I=1, NISOSEP
          READ(9,'(A7)') OLDNAME(I)
          READ(9,*) OLDID(I), OLDISO(I)
          READ(9,'(A7)') NEWNAME(I)
          READ(9,*) NEWID(I), NEWISO(I), NEWMASS(I), NEWMODE(I), NEWTDEP(I), ISOSCALE(I)
          !WRITE(0,*) NEWID(I), NEWISO(I), NEWMASS(I), NEWMODE(I), NEWTDEP(I), ISOSCALE(I)
          DO J=1, MAXVIBVALS
              DO K=1, 2
                  NEWIVIB(K,J,I) = 0
              ENDDO
          ENDDO
          !write(0,*) i, newmode(i)
          IF( NEWMODE(I) .GT. 0 )THEN
             READ(9,*) NEWIVIB(:2,:NEWMODE(I),I)
             READ(9,*) NEWVMR(:NISOVMR,I)
          ENDIF
          NHIISO(NEWID(I))  = 0
          XMASS(:,NEWID(I)) = 0
      ENDDO

      CALL FILECLOSE( 9, 2 )

      WRITE(LUN,31) NISOSEP
      WRITE(LUN,32)

      DO I=1, NISOSEP
          IF( TRIM(NAME(NEWID(I))) .NE. 'OTHER' )THEN
              DO J=1, I-1
                 IF( TRIM(NAME(NEWID(J))) .EQ. TRIM(NAME(NEWID(I))) )THEN
                    WRITE(LUN,37) 'ADDING ISOTOPE ', NEWID(I), NEWISO(I), ' TO MOLECULE ', NEWID(J), NEWISO(J)
                    NEWVMR(:NISOVMR,I) = NEWVMR(:NISOVMR,J)
                    NEWIVIB(:2,:NEWMODE(J),I) = NEWIVIB(:2,:NEWMODE(J),J)
                    ! --- check that this name has one id
                    IF( NEWID(I) .NE. NEWID(J) )THEN
                       WRITE(LUN, 40) NAME(NEWID(J)), NEWID(I), NEWID(J)
                       WRITE(  0, 40) NAME(NEWID(J)), NEWID(I), NEWID(J)
                       CALL SHUTDOWN
                       STOP 1
                    ENDIF
                    ! --- check this same name & id has a different iso
                    IF( NEWISO(I) .EQ. NEWISO(J) )THEN
                       WRITE(LUN, 40) NAME(NEWID(J)), NEWID(I), NEWID(J)
                       WRITE(  0, 40) NAME(NEWID(J)), NEWID(I), NEWID(J)
                       CALL SHUTDOWN
                       STOP 1
                    ENDIF
                    EXIT
                 ENDIF
             ENDDO
             ! --- check that name does not clobber a predefined name
             IF( J .GT. I-1 )THEN
                WRITE(LUN, 41) NEWID(I), NAME(NEWID(I))
                WRITE(  0, 41) NEWID(I), NAME(NEWID(I))
                CALL SHUTDOWN
                STOP 1
             ENDIF
          ENDIF
          ! ---check we have no duplicate names
          DO J=1, MOLTOTAL
             IF( TRIM(NAME(NEWID(I))) .EQ. TRIM(NAME(J)) .AND. TRIM(NAME(J)) .NE. 'OTHER' )THEN
                IF( J .EQ. NEWID(I) )CYCLE
                WRITE(LUN, 42) NAME(J), J, NEWID(I)
                WRITE(  0, 42) NAME(J), J, NEWID(I)
                CALL SHUTDOWN
                STOP 1
             ENDIF
          ENDDO
          ! --- save new data
          NAME(NEWID(I))            = NEWNAME(I)
          NHIISO(NEWID(I))          = NHIISO(NEWID(I)) +1
          XMASS(NEWISO(I),NEWID(I)) = NEWMASS(I)
          IF( NEWMODE(I) .GT. 0 )THEN
             NMODE(NEWID(I))        = NEWMODE(I)
             TDEP(NEWID(I))         = NEWTDEP(I)
             DO J=1, MAXVIBVALS
                 DO K=1, 2
                     IVIB(K,J,NEWID(I)) = NEWIVIB(K,J,I)
                 ENDDO
             ENDDO
             WRITE(LUN,33) OLDNAME(I), OLDID(I), OLDISO(I), NAME(NEWID(I)), &
                          NEWID(I), NHIISO(NEWID(I)),  TDEP(NEWID(I)), ISOSCALE(I)
             WRITE(LUN,36) (IVIB(1:2,J,NEWID(I)),J=1,NMODE(NEWID(I)))
          ELSE
             WRITE(LUN,33) OLDNAME(I), OLDID(I), OLDISO(I), NAME(NEWID(I)), &
                          NEWID(I), NHIISO(NEWID(I))
          ENDIF
      ENDDO

      RETURN
  30  FORMAT(/,"ISOTOPE : NO ISOTOPE SEPARATION PARAMETER FILE, LUN : ", A)
  31  FORMAT(/,"NUMBER OF ISOTOPES TO SEPARATE =", I2)
  32  FORMAT(" OLD NAME, ID, ISO     NEW NAME, ID, ISO    TDEP   S-SCALE")
  33  FORMAT(2(1X,A7,3x,i2,3x,i2,4x),F5.2,2X,F8.6)
!  34  FORMAT(/,"ISOTOPE - COULD NOT OPEN ISOTOPE SEPARATION FILE :",A)
  35  FORMAT(/,"ISOTOPE - NO ISOTOPE SEPARATION SELECTED")
  36  FORMAT(6(I8,',',I2))
  37  FORMAT( 2(A, I3, '/', I2 ))
  40  FORMAT(/,"NEW MOLECULE ID'S FOR GAS : ", A, ' MUST MATCH : ', 2I5 )
  41  FORMAT(/,'MOLECULE ', I5, ' IS ALREADY IN USE AS ',A, ' CANNOT ASSIGN ISOTOPE TO EXISTING MOLECULE')
  42  FORMAT(/,'DUPLICATE MOLECULE NAMES : ', A, 2X, 2I5 )

      END SUBROUTINE RDISOFILE

      END MODULE ISOTOPE
