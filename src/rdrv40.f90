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

      PROGRAM RDRV

      USE PARAMS
      USE SFIT4
      USE DATAFILES
      USE BINPUT_4_0

      IMPLICIT NONE

      CHARACTER (LEN=10)   :: ZTIME='          ', ZONE='          '
      CHARACTER (LEN=8)    :: CDATE='        '
      REAL(8)              :: START_TIME=0.0, END_TIME=0.0

      CALL FILESETUP

! --- DETAIL FILE
      OPEN(UNIT=16, FILE=TFILE(16), STATUS='UNKNOWN', ERR=101, POSITION='ASIS')

      CALL CPU_TIME (START_TIME)

      !WRITE (16, *) VERSION
      !WRITE (*, *) VERSION

      CALL DATE_AND_TIME (CDATE, ZTIME, ZONE)
      WRITE (TAG,*) TRIM(VERSION), ' RUNTIME:', CDATE(1:8), '-', ZTIME(1:2), ':', ZTIME(3:4), ':', ZTIME(5:6)
      !WRITE (TAG,*) TRIM(VERSION), TRIM(BUILDDATE),' RUNTIME:', CDATE(1:8), '-', ZTIME(1:2), ':', ZTIME(3:4), ':', ZTIME(5:6)
      !WRITE (TAG,*) TRIM(VERSION), TRIM(BUILDDATE), ' RUNTIME:', CDATE(5:6), '/', CDATE(7:8), '/', &
      !     & CDATE(1:4), '-', ZTIME(1:2), ':', ZTIME(3:4), ':', ZTIME(5:6)

      WRITE (16, *) TRIM(TAG)
      WRITE ( 6, *) TRIM(TAG)

      CALL READ_BINPUT('sfit4.ctl')

      CALL SFIT ( )

      CALL CPU_TIME (END_TIME)
      WRITE(16,*)''
      WRITE(16,*) 'RDRV: DONE. ELAPSED TIME = ', END_TIME - START_TIME
      CLOSE(16)

      PRINT *, 'RDRV: DONE. ELAPSED TIME = ', END_TIME - START_TIME

      STOP 0

 101  CONTINUE
      WRITE (0, 201) TFILE(16)
      CLOSE(16)
      STOP 1

 201  FORMAT(/,' DETAIL OUTPUT FILE OPEN ERROR-UNIT 16'/,' FILENAME: "',A,'"')
!   10 FORMAT(A)
!   20 FORMAT(A7)

      END PROGRAM RDRV
