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

      MODULE DIAGNOSTIC

      USE OPT
      USE RETVPARAM

      IMPLICIT NONE
      REAL(DOUBLE) :: DOF(3)
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: G, SM, A

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE DOFS( M, N, NLD, NLEV )

      INTEGER, INTENT(IN)                :: N, M, NLD, NLEV
      INTEGER                            :: I, J
      REAL(DOUBLE), DIMENSION(N,M)       :: NM
      REAL(DOUBLE), DIMENSION(N,N)       :: IDEN

      !INTEGER                            :: INFO
      !REAL(DOUBLE), DIMENSION(4*NLEV)    :: WORK
      !REAL(DOUBLE), DIMENSION(NLEV)      :: WR, WI
      !REAL(DOUBLE), DIMENSION(NLEV,NLEV) :: H, Q


!  --- CALCULATE AVERAGING KERNEL NXN MATRIX
!  --- KSK IS AN NMAX*NMAX ARRAY

      ALLOCATE(SM(N, N), A(N,N), G(N,M))

      IF( F_LM )THEN
         CALL MULT( G_LM, KHAT, A, N, M, N)
         G(:N,:M) = G_LM(:N, :M)
      ELSE
         CALL MULT( SHAT(:N,:N), KSK, A, N, N, N)
         G(:N,:M) = MATMUL( SHAT(:N,:N), KS(:N,:M) )
      END IF

      if (F_WRTG) then
         CALL FILEOPEN( 93, 2 )
         WRITE(93,'(A,A)' )TRIM(TAG), ' GAIN MATRIX FOR FULL STATEVECTOR: '
         WRITE(93,*) N, M, ISMIX, NLEV
         WRITE(93,10) ADJUSTR(PNAME(:N))
         DO I=1,N
            WRITE(93,'(10000ES26.18)') ( G(I,J), J=1, M )
         ENDDO
         CALL FILECLOSE( 93, 1 )
      end if

!  --- CALCULATE TRACE OF AK KERNEL FOR ALL PARMS
      DOF(:) = 0.D0
      DO I = 1, N
         DOF(1) = DOF(1) + A(I,I)
      END DO

!  --- CALCULATE TRACE OF AK KERNEL FOR TARGET ONLY
      IF (IFPRF(1)) THEN
         DO I = NLD+1, NLD+NLEV
            DOF(2) = DOF(2) + A(I,I)
         END DO
      ELSE
         ! FOR COLUMN RETRIEVAL THE DOFS IS ALWAYS 1.0
         DOF(2) = 1.0
      END IF

!  --- CALCULATE TRACE OF AK KERNEL FOR TEMPERATURE IF RETRIEVED
      IF( IFTEMP )THEN
         DO I = NTEMP1, NTEMP1+NLEV
            DOF(3) = DOF(3) + A(I,I)
         END DO
      ENDIF

!  --- WRITE AK FOR TARGET
      ! IF COLUMN RETRIEVAL ONLY, THIS MATRIX IS USELESS.
      IF( F_WRTAK .AND. IFPRF(1) )THEN

         CALL FILEOPEN( 81, 1 )
         WRITE(81,'(A,A,A)' )TRIM(TAG), ' AVERAGING KERNELS NLEV X NLEV MATRIX FOR TARGET ONLY : ', GAS(1)
         WRITE(81,* )NLEV, NLEV
         DO I=1, NLEV
             WRITE(81,'(10000ES26.18)') ( A(I+NLD,J+NLD), J=1, NLEV )
         ENDDO
         CALL FILECLOSE( 81, 1 )

      ENDIF

!  --- CALCULATE MEASUREMENT ERROR : Se IS DIAGONAL
!  --- GAIN FOR COMPLETE STATEVECTOR
      CALL MULTDIAG( G(:N,:M), SE, NM, N, M )
      SM(:N,:N) = MATMUL( NM(:N,:M), TRANSPOSE( G(:N,:M) ))

      IF( F_WRTSMEAS .AND. IFPRF(1) )THEN

         CALL FILEOPEN( 82, 2 )
         WRITE(82,'(A,A,A)' )TRIM(TAG), ' MEASUREMENT ERROR WITH RETRIEVAL SE NLEV X NLEV MATRIX FOR TARGET ONLY : ', GAS(1)
         WRITE(82,*) NLEV, NLEV
         DO I=1,NLEV
            WRITE(82,'(10000ES26.18)') ( SM(I+NLD,J+NLD), J=1, NLEV )
         ENDDO
         CALL FILECLOSE( 82, 1 )

      ENDIF

!  --- CALCULATE SMOOTHING ERROR FOR TARGET GAS
!  --- FILL DIAGONAL ELEMENTS OF IDEN, CREATE IDENTITY MATRIX
!  --- REUSE MATRIX SM!

      IF ( F_WRTSMTH  .AND. IFPRF(1) )THEN

         IDEN = 0.0D0
         FORALL( I=1 : N ) IDEN(I,I) = 1.D0

         SM = MATMUL( (A-IDEN), SA(:N,:N) )
         SM = MATMUL( SM, TRANSPOSE( A-IDEN ))

         CALL FILEOPEN( 83, 2 )
         WRITE(83,'(A,A,A)' )TRIM(TAG), ' SMOOTHING ERROR (A-I)#SA#(A-I)T NLEV X NLEV MATRIX FOR TARGET ONLY : ', GAS(1)
         WRITE(83,*) NLEV, NLEV
         DO I=1,NLEV
            WRITE(83,'(10000ES26.18)') ( SM(I+NLD,J+NLD), J=1, NLEV )
         ENDDO
         CALL FILECLOSE( 83, 1 )

      ENDIF

      RETURN

!!  --- AK - SMOOTH
!
!      IF (F_AKSMOOTH .and. IFPRF(1) ) THEN
!         CALL FILEOPEN( 84, 2 )
!         WRITE(84,*) NLEV
!         DO I=1,NLEV
!            WRITE(84,'(10000ES26.18)') ( A(I+NLD,J+NLD)-SM(I+NLD,J+NLD), J=1, NLEV )
!         ENDDO
!         CALL FILECLOSE( 84, 1 )
!      END IF
!! --- COMMENT OUT TO USE LAPACK AND MAKE EIGENVECTORS
!      RETURN
!
!!  --- CALCULATE THE EIGENVALUES/VECTORS OF THE AVERAGING KERNEL
!!  --- REQUIRES COMPILING WITH LAPACK
!!  --- COPY A SINCE DGEEV WILL DESTROY CONTENTS
!!  --- EIGENVALUE DECOMPOSITION
!      H = A(NLD+1:NLEV,NLD+1:NLEV)
!!      CALL DGEEV('N','V',NLEV, H, NLEV, WR, WI, 0, 1, Q, NLEV, WORK, 4*NLEV, INFO)
!      !PRINT *, 'INFO :', INFO
!
!!  --- EIGENVALUES OF AK
!      CALL FILEOPEN( 85, 2 )
!      WRITE(85,*) NLEV, 'first 2 rows are real & imag eigenvalues'
!      WRITE(85,'(10000ES26.18)') ( WR(J), J=1, NLEV )      ! REAL E-VALUE
!      WRITE(85,'(10000ES26.18)') ( WI(J), J=1, NLEV )      ! IMAGINARY #-VALUE
!      DO I=1,NLEV
!          WRITE(85,'(10000ES26.18)') ( Q(I,J), J=1, NLEV )
!      ENDDO
!      CALL FILECLOSE( 85, 1 )

 10   FORMAT( 2000( 12X, A14 ))

      END SUBROUTINE DOFS

      SUBROUTINE RELEASE_MEM_DIA

      IF( ALLOCATED( G )) DEALLOCATE( G )
      IF( ALLOCATED( SM )) DEALLOCATE( SM )
      IF( ALLOCATED( A )) DEALLOCATE( A )


      END SUBROUTINE RELEASE_MEM_DIA

      END MODULE DIAGNOSTIC

