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

      MODULE MATRIX

      USE params
      USE DATAFILES

      IMPLICIT NONE

      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE ADDSUB(A, B, C, IM, JM, ISGN)

      INTEGER, INTENT(IN)        :: IM, JM, ISGN
      REAL(DOUBLE), INTENT(IN)   :: A(IM,JM), B(IM,JM)
      REAL(DOUBLE), INTENT(OUT)  :: C(IM,JM)

      !INTEGER :: I, J

      !print *, im, jm
      !do i=1,im
      !   write(0,'(10es16.4)')( a(i,j), b(i,j), j=1, jm )
      !   enddo

      C = A + B * ISGN

      RETURN

      END SUBROUTINE ADDSUB


!----------------------------------------------------------------------
      SUBROUTINE MULT(A, B, C, IM, JM, KM)

      INTEGER, INTENT(IN)        :: IM, JM, KM
      REAL(DOUBLE), INTENT(IN)   :: A(IM,JM), B(JM,KM)
      REAL(DOUBLE), INTENT(OUT)  :: C(IM,KM)

!      INTEGER :: I, K

!      DO I = 1, IM
!         DO K = 1, KM
!            !C(I,K) = 0.D0
!            !C(I,K) = C(I,K) + SUM(A(I,:)*B(:,K))
!            C(I,K) = SUM( A(I,:) * B(:,K) )
!         END DO
!      END DO

      !write(0,*) im, jm, km
      !stop
      C = MATMUL( A, B )

      RETURN

      END SUBROUTINE MULT


!----------------------------------------------------------------------
     SUBROUTINE MULTDIAG(A, V, AV, N, M)

      INTEGER, INTENT(IN)        :: N, M
      REAL(DOUBLE), INTENT(IN)   :: A(N,M), V(M)
      REAL(DOUBLE), INTENT(OUT)  :: AV(N,M)

      INTEGER :: J

      DO J = 1, M
         AV(:,J) = A(:,J)*V(J)
      END DO

      RETURN

      END SUBROUTINE MULTDIAG


!----------------------------------------------------------------------
     SUBROUTINE TRNSPS(A, AT, IM, JM)

      INTEGER, INTENT(IN)        :: IM, JM
      REAL(DOUBLE), INTENT(IN)   :: A(JM,IM)
      REAL(DOUBLE), INTENT(OUT)  :: AT(IM,JM)

      AT = TRANSPOSE(A)

      RETURN

      END SUBROUTINE TRNSPS


!----------------------------------------------------------------------
      SUBROUTINE INVRT(A, AM1, IM)

      INTEGER, INTENT(IN)         :: IM
      REAL(DOUBLE), INTENT(IN)    :: A(IM,IM)
      REAL(DOUBLE), INTENT(OUT)   :: AM1(IM,IM)

      INTEGER, DIMENSION(IM)      :: IWK1
      INTEGER, DIMENSION(IM*2)    :: IWK2
!      INTEGER                     :: INFO
!      REAL(DOUBLE)                :: ANORM, RCOND, DLANGE, SCOND, DET
!      REAL(DOUBLE)                :: WORK(IM*IM), ATEST( 41, 41 )

! --- UNCOMMENT TO USE LOCAL MATINV
      INTEGER :: IZERO, IDET, IDUM
      REAL(DOUBLE), DIMENSION(IM,IM) :: DUM
      REAL(DOUBLE) :: DET

      !INTEGER :: I, J

!      do i=1,im
!      write(88,'(1000ES26.18)') (a(i,j), j=1,im)
!      enddo

      AM1(:IM,:IM) = A(:IM,:IM)

!      WRITE(33,*) 'Input matrix'
!      DO I=1,IM
!      WRITE(33,'(1000ES26.18)') (AM1(I,J), J=1,IM)
!      ENDDO

! --- UNCOMMENT TO USE LOCAL MATINV
      IZERO = 0
      IDET = 1
      IDUM = 0
      DUM(:IM,:IM) = 0.0D0
      DET = 0.0D0

      CALL MATINV (IM, IM, AM1, IZERO, DUM, IDET, DET, IDUM, IWK1, IWK2)

      IF (DET == 0.D0) THEN
         WRITE (16, *) ' INVRT: ERROR... MATRIX IS SINGULAR...EXITING'
         WRITE ( 0, *) ' INVRT: ERROR... MATRIX IS SINGULAR...EXITING'
         CALL SHUTDOWN
         STOP 3
      ENDIF

      RETURN

! USE DIRECT CALLS TO LAPACK

!  --- COMPUTE NORM OF MATRIX
!      ANORM = DLANGE( "M", IM, IM, AM1, IM, WORK )
!      WRITE(*,10) "Max value : ", ANORM

!      ANORM = DLANGE( "F", IM, IM, AM1, IM, WORK )
!      WRITE(*,10) "Frobenius norm : ", ANORM

!      ANORM = DLANGE( "1", IM, IM, AM1, IM, WORK )
!      WRITE(*,10) "1 norm : ", ANORM

!      ANORM = DLANGE( "I", IM, IM, AM1, IM, WORK )
      !WRITE(*,10) "infinity norm : ", ANORM

!      AM1(:IM,:IM) = AM1(:IM,:IM) / ANORM


!  --- CHECK SCALING
!		CALL DPOEQU( IM, AM1, IM, WORK, SCOND, AMAX, INFO )
!      WRITE(*,10) " DPOEQU : INFO : ", SCOND, AMAX, INFO
!      write(*,*) (work(i),i=1,im)


!  --- COMPUTE CHOLESKY DECOMPOSITION FACTORS - symmetric positive definite
      ! IWK1 is ipiv
!      CALL DPOTRF( "U", IM, AM1, IM, INFO )
      !WRITE(*,11) iwk1
!      WRITE(*,10) " DPOTRF : INFO : ", INFO

!      WRITE(33,*) 'After Cholesky decomp'
!      DO I=1,IM
!      WRITE(33,'(1000ES26.18)') (AM1(I,J), J=1,IM)
!      ENDDO

!  --- COMPUTE THE INVERSE OF A
!      CALL DPOTRI( "U", IM, AM1, IM, INFO )
!      WRITE(*,10) " DPOTRI : INFO : ", info

!  ---  COPY UPPER TRIANGLE MATRIX VALUES TO LOWER PART

!      DO I=1, IM
!         DO J=1, I
!            AM1(I,J) = AM1(J,I)
!         ENDDO
!      ENDDO

!      WRITE(33,*) 'Inverse'
!      DO I=1,IM
!      WRITE(33,'(1000ES26.18)') (AM1(I,J), J=1,IM)
!      ENDDO

!      AM1(:IM,:IM) = A(:IM,:IM)

!  --- CHECK SCALING
!		CALL DGEEQU( IM, AM1, IM, WORK, SCOND, AMAX, INFO )
!      WRITE(*,11) " DGEEQU : INFO : ", SCOND, AMAX, INFO
!      write(*,*) (work(i),i=1,im)


!  --- COMPUTE LU DECOMPOSITION FACTORS - general
      ! IWK1 is ipiv
!      CALL DGETRF( IM, IM, AM1, IM, IWK1, INFO )
      !WRITE(*,11) iwk1
      !WRITE(*,12) " DGETRF : INFO : ", INFO

!  --- COMPUTE RECIPRICAL OF THE CONDITION NUMBER
!      CALL DGECON( "I", IM, AM1, IM, ANORM, RCOND, WORK, IWK2, INFO )
      !WRITE(*,13) "DGECON : INFO, RCOND : ", INFO, RCOND

!  --- COMPUTE THE INVERSE OF A
!      CALL DGETRI( IM, AM1, IM, IWK1, WORK, IM*IM, INFO )
      !WRITE(*,12) " DGETRI : INFO : ", INFO

!      WRITE(33,*) 'inverse'
!      DO I=1,IM
!      WRITE(33,'(1000ES26.18)') (AM1(I,J), J=1,IM)
!      ENDDO

!      ATEST(:,:) = A( 3:43, 3:43 )

!      IM = 41
!  --- COMPUTE LU DECOMPOSITION FACTORS
      ! IWK1 is ipiv
!      CALL DGETRF( IM, IM, ATEST, IM, IWK1, INFO )
      !WRITE(*,11) iwk1
!      WRITE(*,10) " DGETRF : INFO : ", INFO

!  --- COMPUTE RECIPRICAL OF THE CONDITION NUMBER
!      CALL DGECON( "I", IM, ATEST, IM, ANORM, RCOND, WORK, IWK2, INFO )
!      WRITE(*,10) "DGECON : INFO, RCOND : ", INFO, RCOND

!  --- COMPUTE THE INVERSE OF A
!      CALL DGETRI( IM, ATEST, IM, IWK1, WORK, IM*IM, INFO )
!      WRITE(*,10) " DGETRI : INFO : ", info


!      WRITE(33,*) 'funny inverse'
!      DO I=1,IM
!      WRITE(33,'(1000ES26.18)') (ATEST(I,J), J=1,IM)
!      ENDDO

!      STOP


      RETURN

! 10   FORMAT( A32, 5ES17.6 )
! 11   FORMAT( A32, 2ES17.6, I6 )
! 12   FORMAT( A32, 6I7 )
! 13   FORMAT( A32, I7, 2ES17.6 )

      END SUBROUTINE INVRT


!----------------------------------------------------------------------
      SUBROUTINE MATINV(MAX, N, A, M, B, IOP, DETERM, ISCALE, IPIVOT, IWK)

!     F1.3
!***********************************************************************
!
!     PURPOSE - MATINV INVERTS A REAL SQUARE MATRIX A.
!               IN ADDITION THE ROUTINE SOLVES THE MATRIX
!               EQUATION AX=B,WHERE B IS A MATRIX OF CONSTANT
!               VECTORS. THERE IS ALSO AN OPTION TO HAVE THE
!               DETERMINANT EVALUATED. IF THE INVERSE IS NOT
!               NEEDED, USE GELIM TO SOLVE A SYSTEM OF SIMULTANEOUS
!               EQUATIONS AND DETFAC TO EVALUATE A DETERMINANT
!               FOR SAVING TIME AND STORAGE.
!
!     USE     - CALL MATINV(MAX,N,A,M,B,IOP,DETERM,ISCALE,IPIVOT,IWK)
!
!                 MAX - THE MAXIMUM ORDER OF A AS STATED IN THE
!                       DIMENSION STATEMENT OF THE CALLING PROGRAM.
!
!                 N   - THE ORDER OF A, 1.LE.N.LE.MAX.
!
!                 A   - A TWO-DIMENSIONAL ARRAY OF THE COEFFICIENTS.
!                       ON RETURN TO THE CALLING PROGRAM, A INVERSE
!                       IS STORED IN A.
!                       A MUST BE DIMENSIONED IN THE CALLING PROGRAM
!                       WITH FIRST DIMENSION MAX AND SECOND DIMENSION
!                       AT LEAST N.
!
!                 M   - THE NUMBER OF COLUMN VECTORS IN B.
!                       M=0 SIGNALS THAT THE SUBROUTINE IS
!                       USED SOLELY FOR INVERSION,HOWEVER,
!                       IN THE CALL STATEMENT AN ENTRY CORRE-
!                       SPONDING TO B MUST BE PRESENT.
!
!                 B   - A TWO-DIMENSIONAL ARRAY OF THE CONSTANT
!                       VECTOR B. ON RETURN TO CALLING PROGRAM,
!                       X IS STORED IN B. B SHOULD HAVE ITS FIRST
!                       DIMENSION MAX AND ITS SECOND AT LEAST M.
!
!                 IOP - COMPUTE DETERMINANT OPTION.
!                        IOP=0 COMPUTES THE MATRIX INVERSE AND
!                              DETERMINANT.
!                        IOP=1 COMPUTES THE MATRIX INVERSE ONLY.
!
!               DETERM- FOR IOP=0-IN CONJUNCTION WITH ISCALE
!                       REPRESENTS THE VALUE OF THE DETERMINANT
!                       OF A, DET(A),AS FOLLOWS.
!                        DET(A)=(DETERM)(10**100(ISCALE))
!                       THE COMPUTATION DET(A) SHOULD NOT BE
!                       ATTEMPTED IN THE USER PROGRAM SINCE IF
!                       THE ORDER OF A IS LARGER AND/OR THE
!                       MAGNITUDE OF ITS ELEMENTS ARE LARGE(SMALL),
!                       THE DET(A) CALCULATION MAY CAUSE OVERFLOW
!                     (UNDERFLOW). DETERM SET TO ZERO FOR
!                     SINGULAR MATRIX CONDITION, FOR EITHER
!                     I0P=1,OR 0. SHOULD BE CHECKED BY PROGRAMER
!                     ON RETURN TO MAIN PROGRAM.
!
!             ISCALE  - A SCALE FACTOR COMPUTED BY THE
!                       SUBROUTINE TO AVOID OVERFLOW OR
!                       UNDERFLOW IN THE COMPUTATION OF
!                       THE QUANTITY,DETERM.
!
!             IPIVOT  - A ONE DIMENSIONAL INTEGER ARRAY
!                       USED BY THE SUBPROGRAM TO STORE
!                       PIVOTOL INFORMATION. IT SHOULD BE
!                       DIMENSIONED AT LEAST N. IN GENERAL
!                       THE USER DOES NOT NEED TO MAKE USE
!                       OF THIS ARRAY.
!
!                IWK  - A TWO-DIMENSIONAL INTEGER ARRAY OF
!                       TEMPORARY STORAGE USED BY THE ROUTINE.
!                       IWK SHOULD HAVE ITS FIRST DIMENSION
!                       MAX, AND ITS SECOND 2.
!
!     REQUIRED ROUTINES-
!
!     REFERENCE        -FOX,L, AN INTRODUCTION TO NUMERICAL
!                              LINEAR ALGEBRA
!
!     STORAGE          - 542 OCTAL LOCATIONS
!
!     LANGUAGE         -FORTRAN
!     LIBRARY FUNCTIONS -ABS
!
!     RELEASED          - JULY 1973
!
!     LATEST REVISION   - JULY 29, 1981
!                         COMPUTER SCIECES CORPORATION
!                         HAMPTON, VA
!***********************************************************************

      INTEGER :: MAX
      INTEGER :: N
      INTEGER :: M
      INTEGER :: IOP
      INTEGER :: ISCALE
      REAL(DOUBLE) :: DETERM
      INTEGER :: IPIVOT(N)
      INTEGER :: IWK(MAX,2)
      REAL(DOUBLE) :: A(MAX,N)
      REAL(DOUBLE) :: B(MAX,N)

      INTEGER :: JCOLUM, ICOLUM, JROW, IROW, J, I, K, L, L1
      REAL(DOUBLE) :: SWAP, T, AMAX, R1, R2, TMAX, PIVOT, PIVOTI

      EQUIVALENCE (IROW, JROW), (ICOLUM, JCOLUM), (AMAX, T, SWAP)
!
!     INITIALIZATION
!
      ISCALE = 0
!      R1=10.0**100.D0
      R1 = 10.0D0**100.D0
!      R2=1.0/R1
      R2 = 1.0D0/R1
      DETERM = 1.0D0
      IPIVOT = 0
      DO I = 1, N
!
!       SEARCH FOR PIVOT ELEMENT
!
         AMAX = 0.0D0
         DO J = 1, N
            IF (IPIVOT(J) - 1 == 0) CYCLE
            DO K = 1, N
               IF (IPIVOT(K) - 1 > 0) GO TO 740
               IF (IPIVOT(K) - 1 == 0) CYCLE
               TMAX = ABS(A(J,K))
               IF (AMAX - TMAX >= 0.D0) CYCLE
               IROW = J
               ICOLUM = K
               AMAX = TMAX
            END DO
         END DO
         IF (AMAX < 0.D0) GO TO 740
         IF (AMAX <= 0.D0) THEN
            DETERM = 0.0D0
            ISCALE = 0
            GO TO 740
         ENDIF
         IPIVOT(ICOLUM) = 1
!
!       INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
!
         IF (IROW - ICOLUM /= 0) THEN
            DETERM = -DETERM
            DO L = 1, N
               SWAP = A(IROW,L)
               A(IROW,L) = A(ICOLUM,L)
               A(ICOLUM,L) = SWAP
            END DO
            IF (M > 0) THEN
               DO L = 1, M
                  SWAP = B(IROW,L)
                  B(IROW,L) = B(ICOLUM,L)
                  B(ICOLUM,L) = SWAP
               END DO
            ENDIF
         ENDIF
         IWK(I,1) = IROW
         IWK(I,2) = ICOLUM
         PIVOT = A(ICOLUM,ICOLUM)
         IF (IOP < 0) GO TO 740
!
!       SCALE THE DETERMINANT
!
         IF (IOP <= 0) THEN
            PIVOTI = PIVOT
            IF (ABS(DETERM) - R1 >= 0.D0) THEN
               DETERM = DETERM/R1
               ISCALE = ISCALE + 1
               IF (ABS(DETERM) - R1 < 0.D0) GO TO 1060
               DETERM = DETERM/R1
               ISCALE = ISCALE + 1
               GO TO 1060
            ENDIF
            IF (ABS(DETERM) - R2 <= 0.D0) THEN
               DETERM = DETERM*R1
               ISCALE = ISCALE - 1
               IF (ABS(DETERM) - R2 <= 0.D0) THEN
                  DETERM = DETERM*R1
                  ISCALE = ISCALE - 1
               ENDIF
            ENDIF
 1060       CONTINUE
            IF (ABS(PIVOTI) - R1 >= 0.D0) THEN
               PIVOTI = PIVOTI/R1
               ISCALE = ISCALE + 1
               IF (ABS(PIVOTI) - R1 < 0.D0) GO TO 320
               PIVOTI = PIVOTI/R1
               ISCALE = ISCALE + 1
               GO TO 320
            ENDIF
            IF (ABS(PIVOTI) - R2 <= 0.D0) THEN
               PIVOTI = PIVOTI*R1
               ISCALE = ISCALE - 1
               IF (ABS(PIVOTI) - R2 <= 0.D0) THEN
                  PIVOTI = PIVOTI*R1
                  ISCALE = ISCALE - 1
               ENDIF
            ENDIF
  320       CONTINUE
            DETERM = DETERM*PIVOTI
         ENDIF
!
!       DIVIDE PIVOT ROW BY PIVOT ELEMENT
!
         A(ICOLUM,ICOLUM) = 1.0D0
         A(ICOLUM,:) = A(ICOLUM,:)/PIVOT
         IF (M > 0) THEN
            B(ICOLUM,:M) = B(ICOLUM,:M)/PIVOT
         ENDIF
!
!       REDUCE NON-PIVOT ROWS
!
         DO L1 = 1, N
            IF (L1 - ICOLUM == 0) CYCLE
            T = A(L1,ICOLUM)
            A(L1,ICOLUM) = 0.0D0
            A(L1,:) = A(L1,:) - A(ICOLUM,:)*T
            IF (M <= 0) CYCLE
            B(L1,:M) = B(L1,:M) - B(ICOLUM,:M)*T
         END DO
      END DO
!
!     INTERCHANGE COLUMNS
!
      DO I = 1, N
         L = N + 1 - I
         IF (IWK(L,1) - IWK(L,2) == 0) CYCLE
         JROW = IWK(L,1)
         JCOLUM = IWK(L,2)
         DO K = 1, N
            SWAP = A(K,JROW)
            A(K,JROW) = A(K,JCOLUM)
            A(K,JCOLUM) = SWAP
         END DO
      END DO
  740 CONTINUE
      RETURN
      END SUBROUTINE MATINV

      END MODULE MATRIX
