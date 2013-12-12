! ----------------------------------------------------------------
!
      REAL (8) FUNCTION ANDEX (H,SH,GAMMA)
!
!     DOUBLE PRECISION VERSION OF ANDEX - NEEDED FOR IMPROVED GEOMETRY
!
!     *****************************************************************
!     COMPUTES THE INDEX OF REFRACTION AT HEIGHT H, SH IS THE
!     SCALE HEIGHT, GAMMA IS THE VALUE AT H=0 OF THE REFRACTIVITY =
!     INDEX-1
!     *****************************************************************
!
      REAL (8) :: H, SH, GAMMA

      IF (SH.EQ.0.0) THEN
         ANDEX = 1.0D0 + GAMMA
      ELSE
         ANDEX = 1.0D0 + GAMMA*EXP(-H/SH)
      ENDIF
!
      RETURN
!
      END FUNCTION ANDEX
