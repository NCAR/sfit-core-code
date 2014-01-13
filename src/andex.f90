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
