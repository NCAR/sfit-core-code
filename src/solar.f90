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

module solar

! Oct 2005
! Added code so we are equivalent to Hase solar 4.3
! Aug 2005
! Added code so we are equivalent to Hase solar 4.2
!
! July 2005
! Added ABS() to active window calculation so the solar emission line
!   at 811.57 will be seen.
!
! January 2005
! added FP underflow check line 490
!
! September 2004
!  - Moved parameter module globsol3 into this module and added CONTAIN
!     statement
!  - Moved contents of solar_dat_M to this file
!  - Update calculation to FH ver4 changed subroutine name to solarfh
!
!     Last change:  BLA  28 May 2004    3:17 pm
!
!
! Model of the infrared solar transmittance spectrum
!
! by
!
! F. Hase(1), P. Demoulin(2), A. Goldman(3), J. W. Hannigan(4),
! C. P. Rinsland(5), and G. C. Toon(6)
!     1: IMK, Forschungszentrum und Universitaet Karlsruhe, Germany,
!        <frank.hase@imk.fzk.de>
!     2: Institute de Astrophysique, Univ. Liege, Belgium
!     3: Department of Physics, University of Denver, Denver, CO, USA
!     4: National Center for Atmospheric Research, Boulder, CO, USA
!     5: NASA Langley Research Center, Hampton, VA
!     6: Jet Propulsion Laboratory, Pasadena, CA, USA

USE PARAMS
!USE DATAFILES
USE XSECTIONS
USE WRITEOUT

IMPLICIT NONE

LOGICAL                         :: IFCO
LOGICAL, DIMENSION(5)           :: F_RTSOL
LOGICAL                         :: DEC_DOPPLER
CHARACTER*(IFLNMSZ)             :: LINDATEI

INTEGER                         :: NPTSMW, NCOLNS=0

REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: TCO

REAL(DOUBLE), DIMENSION(LNMXCO) :: WSOLAR
REAL(DOUBLE), DIMENSION(LNMXCO) :: ADOP
REAL(DOUBLE), DIMENSION(LNMXCO) :: STR
REAL(DOUBLE), DIMENSION(5)      :: CPARM, SCPARM, CIPARM
REAL(DOUBLE), DIMENSION(MAXBND) :: SEMIEXTFOV, FOVDIA, DNUE, NUESKALCENTER, PHISOL, FOVRATIO

REAL(DOUBLE)                    :: FIRSTNUE, LONGRAD, LATRAD, JULD
REAL(DOUBLE)                    :: NUESKAL_INP

TYPE LOCATION
    SEQUENCE
    INTEGER :: YYYY, MO, DD, HH, MI
    REAL(8) :: SECS, NLAT, ELON
END TYPE

TYPE(LOCATION), DIMENSION(MAXBND,MAXSPE) :: LOCS

TYPE LINIENDATEN
    SEQUENCE
    REAL(8) :: NUE
    REAL(8) :: STRENGTH
    REAL(8) :: GAUSSWIDTH
    REAL(8) :: SKLOR
    REAL(8) :: VARISTRENGTH
    REAL(8) :: VARIWIDTH
    INTEGER :: ILEFT,IRIGHT
END TYPE

CONTAINS

subroutine solarfh( flag )

! CHANGES
!
! 7 JUN 2004
!  - converted to subroutine to sfit2.
!  - the plan is to loop this program (and all its subs) over NBAND and store
!    the transmission spectra consecutivly in TCO.
!  - added fit of shift parameter

!USE params
USE transmis
USE bandparam
!USE datafiles

implicit none

integer :: i, j, nfovmax !,zaehler
integer :: iband, kband, mxlnperbnd, flag

integer, dimension(MAXBND),save :: nlnread
integer, dimension(MAXBND) :: nlnfst, nlnlst
real(8) :: nueskal, radius, xpos, zpos, xznorm, weight, sumweight
real(8), dimension(:), allocatable :: soltrm_mean, soltrm

type(liniendaten),dimension(:,:), pointer, save :: lines
!type(liniendaten),dimension(:,:), allocatable, save :: lines

DEC_DOPPLER = .FALSE.
NUESKAL_INP = 1.0D0;

!  --- CLEAN UP AND EXIT
IF( FLAG == 2 ) THEN
    IF( ASSOCIATED(LINES)) DEALLOCATE (LINES)
    GOTO 100
ENDIF

! INITIALIZE SOLAR X-SECTION ARRAY
!print *, 'shape tco : ', shape(tco)
TCO(:NCROSS) = 0.D0              ! USE OLD SOLAR TRANSMISSION SPECTRUM VARIABLE
LINDATEI     = TFILE(10)         ! USE OLD LINE LIST FILESPEC

! do a quick trip through the solar line file and determine the # of lines in each band

IF ( FLAG == 0 ) THEN

    !JULD    = 2452000.0D0
    !LATRAD  = 0.0D0
    !LONGRAD = 0.0D0
    !WRITE(16,*)'SOLARFH - SETTING DEFAULT LAT, LON, DATE'

! --- determine relevant astronomical values
!       nueskalcenter : scaling of frequency axis due to Earth's orbital motion and rotation
!                    on center of solar disk
!       fovratio      : ratio FOV diam./solar diam.
!       phisol        : angle of solar axis towards earth
   WRITE(16,'(/,A)')' CALCULATION OF ASTRONOMICAL QUANTITIES FOR SOLAR SPECTRUM...'
   PRINT *,'   CALCULATION OF ASTRONOMICAL QUANTITIES ...'

   DO I = 1, NBAND ! FIT WINDOWS
      JULD    = 0.0D0
      LATRAD  = 0.0D0
      LONGRAD = 0.0D0

! --- WE ONLY CALCULATE ONE SOLAR SPECTRUM FOR A BAND
! --- SO WE ARE DISREGARDING THAT DIFFERENT SCANS( TIME ) ARE USING THE SAME SOLAR SPECTRUM
! --- AND TAKE AVERAGE UTC AND LOCATION OF SPECTRA FOR A BAND

      KBAND = NSCAN(I)
      DO J = 1, KBAND ! SOLAR ZENITH ANGLES
         WRITE(16,115) I, J, LOCS(I,J)%YYYY, LOCS(I,J)%MO, LOCS(I,J)%DD, LOCS(I,J)%HH, LOCS(I,J)%MI, &
                       LOCS(I,J)%SECS, LOCS(I,J)%NLAT, LOCS(I,J)%ELON
         JULD    = JULD + JULDAT( LOCS(I,J)%YYYY, LOCS(I,J)%MO, LOCS(I,J)%DD, LOCS(I,J)%HH, LOCS(I,J)%MI, LOCS(I,J)%SECS )
         LATRAD  = LATRAD + LOCS(I,J)%NLAT * PI /180.0D0
         LONGRAD = LONGRAD + LOCS(I,J)%ELON * PI /180.0D0
      ENDDO

      JULD    = JULD / REAL(KBAND,8)
      LATRAD  = LATRAD / REAL(KBAND,8)
      LONGRAD = LONGRAD / REAL(KBAND,8)

      DNUE(I)       = DN(I)                  ! MONOCHROMATIC POINT SPACING
      SEMIEXTFOV(I) = FOVDIA(I)/2000.        ! FOV HALF RADIUS IN RADIANS

      CALL EPHEMERID(JULD,LATRAD,LONGRAD,SEMIEXTFOV(I),NUESKALCENTER(I),FOVRATIO(I),PHISOL(I))

     WRITE(16,102) I, JULD, LATRAD*180.0D0/PI, LONGRAD*180.0D0/PI
     WRITE(16,103) 2.99792D5 * (1.0D0 - NUESKALCENTER(I))
     WRITE(16,104) NUESKALCENTER(I)
     WRITE(16,105) FOVRATIO(I)
     WRITE(16,106) 57.2958 * PHISOL(I)

   ENDDO

! --- for rest system
   IF (DEC_DOPPLER) NUESKALCENTER = NUESKAL_INP

      NCOLNS     = 0
      MXLNPERBND = 0

      DO IBAND = 1, NBAND

         FIRSTNUE = WMON(IBAND)                                ! INCLUDES 10 RESO UNITS FOR SHIFTS
         NPTSMW   = NM(IBAND)

         CALL LOOKATLINES( DNUE(IBAND), NLNFST(IBAND), NLNLST(IBAND) )

         NLNREAD(IBAND) = NLNLST(IBAND) - NLNFST(IBAND) + 1    ! LINE IN THIS BAND
         WRITE(*, 108) IBAND,NLNREAD(IBAND)
         WRITE(16,107) IBAND,NLNREAD(IBAND)
         NCOLNS     = NCOLNS + NLNREAD(IBAND)                  ! SUM OF ALL SOLAR LINES
         MXLNPERBND = MAX( MXLNPERBND, NLNREAD(IBAND) )

      END DO

      !PRINT *,'    MAX NUMBER OF LINES PER BAND: ',MXLNPERBND
      IF( NCOLNS == 0 )GOTO 100
      ALLOCATE( LINES( 1:MXLNPERBND, 1:NBAND ))

ENDIF ! FLAG = 0

!print *, '     size lines', size(lines)
!print *, '    someline   ', lines(1,1)%nue
KBAND = 1
BANDS: DO IBAND = 1, NBAND

   !print *, '    Solar spectrum calculation band: ',iband
   FIRSTNUE = WMON(IBAND)
   NPTSMW   = NM(IBAND)

!  --- Initialize spectral arrays
! (in soltrm the spectrum for a certain position on the solar disk will be written)
! (in soltrm_mean the resulting solar spectrum for given FOV will be written)
   !print *,'Allocating spectral arrays ...'

   ALLOCATE (SOLTRM_MEAN(1:NPTSMW),SOLTRM(1:NPTSMW))
   SOLTRM_MEAN = 0.0D0
   SOLTRM      = 0.0D0

   IF ( FLAG == 0 ) THEN
      CALL READLINES( NLNFST(IBAND), NLNLST(IBAND), LINES( :, IBAND ))
      CALL LINEBOUNDS( DNUE(IBAND), NLNFST(IBAND), NLNLST(IBAND), NUESKALCENTER(IBAND), LINES( :, IBAND ))
   ENDIF

!  --- FOV integration
! nfovmax determines FOV gridpoints: -nfovmax..+nfovmax x -nfovmax..+nfovmax
! xpos,zpos: determine position on solar disk in relative coordinates (right
! handed coordinate system, z-direction perpendicular to observer's LOS, aligned
! with projection of solar rotation axis on celestial sphere, y-axis along LOS

!print *,'      Performing FOV-integration ...'
!  --- Simplified FOV integration, assumes inclination of sun's equator to be zero

   !zaehler = 0
   NFOVMAX   = NINT(20.0D0 * FOVRATIO(IBAND))
   XZNORM    = FOVRATIO(IBAND) / (NFOVMAX + 0.5D0)
   SUMWEIGHT = 0.0D0
   DO I = -NFOVMAX, NFOVMAX

        XPOS    = I * XZNORM
        ZPOS    = 0.5D0 * XZNORM * SQRT(REAL(NFOVMAX * NFOVMAX - I * I,8))
        RADIUS  = SQRT(XPOS * XPOS + ZPOS * ZPOS)

        NUESKAL = NUESKALCENTER(IBAND) * NUESKALDISK(XPOS,ZPOS,0.0D0)
        !NUESKAL = NUESKALCENTER * NUESKALDISK(XPOS,ZPOS,PHISOL)

        CALL MAKESOLSPEC( DNUE(IBAND), RADIUS, NUESKAL, NLNREAD(IBAND), LINES( :, IBAND), SOLTRM )

        WEIGHT      = SQRT(1.0D0 - XPOS * XPOS)
        SOLTRM_MEAN = SOLTRM_MEAN + WEIGHT * SOLTRM
        SUMWEIGHT   = SUMWEIGHT + WEIGHT

   END DO

   SOLTRM_MEAN = SOLTRM_MEAN / SUMWEIGHT

   IF (F_WRTSOLSPEC) THEN
      CALL TOFILE_SPEC( IBAND, NPTSMW, FIRSTNUE, DNUE(IBAND), SOLTRM_MEAN )
   END IF
   !print *, 'size soltrm', size(soltrm_mean)
   !print *, kband, nm(iband), kband, nm(iband)+kband-1

    TCO( KBAND : KBAND + NM(IBAND) - 1) = SOLTRM_MEAN
    KBAND = KBAND + NM(IBAND)

   !print *, '     size lines', size(lines)
   !Deallocierung der Felder
   deallocate (soltrm_mean,soltrm)

END DO BANDS

!IF( FLAG == 0 ) PRINT *,' SOLAR SPECTRUM CALCULATION COMPLETED.'

100 CONTINUE

RETURN

!101 FORMAT( A, ES21.6 )
!101 FORMAT( A, ES21.6 )
!102 FORMAT( A, I4, 3f16.5 )
!103 FORMAT( A, I4, I16 )

102 FORMAT('   BAND, JULIAN DATE, N-LAT, E_LON    : ', I4, 3F17.5 )
103 FORMAT('   RADIAL VEL. (POS. RECEDING) [KM/S] : ', ES21.6 )
104 FORMAT('   DOPPLER EFFECT SCALE FACTOR        : ', ES21.6 )
105 FORMAT('   FOV IN RELATIVE TO SOLAR DIA.      : ', ES21.6 )
106 FORMAT('   SOLAR AXIS INCLINATION [DEG.]      : ', ES21.6 )
107 FORMAT('   BAND, # SOLAR LINES FOUND          : ', I4, I8 )
108 FORMAT('   BAND, # SOLAR LINES FOUND                                  :', I3, I7 )
115 FORMAT('   IBAND, ISPEC, OBSERVATION UTC/LOC  : ',2I4,2x,I4,'/',I2.2,'/',I2.2,2X,I2.2,':',I2.2,':',F04.1,F7.3,'N',F8.3,'W' )

END SUBROUTINE SOLARFH


!====================================================================
!  ephemerid: calculates Doppler scaling of solar spectrum due to
!                 orbital motion and rotation of Earth
!====================================================================
subroutine ephemerid(juld,latrad,longrad,semiextfov,nueskalcenter,fovratio,phisol)

implicit none

real(8),intent(in) :: juld,latrad,longrad,semiextfov
real(8),intent(out) :: nueskalcenter,fovratio,phisol

real(8),parameter :: pi = 3.141592653589793d0
real(8) :: zeit,anom,danom_dt,center,dcenter_dt,vrad_orbit,vrad_diurn &
  ,sideral,laen_true,alpha,delta,kwert,lnull


!print *, juld,latrad,longrad,semiextfov,nueskalcenter,fovratio,phisol


! Zeit vs. 1.Jan. 2000 1.5 TD in Jul. Jahrhunderten
zeit = (juld - 2.451545d6) / 3.6525d4

!====================================================================
! Berechnung der Orbitalkomponente der Radialgeschwindigkeit

! Mittlere Anomalie in rad (und zeitl. Ableitung der m. A. in rad/s)
anom = 6.2400600d0 + zeit * (628.3019553d0 - zeit * (2.72100d-6 &
  + 8.38d-9 * zeit))
danom_dt = 1.99097d-7

! Mittelpunktsgleichung
center = (3.3416073d-2 - zeit * (8.4073d-5 + 2.44d-7 * zeit)) &
  * sin(anom) + (3.48944d-4 - 1.763d-6 * zeit) * sin(2.0d0 * anom) &
  + 5.061d-6 * sin(3.0d0 * anom)
dcenter_dt = danom_dt * (3.3416d-2 * cos(anom) + 3.489d-4 * cos(2.0d0 * anom))

! Radialgeschwindigkeit durch Erdbewegung (astr. Konvention!), m/s
vrad_orbit = 2.50d9 * (danom_dt + dcenter_dt) * sin(anom + center)

!====================================================================
! Berechnung der taeglichen Komponente der Radialgeschwindigkeit

! Mittlere Sternzeit Greenwich (bezogen auf mittl Frühlingspunkt) in rad
sideral = modulo(4.894961212d0 + 6.30038809781d0 * (juld - 2.451545d6) &
+ 6.770708d0 * zeit * zeit - 4.5087d-10 * zeit * zeit * zeit,2.0d0 * pi)

! wahre Laenge der Sonne in rad (Summe aus center und mittl. Laenge)
laen_true = center + 4.89506300d0 + 6.283319668d2 * zeit + 5.2918d-6 * zeit * zeit

! Rektaszension und Deklination der Sonne in rad
alpha = atan2(0.917523d0 * sin(laen_true),cos(laen_true))
delta = asin(0.39768d0 * sin(laen_true))

! Radialgeschwindigkeit durch Erddrehung (astr. Konvention!), m/s
vrad_diurn = 4.63d2 * cos(latrad) * cos(delta) * sin(sideral - longrad - alpha)

!====================================================================
! Berechnung des resultierenden Skalierungsfaktors

nueskalcenter = 1.0d0 - 3.33564d-9 * (vrad_diurn + vrad_orbit)

!====================================================================
! Berechnung des relativen Gesichtsfeldes

fovratio = 2.1486d2 * semiextfov * (1.0d0 + 1.671d-2 * cos(anom + center))
if (fovratio .gt. 1.0d0) fovratio = 1.0d0

!====================================================================
! Berechnung der Neigung des Sonnenaequators (genaehert)

! geometrische mittlere Laenge der Sonne
lnull = 4.8950630d0 + 6.283319668d2 * zeit - 5.2918d-6 * zeit * zeit
kwert = 1.28573d0 + 6.6699211d-7 * (juld - 2396758.0d0)
phisol = 0.1265d0 * sin(lnull + center - kwert)

end subroutine ephemerid



!====================================================================
!  gonext: Einlesen bis zum naechsten $ Zeichen
!====================================================================
subroutine gonext(ifile)

implicit none
integer,intent(in) :: ifile
character(1) :: nextchar

nextchar='x'
do while (nextchar /= '$')
    read(ifile,'(a1)') nextchar
end do

end subroutine gonext



!====================================================================
!  juldat: calculates Julian Date
!====================================================================
real(8) function juldat (jahr_in,monat_in,tag,stunde,minute,sekunde)

implicit none

integer,intent(in) :: jahr_in,monat_in,tag,stunde,minute
real(8),intent(in) :: sekunde

integer :: jahr,monat,iwerta,iwertb
real(8) :: ftag

ftag = tag + stunde / 24.0d0 + minute / (1.440d3) + sekunde / (8.6400d4)

if (monat_in .le. 2) then
    jahr = jahr_in - 1
    monat = monat_in + 12
else
    jahr = jahr_in
    monat = monat_in
end if

iwerta = int(0.01 * jahr)
iwertb = 2 - iwerta + int(0.25d0 * iwerta)

juldat = int(3.6525d2 * (jahr + 4716)) + int(3.06001d1 * (monat + 1)) &
 + ftag + iwertb - 1.5245d3

end function juldat



!====================================================================
!  LineBounds: determines actve window for line calculation
!====================================================================
subroutine LineBounds(dnu,nfirstline,nlastline,nueskal,lines)


implicit none

integer,intent(in) :: nfirstline,nlastline
type(liniendaten),dimension(nlastline-nfirstline+1),intent(inout) :: lines
real(8),intent(in) :: nueskal, dnu

integer :: i,nfirst,nlast
real(8) :: nueradius

do i = 1,nlastline-nfirstline+1
    nueradius = 4.0d0 * (abs(lines(i)%strength) + 0.25d0) &
      * (9.0d0 * lines(i)%gausswidth) * (1.0d0 + 5.0d0 * abs(lines(i)%sklor))
    nfirst = nint((lines(i)%nue * nueskal - firstnue - nueradius) / dnu)
    nlast  = nint((lines(i)%nue * nueskal - firstnue + nueradius) / dnu)
    lines(i)%ileft = nfirst
    lines(i)%iright = nlast
    if (lines(i)%ileft .lt. 1) then
        lines(i)%ileft = 1
    else
        if (lines(i)%ileft > nptsmw) lines(i)%ileft = nptsmw
    end if
    if (lines(i)%iright > nptsmw) then
        lines(i)%iright = nptsmw
    else
        if (lines(i)%iright < 1) lines(i)%iright = 1
    end if
end do

end subroutine LineBounds



!====================================================================
!  LookAtLines: determines number of lines in database
!====================================================================
subroutine LookAtLines( dnu, nfirstline, nlastline )


implicit none

integer,intent(out) :: nfirstline,nlastline

!character(len=85) :: dumchar
integer :: zaehler,status
real(8),parameter :: nueextradius = 20.0d0  ! Erweiterungsbereich fuer Linienauswahl
real(8) :: nue, dnu

!open (10,file=trim(lindatei), iostat = status, status='old',action='read')
!if (status /= 0) print *,'Cannot find solar linedata file'

WRITE (16, '(/A,A)') ' SOLAR LINELIST FILE :', TFILE(10)
CALL FILEOPEN( 10, 3 )

call gonext(10)

zaehler = 0
outer: do
    read (10,'(F12.6,A76)',end = 100,iostat = status) nue !,dumchar
    if (status /= 0) print *,'Sub lookatlines: Cannot read solar lines!'
        zaehler = zaehler + 1
        if (firstnue - nueextradius .lt. nue) then
            nfirstline = zaehler
            do
                read (10,'(F12.6,A76)',end = 100,iostat = status) nue !,dumchar
                if (status /= 0) print *,'Sub lookatlines: Cannot read solar lines!'
                zaehler = zaehler + 1
                if (firstnue + dnu * (nptsmw-1) + nueextradius .lt. nue) then
                    zaehler = zaehler - 1
                    exit outer
                end if
            end do
        end if
end do outer
100 continue

CALL FILECLOSE( 10, 2 )

nlastline = zaehler

end subroutine LookAtLines



!====================================================================
!  MakeSolspec: calculates solar spectrum for specified position
!                on solar disk
!====================================================================
subroutine MakeSolspec( dnu, radius, nueskal, noflines, lines, soltrm )

!USE solar_dat_M

implicit none

integer,intent(in) :: noflines
real(8),intent(in) :: radius,nueskal
type(liniendaten),dimension(noflines),intent(in) :: lines

real(8),dimension(nptsmw),intent(out) :: soltrm

integer :: i,j
real(8) :: linstrength,width,widthqu,widthkub,widthquad,deltanue,restwert
real(8) :: linstr, wshift, nue, tup, tlow, exarg, dnu

! DECODE FITTING PARAMETERS
WSHIFT = CPARM(4) !-1.0d0
LINSTR = CPARM(5) !-1.0d0
!print*, wshift, linstr

soltrm = 0.0d0

! Gesamte optische Dicke aus Beitraegen aller Linien bestimmen
do i = 1,noflines
    if (lines(i)%iright - lines(i)%ileft .gt. 0) then
        linstrength = lines(i)%strength * (1.0d0 + linstr + lines(i)%variStrength * radius * radius)
        width = lines(i)%gausswidth * (1.0d0 + lines(i)%variWidth * radius * radius)
        widthqu = width * width
        widthkub = widthqu * width
        widthquad =  widthkub * width
!print*, linstrength, widthquad
        do j = lines(i)%ileft,lines(i)%iright
            deltanue = abs(firstnue + dnu * (j-1) - nueskal * ( lines(i)%nue + wshift ))
            exarg = deltanue * deltanue / sqrt(widthquad + abs(lines(i)%sklor) * (-0.54d0 * widthquad + deltanue *  &
                    (0.33d0 * widthkub + deltanue * (0.12d0 * widthqu + 0.342d0 * deltanue * width))))
            exarg = MIN(exarg, 664.0d0)
            soltrm(j) = soltrm(j) + linstrength * exp( -exarg )
        end do
    end if
end do

! Uebergang auf Transmissionsspektrum, Minnaert-Offset
! Minnaert-Offset for extended wave number range (valid from 700 cm-1 to 25000 cm-1)

do i = 1,nptsmw
    nue = firstnue + dnu * (i - 1)
    tup = 5700.0d0 + 0.28d0 * (nue - 2000.0d0)
    tlow = 4000.0d0
    restwert =  (exp(1.4387d0 * nue / tup) - 1.0d0) / (exp(1.4387d0 * nue / tlow) - 1.0d0)
    soltrm(i) = restwert + (1.0d0 - restwert) * exp(-soltrm(i))
end do

return

end subroutine MakeSolspec



!====================================================================
!  NueSkalDisk:
!====================================================================
real(8) function NueSkalDisk(xpos,zpos,phisol)

implicit none

real(8),intent(in) :: xpos,zpos,phisol

real(8),parameter :: pi = 3.141592653589793d0
real(8) :: latsol,omega

if (xpos * xpos + zpos * zpos .gt. 1.0d0) then
    NueSkalDisk = 1.0d0
else
    ! determine heliographic latitude associated with xpos,zpos:
    latsol = 0.5d0 * pi - acos(zpos * cos(phisol) - sin(phisol) * &
      sqrt(1.0d0 - xpos * xpos - zpos * zpos))

    ! determine solar rotation frequency to this latitude [rad/s]
    omega = 2.0d0 * pi * (4.2711d-7 + &
      (-9.599d-8 + 1.276d-8 * latsol * latsol) * latsol * latsol)

    ! determine scaling factor
    NueSkalDisk = (1.0d0 - 2.3225d0 * omega * cos(phisol) * xpos)
end if

end function NueSkalDisk



!====================================================================
!  ReadInput: Einlesen der Linien
!====================================================================
subroutine ReadInput(dateiname)


implicit none

character(len=*),intent(in) :: dateiname

open (10,file = trim(dateiname),status = 'old')
call gonext(10)
read (10,*) lindatei

call gonext(10)
read (10,*) firstnue
read (10,*) dnue
read (10,*) nptsmw

call gonext(10)
read (10,*) juld
read (10,*) latrad
read (10,*) longrad
read (10,*) dec_doppler
if (dec_doppler) read (10,*) nueskal_inp

call gonext(10)
read (10,*) semiextfov

close (10)

end subroutine ReadInput



!====================================================================
!  Read solar lines
!====================================================================
subroutine ReadLines( nfirstline, nlastline, lines )


implicit none

integer,intent(in) :: nfirstline,nlastline
type(liniendaten),dimension(1:nlastline-nfirstline+1),intent(inout) :: lines

integer :: i,status
!real(8) :: dumreal

open (10,file=trim(lindatei), iostat = status, status='old',action='read')
if (status /= 0) print *,'Sub readsolarlines: Error in opening file'

do i = 1,nfirstline-1
    read (10,*)
end do

do i = 1, nlastline - nfirstline + 1
    read (10,'(F12.6,ES11.3,ES11.3,ES11.3,ES11.3,ES11.3)',iostat = status) &
       & lines(i)%nue,lines(i)%strength,lines(i)%gausswidth, &
       & lines(i)%sklor,lines(i)%variStrength,lines(i)%variWidth
    if (status /= 0) print *,'Sub readsolarlines: Cannot read solar line!'
end do
 close(10)

end subroutine readlines



!====================================================================
!  tofile_specs: Rausschreiben des Spektrums (Fuer Test-Output)
!====================================================================
SUBROUTINE TOFILE_SPEC(IBAND,NPTS,FIRSTNUE,DNUE,SPEC)

      IMPLICIT NONE

      !CHARACTER(*),INTENT(IN) :: DATEINAME
      INTEGER,INTENT(IN) :: NPTS, IBAND
      REAL(8),INTENT(IN) :: FIRSTNUE,DNUE
      REAL(8),DIMENSION(NPTS),INTENT(IN) :: SPEC

      INTEGER :: I

      IF( IBAND .EQ. 1 ) THEN
         CALL FILEOPEN( 11, 2 )
         WRITE(11,*) TRIM(TAG), ' CALCULATED SOLAR SPECTRUM'
      ELSE
         CALL FILEOPEN( 11, 4 )
      ENDIF

      WRITE(11,'(I10,2F16.6)') NPTS, FIRSTNUE, DNUE
      DO I = 1,NPTS
         WRITE (11,'(F12.6,1X,F7.4)') FIRSTNUE + DNUE * (I - 1), SPEC(I)
      END DO

      CALL FILECLOSE( 11, 1 )

END SUBROUTINE TOFILE_SPEC


end module solar

