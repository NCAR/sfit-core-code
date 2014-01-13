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


module hitran

   use params
   use datafiles
   use bandparam
   use lineparam
   use isotope
   use binput_4_0
   use binput_parse_4_0

   implicit none

   integer, parameter     :: nhit=99, ngal=2, ncia=2, nglines=100000, flagoff=280
   integer, parameter     :: nlmx=1, nsdv=1, nlmlines=100000, nsdlines=100000
   real(8), parameter     :: weps = 1.0d-6
   integer                :: hnml, gnml, lnml, snml
   integer                :: stlun, nfiles, map(nhit)

   type, public :: hitranfile
      integer             :: mo        ! molecule id from subdir name for this file --- replaces map()
      integer             :: flag      ! 0-hitran, 1-cia f&s
      integer             :: lun       ! unit # for this record
      character (len=300) :: buf       ! read buffer
   end type hitranfile

   type, public :: galatrydata
      integer            :: n          ! number of data lines / file(gas)
      integer            :: lun        ! unit for this file
      character (len=255) :: buf        ! read buffer
      integer            :: mo(nglines)    ! mol id
      integer            :: is(nglines)    ! isotope id #
!      real(8)            :: nu(nglines)    ! wavenumber
      real(8)            :: bt(nglines)    ! intensity [cm-1/(molec/cm-2)]
      character(len=60)  :: qa(nglines)    ! quanta data
      ! The format of the quanta fields depends on the species. Refer to
      ! Rothman, L. S. et.al.
      ! The HITRAN 2004 molecular spectroscopic database
      ! Journal Of Quantitative Spectroscopy & Radiative Transfer, 2005, 96, 139-204
      real(8)            :: g0_air(nglines)
      real(8)            :: td_g0_air(nglines)
      real(8)            :: beta(nglines)
      real(4)            :: g0_self(nglines)
      real(4)            :: s_air(nglines)
      real(4)            :: g2_air(nglines)
      real(4)            :: ts_air(nglines)
      real(4)            :: lm_air(nglines) ! line mixing coefficients
      real(4)            :: lm_t1(nglines)  ! extra parameters of F. Hase to model
      real(4)            :: lm_t2(nglines)  ! temperature dependency
   end type galatrydata

   type, public :: linemixfile
      integer            :: n          ! number of data lines / file(gas)
      integer            :: ist        ! record start
      integer            :: lun        ! unit for this file
      character (len=64) :: buf        ! read buffer
      integer            :: mo         ! mol id
   end type linemixfile

   type :: linemixdata
      integer            :: mo         ! mol id
      real(8)            :: nu         ! wavenumber
      real(4)            :: dt(3)      ! intensity [cm-1/(molec/cm-2)]
   end type linemixdata

! this map is obsolete
! the linelist directory structure is the key to the gas names and the molecule id numbers
! those id's and names must be the same in the reference.prf file
! eg a files containing hitran lines is read from one subdir in linelist then the molid will be changed
! to the 2digit integer 0NN of the subdir name and assumed to be for gas 0NN_abcdef
! --- CODE NUMBERS TO CONVERT FROM HITRAN TO ATMOS / SFIT
!     INDEX IS HITRAN ID NUMBER, VALUE IS SFIT
!     SFIT     HITRAN         SFIT        PSEUDOLINES FILE MOLID/ISO
!     #        # NAME         # NAME
      DATA MAP / &
      1,    &! 1 H2O
      2,    &! 2 CO2
      3,    &! 3 O3
      4,    &! 4 N2O
      5,    &! 5 CO
      6,    &! 6 CH4
      7,    &! 7 O2
      8,    &! 8 NO
      9,    &! 9 SO2
      10,   &! 10 NO2
      11,   &! 11 NH3
      12,   &! 12 HNO3
      13,   &! 13 OH
      14,   &! 14 HF
      15,   &! 15 HCL
      16,   &! 16 HBR
      17,   &! 17 HI
      18,   &! 18 CLO
      19,   &! 19 OCS
      20,   &! 20 H2CO
      21,   &! 21 HOCL
      41,   &! 22 N2          HO2
      28,   &! 23 HCN         H202
      30,   &! 24 CH3CL       HONO
      23,   &! 25 H2O2        HO2NO2
      40,   &! 26 C2H2        N2O5
      38,   &! 27 C2H6        CLONO2
      0,    &! 28 PH3         HCN
      36,   &! 29 COF2        CH3F
      50,   &! 30 SF6         CH3CL
      47,   &! 31 H2S         CF4
      46,   &! 32 HCOOH       CCL2F2
      22,   &! 33 HO2         CCL3F3
      0,    &! 34 O           CH3CCL3
      35,   &! 35 CLONO2      CCL4          PS 35/1
      0,    &! 36 NO+         COF2
      0,    &! 37 HOBR        COCLF
      39,   &! 38 C2H4        C2H6
      64,   &! 39 CH3OH*      C2H4
      44,   &! 40 CH3BR*      C2H2
      69,   &! 41 CH3CN*      N2
      31,   &! 42 CF4*        CHF2CL
      43,   &! 43             COCL2
      44,   &! 44             CH3BR
      45,   &! 45             CH3I
      46,   &! 46             HCOOH
      47,   &! 47             H2S
      48,   &! 48             CHCL2F
      49,   &! 49             O2CIA
      50,   &! 50             SF6
      51,   &! 51             NF3
      52,   &! 52             OTHER
      53,   &! 53             OTHER
      54,   &! 54             OTHER
      55,   &! 55             OTHER
      56,   &! 56             OTHER
      57,   &! 57             OTHER
      58,   &! 58             OCLO
      59,   &! 59             F134A
      60,   &! 60             C3H8
      61,   &! 61             F142B
      62,   &! 62             CFC113
      63,   &! 63             F141B
      64,   &! 64             CH3OH
      65,   &! 65             CH3CNPL
      66,   &! 66             C2H6PL
      67,   &! 67             PAN
      68,   &! 68             CH3CHO
      69,   &! 69             CH3CN
      70,   &! 70             OTHER
      71,   &! 71             OTHER
      72,   &! 72             OTHER
      73,   &! 73             OTHER
      74,   &! 74             OTHER
      75,   &! 75             OTHER
      76,   &! 76             OTHER
      77,   &! 77             OTHER
      78,   &! 78             OTHER
      79,   &! 79             OTHER
      80,   &! 80             OTHER
      81,   &! 81             OTHER
      82,   &! 82             OTHER
      83,   &! 83             OTHER
      84,   &! 84             OTHER
      85,   &! 85             OTHER
      86,   &! 86             OTHER
      87,   &! 87             OTHER
      88,   &! 88             OTHER
      89,   &! 89             OTHER
      90,   &! 90             OTHER
      91,   &! 91             OTHER
      92,   &! 92             OTHER
      93,   &! 93             OTHER
      94,   &! 94             OTHER
      95,   &! 95             OTHER
      96,   &! 96             OTHER
      97,   &! 97             OTHER
      98,   &! 98             OTHER
      99    &! 99             OTHER
      /

! map to 1 record in the ascii linelist file made by hbin
! byte range  quantity
! 1   - 160   hitran - see HITRAN docs
! 161 - 172   galatry beta
! 173 - 184   sdv gamma0
! 185 - 196   sdv gamma2
! 197 - 208   sdv eta2
! 209 - 220   lmx ltk1
! 221 - 232   lmx ltk2
! 233 - 244   lmx ylm
! 245 - 256   ---
! 257 - 268   ---
! 269 - 280   ---
! 281 - 288   logical flags for: galatry, fcia, scia, sdv, lmix, 6, 7, 8
! 289 - 300   <nothing>

! from lineparam.f90: GALATRY_FLAG=1, FCIA_FLAG=2, SCIA_FLAG=3, SDV_FLAG=4, LM_FLAG=5

end module hitran

program hbin

! program to create a binary linelist file from initial HITRAN formatted files
!  for use by sfit4

   use hitran

   implicit none

   integer              :: ldx, nl, i, j, ifl
   integer              :: hblun=7, halun=8, istat, iband, inxt(1)
!   integer(long_log)    :: nulm, nuht
   real(double)         :: wavnum, wstr, wstp
   character (len=30)   :: hbfile, hafile
   character (len=200)  :: nam
   character (len=1)    :: pos
   character (len=255)  :: buf
   logical              :: oped, hasc=.TRUE.
   integer              :: iost, dum, ind
   character (len=7), dimension(6) :: sdv_params
   logical :: qu_equal
   character (len=10)   :: ztime='          ', zone='          '
   character (len=8)    :: cdate='        '


! --- we have beta data for 2 gases and not too many lines (sfit4 v0.9)

   type (hitrandata),  dimension(nhit+ncia)   :: hlp
   type (hitranfile),  dimension(nhit+ncia)   :: hfl

!   type (linemixfile), dimension(nlmx)        :: lfl
!   type (linemixdata), dimension(nlmlines)    :: lmx

   type (galatrydata), dimension(ngal)        :: glp
   type (galatrydata), dimension(nlmx)        :: lmx
   type (galatrydata), dimension(nsdv)        :: sdv

   call date_and_time (cdate, ztime, zone)
   write (tag,*) trim(version), ' runtime:', cdate(1:8), '-', ztime(1:2), ':', ztime(3:4), ':', ztime(5:6)
   write ( 6, *) trim(tag)
   write(6,*) ' This version uses quanta from HITRAN to attribute extra parameters to each transition record.'

   print *, ' hbin v0.9.5.0'


   ! --- read in band, isotope info from sfit4.ctl file fr this fit
   call read_ctrl

   ! --- read in paths to HITRAN files
   call read_input( hasc, wave5(1), wave6(nband), HFL, GLP, LMX, SDV )

   ! --- see if we need to separate out isotopes
   !print *, useiso
   if( useiso )call rdisofile( 6 )

   ! --- define & open binary (.hbin) & ascii (.hasc) output files
   write(hbfile,101) wave5(1), wave6(nband)
   do i=1, 30
      if( iachar(hbfile(i:i)) .eq. 32 ) hbfile(i:i) = '0'
   enddo
   hafile = hbfile
   write(hafile(28:30),105) 'asc'
   write(6,102) 'Saving binary line data to : ', hbfile
   if( hasc )write(6,113) 'Saving ascii line data to  : ', hafile

   open( hblun, file=hbfile, status='unknown', iostat=istat, form='unformatted' )
   if( hasc )open( halun, file=hafile, status='unknown', iostat=istat )

   write(hblun) nband
   if( hasc )write(halun,*) nband
   do i=1, nband
      write(hblun) i, wave5(i), wave6(i)
      if( hasc )write(halun,*) i, wave5(i), wave6(i)
   enddo

   ! --- initial fill of HITRAN + CIA line parameters with line from each file
   write(6,102) 'HITRAN file for input:'
   write(6,113) 'index molid  lun  wavnum    name'
   do ifl = 1, hnml

      ! --- get file info
      inquire( hfl(ifl)%lun, position=pos, opened=oped, iostat=iost, name=nam )
      if( iost .ne. 0 )then
         print *,''
         print *, 'HITRAN file error : ', trim(nam)
         print *, hfl(ifl)%lun, pos, oped, iost, trim(nam)
         stop 'hbin error'
      endif

      ! --- read first line info from buffer
      !read( hbuf(ifl), 107 ) hlp(ifl)%mo, hlp(ifl)%is, hlp(ifl)%nu, hlp(ifl)%sl
      call filh( hlp(ifl), hfl(ifl) )

      ! --- print initial line
      write(6, 111) ifl, hlp(ifl)%mo, hfl(ifl)%lun, hlp(ifl)%nu, trim(nam)

   enddo

   ! --- fill Galatry line parameters struct with all line data from each file

   do ifl = 1, gnml

      ! --- first buf already read
      ! glp(ifl)%n = 1
      ! ! HITRAN 2012
      ! print*, glp(ifl)%buf
      ! read( glp(ifl)%buf, 108 ) glp(ifl)%mo(1), glp(ifl)%is(1), &
      !      glp(ifl)%v1(1), glp(ifl)%v2(1), glp(ifl)%branch(1), glp(ifl)%j(1), &
      !      glp(ifl)%g0_air(1), glp(ifl)%beta(1)
      ! HITRAN 2008
      ! read( glp(ifl)%buf, 108 ) glp(ifl)%mo(1), glp(ifl)%is(1), glp(ifl)%nu(1), glp(ifl)%bt(1)
      !print *, ifl, 1, glp(ifl)%mo(1), glp(ifl)%is(1), glp(ifl)%nu(1), glp(ifl)%bt(1)

      ! --- loop over files
      glp(ifl)%n = 1
      do i = 1, nglines
         read( glp(ifl)%lun, 109, end=5 ) buf
         !print *, buf
         ! HITRAN 2012
         ind = glp(ifl)%n
         glp(ifl)%beta(ind) = 0.0D0
         read( buf, 108 ) glp(ifl)%mo(ind), glp(ifl)%is(ind), glp(ifl)%qa(ind), glp(ifl)%g0_air(ind), glp(ifl)%beta(ind)
         if (glp(ifl)%beta(ind) .le. tiny(0.0D0)) cycle
         !         print *, glp(ifl)%v1(ind), glp(ifl)%v2(ind), glp(ifl)%branch(ind), glp(ifl)%j(ind)
         glp(ifl)%n = glp(ifl)%n + 1

      ! ad hoc Galatry files from AG 2008
      !         read( glp(ifl)%buf,108  ) glp(ifl)%mo(i), glp(ifl)%is(i), glp(ifl)%nu(i), glp(ifl)%bt(i)
         !print *, ifl, i, glp(ifl)%mo(i), glp(ifl)%is(i), glp(ifl)%nu(i), glp(ifl)%bt(i)
         goto 6
 5       close( glp(ifl)%lun )
         exit
 6       continue

      enddo

      write(6,117) glp(ifl)%n, ' lines read in Galatry file : ', ifl

   enddo

   ! --- fill of Line Mixing data parameters from each file complete into structures
   if( lnml.gt.0 )then
      write(6,102) 'Line Mixing file for input:'
   endif
   do ifl = 1, lnml

      ! --- get file info
      inquire( lmx(ifl)%lun, position=pos, opened=oped, iostat=iost, name=nam )
      if( iost .ne. 0 )then
         print *,''
         print *, 'Line Mix file error : ', trim(nam)
         print *, lmx(ifl)%lun, pos, oped, iost, trim(nam)
         stop 'hbin error'
      endif

      lmx(ifl)%n = 1
      do i = 1, nglines
         read( lmx(ifl)%lun, 109, end=15 ) buf
         if (len_trim(buf).eq.0) cycle
         ! Read in auxiliary file in HITRAN 2012 format
         ! valid for CO2 line mixing file with F. Hase data
         ind = lmx(ifl)%n
         lmx(ifl)%g2_air(ind) = 0.0D0
         read( buf, 119 ) lmx(ifl)%mo(ind), lmx(ifl)%is(ind), lmx(ifl)%qa(ind)
         read (buf(63:), *) lmx(ifl)%lm_t1(ind), lmx(ifl)%lm_t2(ind), lmx(ifl)%lm_air(ind)

         lmx(ifl)%n = lmx(ifl)%n + 1

         goto 16

15       close( lmx(ifl)%lun )
         exit
16       continue

      end do
      write(6,117) lmx(ifl)%n, ' lines read in LM file : ', ifl

   enddo

   ! --- fill of Speed - Dependent Voigt data parameters from each file complete into structures
   if( snml.gt.0 )then
      write(6,102) 'Speed-Dependent Voigt file for input:'
   endif
   do ifl = 1, snml

      ! --- get file info
      inquire( sdv(ifl)%lun, position=pos, opened=oped, iostat=iost, name=nam )
      if( iost .ne. 0 )then
         print *,''
         print *, 'S-D Voigt file error : ', trim(nam)
         print *, sdv(ifl)%lun, pos, oped, iost, trim(nam)
         stop 'hbin error'
      endif

      sdv(ifl)%n = 1
      do i = 1, nglines
         read( sdv(ifl)%lun, 109, end=25 ) buf
         if (len_trim(buf).eq.0) cycle
         ! Read in auxiliary file in HITRAN 2012 format
         ind = sdv(ifl)%n
         sdv(ifl)%g2_air(ind) = 0.0D0
         read( buf, 119 ) sdv(ifl)%mo(ind), sdv(ifl)%is(ind), &
              sdv(ifl)%qa(ind)
         read (buf(64:), '(f5.4,f5.3,f3.2,f8.6)') sdv(ifl)%g0_air(ind), sdv(ifl)%g0_self(ind),&
              sdv(ifl)%td_g0_air(ind),sdv(ifl)%s_air(ind)
         read (buf, '(tr84,a7,a8,a8,a9,a8,a8)') (sdv_params(j), j=1,6)
         if (len_trim(sdv_params(1)).gt.0) read(sdv_params(1), *) sdv(ifl)%g2_air(ind)
         if (len_trim(sdv_params(3)).gt.0) then
            read(sdv_params(3), *) sdv(ifl)%lm_air(ind)
         end if

         if (sdv(ifl)%g2_air(ind) .le. tiny(0.0D0)) cycle

         sdv(ifl)%n = sdv(ifl)%n + 1

         goto 26
25       close( sdv(ifl)%lun )
         exit
26       continue
      end do
      write(6,117) sdv(ifl)%n, ' lines read in SDV file : ', ifl
   enddo



   nl = 0
   ! --- loop over bands
   do iband = 1, nband

      if( wave3(iband) .ge. wave4(iband) )stop 'sfit4.ctrl: bandpass limits out of order'

      write(6,100) iband, wave5(iband), wave3(iband), wave4(iband), wave6(iband)

      wstr   = wave5(iband)
      wavnum = wstr
      wstp   = wave6(iband)

      do while(( wavnum .le. wstp ) ) !.and. ( wavnum .ge. wstr ))

         ! --- find lowest wavenumber
         !print *, hfl(1:hnml)%lun
         inxt = minloc( hlp(1:hnml)%nu )
         ldx = inxt(1)
         !print *, 'low  ', hfl(ldx)%buf(1:30), wavnum, ldx, hfl(ldx)%lun

         ! --- if in band : write out that data to the hbin file
         if( wavnum .ge. wstr )then


            ! --- check if a Galatry beta can be appended
            do ifl = 1, gnml
               do i = 1, glp(ifl)%n
                  if ( glp(ifl)%mo(i) .eq. hlp(ldx)%mo .and. &
                       glp(ifl)%is(i) .eq. hlp(ldx)%is ) then
                     if (qu_equal(hlp(ldx)%qa, glp(ifl)%qa(i))) then
                        write( hfl(ldx)%buf(161:172), 110 ) glp(ifl)%beta(i)
                        hlp(ldx)%bt = real(glp(ifl)%beta(i),4)
                        hlp(ldx)%flag(GALATRY_FLAG) = .TRUE.
                        dum = flagoff + GALATRY_FLAG
                        write( hfl(ldx)%buf(dum:dum), '(l1)' ) .TRUE.
                     endif
                  end if
               enddo
            enddo

            ! --- check if this line has linemixing data to be attached
            do ifl = 1, lnml
               do i=1, lmx(ifl)%n
                  if ( lmx(ifl)%mo(i) .eq. hlp(ldx)%mo .and. &
                       lmx(ifl)%is(i) .eq. hlp(ldx)%is ) then
                     if (qu_equal(hlp(ldx)%qa, lmx(ifl)%qa(i))) then
                        hlp(ldx)%ylm = lmx(ifl)%lm_air(i)
                        hlp(ldx)%lmtk1 = lmx(ifl)%lm_t1(i)
                        hlp(ldx)%lmtk2 = lmx(ifl)%lm_t2(i)
                        write( hfl(ldx)%buf(220:280), 1121 ) hlp(ldx)%lmtk1, hlp(ldx)%lmtk2, hlp(ldx)%ylm
                        hlp(ldx)%flag(LM_FLAG) = .TRUE.
                        dum = flagoff + LM_FLAG
                        write( hfl(ldx)%buf(dum:dum), '(l1)' ) .TRUE.
                        exit
                     endif ! right wavenumber
                  endif ! right molecule
               enddo ! i, records in file
            enddo ! line mix files


            ! --- check if this line has SDV data to be attached
            do ifl = 1, snml
               do i=1, sdv(ifl)%n
                  if ( sdv(ifl)%mo(i) .eq. hlp(ldx)%mo .and. &
                       sdv(ifl)%is(i) .eq. hlp(ldx)%is ) then
                     if (qu_equal(hlp(ldx)%qa, sdv(ifl)%qa(i))) then
                        hlp(ldx)%gamma0  = real(sdv(ifl)%g0_air(i))          ! gam0 for SDV
                        hlp(ldx)%gamma2  = real(sdv(ifl)%g2_air(i))          ! gam2 for SDV
                        hlp(ldx)%shift0  = real(sdv(ifl)%s_air(i))            ! shift0 for SDV
                        hlp(ldx)%shift2  = real(sdv(ifl)%ts_air(i))            ! td shift0 for SDV
                        hlp(ldx)%ylm = 0.0
                        hlp(ldx)%lmtk1 = 0.0
                        hlp(ldx)%lmtk2 = 0.0
                        if (sdv(ifl)%lm_air(i).gt.tiny(0.0d0)) then
                           hlp(ldx)%ylm   = real(sdv(ifl)%lm_air(i))
                           hlp(ldx)%lmtk1   = 1.0d0
                           hlp(ldx)%lmtk2   = 0.0d0
                           hlp(ldx)%flag(LM_FLAG) = .TRUE.
                           dum = flagoff + LM_FLAG
                           write( hfl(ldx)%buf(dum:dum), '(l1)' ) .TRUE.
                        end if
                        write( hfl(ldx)%buf(172:280), 112 ) hlp(ldx)%gamma0, hlp(ldx)%gamma2, &
                             hlp(ldx)%shift0, hlp(ldx)%shift2, hlp(ldx)%lmtk1, hlp(ldx)%lmtk2, hlp(ldx)%ylm
                        hlp(ldx)%flag(SDV_FLAG) = .TRUE.
                        dum = flagoff + SDV_FLAG
                        write( hfl(ldx)%buf(dum:dum), '(l1)' ) .TRUE.
                        exit
                     endif ! right wavenumber
                  endif ! right molecule
               enddo ! i, records in file
            enddo ! line mix files


            ! --- check if this is an isotope that is to be separated out
            ! --- because of identification of extra line parameters via mo id and quantum numbers,
            !     renumbering of the isotopes has to come last.
            do i=1, nisosep
               !print *, i, nisosep, hlp(lun)%mo, hlp(lun)%is, oldid(i), oldiso(i), newid(i), newiso(i)
               if((hlp(ldx)%mo .eq. oldid(i)) .and. (hlp(ldx)%is .eq. oldiso(i)))then
                  hlp(ldx)%sl = hlp(ldx)%sl / isoscale(i)
                  !write(6,*) i, hlp(ldx)%mo, hlp(ldx)%is, newid(i), newiso(i)
                  hlp(ldx)%mo = newid(i)
                  hlp(ldx)%is = newiso(i)
                  write( hfl(ldx)%buf(1:25), 107 ) newid(i), newiso(i), hlp(ldx)%nu, hlp(ldx)%sl
                  exit
               endif
            enddo

            ! --- save this line
            nl = nl +1
            !print *, 'writing ', hbuf(ldx)(1:30), ldx, hlun(ldx)
            if( hasc )write(halun,'(a)') hfl(ldx)%buf
            write(hblun) hlp(ldx)

         endif

         ! --- re-fill hbuf/hlp with next line from that file
         !print *, ldx, hfl(ldx)%lun
         read( hfl(ldx)%lun, 106, end=10 ) hfl(ldx)%buf
!         read( hlp(ldx)%buf, 107 ) hlp(ldx)%mo, hlp(ldx)%is, hlp(ldx)%nu, hlp(ldx)%sl
         call filh( hlp(ldx), hfl(ldx) )

         !print *,''
         !print *, 'read ', hbuf(ldx)(1:30), ldx, hlun(ldx)

         if( hlp(ldx)%nu .le. wave6(nband) )goto 11

         ! --- end of a file or last read line is too high a wavenumber
      10 continue
         inquire( hfl(ldx)%lun, name=nam )
         write(6,120) '  closing file : ', hfl(ldx)%lun, trim(nam)

         close( hfl(ldx)%lun )
         do i = ldx, hnml-1
            hlp(i)    = hlp(i+1)
            hfl(i)    = hfl(i+1)
         enddo
         hnml = hnml -1
         if( hnml .eq. 0 )goto 20

      11 continue
         wavnum = minval( hlp(1:hnml)%nu )
         !print *, wavnum, wstp, hnml

      enddo ! while

   enddo ! nband

   ! --- closed last hitran files - no more lines
   20 continue

   do ifl=1, hnml
      print *, 'closing file unit : ', hfl(ifl)%lun
      close( hfl(ifl)%lun )
   enddo

   close( halun )
   if( hasc )close( hblun )

   write(6,103) 'Lines saved to output hbin file : ', nl

!   deallocate( hfl, hlp, glp )

stop

100 format( /,'Band: ',i5,4f14.5,f12.3,f12.6 )
101 format( f12.6, '-', f12.6, '.hbin' )
102 format( /, a, a )
103 format( /,a, i10 )
105 format( a3 )
106 format( a200 )
107 format( i2,i1,f12.6,1p,e10.3,10x,0p,f5.4,f5.4,f10.4,f4.2,f8.6,f7.4)
108 format( i2,i1,a60,f7.4,f7.4 )
!HITRAN 2008 108 format( i2,i1,f12.6,f10.0)
109 format( a255 )
110 format( f12.5 )
111 format( 3i5,f12.5,2x,a)
112 format( 8e12.4 )
1121 format( 4e12.4 )
113 format( a, a )
!115 format(a, 2i4, 2(f14.6, 2i4))
!116 format( a,4i4,f8.4,2f14.6)
117 format(/, i5, a, i5)
!118 format( 3i5,i10,f12.5,2x,a)
119 format( i2,i1,a60 )
120 format( a,i4,3x,a )

end program hbin

! --- fill a hitran record from its buffer
subroutine filh( hd, hf )

   use hitran

   type (hitrandata), intent(out)    :: hd
   type (hitranfile), intent(inout)  :: hf

   hd%flag(1:8) = .FALSE.

! --- fill defaults for beta, gamma, eta, lm
   hd%bt     = 0.0            ! beta for galatry
   hd%gamma0 = 0.0            ! gamma 0 for sdv
   hd%gamma2 = 0.0            ! gamma 2 for sdv
   hd%shift0   = 0.0          ! shift 0 for sdv
   hd%shift2   = 0.0          ! shift 2 for sdv
   hd%lmtk1  = 0.0            ! lmtk1 for line mixing
   hd%lmtk2  = 0.0            ! lmtk2 for line mixing
   hd%ylm    = 0.0            ! ylm for line mixing
   write( hf%buf(161:244), 110 ) 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

   ! flag only 0, 1, 2 so far - maybe don't need these anyway...
   select case ( hf%flag )
   case (0)       ! HITRAN line

      ! --- read parameters
      read( hf%buf, 107) hd%mo, hd%is, hd%nu, hd%sl, hd%ea, hd%ah, hd%sh, hd%el, hd%tx, hd%ps, &
                         hd%qa, hd%er, hd%lm, hd%uw, hd%lw, hd%bt

      ! --- map hitran molecule id to sfit id
      hd%mo = hf%mo

   case (1)        ! O2CIA, id 49/0,1 -> 1,2

      ! --- read parameters
      read( hf%buf, 107) hd%mo, hd%is, hd%nu, hd%sl, hd%ea, hd%ah, hd%sh, hd%el, hd%tx, hd%ps

      ! --- map cia molecule id to sfit id
      ! --- in cia files ids for F & C are both 70
      ! --- make first file iso 1 and second file iso 2 - these are the flag ids
      ! --- files in alphabetical order in hbin.input file
      hd%mo = 49

      ! --- map cia molecule iso to sfit iso file iso's are 0,1 (f,s) change to 1,2
      hd%is = hd%is + 1

      if(hd%is .eq. 1) hd%flag(fcia_flag) = .TRUE.      ! fcia o2
      if(hd%is .eq. 2) hd%flag(scia_flag) = .TRUE.      ! scia o2

   case default
   print*, hf%flag
   print*, hf%buf
      stop ' hitran flag out of bounds'

   end select

   write( hf%buf(1:2), '(i2)' ) hd%mo
   write( hf%buf(3:3), '(i1)' ) hd%is
   write( hf%buf(flagoff+1:flagoff+8), '(8l1)' ) hd%flag(1:8)

   return
! 1-160 hitan
! 161 - 172 galatry beta
! 173 - 184 sdv gam0
! 185 - 196 sdv gam2
! 197 - 208 sdv shift0
! 209 - 220 lmx ltk1
! 221 - 232 lmx ltk2
! 233 - 244 lmx ylm
! 245 - 256 ---
! 257 - 268 ---
! 269 - 280 ---
! 281 - 288 lmx flags

107 format( i2, i1, f12.6, 1p, e10.3, e10.0, 0p, 2(f5.4), f10.4, f4.2, f8.6, &
            a60, a18, a1, 2f7.0, f10.0 )
110 format( 7f12.5 )



end subroutine filh


subroutine read_input( hasc, wstr, wstp, HFL, GLP, LFL, SDV )

   use hitran

   implicit none

   logical, intent(inout)  :: hasc
   real(double)            :: wstr, wstp, wavnum
   integer                 :: lun, mo, iso
   character (len=10)      :: ifilename = 'hbin.input'
   integer                 :: j, i, n, istat, ilun=9
   logical                 :: fexist
   character (len=160)     :: buffer, linebuffer, path, filename

   TYPE (GALATRYDATA), intent(inout)   :: GLP(ngal)
   TYPE (HITRANFILE),  intent(inout)   :: HFL(nhit+ncia)
   TYPE (GALATRYDATA), intent(inout)   :: LFL(nlmx)
   TYPE (GALATRYDATA), intent(inout)   :: SDV(nsdv)

   ! --- open hbin.input file if its here
   inquire( file=trim(ifilename), exist = fexist)
   if( .not. fexist )then
      write(*,*) 'File "', trim(ifilename), '" does not exist.'
      stop
   endif
   open( ilun, file=ifilename, status='old', iostat=istat )

   ! --- read in ascii output flag
   call nextbuf( ilun, buffer )
   read(buffer,'(l10)') hasc
   !print*, hasc

   ! --- read path to hitran files
   call nextbuf( ilun, buffer )
   path = trim( buffer )
   write(6,112) 'Linelist path : ', path

   ! --- read number of expected hitran files (max=99)
   call nextbuf( ilun, buffer )
   read(buffer,*) nfiles
   write(6,111) 'Number of HITRAN files to search : ', nfiles
   if( nfiles .gt. nhit )stop 'too many hitran files'

   ! --- loop over files and find valid files with data
   !     position file pointer at starting wavenumber
   !     save its lun in hlp(;)%un 1 to hnml

   ! --- read paths to (up to) 99 hitran formatted files
   ! --- accommodates partial hitran files like the cia pseudo lines since know mol id
   hnml = 0
   do i = 1, nfiles

      ! --- find filename
      call nextbuf( ilun, linebuffer )
      filename = trim(path) // trim(linebuffer)
      n = len_trim(filename)
      ! catch eg 065_CH3CNPL/ 2007.sudo.ch3cn
      do j=1, n
        if( filename(j:j) .eq. ' ' )then
         !print*, i, j, filename(j-1:j-1)
         n = j-1
         exit
        endif
      enddo
      ! catch eg 065_CH3CNPL/
      if( filename(n:n) .eq. '/' )cycle

      !write(6,110) 'Found HITRAN line file : ', trim(filename)

      ! --- open hitran file if its here
      inquire( file=trim(filename), exist = fexist)
      if( .not. fexist )then
         write(*,*) 'HITRAN file "', trim(filename), '" does not exist.'
         stop
      endif

      lun = hnml +10
      open( lun, file=filename, status='old', iostat=istat )

      ! --- find starting wavenumber in file
      do
         read( lun, 100, end=10 ) buffer
         read( buffer, 107 ) mo, iso, wavnum
         if( wavnum .ge. wstr )exit
      enddo
      if( wavnum .lt. wstr )goto 10
      if( wavnum .gt. wstp )goto 10
      goto 11

      ! --- no lines in this region
   10 close( lun )
      goto 12

      ! --- save this line
   11 hnml           = hnml +1
      hfl(hnml)%buf  = buffer
      hfl(hnml)%lun  = lun
      read( linebuffer(1:3), '(i3)' ) hfl(hnml)%mo
      hfl(hnml)%flag = 0
      if( (hfl(hnml)%mo .eq. 49) .and. (mo .eq. 7) ) hfl(hnml)%flag = 1 ! o2cia
      write(6,*)''
      write(6,110) ' File : ', trim(filename)
      write(6,113) wstr, wstp, hnml, hfl(hnml)%lun, hfl(hnml)%mo, mo, wavnum

   12 continue

   enddo
   write(6,114) ' Number of HITRAN molecules/files found : ', hnml


   ! --- Galatry data files - block 2 in hbin.input
   ! --- read number of expected Galatry files (max=2)
   ! --- Galatry files are unique format from hitran
   call nextbuf( ilun, buffer )
   read(buffer,*) nfiles
   write(6,114) 'Number of Galatry files to search : ', nfiles
   if( nfiles .gt. ngal )stop 'too many hitran files'

   ! --- read paths to (up to) 99 hitran files
   stlun = hfl(hnml)%lun
   lun   = stlun
   gnml  = 0
   snml  = 0
   do i = 1, nfiles

      ! --- find the name of the next galatry input file
      call nextbuf( ilun, linebuffer )
      filename = trim(path) // trim(linebuffer)
      n = len_trim(filename)
      if( filename(n:n) .eq. '/' )cycle

      gnml = gnml + 1
      write(6,110) 'Found Galatry line file : ', trim(filename)

      ! --- open Galatry file if its here
      inquire( file=trim(filename), exist = fexist)
      if( .not. fexist )then
         write(*,*) 'Galatry file "', trim(filename), '" does not exist.'
         stop
      endif

      lun = stlun + gnml
      open( lun, file=filename, status='old', iostat=istat )

      ! --- find starting wavenumber in file
!      do
         read( lun, 100, end=20 ) buffer
!         read( buffer, 107 ) mo, iso, wavnum
!         if( wavnum .ge. wstr )exit
!      enddo
 !     if( wavnum .lt. wstr )goto 20
 !     if( wavnum .gt. wstp )goto 20
      goto 21

      ! --- no lines in this region
   20 close( lun )
      gnml   = gnml - 1
      goto 22

      ! --- save this line
   21 glp(gnml)%buf = buffer(1:64)
      glp(gnml)%lun = lun
      read( linebuffer(1:3), '(i3)' ) glp(gnml)%mo(1)
      write(6,110) ' File : ', trim(filename)
      write(6,113) wstr, wstp, gnml, glp(gnml)%lun, mo, glp(gnml)%mo(1)

   22 continue

   enddo
   write(6,114) '   Number of Galatry molecules found : ', gnml
   stlun = stlun + gnml

   ! --- Line mixing data files - block 3 in hbin.input
   ! --- read number of expected LM files (max=2)
   ! --- LM files are unique format from hitran, Galatry...
   ! --- save position in file and read as they are written (like hitran base files)
   ! only accounting for CO2 so far
   call nextbuf( ilun, buffer )
   read(buffer,*) nfiles
   write(6,114) 'Number of LM files to search : ', nfiles
   if( nfiles .gt. nlmx )stop 'too many LM files'

   ! --- read paths to (up to) 99 hitran files
   lnml = 0
   lun = stlun
   do i = 1, nfiles

      ! --- find the name of the next galatry input file
      call nextbuf( ilun, linebuffer )
      filename = trim(path) // trim(linebuffer)
      n = len_trim(filename)
      if( filename(n:n) .eq. '/' )cycle

      lnml = lnml + 1
      write(6,110) 'Found LineMix line file : ', trim(filename)

      ! --- open LineMix file if its here
      inquire( file=trim(filename), exist = fexist )
      if( .not. fexist )then
         write(*,*) 'LM file "', trim(filename), '" does not exist.'
         stop
      endif

      lun = lun + lnml
      !print*, lun
      open( lun, file=filename, status='old', iostat=istat )
      read( lun, 100, end=30 ) buffer
      goto 31
30    close( lun )
      lnml = lnml - 1
      goto 32

      ! --- save this line
31    lfl(lnml)%buf  = buffer(1:160)
      lfl(lnml)%lun  = lun
      read( linebuffer(1:3), '(i3)' ) lfl(lnml)%mo(1)
      write(6,*)''
      write(6,110) ' File : ', trim(filename)
      write(6,113) wstr, wstp, lnml, lfl(lnml)%lun, lfl(lnml)%mo(1)

32    continue

   enddo

   write(6,114) '   Number of LM molecules found : ', lnml
   stlun = stlun + lnml


   ! --- Speed Dependent Voigt data files - block 3 in hbin.input
   ! --- read number of expected SDV files (max=2)
   ! --- From HITRAN 2012 SDV files are simlar format than GALATRY
   ! --- !!! SDV files are unique format from hitran, Galatry, but same as Line Mixing lists...
   ! --- save position in file and read as they are written (like hitran base files)
   ! only accounting for CO2 so far
   call nextbuf( ilun, buffer )
   read(buffer,*) nfiles
   write(6,114) 'Number of SDV files to search : ', nfiles
   if( nfiles .gt. nsdv )stop 'too many SDV files'

   ! --- read paths to (up to) 99 hitran files
   snml = 0
   lun = stlun
   do i = 1, nfiles

      ! --- find the name of the next galatry input file
      call nextbuf( ilun, linebuffer )
      filename = trim(path) // trim(linebuffer)
      n = len_trim(filename)
      if( filename(n:n) .eq. '/' )cycle

      snml = snml + 1
      write(6,110) 'Found SDVoigt line file : ', trim(filename)

      ! --- open SDV file if its here
      inquire( file=trim(filename), exist = fexist )
      if( .not. fexist )then
         write(*,*) 'SDV file "', trim(filename), '" does not exist.'
         stop
      endif

      lun = lun +snml
      !print*, lun
      open( lun, file=filename, status='old', iostat=istat )

      ! --- find starting wavenumber in file
!      do
         read( lun, 100, end=40 ) buffer
      !    read( buffer, 107 ) mo, iso, wavnum
      !    if( wavnum .ge. wstr )exit
      ! enddo
      ! if( wavnum .lt. wstr )goto 40
      ! if( wavnum .gt. wstp )goto 40
      goto 41

      ! --- no lines in this region
   40 close( lun )
      snml = snml - 1
      goto 42

      ! --- save this line
41    sdv(snml)%buf  = buffer(1:160)
      sdv(snml)%lun  = lun
      read( linebuffer(1:3), '(i3)' ) sdv(snml)%mo(1)
      write(6,*)''
      write(6,110) ' File : ', trim(filename)
      write(6,113) wstr, wstp, snml, sdv(snml)%lun, sdv(snml)%mo(1)

   42 continue

   enddo
   write(6,114) '   Number of SDV molecules found : ', snml

   close(ilun)

!do i=1, hnml
!   inquire( hlun(i), position=pos, opened=oped, iostat=iost, name=nam )
!   print *,''
!   print *, i, hlun(i), pos, oped, iost, trim(nam)
!enddo


return

100 format( a160 )
107 format(i2,i1,f12.6,1p,e10.3,10x,0p,f5.4,f5.4,f10.4,f4.2,f8.6,f7.4)
110 format( a, a )
111 format( a, i10 )
112 format( /, a, a )
113 format( 2f14.6, 4i6, 2x, f14.6 )
114 format( /, a, i10 )

end subroutine read_input


! --- read in bands from sfit.ctl file
subroutine read_ctrl

   use hitran

   implicit none
   character (len=9) :: filename = 'sfit4.ctl'
   integer           :: istat, nr_keys, iband, nextra
   logical           :: fexist
   real(double)      :: dwave

   ! --- open sfit4.ctl file if its here
   inquire( file=filename, exist = fexist )
   if( .not. fexist ) then
      write(*,*) 'file ', trim(filename), ' does not exist'
      stop
   endif
   open( bp_nr, file=filename, status='old', iostat=istat)

   ! --- read in band parameters
   do
      call read_line_binput( keyword, nr_keys, value, istat )
      if(( istat .lt. 0 ) .and. ( nr_keys .eq. 0 ))exit
      if( nr_keys .eq. 0 )then
         cycle
      endif
      select case ( trim( adjustl( keyword( 1 ))))
         case ( 'band' )
            call read_band_section( keyword, value )
         case ('file')
            call read_file_section(keyword, value)
         case ('fw')
            call read_fw_section(keyword, value)
      end select
   end do
   write(6,101) 'Found number of bands : ', nband

   ! --- Loop over bands
   do iband = 1, nband
      ! --- 10 res units to account for shifts - from initialize.f90:setup
      dwave  = 10.d0/pmax(iband)
      nextra = nint( dwave/dn(iband))
      !  --- interval for input of line data, dlines accounts for out of band absorption
      wave5(iband) = wave3(iband) - nextra*dn(iband) - dlines - 0.5/pmax(iband)
      wave6(iband) = wave4(iband) + nextra*dn(iband) + dlines + 0.5/pmax(iband)
      write(6,100) iband, wave3(iband), wave5(iband), wave4(iband), wave6(iband), pmax(iband), dn(iband)
   enddo

   close( bp_nr )

return

100 format( i5,4f14.5,f12.3,f12.6 )
101 format( /, a, i10 )

end subroutine read_ctrl


!--------------------------------------------------------------------------------------------
subroutine nextbuf(ifile, buffer)

   implicit none

   integer,       intent(in)   :: ifile
   character(160),intent(out)  :: buffer

   buffer(1:1) = '#'
   do while( adjustl( buffer(1:1)) .eq. '#' )
      read(ifile,'(a160)', err=20, end=21) buffer
      !print *, buffer
      if( len_trim(buffer) .eq. 0 )buffer(1:1) = '#'
   end do

   return

 20 stop 'nextbuf : error reading input file'
 21 stop 'nextbuf : end of input file'

end subroutine nextbuf


!--------------------------------------------------------------------------------------------
logical function qu_equal(quanta1, quanta2)
  ! compares quanta1 and quanta2 field in HITRAN 2004 format, retruns T if they are equal
  ! until now, only string compare, may get more complicated though
  implicit none
  character (len=*), intent(in):: quanta1, quanta2
  !logical :: qu_equal

  qu_equal = .false.
  if (quanta1.eq.quanta2) qu_equal = .true.
  return

end function qu_equal

