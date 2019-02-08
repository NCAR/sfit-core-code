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


program hbin

! program to create a binary linelist file from initial HITRAN formatted files
!  for use by sfit4

   use hitran

   implicit none

   integer              :: ldx, nl, i, j, ifl
   integer              :: hblun=7, halun=8, istat, iband, inxt(1)
!   integer(long_log)    :: nulm, nuht
   real(double)         :: wavnum, wstr, wstp, tmpreal
   character (len=30)   :: hbfile, hafile
   character (len=200)  :: nam
   character (len=1)    :: pos
   character (len=255)  :: buf
   logical              :: oped, hasc
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
   type (galatrydata), dimension(ncorr)        :: elp

   call date_and_time (cdate, ztime, zone)
   write (tag,*) trim(version), ' runtime:', cdate(1:8), '-', ztime(1:2), ':', ztime(3:4), ':', ztime(5:6)
   write ( 6, *) trim(tag)
   write(6,*) ' This hbin version uses quanta from HITRAN to attribute extra parameters to each transition record.'
   !write(6,*) ' hbin v0.9.6.0'


   ! --- read in band, isotope info from sfit4.ctl file fr this fit
   call read_ctrl

   ! --- read in paths to HITRAN files
   call read_input( hasc, wave5(1), wave6(nband), HFL, GLP, LMX, SDV, ELP )

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
         lmx(ifl)%lm_air(ind) = 0.0d0
         lmx(ifl)%lm_t1(ind) = 0.0d0
         lmx(ifl)%lm_t2(ind) = 0.0d0
         ! test for linemixing parameter
         read (buf(64:), *, end=16) lmx(ifl)%lm_air(ind), lmx(ifl)%lm_t1(ind), lmx(ifl)%lm_t2(ind)

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

      sdv(ifl)%n = 0
      do i = 1, nglines
         read( sdv(ifl)%lun, 109, end=25 ) buf
         if (len_trim(buf).eq.0) cycle
         ! Read in auxiliary file in HITRAN 2012 format
         ind = sdv(ifl)%n
         sdv(ifl)%g2_air(ind) = 0.0D0
         read( buf, 119 ) sdv(ifl)%mo(ind), sdv(ifl)%is(ind), &
              sdv(ifl)%qa(ind)
         call flush()
         read (buf(64:), 701) sdv(ifl)%g0_air(ind), tmpreal, tmpreal, tmpreal, sdv(ifl)%g2_air(ind)
!         read (buf(63:), 701, err=26) sdv(ifl)%g0_air(ind), tmpreal, tmpreal, tmpreal, sdv(ifl)%g2_air(ind)
         if (sdv(ifl)%g2_air(ind).lt.tiny(0.0d0)) then
            goto 26
         end if
!              sdv(ifl)%g2_self(ind), sdv(ifl)%s0_self(ind)
!         if (tiny(0.0E0).ge.abs(0.0749D0 - sdv(ifl)%g0_air(ind))) then
!            print *, 'bla2 ', sdv(ifl)%g0_air(ind), sdv(ifl)%s0_self(ind), sdv(ifl)%qa(ind)
!            call flush()
!         end if

!         read (buf(86:), *, end=26) sdv(ifl)%g2_air(ind)

!         if ((sdv(ifl)%g2_air(ind) .le. tiny(0.0D0)).and.(sdv(ifl)%s0_self(ind) .le. tiny(0.0D0))) cycle

         sdv(ifl)%n = sdv(ifl)%n + 1

         goto 26
25       close( sdv(ifl)%lun )
         exit
26       continue
      end do
      write(6,117) sdv(ifl)%n, ' lines read in SDV file : ', ifl
   enddo

!701 format(f6.5,f4.3,f4.2,f8.6,1x,g15.6,1x,g15.6,1x,g15.6)
701 format(f6.5,f5.4,f3.2,f8.6,1x,g15.6)

   ! --- fill CORRELATION line parameters struct with all line data from each file
   do ifl = 1, enml

      ! --- loop over files
      elp(ifl)%n = 1
      do i = 1, nglines
         read( elp(ifl)%lun, 109, end=35 ) buf
         !print *, buf
         ! HITRAN 2012
         ind = elp(ifl)%n
         elp(ifl)%eta(ind) = 0.0D0
         read( buf, 108 ) elp(ifl)%mo(ind), elp(ifl)%is(ind), elp(ifl)%qa(ind), elp(ifl)%g0_air(ind), elp(ifl)%eta(ind)
         if (elp(ifl)%eta(ind) .le. tiny(0.0D0)) cycle
         elp(ifl)%n = elp(ifl)%n + 1

         goto 36
35       close( elp(ifl)%lun )
         exit
36       continue

      enddo
      write(6,117) elp(ifl)%n, ' lines read in CORR file : ', ifl
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

!print*, 'gal'
            ! --- check if a Galatry beta can be appended
            do ifl = 1, gnml
               do i = 1, glp(ifl)%n
                  if ( glp(ifl)%mo(i) .eq. hlp(ldx)%mo .and. &
                       glp(ifl)%is(i) .eq. hlp(ldx)%is ) then
!                       print*, 'gal', hlp(ldx)%qa, glp(ifl)%qa(i)
                     if (qu_equal(hlp(ldx)%qa, glp(ifl)%qa(i))) then
                        write( hfl(ldx)%buf(161:172), 110 ) glp(ifl)%beta(i)
                        hlp(ldx)%bt = real(glp(ifl)%beta(i),4)
                        hlp(ldx)%flag(GALATRY_FLAG) = .TRUE.
                        dum = flagoff + GALATRY_FLAG
                        print*, 'gal',  hlp(ldx)%mo, hlp(ldx)%is, hlp(ldx)%flag(GALATRY_FLAG), hlp(ldx)%qa, glp(ifl)%qa(i)
                        write( hfl(ldx)%buf(dum:dum), '(l1)' ) .TRUE.
                     endif
                  end if
               enddo
            enddo

!print*, 'lnm'
            ! --- check if this line has linemixing data to be attached
            do ifl = 1, lnml
               do i=1, lmx(ifl)%n
                  if ( lmx(ifl)%mo(i) .eq. hlp(ldx)%mo .and. &
                       lmx(ifl)%is(i) .eq. hlp(ldx)%is ) then
                     if (qu_equal(hlp(ldx)%qa, lmx(ifl)%qa(i))) then
                        hlp(ldx)%ylm = lmx(ifl)%lm_air(i)
                        hlp(ldx)%lmtk1 = lmx(ifl)%lm_t1(i)
                        hlp(ldx)%lmtk2 = lmx(ifl)%lm_t2(i)
                        write( hfl(ldx)%buf(209:244), 1121 ) hlp(ldx)%lmtk1, hlp(ldx)%lmtk2, hlp(ldx)%ylm
                        hlp(ldx)%flag(LM_FLAG) = .TRUE.
                        dum = flagoff + LM_FLAG
                        print*, 'lnm',  hlp(ldx)%mo, hlp(ldx)%is, hlp(ldx)%flag(LM_FLAG), hlp(ldx)%qa, lmx(ifl)%qa(i)
                        write( hfl(ldx)%buf(dum:dum), '(l1)' ) .TRUE.
                        exit
                     endif ! right wavenumber
                  endif ! right molecule
               enddo ! i, records in file
            enddo ! line mix files

!print*, 'sdv'
            ! --- check if this line has SDV data to be attached
            do ifl = 1, snml
               do i=1, sdv(ifl)%n
                  if ( sdv(ifl)%mo(i) .eq. hlp(ldx)%mo .and. &
                       sdv(ifl)%is(i) .eq. hlp(ldx)%is ) then
                     if (qu_equal(hlp(ldx)%qa, sdv(ifl)%qa(i))) then
                        hlp(ldx)%gamma0  = real(sdv(ifl)%g0_air(i))          ! gam0 for SDV
                        hlp(ldx)%gamma2  = real(sdv(ifl)%g2_air(i))          ! gam2 for SDV
                        hlp(ldx)%shift0  = real(sdv(ifl)%s0_air(i))          ! shift0 for SDV
                        hlp(ldx)%shift2  = real(sdv(ifl)%ts_air(i))          ! shift2 for SDV
                        write( hfl(ldx)%buf(172:209), 112 ) hlp(ldx)%gamma0, hlp(ldx)%gamma2, &
                             hlp(ldx)%shift0
                        hlp(ldx)%flag(SDV_FLAG) = .TRUE.
                        dum = flagoff + SDV_FLAG
                        print*, 'sdv', hlp(ldx)%mo, hlp(ldx)%is, hlp(ldx)%flag(SDV_FLAG), hlp(ldx)%qa, sdv(ifl)%qa(i)
                        write( hfl(ldx)%buf(dum:dum), '(l1)' ) .TRUE.
                        exit
                     endif ! right wavenumber
                  endif ! right molecule
               enddo ! i, records in file
            enddo ! line mix files

            ! --- check if a CORRELATION can be appended
            ! do ifl = 1, enml
            !    do i = 1, elp(ifl)%n
            !       if ( elp(ifl)%mo(i) .eq. hlp(ldx)%mo .and. &
            !            elp(ifl)%is(i) .eq. hlp(ldx)%is ) then
            !          if (qu_equal(hlp(ldx)%qa, elp(ifl)%qa(i))) then
            !             write( hfl(ldx)%buf(220:232), 110 ) elp(ifl)%eta(i)
            !             hlp(ldx)%eta = real(elp(ifl)%eta(i))
            !             hlp(ldx)%flag(CORR_FLAG) = .TRUE.
            !             dum = flagoff + CORR_FLAG
            !             write( hfl(ldx)%buf(dum:dum), '(l1)' ) .TRUE.
            !          endif
            !       end if
            !    enddo
            ! enddo


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
112 format( 3e12.4 )
1121 format( 3e12.4 )
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
   character aiso*1           ! HITRAN 2016 CO2 has isotopes 11 & 12 labeled A and B, we will leave alone
                              ! in the ascii record but changed to 11, 12 in the binary (internal) record

   hd%flag(1:8) = .FALSE.

! --- fill defaults for beta, gamma, eta, lm
   hd%bt     = 0.0            ! beta for galatry
   hd%gamma0 = 0.0            ! gamma 0 for sdv
   hd%gamma2 = 0.0            ! gamma 2 for sdv
   hd%shift0   = 0.0          ! shift 0 for sdv
   hd%shift2   = 0.0          ! shift 2 for sdv
!   hd%eta      = 0.0
   hd%lmtk1  = 0.0            ! lmtk1 for line mixing
   hd%lmtk2  = 0.0            ! lmtk2 for line mixing
   hd%ylm    = 0.0            ! ylm for line mixing
   write( hf%buf(161:244), 110 ) 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

   ! flag only 0, 1, 2 so far - maybe don't need these anyway...
   select case ( hf%flag )
   case (0)       ! HITRAN line

      read( hf%buf, 107) hd%mo
      if( hd%mo .EQ. 2 ) then
         ! --- read parameters
         read( hf%buf, 108) hd%mo, aiso,  hd%nu, hd%sl, hd%ea, hd%ah, hd%sh, hd%el, hd%tx, hd%ps, &
                            hd%qa, hd%er, hd%lm, hd%uw, hd%lw, hd%bt
         if( aiso .EQ. 'A' )then
            hd%is = 11
         else if( aiso .EQ. 'B' ) then
            hd%is = 12
         else if( aiso .EQ. '0' ) then
            hd%is = 10
         else
            read(aiso, '(i1)') hd%is
         endif
      else
         ! --- read parameters
         read( hf%buf, 107) hd%mo, hd%is, hd%nu, hd%sl, hd%ea, hd%ah, hd%sh, hd%el, hd%tx, hd%ps, &
                            hd%qa, hd%er, hd%lm, hd%uw, hd%lw, hd%bt
      endif

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
      write( hf%buf(3:3), '(i1)' ) hd%is

      if(hd%is .eq. 1) hd%flag(fcia_flag) = .TRUE.      ! fcia o2
      if(hd%is .eq. 2) hd%flag(scia_flag) = .TRUE.      ! scia o2


   case (2)        ! N2CIA, id 52/0,1 -> 1,2

      ! --- read parameters
      read( hf%buf, 107) hd%mo, hd%is, hd%nu, hd%sl, hd%ea, hd%ah, hd%sh, hd%el, hd%tx, hd%ps

      ! --- map cia molecule id to sfit id
      ! --- in cia files ids for F & C are both 410
      ! --- make first file iso 1 and second file iso 2 - these are the flag ids
      ! --- files in alphabetical order in hbin.input file
      hd%mo = 52

      ! --- map cia molecule iso to sfit iso file iso's are 0,1 (f,s) change to 1,2
      hd%is = hd%is + 1
      write( hf%buf(3:3), '(i1)' ) hd%is

      if(hd%is .eq. 1) hd%flag(fcia_flag) = .TRUE.      ! fcia n2
      if(hd%is .eq. 2) hd%flag(scia_flag) = .TRUE.      ! scia n2

   case default
   print*, hf%flag
   print*, hf%buf
      stop ' hitran flag out of bounds'

   end select

   write( hf%buf(1:2), '(i2)' ) hd%mo
   write( hf%buf(flagoff+1:flagoff+8), '(8l1)' ) hd%flag(1:8)

   return
! 1-160 hitran
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
108 format( i2, a1, f12.6, 1p, e10.3, e10.0, 0p, 2(f5.4), f10.4, f4.2, f8.6, &
            a60, a18, a1, 2f7.0, f10.0 )
110 format( 7f12.5 )



end subroutine filh


subroutine read_input( hasc, wstr, wstp, HFL, GLP, LFL, SDV, ELP )

   use hitran
   use binput_4_0
   use binput_parse_4_0

   implicit none

   logical, intent(inout)  :: hasc
   real(double)            :: wstr, wstp, wavnum
   integer                 :: lun, mo, iso
   character (len=10)      :: ifilename = 'hbin.input'
   integer                 :: j, i, n, istat, ilun=9
   logical                 :: fexist
   character (len=160)     :: buffer, linebuffer, path, filename

   integer :: ctl_version = 2 ! 1 - original hbin.input version (till v0.9.4.4)
                              !     test for existence of ASC flag in the first valid line
                              ! 2 - tagged hbin.input version

   TYPE (GALATRYDATA), intent(inout)   :: GLP(ngal)
   TYPE (HITRANFILE),  intent(inout)   :: HFL(nhit+ncia)
   TYPE (GALATRYDATA), intent(inout)   :: LFL(nlmx)
   TYPE (GALATRYDATA), intent(inout)   :: SDV(nsdv)
   TYPE (GALATRYDATA), intent(inout)   :: ELP(ncorr)

   ! --- open hbin.input file if its here
   inquire( file=trim(ifilename), exist = fexist)
   if( .not. fexist )then
      write(*,*) 'File "', trim(ifilename), '" does not exist.'
      stop
   endif

   if(ctl_version.eq.2) then
      call read_hbin(ifilename, istat)
      !print *, 'ISTAT', istat
      if (istat.lt.0) goto 5
      goto 6
   end if

5  continue
   print *, 'Not a valid tagged hbin input, assume old input file version'
   ctl_version = 1

6  continue

   ! --- read in ascii output flag
   if (ctl_version.eq.1) then
      open( ilun, file=ifilename, status='old', iostat=istat )
      call nextbuf( ilun, buffer )
      ! if first line fails, assume it is tagged input
      read(buffer,'(l10)') out_ascii
   end if

   !print*, hasc
   hasc = out_ascii

   ! --- read path to hitran files
   if (ctl_version.eq.1) then
      call nextbuf( ilun, buffer )
      linelist_path = trim( buffer )
   end if
   linelist_path = trim( adjustl(linelist_path) )
   write(6,112) 'Path to Linelist : ', trim(linelist_path)

   ! --- read number of expected hitran files (max=99)
   if (ctl_version.eq.1) then
      call nextbuf( ilun, buffer )
      read(buffer,*) nhit_files
   end if
   write(6,111) 'Number of HITRAN files to search : ', nhit_files
   if( nhit_files .gt. nhit )stop 'too many hitran files'

   ! --- loop over files and find valid files with data
   !     position file pointer at starting wavenumber
   !     save its lun in hlp(;)%un 1 to hnml

   ! --- readlinelist_paths to (up to) 99 hitran formatted files
   ! --- accommodates partial hitran files like the cia pseudo lines since know mol id
   hnml = 0
   do i = 1, nhit_files

      ! --- find filename
      if (ctl_version.eq.1) then
         call nextbuf( ilun, linebuffer )
         filename = trim(linelist_path) // trim(linebuffer)
      else
         filename = adjustl(trim(linelist_path) // trim(hitran_files(i)))
         linebuffer = trim(hitran_files(i))
      end if
      n = len_trim(filename)
      !print*, n
      !print*, linebuffer
      !print*, filename
      !print*, nhit_files, i, hitran_files(i)

      ! catch eg 065_CH3CNPL/ 2007.sudo.ch3cn
      do j=1, n
        if( filename(j:j) .eq. ' ' )then
         print*, i, j, filename(j-1:j-1)
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
         read( buffer, 108 ) mo, wavnum
         if( wavnum .ge. wstr )exit
      enddo
      if( wavnum .lt. wstr )goto 10
      if( wavnum .gt. wstp )goto 10
      goto 11

      ! --- no lines in this region
   10 close( lun )
      goto 12

      ! --- save this line
   11 hnml           = hnml +1   ! actual # of hitran files
      hfl(hnml)%buf  = buffer
      hfl(hnml)%lun  = lun
      read( linebuffer(1:3), '(i3)' ) hfl(hnml)%mo
      hfl(hnml)%flag = 0
      if( (hfl(hnml)%mo .eq. 49) .and. (mo .eq.  7) ) hfl(hnml)%flag = 1 ! o2cia
      if( (hfl(hnml)%mo .eq. 52) .and. (mo .eq. 41) ) hfl(hnml)%flag = 2 ! n2cia
      write(6,*)''
      write(6,110) ' File : ', trim(filename)
      write(6,113) wstr, wstp, hnml, hfl(hnml)%lun, hfl(hnml)%mo, mo, wavnum

   12 continue

   enddo
   write(6,114) ' Number of HITRAN molecules/files found : ', hnml


   ! --- Dicke narrowing narrowing or Galatry data files - block 2 in hbin.input
   ! --- read number of expected Galatry files (max=2)
   ! --- Galatry files are unique format from hitran

   if (ctl_version.eq.1) then
      call nextbuf( ilun, buffer )
      read(buffer,*) ngal_files
   end if
   write(6,114) 'Number of Galatry files to search : ', ngal_files
   if( ngal_files .gt. ngal )stop 'too many hitran files'

   ! --- readlinelist_paths to (up to) 99 hitran files
   stlun = hfl(hnml)%lun
   lun   = stlun
   gnml  = 0
   snml  = 0
   do i = 1, ngal_files

      ! --- find the name of the next galatry input file
      if (ctl_version.eq.1) then
         call nextbuf( ilun, linebuffer )
         filename = trim(linelist_path) // trim(linebuffer)
      else
         linebuffer = trim(gal_files(i))
         filename = trim(linelist_path) // linebuffer
      end if

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

      lun = stlun + 1
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
print *, 1
print*, buffer
print *, 2
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
   if (ctl_version.eq.1) then
      call nextbuf( ilun, buffer )
      read(buffer,*) nlm_files
   end if
   write(6,114) 'Number of LM files to search : ', nlm_files
   if( nlm_files .gt. nlmx )stop 'too many LM files'

   ! --- readlinelist_paths to (up to) 99 hitran files
   lnml = 0
   lun = stlun
   do i = 1, nlm_files

      ! --- find the name of the next galatry input file
      if (ctl_version.eq.1) then
         call nextbuf( ilun, linebuffer )
         filename = trim(linelist_path) // trim(linebuffer)
      else
         filename = trim(linelist_path) // trim(lm_files(i))
         n = len_trim(filename)
      end if
      if( filename(n:n) .eq. '/' )cycle

      lnml = lnml + 1
      write(6,110) 'Found LineMix line file : ', trim(filename)

      ! --- open LineMix file if its here
      inquire( file=trim(filename), exist = fexist )
      if( .not. fexist )then
         write(*,*) 'LM file "', trim(filename), '" does not exist.'
         stop
      endif

      lun = lun + 1
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
   if (ctl_version.eq.1) then
      call nextbuf( ilun, buffer )
      read(buffer,*) nsdv_files
   end if
   write(6,114) 'Number of SDV files to search : ', nsdv_files
   if( nsdv_files .gt. nsdv )stop 'too many SDV files'

   ! --- read paths to (up to) 99 hitran files
   snml = 0
   lun = stlun
   do i = 1, nsdv_files

      ! --- find the name of the next galatry input file
      if (ctl_version.eq.1) then
         call nextbuf( ilun, linebuffer )
         filename = trim(linelist_path) // trim(linebuffer)
      else
         filename = adjustl(trim(linelist_path) // trim(sdv_files(i)))
      end if
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

      lun = lun + 1
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
!107 format(i2,i1,f12.6,1p,e10.3,10x,0p,f5.4,f5.4,f10.4,f4.2,f8.6,f7.4)
108 format(i2,1x,f12.6)
110 format( a, a )
111 format( a, i10 )
112 format( /, a, a )
113 format( 2f14.6, 4i6, 2x, f14.6 )
114 format( /, a, i10 )

end subroutine read_input


! --- read in bands from sfit.ctl file
subroutine read_ctrl

   use hitran
   use binput_4_0
   use binput_parse_4_0
   use params
   use datafiles
   use bandparam
   use lineparam
   use isotope


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
   write(6,*) ''

   ! --- Loop over bands
   do iband = 1, nband
      ! --- 10 res units to account for shifts - from initialize.f90:setup
      dwave  = resunits/pmax(iband)
      !dwave  = 10./pmax(iband)
      nextra = nint( dwave/dn(iband))
      !  --- interval for input of line data, dlines accounts for out of band absorption
      wave5(iband) = wave3(iband) - nextra*dn(iband) - dlines - 0.5/pmax(iband)
      wave6(iband) = wave4(iband) + nextra*dn(iband) + dlines + 0.5/pmax(iband)
      write(6,100) iband, wave3(iband), wave5(iband), wave4(iband), wave6(iband), pmax(iband), dn(iband)
   enddo

   write(6,*) ''
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

function qu_equal(quanta1, quanta2)
  ! compares quanta1 and quanta2 field in HITRAN 2004 format, retruns T if they are equal
  ! until now, only string compare, may get more complicated though
  implicit none
  character (len=*), intent(in):: quanta1, quanta2
  logical :: qu_equal

!print *, quanta1
!print *, quanta2


  qu_equal = .false.
  if (quanta1.eq.quanta2) qu_equal = .true.
  return

end function qu_equal
