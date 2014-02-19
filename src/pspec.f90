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

program pspec

! simple program to prepare the ascii spectral file that sfit4 expects
! reads ascii input file "pspec.inp"
! in that are names of binary spectra files and spectral window regions
! output is the 't15asc' file for tan sfit4 run.

! the work is done by spec.f90:
! each binary file is read & the window extracted, the spectral resolution is
! diminished and zero filled if requested, the spectra may be ratioed with
! a low resolution spectra if requested.

! parts previously developed by C.P.Rinsland, A.Goldman, P.Manning, F.Hase

   use spec
   use params
   use bandparam

   implicit none

   interface
      subroutine gonext(ifile, buffer)
         integer,intent(in)         :: ifile
         character(80),intent(out)  :: buffer
      end subroutine gonext
      subroutine read_ctrl
      end subroutine read_ctrl
   end interface

   character(len=80) :: bnrfile, ratfile, buffer
   character(len=1)  :: loc
   real          (8) :: obs_lat, obs_lon, obs_alt, roe, newlo, newhi, zflag
   integer       (4) :: nspec, j, n, m
   integer       (4) :: nterp, rflag, iband, bflag, fflag, oflag, vflag
   character (len=10)   :: ztime='          ', zone='          '
   character (len=8)    :: cdate='        '

   ! blun 7  - input bnr file
   ! clun 8  - cinput
   ! ilun 9  - input file pspec.inp
   ! rlun 10 - ratio file
   ! tlun 15 - t15asc
   ! olun 16 - simple output for fitbn tag & zero
   ! nlun 20 - output bnr file newratio.ispec

   ! bflag = 1 if a block was written to the t15 file

   ! loc is a location specifier for NZ data

! --- options to kpno

   ! fflag - file type
      ! = 0 open fortran unformatted,
      ! = 1 stream or binary or c-type binary file

   ! zflag  - zero offset
     ! = 0 no zero offset,
     ! = 1 try w/ baselincorrect,
     ! 0 < z < 1 use this value,
     ! = 2 use combo 2 + 4 for 10m region

   ! oflag  - output
     ! = 1 output t15asc file
     ! = 2 bnr file
     ! = 3 both

   ! vflag  - verbosity
     ! = 0 no output from baseline correct or s bnr or block output for plotting
     ! = 1 verbose output from bc and zeroed bnr but no blockout
     ! = 2 verbose, zeroed bnr and blockout for plotting

   ! nterp  - interpolation
     ! = 1 no interpolation of points
     ! = N interpolate xN points

   ! rflag   - envelope ratio
     ! = 1 ratio
     ! = 0 not

   call date_and_time (cdate, ztime, zone)
   write (tag,*) trim(version), ' runtime:', cdate(1:8), '-', ztime(1:2), ':', ztime(3:4), ':', ztime(5:6)
   write (6,  *) trim(tag)

   open (unit=ilun, file='pspec.input', status='old', err=668)
   call gonext(ilun, buffer)
   read (buffer,*) obs_lat

   call gonext(ilun, buffer)
   read (buffer,*) obs_lon

   call gonext(ilun, buffer)
   read (buffer,*) obs_alt

   call gonext(ilun, buffer)
   read (buffer,*) oflag, vflag

   call gonext(ilun, buffer)
   read (buffer,*) nsnr

   do j=1, nsnr
      call gonext(ilun, buffer)
      !print*, buffer
      read (buffer,*) snrid(j), psnr(1,j), psnr(2,j)
   enddo

   call gonext(ilun, buffer)
   read (buffer,*) nspec

   write(6,103) 'Observation Lat, Lon, Alt : ', obs_lat, obs_lon, obs_alt
   write(6,102) 'Output flag : ', oflag
   write(6,102) 'Verbose flag : ', vflag

   ! get band data from sfi4.ctl file
   call read_ctrl()

   ! open output file
   open(unit=tlun, file='t15asc.4', status='unknown', err=555)

   if( vflag .gt. 0 )open(unit=vlun,file='pspec_zero.dtl')

   ! loop over (# of spectra) x (# of windows)
   m = 0
   n = 0
   bflag = 0

   ! loop over bnr files
   do j=1, nspec

      call gonext(ilun, buffer)
      read(buffer,'(a80)') bnrfile
      call gonext(ilun, buffer)
      read(buffer,*) roe, nterp, rflag, fflag, zflag

      ! each t15asc block requires radius of earth, solar zenith angle, lat, lon,
      ! date and time
      ! ckopus.c puts most of this in the bnr header except lat, lon and roe
      ! hence these are given here in pspec.inp

      write(6,101) "Opening bnr for input : ", bnrfile(1:len_trim(bnrfile))
      if( fflag .eq. 0 )then
         open(unit=blun, file=bnrfile(1:len_trim(bnrfile)), form='unformatted', status='old', err=666)
      elseif( fflag .eq. 1 )then
         open(unit=blun, file=bnrfile(1:len_trim(bnrfile)), form='unformatted', access='stream', status='old', err=666)
      else
         print*, ' file flag out of range 0 - fortran type or 1 c-type.'
         stop '1'
      endif

      if( rflag .eq. 1) then
         call gonext(ilun, buffer)
         read(buffer,'(a80)') ratfile
         write(6,105) 'Opening ratio file : ', ratfile(1:len_trim(ratfile))
         if( fflag .eq. 0 )then
            open(unit=rlun, file=ratfile(1:len_trim(ratfile)), form='unformatted', status='old', err=666)
         elseif( fflag .eq. 1 )then
            open(unit=rlun, file=ratfile(1:len_trim(bnrfile)), form='unformatted', access='stream', status='old', err=666)
         else
            print*, ' file flag out of range 0 - fortran type or 1 c-type.'
            stop '2'
         endif
      endif

      ! loop over bands for this bnr file
      do iband = 1, nband
         n = n +1
         write(6,104) "Band : ", n
         newlo = wave3(iband)-.1
         newhi = wave4(iband)+.1

         call kpno( pmax(iband), newlo, newhi, roe, obs_lat, obs_lon, nterp, rflag, oflag, zflag, vflag, bflag, loc )

         if( bflag .eq. 0 ) m = m +1
         rewind(blun)
         rewind(rlun)
      enddo

      close(blun)
      close(rlun)

   enddo

   close(tlun)   ! t15asc file
   close(blun)   ! data bnr file
   close(ilun)   ! input ascii file
   close(rlun)   ! ratio bnr file
   if( vflag .gt. 0 )close( vlun ) ! zero correction verbose output

   write(6,102)'Blocks written to t15asc : ', m

   stop 'pspec .done.'

! 1001 format( a )
  101 format(/,a32,a)
!1002 format( * )
  102 format(a32,4i20)
  103 format(a32,5f20.6)
  104 format(/,a32,4i20)
  105 format(a32,a)
!  208 format(2f10.5)
!  209 format(f5.0)
!  212 format(2i5)
!  500 format(a6)
  555 write ( *, 556)
  556 format(' t15asc error from pspec')
      stop 'pspec'
  666 write ( *, 667)
  667 format(' abort-binary file open error called from pspec')
      stop 'pspec'
  668 write(*,669)
  669 format(' abort-input file open error called from pspec')
      stop 'pspec'

! 1002 format(8f10.3)

! 1008 format(f10.3,7i5)
      stop 'pspec'
      end program pspec

! --- read in bands from sfit.ctl file
subroutine read_ctrl

   use params
   use datafiles
   use bandparam
   use lineparam
   use isotope
   use binput_4_0
   use binput_parse_4_0


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
   write(6,101) 'Number of bands from .ctl : ', nband

   ! --- Loop over bands
   do iband = 1, nband
      ! --- 10 res units to account for shifts - from initialize.f90:setup
      dwave  = 10.d0/pmax(iband)
      nextra = nint( dwave/dn(iband))
      !  --- interval for input of line data, dlines accounts for out of band absorption
      wave5(iband) = wave3(iband) - nextra*dn(iband) - dlines
      wave6(iband) = wave4(iband) + nextra*dn(iband) + dlines
      write(6,100) iband, wave3(iband), wave4(iband), pmax(iband)
   enddo

   close( bp_nr )

return

100 format( 32x,i5,4f14.5,f12.3,f12.6 )
101 format( /, a32, i20 )

end subroutine read_ctrl

!--------------------------------------------------------------------------------------------
subroutine gonext(ifile, buffer)

   implicit none
   integer,intent(in)         :: ifile
   character(80),intent(out)  :: buffer

   buffer(1:1) = '#'
   do while (buffer(1:1) .eq. '#')
       read(ifile,'(a80)') buffer
       !print *, buffer
   end do

end subroutine gonext

