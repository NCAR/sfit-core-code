program pspec

! simple program to prepare the ascii spectral file that sfit4 expects
! reads ascii input file "prepspec.inp"
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

character(len=80) :: bnrfile, ratfile, buffer
character(len=1)  :: loc
real          (8) :: obs_lat, obs_lon, obs_alt, roe !, tag
integer       (4) :: nspec, j, n, m !obs_yy, obs_mm, obs_dd, obs_hh, obs_nn, obs_ss
integer       (4) :: nterp, rflag, iband, flag, fflag

! lun 8  - input file - prepspec.inp
! lun 20 - output file - t15asc
! lun 7  - input bnr file
! lun 9  - ratio file (in bnr format)
! lun 15 - output bnr file newratio.ispec (last block only)

! fflag - file open 0 fortran unformatted, 1 stream or binary or c-type
! rflag-  ratio with envelope function

open (unit=8, file='pspec.input', status='old', err=668)
call gonext(8, buffer)
read (buffer,*) obs_lat
call gonext(8, buffer)
read (buffer,*) obs_lon
call gonext(8, buffer)
read (buffer,*) obs_alt
call gonext(8, buffer)
read (buffer,*) nsnr
do j=1, nsnr
   call gonext(8, buffer)
   print*, buffer
   read (buffer,*) snrid(j), psnr(1,j), psnr(2,j)
enddo
call gonext(8, buffer)
read (buffer,*) nspec

! get band data from sfi4.ctl file
call read_ctrl()

! open output file
open(unit=20, file='t15asc.4', status='unknown', err=555)

! loop over (# of spectra) x (# of windows)
m = 0
n = 0
do j=1, nspec

   call gonext(8, buffer)
   read(buffer,'(a80)') bnrfile
   call gonext(8, buffer)
   read(buffer,*) roe, nterp, rflag, fflag

   ! each t15asc block requires radius of earth, solar zenith angle, lat, lon,
   ! date and time
   ! ckopus.c puts most of this in the bnr header except lat, lon and roe
   ! hence these are given here in prepspec.inp

   !call gonext(8, buffer)
   !read(buffer,*) obs_yy, obs_mm, obs_dd, obs_hh, obs_nn, obs_ss

   if( rflag .eq. 1) then
      call gonext(8, buffer)
      read(buffer,'(a80)') ratfile
   endif

   if( fflag .eq. 0 )then
      print *, " Opening bnr for input : ", bnrfile(1:len_trim(bnrfile))
      open(unit=7, file=bnrfile(1:len_trim(bnrfile)), form='unformatted', status='old', err=666)
   elseif( fflag .eq. 1 )then
      print *, " Opening bnr for input : ", bnrfile(1:len_trim(bnrfile))
      open(unit=7, file=bnrfile(1:len_trim(bnrfile)), form='unformatted', access='stream', status='old', err=666)
   else
      print*, ' file flag out of range 0 - fortran type or 1 c-type.'
      stop 1
   endif

   if( rflag .eq. 1) then
      print *, ' Opening ratio file : ', ratfile(1:len_trim(ratfile))
      open(unit=9, file=ratfile(1:len_trim(ratfile)), form='unformatted', status='old', err=666)
   endif

   do iband = 1, nband
      n = n +1
      print *, n, " : Input : ", bnrfile(1:len_trim(bnrfile))
      call kpno( pmax(iband), wave3(iband)-.1, wave4(iband)+.1, roe, obs_lat, obs_lon, nterp, rflag, flag, loc )
      if( flag .eq. 0 ) m = m +1
      rewind(7)
      rewind(9)
   enddo

   close(7)
   close(9)

enddo

close(20)   ! t15asc file
close(7)    ! data bnr file
close(8)    ! input ascii file
close(9)    ! ratio bnr file

print *,'  Blocks written to t15asc : ', m
stop 'prepspec .done.'


! 1001 format( a )
!  101 format(a1)
!1002 format( * )
!  102 format(4i5,a20)
!  103 format( 2(a, i3), 2f12.4)
!  104 format( f16.5 )
!  208 format(2f10.5)
!  209 format(f5.0)
!  212 format(2i5)
!  500 format(a6)
  555 write ( *, 556)
  556 format(' t15asc error from prepspec')
      stop 'prepspec'
  666 write ( *, 667)
  667 format(' abort-binary file open error called from prepspec')
      stop 'prepspec'
  668 write(*,669)
  669 format(' abort-input file open error called from prepspec')
      stop 'prepspec'

! 1002 format(8f10.3)

! 1008 format(f10.3,7i5)
      stop 'prepspec'
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
   write(6,101) 'Found number of bands : ', nband

   ! --- Loop over bands
   do iband = 1, nband
      ! --- 10 res units to account for shifts - from initialize.f90:setup
      dwave  = 10.d0/pmax(iband)
      nextra = nint( dwave/dn(iband))
      !  --- interval for input of line data, dlines accounts for out of band absorption
      wave5(iband) = wave3(iband) - nextra*dn(iband) - dlines
      wave6(iband) = wave4(iband) + nextra*dn(iband) + dlines
      write(6,100) iband, wave3(iband), wave4(iband), pmax(iband), dn(iband)
   enddo

   close( bp_nr )

return

100 format( i5,4f14.5,f12.3,f12.6 )
101 format( /, a, i10 )

end subroutine read_ctrl

!--------------------------------------------------------------------------------------------
subroutine gonext(ifile, buffer)

implicit none
integer,intent(in) :: ifile
character(80),intent(out) :: buffer

buffer(1:1) = '#'
do while (buffer(1:1) .eq. '#')
    read(ifile,'(a80)') buffer
    !print *, buffer
end do

end subroutine gonext

