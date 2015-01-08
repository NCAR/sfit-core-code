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
   use pspec_read_mods

   implicit none

   !------------------
   ! Type declerations
   !------------------
   real(double)                    :: obs_lat, obs_lon, obs_alt, newhi, newlo
   real(double),allocatable        :: roe(:)
   real(8), allocatable            :: zflag(:)
   integer,allocatable             :: nterp(:), rflag(:), fflag(:)
   integer                         :: oflag, vflag, arg_cnt, file_rtncode, nchar
   integer                         :: m, n, j, bflag, n_spec, iband
   character (len=200)             :: infile
   character(len=1)                :: loc
   character (len=10)              :: ztime='          ', zone='          '
   character (len=8)               :: cdate='        '
   logical                         :: fexist

   !-----------------------------------------------------------------------
   ! Deferred length character components only supported in gfortran v >4.9
   ! For now set character len = 200
   !-----------------------------------------------------------------------
   type (char_mat), dimension (:), allocatable :: bnr_fname_mat
   type (char_mat), dimension (:), allocatable :: ratio_fname_mat
   
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

   !-----------------------------
   ! Check command line arguments
   !-----------------------------
   arg_cnt = command_argument_count()

   if ( arg_cnt /= 1 ) then
       print*,'Usage: pspec [input file]'
       stop
   endif

   !-------------------------
   ! Grab input file and open
   !-------------------------
   call get_command_argument(1,infile)

   open (unit=ilun, file=adjustl(trim(infile)), status='old', err=668,iostat=file_rtncode)
   if ( file_rtncode /= 0 ) then
       print*,'Unable to open file: '// trim(infile)
       stop
   endif 

   !---------------------------------------
   ! Determine file type based on extension
   !---------------------------------------
   call to_lower(infile)
   nchar = len(trim(infile))

   !---------------------------------------
   ! Read input file depending of file type
   !---------------------------------------
   if ( infile(nchar-4:nchar) == 'input' ) then             ! Old input style
       call pspec_read_old(ilun,obs_lat,obs_lon,obs_alt,oflag,vflag,nsnr,psnr,roe,nterp,rflag,fflag,zflag,&
                           bnr_fname_mat,ratio_fname_mat)

    else if ( infile(nchar-2:nchar) == 'ctl' ) then       ! .ctl input style
        call read_pspec_ctl(ilun,obs_lat,obs_lon,obs_alt,oflag,vflag,nsnr,psnr,roe,nterp,rflag,fflag,zflag,&
                           bnr_fname_mat,ratio_fname_mat)

    else
        print*,'Unrecognized pspec input file type. File extension should be .ctl or .input!!' 
        stop
    end if

   !--------------------------------
   ! Open output file and write each
   ! ascii block
   !--------------------------------
   open(unit=tlun, file='t15asc.4', status='unknown', err=555)

   if( vflag .gt. 0 )open(unit=vlun,file='pspec_zero.dtl')

   ! loop over (# of spectra) x (# of windows)
   m      = 0
   n      = 0
   bflag  = 0
   n_spec = size(bnr_fname_mat)

   ! loop over bnr files
   do j=1, n_spec

      ! each t15asc block requires radius of earth, solar zenith angle, lat, lon,
      ! date and time
      ! ckopus.c puts most of this in the bnr header except lat, lon and roe
      ! hence these are given here in pspec.inp

      write(6,101) "Opening bnr for input : ", adjustl(trim(bnr_fname_mat(j)%fname))
      inquire( file=adjustl(trim(bnr_fname_mat(j)%fname)), exist = fexist )
      if( .not. fexist ) then
          write(*,*) 'file ', adjustl(trim(bnr_fname_mat(j)%fname)), ' does not exist'
          stop
      endif

      if( fflag(j) .eq. 0 )then
         open(unit=blun, file=adjustl(trim(bnr_fname_mat(j)%fname)), form='unformatted', status='old', err=666)
      elseif( fflag(j) .eq. 1 )then
         open(unit=blun, file=adjustl(trim(bnr_fname_mat(j)%fname)), form='unformatted', access='stream', status='old', err=666)
      else
         print*, ' file flag out of range 0 - fortran type or 1 c-type.'
         stop '1'
      endif

      if( rflag(j) .eq. 1) then
         write(6,105) 'Opening ratio file : ', ratio_fname_mat(j)%fname

         if( fflag(j) .eq. 0 )then
            open(unit=rlun, file=adjustl(trim(ratio_fname_mat(j)%fname)), form='unformatted', status='old', err=666)
         elseif( fflag(j) .eq. 1 )then
            open(unit=rlun, file=adjustl(trim(ratio_fname_mat(j)%fname)), form='unformatted', access='stream', status='old', err=666)
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

         call kpno( pmax(iband), newlo, newhi, roe(j), obs_lat, obs_lon, nterp(j), &
                    rflag(j), oflag, zflag(j), vflag, bflag, loc )

         if( bflag .eq. 0 ) m = m +1
         rewind(blun)
         rewind(rlun)
      enddo

      close(blun)
      close(rlun)

   enddo

   !--------------------------
   ! Screen output information
   !--------------------------
   call date_and_time (cdate, ztime, zone)
   write (tag,*) trim(version), ' runtime:', cdate(1:8), '-', ztime(1:2), ':', ztime(3:4), ':', ztime(5:6)
   write (6,  *) trim(tag)
   
   write(6,103) 'Observation Lat, Lon, Alt : ', obs_lat, obs_lon, obs_alt
   write(6,102) 'Output flag : ', oflag
   write(6,102) 'Verbose flag : ', vflag
   write(6,102) 'Blocks written to t15asc : ', m
   write(6,102) 'Number of bands from .ctl : ', nband

   !--------------
   ! File clean up
   !--------------
   call cleanup(tlun,blun,ilun,rlun,vflag,vlun)
   stop 'pspec .done.'

!--------------
! Write Formats
!--------------
  101 format( ' ', a100 )
  102 format( a40, i20 )
  103 format( a40, 3f20.6)
  104 format(/,a32,4i20)
  105 format(a32,a)
  555 write ( *, 556)
  556 format(' t15asc error from pspec')
      stop 'pspec'
  666 write ( *, 667)
  667 format(' abort-binary file open error called from pspec')
      stop 'pspec'
  668 write(*,669)
  669 format(' abort-input file open error called from pspec')
      stop 'pspec'

end program pspec

