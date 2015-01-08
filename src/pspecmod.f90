MODULE pspec_read_mods

   use params
   use datafiles
   use bandparam
   use lineparam
   use isotope
   use binput_4_0
   use binput_parse_4_0
   use strings

   implicit none

   !------------------------------------
   ! Declare derived data type for array 
   ! file names
   !------------------------------------
   type char_mat
       character(len=200) :: fname
   end type char_mat

   !----------------------
   ! Set member visibility
   !----------------------
   !private            ::   gonext_skip


CONTAINS

!------------------------------------------
! Subroutine for old style pspec input file
!------------------------------------------
subroutine pspec_read_old(ilun,            &   ! Input
                          obs_lat_s,       &   ! Output
                          obs_lon_s,       &   ! Output
                          obs_alt_s,       &   ! Output
                          oflag_s,         &   ! Output
                          vflag_s,         &   ! Output
                          nsnr_s,          &   ! Output
                          psnr_s,          &   ! Output
                          roe_s,           &   ! Output
                          nterp_s,         &   ! Output
                          rflag_s,         &   ! Output
                          fflag_s,         &   ! Output
                          zflag_s,         &   ! Output
                          bnr_fname_mat_s, &   ! Output
                          ratio_fname_mat_s)   ! Output
                          

   implicit none

   character(len=80)                    :: buffer
   real(double),intent(out)             :: obs_lat_s, obs_lon_s, obs_alt_s
   real(double),intent(out),allocatable :: roe_s(:)
   real(8)                              :: psnr_s(2,100)    ! This has to dimensioned the same as psnr in spec.f90
   integer,intent(out),allocatable      :: nterp_s(:), rflag_s(:), fflag_s(:)
   real(8), allocatable                 :: zflag_s(:)
   integer,intent(out)                  :: oflag_s, vflag_s, nsnr_s
   integer                              :: n_spec, j
   character(len=2)                     :: dum
   integer, intent(in)                  :: ilun

   type (char_mat), intent(out), dimension (:), allocatable :: bnr_fname_mat_s
   type (char_mat), intent(out), dimension (:), allocatable :: ratio_fname_mat_s

   call gonext_skip(ilun, buffer)
   read (buffer,*) obs_lat_s

   call gonext_skip(ilun, buffer)
   read (buffer,*) obs_lon_s

   call gonext_skip(ilun, buffer)
   read (buffer,*) obs_alt_s

   call gonext_skip(ilun, buffer)
   read (buffer,*) oflag_s, vflag_s

   call gonext_skip(ilun, buffer)
   read (buffer,*) nsnr_s
 
   do j=1, nsnr_s
      call gonext_skip(ilun, buffer)
      !print*, buffer
      read (buffer,*) dum, psnr_s(1,j), psnr_s(2,j)
   enddo

   call gonext_skip(ilun, buffer)
   read (buffer,*) n_spec

   call read_ctrl()

   allocate (bnr_fname_mat_s(n_spec), ratio_fname_mat_s(n_spec), roe_s(n_spec), nterp_s(n_spec), &
             rflag_s(n_spec), fflag_s(n_spec), zflag_s(n_spec))

   do j=1, n_spec
      call gonext_skip(ilun, buffer)
      read(buffer,'(a80)') bnr_fname_mat_s(j)%fname
      call gonext_skip(ilun, buffer)
      read(buffer,*) roe_s(j), nterp_s(j), rflag_s(j), fflag_s(j), zflag_s(j)

      if( rflag_s(j) .eq. 1) then
         call gonext_skip(ilun, buffer)
         read(buffer,'(a80)') ratio_fname_mat_s(j)%fname
      end if
   end do 


end subroutine pspec_read_old

!-----------------------------------------------
! Subroutine to read .ctl style pspec input file
!-----------------------------------------------
subroutine read_pspec_ctl(ilun,            &   ! Input
                          obs_lat_s,       &   ! Output
                          obs_lon_s,       &   ! Output
                          obs_alt_s,       &   ! Output
                          oflag_s,         &   ! Output
                          vflag_s,         &   ! Output
                          nsnr_s,          &   ! Output
                          psnr_s,          &   ! Output
                          roe_s,           &   ! Output
                          nterp_s,         &   ! Output
                          rflag_s,         &   ! Output
                          fflag_s,         &   ! Output
                          zflag_s,         &   ! Output
                          bnr_fname_mat_s, &   ! Output
                          ratio_fname_mat_s)   ! Output


   implicit none

   integer, intent(in)                      :: ilun
   
   real(double),intent(out)                 :: obs_lat_s, obs_lon_s, obs_alt_s
   real(double),intent(out),allocatable     :: roe_s(:)
   real(8)                                  :: psnr_s(2,100)    ! This has to dimensioned the same as psnr in spec.f90
   integer,intent(out),allocatable          :: nterp_s(:), rflag_s(:), fflag_s(:)
   real(8), allocatable                     :: zflag_s(:)
   integer,intent(out)                      :: oflag_s, vflag_s, nsnr_s

   integer, parameter                       :: StrMax=5, Nmax=20
   character (len=StrMax), dimension(Nmax)  :: args
   character (len=255)                      :: val
   integer                                  :: n_spec, snr_loc, bnr_loc, i
   integer                                  :: istat, nr_keys, iband, nextra, band_flg,nargs
   real(double)                             :: dwave

   type (char_mat), intent(out), dimension (:), allocatable :: bnr_fname_mat_s
   type (char_mat), intent(out), dimension (:), allocatable :: ratio_fname_mat_s

   !----------------------------------------------------
   ! pspec file already opened in calling program (ilun)
   ! Read in input parameters
   !----------------------------------------------------
   band_flg = 0
   snr_loc  = -1
   bnr_loc  = -1

   do
      call read_line_binput( ilun, keyword, nr_keys, value, istat )

      if(( istat .lt. 0 ) .and. ( nr_keys .eq. 0 )) exit
      if( nr_keys .eq. 0 ) cycle
      !print*,'keyword = ', keyword(3)
      !print*,'Value = ',value
      !print*,'nr_keys  = ',nr_keys
      !print*,'istat  = ',istat

      select case ( trim( adjustl( keyword( 1 ))))
         case ( 'band' )
            band_flg = 1
            call read_band_section( keyword, value )

         case ( 'loc' )   ! Location information
             select case (trim(adjustl(keyword(2))))
                 case ('lat')
                     read(value,*) obs_lat_s
                 case ('lon')
                     read(value,*) obs_lon_s
                 case ('alt')
                     read(value,*) obs_alt_s
             end select
   
         case ( 'verb' )  ! Verbosity flags
             select case (trim(adjustl(keyword(2))))
                 case ('oflag')
                     read(value,*) oflag_s
                 case ('vflag')
                     read(value,*) vflag_s
             end select

         case ('snr')   ! SNR windows
             if (len_trim(keyword(2)).eq.0) then
                val = adjustl(trim(value))
                call parse(val,' ',args,nargs)
                nsnr_s = nargs

                if ( nargs == 0 ) then
                    print*,'Pspec error: No SNR windows specified!!'
                    stop
                end if

             else
                 do i = 1,nsnr_s
                   if ( trim(adjustl(keyword(2))) == trim(adjustl(args(i))) ) snr_loc = i
                 end do

                 select case (trim(adjustl(keyword(3))))
                     case ('nu_start')
                        if ( snr_loc == -1) then
                           print*,'Error in snr format in pspec.ctl file!!'
                           stop
                        end if
                        read(value, *) psnr_s(1,snr_loc)
                     case ('nu_stop')
                        read(value, *) psnr_s(2,snr_loc)
                     case default
                        print*,     'Key ', trim(keyword(3)), ' not contained in pspec section : snr'
                        write(16,*) 'Key ', trim(keyword(3)), ' not contained in pspec section : snr'
                 end select
             end if

         case ('bnr')   ! BNR files
             if (len_trim(keyword(2)).eq.0) then
                val = adjustl(trim(value))
                call parse(val,' ',args,nargs)
                n_spec = nargs

                if ( nargs == 0 ) then
                    print*,'Pspec error: No BNR files specified!!'
                    stop
                end if

                allocate (bnr_fname_mat_s(n_spec), ratio_fname_mat_s(n_spec), roe_s(n_spec), nterp_s(n_spec), &
                          rflag_s(n_spec), fflag_s(n_spec), zflag_s(n_spec))
                
             else
                 do i = 1,n_spec
                   if ( trim(adjustl(keyword(2))) == args(i)) bnr_loc = i
                 end do

                 select case (trim(adjustl(keyword(3))))
                     case ('file')
                        if ( bnr_loc == -1) then
                           print*,'Error in bnr format in pspec.ctl file!!'
                           stop
                        end if
                        read(value, *) bnr_fname_mat_s(bnr_loc)%fname 
                     case ('roe')
                        read(value, *) roe_s(bnr_loc)
                     case ('nterp')
                        read(value, *) nterp_s(bnr_loc)
                     case ('rflag')
                        read(value, *) rflag_s(bnr_loc)
                     case ('fflag')
                        read(value, *) fflag_s(bnr_loc)
                     case ('zflag')
                        read(value, *) zflag_s(bnr_loc)
                     case ('ratio')
                        read(value, *) ratio_fname_mat_s(bnr_loc)%fname 
                     case default
                        print*,     'Key ', trim(keyword(3)), ' not contained in pspec section : bnr'
                        write(16,*) 'Key ', trim(keyword(3)), ' not contained in pspec section : bnr'
                 end select
             end if
      end select
   end do

   !---------------------------------------------------------------
   ! This is typically done in read ctl subroutine; however, we are
   ! assuming no ctl file in this case
   !---------------------------------------------------------------
   if ( band_flg == 1) then
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
   end if

   !--------------------------------------------------
   ! If unable to find band information in pspec input
   ! file then try sfit.ctl file
   !--------------------------------------------------
   if ( band_flg == 0 ) then
       print*,'Going here?!?!?!?!'
       call read_ctrl()
       band_flg = 1

   end if

   !---------------------------------------------
   ! Test to see if any band information is found
   !---------------------------------------------
   if ( band_flg == 0 ) then
       print*,'Pspec Error: Unable to find band information in pspec.ctl file or sfit.ctl file...exiting'
       stop
   end if

   write(6,101) 'Number of bands from .ctl : ', nband

   close( bp_nr )

100 format( 32x,i5,4f14.5,f12.3,f12.6 )
101 format( /, a32, i20 )

end subroutine read_pspec_ctl

!-----------------------------------------------
! Subroutine to read in bands from sfit.ctl file
!-----------------------------------------------
subroutine read_ctrl

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
      call read_line_binput( bp_nr, keyword, nr_keys, value, istat )
      if(( istat .lt. 0 ) .and. ( nr_keys .eq. 0 ))exit
      if( nr_keys .eq. 0 )then
         cycle
      endif
      select case ( trim( adjustl( keyword( 1 ))))
         case ( 'band' )
            call read_band_section( keyword, value )
      end select
   end do

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
   100 format( 32x,i5,4f14.5,f12.3,f12.6 )

end subroutine read_ctrl

!-------------------------------------------------
! Subroutine to read a line in old style input 
! file and skip comments '#'
!-------------------------------------------------
subroutine gonext_skip(ifile, buffer)

   implicit none
   integer,intent(in)         :: ifile
   character(80),intent(out)  :: buffer

   buffer(1:1) = '#'
   do while (buffer(1:1) .eq. '#')
       read(ifile,'(a80)') buffer
       !print *, buffer
   end do

end subroutine gonext_skip

!-------------------------------------------
! Subroutine to convert string to lower case
!-------------------------------------------
subroutine To_lower(str)
   
   character(*), intent(in out) :: str
   integer :: i
 
   do i = 1, len(str)
     select case(str(i:i))
       case("A":"Z")
         str(i:i) = achar(iachar(str(i:i))+32)
     end select
   end do  
end subroutine To_Lower

!---------------------
! Clean up  subroutine
!---------------------
subroutine cleanup(tlun,blun,ilun,rlun,vflag,vlun)
  
   implicit none
   integer, intent(in) :: tlun,blun,ilun,rlun,vflag,vlun

   close(tlun)   ! t15asc file
   close(blun)   ! data bnr file
   close(ilun)   ! input ascii file
   close(rlun)   ! ratio bnr file
   if( vflag .gt. 0 )close( vlun ) ! zero correction verbose output

end subroutine cleanup

END MODULE pspec_read_mods
