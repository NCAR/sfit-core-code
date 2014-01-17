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

module binput_parse_4_0

  use params
  use retvparam
  use frwdmdl
  use bandparam
  use lineparam
  use solar
  use channel
  use initialize
  use opt
  use datafiles
  use isotope
  use writeout

  implicit none;
  save

  character (len=255), dimension(5) :: keyword
  character (len=1023) :: value
  character (len=7), dimension(10) :: gas_prf, gas_col
  logical, dimension(10) :: gas_detail=.false.
  logical :: f_gasprf=.false., f_gascol=.false.
  integer, dimension(maxbnd) :: nbeam
  integer                    :: nprf
contains

  subroutine read_file_section(keyword, value)

    implicit none

    character (len=*), dimension(*),intent(in) :: keyword
    character (len=*),intent(in) :: value

    if( trim(adjustl(keyword(2))) .eq. 'in' )then

       select case(trim(adjustl(keyword(3))))
          case ('stalayers')
             tfile(71) = trim(adjustl(value))
          case ('spectrum')
             tfile(15) = trim(adjustl(value))
          case ('modulation_fcn')
             tfile(23) = trim(adjustl(value))
          case ('phase_fcn')
             tfile(24) = trim(adjustl(value))
          case ('sa_matrix')
             tfile(62) = trim(adjustl(value))
          case ('isotope')
             tfile(09) = trim(adjustl(value))
          case ('linelist')
             tfile(14) = trim(adjustl(value))
          case ('solarlines')
             tfile(10) = trim(adjustl(value))
          case ('refprofile')
             tfile(72) = trim(adjustl(value))
          case default
             WRITE(16,*) 'BINPUT_PARSE_4_0:READ_FILE_SECTION: Key ', &
                  trim(keyword(3)), ' not contained in section file.in'
             WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_FILE_SECTION: Key ', &
                  trim(keyword(3)), ' not contained in section file.in'
             CALL SHUTDOWN
             STOP 1
          end select
          
       elseif( trim(adjustl(keyword(2))) .eq. 'out' )then

       select case(trim(adjustl(keyword(3))))
          case ('solarspectrum')
             TFILE(11) = trim(adjustl(value))
          case ('summary')
             tfile(20) = trim(adjustl(value))
          case ('pbpfile')
             tfile(08) = trim(adjustl(value))
          case ('statevec')
             tfile(18) = trim(adjustl(value))
          case ('k_matrix')
             TFILE(66) = trim(adjustl(value))
          case ('g_matrix')
             TFILE(93) = trim(adjustl(value))
!          case ('kout_matrices')
!             TFILE(66) = trim(adjustl(value))
          case ('shat_matrix')
             TFILE(64) = trim(adjustl(value))
          case ('sa_matrix')
             TFILE(63) = trim(adjustl(value))
          case ('retprofiles')
             TFILE(88) = trim(adjustl(value))
          case ('aprprofiles')
             TFILE(87) = trim(adjustl(value))
          case ('ak_matrix')
             TFILE(81) = trim(adjustl(value))
          case ('ab_matrix')
             TFILE(92) = trim(adjustl(value))
          case ('parm_vectors')
             TFILE(89) = trim(adjustl(value))
          case ('seinv_vector')
             TFILE(67) = trim(adjustl(value))
          case ('sainv_matrix')
             TFILE(69) = trim(adjustl(value))
          case ('smeas_matrix')
             TFILE(82) = trim(adjustl(value))
          case ('ssmooth_matrix')
             TFILE(83) = trim(adjustl(value))
          case ('kb_matrix')
             TFILE(90) = trim(adjustl(value))
          case default
             WRITE(16,*) 'BINPUT_PARSE_4_0:READ_FILE_SECTION: Key ', &
                  trim(keyword(3)), ' not contained in section file.out'
             WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_FILE_SECTION: Key ', &
                  trim(keyword(3)), ' not contained in section file.out'
             CALL SHUTDOWN
             STOP 1
          end select
       else
       WRITE(16,*) 'BINPUT_PARSE_4_0:READ_FILE_SECTION: Key ', &
            trim(keyword(2)), ' not contained in section file'
       WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_FILE_SECTION: Key ', &
            trim(keyword(2)), ' not contained in section file'
       CALL SHUTDOWN
       STOP 1
    endif
 
end subroutine read_file_section


  subroutine read_gas_section(keyword, value)

    implicit none

    character (len=*), dimension(*),intent(in) :: keyword
    character (len=*),intent(in) :: value
    character (len=255) :: val
    integer             :: pos=0, nr=0, nr1=0, ncol = 0
    logical             :: flag

    val = trim(adjustl(value))

    if (trim(keyword(2)).eq.'layers')then
           !print*,8
       read(value,*) nlayers
    end if
    !print *, 'nlayers read : ', nlayers

!  --- profile section

    if (trim(keyword(2)).eq.'profile') then
       select case (trim(adjustl(keyword(3))))
       case ('list')
          if(f_gasprf) then
             WRITE(16,*) 'BINPUT_PARSE_4_0:READ_GAS_SECTION: gas.profile.list ALREADY GIVEN'
             WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_GAS_SECTION: gas.profile.list ALREADY GIVEN'
             CALL SHUTDOWN
             STOP 1
          end if
          nprf = 0
          pos = index(adjustl(val),' ')
   !       write(*,*) val, pos
          if (pos.eq.0) write(*,*) 'No gas given in binput?'
          do
             if (len_trim(val).eq.0) exit
             nprf = nprf + 1
             if (pos.gt.0) then
                gas_prf(nprf) = trim(adjustl(val(1:pos)))
             else
                gas_prf(nprf) = trim(adjustl(val(1:len_trim(val))))
                exit
             end if
             val = adjustl(val(pos+1:len(val)))
             pos = index(trim(adjustl(val)),' ')
          end do
          ! In case columns have already been read in, update parameters in NRET arrays
          gas(nprf+1:nret+nprf)    = gas(1:nret)
          gas(1:nprf)              = gas_prf(1:nprf)
          colsf(nprf+1:nret+nprf)  = colsf(1:nret)
          colsf(1:nret)            = 1.0d0
          scolsf(nprf+1:nret+nprf) = scolsf(1:nret)
          nret                     = nret+nprf
          ifprf(1:nprf)            = .true.
          f_gasprf                 = .true.

       end select

       flag = .false.
       do nr=1,nprf
          if (trim(adjustl(gas_prf(nr))).eq.trim(adjustl(keyword(3)))) then
             flag = .true.
             exit
          end if
       end do

       if (.not.flag) then
          return
       end if

       select case (trim(adjustl(keyword(4))))
       case ('correlation')
          if (len_trim(keyword(5)).eq.0) then
             read(value,*) correlate(nr)
          else
             select case (trim(adjustl(keyword(5))))
             case ('type')
                read(value,*) ifoff(nr)
             case ('width')
                read(value,*) zwid(nr)
             case ('minalt')
                read(value,*) zgmin(nr)
             case ('maxalt')
                read(value,*) zgmax(nr)
             case default
                WRITE(16,*) 'BINPUT_PARSE_4_0:READ_GAS_SECTION: Key ', trim(keyword(4)), &
                            ' not contained in section gas...correlation'
                WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_GAS_SECTION: Key ', trim(keyword(4)), &
                            ' not contained in section gas...correlation'
                CALL SHUTDOWN
                STOP 1
             end select
           endif
       case ('scale')
         read(value,*) colsf(nr)
       case ('sigma')
         read(value,*) sig(1:nlayers, nr)
       case ('logstate')
          !print*,2
          read(value,*) log_statev(nr)
          if (log_statev(nr))then
             ilogretrieval(nr) = 1
          else
             ilogretrieval(nr) = 0
          end if
       end select
    endif

!  --- column section

    if (trim(keyword(2)).eq.'column') then

       select case (trim(adjustl(keyword(3))))
       case ('list')
       if(f_gascol) then
          WRITE(16,*) 'BINPUT_PARSE_4_0:READ_GAS_SECTION: Key ', trim(keyword(4)), ' not contained in section gas...correlation'
          WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_GAS_SECTION: Key ', trim(keyword(4)), ' not contained in section gas...correlation'
          CALL SHUTDOWN
          STOP 1
       end if
       pos = index(adjustl(val),' ')
          !       write(*,*) val, pos
          if (pos.eq.0) write(*,*) 'No gas given in binput?'
          do
             if (len_trim(val).eq.0) exit
             ncol = ncol + 1
             if (pos.gt.0) then
                gas_col(ncol) = trim(adjustl(val(1:pos)))
             else
                gas_col(ncol) = trim(adjustl(val(1:len_trim(val))))
                exit
             end if
             val = adjustl(val(pos+1:len(val)))
             pos = index(trim(adjustl(val)),' ')
          end do
          gas(nret+1:nret+ncol) = gas_col(1:ncol)
          nret = nret + ncol
          f_gascol = .true.
       end select

       flag = .false.
       do nr=1,ncol
          if (trim(adjustl(gas_col(nr))).eq.trim(adjustl(keyword(3)))) then
             flag = .true.
             nr1 = nr
             exit
          end if
       end do
       if (.not.flag) then
          return
       end if

       nr1 = nr1 + nprf

       select case (trim(adjustl(keyword(4))))
       case ('sigma')
         read(value,*) scolsf(nr1)
       case ('scale')
         read(value,*) colsf(nr1)
       end select

     endif


  end subroutine read_gas_section

  subroutine read_fw_section(keyword, value)

    implicit none
    !logical :: tflag
    character (len=*), dimension(*),intent(in) :: keyword
    character (len=*),intent(in) :: value
    character (len=1023) :: val
    integer pos

    select case (trim(adjustl(keyword(2))))
    case ('isotope_separation')
       read(value,*) useiso
    case ('delnu')
       read(value,*) delnu
    case('lshapemodel')
       read(value,*) lshapemodel
    case('linemixing')
       if (len_trim(keyword(3)).eq.0) then
          read(value, *) use_lm
       else
          select case (trim(adjustl(keyword(3))))
          case('gas')
             val = adjustl(trim(value))
             nr_lmgas  = 0
             pos = index(adjustl(val),' ')
             do
                if (len_trim(val).eq.0) exit
                nr_lmgas = nr_lmgas + 1
                if (pos.gt.0) then
                   read(val(1:pos),*) lm_gas(nr_lmgas)
                else
                   read(val(1:len_trim(adjustl(val))),*) lm_gas(nr_lmgas)
                   exit
                end if
                val = adjustl(val(pos+1:len(val)))
                pos = index(trim(adjustl(val)),' ')
             end do
          case default
             WRITE(16,*) 'BINPUT_PARSE_4_0:READ_FW_SECTION: Parameter ', trim(keyword(3)), 'not defined for fm.linemixing'
             WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_FW_SECTION: Parameter ', trim(keyword(3)), 'not defined for fm.linemixing'
             CALL SHUTDOWN
             STOP 1
          end select
       end if
    case('solar_spectrum')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) ifco
!       else
!          select case (trim(adjustl(keyword(3))))
!          case ('shift')
!             read(value,*) cparm(4)
!          end select
       end if
    case ('pressure_shift')
       read(value,*) fps
    case ('apod_fcn')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_EAPOD
       else
          select case (trim(adjustl(keyword(3))))
          case ('type')
             read(value,*) ieap
          case ('order')
             read(value,*) neap
          end select
       endif
    case ('phase_fcn')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_EPHASE
       else
          select case (trim(adjustl(keyword(3))))
          case ('type')
             read(value,*) iephs
          case ('order')
             read(value,*) nephs
          end select
       endif
    case ('emission')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) emission
          if (emission) then
             iemission = 1
          else
             iemission = 0
          end if
       else
          select case (trim(adjustl(keyword(3))))
          case ('T_infinity')
             read(value,*) emission_t_back
          case ('object')
             read(value,*) emission_object
          case ('normalized')
             read(value,*) emission_norm
             if (emission_norm) then
                ienorm = 1
             else
                ienorm = 0
             end if
          case default
             WRITE(16,*) 'BINPUT_PARSE_4_0:READ_FW_SECTION: Key ', trim(keyword(3)), ' not contained in section fw.emission'
             WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_FW_SECTION: Key ', trim(keyword(3)), ' not contained in section fw.emission'
             CALL SHUTDOWN
             STOP 1
          end select
       end if
    case ('raytonly')
       read(value,*) raytonly
    case default
       WRITE(16,*) 'BINPUT_PARSE_4_0:READ_FW_SECTION: Key ', trim(keyword(2)), ' not contained in section : fw'
       WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_FW_SECTION: Key ', trim(keyword(2)), ' not contained in section : fw'
       CALL SHUTDOWN
       STOP 1
    end select

  end subroutine read_fw_section


  subroutine read_kb_section(keyword, value)
    implicit none

    character (len=*), dimension(*), intent(in) :: keyword
    character (len=*), intent(in) :: value
    !character (len=255) :: tmpstr

    if (len_trim(keyword(2)).eq.0) then
       read(value,*) f_kb
       return
    end if

    select case (trim(adjustl(keyword(2))))
    case('profile')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) f_kb_profile
          return
       end if
       select case (trim(adjustl(keyword(3))))
       case('gas')
          s_kb_profile_gases = adjustl(trim(value))
       case default
          print*, 'Key ', trim(keyword(3)), ' not contained in section kb.profile'
          write(16,*) 'Key ', trim(keyword(3)), ' not contained in section kb.profile'
       end select
    case('temperature')
       read(value,*) f_kb_temp
    case ('slope')
       read(value,*) f_kb_slope
    case ('curvature')
       read(value,*) f_kb_curvature
    case ('solshft')
       read(value,*) f_kb_solshft
    case ('solstrnth')
       read(value,*) f_kb_solstrnth
    case ('phase')
       read(value,*) f_kb_phase
    case ('dwshift')
       read(value,*) f_kb_ifdiff
    case ('wshift')
       read(value,*) f_kb_wshift
    case ('apod_fcn')
       read(value,*) f_kb_eap
    case ('phase_fcn')
       read(value,*) f_kb_ephs
    case ('zshift')
       read(value,*) f_kb_zshift
    case ('sza')
       read(value,*) f_kb_sza
    case ('omega')
       read(value,*) f_kb_fov
    case ('max_opd')
       read(value,*) f_kb_opd
    case ('line')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) f_kb_line
          return
       end if
       select case (trim(adjustl(keyword(3))))
       case ('type')
          read(value,*) i_kb_line_type
       case('gas')
          S_KB_LINE_GASES = adjustl(trim(value))
       case default
          WRITE(16,*) 'BINPUT_PARSE_4_0:READ_KB_SECTION: Key ', trim(keyword(3)), ' not contained in section kb.line'
          WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_KB_SECTION: Key ', trim(keyword(3)), ' not contained in section kb.line'
          CALL SHUTDOWN
          STOP 1
       end select
    case default
       WRITE(16,*) 'BINPUT_PARSE_4_0:READ_KB_SECTION: Key ', trim(keyword(2)), ' not contained in section : kb'
       WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_KB_SECTION: Key ', trim(keyword(2)), ' not contained in section : kb'
       CALL SHUTDOWN
       STOP 1
    end select

  end subroutine read_kb_section

  subroutine read_rt_section(keyword, value)
    implicit none

    character (len=*), dimension(*), intent(in) :: keyword
    character (len=*), intent(in) :: value

    character (len=255) :: tmpstr
    !integer :: nr
    logical :: tflag

    if (len_trim(keyword(2)).eq.0) then
       read(value,*) retflg
       return
    end if

    select case (trim(adjustl(keyword(2))))
    case ('temperature')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) iftemp
       else
          select case (trim(adjustl(keyword(3))))
          case ('sigma')
             read(value,*) tsigma(1:nlayers)
          end select
       end if
    case ('lm')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_LM
       else
          tmpstr = keyword(3)
          select case (trim(adjustl(keyword(3))))
          case('gamma_start')
             read(value,*) gamma_start
             !          case('stop')
             !             read(value,*) stop_criterion
          case('gamma_inc')
             read(value,*) gamma_inc
          case('gamma_dec')
             read(value,*) gamma_dec
          end select
       end if
    case ('convergence')
       read(value,*) convergence
    case ('tolerance')
       read(value,*) tol
    case ('max_iteration')
       read(value,*) itrmax
    case ('wshift')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_WSHIFT
       else
          select case (trim(adjustl(keyword(3))))
          case('type')
             read(value,*) isparm
          case ('apriori')
             read(value, *) wshft
          case ('sigma')
             read(value, *) swshft
          end select
       end if
!    case ('offset')
!       if (len_trim(keyword(3)).eq.0) then
!          read(value,*) tflag
!          if (tflag) nback = 1
!       else
!          iphase = 1
!          select case (trim(adjustl(keyword(3))))
!          case ('apriori')
!             read(value, *) bckoff
!          case ('sigma')
!             read(value, *) sbckoff
!          end select
!       end if
    case ('slope')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_BACKG
          if( F_BACKG )nback = 2
       else
          select case (trim(adjustl(keyword(3))))
          case ('apriori')
             read(value, *) bcksl
          case ('sigma')
             read(value, *) sbcksl
          end select
       end if
    case ('curvature')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) tflag
          if (tflag)then
             F_BACKG = .TRUE.
             nback = 3
          endif
       else
          select case (trim(adjustl(keyword(3))))
          case ('apriori')
             read(value, *) bckcrv
          case ('sigma')
             read(value, *) sbckcrv
          end select
       end if
    case ('solshift')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_RTSOL(4)
       else
          select case (trim(adjustl(keyword(3))))
          case ('apriori')
             read(value,*) ciparm(4)
          case ('sigma')
             read(value,*) scparm(4)
          end select
       end if
    case ('solstrnth')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_RTSOL(5)
       else
          select case (trim(adjustl(keyword(3))))
          case ('apriori')
             read(value,*) ciparm(5)
          case ('sigma')
             read(value,*) scparm(5)
          end select
       end if
    case ('phase')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) ifphase
       else
          select case (trim(adjustl(keyword(3))))
          case ('apriori')
             read(value, *) phs
          case ('sigma')
             read(value, *) sphs
          end select
       end if
    case ('apod_fcn')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_RTAPOD
          !if (tflag) irteap = 1
       else
          select case (trim(adjustl(keyword(3))))
          case ('apriori')
             read(value, *) eappar
          case ('sigma')
             read(value, *) seappar
          end select
       end if
    case ('phase_fcn')
       if (len_trim(keyword(3)).eq.0) then
          read(value,*) F_RTPHASE
          !if (tflag) irtephs = 1
       else
          select case (trim(adjustl(keyword(3))))
          case ('apriori')
             read(value, *) ephspar
          case ('sigma')
             read(value, *) sephspar
          end select
       end if
    case ('ifcalc_se')
       read(value, *) ifcalcse
    case ('dwshift')
       read(value, *) ifdiff
    case default
       WRITE(16,*) 'BINPUT_PARSE_4_0:READ_RT_SECTION: Key ', trim(keyword(2)), ' not contained in section : rt'
       WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_RT_SECTION: Key ', trim(keyword(2)), ' not contained in section : rt'
       CALL SHUTDOWN
       STOP 1
    end select


  end subroutine read_rt_section

  subroutine read_band_section(keyword, value)

    implicit none
    character (len=*), dimension(*),intent(in) :: keyword
    character (len=*), intent(in) :: value

    integer :: nr, nr_band, nr_band_2, nr_beam, pos
    integer :: nr_beams_2
    integer, dimension(maxbnd) :: nbeams
    !integer :: iostat
    character (len=1023)  :: val
    logical :: flag


    if (len_trim(keyword(2)).eq.0) then
       val = adjustl(trim(value))
       nr_band = 0
       pos = index(adjustl(val),' ')
       do
          if (len_trim(val).eq.0) exit
          nband = nband + 1
          !print *, nband
          if (pos.gt.0) then
             read(val(1:pos),*) nbands(nband)
          else
             read(val(1:len_trim(adjustl(val))),*) nbands(nband)
             exit
          end if
          val = adjustl(val(pos+1:len(val)))
          pos = index(trim(adjustl(val)),' ')
       end do
       return
    else
       read(keyword(2),*) nr_band_2
    end if

    flag = .false.

    do nr = 1,nband
       if (nr_band_2.eq.nbands(nr)) then
          flag = .true.
          nr_band = nr
       end if
    end do

    if (.not.flag) return

    select case (trim(adjustl(keyword(3))))
    case ('tempretb')
       read(value,*) tretb (nr_band)
    case ('wave_factor')
       read(value,*) wavfac(nr_band)
    case ('max_opd')
       read(value,*) pmax(nr_band)
    case ('omega')
       read(value,*) omega(nr_band)
    case ('apodization_code')
       read(value,*) iap(nr_band)
    case ('calc_point_space')
       read(value,*) dn(nr_band)
    case ('nu_start')
       read(value, *) wave3(nr_band)
    case ('nu_stop')
       read(value, *) wave4(nr_band)
    case ('snr')
       read(value,*) scnsnr(nr_band,1)
       scnsnr(nr_band,2:maxspe) = scnsnr(nr_band,1)
    case ('gasb')
       val = value
       pos = index(adjustl(trim(val)),' ')
       if (pos.eq.0) write(*,*) 'No gas given in band ', nr_band, '?'
       nretb(nr_band) = 0
       do
          if (len_trim(val).eq.0) exit
          nretb(nr_band) = nretb(nr_band) + 1
          if (pos.gt.0) then
             gasb(nr_band,nretb(nr_band)) = trim(adjustl(val(1:pos)))
          else
             gasb(nr_band,nretb(nr_band)) = val(1:len_trim(adjustl(val)))
             exit
          end if
          val = adjustl(val(pos+1:len(val)))
          pos = index(trim(adjustl(val)),' ')
          !print *, 'pos ', pos
          !print *, 'nretb ', nretb(nr_band)
          !print *, 'nr_band ',nr_band
          !print *, 'gasb ', gasb(nr_band,nretb(nr_band))
       end do
    case ('zshift')
       if (len_trim(keyword(4)).eq.0) then
          read(value, *) f_zshift(nr_band)
       else
          select case (trim(adjustl(keyword(4))))
          case ('type')
             read(value, *) izero(nr_band)
          case ('apriori')
             read(value, *) zshift(nr_band,1)
          case ('sigma')
             read(value, *) szero(nr_band)
          end select
       end if
    case ('beam')
      if (len_trim(keyword(4)).eq.0) then
         nbeam_of_band(nr_band) = 0
         nbeam(nband) = 0
         val = adjustl(trim(value))
         pos = index(adjustl(val),' ')
         do
!rint*, nband, nbeam(nband)
            if (len_trim(val).eq.0) exit
            nbeam(nband) = nbeam(nband) + 1
            if (pos.gt.0) then
               read(val(1:pos),*) nbeams(nbeam(nband))
            else
               read(val(1:len_trim(adjustl(val))),*) nbeams(nbeam(nband))
               exit
            end if
            val = adjustl(val(pos+1:len(val)))
            pos = index(trim(adjustl(val)),' ')
         end do
         return
      else if (trim(adjustl(keyword(4))).eq.'model') then
         val = adjustl(trim(value))
         channel_model_of_band(nr_band) = val(1:2)
      else
         read(keyword(4),*) nr_beams_2
         flag = .false.

         do nr = 1,nbeam(nband)
!rint*,nr, nbeam(nband)
            if (nr_beams_2.eq.nbeams(nr)) then
               flag = .true.
               nr_beam = nr
            end if
         end do

         if (.not.flag) return
         nbeam_of_band(nr_band) = nbeam_of_band(nr_band) + 1
         select case (trim(adjustl(keyword(5))))
         case ('apriori')
            read(value,*) cciparm(nr_band,nr_beam,:)
         case ('sigma')
            read(value,*) schan_scale(nr_band,nr_beam,:)
         end select
      end if
   case default
       WRITE(16,*) 'BINPUT_PARSE_4_0:READ_BAND_SECTION: Key ', trim(keyword(3)), ' not contained in section : band'
       WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_BAND_SECTION: Key ', trim(keyword(3)), ' not contained in section : band'
       CALL SHUTDOWN
       STOP 1
   end select
 end subroutine read_band_section

  subroutine read_spectrum_section(keyword, value)

    implicit none
    character (len=*), dimension(*),intent(in) :: keyword
    character (len=*), intent(in) :: value
    !real(8) :: snr

    character (len=1023)  :: val
    logical :: flag
    integer :: nr_band_2, nr_snr=0, nr, pos
    integer :: nsnr(maxsnr)

    select case (trim(adjustl(keyword(2))))
    case('snr')
       if (len_trim(keyword(3)).eq.0) then
          val = adjustl(trim(value))
          nr_snr = 0
          pos = index(adjustl(val),' ')
          do
             if (len_trim(val).eq.0) exit
             nstnr = nstnr + 1
             !print *, 'nstnr', nstnr
             if (pos.gt.0) then
                read(val(1:pos),*) nsnr(nstnr)
             else
                read(val(1:len_trim(adjustl(val))),*) nsnr(nstnr)
                exit
             end if
             val = adjustl(val(pos+1:len(val)))
             pos = index(trim(adjustl(val)),' ')
          end do
          return
       else
          read(keyword(3),*) nr_band_2
       end if

       flag = .false.
       do nr = 1, nstnr
          if (nr_band_2.eq.nsnr(nr)) then
             flag = .true.
             nr_snr = nr
          end if
       end do

       if (.not.flag) return

       select case (trim(adjustl(keyword(4))))
       case ('nu_start')
          read(value, *) wwv0(nr_snr)
       case ('nu_stop')
          read(value, *) wwv1(nr_snr)
       case ('snr')
          read(value, *) gstnr(nr_snr)
       end select
       case default
          print*, 'Key ', trim(keyword(2)), ' not contained in section : sp'
          write(16,*) 'Key ', trim(keyword(2)), ' not contained in section : sp'
       end select

     end subroutine read_spectrum_section


     subroutine read_output_section(keyword, value)

       implicit none
       character (len=*), dimension(*),intent(in) :: keyword
       character (len=*), intent(in) :: value

       select case (trim(adjustl(keyword(2))))
       case ('level')
          read(value, *) OUTPUTLEVL
       case ('gas_spectra')
          if (len_trim(keyword(3)).eq.0) then
             read(value,*) F_WRTGASSPC
          else
             select case (trim(adjustl(keyword(3))))
             case ('type')
                read(value,*) GASOUTTYPE
             end select
          end if
       case ('k_matrix')
          read(value,*) f_wrtk
       case ('sa_matrix')
          read(value,*) f_wrtsa
       case ('shat_matrix')
          read(value,*) f_wrtshat
       case ('retprofiles')
          read(value,*)  F_WRTRPRF
       case ('aprprofiles')
          read(value,*)  F_WRTAPRF
       case ('ak_matrix')
          read(value,*)  F_WRTAK
       case ('ab_matrix')
          read(value,*)  F_WRTAB
       case ('summary')
          read(value,*)  F_WRTSUMRY
       case ('pbpfile')
          read(value,*)  F_WRTPBP
!       case ('pbpfile_kb')
!          read(value,*)  F_WRTPBP_KB
       case ('channel')
          read(value,*)  F_WRTCHANNEL
       case ('parm_vectors')
          read(value,*)  F_WRTPARM
       case ('seinv_vector')
          read(value,*)  F_WRTSEINV
       case ('sainv_matrix')
          read(value,*)  F_WRTSAINV
       case ('smeas_matrix')
          read(value,*)  F_WRTSMEAS
       case ('ssmooth_matrix')
          read(value,*)  F_WRTSMTH
       case ('raytrace')
          if (len_trim(keyword(3)).eq.0) then
             read(value,*)  F_WRTRAYTC
          else
             select case (trim(adjustl(keyword(3))))
             case ('type')
                read(value,*) RAYOUTTYPE
             end select
          end if

       case ('solarspectrum')
          read(value,*)  F_WRTSOLSPEC
       case ('levmardet')
          read(value,*)  F_WRTLM
       case ('statevec')
          read(value,*)  F_WRTSTV
       case ('xscdetail')
          read(value,*)  XSC_DETAIL
       case ('g_matrix')
          read(value,*)  F_WRTG
       case default
          WRITE(16,*) 'BINPUT_PARSE_4_0:READ_OUTPUT_SECTION: Key ', trim(keyword(2)), ' not contained in section : output'
          WRITE( 0,*) 'BINPUT_PARSE_4_0:READ_OUTPUT_SECTION: Key ', trim(keyword(2)), ' not contained in section : output'
          CALL SHUTDOWN
          STOP 1
       end select

     end subroutine read_output_section
   end module binput_parse_4_0

