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

module binput_4_0

! comment out write isoflag - read in from old but not used in new sfit4.ctl

  use params
  use binput_parse_4_0
  use retvparam
  use bandparam
  use lineparam
  use solar
  use channel
  use initialize
  use opt
  use writeout

  implicit none
  save
  integer, parameter :: bp_nr = 10
  logical :: sfix

contains

subroutine read_binput(filename)

  character (len=*), intent(in) :: filename
  character (len=255), dimension(10) :: keyword
  character (len=1023) :: value
  integer :: file_stat, nr_keys, nr
  logical :: bp_exist

  nret = 0

  INQUIRE (FILE=FILENAME, EXIST = BP_EXIST)
  IF (.NOT.BP_EXIST) THEN
     WRITE(16,*) 'BINPUT_4_0:READ_BINPUT: FILE ', TRIM(FILENAME), ' DOES NOT EXIST'
     WRITE( 0,*) 'BINPUT_4_0:READ_BINPUT: FILE ', TRIM(FILENAME), ' DOES NOT EXIST'
     CALL SHUTDOWN
     STOP 1
  END IF

  open(bp_nr, file=filename, status='old', iostat = file_stat)

  do
     call read_line_binput(keyword, nr_keys, value, file_stat)

     if ((file_stat.lt.0).and.(nr_keys.eq.0)) exit

     if (nr_keys.eq.0)then
        cycle
     end if

     select case (trim(adjustl(keyword(1))))
     case ('file')
        call read_file_section(keyword, value)
     case ('gas')
        call read_gas_section(keyword, value)
     case ('fw')
        call read_fw_section(keyword, value)
     case ('kb')
        call read_kb_section(keyword, value)
     case ('rt')
        call read_rt_section(keyword, value)
     case ('band')
        call read_band_section(keyword, value)
     case ('sp')
        call read_spectrum_section(keyword, value)
     case ('out')
        call read_output_section(keyword, value)
     case default
        print *, 'Section ', trim(keyword(1)), ' not defined'
        write(16, *) 'Section ', trim(keyword(1)), ' not defined'
     end select

  end do

  do nr = 1,nband
     if ( mod(nbeam_of_band(nr),2).eq.1) then
        write(*,*) 'Error in beam parameters of band ', nr
        write(16,*) 'Error in beam parameters of band ', nr
     else
        nbeam_of_band(nr) = nbeam_of_band(nr) / 2
     end if
  end do

  close(bp_nr)

end subroutine read_binput

subroutine read_line_binput(keyword, nr_keyword, value, file_stat)

    implicit none

    character (len=*), dimension(:), intent(out) :: keyword
    character (len=*), intent(out) :: value
    integer, intent(out) :: file_stat
    integer, intent(out) :: nr_keyword

    character (len=1023) ::  line
    character (len=1023) :: line_complete
    character (len=255) :: kw
    integer pos,nr,error
    logical flag

    line_complete = ''
    flag = .false.
    file_stat = 0

    nr_keyword = 0

    value = ''


    do
       read (bp_nr, fmt='(a)', iostat=error) line
       if (error.lt.0) then
          file_stat = -1
          goto 1001
       end if
       ! replace tab  characters with blanks
       do nr=1,len_trim(line)
          if (line(nr:nr).eq.achar(9)) then
             line(nr:nr)=' '
          end if
       end do
       line = adjustl(line)
       pos = len_trim(line)

       ! remove comment (everything after #)
       pos=index(line,'#')
       if (pos.gt.0)  line(pos:) = ' '

       if (len(trim(line)).eq.0) then
          goto 1001
       end if
       if (line(1:1).eq.'#') then
          goto 1001
       end if
       if (line(1:1).eq.'.') then
          goto 1001
       end if
       pos=index(line,'=')
       if(flag.and.index(line,'=').gt.0) then
          backspace(bp_nr)
          goto 1001
       end if
       flag = .true.
       line_complete = trim(line_complete)//' '//trim(line)
    end do

1001 continue


    if (len_trim(line_complete).eq.0) then
       nr_keyword = 0
       return
    end if


    pos = index(line_complete, '=')
    value = trim(line_complete(pos+1:len_trim(line_complete)))

!    write(*,*) pos, len_trim(line_complete),line_complete
!    write(*,*) value

    kw = trim(line_complete(1:pos-1))
    pos = index(kw,'.');
    if (pos.eq.0) then
       nr = 0
       goto 1000
    end if
    do nr = 1,1000
       keyword(nr) = trim(kw(1:pos-1))
       kw = kw(pos+1:len(kw))
       pos = index(kw,'.')
       if (pos.eq.0) then
          goto 1000
       end if
    end do



1000 continue
    nr = nr + 1
    keyword(nr) = trim(kw)
    nr_keyword = nr
    keyword(nr+1:) = ''

  end subroutine read_line_binput

end module binput_4_0
