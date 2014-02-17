program convert_t15asc_bnr

! Converts t15asc file to a bnr file of toronto style

implicit none

real(4) :: tobs
real(8) :: wlow,whi,spac
INTEGER, PARAMETER :: DOUBLE = SELECTED_REAL_KIND(13)
REAL(DOUBLE)       :: SZA1, ROE1, LAT1, LON1, SECS, BSNR, space
CHARACTER (LEN=80) :: TITLE, file
integer i
INTEGER            :: YYYY, MO, DD, HH, MI
integer(4) :: npfile

write(*,*) 'Filename?'
read(*,*) file

open(unit=15, file=file)
open(unit=16, file=trim(file)//'.bnr', form='unformatted', access='stream',status='unknown')

READ(15, *) SZA1, ROE1, LAT1, LON1, BSNR
READ(15, *) YYYY, MO, DD, HH, MI, SECS
READ(15, '(a80)') TITLE
READ (15, *) WLOW, WHI, SPACE, NPFILE
!06/17/2004 15:09:24UT Z:54.203 A:335.332 D:1443.00 R:0.0035 P:BX F:3.8636mr
! write header, but only formally and the values which are in the file
1 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.2,2(3x,f7.3),3x,f6.4,8x,f6.3)
title = ' '
write(title,1) mo,dd,yyyy,hh,mi,floor(secs),0.0,0.0,0.0,space,45.0
write(16) title
write(16) wlow, whi, space, npfile

DO I = 1, NPFILE
   READ (15, *) TOBS
   write(16) TOBS
end DO

close(15)
close(16)

end program convert_t15asc_bnr
