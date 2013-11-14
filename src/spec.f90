module spec

use matrix

implicit none


integer       (4) :: nsnr
real          (8) :: psnr(2,100)
character(len=2)  :: snrid(100)

contains

!---------------------------------------------------------------------------------
!
! ratio's a spectrum taken with the bruker using filters 1,2,3,4,5,6 or 8
!  with a low res spectra which was also normalized with a BBsolar/BB1000˚C
!  calls a 5 point interpolation function taken from Meeuse
!  The ratio file is a bnr and is already open - LUN=10
!
subroutine ratio( outspec, wstart, dnue, np )

integer, intent(inout)                           :: np
real (kind=4), dimension(np),  intent(inout)     :: outspec
real (kind=8), intent(inout)                     :: dnue, wstart

character (len=80)                            :: rtitl
real      (kind=4), dimension(:), allocatable :: rmp, rwv, wavs
real      (kind=8)                            :: rlo, rhi, rspac, y !interp, y
integer   (kind=4)                            :: i, iil, rpts, startpt, endpt
integer   (kind=4), dimension(1)              :: ilow

!  outspec  - pointer to the spectrum on input and output.
!              on output it is the intersection of the 2 spectra
!  wstart   - starting wavenumber of spectrum on input and output and may change
!  dnue     - point spacing same in input and output
!  np       - number of points may change on output

print *,''
print *, 'Ratio data spectra with low-res ratio spectrum.'

! read in data from bnr - already open and rewound

read(9) rtitl
print *, rtitl

read(9) rlo, rhi, rspac, rpts
print *,''
print *, "Ratio spectra:"
write(*,'(1x,a10,f20.10)')'low wave#  ',rlo
write(*,'(1x,a10,2f20.10)')'high wave#  ',rhi
write(*,'(1x,a10,f20.10)')'spacing  ',rspac
write(*,'(1x,a10,i20)')'# points  ',rpts

allocate( rmp( rpts ), rwv(rpts) )
read(9) rmp

!print *,'enter a new rpts:'
!read(5,*) rpts

! create abscissa array of ratio spectrum for interpolation
do i=1, rpts
   rwv(i)= i-1
   end do
rwv = real( rwv * rspac + rlo, 4 )
!print *, rwv(1)

! create abscissa array of data spectrum for interpolation
allocate( wavs( np ))
do i=1, np
   wavs(i)= i-1
   end do
wavs = real( wavs * dnue + wstart, 4 )

print *,''
print *, "Data spectra:"
write(*,'(1x,a10,f20.10)')'low wave#  ',wavs(1)
write(*,'(1x,a10,2f20.10)')'high wave#  ',wavs(np)
write(*,'(1x,a10,f20.10)')'spacing  ',dnue
write(*,'(1x,a10,i20)')'# points  ',np

print *,''

if(( wavs(1) .GT. rhi ) .OR. ( wavs(np) .LT. rlo )) then
   print *,' ratio and data spectra do not overlap...'
   ! we don't change a thing
   deallocate( wavs, rwv, rmp )
   return
   endif

! interpolate envelope spectra and ratio to atm spectra
startpt = 0
endpt   = np
do i= 1, np

   ilow = minloc(( rwv-wavs(i) ) , mask=(rwv-wavs(i)) > 0.0 )

   if( ilow(1) .EQ. 0 ) then
      print *, ' End of ratio spectrum :', i, wavs(i)
      endpt = i -1
      exit
   endif

   if( startpt .eq. 0 ) then
     if( ilow(1) .LT. 4 ) then
        cycle
     else
        print *,' Setting starting wavenumber for ratio : ', i, rwv(ilow(1)-1), wavs(i)
        startpt = i
     endif
   endif

   iil = ilow(1) - 3
   !print *, ilow(1), iil
   !print *, rmp(iil:iil+4), rwv(ilow(1)-1), wavs(i), rwv(ilow(1))
   !print *, rspac
   y = interp( rmp(iil:iil+4), rspac, rwv(ilow(1)-1), wavs(i) )

   ! shift spectrum to the left if the ratio started at a higher wavenumber
   !  than the data (or start == 1)
   outspec(i-startpt+1) = outspec(i) / real( y, 4 )

   !print *, rwv(ilow(1)-1), wavs(i), y, rmp(ilow(1)-1)
   !print *,''

enddo

! reset output starting wavenumber & number of points
wstart = wavs(startpt)
np = np - startpt +1 - (np-endpt)

deallocate( wavs, rwv, rmp )

end subroutine ratio


!---------------------------------------------------------------------------------
subroutine sincinterp ( inspec, outspec, n, wlow, space, opdmax, nterp )

   ! From F. Hase 2002

   integer (4), intent(in)                 :: nterp
   integer (4), intent(inout)              :: n
   real    (8), intent(inout)              :: wlow, space, opdmax
   real    (4), dimension(n),  intent(in)  :: inspec
   real    (4), dimension(:), allocatable, intent(out) :: outspec

   integer (4) :: nofpts_in, nofpts_out, sincradius, nmaxsinc, ninterpol
   integer (4) :: i, j, npos
   real    (8) :: deltanue_out, deltanue_lu, deltanue_in, firstnue_in, firstnue_out
   real    (8) :: fnpos, remain, normsinc, xwert, ywert, faktor, dfirstnue

   real    (8), parameter                 :: pi  = 3.141592653589793d0
   real    (8), parameter                 :: eps = 1.0d-5
   real    (8), dimension(:), allocatable :: sinc

   nmaxsinc  = 50
   ninterpol = nterp

   firstnue_in = wlow
   !print *, ' firstnue_in          : ', firstnue_in

   deltanue_in = space
   !print *, ' deltanue_in          : ', deltanue_in

   nofpts_in   = n
   !print *, ' nofpts_in            : ', nofpts_in

   if( abs(deltanue_in) .ge. 1.00001d0 / (2.0d0 * opdmax) )then
      call warnout('Input spectrum undersampled!...return input spectrum ',0)
      allocate ( outspec( nofpts_in ))
      outspec = inspec
      return
   end if

   write(6,102) 'OPDMAX                      : ', opdmax
   write(6,101) 'Interpolation factor        : ', ninterpol

   deltanue_lu = 0.5d0 / opdmax
   write(6,102)' Minimal spacing             : ', deltanue_lu

   if( ninterpol .gt. 0 )then
      deltanue_out = deltanue_lu / real(ninterpol,8)
   else
      deltanue_out = deltanue_in
      deltanue_lu = deltanue_in
      ninterpol = 1
   endif
   write(6,102)' Final spacing               : ', deltanue_out

   sincradius = nint(nmaxsinc * ninterpol * deltanue_lu / deltanue_in) + 1
   !print *, ' sincradius           : ', sincradius

   firstnue_out = (int((firstnue_in + sincradius * deltanue_in) / deltanue_out) + 2) * deltanue_out
   write(6,102)' Final first wavenumber      : ', firstnue_out

   nofpts_out = int(real(nofpts_in - 2 * sincradius,8) * deltanue_in / deltanue_out) - 4
   write(6,101)' Number of points            : ', nofpts_out

   dfirstnue = firstnue_out - firstnue_in
   !print *, ' dfirstnue            : ', dfirstnue

   !faktor = ninterpol * deltanue_in / deltanue_lu
   faktor = deltanue_in / deltanue_lu
   !print *, ' faktor               : ', faktor

   allocate ( outspec( nofpts_out ))
   outspec = 0.0

   allocate ( sinc( -sincradius : sincradius ))
   sinc = 0.0d0

   ! loop over new output spectrum
   do i = 1, nofpts_out

      fnpos  = 1.0d0 + (dfirstnue + (i-1) * deltanue_out) / deltanue_in
      npos   = nint(fnpos)
      remain = fnpos - real(npos,8)

      !print *, ' fnpos                        : ', fnpos
      !print *, ' npos                         : ', npos
      !print *, ' remain                       : ', remain

      ! generate sinc function
      normsinc = 0.0d0
      do j = -sincradius, sincradius
         xwert = faktor * abs(real(j,8) - remain) + eps
         sinc(j) = sin(pi * xwert) / xwert * cos (0.5d0 * pi * xwert / (faktor * sincradius))
         normsinc = normsinc + sinc(j)
         !write(44,*)  xwert, sinc(j), normsinc
      end do
      !stop

      !print *, ' normsinc  : ', normsinc

      ! interpolation
      ywert = 0.0d0
      do j = -sincradius, sincradius
         ywert = ywert + sinc(j) * real(inspec(npos+j),8)
      end do

      outspec(i) = real( ywert / normsinc, 4 )

   end do

   wlow  = firstnue_out
   space = deltanue_out
   n     = nofpts_out

   deallocate ( sinc )

   return

 101  format( a32, i20 )
 102  format( a32, 2f20.6)
! 103  format( a32, 2f20.12)
! 104  format( a32, 2e20.12 )

! 10   format( a, f10.3, i5 )
 ! 102  format( a32, f11.4, i2, es15.6, i10 )

end subroutine sincinterp


!---------------------------------------------------------------------------------
subroutine calcsnr( wavs, amps, npfile, wlim1, wlim2, spac, opdmax, nterp, noise )

   ! calculate the snr from a small region near the microwindow wanted
   ! use peak signal inmicrowindow
   ! degrade resolution and or interpolate points first so snr is appropriate
   ! to fitted window

   integer   (4), intent(in)   :: npfile, nterp
   real      (8), intent(in)   :: wavs(npfile), spac, opdmax, wlim1, wlim2
   real      (4), intent(in)   :: amps(npfile)
   real      (4), dimension(:), allocatable :: outspec
   real      (8), dimension(:), allocatable :: A, B
   real      (8), dimension(:,:), allocatable :: X, XIT, XINV
   integer   (4)               :: i, j, k, iil, iih, np, order
   integer   (4), dimension(1) :: ilow, ihi
   real      (8)               :: noise, mind, mean, wstart, dnue, opdm, w1, w2, determ



   ! get snr nearest to our mw
   mind = 10000.
   do i=1, nsnr
      noise = abs( (psnr(1,i)+psnr(2,i))/2. - (wlim1+wlim2)/2.)
      if( noise .lt. mind )then
         mind = noise
         k = i
      endif
   enddo

   ! get the spectra in this region +- 1 wavenumber more
   if (nterp.gt.0) then
      ! no resampling
      w1 = psnr(1,k)-3.
      w2 = psnr(2,k)+3.
   else
      w1 = psnr(1,k)
      w2 = psnr(2,k)
   end if
   ilow = minloc(( wavs-w1 ), mask=((wavs-w1) > 0.0D0))
   ihi  = minloc(( wavs-w2 ), mask=((wavs-w2) > 0.0D0))
   iil = ilow(1) - 1
   iih = ihi(1)
   np = iih-iil+1

   write(6,102) '   snr region                 : ', psnr(:,k)
   write(6,102) '   resample region            : ', wavs(iil), wavs(iih)
   write(6,101) '   points in resample calc.   : ', np

   if( np .le. 2 ) then
       write(*,104) '*** calcsnr : Less than 3 points found(1)'
       stop 1
   end if

   ! resample and / or degrade resolution
   ! nterp =  0 - skip resample & resolution degradation
   ! nterp =  1 - minimally sample at opdmax
   ! nterp >  1 - interpolate nterp-1 points

   if( nterp .eq. 0 )then
      allocate( outspec( np ))
      outspec = amps(iil:iih)
   else
      print*,'Resample snr region...'
      opdm = opdmax
      dnue = spac
      wstart = wavs(iil)
      call sincinterp( amps(iil:iih), outspec, np, wstart, dnue, opdm, nterp )

      ! get back the snr sub-region
      iil = 0
      do i=1, np
         if( (i-1)*dnue + wstart .gt. psnr(1,k) )then
            iil = i -1
            exit
         endif
      enddo
      if( iil .eq. 0 )stop 2
      iih = 0
      do i=iil, np
         if( (i-1)*dnue + wstart .gt. psnr(2,k) )then
            iih = i -1
            exit
         endif
      enddo
      if( iih .eq. 0 )stop 3
      np = iih-iil+1
      allocate( outspec( np ))
      outspec = amps(iil:iih)
   endif




! fit <order> order polynomial to spectrum
! see http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html

!   construct matrix X
   order = 4 ! more than 4 does not to be sensible
   allocate(X(np,order), XIT(order,order), XINV(order,order), A(order), B(np))
   do i = 0,order-1
      X(:, i+1) = (/ (j**i, j=1,np) /)
   end do
!     
   mean = sum(outspec(1:np)) / real( np, 8 )

   B(1:np) = outspec(1:np) - mean
   XIT = matmul(transpose(X),X)
   call invrt(XIT, XINV,order)
   A = matmul(matmul(XINV,transpose(X)),B)

   ! remove polynomial from noise window
   B(:) = B(:) - matmul(X,A)
!   do i = 1,np
!      write(10,*)  wavs(iil)+real(i-1)*spac, outspec(i) - mean, B(i) 
!   end do
!  no polynomial fitted to noise window
!   noise = sqrt(dot_product(outspec(1:np)-mean, outspec(1:np)-mean) / real( np, 8 ) )
   noise = sqrt(dot_product(B(:), B(:)) / real( np, 8 ) )
   deallocate(X,XIT,XINV, A, B)

   write(6,101) '   points in snr region     : ', np
   write(6,102) '   mean signal              : ', mean
   write(6,102) '   rms noise                : ', noise
   write(6,102) '   mean snr                 : ', mean/noise

  return

 101  format( a32, i20 )
 102  format( a32, 2f20.12)
! 103  format( a32, 2e20.12 )
 104  format( a48 )


end subroutine calcsnr


!---------------------------------------------------------------------------------
subroutine kpno( opdmax, wl1, wl2, roe, lat, lon, nterp, rflag, flag, loc )

   character (len=80)  :: title
   character (len=1)   :: loc
   real      (4), dimension(:), allocatable :: amps, outspec
   real      (8), dimension(:), allocatable :: wavs
   real      (8) opdmax, wlow, whi, spac, wlim1, wlim2, wl1, wl2, wstart, dnue, roe, noise, peak
   integer   (4) :: npfile, i, iil, iih, np, nterp, rflag, flag
   integer   (4), dimension(1) :: ilow, ihi
   integer       :: yy, mm, dd, hh, nn, ss
   real          :: sza, azm, dur, res, fov
   real      (8) :: lat, lon !tag

   write(*,102) "Desired limits for this window :", wl1, wl2

   flag = 0
   if( opdmax .lt. tiny( 1.0d0 ))then
      print*,' opdmax is 0'
      stop 1
   endif

   ! bnr has been rewound - get title
   read(7) title
   print*,'BNR header:'
   print '(4x,a)', trim(title)

   res = 0.0
   call parsetitle( title, yy, mm, dd, hh, nn, ss, sza, azm, dur, res, fov, loc, roe )
   !print *, yy, mm, dd, hh, nn, ss, sza, azm, dur, res, fov
   !tag = real(dd) + 100.*real(mm) + 10000.*real(yy) + ((real(ss)/60. + real(nn))/60. + real(hh))/24.

   ! get spectra paameters
   read(7) wlow, whi, spac, npfile
   write(*,102) "Spectral params           : ", wlow, whi, spac, npfile
   print*,''

   ! get amplitude data from bnr
   allocate( amps( npfile ), wavs( npfile ))
   read (7, err = 200) amps

   ! calculate wavenumbers
   do i=1, npfile
       wavs(i)= real( i-1, 8 )
   end do
   wavs = wavs*spac + wlow

   ! widen band for interpolation
   wlim1 = wl1 - 10.0
   wlim2 = wl2 + 10.0

   write(*,102) "Limits for interpolation     : ", wlim1, wlim2

   if( wlim1 .gt. whi .or. wlim2 .lt. wlow )then
       write ( *, 801) wlim1, wlim2
       flag = 1
       return
       !close (7)
       !stop 'kpno'
   endif

   ! reset if we don't have the bandwidth
   if( wlim1 .lt. wlow ) then
       write(*,*) " ***Desired interpolation low limit less than bnr file limit."
       wlim1 = wlow
       write(*,102) " Re-setting to file limit", wlim1
   endif

   if( wlim2 .gt. whi ) then
       write(*,*) " ***Desired interpolation high limit greater than bnr file limit."
       wlim2 = whi
       write(*,102) " Re-setting to file limit", wlim2
   endif

   if( wavs(1)-wlim1 .gt. 0.0 ) then
      print *, wavs(1), wlim1, wavs(1)-wlim1
      write(*,104) '***Requested lower bound out of range'
      goto 99
   end if

   if( wavs(npfile)-wlim2 .le. 0.0 ) then
       write(*,104) '***Requested upper bound out of range'
       goto 99
   end if

   ! calculate snr at nearest interval
   print *, ' Calculate snr...'
   call calcsnr( wavs, amps, npfile, wlim1, wlim2, spac, opdmax, nterp, noise )

   ! back to the fit microwindow
   ! resample and / or degrade resolution
   ilow = minloc(( wavs-wlim1 ), mask=((wavs-wlim1) > 0.0D0))
   ihi  = minloc(( wavs-wlim2 ), mask=((wavs-wlim2) > 0.0D0))

   iil = ilow(1) - 1
   iih = ihi(1)

   if( (iih-iil) .le. 2 ) then
       write(*,104) '***Less than 3 points found(1)'
       goto 99
   end if

   np = iih-iil+1
   wstart = wavs(iil)
   dnue = spac

   deallocate( wavs )

   !print *, opdmax, 0.9/res, nterp
   ! nterp =  0 - skip resample & resolution degradation
   ! nterp =  1 - minimally sample at opdmax
   ! nterp >  1 - interpolate nterp-1 points

   if( nterp .eq. 0 )then
      allocate( outspec( np ))
      outspec = amps(iil:iih)
   else
      print*,'Resample fit microwindow...'
      call sincinterp( amps(iil:iih), outspec, np, wstart, dnue, opdmax, nterp )
   endif

   ! ratio if requested
   if( rflag .eq. 1 ) call ratio( outspec, wstart, dnue, np )

   ! test that we still have a valid band
   allocate( wavs( np ))
   do i=1,np
       wavs(i)= real( i-1, 8 )
   end do

   wavs = wavs*dnue + wstart

   ilow = minloc(( wavs-wl1 ), mask=((wavs-wl1) > 0.0D0 ))
   ihi  = minloc(( wavs-wl2 ), mask=((wavs-wl2) > 0.0D0 ))

   iil = ilow(1) - 1
   iih = ihi(1)

   if( (iih-iil) .le. 2 ) then
       write(*,104) '***Less than 3 points found(2)'
       goto 99
   end if

   wlow = wavs(iil)
   whi  = wavs(iih)
   np   = iih-iil+1

   peak = maxval( outspec(iil:iih) )

   print *,''
   print *,"Final spectra:"
   write(*,102)'  low wavenumber             : ',wlow
   write(*,102)'  high wavenuumber           : ',whi
   write(*,102)'  spacing                    : ',dnue
   write(*,101)'  number of points           : ',np
   write(*,102)'  peak signal                : ',peak
   write(*,102)'  RMS SNR                    : ',peak/noise

   print *,''
   print *, ' writing to t15asc'

   ! write to t15asc
   !write(20, 105) tag
   write(20, 106) sza, roe, lat, lon, peak/noise
   write(20, 107) yy, mm, dd, hh, nn, ss
   write(20, 888) title
   write(20, * ) wlow, whi, dnue, np
   write(20, '(e16.6)' ) (outspec(i), i=iil, iih)

   ! write out a bnr while we are at it
   !bfile(1:12) = 'newratio.bnr'
   !print *, ' writing to bnr    : ', bfile
   !open(unit=15, file=bfile, form = 'unformatted')
   !write(15) title
   !write(15) wlow, whi, dnue, np
   !write(15) outspec(iil:iih)
   !close(15)

   deallocate( amps, wavs, outspec )
   flag = 0

   return

 200  write(*,104) '***Error reading spectral values from bnr file'
 99   deallocate( amps, wavs )

   close (7)
   stop 'kpno'

 101  format( a32, i20 )
 102  format( a32, 2f20.12, es15.6, i10 )
 104  format( a48 )
! 105  format( f16.5 )
 106  format( 10f12.5 )
 107  format( i10, 10i5 )
 801  format(5x,'user range ',2f10.4,' not on tape!')
 888  format(a80)

end subroutine kpno


!---------------------------------------------------------------------------------
subroutine parsetitle( title, yy, mm, dd, hh, nn, ss, sza, azm, dur, res, fov, loc, roe )

character (len=80), intent(out)   :: title
character (len=1), intent(inout)  :: loc
integer, intent(out)              :: yy, mm, dd, hh, nn, ss
real, intent(out)                 :: sza, azm, dur, fov, res
real(8)                           :: roe
integer                           :: m = 0

! from bnr.c
!// 1char key for values in header
!
!// Z  - solar zenith angle
!// A  - solar azimuth
!// D  - scan(s) duration
!// F  - filter # (was field of view - jan 2011)
!// V  - field of view
!// K  - observation altitude
!// UT - universal time
!// R  - resolution cm-1
!// P  - apodization function nemonic 2 or 3 chars
!// T  - latitude
!// N  - longitude positive West
!// E  - ROE in direction of azimuth

!06/17/2004 15:09:24UT Z:54.203 A:335.332 D:1443.00 R:0.0035 P:BX F:3.8636mr
read(title,1,err=11)mm, dd, yy, hh, nn, ss, sza, azm, dur, res, fov
goto 10

11 continue
m = m + 1
!04/06/2007 17:29:08UT Z:70.370 A:013.62 D:0721.2 R:0.0035 P:BX F:03.8636mr
read(title,2,err=12)mm, dd, yy, hh, nn, ss, sza, azm, dur, res, fov
goto 10

12 continue
m = m + 1
!03/05/2002 18:59:41UT Z:84.62 DUR:203.33 RES:0.0035 APD:BX FOV:2.2727mr
read(title,3,err=13)mm, dd, yy, hh, nn, ss, sza, dur, res, fov
azm = -999.
goto 10

13 continue
m = m + 1
! ipy
!04/26/2007 17:17:35UT Z:73.980 A:170.208 D:90863.46 R:0.0035 P:BX F:3.8636mr
read(title,4,err=14)mm, dd, yy, hh, nn, ss, sza, dur, res, fov
goto 10

14 continue
m = m + 1
!04/05/2006 17:31:02UT Z:72.390 A:269.548 D:256.14 R:0.0035 P:BX F:1.1962mr
read(title,5,err=15)mm, dd, yy, hh, nn, ss, sza, dur, res, fov
goto 10

15 continue
m = m + 1
!NIWA Lauder BNR format
!L  2013/06/18 17:31:00 SZA=72.39 RES=250.00 FOV=3.45 DUR=345.0 APO=BOX
read(title,6,err=16)loc, yy, mm, dd, hh, nn, ss, sza, res, fov, dur
azm = -999.0
goto 10
!also, required in pspec.f90 open unformatted needs access='stream'
!either as a pre-compiler option or a flag in pspec.input

16 continue
m = m + 1 !6
!sfit4 v0.9.4.1 (ckopus.c)
!20120901 17:08:39UT Z:76.436 A:266.26 D:0204.7 R:0.0035 P:BX V:01.9139 E:6380
read(title,7,err=21)yy, mm, dd, hh, nn, ss, sza, azm, dur, res, fov, roe
goto 10

21 print*, 'spec:parsetitle: header read', m
print*,yy, mm, dd, hh, nn, ss, sza, azm
stop

10 continue
return

1 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.0,2(3x,f7.0),3x,f6.0,8x,f6.0)
2 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.0,3x,f6.0,3x,f7.0,3x,f6.0,8x,f6.0)
3 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f5.0,5x,f6.0,5x,f6.0,12x,f6.0)
4 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.0,3x,f7.0,3x,f8.0,3x,f6.0,8x,f6.0)
5 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.0,3x,f7.0,3x,f6.0,3x,f6.0,8x,f6.0)
6 format(a1,2x,i4,5(1x,i2),5x,f5.1,5x,f6.4,5x,f4.1,5x,f5.1)
7 format(i4,2i2,1x,2(i2,1x),i2,5x,f7.0,3x,f6.0,3x,f6.0,3x,f6.0,8x,f7.0,3x,f4.0)

end subroutine parsetitle


!---------------------------------------------------------------------------------
subroutine warnout(text,predec)

implicit none

character (len=*) :: text
integer,intent(in) :: predec
integer :: intdec

print *,'warning:'
print *, trim(text)
print *,'to shutdown program: enter 0, to go on: enter 1. '
read *, intdec
if (intdec .eq. 0 .or. predec .eq. 0) stop

end subroutine warnout


!---------------------------------------------------------------------------------
!  5-point interpolation scheme from Meeuse
!
real(8) function interp( rmp, rspac, rwv, wav )

implicit none

real(4), dimension(5) :: rmp
real(8) :: rspac, a,b,c,d,e,f,g,h,j,k,n
real(4) :: rwv, wav

!print *, wav, rwv, rspac, rmp
n = (wav-rwv)/rspac

a = rmp(2) - rmp(1)
b = rmp(3) - rmp(2)
c = rmp(4) - rmp(3)
d = rmp(5) - rmp(4)

e = b - a
f = c - b
g = d - c

h = f - e
j = g - f

k = j - h

interp = rmp(3) + n*((b+c)/2. - (h+j)/12.) + n**2*(f/2. - k/24.) + n**3*(h+j)/12. + n**4*k/24.

return

end function interp

end module spec
