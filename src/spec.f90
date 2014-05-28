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

module spec

use fitting

implicit none

integer       (4) :: nsnr
real          (8) :: psnr(2,100)
character(len=2)  :: snrid(100)

integer :: vrsn, blun, rlun, tlun, clun, nlun, olun, ilun, vlun

 ! blun 7  - input bnr file
 ! clun 8  - cinput
 ! ilun 9  - input file pspec.inp
 ! rlun 10 - ratio file
 ! tlun 15 - output bnr file newratio.ispec
 ! olun 16 - simple output for fitbn tag & zero
 ! nlun 20 - t15asc
 ! (bp_nr = 10 sfit4.ctl)

parameter( vrsn=3, blun=7, clun=8, ilun=9, rlun=11, tlun=15, olun=16, nlun=20, vlun=12 )

contains

real(8) function mktag( yy, mm, dd, hh, nn, ss) result( tag )

   integer (4), intent(in) :: yy, mm, dd, hh, nn, ss
   real    (8)             :: dec

   ! make time tag

   !print *, yy, mm, dd, hh, nn, ss, sza, azm, dur, res, fov
   tag = real(dd,8) + 100.d0*real(mm,8) + 10000.d0*real(yy,8)
   dec = ((real(ss,8)/60.d0 + real(nn,8))/60.d0 + real(hh,8))/24.d0
   !print 107, tag, dec
   tag = tag + dec
   !print *, ((real(ss)/60. + real(nn))/60. + real(hh))/24.

   print 107, tag

   return

107 format( "     File time tag : ", f16.7, 2x, f16.7, g27.17, d15.7 )

end function mktag

!******************************************************************************
real(8) function bc4( sp, wv, n, wmid, n10m, old10m, vflag ) result(zero)

! use this method just for 10micron O3 region


! another take completely
! looks at all points near zero
! tries to eliminate those not randomly around 'zero'
! then makes a sraight line through them

! computes straight line through selected (possible) fully absorbed regions
! in the bandpass 750 - 1350 cm-1
! returns the zero offset computed at wavenumber wmid
! returns sp the entire spectra with the baseline corrected

      implicit none

      logical verbose
      integer(4),    intent (in)    :: n, n10m, vflag
      real   (8),    intent (in)    :: wv(:), sp(:)
      real   (8),    intent (in)    :: wmid
      real   (8),    intent (out)   :: old10m(2,n10m)

      !real (8), parameter :: eps = 0.001d0

      integer :: i, j, cnt, cn2, pos, neg, ndif
      real(8) :: avgsw, factor, maxsw, minsw, nstd, pstd, nrstd, prstd !,mdnsw
      real(8), allocatable  :: wavew(:), specw(:), zerod(:), ptwnd(:), newsp(:), lstwv(:), defwv(:), defsp(:)
      real(8), dimension(3) :: curve

      ndif     = 1
      verbose = .false.
      if( vflag .gt. 0 )verbose = .true.

      if( verbose )write(vlun,303) 'bc4: midpoint spectral value : ', wmid
      zero = 0.0d0
! --- quick check that we are in the right region
      if( wmid .lt. 1000. .or. wmid .gt. 1060. )then
         print *," bc4: Spectra out of range - no 10µ baseline corection applied not!"
        ! return
      endif

! --- get 20% of lowest valued points

      factor = 0.25d0
      maxsw  = maxval(sp)
      if( verbose )write(vlun,303) ' initial max spectral value : ', maxsw

      cn2 = 0
      do i=1, n
         if( sp(i) .lt. maxsw )cn2 = cn2 +1
      enddo
      if( verbose )write(vlun,306) ' initial # points in window : ', cn2
      if( cn2 .lt. 1000 ) call exit(5)

      allocate( defwv( cn2 ), defsp( cn2 ))

      cn2 = 0
      do i=1, n
         if( sp(i) .lt. maxsw )then
            cn2 = cn2 +1
            defsp(cn2) = sp(i)
            defwv(cn2) = wv(i)
         endif
      enddo

      cnt = cn2

      do j=1, 20

         if( verbose )write(vlun,*)''
         cn2 = 0
         maxsw = maxval(defsp)
         !print *, maxsw
         maxsw = factor * maxsw
         if( verbose )write(vlun,306) ' iter,  factor, maxsp : ', j, factor, maxsw
         do i=1, cnt
            if( defsp(i) .lt. maxsw )cn2 = cn2 +1
         enddo
         if( verbose )write(vlun,306) ' # points in window : ', cn2
         if( cn2 .lt. 100 )exit

         if( allocated( wavew ))deallocate( wavew )
         if( allocated( specw ))deallocate( specw )
         if( allocated( ptwnd ))deallocate( ptwnd )
         if( allocated( newsp ))deallocate( newsp )
         allocate( newsp( cn2 ))
         allocate( wavew( cn2 ))
         allocate( specw( cn2 ))
         allocate( ptwnd( cn2 ))
         newsp = 0.0d0
         wavew = 0.0d0
         specw = 0.0d0
         ptwnd = 0.0d0

         cn2 = 0
         do i=1, cnt
            if( defsp(i) .lt. maxsw )then
               cn2 = cn2 +1
               specw(cn2) = defsp(i)
               wavew(cn2) = defwv(i)
            endif
         enddo

         if( verbose )write(vlun,303) ' min & max wavenumber : ', minval(wavew), maxval(wavew)
         maxsw = maxval(specw)
         if( verbose )write(vlun,303) ' max spectral value : ', maxsw
         minsw = minval(specw)
         if( verbose )write(vlun,303) ' min spectral value : ', minsw
         avgsw = sum(specw(:)) / real(cn2,8)
         if( verbose )write(vlun,303) ' mean spectral value : ', avgsw

! --- fit line through
         ptwnd      = wavew - wavew(1)
         curve(1:2) = polyfit(ptwnd(:), specw(:), cn2, 1)

         if( allocated( zerod ))deallocate( zerod )
         allocate( zerod( cn2 ))
         if( allocated( lstwv ))deallocate( lstwv )
         allocate( lstwv( cn2 ))
         lstwv = wavew

! --- subtract zero line from spectra
         newsp = curve(1) + curve(2)*ptwnd
         zerod = specw - newsp
         if( verbose )write(vlun,303) ' offset & slope : ', curve(1), curve(2)

!         do i=1, cn2, 100
!            print 307, i, wavew(i), specw(i), newsp(i), zerod(i)
!         enddo

         maxsw = maxval(zerod)
         if( verbose )write(vlun,303) ' max zeroed value : ', maxsw
         minsw = minval(zerod)
         if( verbose )write(vlun,303) ' min zeroed value : ', minsw
         avgsw = sum(zerod(:)) / real(cn2,8)
         if( verbose )write(vlun,303) ' mean zeroed value : ', avgsw

         pos   = 0
         neg   = 0
         nstd  = 0.0d0
         pstd  = 0.0d0
         nrstd = 0.0d0
         prstd = 0.0d0
         do i=1, cn2
            if( zerod(i) .gt. 0.0 )then
               pos = pos +1
               pstd = pstd + zerod(i)
            else
               neg = neg +1
               nstd = nstd + zerod(i)
            endif
         enddo

         pstd = pstd / real(pos,8)
         nstd = nstd / real(neg,8)

         do i=1, cn2
            if( zerod(i) .gt. 0.0 )then
               prstd = prstd + (pstd - zerod(i))**2
            else
               nrstd = nrstd + (nstd - zerod(i))**2
            endif
         enddo

         nstd = sqrt(nrstd)
         pstd = sqrt(prstd)

         if( verbose )write(vlun,306) ' zeroed points greater then 0 : ', pos, pstd
         if( verbose )write(vlun,306) ' zeroed points less then 0 : ', neg, nstd

         if( abs(real( pos - ndif*neg,4)) .le. 0.01*cn2 )exit

         if( abs(real( pos - ndif*neg,4)) .le. 0.03*cn2 )then
            if( pos .gt. ndif*neg ) factor = factor * 1.03d0
            if( pos .lt. ndif*neg ) factor = factor * 0.97d0
         else
            if( pos .gt. ndif*neg ) factor = factor * 1.15d0
            if( pos .lt. ndif*neg ) factor = factor * 0.8d0
         endif

! --- actual zero offset for this band
         zero = curve(1) + curve(2) * (wmid - wavew(1))
         if( verbose )write(vlun,303) 'Zero at midpt : ', zero, wmid

      enddo

      !if( abs( mdnsw - avgsw ) .gt. eps )goto 25

      ptwnd(:)   = wavew(:) - wavew(1)
      curve(1:2) = polyfit(ptwnd(:), specw(:), cn2, 1)

      write(6,302) '10µ Best fit = a + bx :', curve(1), curve(2)

! actual zero offset for this band
      zero = curve(1) + curve(2) * (wmid-wavew(1))
      write(6,303) 'Zero at midpt : ', zero, wmid

 ! generate this line
      old10m(1,1) = wv(1) - (wv(n)-wv(1))/(n10m-1)
      do i=1, n10m
         old10m(1,i) = old10m(1,1) + i * (wv(n)-wv(1))/(n10m-1)
         old10m(2,i) = curve(1) + curve(2) * (old10m(1,i)-old10m(1,1))
      enddo

      if( allocated( wavew )) deallocate( wavew )
      if( allocated( specw )) deallocate( specw )
      if( allocated( ptwnd )) deallocate( ptwnd )
      if( allocated( newsp )) deallocate( newsp )
      if( allocated( zerod )) deallocate( zerod )
      if( allocated( lstwv )) deallocate( lstwv )
      if( allocated( defwv )) deallocate( defwv )
      if( allocated( defsp )) deallocate( defsp )

      return

! 25   print *, ' bc4 no convergence.'
!      call exit(8)

      return

! 10   format(3(a,g12.5))
! 11   format( a, f0.1)
 302  format( /, a32, 6f20.10 )
 303  format( a32, 6f20.10 )
! 305  format( a32, 4f12.9 )
 306  format( a32, i12, 10f20.10 )
! 307  format( i4, f10.4, 3f12.4 )

end function bc4


!---------------------------------------------------------------------------------

real(8) function bc2( sp, wavelength, n, wmid, vflag, noise ) result (zero)


! second improved but similar implmentation to orig
! added more points near 1000 for O3
! changed criteria for saturated window
! still not always working but close

! computes straight line thruough selected (possible) fully absorbed regions
! in the bandpass 750 - 1350 cm-1
! returns the zero offset comupted at wavenumber wmid
! returns sp the entire spectra with the baseline corrected

      implicit none

      logical verbose, blockout

      integer(4),    intent (in)    :: n, vflag
      real   (8),    intent (in)    :: wavelength(:)
      real   (4),    intent (inout) :: sp(:)
      real   (8),    intent (in)    :: wmid
      real   (8),    intent (out)   :: noise

      integer, parameter :: nsat = 41  ! # of saturated regions
      integer, parameter :: n10m = 500 ! # of points in 10 micron region

      integer                    :: i, k, l, mm, count, count2, count3, stdcount, iswap1
      integer                    :: iih, iil, above, below
      integer(4), dimension(1)   :: ihi, ilow, iswap
      integer, dimension(4,nsat) :: inband
      real(8), dimension(2,nsat) :: satarr
      real(8), dimension(nsat)   :: stdarr
      real(8)                    :: temp, initmax, dstncmax, mean
      real(8)                    :: stdev, meansw, runningsum, meanstd, runningmeanstd, distnc, mdwav, mdpnt, azer
      real(8), allocatable       :: wavewindow(:), specwindow(:), zeroed(:), ptwnd(:), newsp(:)
      real(8), allocatable       :: allsatwave(:), allsatspec(:)
      real(8), dimension(2,n10m) :: old10m

      real(8), dimension(3)      :: curve = 0.0d0

      verbose  = .false.
      blockout = .false.
      if( vflag .gt. 0 )verbose  = .true.
      if( vflag .gt. 1 )blockout = .true.
print*,blockout
      dstncmax = 50.0d0
      zero = 0.0d0
! quick check that we are in the right region
      if( wmid .lt. 760. .or. wmid .gt. 1340. )then
         write(6,301)"bc2 : Spectra out of range - no baseline corection applied."
         zero = -999.d0
         return
      endif

     ! the saturated regions - narrow
     satarr(1,1)  = (751.28d0)
     satarr(2,1)  = (751.41d0)
     satarr(1,37) = (752.75d0)
     satarr(2,37) = (752.90d0)
     satarr(1,2)  = (754.20d0)
     satarr(2,2)  = (754.40d0)
     satarr(1,38) = (754.56d0)
     satarr(2,38) = (754.73d0)
     satarr(1,39) = (755.71d0)
     satarr(2,39) = (755.83d0)
     satarr(1,40) = (757.22d0)
     satarr(2,40) = (757.38d0)
     satarr(1,41) = (758.75d0)
     satarr(2,41) = (758.82d0)
     satarr(1,34) = (760.26d0)
     satarr(2,34) = (760.30d0)
     satarr(1,3)  = (776.9300d0)
     satarr(2,3)  = (777.0200d0)
     satarr(1,4)  = (784.4400d0)
     satarr(2,4)  = (784.4800d0)
     satarr(1,5)  = (795.8600d0)
     satarr(2,5)  = (795.9200d0)
     satarr(1,6)  = (798.5300d0)
     satarr(2,6)  = (798.5800d0)
     satarr(1,7)  = (803.5350d0)
     satarr(2,7)  = (803.5650d0)
     satarr(1,35) = (827.6500d0)
     satarr(2,35) = (827.7500d0)
     satarr(1,36) = (839.80d0)
     satarr(2,36) = (839.95d0)
     satarr(1,8)  = (849.5700d0)
     satarr(2,8)  = (849.6000d0)
     satarr(1,9)  = (852.4100d0)
     satarr(2,9)  = (852.4500d0)

     satarr(1,37) = (1001.0d0)
     satarr(2,37) = (1070.0d0)

!     satarr(1,42) = (992.0)
!     satarr(2,42) = (1000.0)

!     satarr(1,37) = (00.0)
!     satarr(2,37) = (00.0)

     satarr(1,10) = (1106.7150d0)
     satarr(2,10) = (1106.7400d0)
     satarr(1,11) = (1135.7100d0)
     satarr(2,11) = (1135.7800d0)
     satarr(1,12) = (1174.4500d0)
     satarr(2,12) = (1174.6000d0)
     satarr(1,13) = (1186.9500d0)
     satarr(2,13) = (1187.0800d0)
     satarr(1,14) = (1198.1400d0)
     satarr(2,14) = (1198.2200d0)
     satarr(1,15) = (1212.1800d0)
     satarr(2,15) = (1212.3200d0)
     satarr(1,16) = (1244.0500d0)
     satarr(2,16) = (1244.2000d0)
     satarr(1,17) = (1260.1000d0)
     satarr(2,17) = (1260.50d0)
     satarr(1,18) = (1268.2800d0)
     satarr(2,18) = (1268.5000d0)
     satarr(1,19) = (1271.40d0)
     satarr(2,19) = (1272.00d0)
     satarr(1,20) = (1287.2000d0)
     satarr(2,20) = (1287.5000d0)
     satarr(1,21) = (1287.7500d0)
     satarr(2,21) = (1287.8500d0)
     satarr(1,22) = (1288.1500d0)
     satarr(2,22) = (1288.5000d0)
     satarr(1,23) = (1296.4500d0)
     satarr(2,23) = (1296.5100d0)
     satarr(1,24) = (1296.6500d0)
     satarr(2,24) = (1296.7700d0)
     satarr(1,25)=(1305.30d0)
     satarr(2,25)=(1306.43d0)
     satarr(1,26)=(1319.0000d0)
     satarr(2,26)=(1320.50d0)
     satarr(1,27)=(1336.7500d0)
     satarr(2,27)=(1337.07d0)
     satarr(1,28)=(1337.55d0)
     satarr(2,28)=(1337.89d0)
     satarr(1,29)=(1340.00d0)
     satarr(2,29)=(1340.70d0)
     satarr(1,30)=(1349.20d0)
     satarr(2,30)=(1349.52d0)
     satarr(1,31)=(1346.97d0)
     satarr(2,31)=(1347.11d0)
     satarr(1,32)=(1347.86d0)
     satarr(2,32)=(1348.09d0)
     satarr(1,33)=(1349.20d0)
     satarr(2,33)=(1349.52d0)

! Now put in increasing wavenumber order
   do i=1, nsat-1
      iswap = minloc(satarr(1,i:nsat))
      iswap1 = iswap(1)+i-1
      if( iswap1 .ne. i )then
         temp = satarr(1,i)
         satarr(1,i) = satarr(1,iswap1)
         satarr(1,iswap1) = temp
         temp = satarr(2,i)
         satarr(2,i) = satarr(2,iswap1)
         satarr(2,iswap1) = temp
      endif
   enddo

! --- start with 'normalized' spectrum
   initmax = real( maxval(sp), 8 )
   sp(:)   = sp(:)/real( initmax, 4 )


! --- pick the saturated regions
! --- first count to find out how big the "allsat" arrays should be
      count2 = 0
      count3 = 0
      above  = 0
      below  = 0
      distnc = 1000.0d0
      l = 0

      do i = 1, nsat

         if (wavelength(1) .lt. satarr(1,i) .and. wavelength(n) .gt. satarr(2,i)) then
            l = l +1
            count3 = count3 +1
            inband(1,l) = i

            ilow = minloc(( wavelength-satarr(1,i) ), mask=(wavelength-satarr(1,i)) > 0.0 )

            ihi  = minloc(( wavelength-satarr(2,i) ), mask=(wavelength-satarr(2,i)) > 0.0 )

            iil  = ilow(1)
            iih  = ihi(1) -1

            inband(2,l) = iih - iil +1
            ! special case for 10 mic - uses bc4
            !if( satarr(1,i) .eq. 1001.0 ) inband(2,l) = n10µ
            inband(3,l) = iil
            inband(4,l) = iih

            meansw = (satarr(2,i) + satarr(1,i))/2.
            if( meansw .lt. wmid )below = below +1
            if( meansw .gt. wmid )above = above +1
            if( abs(meansw-wmid) .lt. distnc )distnc = abs( meansw - wmid)

            count2 = count2 + inband(2,l)
         endif

      enddo

      if( count2 .eq. 0 )go to 667
      if( verbose )write(vlun,*) ''
      if( verbose )write(vlun,304) ' ? sat points in spec : ', count2
      if( verbose )write(vlun,304) ' ? sat bands in spec : ', count3
      if( verbose )write(vlun,304) ' ? sat bands below midpt : ', below
      if( verbose )write(vlun,304) ' ? sat bands above midpt : ', above
      if( verbose )write(vlun,310) ' Closest sat band to midpt : ', distnc

!0test
      if( (above .eq. 0 .or. below .eq. 0) .and. distnc .gt. dstncmax )return

      allocate(allsatwave(count2))
      allocate(allsatspec(count2))

! --- ouput region by region for plotting
      if( blockout )open(unit=98,file='pspec_zero_blocks.dtl')
      if( blockout )write(98,'(a8,i5)') 'nblocks', count3

! --- now loop through each saturated feature and calculate for flatness
      count  = 1
      above  = 0
      below  = 0
      distnc = 1000.0d0

      do l = 1, count3 !through each region

         i   = inband(1,l)    ! satarr index
         k   = inband(2,l)    ! # spectral points in this sat band
         iil = inband(3,l)
         iih = inband(4,l)

! --- reset the wavewindow and specwindow which will be used each time
         if( allocated(wavewindow) ) deallocate( wavewindow )
         if( allocated(specwindow) ) deallocate( specwindow )
         if( allocated(zeroed) ) deallocate( zeroed )
         if( allocated(newsp) ) deallocate( newsp )
         if( allocated(ptwnd) ) deallocate( ptwnd )

         allocate(wavewindow(k))
         allocate(specwindow(k))
         allocate(zeroed(k))
         allocate(ptwnd(k), newsp(k))

         wavewindow(:) = wavelength(iil:iih)
         specwindow(:) = real(sp(iil:iih),8)
         mdwav         = (wavelength(iil) + wavelength(iih)) /2.0
         mdpnt         = mdwav - wavewindow(1)

         if( satarr(1,i) .ne. 1001.0d0 .and. satarr(1,i) .ne. 992.0d0 )then

            ptwnd(:)      = wavewindow(:) - wavewindow(1)

            ! --- fit a second degree polynomial to the region
            ! --- curve = a + bx + cx^2
            curve(1:3) = polyfit( ptwnd(:), specwindow(:), k, 2 )
!            print*, polyfit( ptwnd(:), specwindow(:), k, 2 )


            ! --- zero offset at midpoint
            zero = curve(1) + curve(2)*mdpnt + curve(3)*mdpnt*mdpnt

            ! --- subtract zero line from spectra
            newsp  = curve(1) + curve(2)*ptwnd + curve(3)*ptwnd**2
            zeroed =  specwindow - newsp

            ! --- calculate the standard deviation in that region too

            ! --- average value of the zeroed spectrum
            azer = sum(zeroed(:)) / max(1.0d0,real(k,8))

            ! --- calculate the standard deviation, intermediate step required here
            runningsum = 0.0d0
            do mm = 1, k
               runningsum=runningsum+((zeroed(mm)-azer)**2)
            enddo
            stdev = sqrt(runningsum/real(k,8))

            ! --- for plotting blocks
            if( blockout )then
               write(98,'(a8,2i5,10f14.4)') 'block', l, k, curve(1:3), zero, stdev, wavelength(iih)-wavelength(iil), mdwav
               do i=1, k
                  write(98,307) i, ptwnd(i), wavewindow(i), specwindow(i), newsp(i), zeroed(i)
               enddo
            endif
            if( verbose )write(vlun,306) 'npts mdpt mzer azer std off slp crv : ', k, mdwav, zero, azer, stdev, curve(1:3)

            ! --- if they pass the polynomial fit, add both the region and the std dev to arrays
            if ( curve(3)/(wavelength(iih)-wavelength(iil)) .lt. 50. .and. zero .lt. 0.25d0 .and. stdev .lt. 0.05 )then

               allsatwave(count:count+k-1)=wavewindow(:)
               allsatspec(count:count+k-1)=specwindow(:)
               stdarr(l) = stdev
               count     = count + k

               if( mdwav .lt. wmid )below = below +1
               if( mdwav .gt. wmid )above = above +1
               if( abs(mdwav-wmid) .lt. distnc )distnc = abs( mdwav - wmid)

               if( verbose )write(vlun,306) 'good region : npts mdpt mzer stdv : ', k, mdwav, zero, stdev

            endif

         else
            ! --- 60cm-1 around 10 micron O3 treat this as a large band & try to get 50 (n10µ)
            !     zeroed data points from it
            zero = bc4( specwindow, wavewindow, k, mdwav, n10m, old10m, vflag )

            allsatwave(count:count+n10m-1) = old10m(1,:)
            allsatspec(count:count+n10m-1) = old10m(2,:)
            stdarr(l) = 0.05d0 ! did not calc yet   stdev
            count     = count + n10m

            if( mdwav .lt. wmid )below = below +1
            if( mdwav .gt. wmid )above = above +1
            if( abs(mdwav-wmid) .lt. distnc )distnc = abs( mdwav - wmid)

            if( verbose )write(vlun,306) 'good region : npts mdpt mzer stdv : ', n10m, mdwav, zero, stdev

         endif

      end do ! finished that region
      count = count - 1

      if( blockout )close(98)
      if( verbose )write(vlun,304) ' Found sat bands below midpt : ', below
      if( verbose )write(vlun,304) ' Found sat bands above midpt : ', above
      if( verbose )write(vlun,310) ' Closest sat band to midpt : ', distnc
      if( verbose )write(vlun,304) " Points in allsat vector : ", count

! --- calculate a 2nd order poly using all points
     curve(1:3) = polyfit(allsatwave(1:count), allsatspec(1:count), count, 2)

! --- calculate the mean std
      runningmeanstd = 0.0d0
      stdcount       = 0
      do i = 1, nsat
       if( abs(stdarr(i)) .gt. tiny( 0.0d0 ))then
         runningmeanstd = runningmeanstd + stdarr(i)
         stdcount       = stdcount + 1
       endif
      enddo
      meanstd = runningmeanstd/stdcount

      write(6,303) 'Best fit= a + bx + cx^2, std : ', curve(1), curve(2), curve(3), meanstd

! --- calculate the rms noise of this zero vector
      noise = 0.0d0
      !do i=1, count
      !   noise = noise + &
      !     ( allsatspec(i) - ( curve(1) + ( curve(2) + curve(3)*allsatwave(i) )*allsatwave(i)) ) * &
      !     ( allsatspec(i) - ( curve(1) + ( curve(2) + curve(3)*allsatwave(i) )*allsatwave(i)) )
      !   write(66,*) allsatwave(i), allsatspec(i), (curve(1) + (curve(2) + curve(3)*allsatwave(i) ) * allsatwave(i)), &
      !               allsatspec(i) - (curve(1) + (curve(2) + curve(3)*allsatwave(i) ) * allsatwave(i))
      !enddo
      !noise = sqrt( noise / real( count, 8 ))

! --- assume horizontal band - these points are random around zero by definition if we got this far
      mean = sum(allsatspec(1:count)) / real( count, 8 )
      noise = sqrt(dot_product(allsatspec(1:count)-mean, allsatspec(1:count)-mean) / real( count, 8 ) )

! --- actual zero offset for this band
      zero = curve(1) + curve(2)*wmid + curve(3)*wmid*wmid

! --- don't go too far afield !!!
! --- 40 gets o3 at 1002 and mid pt of 10µ region ~1025
      if(( above .eq. 0 .or. below .eq. 0 ) .and. distnc .gt. dstncmax )then
         write(6,302) 'Zero found for this region : ', zero
         zero = 0.0d0
         sp(:) = real(initmax,4) * sp(:)
         print *, 'No zero offset applied...return now.'
      else
         print *, 'Zero correcting this spectrum.'
! --- send back the entire spectrum zero level adjusted
!     and scaled back to original peak value
         do i=1, n
           temp  = real(wavelength(i),8)
           sp(i) = real(initmax,4) * (sp(i) - real(curve(1) + ( curve(2)+ curve(3)*temp )*temp ,4))
         enddo
      endif

      if( allocated( wavewindow) ) deallocate( wavewindow )
      if( allocated( specwindow) ) deallocate( specwindow )
      if( allocated( zeroed) )     deallocate( zeroed )
      if( allocated( newsp) )      deallocate( newsp )
      if( allocated( ptwnd) )      deallocate( ptwnd )
      if( allocated( allsatwave) ) deallocate( allsatwave )
      if( allocated( allsatspec) ) deallocate( allsatspec )

      return

667   continue
      print *, ' None of the saturated regions are in this file'
      print *, ' Spectrum is probably not filter 6...'
      zero = 0.0d0

      return

! 10   format(3(a,g12.5))
 301  format( /, a )
 302  format( a32, 6f12.5 )
 303  format( a32, 6f20.10 )
 304  format( a32, 6i8 )
! 305  format( a32, 4f12.9 )
 306  format( a32, i12, 10f12.4 )
 307  format( i4, f10.4, 5f12.5 )
 310   format( a32, f12.1)

end function bc2


!---------------------------------------------------------------------------------
!
! ratio's a spectrum taken with the bruker using filters 1,2,3,4,5,6 or 8
!  with a low res spectra which was also normalized with a BBsolar/BB1000˚C
!  calls a 5 point interpolation function taken from Meeuse
!  The ratio file is a bnr and is already open - LUN=10
!
subroutine ratio( outspec, wstart, dnue, np )

   integer, intent(inout)                         :: np
   real (kind=4), dimension(np), intent(inout)    :: outspec
   real (kind=8), intent(inout)                   :: dnue, wstart

   character (len=80)                             :: rtitl
   real      (kind=4), dimension(:), allocatable  :: rmp, rwv, wavs
   real      (kind=8)                             :: rlo, rhi, rspac, y !interp, y
   integer   (kind=4)                             :: i, iil, rpts, startpt, endpt
   integer   (kind=4), dimension(1)               :: ilow

   !  outspec  - pointer to the spectrum on input and output.
   !              on output it is the intersection of the 2 spectra
   !  wstart   - starting wavenumber of spectrum on input and output and may change
   !  dnue     - point spacing same in input and output
   !  np       - number of points may change on output

   write(6,111) 'Ratio data spectra with low-res ratio spectrum.'

   ! read in data from bnr - already open and rewound

   read(rlun) rtitl
   write(6,112) rtitl

   read(rlun) rlo, rhi, rspac, rpts
   write(6,102)'Low wavenumber : ', rlo
   write(6,102)'High wavenumber : ', rhi
   write(6,102)'Spacing : ', rspac
   write(6,101)'Number of points : ', rpts

   allocate( rmp( rpts ), rwv(rpts) )
   read(rlun) rmp

   !print *,'enter a new rpts:'
   !read(5,*) rpts

   ! create abscissa array of ratio spectrum for interpolation
   do i=1, rpts
      rwv(i)= real(i-1,4)
      end do
   rwv = real( rwv * rspac + rlo, 4 )
   !print *, rwv(1)

   ! create abscissa array of data spectrum for interpolation
   allocate( wavs( np ))
   do i=1, np
      wavs(i)= real(i-1,4)
      end do
   wavs = real( wavs * dnue + wstart, 4 )

   write(6,111) "Data spectra :"
   write(6,102)'low wavenumber : ', wavs(1)
   write(6,102)'high wavenumber : ', wavs(np)
   write(6,102)'spacing : ', dnue
   write(6,101)'# points : ', np
   write(6,*)''

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
           write(6,102)'Starting wavenumber for ratio : ', rwv(ilow(1)-1), wavs(i)
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
   wstart = real(wavs(startpt),8)
   np = np - startpt +1 - (np-endpt)

   deallocate( wavs, rwv, rmp )

return

 101  format( a32, i20 )
 102  format( a32, 2f20.12, f20.12, i10 )
 111  format( /, a )
 112  format( 4x, a )

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
      call warnout('Input spectrum undersampled!...return input spectrum ')
      allocate ( outspec( nofpts_in ))
      outspec = inspec
      return
   end if

   write(6,102) 'OPDMAX : ', opdmax
   write(6,101) 'Interpolation factor : ', ninterpol

   deltanue_lu = 0.5d0 / opdmax
   write(6,102)'Minimal spacing : ', deltanue_lu

   if( ninterpol .gt. 0 )then
      deltanue_out = deltanue_lu / real(ninterpol,8)
   else
      deltanue_out = deltanue_in
      deltanue_lu = deltanue_in
      ninterpol = 1
   endif
   write(6,102)'Final spacing : ', deltanue_out

   sincradius = nint(nmaxsinc * ninterpol * deltanue_lu / deltanue_in) + 1
   !print *, ' sincradius           : ', sincradius

   firstnue_out = real((int((firstnue_in + sincradius * deltanue_in) / deltanue_out) + 2), 8) * deltanue_out
   write(6,102)'Final first wavenumber : ', firstnue_out

   nofpts_out = int(real(nofpts_in - 2 * sincradius,8) * deltanue_in / deltanue_out) - 4
   write(6,101)'Number of points : ', nofpts_out

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

      !print *, ' normsinc  : ', i, normsinc

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
! 105  format( a, f10.3, i5 )
! 106  format( a32, f11.4, i2, es15.6, i10 )

end subroutine sincinterp


!---------------------------------------------------------------------------------
subroutine kpno( opdmax, wl1, wl2, roe, lat, lon, nterp, rflag, oflag, zflag, vflag, bflag, loc )

   character (len=80)  :: title, bfile
   character (len=1)   :: loc
   real      (4), dimension(:), allocatable :: amps, outspec
   real      (8), dimension(:), allocatable :: wavs
   integer   (4), dimension(1)              :: ilow, ihi
   real      (8) :: opdmax, wlow, whi, spac, wlim1, wlim2, wl1, wl2, wstart, dnue, roe, noise, peak
   logical       :: writezero
   integer   (4) :: npfile, i, iil, iih, np, nterp, rflag, bflag, oflag, vflag
   integer       :: yy, mm, dd, hh, nn, ss
   real      (4) :: sza, azm, dur, res, fov
   real      (8) :: lat, lon, pspc, tag, wmid, zero, zflag

   writezero = .false.
   if( vflag .gt. 0 )writezero = .true.

   write(6,101) 'Interpolation flag : ', nterp
   write(6,101) 'Ratio flag : ', rflag
   write(6,102) 'Zero flag : ', zflag

   noise = 0.0d0
   zero  = 0.0d0

   ! wl1   - wanted low wv #
   ! wl2   - wanted high wv #
   ! wmid  - mid point in this region region

   ! wlow  - low wv # in input file
   ! whi   - high wv # in input file

   ! wlim1 - extended low wv # for interpolation & ratio
   ! wlim2 - extended high wv # for interpolation & ratio

   wmid = (wl2+wl1)/2.
   write(6,102) "Limits for this window : ", wl1, wl2, wmid

   bflag = 0
   if( opdmax .lt. tiny( 0.0d0 ))then
      print*,' opdmax is 0. ?!'
      stop '1'
   endif

   ! bnr has been rewound - get title
   read(blun) title
   write(6,111) 'BNR header:'
   write(6,112) trim(title)

   res = 0.0
   call parsetitle( title, yy, mm, dd, hh, nn, ss, sza, azm, dur, res, fov, loc, roe )
   !print *, yy, mm, dd, hh, nn, ss, sza, azm, dur, res, fov

   ! get spectra parameters
   read(blun) wlow, whi, spac, npfile
   write(6,102)'Low wavenumber : ', wlow
   write(6,102)'High wavenumber : ', whi
   write(6,102)'Spacing : ', spac
   write(6,101)'Number of points : ', npfile

   ! get amplitude data from bnr
   if( allocated( amps ) )deallocate( amps )
   if( allocated( wavs ) )deallocate( wavs )
   allocate( amps( npfile ), wavs( npfile ))
   read (blun, err = 200) amps

   ! calculate wavenumbers
   do i=1, npfile
       wavs(i)= real( i-1, 8 )
   end do
   wavs = wavs*spac + wlow



! --- Step 1 : zero correct the whole spectra
   ! return spectrum (amps) in place
   ! return value zero is for the mid point of the wanted window
   ! zflag  :
   !     = 0 no zero offset,
   !     = 0 < zflag < 1, use this value for baselincorrect,
   !     = 2 use combo 2 + 4 for 10µ region
   ! noise is the rms noise from the zero (allsatspec) vector - by definition random around zoro

   if( idint(zflag) .eq. 0 )then
      zero  = 0.0d0
      noise = 0.0d0

   else if( abs(zflag) .gt. 0. .and. abs(zflag) .lt. 1. )then
      zero    = zflag
      noise   = 0.0d0
      amps(:) = amps(:) - real(zero,4)

   else if( idint(zflag) .gt. 1 .and. idint(zflag) .lt. 3 )then

      if( idint(zflag) .eq. 2 ) zero = bc2( amps, wavs, npfile, wmid, vflag, noise )

      if( zero .ne. -999.d0 )then
         write(6, 102) 'Zero offset determined : ', zero
         write(6, 102) ' at wavenumber : ', wmid
         write(6, 102) 'RMS Noise from zero : ', noise

         if( writezero )then
            open( unit=33, file='zeroed.bnr',status='unknown', form='unformatted' )
            write(33) title
            write(33) wlow, whi, spac, npfile
            write(33) amps  !/maxval(amps)
            close(33)
         endif
      endif
   else
      print *, "!! zflag not a valid option : ", zflag
      print *, "Valid zflag options: 0, 0 < zflag < 1, 2"
      call exit(2)

   endif



! --- Step 2 : Extend micro window region for interpolation
   ! Choose a region around the desired band for ratio and interpolation
   ! take any interpolation into account
   pspc = 50.*nterp/opdmax
   write(6,*)''
   if( vflag .gt. 0 )write(6,102) 'Interpolation extension : ', pspc
   wlim1 = real(wl1 - pspc,8)
   wlim2 = real(wl2 + pspc,8)

   write(6,102) "Extended limits : ", wlim1, wlim2
   write(6,102) "File limits : ", wlow, whi

   if( wlim1 .gt. whi .or. wlim2 .lt. wlow )then
       write ( *, 801) wlim1, wlim2
       bflag = 1
       return
   endif

   ! reset if we don't have the bandwidth
   if( wlim1 .lt. wlow ) then
       write(*,*) "***Desired interpolation low limit less than bnr file limit."
       wlim1 = wlow
       write(*,102) "Re-setting to file limit", wlim1
   endif

   if( wlim2 .gt. whi ) then
       write(*,*) "***Desired interpolation high limit greater than bnr file limit."
       wlim2 = whi
       write(*,102) "Re-setting to file limit", wlim2
   endif

   if( wavs(1)-wlim1 .gt. 0.0 ) then
      print *, wavs(1), wlim1, wavs(1)-wlim1
      write(*,104) '***Requested lower bound out of range'
      goto 99
   end if

   if( wavs(npfile)-wlim2 .lt. 0.0 ) then
       write(*,104) '***Requested upper bound out of range'
       goto 99
   end if



! --- Step 3 : Calculate SNR
   ! calculate snr at nearest interval
   write(6,111) 'Calculate noise...'
   !noise=0.0 !0test
   call calcsnr( wavs, amps, npfile, wlim1, wlim2, spac, opdmax, nterp, noise, vflag )



! --- Step 4 : Interpolate if requested
   ! back to the fit microwindow
   ! resample and / or degrade resolution
   ilow = minloc(( wavs-wlim1 ), mask=((wavs-wlim1) > 0.0D0))
   ihi  = minloc(( wavs-wlim2 ), mask=((wavs-wlim2) > 0.0D0))

   iil = ilow(1) - 1
   iih = ihi(1)
   if( vflag .gt. 0 )write(6,109) 'Spectra segment before ratio : ',iil, wavs(iil), iih, wavs(iih), iih-iil

   if( (iih-iil) .le. 2 ) then
       write(*,104) '***Less than 3 points found(1)'
       goto 99
   end if

   np     = iih-iil+1
   wstart = wavs(iil)
   dnue   = spac
   deallocate( wavs )

   ! nterp =  0 - skip resample & resolution degradation
   ! nterp =  1 - minimally sample at opdmax
   ! nterp >  1 - interpolate nterp-1 points at opdmax
   if( nterp .eq. 0 )then
      allocate( outspec( np ))
      outspec = amps(iil:iih)
   else
      write(6,111) 'Resample fit microwindow...'
      if( vflag .gt. 0 )write(6,*)"to interp : ", iil, iih, np, wstart, dnue
      call sincinterp( amps(iil:iih), outspec, np, wstart, dnue, opdmax, nterp )
      if( vflag .gt. 0 )write(6,*)'from interp : ',np, dnue, wstart
   endif

   peak = real( maxval( outspec(:) ), 8 )

! --- Step 5: Ratio if requested
   if( rflag .eq. 1 ) call ratio( outspec, wstart, dnue, np )


   ! test that we still have a valid band
   allocate( wavs( np ))
   do i=1,np
       wavs(i)= real( i-1, 8 )
   end do
   wavs = wavs*dnue + wstart

   ilow = minloc(( wavs-wl1 ), mask=((wavs-wl1) > 0.0D0 ))
   ihi  = minloc(( wavs-wl2 ), mask=((wavs-wl2) > 0.0D0 ))

   iil = ilow(1)
   if(iil .ne. 1 )iil = ilow(1) - 1
   iih = ihi(1)

   if( (iih-iil) .le. 2 ) then
       write(*,104) '***Less than 3 points found(2)'
       goto 99
   end if

   !print*, iil, iih
   !print*,wstart, wavs(iil), wavs(iih), wl1

   wlow = wavs(iil)
   whi  = wavs(iih)
   np   = iih-iil+1


! --- write out final spectrum
   print *,''
   write(6,111)'Final spectra : '
   write(*,102)'Low wavenumber : ',wlow
   write(*,102)'High wavenuumber : ',whi
   write(*,102)'Spacing : ',dnue
   write(*,101)'Number of points : ',np
   write(*,102)'Peak signal : ',peak
   write(*,102)'RMS SNR in fit region : ',peak/noise

   print *,''

   if( oflag .eq. 1 .or. oflag .eq. 3 )then

      ! write to t15asc
      write(6, 103) 'Writing to : ', 't15asc.4'
      write(tlun, 106) sza, roe, lat, lon, peak/noise
      write(tlun, 107) yy, mm, dd, hh, nn, ss
      write(tlun, 888) title
      write(tlun, 108) wlow, whi, dnue, np
      !if( oflag .eq. 3 ) outspec(iil:iih) = outspec(iil:iih)/maxval(outspec(iil:iih))
      write(tlun, '(e16.6)' ) (outspec(i), i=iil, iih)
   endif

   if( oflag .eq. 2 .or. oflag .eq. 3 )then

      ! write out a bnr while we are at it
      write(bfile, 105) tag, ".bnr"
      write(bfile, 110) "t15.bnr"
      write(6, 103) 'Writing to bnr : ', trim(bfile)
      open(unit=nlun, file=trim(bfile), form = 'unformatted')
      write(nlun) title
      write(nlun) wlow, whi, dnue, np
      !if( oflag .eq. 4 ) outspec(iil:iih) = outspec(iil:iih)/maxval(outspec(iil:iih))
      write(nlun) outspec(iil:iih)
      close(nlun)

   endif

   deallocate( amps, wavs, outspec )
   !flag = 0

   return

 200  write(*,104) '***Error reading spectral values from bnr file'
 99   deallocate( amps, wavs )

   close (7)
   stop 'kpno'

 101  format( a32, i20 )
 102  format( a32, 2f20.12, f20.12, i10 )
 104  format( a48 )
 103  format( a32, a )
 105  format( f16.5, a4 )
 106  format( 10f12.5 )
 107  format( i10, 10i5 )
 108  format( 3f22.12, i15 )
 109  format( /, a32, 3(i20,f20.12))
 110  format( a7 )
 111  format( /, a )
 112  format( 4x, a )
 !113  format( a )
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
real(4)                           :: opd
character (len=3)                 :: apd
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
!// O  - optical path difference

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

16 continue
m = m + 1 !6
!sfit4 v0.9.4.1 (ckopus.c)
!20120901 17:08:39UT Z:76.436 A:266.26 D:0204.7 R:0.0035 P:BX V:01.9139 E:6380
read(title,7,err=17)yy, mm, dd, hh, nn, ss, sza, azm, dur, res, fov, roe
goto 10

17 continue
m = m + 1
!BIRA type of spectra headers
!20110125 04:04:55  102s ZT=04:05:45 OPD= 81.97 FOV= 6.36 APF=BX aS 61.553  28.45
!ZT is the median of the recording time... and was used to get the ZA
!print *,trim(title)
read(title,8,err=18)yy, mm, dd , dur, hh, nn, ss, fov, sza
!print *,yy,mm,dd,hh,nn,ss,dur,fov,sza
goto 10

18 continue
m = m + 1
! zephyr2
!20090602 09:54:13  205s ZT=09:55:55 OPD=257.14 FOV= 2.75 APF=BX aS 60.531 999.99
read(title,9,err=21)yy, mm, dd, dur, hh, nn, ss, fov, sza, azm
goto 10

21 print*, 'spec:parsetitle: header read', m
print*,yy, mm, dd, hh, nn, ss, sza, azm
stop '4'

10 continue
return

1 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.0,2(3x,f7.0),3x,f6.0,8x,f6.0)
2 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.0,3x,f6.0,3x,f7.0,3x,f6.0,8x,f6.0)
3 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f5.0,5x,f6.0,5x,f6.0,12x,f6.0)
4 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.0,3x,f7.0,3x,f8.0,3x,f6.0,8x,f6.0)
5 format(2(i2,1x),i4,1x,2(i2,1x),i2,5x,f6.0,3x,f7.0,3x,f6.0,3x,f6.0,8x,f6.0)
6 format(a1,2x,i4,5(1x,i2),5x,f5.1,5x,f6.4,5x,f4.1,5x,f5.1)
7 format(i4,2i2,1x,2(i2,1x),i2,5x,f7.0,3x,f6.0,3x,f6.0,3x,f6.0,8x,f7.0,3x,f4.0)
8 format(i4,2(i2),11x,f3.0,5x,3(i2,1x),15x,f6.0,10x,f6.0)
9 format(i4,i2,i2,11x,f3.0,5x,3(i2,1x),16x,f4.0,11x,f6.0,1x,f6.0)

end subroutine parsetitle


!---------------------------------------------------------------------------------
subroutine warnout(text)

   implicit none

   character (len=*) :: text
   integer :: intdec

   print *,'warning:'
   print *, trim(text)
   print *,'to shutdown program: enter 0, to go on: enter 1. '
   read *, intdec
   !if (intdec .eq. 0 .or. predec .eq. 0) stop
   if (intdec .eq. 0 ) call exit (6)

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

   a = real( rmp(2) - rmp(1), 8)
   b = real( rmp(3) - rmp(2), 8)
   c = real( rmp(4) - rmp(3), 8)
   d = real( rmp(5) - rmp(4), 8)

   e = b - a
   f = c - b
   g = d - c

   h = f - e
   j = g - f

   k = j - h

   interp = rmp(3) + n*((b+c)/2. - (h+j)/12.) + n**2*(f/2. - k/24.) + n**3*(h+j)/12. + n**4*k/24.

return

end function interp

!---------------------------------------------------------------------------------
subroutine calcsnr( wavs, amps, npfile, wlim1, wlim2, spac, opdmax, nterp, noise, vflag )

   ! calculate the snr from a small region near the microwindow wanted
   ! use peak signal inmicrowindow
   ! degrade resolution and or interpolate points first so snr is appropriate
   ! to fitted window

   integer   (4), intent(in)    :: npfile, nterp, vflag
   real      (8), intent(in)    :: wavs(npfile), spac, opdmax, wlim1, wlim2
   real      (4), intent(in)    :: amps(npfile)
   real      (8), intent(inout) :: noise
   real      (8), dimension(:), allocatable :: x, y, z, curve
   real      (4), dimension(:), allocatable :: outspec
   integer   (4)                :: i, k, iil, iih, np
   integer   (4), dimension(1)  :: ilow, ihi
   real      (8)                :: mind, mean, wstart, dnue, opdm, w1, w2

!print*, wavs

   ! if we already calculated noise in the 10µ region adjust for reduced resolution if needed
   ! noise from zero is calculated with the initial spectrum
   if( noise .gt. tiny( 0.0d0 ))then
      if( nterp .eq. 0 )return
      opdm = 0.5d0 / spac
      if( opdmax .lt. opdm )then
         noise = noise * sqrt( opdmax / opdm ) * real(nterp,8)
         write(6,102) 'Noise (10µ) after resample : ', noise
         return
      endif
   endif

   ! get snr nearest to our mw
   k    = 0
   mind = 10000.d0
   do i=1, nsnr
      !print*, i, psnr(:,i), wlim1, wlim2
      noise = abs( (psnr(1,i)+psnr(2,i))/2. - (wlim1+wlim2)/2.)
      if( noise .lt. mind )then
         mind = noise
         k = i
      endif
   enddo

  if( vflag .gt. 1 )then
      open(66,file='noisefit.txt')
      write(66,*)'Nearest exact noise region in raw spectrum'
      w1 = psnr(1,k)
      w2 = psnr(2,k)
      ilow = minloc(( wavs-w1 ), mask=((wavs-w1) > 0.0D0))
      ihi  = minloc(( wavs-w2 ), mask=((wavs-w2) > 0.0D0))
      iil = ilow(1) - 1
      iih = ihi(1)
      write(66,*) 1, iih - iil + 1
      do i=iil, iih
         write(66,*) wavs(i), amps(i)
      enddo
   endif

   ! get the spectra in this region +- wavenumber buffer
   w1 = psnr(1,k)-50.*nterp/opdmax
   w2 = psnr(2,k)+50.*nterp/opdmax
   ilow = minloc(( wavs-w1 ), mask=((wavs-w1) > 0.0D0))
   ihi  = minloc(( wavs-w2 ), mask=((wavs-w2) > 0.0D0))
   iil = ilow(1) - 1
   iih = ihi(1)
   np = iih-iil+1

   write(6,102) 'SNR region : ', psnr(1,k), psnr(2,k)
   write(6,102) 'Resample region : ', wavs(iil), wavs(iih)
   write(6,101) 'Points in resample : ', np

   if( np .le. 2 ) then
       write(*,104) '*** calcsnr : Less than 3 points found(1)'
       stop '1'
   end if

   if( vflag .gt. 1 )then
      write(66,*)'Extended noise region in raw spectrum'
      write(66,*) 2, iih - iil + 1
      do i=iil, iih
         write(66,*) wavs(i), amps(i)
      enddo
   endif

   ! resample and / or degrade resolution
   ! nterp =  0 - skip resample & resolution degradation
   ! nterp =  1 - minimally sample at opdmax
   ! nterp >  1 - interpolate nterp-1 points

   if( nterp .eq. 0 )then
      allocate( outspec( np ))
      outspec = amps(iil:iih)
      dnue    = spac
      wstart  = wavs(iil)
   else
      write(6,105)'Resample snr region...'
      opdm    = opdmax
      dnue    = spac
      wstart  = wavs(iil)
      call sincinterp( amps(iil:iih), outspec, np, wstart, dnue, opdm, nterp )
   endif

   if( vflag .gt. 1 )then
      write(66,*)'Extended noise region in resampled spectrum'
      write(66,*) 3, np
      do i=1, np
         write(66,*) (i-1)*dnue + wstart, outspec(i)
      enddo
   endif

   ! get back the snr sub-region
   allocate( z(np) )
   do i=1, np
      z(i) = (i-1)*dnue + wstart
   enddo

   iil = 0
   do i=1, np
      !print*, dnue, wstart, psnr(1,k)
      z(i) = (i-1)*dnue + wstart
      if(  z(i) .gt. psnr(1,k) )then
         iil = i -1
         exit
      endif
   enddo
   if( iil .eq. 0 )stop '2'
   iih = 0
   do i=iil, np
      if( z(i) .gt. psnr(2,k) )then
         iih = i -1
         exit
      endif
   enddo
   if( iih .eq. 0 )stop '3'

   np = iih - iil +1
   mean = real(sum(outspec(iil:iih)), 8) / real( np, 8 )

   ! assume horizontal band
   !noise = sqrt(dot_product(outspec(iil:iih)-mean, outspec(iil:iih)-mean) / real( np, 8 ) )

   ! fit a parabola
   allocate( x(np), y(np), curve(3) )
   do i=1, np
      x(i) = real(i,8)
   enddo
   !print*, size(x), size( outspec(iil:iih))

   curve(1:3) = polyfit( x, real( outspec(iil:iih), 8 ), np, 2 )
   do i=1, np
      y(i) = outspec(iil+i-1) - (curve(1) + (curve(2) + curve(3)*x(i) ) * x(i))
   enddo

   if( vflag .gt. 1 )then
      write(66,*)'Exact noise region in resampled spectrum, i, w#, spec, fit, diff'
      write(66,*) 4, np, iil*dnue + wstart, dnue
      do i=1, np
         write(66,*) x(i), z(iil+i-1), outspec(iil+i-1), (curve(1) + (curve(2) + curve(3)*x(i) ) * x(i)), y(i)
      enddo
      close(66)
   endif

   noise = 0.0d0
   do i=1, np
      noise = noise + y(i) * y(i)
   enddo

   noise = sqrt( noise / real(np,8))

   write(6,101) 'Points in snr region : ', np
   write(6,102) 'Mean signal : ', mean
   write(6,102) 'RMS noise : ', noise
   write(6,102) 'Mean SNR in snr region : ', mean/noise

   deallocate( x, y, z, curve )

   return

 101  format( a32, i20 )
 102  format( a32, 2f20.12)
 !103  format(/, a32, 2e20.12 )
 104  format( a48 )
 105  format( /, a )


end subroutine calcsnr

end module spec