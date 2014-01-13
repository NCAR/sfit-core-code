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

  module voigt_sdv_lm
    use params

    implicit none

! Mathias Palm September 2012
! this module has been adapted from ACE-FTS_linemixing.for described in
! Boone, C. D.; Walker, K. A. & Bernath, P. F.
! An efficient analytical approach for calculating line mixing in atmospheric remote sensing applications
! J. Quant. Spectrosc. Radiat. Transfer, 2011, 112, 980 - 989

!    double precision ll_freq,temper,pres,ll_ntpwid,fitting(5),&
!         voigt(10000),sdvoigt(10000),dnu
    ! c Wavenumber
    ! ll_freq = 1218.52
    ! c Temperature
    ! temper = 220.d0
    ! c Pressure
    ! pres = 0.283d0
    ! c Air broadening exponent
    ! ll_ntpwid = 0.59d0
    ! c gamma_2
    ! fitting(1) = 0.04d0
    ! c gamma_0
    ! fitting(2) = 0.12d0
    ! c eta_0 (pressure shift)
    ! fitting(3) = -0.006d0
    ! c eta_2
    ! fitting(4) = 0.0096d0

    ! c Choose a wavenumber spacing for the calculations
    ! dnu = 0.0005d0

    ! call sdv(ll_freq,dnu,temper,pres,ll_ntpwid,
    ! *          fitting,voigt,sdvoigt)
    ! stop

    PRIVATE

    real(double) :: fac_dx, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,adop
    real(double) :: dnu,doppler, ds2,g2,eta0,eta2,rtln2,rtpinv, z1,y_keep
    real(double) :: mixval, sqrtdelta, alphadelta,pshift, normalize

    PUBLIC :: sdv_misc, sdvmix, voigtmix


  contains


    subroutine sdv_misc(adop_tmp,tratio,pres,ll_ntpwid,fitting)

      real(double) :: adop_tmp, tratio, pres, ll_ntpwid
      real(double) :: delta, lorwidth, xy, xy2, xz, eta0
      real(double), dimension(4) :: fitting
      ! double precision ll_freq,pres,tratio, &
      !      dopwid,ll_ntpwid,x,xy,xy2,y,dx,adop_tmp, &
      !      nu,xb,xz,z1,lorwidth,dvx,dvy,&
      !      fitting(4), delta

      !print *, 'sdv_misc'
      ! Define some constants
      rtln2 = 0.8325546111576977d0
      rtpinv = 1.d0/1.772453851d0

      lorwidth = fitting(2)*tratio**ll_ntpwid
      g2 = fitting(1)*pres*tratio**ll_ntpwid
      eta0 = pres*fitting(3)
      !c Calculate the line mixing portion.  Here, fitting(4) is the
      !c line mixing coefficient Y.  Assume the same form for the P
      !c and T dependence, and (in this example) assume the same
      !c coefficient (n) for T-dependence
      !c
      mixval = pres*fitting(4)
      adop = adop_tmp
      normalize = 1.0d0
!      normalize = adop*rtpinv

      !c xy = gamma_D / 2 / sqrt(ln2) = omega_o * sqrt(k*T/2/m/c/c)
      !c
      xy = 0.5d0/adop
      xy2 = xy*xy

      ! this is used for SDV, if gamma_2 is 0 use voigt instead!
      if (g2.gt.tiny(g2)) then
         delta = xy2/g2/g2
         sqrtdelta = dsqrt(delta)
         !c alphadelta = alpha + delta, and pshift contains the pressure shift term
         alphadelta = lorwidth/g2 - 1.5d0 + delta
         pshift = -eta0/g2

         ds2 = 1.d0/dsqrt(2.d0)

         a1 = g2*g2+eta2*eta2
         a2 = g2/a1
         a3 = eta2/a1
         a4 = a2*a2-a3*a3
         a5 = 2.d0*a2*a3
         
         xy = 0.5d0/adop
         xy2 = xy*xy
         
         a8 = xy2*a4
         a9 = xy2*a5
         
         a6 = a2*lorwidth - 1.5d0 + a8 + eta0*a3
         a7 = lorwidth*a3 - eta0*a2 + a9
         
         xz = -a3
         
         z1 = 0.195d0*dabs(xz)-0.176d0



         z1 = a9
         if(z1 .eq. 0.) then
            z1 = 0.d0
         else
            z1 = a9/dabs(a9)
         endif
         
         a10 = dsqrt(a8*a8+a9*a9)
         a11 = dsqrt(a10+a8)*ds2
         a12 = z1*dsqrt(a10-a8)*ds2
         
         fac_dx = a1
      end if
      y_keep = lorwidth * adop


    end subroutine sdv_misc


    ! c-----------------------------------------------sdvmix
    ! c subroutine sdvmix
    ! c
    ! c Calculates the speed dependent Voigt profile plus line mixing for a
    ! c range of detuning from line center (+/- 5000 points away from line
    ! c center on the user defined calculation grid).  This simplified routine
    ! c assumes a single line in the center of the calculation window.
    ! c In real applications, the subroutine will need to be edited to
    ! c determine where the line falls on the calculation grid.
    ! c

    subroutine sdvmix(nu, sdmix)

      real(double) :: nu, sdmix

      real(double), dimension(2) :: denon1,numer1,nummix1
      real(double), dimension(4) :: denon2,numer2,nummix2
      real(double), dimension(5) :: denon3,numer3,nummix3
      real(double), dimension(14) :: denon4,numer4,nummix4

      real(double) :: sdv1=0.0, sdv2=0.0, sdvimag=0.0, smix1=0.0, smix2=0.0
      real(double) :: xsdv, xz, y2, yy, yy1, yy2, zz

      integer :: ihum=0,ihum1=0,ihum2=0

      !print *, 'sdvmix'
! c SDV inputs are z1 and z2, where
! c z1 = 1/sqrt(2)*sqrt(sqrt[(alpha + delta)^2 + beta^2] + alpha + delta)
! c    + i*sign(beta)/sqrt(2)*sqrt(sqrt[(alpha + delta)^2 + beta^2] - alpha - delta)
! c       - sqrt(delta)
! c z2 = z1 + 2*sqrt(delta)
! c
! c  alpha = gamma_o / gamma_2 - 1.5
! c  beta = (omega - omega_o + eta0) / gamma_2
! c  delta = gamma_D * gamma_D / gamma_2 / gamma_2 / 4 / ln2
! c
! c Note that the pressure shift (eta0) has been included in beta
! c

       xsdv = nu/g2
       ds2 = 1.d0/dsqrt(2.d0)

!c Imaginary term (beta) in SDV calculation (omega - omega_o + eta0) / g2,
!c where eta0 is the pressure shift term
           sdvimag = pshift + xsdv

!c z1 = sign(beta)
           z1 = sdvimag
           if(z1 .eq. 0.) then
              z1 = 0.d0
           else
              z1 = sdvimag/dabs(sdvimag)
           endif

!c a1 = 	sqrt[(alpha + delta)^2 + beta^2], where alpha, beta, and delta
!c            are defined above
!c
           a1 = dsqrt(alphadelta*alphadelta+sdvimag*sdvimag)

!c a2 = 1/sqrt(2)*sqrt(sqrt[(alpha + delta)^2 + beta^2] + alpha + delta)
           a2 = dsqrt(a1+alphadelta)*ds2

!c a3 = sign(beta)/sqrt(2)*sqrt(sqrt[(alpha + delta)^2 + beta^2] - alpha - delta)
           a3 = z1*dsqrt(a1-alphadelta)*ds2
           xz = -a3

       z1 = 0.195d0*dabs(xz)-0.176d0

! c The speed-dependent Voigt function is the difference between two terms.
! c The calculation assumes the form y - i*x, where y is the real component
! c of the input variable, and x is the imaginary component.
! c
! c  Each term has the same imaginary component (a3), but different real
! c  components.  We use the Humlicek approach for calculating the terms.
! c Ensure that both terms are in the same Humlicek interval.

           yy1 = a2 - sqrtdelta
           yy2 = a2 + sqrtdelta

           zz = dabs(xz) + yy1
           if(zz .ge. 15.d0) then
              ihum1 = 1
           else if(zz .lt. 15.d0 .and. zz .ge. 6.d0) then
              ihum1 = 2
           else if(zz .lt. 6.d0 .and. yy1 .ge. z1) then
              ihum1 = 3
           else if(zz .lt. 6.d0 .and. yy1 .lt. z1) then
              ihum1 = 4
           endif
           zz = dabs(xz) + yy2
           if(zz .ge. 15.d0) then
              ihum2 = 1
           else if(zz .lt. 15.d0 .and. zz .ge. 6.d0) then
              ihum2 = 2
           else if(zz .lt. 6.d0 .and. yy2 .ge. z1) then
              ihum2 = 3
           else if(zz .lt. 6.d0 .and. yy2 .lt. z1) then
              ihum2 = 4
           endif
           ihum = ihum1
    if(ihum2 .gt. ihum1) ihum = ihum2



           if(ihum .eq. 1) then
              call newy1(yy1,denon1,numer1)
              call mix1(yy1,nummix1)
              sdv1 = humlicek1(xz,denon1,numer1)
              smix1 = mixing1(xz,denon1,nummix1)
           else if(ihum .eq. 2) then
              call newy2(yy1,denon2,numer2)
              call mix2(yy1,nummix2)
              sdv1 = humlicek2(xz,denon2,numer2)
              smix1 = mixing2(xz,denon2,nummix2)
           else if(ihum .eq. 3) then
              call newy3(yy1,denon3,numer3)
              call mix3(yy1,nummix3)
              sdv1 = humlicek3(xz,denon3,numer3)
              smix1 = mixing3(xz,denon3,nummix3)
           else if(ihum .eq. 4) then
              call newy4(yy1,denon4,numer4)
              call mix4(yy1,nummix4)
              yy = yy1*yy1
              y2 = 2.d0*yy1
              sdv1 = humlicek4(xz,yy,y2,denon4,numer4)
              smix1 = mixing4(xz,yy,y2,denon4,nummix4)
           else
              write(*,*) 'could not find interval for yy1'
           endif

           if(ihum .eq. 1) then
              call newy1(yy2,denon1,numer1)
              call mix1(yy2,nummix1)
              sdv2 = humlicek1(xz,denon1,numer1)
              smix2 = mixing1(xz,denon1,nummix1)
           else if(ihum .eq. 2) then
              call newy2(yy2,denon2,numer2)
              call mix2(yy2,nummix2)
              sdv2 = humlicek2(xz,denon2,numer2)
              smix2 = mixing2(xz,denon2,nummix2)
           else if(ihum .eq. 3) then
              call newy3(yy2,denon3,numer3)
              call mix3(yy2,nummix3)
              sdv2 = humlicek3(xz,denon3,numer3)
              smix2 = mixing3(xz,denon3,nummix3)
           else if(ihum .eq. 4) then
              call newy4(yy2,denon4,numer4)
              call mix4(yy2,nummix4)
              yy = yy2*yy2
              y2 = 2.d0*yy2
              sdv2 = humlicek4(xz,yy,y2,denon4,numer4)
              smix2 = mixing4(xz,yy,y2,denon4,nummix4)
           else
              write(*,*) 'could not find interval for yy2'
           endif

!c The Voigt and the SDV line shapes have the same normalization factor:
!c sqrt(ln2/pi)/gamma_D
!c
           sdmix = normalize*((sdv1-sdv2) + mixval*(smix1-smix2))


     return
   end subroutine sdvmix
!c-----------------------------------------------sdvmix


! c-----------------------------------------------voigtmix
! c subroutine voigtmix
! c
! c Calculates the line profile (Voigt plus line mixing) for a range of
! c detuning from line center (+/- 5000 points away from line center
! c on the user defined calculation grid).  This simplified routine
! c assumes a single line in the center of the calculation window.
! c In real applications, the subroutine will need to be edited to
! c determine where the line falls on the calculation grid, as well
! c as the range of points for which calculations are to be performed.
! c
! c inputs: ll_freq (line center in wavenumbers)
! c         dnu (point spacing on the calculation grid in wavenumbers)
! c         temper (temperature in K)
! c	  pres (pressure in atmospheres)
! c	  ll_ntpwid (air broadening temperature exponent)
! c	  mass (molecular mass of the molecule in question)
! c	  fitting() (an array containing values for calculating the
! c	               line profile)
! c
! c output: vgtmix (the calculated line profile, Voigt plus line mixing,
! c
! for the given set of conditions).

   subroutine voigtmix(nu, vgtmix)
     real(double) :: nu, vgtmix

     real(double), dimension(2) :: denon1,numer1,nummix1
     real(double), dimension(4) :: denon2,numer2,nummix2
     real(double), dimension(5) :: denon3,numer3,nummix3
     real(double), dimension(14) :: denon4,numer4,nummix4

     real(double) :: yy, y2, zv, zz
     !print *, 'voigtmix'

     nu = nu*adop

     z1 = 0.195d0*dabs(nu)-0.176d0
     zz = dabs(nu) + y_keep
     if(zz .ge. 15.d0) then
        call newy1(y_keep,denon1,numer1)
        call mix1(y_keep,nummix1)
        !c Voigt profile
        zv = humlicek1(nu,denon1,numer1)
        !c Voigt with line mixing
        vgtmix = zv + mixval*mixing1(nu,denon1,nummix1)
     endif
     if(zz .lt. 15.d0 .and. zz .ge. 5.5) then
        call newy2(y_keep,denon2,numer2)
        call mix2(y_keep,nummix2)
        !c Voigt profile
        zv = humlicek2(nu,denon2,numer2)
        !c Voigt with line mixing
        vgtmix = zv + mixval*mixing2(nu,denon2,nummix2)
     endif
     if(zz .lt. 5.5d0 .and. y_keep .ge. z1) then
        call newy3(y_keep,denon3,numer3)
        call mix3(y_keep,nummix3)
        !c Voigt profile
        zv = humlicek3(nu,denon3,numer3)
        !c Voigt with line mixing
        vgtmix = zv + mixval*mixing3(nu,denon3,nummix3)
     endif
     if(zz .lt. 5.5d0 .and. y_keep .lt. z1) then
        call newy4(y_keep,denon4,numer4)
        call mix4(y_keep,nummix4)
        yy = y_keep*y_keep
        y2 = 2.d0*y_keep
        !c Voigt profile
        zv = humlicek4(nu,yy,y2,denon4,numer4)
        !c Voigt with line mixing
        vgtmix = zv + mixval*mixing4(nu,yy,y2,denon4,nummix4)
     endif

     vgtmix = normalize*vgtmix
     return
   end subroutine voigtmix
!c-----------------------------------------------voigtmix


! c----------------------------------------humlicek1
! c function humlicek1
! c
! c Calculates the real component K(x,y) of the Humlicek approximate
! c solution for the complex probability function in "region 1."
! c
! c inputs: x (= sqrt(ln2)* [omega - omega_o + eta] / [Doppler width])
! c                  where omega - omega_o is the detuning from the center of
! c                  the line in question (located at position omega_o), and
! c                  eta is the pressure shift at the given pressure and temperature.
! c         denon (coefficients in the polynomial for the denominator in the
! c                    Humlicek algorithm - previously calculated in the
! c                    subroutine newy1)
! c         numer (coefficients in the polynomial for the numerator in the
! c                    Humlicek algorithm - previously calculated in the
! c                     subroutine newy1)
! c
! c output: humlicek1 (approximate solution for K(x,y), the real component of the
! c                    complex probability function for Humlicek's region 1)
! c
      double precision function humlicek1(x,denon,numer)
      double precision x,xx,denon(2),numer(2)
      xx = x*x
      humlicek1 = (numer(1) + numer(2)*xx)/(denon(1)+xx*(denon(2)+xx))
      return
    end function humlicek1
!c----------------------------------------humlicek1

!c----------------------------------------mixing1
 ! function mixing1

 ! Calculates the imaginary component L(x,y) of the Humlicek approximate
 ! solution for the complex probability function in "region 1."

 ! inputs: x (= sqrt(ln2)* [omega - omega_o + eta] / [Doppler width])
 !                  where omega - omega_o is the detuning from the center of
 !                  the line in question (located at position omega_o), and
 !                  eta is the pressure shift at the given pressure and temperature.
 !         denon (coefficients in the polynomial for the denominator in the
 !                    Humlicek algorithm - previously calculated in the
 !                    subroutine newy1)
 !         numer (coefficients in the polynomial for the numerator in the
 !                    Humlicek algorithm - previously calculated in the
 !                     subroutine mix1)

 ! output: mixing1 (approximate solution for L(x,y), the imaginary component of
 !                    the complex probability function for Humlicek's region 1)

      double precision function mixing1(x,denon,numer)
      double precision x,xx,denon(2),numer(2)
      xx = x*x
      mixing1 = x*(numer(1) + numer(2)*xx)/(denon(1)+xx*(denon(2)+xx))
      return
    end function mixing1
!c----------------------------------------mixing1


! c----------------------------------------newy1
! c subroutine newy1
! c
! c Generates expressions required to calculate the numerator and
! c denominator for the real component K(x,y) of the Humlicek
! c approximate solution for the complex probability function in
! c "region 1."
! c
! c input: y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c
! c outputs: numer (coefficients in the polynomial for the numerator in
! c                   the real part of the Humlicek algorithm)
! c          denon (coefficients in the polynomial for the denominator in
! c                    the Humlicek algorithm)
! c
      subroutine newy1(y,denon,numer)
      double precision denon(2),numer(2),y,yy
      yy = y*y
      numer(1) = y*(0.2820948d0+0.5641896d0*yy)
      numer(2) = 0.5641896d0*y
      denon(1) = 0.25d0 + yy + yy*yy
      denon(2) = 2.d0*yy - 1.d0
      return
    end subroutine newy1
!c----------------------------------------newy1


! c----------------------------------------newy2
! c subroutine newy2
! c
! c Generates expressions required to calculate the numerator and
! c denominator for the real component K(x,y) of the Humlicek
! c approximate solution for the complex probability function in
! c "region 2."
! c
! c input: y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c
! c outputs: numer (coefficients in the polynomial for the numerator in
! c                   the real part of the Humlicek algorithm)
! c          denon (coefficients in the polynomial for the denominator in
! c                    the Humlicek algorithm)
! c
      subroutine newy2(y,denon,numer)
      double precision denon(4),numer(4),y,yy
      yy = y*y
      numer(1) = y*(1.0578555d0+yy*(4.6545642d0 + yy*(3.1030428d0 + &
           yy*0.5641896d0)))
      numer(2) = y*(2.9619954d0+yy*(0.5641896d0 + yy*1.6925688d0))
      numer(3) = -y*(2.5388532d0 - yy*1.6925688d0)
      numer(4) = y*0.5641896d0
      denon(1) = 0.5625d0 + yy*(4.5d0 + yy*(10.5d0 + yy*(6.d0 + yy)))
      denon(2) = yy*(9.d0 + yy*(6.d0 + yy*4.d0)) - 4.5d0
      denon(3) = 10.5d0 + 6.d0*yy*(yy - 1.d0)
      denon(4) = 4.d0*yy - 6.d0
      return
    end subroutine newy2
!c----------------------------------------newy2


! c----------------------------------------mix1
! c subroutine mix1
! c
! c Generates expressions required to calculate the numerator for the
! c imaginary component L(x,y) of the Humlicek approximate solution for
! c the complex probability function in "region 1."
! c
! c input: y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c
! c output: nummix (coefficients in the polynomial for the numerator in
! c                      the imaginary part of the Humlicek algorithm)
! c
      subroutine mix1(y,nummix)
      double precision y, nummix(2)
      nummix(1) = .5641896d0*y*y - .2820948d0
      nummix(2) = .5641896d0
      return
    end subroutine mix1
!c----------------------------------------mix1




! c----------------------------------------mix2
! c subroutine mix2
! c
! c Generates expressions required to calculate the numerator for the
! c imaginary component L(x,y) of the Humlicek approximate solution for
! c the complex probability function in "region 2."
! c
! c input: y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c
! c output: nummix (coefficients in the polynomial for the numerator in
! c                      the imaginary part of the Humlicek algorithm)
! c
      subroutine mix2(y,nummix)
      double precision y,y2,nummix(4)
      y2 = y*y
      nummix(1) = -1.0578555 + y2*(2.9619954 + y2*(2.5388532 + y2*0.5641896))
      nummix(2) = 4.6545642 +y2*(-0.5641896 + y2*1.6925688)
      nummix(3) = -3.1030428 + 1.6925688*y2
      nummix(4) = 0.5641896
      return
    end subroutine mix2
!c----------------------------------------mix2

! c----------------------------------------humlicek2
! c function humlicek2
! c
! c Calculates the real component K(x,y) of the Humlicek approximate
! c solution for the complex probability function in "region 2."
! c
! c inputs: x (= sqrt(ln2)* [omega - omega_o + eta] / [Doppler width])
! c                  where omega - omega_o is the detuning from the center of
! c                  the line in question (located at position omega_o), and
! c                  eta is the pressure shift at the given pressure and temperature.
! c         denon (coefficients in the polynomial for the denominator in the
! c                    Humlicek algorithm - previously calculated in the
! c                    subroutine newy2)
! c         numer (coefficients in the polynomial for the numerator in the
! c                    Humlicek algorithm - previously calculated in the
! c                     subroutine newy2)
! c
! c output: humlicek2 (approximate solution for K(x,y), the real component of the
! c                    complex probability function for Humlicek's region 2)
! c
      double precision function humlicek2(x,denon,numer)
      double precision x,xx,denon(4),numer(4)
      xx = x*x
      humlicek2 = (numer(1) + xx*(numer(2)+xx*(numer(3)+xx*numer(4))))&
           /(denon(1)+xx*(denon(2)+xx*(denon(3)+xx*(denon(4)+xx))))
      return
    end function humlicek2
!c----------------------------------------humlicek2

!c----------------------------------------mixing2
! c function mixing2
! c
! c Calculates the imaginary component L(x,y) of the Humlicek approximate
! c solution for the complex probability function in "region 2."
! c
! c inputs: x (= sqrt(ln2)* [omega - omega_o + eta] / [Doppler width])
! c                  where omega - omega_o is the detuning from the center of
! c                  the line in question (located at position omega_o), and
! c                  eta is the pressure shift at the given pressure and temperature.
! c         denon (coefficients in the polynomial for the denominator in the
! c                    Humlicek algorithm - previously calculated in the
! c                    subroutine newy2)
! c         numer (coefficients in the polynomial for the numerator in the
! c                    Humlicek algorithm - previously calculated in the
! c                     subroutine mix2)
! c
! c output: mixing2 (approximate solution for L(x,y), the imaginary component of
! c                    the complex probability function for Humlicek's region 2)
! c
      double precision function mixing2(x,denon,numer)
      double precision x,xx,denon(4),numer(4)
      xx = x*x
      mixing2 = x*(numer(1)+xx*(numer(2)+xx*(numer(3)+xx*numer(4))))&
           /(denon(1)+xx*(denon(2)+xx*(denon(3)+xx*(denon(4)+xx))))
      return
    end function mixing2
!c----------------------------------------mixing2

! c----------------------------------------humlicek3
! c function humlicek3
! c
! c Calculates the real component K(x,y) of the Humlicek approximate
! c solution for the complex probability function in "region 3."
! c
! c inputs: x (= sqrt(ln2)* [omega - omega_o + eta] / [Doppler width])
! c                  where omega - omega_o is the detuning from the center of
! c                  the line in question (located at position omega_o), and
! c                  eta is the pressure shift at the given pressure and temperature.
! c         denon (coefficients in the polynomial for the denominator in the
! c                    Humlicek algorithm - previously calculated in the
! c                    subroutine newy3)
! c         numer (coefficients in the polynomial for the numerator in the
! c                    Humlicek algorithm - previously calculated in the
! c                     subroutine newy3)
! c
! c output: humlicek3 (approximate solution for K(x,y), the real component of the
! c                    complex probability function for Humlicek's region 3)
! c
      double precision function humlicek3(x,denon,numer)
      double precision x,xx,denon(5),numer(5)
      xx = x*x
      humlicek3 = (numer(1) + xx*(numer(2)+xx*(numer(3)+xx*  &
           (numer(4)+xx*numer(5)))))/(denon(1)+xx*(denon(2)+ &
           xx*(denon(3)+xx*(denon(4)+xx*(denon(5)+xx)))))
      return
    end function humlicek3
!c----------------------------------------humlicek3

! c----------------------------------------mixing3
! c function mixing3
! c
! c Calculates the imaginary component L(x,y) of the Humlicek approximate
! c solution for the complex probability function in "region 3."
! c
! c inputs: x (= sqrt(ln2)* [omega - omega_o + eta] / [Doppler width])
! c                  where omega - omega_o is the detuning from the center of
! c                  the line in question (located at position omega_o), and
! c                  eta is the pressure shift at the given pressure and temperature.
! c         denon (coefficients in the polynomial for the denominator in the
! c                    Humlicek algorithm - previously calculated in the
! c                    subroutine newy3)
! c         numer (coefficients in the polynomial for the numerator in the
! c                    Humlicek algorithm - previously calculated in the
! c                     subroutine mix3)
! c
! c output: mixing3 (approximate solution for L(x,y), the imaginary component of
! c                    the complex probability function for Humlicek's region 3)
! c
      double precision function mixing3(x,denon,numer)
      double precision x,xx,denon(5),numer(5)
      xx = x*x
      mixing3 = x*(numer(1) + xx*(numer(2)+xx*(numer(3)+xx* &
           (numer(4)+xx*numer(5)))))/(denon(1)+xx*(denon(2)+ &
           xx*(denon(3)+xx*(denon(4)+xx*(denon(5)+xx)))))
      return
    end function mixing3
!c----------------------------------------mixing3

!c----------------------------------------newy3
!c subroutine newy3
!c
!c Generates expressions required to calculate the numerator and
!c denominator for the real component K(x,y) of the Humlicek
!c approximate solution for the complex probability function in
!c "region 3."
! c
! c input: y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c
! c outputs: numer (coefficients in the polynomial for the numerator in
! c                   the real part of the Humlicek algorithm)
! c          denon (coefficients in the polynomial for the denominator in
! c                    the Humlicek algorithm)
! c
      subroutine newy3(y,denon,numer)
      double precision denon(5),numer(5),y
      numer(1) = 272.102d0 + y*(973.778d0+y*(1629.76d0+y*(1678.33d0+y*&
           (1174.8d0+y*(581.746d0+y*(204.51d0+y*(49.5213d0+y*(7.55895d0+&
           y*0.564224d0))))))))
      numer(2) = y*(-2.34403d0+y*(220.843d0+y*(336.364d0+y*(247.198d0+ &
           y*(100.705d0+y*(22.6778d0+y*2.25689d0)))))) - 60.5644d0
      numer(3) = 4.58029d0 + y*(18.546d0+y*(42.5683d0+y*(52.8454d0+ &
           y*(22.6798d0+y*3.38534d0))))
      numer(4) = y*(1.66203d0+y*(7.56186d0+y*2.25689d0)) - 0.128922d0
      numer(5) = 0.000971457d0 + 0.564224d0*y
      denon(1) = 272.102d0+y*(1280.83d0+y*(2802.87d0+y*(3764.97d0+y* &
           (3447.63d0+y*(2256.98d0+y*(1074.41d0+y*(369.199d0+y*(88.2674d0 &
           +y*(13.3988 + y)))))))))
      denon(2) = 211.678d0+y*(902.306d0+y*(1758.34d0+y*(2037.31d0+y* &
           (1549.68d0+y*(793.427d0+y*(266.299d0+y*(53.5952d0+5.d0*y)))))))
      denon(3) = 78.866d0+y*(308.186d0+y*(497.302d0+y*(479.258d0+y* &
           (269.292d0+y*(80.3928d0+y*10.d0)))))
      denon(4) = 22.0353d0+y*(55.0293d0+y*(92.7586d0+y*(53.5952d0+ &
           y*10.d0)))
      denon(5) = 1.49645d0+y*(13.3988d0+y*5.d0)
      return
    end subroutine newy3
!c----------------------------------------newy3



! c----------------------------------------mix3
! c subroutine mix3
! c
! c Generates expressions required to calculate the numerator for the
! c imaginary component L(x,y) of the Humlicek approximate solution for
! c the complex probability function in "region 3."
! c
! c input: y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c
! c output: nummix (coefficients in the polynomial for the numerator in
! c                      the imaginary part of the Humlicek algorithm)
! c
      subroutine mix3(y,nummix)
      double precision y,nummix(5)
      nummix(1) = 307.0522 + y*(900.8651 + y*(1215.616 + y*(&
           988.17457 + y*(534.07725 + y*(196.836 + y*(&
           48.97184 + y*(7.557974 + y*0.5642236)))))))
      nummix(2) = 33.630785 + y*(178.55242 + y*(284.30395 + y*(&
           231.99739 + y*(99.056246 + y*(22.673922 +     &
           y*2.2568944)))))
      nummix(3) = 14.1547 + y*(35.161375 + y*(51.196966 + y*( &
           22.673922 + y*3.3853416)))
      nummix(4) = 1.1125621 + y*(7.557974 + y*2.2568944)
      nummix(5) = 0.5642236
      return
    end subroutine mix3
!c----------------------------------------mix3

! c----------------------------------------newy4
! c subroutine newy4
! c
! c Generates expressions required to calculate the numerator and
! c denominator for the real component K(x,y) of the Humlicek
! c approximate solution for the complex probability function in
! c "region 4."
! c
! c input: y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c
! c outputs: numer (coefficients in the polynomial for the numerator in
! c                   the real part of the Humlicek algorithm)
! c          denon (coefficients in the polynomial for the denominator in
! c                    the Humlicek algorithm)
! c
      subroutine newy4(y,denon,numer)
      double precision denon(14),numer(14),y,yy
      yy = y*y
      numer(1) = y*(1.16028d9-yy*(9.86606d8-yy*(4.56662d8-yy*(1.53575d8 &
           -yy*(4.08168d7-yy*(9.69464d6-yy*(1.6841d6-yy*(320772.d0-yy*     &
           (40649.2d0-yy*(5860.68d0-yy*(571.688d0-yy*(72.9359d0-yy*        &
           (2.35944d0-yy*0.56419d0)))))))))))))
      numer(2) = y*(-5.60506d8-yy*(9.85386d8-yy*(8.06985d8-yy*           &
           (2.91877d8-yy*(8.64829d7-yy*(7.72359d6-yy*(3.59915d6-yy*(        &
           234416.d0-yy*(45251.3d0-yy*(2269.19d0+yy*(234.143d0-yy*          &
           (23.0312d0-yy*7.33447d0))))))))))))
      numer(3) = y*(-6.51523d8+yy*(2.47157d8+yy*(2.94262d8-yy*(2.04467d8 &
           -yy*(2.29303d7-yy*(2.3818d7-yy*(576054.d0+yy*(98079.1d0-yy*(     &
           25338.3d0-yy*(1097.77d0+yy*(97.6203d0-yy*44.0068d0)))))))))))
      numer(4) = y*(-2.63984d8+yy*(2.70166d8-yy*(9.96224d7+yy*(4.15013d7  &
           -yy*(3.83112d7+yy*(2.2404d6-yy*(303569.d0+yy*(66431.3d0-yy*(   &
           8381.97d0+yy*(228.563d0-yy*161.358d0))))))))))
      numer(5) = y*(-6.31771d7+yy*(1.40676d8+yy*(5.56966d6+yy*(2.46201d7  &
           +yy*(468141.d0-yy*(1.003d6+yy*(66212.1d0-yy*(23507.6d0+yy*(    &
           296.38d0-yy*403.396d0)))))))))
      numer(6) = y*(-1.69846d7+yy*(4.07381d6-yy*(3.32896d7+yy*(1.93114d6  &
           +yy*(934717.d0-yy*(8821.d0+yy*(37544.8d0+yy*(125.589d0-        &
           yy*726.112d0))))))))
      numer(7) = y*(-1.23164d6+yy*(7.52883d6-yy*(900011.d0+yy*(186682.d0  &
           -yy*(79902.7d0+yy*(37371.9d0-yy*(260.198d0+yy*968.15d0)))))))
      numer(8) = y*(-610621.d0+yy*(86407.7d0+yy*(153468.d0+yy*(72521.d0   &
           +yy*(23137.1d0-yy*(571.645d0+yy*968.15d0))))))
      numer(9) = y*(-23586.5d0+yy*(49883.8d0+yy*(26538.6d0+yy*(8073.14d0  &
           -yy*(575.165d0+yy*726.112d0)))))
      numer(10) = y*(-8009.1d0+yy*(2198.85d0+yy*(953.65-yy*(352.467d0     &
           +yy*403.395d0))))
      numer(11) = y*(-622.058d0-yy*(271.203d0+yy*(134.792d0+              &
           yy*161.358d0)))
      numer(12) = -y*(77.0536d0+yy*(29.7896d0+yy*44.0068d0))
      numer(13) = -y*(2.92264d0+yy*7.33447d0)
      numer(14) = -0.56419d0*y

      denon(1) = 1.02827d9-yy*(1.5599d9-yy*(1.17022d9-yy*(5.79099d8-yy*    &
           (2.11107d8-yy*(6.11148d7-yy*(1.44648d7-yy*(2.85721d6-yy*        &
           (483737.d0-yy*(70946.1d0-yy*(9504.65d0-yy*(955.194d0-yy*        &
           (126.532d0-yy*(3.68288d0-yy)))))))))))))
      denon(2) = 1.5599d9-yy*(2.28855d9-yy*(1.66421d9-yy*(7.5383d8-yy*     &
           (2.89676d8-yy*(7.01358d7-yy*(1.39465d7-yy*(2.84954d6-yy*        &
           (498334.d0-yy*(55600.d0-yy*(3058.26d0+yy*(533.254d0-yy*         &
           (40.5117d0-yy*14.d0))))))))))))
      denon(3) = 1.17022d9-yy*(1.66421d9-yy*(1.06002d9-yy*(6.60077d8       &
           -yy*(6.33497d7-yy*(4.60396d7-yy*(1.4841d7-yy*(1.06352d6+yy*     &
           (217801.d0-yy*(48153.3d0-yy*(1500.17d0+yy*(198.875d0-           &
           yy*91.d0)))))))))))
      denon(4) = 5.79099d8-yy*(7.5383d8-yy*(6.60077d8+yy*(5.40371d7+       &
           yy*(1.99846d8-yy*(6.87655d6+yy*(6.89002d6-yy*(280427.d0+yy*     &
           (161461.d0-yy*(16493.7d0+yy*(567.163d0-yy*364.d0))))))))))
      denon(5) = 2.11107d8-yy*(2.89676d8-yy*(6.33497d7-yy*(1.99846d8        &
           +yy*(5.01017d7+yy*(5.25722d6-yy*(1.9547d6+yy*(240373.d0-yy*     &
           (55582.d0+yy*(1012.79d0-yy*1001.d0)))))))))
      denon(6) = 6.11148d7-yy*(7.01358d7-yy*(4.60396d7-yy*(6.87655d6       &
           -yy*(5.25722d6+yy*(3.04316d6+yy*(123052.d0-yy*(106663.d0+yy*    &
           (1093.81-yy*2002.d0))))))))
      denon(7) = 1.44648d7-yy*(1.39465d7-yy*(1.4841d7+yy*(6.89002d6        &
           +yy*(1.9547d6-yy*(123052.d0+yy*(131336.5d0+yy*(486.14d0         &
           -yy*3003.d0)))))))
      denon(8) = 2.85721d6-yy*(2.84955d6-yy*(1.06352d6+yy*(280427.d0       &
           -yy*(240373.d0+yy*(106663.d0-yy*(486.14d0+yy*3432.d0))))))
      denon(9) = 483737.d0-yy*(498334.d0+yy*(217801.d0+yy*(161461.d0       &
           +yy*(55582.d0-yy*(1093.81d0+yy*3003.d0)))))
      denon(10) = 70946.1d0-yy*(55600.d0+yy*(48153.3d0+yy*(16493.7d0       &
           -yy*(1012.79d0+yy*2002.d0))))
      denon(11) = 9504.65d0-yy*(3058.26d0+yy*(1500.17d0-yy*(567.163d0      &
           +yy*1001.d0)))
      denon(12) = 955.194d0+yy*(533.254d0+yy*(198.875d0+yy*364.d0))
      denon(13) = 126.532d0+yy*(40.5117d0+yy*91.d0)
      denon(14) = 3.68288d0+14.d0*yy
      return
    end subroutine newy4
!c----------------------------------------newy4


! c----------------------------------------humlicek4
! c function humlicek4
! c
! c Calculates the real component K(x,y) of the Humlicek approximate
! c solution for the complex probability function in "region 4."
! c
! c inputs: x (= sqrt(ln2)* [omega - omega_o + eta] / [Doppler width])
! c                  where omega - omega_o is the detuning from the center of
! c                  the line in question (located at position omega_o), and
! c                  eta is the pressure shift at the given pressure and temperature.
! c         y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c         yy = y * y
! c         y2 = 2.0*y
! c         denon (coefficients in the polynomial for the denominator in the
! c                    Humlicek algorithm - previously calculated in the
! c                    subroutine newy4)
! c         numer (coefficients in the polynomial for the numerator in the
! c                    Humlicek algorithm - previously calculated in the
! c                     subroutine newy4)
! c
! c output: humlicek4 (approximate solution for K(x,y), the real component of the
! c                    complex probability function for Humlicek's region 4)
! c
      double precision function humlicek4(x,yy,y2,denon,numer)
      double precision x,xx,yy,y2,denon(14),numer(14)
      xx = x*x
      humlicek4 = (numer(1) + xx*(numer(2)+xx*(numer(3)+xx*               &
           (numer(4)+xx*(numer(5)+xx*(numer(6)+xx*(numer(7)+xx*(numer(8)+ &
           xx*(numer(9)+xx*(numer(10)+xx*(numer(11)+xx*(numer(12)+        &
           xx*(numer(13)+xx*numer(14))))))))))))))/                       &
           (denon(1)+xx*(denon(2)+xx*(denon(3)+xx*(denon(4)+xx*           &
           (denon(5)+xx*(denon(6)+xx*(denon(7)+xx*(denon(8)+xx*           &
           (denon(9)+xx*(denon(10)+xx*(denon(11)+xx*(denon(12)+xx*        &
           (denon(13)+xx*(denon(14)+xx))))))))))))))
      humlicek4 = dexp(yy-xx)*cos(y2*x) - humlicek4
      return
    end function humlicek4
!c----------------------------------------humlicek4

! c----------------------------------------mixing4
! c function mixing4
! c
! c Calculates the imaginary component L(x,y) of the Humlicek approximate
! c solution for the complex probability function in "region 4."
! c
! c inputs: x (= sqrt(ln2)* [omega - omega_o + eta] / [Doppler width])
! c                  where omega - omega_o is the detuning from the center of
! c                  the line in question (located at position omega_o), and
! c                  eta is the pressure shift at the given pressure and temperature.
! c         denon (coefficients in the polynomial for the denominator in the
! c                    Humlicek algorithm - previously calculated in the
! c                    subroutine newy4)
! c         numer (coefficients in the polynomial for the numerator in the
! c                    Humlicek algorithm - previously calculated in the
! c                     subroutine mix4)
! c
! c output: mixing4 (approximate solution for L(x,y), the imaginary component of
! c                    the complex probability function for Humlicek's region 4)
! c
      double precision function mixing4(x,yy,y2,denon,numer)
      double precision x,xx,yy,y2,denon(14),numer(14)
      xx = x*x
      mixing4 = x*(numer(1) + xx*(numer(2)+xx*(numer(3)+xx*                &
           (numer(4)+xx*(numer(5)+xx*(numer(6)+xx*(numer(7)+xx*(numer(8)+  &
           xx*(numer(9)+xx*(numer(10)+xx*(numer(11)+xx*(numer(12)+         &
           xx*(numer(13)+xx*numer(14))))))))))))))/                        &
           (denon(1)+xx*(denon(2)+xx*(denon(3)+xx*(denon(4)+xx*            &
           (denon(5)+xx*(denon(6)+xx*(denon(7)+xx*(denon(8)+xx*            &
           (denon(9)+xx*(denon(10)+xx*(denon(11)+xx*(denon(12)+xx*         &
           (denon(13)+xx*(denon(14)+xx))))))))))))))
      mixing4 = dexp(yy-xx)*dsin(y2*x) - mixing4
      return
    end function mixing4
!c----------------------------------------mixing4

! c----------------------------------------mix4
! c subroutine mix4
! c
! c Generates expressions required to calculate the numerator for the
! c imaginary component L(x,y) of the Humlicek approximate solution for
! c the complex probability function in "region 4."
! c
! c input: y (= sqrt(ln2)* [pressure broadened width] / [Doppler width])
! c
! c output: nummix (coefficients in the polynomial for the numerator in
! c                      the imaginary part of the Humlicek algorithm)




    subroutine mix4(y,nummix)
      double precision y,y2,nummix(14)
      y2 = y*y
      nummix(1) = -1.1602757d9 + y2*(-5.6050464d8 + y2*(6.5152333d8 &
           + y2*(-2.6389408d8 + y2*(6.31771d7 + y2*(-1.6984609d7     &
           + y2*(1.2316455d6 +y2*(-6.10622d5 + y2*(23586.527         &
           + y2*(-8009.0976 + y2*(622.05559 + y2*(-77.05345          &
           + y2*(2.922644 - y2*0.56419))))))))))))
      nummix(2) = -9.8660434d8 + y2*(9.85361d8 + y2*(2.4715686d8    &
           + y2*(-2.7016652d8 + y2*(1.4067656d8 + y2*(-4.0738161d6   &
           + y2*(7.5288304d6 + y2*(-84607.61 + y2*(49883.805         &
           + y2*(-2198.8554 +y2*(-271.20157 + y2*(29.789641          &
           - y2*7.33447)))))))))))
      nummix(3) = -4.5666204d8 + y2*(8.0698521d8 + y2*(-2.942619d8   &
           + y2*(-9.9622353d7 + y2*(-5.569648d6 + y2*(-3.3289558d7   &
           + y2*(9.000102d5 + y2*(1.5346779d5 + y2*(-26538.547       &
           + y2*(953.655 + y2*(134.79153 - y2*44.00682))))))))))
      nummix(4) = -1.5357482d8 + y2*(2.9187633d8 + y2*(-2.0446691d8  &
           + y2*(4.1501295d7 + y2*(2.4620116d7 + y2*(1.9311433d6       &
           + y2*(-1.866818d5 + y2*(-72520.91 + y2*(8073.149            &
           + y2*(352.4668 - y2*161.3582)))))))))
      nummix(5) = -4.0816819d7 + y2*(8.6482898d7 + y2*(-2.293025d7   &
           + y2*(3.83112d7 + y2*(-4.681421d5 + y2*(-9.347168d5         &
           + y2*(-79902.49 + y2*(23137.06 + y2*(575.167              &
           - y2*403.395))))))))
      nummix(6) = -9.6946319d6 + y2*(7.7235907d6 +y2*(-2.3818001d7   &
           + y2*(-2.2404003d6 +y2*(-1.002996d6 +y2*(-8821.             &
           + y2*(37371.88 + y2*(571.647 +y2*(-726.113))))))))
      nummix(7) = -1.6840975d6 + y2*(3.5991546d6 +y2*(-5.7605388d5   &
           + y2*(-3.0356892d5 + y2*(66212.07 + y2*(37544.8             &
           + y2*(260.197 - y2*968.15))))))
      nummix(8) = -3.2077208d5 + y2*(2.3441652d5 + y2*(98079.06       &
           + y2*(66431.22 + y2*(23507.63 + y2*(-125.59                  &
           - y2*968.15)))))
      nummix(9) = -40649.182 + y2*(45251.326 + y2*(25338.301           &
           + y2*(8381.969 + y2*(-296.38 - y2*726.113))))
      nummix(10) = -5860.6757 + y2*(2269.1863 + y2*(1097.7707          &
           +y2*(-228.5628 - y2*403.3958)))
      nummix(11) = -571.6872 + y2*(-234.14331 + y2*(-97.62033          &
           - y2*161.35834))
      nummix(12) = -72.935866 - y2*(23.03124 + y2*44.00682)
      nummix(13) = -2.359444 - 7.33447*y2
      nummix(14) = -0.56419
      return
    end subroutine mix4
!c----------------------------------------mix4



  end module voigt_sdv_lm
