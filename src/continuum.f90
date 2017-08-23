module continuum

  use params
  use xsections

  implicit none
  
  integer, parameter :: cont_poly_max = 10

  integer :: ncont = 0, n_contabs
  logical :: f_continuum = .false., f_mtckd = .false.
  integer :: abscont_type, abscont_order
  real(double), dimension(cont_poly_max) :: abscont_param, abscont_sparam
  real(double), dimension(:), allocatable :: cont_param
  real(double), dimension(:,:,:), allocatable :: mtckd
  real(double) :: cont_z_abs, cont_alpha
  
  
contains

  subroutine calc_continuum(param)
    ! wrapper for calculation of the continua. The continuums absorption gets
    ! inserted into the cross sections after the gas crosssections

    implicit none

    real(double), dimension(:), intent(in) :: param
    

    integer :: iband, k, j, l
    integer :: mone, mxne
    real(double) :: wone, wxne, wmid
    real(double) :: polynom

    if (f_continuum) then
    
       select case (abscont_type) 
          !     case (0)
          !        ! Offset only
          !        !  --- Loop over layers
          !        DO K = 1, KMAX
          !           CROSS(nret+2,K,:NCROSS) = param(1)*(P(k)/P(kmax))
          !        end DO
          !     case (0)
          !        ! Offset and slope
          !        mone = 1
          !        DO IBAND = 1, NBAND
          !           mxne = mone + nm(iband) - 1
          !           wone = wstart(iband)
          !           wxne = wstart(iband) + dn(iband)*nm(iband)
          !           wmid = (wone + wxne)/2.0d0
          ! !          print *, mone, mxne, wone, wxne, wmid, dn(iband)
          ! !          print *, param(1), param(2), abscont_strength, abscont_tilt
          !           DO K = 1, KMAX     
          !              CROSS(nret+2,K,mone:mxne) = param(1)*(P(k)/P(kmax))
          !              do j = mone,mxne
          !                 CROSS(nret+2,K,j) = CROSS(nret+2,K,j) * ( 1.0D0 + param(2)*((wone+dble(j)*dn(iband))-wmid)/(wxne-wone) )
          !              end do
          !           end DO
          !           mone = mone + nm(iband)
          !        end DO
       case (1)
          ! n-th order polynomial
          mone = 1
          CROSS(nret+2,:KMAX,:ncross) = 0.0d0
          DO IBAND = 1, NBAND
             mxne = mone + nm(iband) - 1
             wone = wstart(iband)
             wxne = wstart(iband) + dn(iband)*nm(iband)
             wmid = (wone + wxne)/2.0d0
             !          print *, mone, mxne, wone, wxne, wmid, dn(iband)
             !          print *, param(1), param(2), abscont_strength, abscont_tilt
             DO K = 1, KMAX     
                do j = mone,mxne
                   polynom = param(1)
                   do l = 1,n_contabs-1
                      polynom = polynom + param(l+1)*(((wone+dble(j)*dn(iband))-wmid)/(wxne-wone))**l
                   end do
                   CROSS(nret+2,K,j) = CROSS(nret+2,K,j) + polynom*(P(k)/P(kmax))
                   !                print *,  CROSS(nret+2,K,j)
                end do
             end DO
             mone = mone + nm(iband)
          end DO
       case(2)
          ! an absorbing layer modeled at the altitude
          ! z_abs with the absorbing strength cont_alpha (can be retrieved)
          mone = 1
          CROSS(nret+2,:KMAX,:ncross) = 0.0d0
          cont_alpha = param(1)
          DO IBAND = 1, NBAND
             mxne = mone + nm(iband) - 1
             wone = wstart(iband)
             wxne = wstart(iband) + dn(iband)*nm(iband)
             wmid = (wone + wxne)/2.0d0
             DO K = 1, KMAX-1
                IF ((z(k).ge.cont_z_abs).and.(z(k+1).le.cont_z_abs)) then
                   do j = mone,mxne
                      polynom = param(1)
                      do l = 1,n_contabs-1
                         polynom = polynom + param(l+1)*(((wone+dble(j)*dn(iband))-wmid)/(wxne-wone))**l
                      end do
                      CROSS(nret+2,K,j) = CROSS(nret+2,K,j) - polynom
                   end do
                end IF
             end DO
             mone = mone + nm(iband)
          end DO
       end select
    end if
  end subroutine calc_continuum

  subroutine h2o_continuum()
    ! calculates the continuum absorption for a give atmoshere
    ! It is a wrapper for the MT-CKD continuum and sets up the variables for

    implicit none
    integer iband,k, kk, mxone, nmon, kvert, nptabs
    real (double) :: wtot, wa, wn2, vmrh2o,w_dry
    real (double) :: V1ABS,V2ABS,DVABS,ABSRB,xlength
    real (double) :: XID,SECANT,PAVE,TAVE,HMOLID,XALTZ
    real (double) :: WK,PZL,PZU,TZL,TZU,WBROAD,DV,V1 ,V2 ,TBOUND
    real (double) :: EMISIV,FSCDID,NMOL_C,LAYER ,YI1,YID,LSTWDF
    real (double) :: XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL 
    COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
    COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
         &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
         &                EMISIV,FSCDID(17),NMOL_C,LAYER ,YI1,YID(10),LSTWDF
      common /cntscl/ XSELF(7),XFRGN(7),XCO2C,XO3CN,XO2CN,XN2CN,XRAYL 

    kvert = nspec + 1
    mxone = 1
    do iband = 1, nband
       do k = 1, kmax-1
          pave = p(k)*1013.15d0 ! convert in mbar
          tave = t(k)
          xlength = abs((z(k+1) - z(k))*10000.0d0) !thickness needed in cm

          vmrh2o = xgas(1,k)
          wtot = xgas(1,k)*ccc(kvert,k)
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('H2O')) then
                vmrh2o = x(kk,k)
                wtot = x(kk,k)*ccc(kvert,k)
             end if
          end do
          
          w_dry = wtot * (1.-vmrh2o)
          wk(1) = vmrh2o * w_dry
          !ARGON
          WA     = 0.009     * W_dry 
          
          !NITROGEN
          wn2 = xgas(41,k) * ccc(kvert,k)
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('N2')) then
                wn2 = x(kk,k)*ccc(kvert,k)
             end if
          end do

          ! CO2
          wk(2) = xgas(2,k)*ccc(kvert,k)
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('CO2')) then
                wk(2) = x(kk,k)*ccc(kvert,k)
             end if
          end do

          ! Ozone
          wk(3) = xgas(3,k)*ccc(kvert,k)
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('O3')) then
                wk(3) = x(kk,k)*ccc(kvert,k)
             end if
          end do

          ! Oxygen
          wk(7) = xgas(7,k)*ccc(kvert,k)
          
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('O2')) then
                wk(7) = x(kk,k)*ccc(kvert,k)
             end if
          end do

          wbroad=wn2+wa
          nmol_c = 7
             
          nmon = nm(iband)
          v1abs = wstart(iband)
          v2abs = wstart(iband) + nmon*spac(iband)
          dvabs = spac(iband)
          nptabs = nmon
          v1 = v1abs
          v2 = v2abs
          xself(1:7) = 1.0d0
!          print *, pave, tave, xlength, wk(1), wk(2), wk(3), wk(7)
          call contnm(1)

!          print *, absrb(1:nmon)

          mtckd(1, k, mxone:mxone+nmon) = absrb(1:nmon)
          mxone = mxone + nmon
       end do
    end do
    
  end subroutine h2o_continuum
  
end module continuum
