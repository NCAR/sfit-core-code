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

  
end module continuum
