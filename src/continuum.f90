module continuum

  use params
  use xsections

  implicit none
  
  integer, parameter :: cont_poly_max = 10

  integer :: ncont = 0, n_contabs
  logical :: f_contabs = .false.
  integer :: abscont_type, abscont_order
  real(double), dimension(cont_poly_max) :: abscont_param, abscont_sparam
  real(double), dimension(:), allocatable :: cont_param
  
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

    abscont_type = 2
    select case (abscont_type) 
    case (0)
       ! Offset only
       !  --- Loop over layers
       DO K = 1, KMAX
          CROSS(nret+2,K,:NCROSS) = param(1)*(P(k)/P(kmax))
       end DO
    case (1)
       ! Offset and slope
       mone = 1
       DO IBAND = 1, NBAND
          mxne = mone + nm(iband) - 1
          wone = wstart(iband)
          wxne = wstart(iband) + dn(iband)*nm(iband)
          wmid = (wone + wxne)/2.0d0
!          print *, mone, mxne, wone, wxne, wmid, dn(iband)
!          print *, param(1), param(2), abscont_strength, abscont_tilt
          DO K = 1, KMAX     
             CROSS(nret+2,K,mone:mxne) = param(1)*(P(k)/P(kmax))
             do j = mone,mxne
                CROSS(nret+2,K,j) = CROSS(nret+2,K,j) * ( 1.0D0 + param(2)*((wone+dble(j)*dn(iband))-wmid)/(wxne-wone) )
             end do
          end DO
          mone = mone + nm(iband)
       end DO
    case (2)
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
                polynom = 0.0d0
                do l = 0,n_contabs-1
                   polynom = polynom + param(l+1)*(((wone+dble(j)*dn(iband))-wmid)/(wxne-wone))**l
                end do
                CROSS(nret+2,K,j) = CROSS(nret+2,K,j) + polynom*(P(k)/P(kmax))
             end do
          end DO
          mone = mone + nm(iband)
       end DO       
    end select
  end subroutine calc_continuum
  
end module continuum
