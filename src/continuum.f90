module continuum

  use params
  use xsections

  implicit none
  
  integer :: ncont = 0, n_contabs
  logical :: f_contabs = .false.
  integer :: abscont_type
  real(double), dimension(2) :: abscont_param
  real(double), dimension(:), allocatable :: cont_param
  
  contains

  subroutine calc_continuum(param)
    ! wrapper for calculation of the continua. The continuums absorption gets
    ! inserted into the cross sections after the gas crosssections

    implicit none

    real(double), dimension(:), intent(in) :: param
    

    integer :: iband, k, j
    integer :: mone, mxne
    real(double) :: wone, wxne, wmid
    

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
                CROSS(nret+2,K,j) = CROSS(nret+2,K,j) * ( 1 + param(2)*((wone+dble(j)*dn(iband))-wmid)/(wxne-wone) )
             end do
          end DO
          mone = mone + nm(iband)
       end DO
    end select
  end subroutine calc_continuum
  
end module continuum
