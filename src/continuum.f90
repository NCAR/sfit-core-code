module continuum

  use params
  use xsections

  implicit none
  
  integer, parameter :: cont_poly_max = 10

  integer :: ncont = 0, n_contabs
  logical :: f_continuum = .false.
  integer :: abscont_type, abscont_order
  real(double), dimension(cont_poly_max) :: abscont_param, abscont_sparam
  real(double), dimension(:), allocatable :: cont_param
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
       case(1)
          ! n-th order polynomial
          mone = 1
          CROSS(nret+2,:KMAX,:ncross) = 0.0d0
          DO IBAND = 1, NBAND
             mxne = mone + nm(iband) - 1
             wone = wstart(iband)
             wxne = wstart(iband) + dn(iband)*nm(iband)
             wmid = (wone + wxne)/2.0d0
             DO K = 1, KMAX
                IF ((abscont_type.ne.2).or.((k.lt.kmax).and.(z(k).ge.cont_z_abs).and.(z(k+1).le.cont_z_abs))) then
                   print *, k, param(1)
                   do j = mone,mxne
                      polynom = param(1)
                      do l = 1,n_contabs-1
                         polynom = polynom + param(l+1)*(((wone+dble(j)*dn(iband))-wmid)/(wxne-wone))**l
                      end do
                      
                      CROSS(nret+2,K,j) = CROSS(nret+2,K,j) + polynom*(P(k)/P(kmax))
                   end do
                end IF
             end DO
          end DO
       case(2)
          ! an absorbing layer modeled at the altitude
          ! z_abs with the absorbing strength cont_alpha (can be retrieved)
          mone = 1
          CROSS(nret+2,:KMAX,:ncross) = 0.0d0
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
                      CROSS(nret+2,K,j) = CROSS(nret+2,K,j) + polynom
                   end do
                end IF
             end DO
             mone = mone + nm(iband)
          end DO
       end select
    end if
  end subroutine calc_continuum

  
end module continuum
