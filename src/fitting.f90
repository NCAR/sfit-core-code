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

module fitting

    !-------------
    ! Modules used
    !-------------
    use matrix       ! Uncomment this to use matrix inversion procedure in matrix.f90

    implicit none

    !--------------------------------
    ! Set kind type for use in module
    !--------------------------------
    integer,parameter :: rd = selected_real_kind(15)
    integer,parameter :: id = selected_int_kind(15)
    integer,parameter :: is = selected_int_kind(8)


    !----------------------
    ! Set member visibility
    !----------------------
    private     ! Default all to private
    private      :: arith_mean
    private      :: arith_mean_pure
    private      :: lst_sqrs_poly_fun
    private      :: lst_sqrs_poly_sub
    public       :: polyfit

contains

   function polyfit(vx, vy, npnts, norder)

      implicit none

      !integer, parameter                    :: dp = selected_real_kind(15, 307)
      integer, intent(in)                 :: norder, npnts ! polynomial order, there are +1 coefficients including constant
      real(8), dimension(:), intent(in)   :: vx, vy
      real(8), dimension(norder+1)        :: polyfit

      integer                             :: ncoeff, mode
      real(8)                             :: r_sqrd, rmse !, chisq
      real(8), dimension(:), allocatable  :: sigmay, a

      if( npnts .ne. size(vx) .and. npnts .ne. size(vy) )stop ' npnts & size of arrays vx & vy need to be the same.'

      mode   = 0
      ncoeff = norder + 1

      allocate( sigmay(npnts), a(ncoeff) )

      a      = 0.0d0
      sigmay = 0.0d0

      call lst_sqrs_poly_sub( vx, vy, npnts, ncoeff, a, r_sqrd, rmse )

      !call polfit( vx, vy, sigmay, npts, nterms, mode, a, chisq )

      polyfit(:ncoeff) = a(:ncoeff)

      !polyfit = polyfitlp(vx, vy, d)

      deallocate( sigmay, a )

      return

   end function polyfit


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Name:
!       Math_functions
!
! Purpose:
!       Contains simple math functions that are not intrinsic to fortran
!
!
!
! Calling Sequence:
!       USE math_functions
!
!
! Modules:
!       kind_type:          Module containing definitions for variable kind type (NOT USED ANYMORE)
!
! Contains:
!       arith_mean
!       -------------
!       Function to calculate the Arithmetic mean of a rank 1 array
!
!       lst_sqrs_poly_fun
!       -------------
!       Function to calculate the coefficients of a polynomial fit using least squares
!
!       lst_sqrs_poly_sub
!       -------------
!       Subroutine to calculate the coefficients of a polynomial fit using least squares
!
!       inverse (Private)
!       -------------
!       Subroutine to find the inverse of a square nxn matrix
!
!
!
! Notes:
!       1) The least squares polynomial fit can either use the matrix inversion procedure
!          included in this module or the matrix procedure included in matrix.f90 depending
!          on what is commented out.
!
!
! Version History:
!   --September, 2013  => Commented out nearestpoint_v1 and linear_interp1D as these depend on intrinsic ifort compiler
!                         routine. Added lst_sqrs_fit which does not call LAPACK.
!   --April, 2011 => Added procedure overloading to allow for scalar or array inputs
!   --March, 2011 => Created by Eric Nussbaumer (ebaumer@atmos.umd.edu)
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!--------------------------------------------------------------------------------------------------------------
! Name:
!       lst_sqrs_poly_fun  (Function version)
!
! Purpose:
!       This program fits a polynomial of one variable given discrete data
!       Y = coeff(1) + coeff(2) * X + coeff(3) * X^2 + ... + coeff(norder+1) * X^norder
!
!
! Calling Sequence:
!    coeff = lst_sqrs_poly( x,             & ! Input
!                           y,             & ! Input
!                           norder)        & ! Input
!
!
! Input Arguments:
!
!          x:               x values of data points
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-1 array
!                           ATTRIBUTES:  Intent( in )
!
!          y:               y values of data points
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-1 array
!                           ATTRIBUTES:  Intent( in )
!
!      nordr:               The order of the fit
!                           UNITS:       N/A
!                           TYPE:        int( rd )
!                           DIMENSION:   Scalar
!                           ATTRIBUTES:  Intent( in )
!
!
!
!
! Output arguments (Function Result):
!
!      coeff:               Calculated coefficients
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-1 array
!                           ATTRIBUTES:  Intent( out )
!
!
!
!
! Calls:
!       inverse -- Procedure to calculate the inverse of a matrix (local procedure)
!       INVRT   -- Procedure to calculate the inverse of a matrix (from matrix.f90)
!
! Notes:
!       1) The least squares polynomial fit can either use the matrix inversion procedure
!          included in this module or the matrix procedure included in matrix.f90 depending
!          on what is commented out.
!
!
! Procedure:
!
!
! Reference:
!
!
! Version History:
!   --September, 2013 => Created by Eric Nussbaumer (ebaumer@ucar.edu)
!
!-----------------------------------------------------------------------------------------------------------------
    function lst_sqrs_poly_fun( x,             & ! Input
                                y,             & ! Input
                                norder)        & ! Input
                                result(coeff)    ! Output


    !--------------------------------------------------------------------------
    !                         -- TYPE DECLARATIONS --
    !--------------------------------------------------------------------------
    implicit none

    ! Input variables
    real(rd),    intent(in)               :: x(:), y(:)
    integer(id), intent(in)               :: norder

    ! Output variables
    real(rd),dimension(:)                 :: coeff(norder+1)

    ! Local Parameters
    real(rd),dimension(:,:),allocatable   :: A, AT, ATA, ATA_inv
    integer(id)                           :: npnts, i, j
    !integer(id)                          :: ncoeff
    integer(is)                           :: ncoeff


    !--------------------------------------------------------------------------
    !                         -- INITIALIZE VARS --
    !--------------------------------------------------------------------------
    npnts  = size(x)
    ncoeff = size(coeff)

    if ( size(x) /= size(y) ) then
        write(0,*) 'Input vectors x and y are not the same size!!!'
        stop
    endif

    allocate( A(npnts,ncoeff)        )
    allocate( AT(ncoeff,npnts)       )
    allocate( ATA(ncoeff,ncoeff)     )
    allocate( ATA_inv(ncoeff,ncoeff) )

    coeff   = 0.0_rd


    !--------------------------------------------------------------------------
    !                           -- Calculations --
    !--------------------------------------------------------------------------
    forall ( i = 0:norder, j = 1:npnts )
        A(j,i+1) = x(j)**i
    end forall

    AT  = transpose(A)
    ATA = matmul( AT, A )


    !--------------------
    ! Find inverse of ATA
    !--------------------
    !call inverse( ATA, ncoeff, ATA_inv  )    ! This is matrix inversion call for local procedure

    call INVRT( ATA, ATA_inv, ncoeff )        ! This is matrix inversion call for procedure in matrix.f90


    !----------------------
    ! Find fit coefficients
    !----------------------
    coeff = matmul( matmul( ATA_inv, AT ), y )

    !---------------
    ! Release memory
    !---------------
    deallocate(ATA)
    deallocate(AT)
    deallocate(A)

    end function lst_sqrs_poly_fun


!--------------------------------------------------------------------------------------------------------------
! Name:
!       lst_sqrs_poly_sub  (Subroutine version)
!
! Purpose:
!       This program fits a polynomial of one variable given discrete data
!       Y = coeff(1) + coeff(2) * X + coeff(3) * X^2 + ... + coeff(norder+1) * X^norder
!
!
! Calling Sequence:
!    call  subroutine lst_sqrs_poly_sub( x,             & ! Input
!                                        y,             & ! Input
!                                        norder,        & ! Input
!                                        coeff,         & ! Output
!                                        r_sqrd,        & ! Output
!                                        rmse )           ! Output
!
!
! Input Arguments:
!
!          x:               x values of data points
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-1 array
!                           ATTRIBUTES:  Intent( in )
!
!          y:               y values of data points
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-1 array
!                           ATTRIBUTES:  Intent( in )
!
!       npnts:              The number of data points in x & y
!                           UNITS:       N/A
!                           TYPE:        int( is )
!                           DIMENSION:   Scalar
!                           ATTRIBUTES:  Intent( in )
!
!      ncoeff:              The polynomial order + 1 of the fit
!                           UNITS:       N/A
!                           TYPE:        int( is )
!                           DIMENSION:   Scalar
!                           ATTRIBUTES:  Intent( in )
!
!
!
!
! Output arguments:
!
!      coeff:               Calculated coefficients
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-1 array
!                           ATTRIBUTES:  Intent( out )
!
!      r_sqrd:              R^2 value of fit
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Scalar
!                           ATTRIBUTES:  Intent( out )
!
!      rmse:                Root mean square error of fit
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Scalar
!                           ATTRIBUTES:  Intent( out )
!
!
!
! Calls:
!       inverse -- Procedure to calculate the inverse of a matrix (local procedure)
!       INVRT   -- Procedure to calculate the inverse of a matrix (from matrix.f90)
!
!
! Notes:
!       1) The least squares polynomial fit can either use the matrix inversion procedure
!          included in this module or the matrix procedure included in matrix.f90 depending
!          on what is commented out.
!
!
! Procedure:
!
!
! Reference:
!
!
! Version History:
!   --September, 2013 => Created by Eric Nussbaumer (ebaumer@ucar.edu)
!
!-----------------------------------------------------------------------------------------------------------------
    subroutine lst_sqrs_poly_sub( x,             & ! Input
                                  y,             & ! Input
                                  npnts,         & ! Input size of x & y
                                  ncoeff,        & ! Input order of polynomial
                                  coeff,         & ! Output polynomial coeficients
                                  r_sqrd,        & ! Output
                                  rmse )           ! Output


    !--------------------------------------------------------------------------
    !                         -- TYPE DECLARATIONS --
    !--------------------------------------------------------------------------
    implicit none

    ! Input variables
    real(rd),    intent(in)               :: x(:), y(:)
    integer(is), intent(in)               :: npnts, ncoeff

    ! Output variables
    real(rd),dimension(:),intent(out)     :: coeff(ncoeff)
    real(rd), intent(out)                 :: r_sqrd, rmse

    ! Local Parameters
    real(rd)                              :: y_bar
    real(rd),dimension(:,:),allocatable   :: A, AT, ATA, ATA_inv
    real(rd),dimension(:),  allocatable   :: y_calc, s_tot, s_res
    integer(id)                           :: i, j
    !integer(id)                          :: ncoeff

    !--------------------------------------------------------------------------
    !                         -- INITIALIZE VARS --
    !--------------------------------------------------------------------------
    !npnts  = size(x)
    !ncoeff = size(coeff)

    if ( size(x) /= size(y) ) then
        write(0,*) 'Input vectors x and y are not the same size!!!'
        stop
    endif

    allocate( A(npnts,ncoeff)        )
    allocate( AT(ncoeff,npnts)       )
    allocate( ATA(ncoeff,ncoeff)     )
    allocate( ATA_inv(ncoeff,ncoeff) )

    allocate( s_tot(npnts) )
    allocate( s_res(npnts) )
    allocate( y_calc(npnts))

    coeff   = 0.0_rd
    r_sqrd  = 0.0_rd
    y_bar   = 0.0_rd
    y_calc  = 0.0_rd


    !--------------------------------------------------------------------------
    !                           -- Calculations --
    !--------------------------------------------------------------------------
    forall ( i = 0:ncoeff-1, j = 1:npnts )
        A(j,i+1) = x(j)**i
    end forall

    AT  = transpose(A)
    ATA = matmul( AT, A )


    !--------------------
    ! Find inverse of ATA
    !--------------------
    !call inverse( ATA, ncoeff, ATA_inv  )          ! This is matrix inversion call for local procedure

    call INVRT( ATA, ATA_inv, ncoeff )              ! This is matrix inversion call for procedure in matrix.f90

    !----------------------
    ! Find fit coefficients
    !----------------------
    coeff = matmul( matmul( ATA_inv, AT ), y )


    !----------------------------------------------------
    ! Calculate R-squared or Coefficient of Determination
    !----------------------------------------------------
    y_bar = arith_mean( y )

    do i = 1,npnts
        do j = 0,ncoeff-1
            y_calc(i) = y_calc(i) + ( coeff(j+1) * x(i)**j )
        end do
        s_tot(i) = ( y(i) - y_bar )**2
        s_res(i) = ( y_calc(i) - y(i) )**2
    end do

    r_sqrd = 1.0_rd - ( sum(s_res) / sum(s_tot) )


    !---------------------------------
    ! Calculate Root mean square error
    !---------------------------------
    rmse = sqrt( sum(s_res) / npnts )


    !---------------
    ! Release memory
    !---------------
    if( allocated( ATA ))deallocate(ATA)
    if( allocated( AT ))deallocate(AT)
    if( allocated( A ))deallocate(A)
    if( allocated( ATA_inv ))deallocate(ATA_inv)

    if( allocated( s_tot ))deallocate(s_tot)
    if( allocated( s_res ))deallocate(s_res)
    if( allocated( y_calc ))deallocate(y_calc)

    end subroutine lst_sqrs_poly_sub



!--------------------------------------------------------------------------------------------------------------
! Name:
!       arith_mean
!
! Purpose:
!       Calculate the weighted arithmetic mean of a rank 1 array. If weights are one or not included then the
!       the straight arithmetic mean is calculated
!
!
! Calling Sequence:
!       x_bar =             arith_mean( x,             & ! Input
!                                       weights,       & ! Optional Input
!                                       message_log  )   ! Error Messaging (Optional)
!
!
! Input Arguments:
!       x:                  values to be averaged
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-1 array
!                           ATTRIBUTES:  Intent( in )
!
!       weights:            Weights representing the reliability of the influence upone the mean by
!                           each respective values
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-1 array
!                           ATTRIBUTES:  Intent( in ), Optional
!
!
!
! Output arguments (Function Result):
!       x_bar:             The arithmatic mean of x
!                          UNITS:      Same as X
!                          TYPE:       real( rd )
!                          DIMENSION:  Scalar
!
! Calls:
!
! Notes:
!   -- There are two versions of this function 1) arith_mean_pure which is declared as a pure function for use in
!                                                 FORALL constructs
!                                              2) arith_mean which is the regular version.
!
!
! Procedure:
!       The weighted arithmetic mean is defined as:
!
!
!
!
!       When weights are not present or = 1 then the arithmetic mean is calculated by:
!
!           _      1
!           X  =  --- * sum( x_i, i=1,N )
!                  N
!
!
!
! Version History:
!   --March, 2011 => Created by Eric Nussbaumer (ebaumer@atmos.umd.edu)
!
!
!-----------------------------------------------------------------------------------------------------------------
!====================================================================================
!                           Pure Version -for use in forall constructs
!====================================================================================

    pure function arith_mean_pure( x,             & ! Input
                                   weights)       & ! Input,           Optional
                                   result(x_bar)    ! Output

    !#--------------------------------------------------------------------------#
    !#                         -- TYPE DECLARATIONS --                          #
    !#--------------------------------------------------------------------------#

    implicit none
    ! Input variables
    real(rd),intent(in)               :: x(:)
    real(rd),optional,intent(in)      :: weights(:)

    ! Function results
    real(rd)                          :: x_bar

    ! Local Parameters
    integer(id)                       :: N


    !#--------------------------------------------------------------------------#
    !#                     -- INITIALIZE RETURN VALUE --                        #
    !#--------------------------------------------------------------------------#

    x_bar = 0.0_rd


    !#--------------------------------------------------------------------------#
    !#                           -- CHECK INPUT --                              #
    !#--------------------------------------------------------------------------#

    N = size(x)


    !#--------------------------------------------------------------------------#
    !#                          -- Calculate Mean --                            #
    !#--------------------------------------------------------------------------#

    if ( present( weights ) ) then

        x_bar = sum( weights * x ) / sum( weights )

    else

        x_bar = sum( x ) / N

    endif
    end function arith_mean_pure



!====================================================================================
!                           Regular Version
!====================================================================================

    function arith_mean( x,             & ! Input
                         weights)       & ! Input,           Optional
                         result(x_bar)    ! Output

    !#--------------------------------------------------------------------------#
    !#                         -- TYPE DECLARATIONS --                          #
    !#--------------------------------------------------------------------------#

    implicit none
    ! Input variables
    real(rd),intent(in)               :: x(:)
    real(rd),optional,intent(in)      :: weights(:)

    ! Function results
    real(rd)                          :: x_bar

    ! Local Parameters
    integer(id)                       :: N



    !#--------------------------------------------------------------------------#
    !#                     -- INITIALIZE RETURN VALUE --                        #
    !#--------------------------------------------------------------------------#

    x_bar = 0.0_rd


    !#--------------------------------------------------------------------------#
    !#                           -- CHECK INPUT --                              #
    !#--------------------------------------------------------------------------#

    N = size(x)

!     if ( N <= 1 ) then
! 		print*, 'In arith_mean: Must have > 1 inputs to calculate mean'
!         return
!
!     elseif ( size( x ) /= size( weights ) ) then
! 		print*, 'In arith_mean: X and weights must have the same size'
!         return
!
!     endif


    !#--------------------------------------------------------------------------#
    !#                          -- Calculate Mean --                            #
    !#--------------------------------------------------------------------------#

    if ( present( weights ) ) then

        x_bar = sum( weights * x ) / sum( x )

    else

        x_bar = sum( x ) / N

    endif
    end function arith_mean


!--------------------------------------------------------------------------------------------------------------
! Name:
!       inverse
!
! Purpose:
!       Calculates the inverse of a nxn matrix based on Doolittle LU factorization for Ax=b
!
!
! Calling Sequence:
!       call    inverse( a, &     ! Input
!                        n, &     ! Input
!                        a_inv  ) ! Output
!
!
! Input Arguments:
!       a:                  Matrix (nxn) to be inverted
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   Rank-2 matrix
!                           ATTRIBUTES:  Intent( in )
!
!       n:                  Dimension of matrix to be inverted
!                           UNITS:       N/A
!                           TYPE:        real( rd )
!                           DIMENSION:   scalar
!                           ATTRIBUTES:  Intent( in )
!
!
!
! Output arguments (Function Result):
!       a_inv:             Inverse of matrix A
!                          UNITS:      N/A
!                          TYPE:       real( rd )
!                          DIMENSION:  Rank-2 matrix
!                          ATTRIBUTES:  Intent( Out )
!
! Calls:
!
! Notes:
!   Original matrix a is destroyed during calculation
!
!
! Procedure:
!
!                  N
!
!
!
! Version History:
!   --September, 2013 => Created by Eric Nussbaumer (ebaumer@ucar.edu) based on code from
!                        Alex G. (http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90)
!
!-----------------------------------------------------------------------------------------------------------------
  subroutine inverse( a, &     ! Input
                      n, &     ! Input
                      a_inv  ) ! Output


    !--------------------------------------------------------------------------
    !                         -- TYPE DECLARATIONS --
    !--------------------------------------------------------------------------
    implicit none
    integer(is),intent(in)   :: n
    real(rd),intent(out)     :: a_inv(n,n)
    real(rd)                 :: a(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
    real(rd)                 :: coeff
    integer(id)              :: i, j, k

    !--------------------------------------------------
    ! step 0: initialization for matrices L and U and b
    !--------------------------------------------------
    L=0.0_rd
    U=0.0_rd
    b=0.0_rd

    !----------------------------
    ! step 1: forward elimination
    !----------------------------
    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
        end do
    end do

    !----------------------------------------------------
    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    !----------------------------------------------------
    do i=1,n
        L(i,i) = 1.0_rd
    end do

    !-------------------------------------------
    ! U matrix is the upper triangular part of A
    !-------------------------------------------
    do j=1,n
        do i=1,j
            U(i,j) = a(i,j)
        end do
    end do

    !------------------------------------------------
    ! Step 3: compute columns of the inverse matrix C
    !------------------------------------------------
    do k=1,n
        b(k)=1.0_rd
        d(1) = b(1)
        !---------------------------------------------------
        ! Step 3a: Solve Ld=b using the forward substitution
        !---------------------------------------------------
        do i=2,n
            d(i)=b(i)
            do j=1,i-1
                d(i) = d(i) - L(i,j)*d(j)
            end do
        end do
        !------------------------------------------------
        ! Step 3b: Solve Ux=d using the back substitution
        !------------------------------------------------
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
                x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
        end do
        !----------------------------------------------------
        ! Step 3c: fill the solutions x(n) into column k of C
        !----------------------------------------------------
        do i=1,n
            a_inv(i,k) = x(i)
        end do
        b(k)=0.0_rd
    end do

    end subroutine inverse

end module fitting

