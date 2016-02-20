program convert

  ! corrects matrices to have one row in the matrix in one row of the file
  
  implicit none
  character (len=10000) :: tmp, read_format, write_format
  integer :: ind, nr_row, nr_column, itmp, ind2, itmp2
  real(kind=8), allocatable, dimension(:) :: row 

  
!  call system ('mv g.out g.out.orig')
  open(10, file = 'g.out.orig', status='old')
  open(20, file = 'g.out.new', status='unknown')

  read(10, *) tmp! skip first line
  write(20,'(A)') trim(tmp)
  read(10, *), nr_row, nr_column, itmp, itmp2
  write(20, *) nr_row, nr_column, itmp, itmp2
  read(10, '(A)') tmp! skip third line
  write(20, *) trim(tmp)! skip third line

  allocate(row(nr_column))
  
  write(write_format, '(a,i5,a)') '(', nr_column, 'ES26.18)'  
  do ind = 1,nr_row
     read(10,*) row
     write(20,trim(write_format)) row
  end do

  deallocate(row)
  close(10)
  close(20)
  
end program convert
