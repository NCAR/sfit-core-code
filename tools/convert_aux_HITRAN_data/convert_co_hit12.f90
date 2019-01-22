program convert_devy

! Convert CO line mixing and sdv data as given by Devy 2013 to the format
! read by hbin. In particular, convert the speed dependency value by division
! by gamma_0

  implicit none



  character(len=255) :: line
  character(len=60) :: quanta
  integer :: mo, is 
  real :: g0_air,g0_self,td_g0_air,s_air, g2, lm1st

  open(10, file = '05_hit12_SDV.orig')
  open(20, file = '05_hit12_LM1ST.txt')
  open(30, file = '05_hit12_SDV.txt')

  read (10,'(a255)',end=10) line
  do 
     read (10,'(a255)',end=10) line
     read( line, '(i2,i1,a60)' ) mo,is, quanta
     read (line(64:), '(f6.5,f4.3,f3.2,f8.6,f8.6)') g0_air,g0_self,td_g0_air,s_air
     read (line(85:), *, end=26) g2
     write(30, '(i2,i1,a60,f6.5,f5.4,f3.2,f8.6,1x,g12.6)') mo,is,quanta,g0_air,g0_self,td_g0_air,s_air,g2*g0_air
26   continue
     read (line(100:), *, end=27) lm1st
     write(20, '(i2,i1,a60,1x,e10.3,1x,e10.3,1x,e10.3)')mo,is, quanta,0.0d0, 0.0d0, lm1st
27   continue
  end do
  
10 continue
  close(10)
  close(20)
  close(30)

end program convert_devy
