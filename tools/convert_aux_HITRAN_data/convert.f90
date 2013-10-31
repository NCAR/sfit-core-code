program convert

  implicit none



  logical :: declm, decshape
  integer :: hitspecies
  real :: hitnue,hitkappacm,hitdummy,hitlorwidthf,hitlorwidths,hitergniveau,hitlortdepend,hitpshift
  real :: auxbeta, auxbetatdepend, auxgam2, auxgam2tdepend, auxylm, auxlmtk1, auxlmtk2
  character(len=15) :: quanta1, quanta2, quanta3, quanta4
  character(len=160):: hitran_line, tmp


  open(10, file = '02_hit08_f53+LM_s2.par')
  open(20, file = '02_hase_lm.dat')

  do 
     read (10,'(a160)',end = 10) hitran_line
     read(hitran_line, 801) hitspecies,hitnue,hitkappacm,hitdummy &
          ,hitlorwidthf,hitlorwidths,hitergniveau,hitlortdepend,hitpshift,quanta1, quanta2, quanta3, quanta4, tmp
    
     if (declm) then
        write(20,701) hitspecies, quanta1, quanta2, quanta3, quanta4, auxlmtk1, auxlmtk2, auxylm
        declm = .false.
     end if
     
     if (hitspecies .eq. 997) then  ! 997 line mixing
        decLM = .true.
        auxLMTK1 = hitkappacm
        auxLMTK2 = hitdummy
        auxYLM = hitergniveau
     else if (hitspecies .eq. 998) then  ! 998 SDV
        decshape = .true.
        auxgam2 = hitkappacm
        auxgam2tdepend = hitdummy
     else if (hitspecies .eq. 999) then  ! 999 Galatry
        decshape = .true.
        auxbeta = hitkappacm
        auxbetatdepend = hitdummy
     end if
     
  end do
  
10 continue
  close(10)
  close(20)
  
701 format(i3,a15,a15,a15,a15,e10.4,1x,e10.4,1x,e10.4)
801 format(I3,F12.6,E10.3,E10.3,F6.5,F4.3,F10.4,F4.2,F8.6,a15,a15,a15,a15,a160)  ! HIT 08 format
end program convert
