program convert

  ! convert linelists in proffit format to sfit4 format
  
  implicit none



  logical :: declm, decLM2, decshape
  integer :: hitspecies
  real :: hitnue,hitkappacm,hitdummy,hitlorwidthf,hitlorwidths,hitergniveau,hitlortdepend,hitpshift
  real :: auxbeta, auxbetatdepend, auxgam2, auxgam2tdepend, auxylm, auxlmtk1, auxlmtk2
  real :: auxylm2, auxlm2tk1, auxlm2tk2
  real :: LMY0f,LMY1f,LMY2f,LMY0s,LMY1s,LMY2s,LMI0f,LMI1f,LMI0s,LMI1s,LMD0f,LMD1f,LMD0s,LMD1s
  character(len=15) :: quanta1, quanta2, quanta3, quanta4
  character(len=160):: hitran_line, tmp


  open(10, file = '06_kit15_2015-12-08_LineFit.par')
  open(20, file = 'profit_HITRAN.dat', status='unknown')
  open(30, file = 'profit_LM1ST.dat', status='unknown')
  open(40, file = 'profit_LM2RD.dat', status='unknown')

  declm=.false.
  declm2=.false.
  do 
     read (10,'(a160)',end = 10) hitran_line
     read (hitran_line, '(I3,A157)') hitspecies,tmp
     if (hitspecies.eq.995) then !Rosenkrantz line mixing
        read (hitran_line,803) hitspecies,hitnue,LMY0f,LMY1f,LMY2f,LMY0s,LMY1s,LMY2s
        decLM = .true.
        auxLMTK1 = LMY1f/LMY0f
        auxLMTK2 = LMY2f/LMY0f
        auxYLM = LMY0f
     end if
     if (hitspecies.eq.994) then !Smith line mixing
        read (hitran_line,803) hitspecies,hitnue,LMI0f,LMI1f,LMI0s,LMI1s,LMD0f,LMD1f,LMD0s,LMD1s
        decLM2 = .true.
        auxLM2TK1 = LMI1f
        auxLM2TK2 = LMI0s
        auxYLM2 = LMI0f
     end if
     if (hitspecies.lt.900) then
        read(hitran_line, 801) hitspecies,hitnue,hitkappacm,hitdummy &
             ,hitlorwidthf,hitlorwidths,hitergniveau,hitlortdepend,hitpshift,quanta1, quanta2, quanta3, quanta4, tmp
        if (declm) then
           write(30,701) hitspecies, quanta1, quanta2, quanta3, quanta4, auxlmtk1, auxlmtk2, auxylm
           declm = .false.
        end if
        if (declm2) then
           write(40,701) hitspecies, quanta1, quanta2, quanta3, quanta4, auxlm2tk1, auxlm2tk2, auxylm2
           declm2 = .false.
        end if
        write(20,'(a160)') hitran_line
     end if
  end do
  
10 continue
  close (10)
  close (20)
  close (30)
  close (40)
701 format(i3,a15,a15,a15,a15,E10.3,1x,E10.3,1x,E10.3,1x,E10.3,1x,E10.3,1x,E10.3,1x,E10.3)
801 format(I3,F12.6,E10.3,E10.3,F6.5,F4.3,F10.4,F4.2,F8.6,a15,a15,a15,a15,a160)  ! HIT 08 format
803 format(I3,F12.6,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3,E10.3)! AUX format in PROFFIT linelist
end program convert
