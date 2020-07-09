      integer*4
     & lunr,    ! Logical Unit Number for opening "molparam.dat"
     & kgas,    ! Gas number (H2O=1, CO2=2, O3=3, ....)
     & iso,     ! isotopomer number (most abundant=1,...)
     & mspeci,  ! Maximum number of different species
     & nspeci,  ! Actual number of different species
     & jspeci,  ! Species index
     & nmode,   ! Number of different vibrational modes (ie 3N-6) of each gas
     & mmode,   ! Maximum possible value of NMODE
     & j        ! scratch variable

      parameter (mmode=20,mspeci=110,lunr=14)

      integer*4
     & dgen(mmode)     ! Degeneracy of vibrational mode.

      real*4
     & fia,            ! Fractional Isotopomeric Abundnce
     & mol_wt(mspeci), ! Molecular Weight
     & tdrpf(mspeci),  ! Temperature Dependence of Rotational Partition Function
     & vibfrq(mmode)   ! Vibrational frequencies

      character
     & gas*8,            ! Name of the gas (e.g. H2O, CO2, O3, ....)
     & speci(mspeci)*24  ! the chemical formula of the species

      open(unit=lunr,file='isotopomers.dat',status='old')
      do jspeci=1,mspeci
        read(lunr,77,end=99)kgas,iso,gas,speci(jspeci),fia,
     &  mol_wt(jspeci),tdrpf(jspeci),nmode,(vibfrq(j),dgen(j),j=1,nmode)
 77     format(2i2,1x,a8,a24,f10.8,i3,f5.2,i3,20(f6.0,i2))
        if(nmode.gt.mmode) stop ' Increase parameter MMODE'
      end do
      write(*,*)'Increase parameter mspeci'
 99   close(lunr)
      nspeci=jspeci-1
      write(*,*)nspeci
      stop
      end
