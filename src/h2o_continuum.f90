module h2o_continuum

  use params
  use xsections

  implicit none
  
  real(double), dimension(:,:,:), allocatable :: mtckd
  logical :: f_mtckd = .false.
  
contains


  subroutine calc_h2o_continuum()
    ! calculates the continuum absorption for a give atmoshere
    ! It is a wrapper for the MT-CKD continuum and sets up the variables for
    
    
    
    
    real(double) :: V1ABS,V2ABS,DVABS, ABSRB
    integer :: NPTABS, NMOL_C,LAYER ,LSTWDF,NPTC,NPTh
    integer :: IRD,IPRcnt,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL
    integer :: NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL  
    integer :: NLTEFL,LNFIL4,LNGTH4
    integer :: i,k,kk,icflg,jrad,mxone, nmon, kvert, ksmax2,iband
    
    
    real(double) :: PAVE,TAVE
    real(double) :: WK,PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND
    real(double) :: EMISIV,FSCDID,YI1
    real(double) :: V1C,V2C,DVC,C
    real(double) :: V1h,V2h,DVh,Ch,csh2o,cfh2o
    real(double) :: XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
    real(double) :: vi,vmrh2o,W_dry,wa,wn2,wtot,xcnt,xlength
    
    CHARACTER*18 HNAMCNT,HVRCNT
    !                                                                         F00100
    CHARACTER*8      XID,       HMOLID,      YID 
    REAL*8               SECANT,       XALTZ
    
    COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)              ! 500060
    !
    COMMON /CVRCNT/ HNAMCNT,HVRCNT
    !
    COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),  &    ! F00130
         WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &  ! F00140
         EMISIV,FSCDID(17),NMOL_C,LAYER ,YI1,YID(10),LSTWDF   ! F00150
    !
    !
    COMMON /XCONT/  V1C,V2C,DVC,NPTC,C(IPTS2) 
    !
    !********************************************
    COMMON /cnth2o/ V1h,V2h,DVh,NPTh,Ch(n_absrb),csh2o(n_absrb),cfh2o(n_absrb)
    !********************************************
    !
    COMMON /IFIL/ IRD,IPRcnt,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, &
         NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,   &    ! F00180
         NLTEFL,LNFIL4,LNGTH4                                 ! F00190
    
    common /cntscl/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
    
    icflg = -999
    !
    
    xself = 1.0d0
    XFRGN = 1.0d0
    XCO2C = 1.0d0
    XO3CN = 1.0d0
    XO2CN = 1.0d0
    XN2CN = 1.0d0
    XRAYL = 1.0d0
    
    
    kvert = nspec + 1
    mxone = 1
    
    ksmax2 = kztan(iscan(1,1))
    
    do iband = 1, nband
       do k = 1, ksmax2-1
          pave = p(k)*1013.15d0 ! convert in mbar
          tave = t(k)
          xlength = abs((z(k+1) - z(k)))*10000.0d0 !thickness needed in cm
          
          vmrh2o = xgas(1,k)
          wtot = xgas(1,k)*ccc(kvert,k)
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('H2O')) then
                vmrh2o = x(kk,k)
                wtot = x(kk,k)*ccc(kvert,k)
             end if
          end do
          
          w_dry = wtot * (1.-vmrh2o)
          wk(1) = vmrh2o * w_dry
          !ARGON
          WA     = 0.009     * W_dry 
          
          !NITROGEN
          wn2 = xgas(41,k) * ccc(kvert,k)
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('N2')) then
                wn2 = x(kk,k)*ccc(kvert,k)
             end if
          end do
          
          ! CO2
          wk(2) = xgas(2,k)*ccc(kvert,k)
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('CO2')) then
                wk(2) = x(kk,k)*ccc(kvert,k)
             end if
          end do
          
          ! Ozone
          wk(3) = xgas(3,k)*ccc(kvert,k)
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('O3')) then
                wk(3) = x(kk,k)*ccc(kvert,k)
             end if
          end do
          
          ! Oxygen
          wk(7) = xgas(7,k)*ccc(kvert,k)
          
          do kk = 1, nret
             if( trim(name(igas(kk))) .eq. trim('O2')) then
                wk(7) = x(kk,k)*ccc(kvert,k)
             end if
          end do
          
          
          wtot = 0.0d0
          w_dry = 0.0d0
          wa = 0.0d0
          wbroad=wn2+wa
          nmol_c = 7
          
          nmon = nm(iband)
          v1abs = wstart(iband)
          v2abs = wstart(iband) + nm(iband)*dn(iband)
          dvabs = dn(iband)
          nptabs = nm(iband)
          v1 = v1abs
          v2 = v2abs
          !          print *, v1abs, v2abs, dvabs, pave, tave, z(k+1), z(k), xlength, wbroad, wk(1:7)
                    
          absrb(1:nptabs)=0.0
          
          call contnm(1)
          
          open(10, file='h2ocont')
          DO I=1,NPTABS
             VI=V1ABS+dble(I-1)*DVABS
             WRITE (10, 910) VI, ABSRB(I)
          end DO
910       FORMAT(F10.3,1P,E12.3)
          close(10)

          mtckd(1, k, mxone:mxone+nm(iband)-1) = absrb(1:nm(iband))
          !          print *, k, mxone,mxone+nmon
       end do
    end do
    mxone = mxone + nmon
  end subroutine calc_h2o_continuum
  
end module h2o_continuum
