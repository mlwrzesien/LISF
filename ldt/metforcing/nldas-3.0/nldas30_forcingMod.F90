!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nldas30_forcingMod
!BOP
! !MODULE: nldas30_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the NLDAS-3 forcing data.
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt nldas30\_struc}
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the ECMWF data
!  \item[nldas30time1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[nldas30time2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[nldas30dir]
!    Directory containing the input data
!  \item[nldas30hgt_file]
!    File with the terrain height definition for the input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for bilinear interpolation.
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for n. neighbor interpolation.
!  \item[findtime1, findtime2]
!    boolean flags to indicate which time is to be read for
!    temporal interpolation.
!  \end{description}
!
! !USES:
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_nldas30      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nldas30_struc

!EOP
  type, public ::  nldas30_type_dec
     real         :: ts
     integer      :: nc, nr
     character(len=LDT_CONST_PATH_LEN) :: nldas30dir   !MERRA2 Forcing Directory
     real*8       :: nldas30time1, nldas30time2, ringtime
     logical      :: reset_flag

     integer                :: mi
     real                   :: gridDesc(20)
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)

     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)
     integer, allocatable   :: n113(:)
     integer                :: findtime1, findtime2
     logical                :: startFlag, dayFlag
     real, allocatable      :: nldasforc1(:,:,:), nldasforc2(:,:,:)

     integer                :: uselml
     integer                :: usecorr

     character*140          :: nldas30hgt_file

  end type nldas30_type_dec

  type(nldas30_type_dec), allocatable :: nldas30_struc(:)

contains

!BOP
!
! !ROUTINE: init_nldas30
! \label{init_nldas30}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
! 23 Oct 2025: M. Wrzesien, code for NLDAS-3 based on MERRA2
!
! !INTERFACE:
  subroutine init_nldas30(findex)

! !USES:
    use LDT_coreMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for MERRA2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_nldas30](\ref{readcrd_nldas30}) \newline
!     reads the runtime options specified for MERRA2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    integer :: n
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt

    allocate(nldas30_struc(LDT_rc%nnest))
    write(LDT_logunit,fmt=*)"[INFO] Initializing NLDAS-3 forcing grid ... "

 !- Read in config file entries:
    call readcrd_nldas30()

    do n=1, LDT_rc%nnest
       nldas30_struc(n)%ts = 3600  
       call LDT_update_timestep(LDT_rc, n, nldas30_struc(n)%ts)
    enddo

    nldas30_struc%reset_flag = .false.

    LDT_rc%met_nf(findex) = 8
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true.

    nldas30_struc%nc = 11700
    nldas30_struc%nr = 6500
    LDT_rc%met_nc(findex) = nldas30_struc(1)%nc
    LDT_rc%met_nr(findex) = nldas30_struc(1)%nr

    do n=1,LDT_rc%nnest
       nldas30_struc(n)%nc = 11700
       nldas30_struc(n)%nr = 6500

       nldas30_struc(n)%gridDesc(:) = 0.

       nldas30_struc(n)%gridDesc(1) = 0
       nldas30_struc(n)%gridDesc(2) = nldas30_struc(n)%nc
       nldas30_struc(n)%gridDesc(3) = nldas30_struc(n)%nr
       nldas30_struc(n)%gridDesc(4) = 7.005
       nldas30_struc(n)%gridDesc(5) = -168.995
       nldas30_struc(n)%gridDesc(6) = 128
       nldas30_struc(n)%gridDesc(7) = 71.995
       nldas30_struc(n)%gridDesc(8) = -52.005
       nldas30_struc(n)%gridDesc(9) = 0.01
       nldas30_struc(n)%gridDesc(10) = 0.01
       nldas30_struc(n)%gridDesc(20) = 0

       LDT_rc%met_gridDesc(findex,1:20) = nldas30_struc(n)%gridDesc(1:20)

       nldas30_struc(n)%mi = nldas30_struc(n)%nc*nldas30_struc(n)%nr

       ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )

        case( "bilinear" )
          allocate(nldas30_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, nldas30_struc(n)%gridDesc(:),&
               nldas30_struc(n)%n111,nldas30_struc(n)%n121,&
               nldas30_struc(n)%n211,nldas30_struc(n)%n221,&
               nldas30_struc(n)%w111,nldas30_struc(n)%w121,&
               nldas30_struc(n)%w211,nldas30_struc(n)%w221)

        case( "budget-bilinear" )
          allocate(nldas30_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas30_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, nldas30_struc(n)%gridDesc(:),&
               nldas30_struc(n)%n111,nldas30_struc(n)%n121,&
               nldas30_struc(n)%n211,nldas30_struc(n)%n221,&
               nldas30_struc(n)%w111,nldas30_struc(n)%w121,&
               nldas30_struc(n)%w211,nldas30_struc(n)%w221)

          allocate(nldas30_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas30_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas30_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas30_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas30_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas30_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas30_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas30_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n, nldas30_struc(n)%gridDesc(:),&
               nldas30_struc(n)%n112,nldas30_struc(n)%n122,&
               nldas30_struc(n)%n212,nldas30_struc(n)%n222,&
               nldas30_struc(n)%w112,nldas30_struc(n)%w122,&
               nldas30_struc(n)%w212,nldas30_struc(n)%w222)

        case( "neighbor" )
          allocate(nldas30_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n, nldas30_struc(n)%gridDesc(:),&
               nldas30_struc(n)%n113)

        case( "none" )
          write(LDT_logunit,*) "[INFO] No interpolation applied for NLDAS-3 ..."
          write(LDT_logunit,*) " -- Reading native grid to be written out to -- " 

        case default
          write(LDT_logunit,*) '[ERR] Interpolation option '// &
               trim(LDT_rc%met_gridtransform(findex))//&
               ' for NLDAS-3 forcing is not supported'
          call LDT_endrun()
      end select

      call LDT_registerAlarm("NLDAS3 forcing alarm",&
           86400.0,86400.0)
      nldas30_struc(n)%startFlag = .true.
      nldas30_struc(n)%dayFlag = .true.

      if (trim(LDT_rc%runmode) == "Metforce processing" .or. &
          trim(LDT_rc%runmode) == "Metforce temporal downscaling" .or. &
          trim(LDT_rc%runmode) == "Statistical downscaling of met forcing") then
         allocate(nldas30_struc(n)%nldasforc1(&
           LDT_rc%met_nf(findex), 24, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n)))
         allocate(nldas30_struc(n)%nldasforc2(&
           LDT_rc%met_nf(findex), 24, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n)))

         nldas30_struc(n)%nldasforc1 = LDT_rc%udef
         nldas30_struc(n)%nldasforc2 = LDT_rc%udef
      endif
    enddo

    if (trim(LDT_rc%runmode) == "Metforce processing" .or. &
        trim(LDT_rc%runmode) == "Metforce temporal downscaling" .or. &
        trim(LDT_rc%runmode) == "Statistical downscaling of met forcing") then
       write(LDT_logunit,*)"[INFO] NLDAS-3 time interp option :: ",&
             trim(LDT_rc%met_tinterp(findex))
    endif

  end subroutine init_nldas30
end module nldas30_forcingMod

