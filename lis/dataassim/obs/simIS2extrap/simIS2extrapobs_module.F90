!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: simIS2extrapobs_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle synthetic snow water equivalent (SNWD) observations. 
!   
! !REVISION HISTORY: 
!  27Feb05    Sujay Kumar;   Initial Specification
! 
! 
module simIS2extrapobs_module
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: simIS2extrapobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: IS2extrap_struc
!EOP
  type, public:: IS2extrap_dec

     logical                :: startMode
     integer                :: nc
     integer                :: nr
     integer                :: mi
     real                   :: datares
     real                   :: ssdev_inp
!     type(proj_info)        :: proj
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real,  allocatable     :: w11(:)
     real,  allocatable     :: w12(:)
     real,  allocatable     :: w21(:)
     real,  allocatable     :: w22(:)
     real, allocatable      :: snwd(:)
!     real  :: gridDesci(50)
  end type IS2extrap_dec

  type(IS2extrap_dec),allocatable :: IS2extrap_struc(:)


contains
!BOP
! 
! !ROUTINE: simIS2extrapobs_setup
! \label{simIS2extrapobs_setup}
! 
! !INTERFACE: 
  subroutine simIS2extrapobs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION:
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for SNWD assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n 
    integer                ::  ftn
    integer                ::  i
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  synsnwdobsdir
    character*100          ::  temp
    real, allocatable          :: ssdev(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real,         allocatable  ::  varmin(:)
    real,         allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    real                   :: gridDesci(20)

!    print *,'before allocate'
    allocate(IS2extrap_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"simulated IS2extrap Snow Depth data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synsnwdobsdir,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            synsnwdobsdir, rc=status)
       call LIS_verify(status)
    enddo

    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(OBS_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Update Time",&
            -99.0, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
            LIS_rc%obs_ngrid(k),rc=status)
       call LIS_verify(status)

    enddo

    write(LIS_logunit,*)'read simulated IS2extrap Snow Depth data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%obs_ngrid(k). 
!   
!----------------------------------------------------------------------------
    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_obsvecGrid(n,k),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)
       
!Perturbations State
       write(LIS_logunit,*) 'Opening attributes for observations ',&
            trim(LIS_rc%obsattribfile(k))
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
       read(ftn,*)
       read(ftn,*) LIS_rc%nobtypes(k)
       read(ftn,*)
    
       allocate(vname(LIS_rc%nobtypes(k)))
       allocate(varmax(LIS_rc%nobtypes(k)))
       allocate(varmin(LIS_rc%nobtypes(k)))
       
       do i=1,LIS_rc%nobtypes(k)
          read(ftn,fmt='(a40)') vname(i)
          read(ftn,*) varmin(i),varmax(i)
          write(LIS_logunit,*) vname(i),varmin(i),varmax(i)
       enddo
       call LIS_releaseUnitNumber(ftn)  
       
 !      call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
 !      call LIS_verify(status)

       allocate(ssdev(LIS_rc%obs_ngrid(k)))

       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 

          allocate(obs_pert%vname(1))
          allocate(obs_pert%perttype(1))
          allocate(obs_pert%ssdev(1))
          allocate(obs_pert%stdmax(1))
          allocate(obs_pert%zeromean(1))
          allocate(obs_pert%tcorr(1))
          allocate(obs_pert%xcorr(1))
          allocate(obs_pert%ycorr(1))
          allocate(obs_pert%ccorr(1,1))

          call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
               obs_pert)

          ssdev = obs_pert%ssdev(1)

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 

          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1), rc=status)
          call LIS_verify(status)
          
          if(LIS_rc%obs_ngrid(k).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
             call LIS_verify(status)
          endif

          call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
               obs_pert%stdmax(1), rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
               obs_pert%zeromean(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
               obs_pert%tcorr(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
               obs_pert%xcorr(1),rc=status)
          
          call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
               obs_pert%ycorr(1),rc=status)

          call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
               obs_pert%ccorr(1,:),itemCount=1,rc=status)

          call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
          call LIS_verify(status)
       endif
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)
    enddo

    write(LIS_logunit,*) 'Created the States to hold the observations data'


    do n=1,LIS_rc%nnest
       IS2extrap_struc(n)%nc = 60
       IS2extrap_struc(n)%nr = 59

!       call map_set(PROJ_LATLON, 24.94958,-124.73375,&
!            0.0, 0.00833,0.00833, 0.0,&
!            IS2extrap_struc(n)%nc,&
!            IS2extrap_struc(n)%nr,&
!            IS2extrap_struc(n)%proj)

!       gridDesci = 0
!       IS2extrap_struc(n)%gridDesci(1) = 3
!       IS2extrap_struc(n)%gridDesci(2) = IS2extrap_struc(n)%nc
!       IS2extrap_struc(n)%gridDesci(3) = IS2extrap_struc(n)%nr
!       IS2extrap_struc(n)%gridDesci(4) = 35.100
!       IS2extrap_struc(n)%gridDesci(5) = -112.497
!       IS2extrap_struc(n)%gridDesci(6) = 8
!       IS2extrap_struc(n)%gridDesci(7) = 46.0
!       IS2extrap_struc(n)%gridDesci(8) = 0.50
!       IS2extrap_struc(n)%gridDesci(9) = 0.50
!       IS2extrap_struc(n)%gridDesci(10) = 39.0
!       IS2extrap_struc(n)%gridDesci(11) = -100.0
!       IS2extrap_struc(n)%gridDesci(20) = 0

       gridDesci = 0
       gridDesci(1) = 3
       gridDesci(2) = IS2extrap_struc(n)%nc
       gridDesci(3) = IS2extrap_struc(n)%nr
       gridDesci(4) = 37.71577
       gridDesci(5) = -119.84483
       gridDesci(6) = 8
       gridDesci(7) = 38.00
       gridDesci(8) = 1.0
       gridDesci(9) = 1.0
       gridDesci(10) = 37.00
       gridDesci(11) = -120.0
       gridDesci(20) = 0

!       IS2extrap_struc(n)%datares = 0.00833

       allocate(IS2extrap_struc(n)%n11(&
            IS2extrap_struc(n)%nc*IS2extrap_struc(n)%nr))

       IS2extrap_struc(n)%mi = IS2extrap_struc(n)%nc*IS2extrap_struc(n)%nr
!       allocate(IS2extrap_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IS2extrap_struc(n)%n12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IS2extrap_struc(n)%n21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IS2extrap_struc(n)%n22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       allocate(IS2extrap_struc(n)%w11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IS2extrap_struc(n)%w12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IS2extrap_struc(n)%w21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(IS2extrap_struc(n)%w22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       allocate(IS2extrap_struc(n)%snwd(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       IS2extrap_struc(n)%snwd = LIS_rc%udef

!       call bilinear_interp_input(n, gridDesci(:),& 
!            IS2extrap_struc(n)%n11,IS2extrap_struc(n)%n12,&
!            IS2extrap_struc(n)%n21,IS2extrap_struc(n)%n22,&
!            IS2extrap_struc(n)%w11,IS2extrap_struc(n)%w12,&
!            IS2extrap_struc(n)%w21,IS2extrap_struc(n)%w22)


!       print *,'before upscale in simModule'

       call upscaleByAveraging_input(&
            gridDesci(:),&
            LIS_rc%obs_gridDesc(k,:),&
            IS2extrap_struc(n)%nc*IS2extrap_struc(n)%nr,&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
            IS2extrap_struc(n)%n11)

!       print *,'after upscaline in simModule'

!       call LIS_registerAlarm("IS2extrap read alarm",&
!            86400.0, 86400.0)

       IS2extrap_struc(n)%startMode = .true.

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

!       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
!       call LIS_verify(status)
      
!       print *,'test end of simModule'
    enddo

    
  end subroutine simIS2extrapobs_setup
  
end module simIS2extrapobs_module
