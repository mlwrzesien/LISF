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
! !ROUTINE: noahmp401_map_snow
! \label{noahmp401_map_snow}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 12 Jun 2013, May 31 2014: Yuqiong Liu; Modified to combine different DI approaches,
!   Please refer to Rodell & House (2004, JHM), De Lannoy (2012, WRR) and Liu et al. (2013, AWR)
!
! !INTERFACE:
subroutine noahmp401_map_snow(n,k,OBS_State,LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_constantsMod, only  : LIS_CONST_TKFRZ
  use LIS_logMod,   only  : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_lsmMod
  use noahmp401_lsmMod

!  use LIS_topoMod, only : LIS_topo

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_Incr_State
! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding 
!  variables in the LSM state for SCA data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
  !EOP

  real, parameter          :: nominal_swe = 10.0
  type(ESMF_Field)         :: sweincrField
  type(ESMF_Field)         :: obs_sca_field
  real, pointer            :: sweincr(:)
  type(ESMF_Field)         :: snodincrField
  real, pointer            :: snodincr(:)
  real, pointer            :: scaobs(:)
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  character*100,allocatable    :: obs_state_objs(:)
  real                     :: tempc, dsnew, hnewc, snowhc
  real                     :: snowh, newsn, newsnc
  integer                  :: st_id, en_id
  real, allocatable            :: swe(:)
  real, allocatable            :: snod(:)

  allocate(swe(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(snod(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  call ESMF_StateGet(LSM_Incr_State,"SWE",sweincrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweincrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodincrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(snodincrField,localDE=0,farrayPtr=snodincr,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))

  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_sca_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_sca_field,localDE=0,farrayPtr=scaobs,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     swe(t)  = noahmp401_struc(n)%noahmp401(t)%sneqv
     snod(t) = noahmp401_struc(n)%noahmp401(t)%snowh
  enddo
! Based on SCA, we update the SWE based on the rule

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
! remove the snow if observation has no snow
     call LIS_lsm_DAmapTileSpaceToObsSpace(n,k,t,st_id,en_id)
     if(scaobs(st_id).ne.-9999.0) then
        if(scaobs(st_id).lt.0.01) then
           swe(t) = 0.0
           snod(t) = 0.0
        endif
     endif
  enddo
!     if(scaobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index).ne.-9999.0) then 
!        if(scaobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index).lt.0.01) then 
!           swe(t) = 0.0
!           snod(t) = 0.0
!        endif        
!     endif
 
 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     sweincr(t)  = swe(t)  - noahmp401_struc(n)%noahmp401(t)%sneqv
     snodincr(t) = snod(t) - noahmp401_struc(n)%noahmp401(t)%snowh
!     if(LIS_localPet.eq.1.and.t.eq.122) then
!        print*, LIS_localPet,t, scaobs(t),swe(t), noahmp401_struc(n)%noahmp40!(t)%sneqv,sweincr(t)
!     endif
  enddo

  deallocate(obs_state_objs)
  deallocate(swe)
  deallocate(snod)
  
end subroutine noahmp401_map_snow
