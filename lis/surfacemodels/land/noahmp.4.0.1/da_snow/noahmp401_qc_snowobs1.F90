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
! !ROUTINE: noahmp401_qc_snowobs
! \label{noahmp401_qc_snowobs}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!  30 Jan 2015: Yuqiong Liu; added additional QC
!  03 Oct 2018: Yeosang Yoon; Modified for NoahMP 3.6
!  14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1
! 15 May 2019: Yeosang Yoon; Modified for NoahMP 4.0.1 and LDTSI
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with SNOW
!
! !INTERFACE:
subroutine noahmp401_qc_snowobs1(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the snow observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  ground is fully or partially covered with snow. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_snow_field

  real, pointer            :: snowobs(:)
  integer                  :: t,iz
  integer                  :: c,r,i,j,c1,r1,c2,r2
  integer                  :: gid
  logical                  :: found
  integer                  :: status
  integer                  :: rad
  real                     :: stc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: vegt(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: fveg_obs(LIS_rc%obs_ngrid(k))
  real                     :: tv_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))
  real                     :: mswe_obs(LIS_rc%obs_ngrid(k))
  real                     :: veg2d(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                     :: obs2d(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                     :: mswe2d(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                     :: snowliqtot(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: snowliqtot_obs(LIS_rc%obs_ngrid(k))

  call ESMF_StateGet(OBS_State,"Observation01",obs_snow_field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in noahmp401_qc_snowobs")
  call ESMF_FieldGet(obs_snow_field,localDE=0,farrayPtr=snowobs,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet failed in noahmp401_qc_snowobs")
  
  snowliqtot = 0.0
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     !stc1(t) = noahmp401_struc(n)%noahmp401(t)%sstc(1) ! get snow/veg temp.
     stc1(t) = noahmp401_struc(n)%noahmp401(t)%tslb(1) ! get snow/veg temp.
     vegt(t) = LIS_surface(n,1)%tile(t)%vegt

     do iz=1, NOAHMP401_struc(n)%nsnow
        snowliqtot(t) = snowliqtot(t) + NOAHMP401_struc(n)%noahmp401(t)%snowliq(iz)
     enddo
  enddo

  call LIS_convertPatchSpaceToObsSpace(n,k,&       
       LIS_rc%lsm_index, noahmp401_struc(n)%noahmp401(:)%tv,tv_obs) !tv: vegetation temperature. unit: K 
  call LIS_convertPatchSpaceToObsSpace(n,k,LIS_rc%lsm_index, &    !fveg: green vegetation fraction. unit: - 
       noahmp401_struc(n)%noahmp401(:)%fveg,fveg_obs)
  

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index,stc1,stc1_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index,noahmp401_struc(n)%noahmp401(:)%sneqv,mswe_obs)
  
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index,vegt,vegt_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index,snowliqtot,snowliqtot_obs)


  veg2d = 0.0
  mswe2d = 0.0
  obs2d = -9999.0
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
!        print *, LIS_localPet,c,r
        gid  = LIS_obs_domain(n,k)%gindex(c,r)
!        print *, snowliqtot_obs(gid),c,r,gid
        if(gid.ne.-1) then
           obs2d(c,r) = snowobs(gid)
           if(vegt_obs(gid).le.4.or.fveg_obs(gid).gt.0.7) then 
              veg2d(c,r) = 1.0
           else
              veg2d(c,r) = 0.0
           endif
           mswe2d(c,r) = mswe_obs(gid)
        endif
     enddo
  enddo

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        rad = 1
        found = .false.
        do while(.not.found) 
           if(veg2d(c,r).eq.1) then
              c1 = max(1,c-rad)
              c2 = min(LIS_rc%obs_lnc(k),c+rad)
              r1 = max(1,r-rad)
              r2 = min(LIS_rc%obs_lnr(k),r+rad)
              
              do j=r1,r2
                 do i=c1,c2
                    if(obs2d(i,j).ne.-9999.0.and.veg2d(c,r).eq.0) then
                       if(mswe2d(i,j).gt.0) then 
                          obs2d(c,r) = obs2d(i,j)*mswe2d(c,r)/mswe2d(i,j)
                       else
                          obs2d(c,r) = 0.0
                       endif
                       found = .true.
                       exit
                    endif
                 enddo
                 if(found) exit              
              enddo
              if(found) exit                            
           endif
           if(.not.found) rad = rad+1
           if(rad.gt.3) then
              found = .true. 
              exit
           endif
        enddo
     enddo
  enddo
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        gid = LIS_obs_domain(n,k)%gindex(c,r)
        if(gid.ne.-1) then        
           snowobs(gid) = obs2d(c,r)
        endif 
     enddo
  enddo 
   
  do t=1,LIS_rc%obs_ngrid(k)
     if(snowobs(t).ne.LIS_rc%udef) then
!        if(fveg_obs(t).gt.0.7) then
!           snowobs(t) = LIS_rc%udef
!        elseif(vegt_obs(t).le.4) then !forest types
!           snowobs(t) = LIS_rc%udef
        if(vegt_obs(t).eq.LIS_rc%glacierclass) then !TML: Eliminate Glaciers
           snowobs(t) = LIS_rc%udef
!assume that snow will not form at 5 deg. celcius or higher ground temp. 
        elseif(tv_obs(t).ge.278.15) then
           snowobs(t) = LIS_rc%udef
        elseif(stc1_obs(t).ge.278.15) then
           snowobs(t) = LIS_rc%udef
        elseif(snowliqtot(t).ge.1.0) then
           snowobs(t) = LIS_rc%udef
        endif
     endif
  enddo

  
end subroutine noahmp401_qc_snowobs1

