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
! !ROUTINE: noahmp401_updateswevars
! \label{noahmp401_updateswevars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1 and SNODEP
! 15 May 2019: Yeosang Yoon; Modified for NoahMP 4.0.1 and LDTSI
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with SWE
! 09 Jan 2020: Yeosang Yoon; Updated QC
!
! !INTERFACE:
subroutine noahmp401_updateswevars(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use noahmp401_lsmMod
  use LIS_logMod,   only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Returns the swe related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
!
!EOP

  type(ESMF_Field)       :: sweField, sweIncrField

  integer                :: t, gid
  integer                :: status
  real, pointer          :: swe(:), sweincr(:)
  real                   :: swetmp, sndens
  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: swemean(LIS_rc%ngrid(n))
  integer                :: nswemean(LIS_rc%ngrid(n))
 
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)


  update_flag    = .true.
  perc_violation = 0.0
  swemean       = 0.0
  nswemean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     swetmp = swe(t) + sweincr(t)

     if((swetmp.lt.0 )) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) +1
     endif

  enddo

  do gid=1,LIS_rc%ngrid(n)
     perc_violation(gid) = perc_violation(gid) / real(LIS_rc%nensem(n))
  enddo

! For ensembles that are unphysical, compute the ensemble average after excluding them. This
! is done only if the majority of the ensemble members are good (>80%)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     if(.not.update_flag(gid)) then         ! false
        if(perc_violation(gid).lt.0.2) then
           if(swe(t)+sweincr(t).ge.0) then
              swemean(gid) = swemean(gid) + swe(t)+sweincr(t)
              nswemean(gid) = nswemean(gid) + 1
           else
             swemean(gid) = 0.0
           endif
        endif
     endif
  enddo

  do gid=1,LIS_rc%ngrid(n)
     if(nswemean(gid).gt.0) then
        swemean(gid) = swemean(gid) / real(nswemean(gid))
     endif
  enddo

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not update.

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)


     swetmp = swe(t) + sweincr(t)

!Use the model's snow density from the previous timestep
     sndens = 0.0
     if(noahmp401_struc(n)%noahmp401(t)%snowh.gt.0) then
        sndens = noahmp401_struc(n)%noahmp401(t)%sneqv/&
             noahmp401_struc(n)%noahmp401(t)%snowh
     endif

     if(update_flag(gid)) then
        swe(t) = swetmp
     elseif(perc_violation(gid).lt.0.2) then
       if(swetmp.lt.0.0) then  ! average of the good ensemble members
          swe(t) = swemean(gid)
       else
          swe(t) = swetmp
       endif
     else            ! do not update
       swe(t) = noahmp401_struc(n)%noahmp401(t)%sneqv
     end if

  enddo

end subroutine noahmp401_updateswevars

