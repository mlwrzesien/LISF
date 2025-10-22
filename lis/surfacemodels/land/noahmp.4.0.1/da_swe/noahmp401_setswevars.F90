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
! !ROUTINE: noahmp401_setswevars
! \label{noahmp401_setswevars}
!
! !REVISION HISTORY:
! 15 Aug 2017: Sujay Kumar; Initial Specification
! 03 Oct 2018: Yeosang Yoon; Modified for NoahMP 3.6
!
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 21 Jul 2011: James Geiger; Modified for Noah 3.2
! 03 Oct 2018: Yeosang Yoon; Modified for NoahMP 3.6
! 14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1 and SWEEP
! 15 May 2019: Yeosang Yoon; Modified for NoahMP 4.0.1 and LDTSI
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with SWE
!
! !INTERFACE:
subroutine noahmp401_setswevars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  use LIS_logMod, only : LIS_logunit, LIS_verify, LIS_endrun
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the snow progognostic variables to noah's
!  model space. The state vector consists of total SWE and snow depth. 
!  This routine also updates other model prognostics (snice, snliq,
!  snow thickness, snow temperature) based on the update. 
! 
!EOP
  type(ESMF_Field)       :: sweField
  real, pointer          :: swe(:)
  real                   :: dsneqv,dsnowh,snoden
  integer                :: t
  integer                :: status
  
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(noahmp401_struc(n)%noahmp401(t)%sneqv.gt.0) then

        dsneqv = swe(t) - noahmp401_struc(n)%noahmp401(t)%sneqv  !in m

        if(noahmp401_struc(n)%noahmp401(t)%tslb(1).ge.278.15.and.& !MLWdebugging
             dsneqv.gt.0) then
           dsneqv = 0
           dsnowh = 0
        endif

        if(noahmp401_struc(n)%noahmp401(t)%snowh.ne.0) then 
           snoden = noahmp401_struc(n)%noahmp401(t)%sneqv/&
                noahmp401_struc(n)%noahmp401(t)%snowh           
           dsnowh = dsneqv/snoden
        else
           dsnowh = 0
           dsneqv = 0
        endif

        if(abs(dsneqv).ge.0.010) then           
           ! update
           call noahmp401_snow_update(n, t, dsneqv, dsnowh)
        endif
     endif
     
  enddo
end subroutine noahmp401_setswevars


