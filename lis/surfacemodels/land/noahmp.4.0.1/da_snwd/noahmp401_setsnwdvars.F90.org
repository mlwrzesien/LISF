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
! !ROUTINE: noahmp401_setsnwdvars
! \label{noahmp401_setsnwdvars}
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
! 14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1 and SNODEP
! 15 May 2019: Yeosang Yoon; Modified for NoahMP 4.0.1 and LDTSI
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with SNWD
!
! !INTERFACE:
subroutine noahmp401_setsnwdvars(n, LSM_State)
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
  type(ESMF_Field)       :: snodField
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real                   :: dsneqv,dsnowh,snoden
  integer                :: t
  integer                :: status
  
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(noahmp401_struc(n)%noahmp401(t)%sneqv.gt.0) then

        dsnowh = snod(t) - noahmp401_struc(n)%noahmp401(t)%snowh  !in m

        if(noahmp401_struc(n)%noahmp401(t)%snowh.ne.0) then 
           snoden = noahmp401_struc(n)%noahmp401(t)%sneqv/&
                noahmp401_struc(n)%noahmp401(t)%snowh           
           dsneqv = dsnowh*snoden           
        else
           dsnowh = 0
           dsneqv = 0
        endif
        ! update
        call noahmp401_snow_update(n, t, dsneqv, dsnowh)
     endif
  enddo
end subroutine noahmp401_setsnwdvars


