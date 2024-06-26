!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_transform_usafsi
! \label{noahmp401_transform_usafsi}
!
! !REVISION HISTORY:
! 25Jun2006: Sujay Kumar: Initial Specification
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!  03 Oct 2018; Yeosang Yoon; Modified for NoahMP 3.6
!  14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1
! 15 May 2019: Yeosang Yoon; Modified for NoahMP 4.0.1 and LDTSI
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with USAFSI
!
! !INTERFACE:
subroutine noahmp401_transform_usafsi(n,OBS_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use noahmp401_lsmMod
!EOP
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the USAFSI state
!  (mm) to the lsm state
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP

  ! Since USAFSI is already in meters, no work is needed here.

end subroutine noahmp401_transform_usafsi
