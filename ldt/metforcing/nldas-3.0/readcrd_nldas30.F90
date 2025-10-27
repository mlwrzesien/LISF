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
! !ROUTINE: readcrd_nldas30
! \label{readcrd_nldas30}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
!
! !INTERFACE:    
subroutine readcrd_nldas30()
! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod
  use nldas30_forcingMod, only : nldas30_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NLDAS-3 forcing
!  from the LDT configuration file. 
!   MLW: still need to add config option for reading in elevation for downscaling
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LDT_config,                                &
                            "NLDAS-3 forcing directory:",rc=rc)
  do n = 1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,                          &
                                  nldas30_struc(n)%nldas30dir,rc=rc)
     call LDT_verify(rc,'NLDAS-3 forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS-3 forcing terrain height file:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,nldas30_struc(n)%nldas30hgt_file,&
          rc=rc)
     call LDT_verify(rc,&
          'NLDAS-3 forcing terrain height file: not defined')
  enddo


  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) '[INFO] Using NLDAS-3 forcing'
     write(LDT_logunit,*) '[INFO] NLDAS-3 forcing directory: ',&
           trim(nldas30_struc(n)%nldas30DIR)

     nldas30_struc(n)%nldas30time1 = 3000.0
     nldas30_struc(n)%nldas30time2 = 0.0

  enddo
end subroutine readcrd_nldas30
