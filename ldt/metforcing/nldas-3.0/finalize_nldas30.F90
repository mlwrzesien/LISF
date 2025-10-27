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
! !ROUTINE: finalize_nldas30
! \label{finalize_nldas30}
!
! !REVISION HISTORY: 
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
! 23 Oct 2025: M Wrzesien, based on MERRA2 code
! 
! !INTERFACE:

subroutine finalize_nldas30(findex)

! !USES:
  use LDT_coreMod,       only : LDT_rc
  use nldas30_forcingMod, only : nldas30_struc
!
! !DESCRIPTION:
!  Routine to cleanup NLDAS-3 forcing related memory allocations.   
! 
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n=1,LDT_rc%nnest
    select case( LDT_rc%met_gridtransform(findex) )

     case( "bilinear" )
       deallocate(nldas30_struc(n)%n111)
       deallocate(nldas30_struc(n)%n121)
       deallocate(nldas30_struc(n)%n211)
       deallocate(nldas30_struc(n)%n221)
       deallocate(nldas30_struc(n)%w111)
       deallocate(nldas30_struc(n)%w121)
       deallocate(nldas30_struc(n)%w211)
       deallocate(nldas30_struc(n)%w221)

     case( "budget-bilinear" )
       deallocate(nldas30_struc(n)%n111)
       deallocate(nldas30_struc(n)%n121)
       deallocate(nldas30_struc(n)%n211)
       deallocate(nldas30_struc(n)%n221)
       deallocate(nldas30_struc(n)%w111)
       deallocate(nldas30_struc(n)%w121)
       deallocate(nldas30_struc(n)%w211)
       deallocate(nldas30_struc(n)%w221)
       deallocate(nldas30_struc(n)%n112)
       deallocate(nldas30_struc(n)%n122)
       deallocate(nldas30_struc(n)%n212)
       deallocate(nldas30_struc(n)%n222)
       deallocate(nldas30_struc(n)%w112)
       deallocate(nldas30_struc(n)%w122)
       deallocate(nldas30_struc(n)%w212)
       deallocate(nldas30_struc(n)%w222)

     case( "neighbor" )
       deallocate(nldas30_struc(n)%n113)
    end select

 enddo
 deallocate(nldas30_struc)

end subroutine finalize_nldas30
