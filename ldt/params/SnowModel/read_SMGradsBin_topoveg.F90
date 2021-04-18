!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_SMGradsBin_topoveg
! \label{read_SMGradsBin_topoveg}
!
! !REVISION HISTORY:
!  16Jul2020: Kristi Arsenault; Added SnowModel topo-veg reader
!
! !INTERFACE:
subroutine read_SMGradsBin_topoveg(n, array1, array2)

! !USES:
  use LDT_coreMod,        only : LDT_rc
  use LDT_logMod,         only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use SnowModel_parmsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array1(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real, intent(inout) :: array2(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves SnowModel's topo-veg data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array1]
!   output field with the retrieved topo data
!  \item[array2]
!   output field with the retrieved vege data
!  \end{description}
!EOP

  integer :: ftn
  integer :: c, r
  logical :: file_exists
! ____________________________

  array1 = LDT_rc%udef
  array2 = LDT_rc%udef

  inquire(file=trim(SnowModel_struc(n)%topoveg_file), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "SnowModel topo-veg map, ",&
           trim(SnowModel_struc(n)%topoveg_file),", not found."
     call LDT_endrun
  endif
  select case ( SnowModel_struc(n)%topoveg_gridtransform )
    case( "none", "neighbor" ) 
      write(LDT_logunit,*) "[INFO] Reading Grads_binary topoveg file: ",&
            trim(SnowModel_struc(n)%topoveg_file)
  case default
     write(LDT_logunit,*) "[ERR] Since the Topo-veg field involves discrete data values,"
     write(LDT_logunit,*) "  only 'neighbor' is currently supported spatial"
     write(LDT_logunit,*) "  transform types.  Please check your entries for this parameter."
     call LDT_endrun
  end select

  ftn = LDT_getNextUnitNumber()
  open( ftn, file=SnowModel_struc(n)%topoveg_file, &
        form="unformatted", access='direct',status='old', &
        recl=4*LDT_rc%lnc(n)*LDT_rc%lnr(n) )

!     Original code from Snowmodel's preprocess.f
!       open (unit=37,file=topoveg_fname, &
!             form='unformatted',access='direct',recl=4*nx*ny)
!       read (37,rec=1) ((topo_land(i,j),i=1,nx),j=1,ny)

  ! Read in topographic map
  read(ftn,rec=1) ((array1(c,r,1),c=1,LDT_rc%lnc(n)),r=1,LDT_rc%lnr(n))
  ! Read in vegetation type map
  read(ftn,rec=2) ((array2(c,r,1),c=1,LDT_rc%lnc(n)),r=1,LDT_rc%lnr(n))
 
  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading Grads_binary topo-veg file"

end subroutine read_SMGradsBin_topoveg
