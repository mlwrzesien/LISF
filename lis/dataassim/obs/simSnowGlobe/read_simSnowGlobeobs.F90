!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_simSnowGlobeobs
!  \label{read_simSnowGlobeobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_simSnowGlobeobs(n, k, OBS_State, OBS_Pert_State) 
  ! !USES:

  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use LIS_pluginIndices
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the synthetic SWE observations produced from a LIS control run. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: sweField

  real,    pointer    :: obsl(:)
  real                :: sweobs(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real, allocatable       :: dummy(:)

  character(len=LIS_CONST_PATH_LEN) :: sweobsdir
  character(len=LIS_CONST_PATH_LEN) :: name

  logical             :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  logical             :: file_exists

  logical             :: readflag
  integer             :: sweid, status
  integer             :: fnd
  integer             :: c,r
  integer             :: ftn,t


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sweobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call simSnowGlobeSWE_filename(name,sweobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif

  if (readflag) then 
     allocate(dummy(LIS_rc%obs_ngrid(k)))
     write(LIS_logunit,*)  'Reading syn data ',trim(name)
     
     call ESMF_StateGet(OBS_State,"Observation01",sweField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(sweField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     

!     open(90, file=trim(name),form='unformatted')
!     do t=1,1
!        if(t==1) then 
!           call LIS_readvar_gridded(90,n,obsl)
!        else 
!           call LIS_readvar_gridded(90,n,dummy)
!        endif
!     end do
!     close(90)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     call LIS_verify(nf90_open(path=trim(name),mode=NF90_NOWRITE,ncid=ftn),&
          'Error opening file '//trim(name))
     call LIS_verify(nf90_inq_varid(ftn,'SWE_tavg',sweid),&
          'Error nf90_inq_varid: SWE_tavg')

     call LIS_verify(nf90_get_var(ftn,sweid,sweobs, &
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/), &
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'Error in nf90_get_var')
     
     call LIS_verify(nf90_close(ftn))

     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = &
                   sweobs(c,r)
           end if
        end do
     end do
     
#endif
     !-------------------------------------------------------------------------
     !  Apply LSM based QC and screening of observations
     !-------------------------------------------------------------------------
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+" &
          //trim(LIS_SnowGlobeSWEId)//char(0), n, k, OBS_state)

     call LIS_checkForValidObs(n, k, obsl, fnd, sweobs)     

     if (fnd .eq. 0) then
        data_upd_flag_local = .false.
     else
        data_upd_flag_local = .true.
     endif
     
#if (defined SPMD)
     call MPI_ALLGATHER(data_upd_flag_local, 1, &
          MPI_LOGICAL, data_upd_flag(:), &
          1, MPI_LOGICAL, LIS_mpi_comm, status)
     data_upd = any(data_upd_flag)
#else
     data_upd = data_upd_flag_local
#endif

     if (data_upd) then          
        
        do t=1,LIS_rc%obs_ngrid(k)
           gid(t) = t
           if(obsl(t).ne.-9999.0) then 
              assimflag(t) = 1 
           else
              assimflag(t) = 0 
           endif
        enddo
        
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .true., rc=status)
        call LIS_verify(status)
        
        if (LIS_rc%obs_ngrid(k) .gt. 0) then
           call ESMF_AttributeSet(sweField,"Grid Number",&
                gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(sweField,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
        endif
        
     else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .false., rc=status)
        call LIS_verify(status)
        return
     end if
  
     do t=1,LIS_rc%obs_ngrid(k)
        if(obsl(t).ne.-9999.0) then 
           if(obsl(t).lt.0.0) obsl(t) = 0.0
!        if(obsl(t).gt.200.0 ) obsl(t) = 200.0
        endif
     enddo
     deallocate(dummy)
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
  end if

end subroutine read_simSnowGlobeobs

subroutine simSnowGlobeSWE_filename(name, ndir, yr, mo,da,hr,mn)

  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  character(len=LIS_CONST_PATH_LEN) :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  

  name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/SimObs_'//&
       trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.nc'  

end subroutine simSnowGlobeSWE_filename



