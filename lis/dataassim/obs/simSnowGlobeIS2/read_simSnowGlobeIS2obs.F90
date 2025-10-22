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
! !ROUTINE: read_simSnowGlobeIS2obs
!  \label{read_simSnowGlobeIS2obs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_simSnowGlobeIS2obs(n, k, OBS_State, OBS_Pert_State) 
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
  use simSnowGlobeIS2obs_module, only : SnowGlobeIS2_struc
  
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
!  reads the synthetic SNWD observations produced from a LIS control run. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: snwdField

  real,    pointer    :: obsl(:)
  real                :: snwdobs_out(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                :: snwdobs(SnowGlobeIS2_struc(n)%nc,SnowGlobeIS2_struc(n)%nr)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real, allocatable       :: dummy(:)

  character(len=LIS_CONST_PATH_LEN) :: snwdobsdir
  character(len=LIS_CONST_PATH_LEN) :: name

  logical             :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd

  logical*1, allocatable :: snwd_data_b(:)
  logical*1              :: snwdobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  logical             :: file_exists
  integer*2, allocatable :: var(:,:)
  real, allocatable      :: snwd1d(:)

  logical             :: readflag
  integer             :: snwdid, status
  integer             :: fnd
  integer             :: c,r
  integer             :: ftn,t


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       snwdobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call simSnowGlobeIS2SNWD_filename(name,snwdobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)
!  print *,'snwdobsdir=',snwdobsdir
!  print *,'name=',name

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true.
!     print *,'reading SimObs file' 
  else 
     readflag = .false.
  endif

  if (readflag) then 
     allocate(dummy(LIS_rc%obs_ngrid(k)))
     allocate(snwd_data_b(SnowGlobeIS2_struc(n)%nc*SnowGlobeIS2_struc(n)%nr))
     allocate(var(SnowGlobeIS2_struc(n)%nc,SnowGlobeIS2_struc(n)%nr))
     allocate(snwd1d(SnowGlobeIS2_struc(n)%nc*SnowGlobeIS2_struc(n)%nr))

     call ESMF_StateGet(OBS_State,"Observation01",snwdField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(snwdField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
  
     write(LIS_logunit,*)  'Reading syn data ',trim(name)
     
     call ESMF_StateGet(OBS_State,"Observation01",snwdField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(snwdField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     call LIS_verify(nf90_open(path=trim(name),mode=NF90_NOWRITE,ncid=ftn),&
          'Error opening file '//trim(name))
     call LIS_verify(nf90_inq_varid(ftn,'SnowDepth_tavg',snwdid),&
          'Error nf90_inq_varid: SnowDepth_tavg')

     call LIS_verify(nf90_get_var(ftn,snwdid,snwdobs, &
          count=(/SnowGlobeIS2_struc(n)%nc,SnowGlobeIS2_struc(n)%nr,1/),&
          start=(/1,1,1/)),&
          'Error in nf90_get_var')
 !    print *,'within block'
     call LIS_verify(nf90_close(ftn))


!--------------------------------------------------------------------------
! Interpolate to the observation grid
!-------------------------------------------------------------------------- 
        snwd1d = LIS_rc%udef
        snwd_data_b = .false.

        do r=1, SnowGlobeIS2_struc(n)%nr
           do c=1, SnowGlobeIS2_struc(n)%nc
              if(snwdobs(c,SnowGlobeIS2_struc(n)%nr-r+1).ge.0) then
                 snwd1d(c+(r-1)*SnowGlobeIS2_struc(n)%nc) = &
                      snwdobs(c,SnowGlobeIS2_struc(n)%nr-r+1)
                 snwd_data_b(c+(r-1)*SnowGlobeIS2_struc(n)%nc)=.true.
              endif
           enddo
        enddo
        deallocate(var)

        call upscaleByAveraging(&
             SnowGlobeIS2_struc(n)%nc*SnowGlobeIS2_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
             LIS_rc%udef,&
             SnowGlobeIS2_struc(n)%n11,&
             snwd_data_b,snwd1d,&
             snwdobs_b_ip,snwdobs_out)

        deallocate(snwd1d)
        deallocate(snwd_data_b)

!end MLW edits





     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = &
                   snwdobs_out(c+(r-1)*LIS_rc%obs_lnc(k))
           end if
        end do
     end do
     
#endif
     !-------------------------------------------------------------------------
     !  Apply LSM based QC and screening of observations
     !-------------------------------------------------------------------------
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+" &
          //trim(LIS_SnowGlobeIS2SNWDId)//char(0), n, k, OBS_state)

     call LIS_checkForValidObs(n, k, obsl, fnd, snwdobs)     

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
           call ESMF_AttributeSet(snwdField,"Grid Number",&
                gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(snwdField,"Assimilation Flag",&
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

end subroutine read_simSnowGlobeIS2obs

subroutine simSnowGlobeIS2SNWD_filename(name, ndir, yr, mo,da,hr,mn)

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

end subroutine simSnowGlobeIS2SNWD_filename



