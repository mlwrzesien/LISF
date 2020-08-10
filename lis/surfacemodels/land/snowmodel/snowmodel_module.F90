!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module snowmodel_module
!BOP
!
! !MODULE: snowmodel_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the 
!  data structure containing the SnowModel 1-d variables.
!  The variables specified in the data structure include: 
!
!  \begin{description}
!   \item[slopetype]
!    slope type (integer index)
!   \item[nroot] 
!    Number of root layers, depends on the vegetation type
!   \end{description}
!
! !REVISION HISTORY:
!  14 Apr 2020: K. Arsenault; Added G. Liston's SnowModel 
! 
!EOP
  implicit none

  PRIVATE
  type, public :: snowmodeldec

     integer :: nroot
     integer :: vegt

     real, allocatable :: lyrthk(:) !used only if distributed depths are defined

!-------------------------------------------------------------------------
! SnowModel State Variables
!-------------------------------------------------------------------------
     ! States that can be 0'd out at the start of the next water year
     real :: swe_depth      ! swe_depth, swed: Snow water equivalent depth (m)
     real :: snow_depth     ! snod: snow depth (m)  
     real :: ro_snow_grid   ! gets initialized as ro_snow
     real :: canopy_int_old
     real :: swe_depth_old

     ! Multilayer states
     real :: tslsnowfall    ! Initialized with tsls_threshold
     real :: snod_layer     ! (i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
     real :: swed_layer     ! (i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
     real :: ro_layer       ! (i,j,k) = ro_snow, initialized
     real :: T_old_layer    ! (i,j,k) = 273.15, initialized
     real :: KK             ! (i,j) = 0; this is the total number of layers set to

!       1) ((snow_d(i,j),i=1,nx),j=1,ny)
!       2) ((snow_depth(i,j),i=1,nx),j=1,ny)
!       3) ((canopy_int(i,j),i=1,nx),j=1,ny)
!       4) ((soft_snow_d(i,j),i=1,nx),j=1,ny)
!       5) ((ro_snow_grid(i,j),i=1,nx),j=1,ny)
!       6) ((swe_depth(i,j),i=1,nx),j=1,ny)
!       7) ((ro_soft_snow_old(i,j),i=1,nx),j=1,ny)
!       8) ((snow_d_init(i,j),i=1,nx),j=1,ny)
!       9) ((swe_depth_old(i,j),i=1,nx),j=1,ny)
!       10) ((canopy_int_old(i,j),i=1,nx),j=1,ny)
!       11) ((topo(i,j),i=1,nx),j=1,ny)
!       12) ((sum_sprec(i,j),i=1,nx),j=1,ny)


!-------------------------------------------------------------------------
! SnowModel Forcing Variables
!-------------------------------------------------------------------------
     real :: tair            ! tair_grid: Air temperature (K)
     real :: rh              ! rh_grid: Relative humidity (%)
     real :: windspd         ! windspd_grid, wspd: Wind speed (m/2)
     real :: winddir         ! winddir_grid, wdir: Wind direction (deg, true N)
     real :: uwind           ! uwind_grid
     real :: vwind           ! vwind_grid
     real :: totprecip       ! prec_grid: Water-equivalent precip (m/dt)
     real :: rainf           ! rpre or rain: Liquid precip, rain (m/dt)
     real :: snowf           ! spre or sprec: Solid precip, snowfall (m/dt)
     real :: swdown          ! Qsi_grid, qsin: incoming solar radiation (W/m2)
     real :: lwdown          ! Qli_grid, qlin: incoming lw radiation (W/m2)
     real :: cloud_fraction  ! cldf: cloud fraction (grid) (0-1)
     real :: qair            ! Specific humidity
     real :: psurf
     real :: fheight

!-------------------------------------------------------------------------
! SnowModel Parameters
!-------------------------------------------------------------------------
     real :: smtopo          ! SnowModel topography field
     real :: smvege          ! SnowModel vegetation class map

!-----------------------------------------------------------------------
! SnowModel Output variables
!-----------------------------------------------------------------------
     ! EBal:
     real :: tsfc        ! tsfc: Surface (skin) temperature (deg C)
     real :: qle         ! qlem: Emitted longwave radiation (W/m2)
     real :: qh          ! Sensible heat flux (W/m2)
     real :: qe          ! Latent heat flux (W/m2)
     real :: qc          ! Conductive heat flux (W/m2)
     real :: qm          ! Melt energy flux (W/m2)
     real :: albedo      ! Albedo (0-1), albd
     real :: e_balance   ! Energy balance error (W/m2)

     real :: swnet       ! Net shortwave (W/m2)
     real :: lwnet       ! Net longwave (W/m2)
   
     ! Snowpack:
     real :: sden        ! sden, xro_snow: snow density (kg/m3)
     real :: runoff      ! runoff, roff: runoff from snowpack base (m/dt)
     real :: glmt        ! glacier_melt, glmt: SWE melt from glacier ice (m)
     real :: qcs         ! qcs: Canopy sublimation (m/dt)
     real :: canopint    ! canopy_int:  Canopy interception store (m)
     real :: snowcover   

     real :: sublim      ! swesublim, ssub: static-surface sublimation (m)
     real :: swemelt     ! swemelt, smlt: SWE melt (m)

     real :: sumsublim   ! Summed static-surface sublimation (m)
     real :: sumqcs      ! Summed canopy sublim during year (m)
     real :: sumprec     ! Summed precipitation during year (m)
     real :: sumsprec    ! sspr, Summed snow precip during year (m)
     real :: sumroff     ! Summed runoff during the year (m)
     real :: sumswemelt  ! ssmt, Summed snow-water-equivalent melt (m)
     real :: sumunload   ! Sum_unload: Summed canopy unloading during the year (m)
     real :: w_balance   ! wbal: Summed water balance error durign year (m)

! SnowTran:
     real :: snow_d
     real :: wbal_qsubl
     real :: wbal_salt
     real :: wbal_susp
     real :: wbal_subgrid
     real :: sum_qsubl
     real :: sum_trans

! Multi-layer snow fields:
     real :: diam_layer  ! (i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
     real :: flux_layer  ! (i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
     real :: Told_layer  ! (i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
     real :: gamma_layer ! (i,j,k),i=1,nx),j=1,ny),k=1,nz_max)

  end type snowmodeldec

  real, public, allocatable :: glb_topoland(:,:)
  real, public, allocatable :: glb_vegtype(:,:)

end module snowmodel_module
