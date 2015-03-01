module ocean_model_mod

use   mpp_domains_mod, only: domain2d, mpp_define_layout, mpp_define_domains, mpp_get_compute_domain

use           fms_mod, only: error_mesg, FATAL, write_version_number, mpp_npes, read_data
use           fms_mod, only: field_size, field_exist, get_mosaic_tile_grid

use  time_manager_mod, only: time_type

use coupler_types_mod, only: coupler_2d_bc_type

use        mosaic_mod, only: get_mosaic_ntiles, get_mosaic_grid_sizes, get_mosaic_xgrid
use        mosaic_mod, only: get_mosaic_xgrid_size, calc_mosaic_grid_area

use     constants_mod, only: PI, RADIUS

implicit none
private

public :: ocean_model_init, ocean_model_end, update_ocean_model, &
          ice_ocean_boundary_type, ocean_grids_type, &
          ocean_model_flux_init, ocean_model_init_sfc, &
          ocean_stock_pe, ocean_model_restart, &
          ice_ocn_bnd_type_chksum, ocean_public_type_chksum

public    ocean_model_data_get
interface ocean_model_data_get
   module procedure ocean_model_data1D_get 
   module procedure ocean_model_data2D_get 
end interface

!-----------------------------------------------------------------------

type ice_ocean_boundary_type
  real, dimension(:,:), pointer :: u_flux =>NULL(), &
                                   v_flux =>NULL(), &
                                   t_flux =>NULL(), &
                                   q_flux =>NULL(), &
                                   salt_flux =>NULL(), &
                                   lw_flux =>NULL(), &
                                   sw_flux =>NULL(), &
                                   sw_flux_nir_dir =>NULL(), &
                                   sw_flux_nir_dif =>NULL(), &
                                   sw_flux_vis_dir =>NULL(), &
                                   sw_flux_vis_dif =>NULL(), &
                                   lprec =>NULL(), &
                                   fprec  =>NULL()
  real, dimension(:,:), pointer :: runoff =>NULL(), &
                                   calving  =>NULL()
  real, pointer, dimension(:,:) :: runoff_hflx     =>NULL() ! heat flux of liquid runoff (kg/m2/s) 
  real, pointer, dimension(:,:) :: calving_hflx    =>NULL() ! heat flux of frozen runoff (kg/m2/s) 
  real, dimension(:,:), pointer :: p  =>NULL()
  real, dimension(:,:,:), pointer :: data  =>NULL()
  integer :: xtype
  type(coupler_2d_bc_type)      :: fluxes
end type ice_ocean_boundary_type

!-----------------------------------------------------------------------

 type ocean_grids_type
    real,    pointer, dimension(:)   :: lon_bnd =>NULL(), lat_bnd =>NULL()
    real,    pointer, dimension(:,:) :: lon =>NULL(), lat =>NULL()
    logical, pointer, dimension(:,:) :: mask  =>NULL()
 end type

!-----------------------------------------------------------------------

type, public :: ocean_public_type
   type (domain2d)               :: Domain
   type (ocean_grids_type)       :: Global, Data
   real, pointer, dimension(:,:) :: t_surf =>NULL() , &
                                    frazil =>NULL() , &
                                    u_surf =>NULL() , &
                                    v_surf =>NULL() , &
                                    s_surf =>NULL() , &
                                    area   =>NULL() , &
                                    sea_lev=>NULL()
   logical, pointer, dimension(:,:) :: maskmap =>NULL()
   logical :: is_ocean_pe = .false.
   integer, pointer :: pelist(:) =>NULL()
   integer, dimension(3)            :: axes    
   type(coupler_2d_bc_type)         :: fields
end type ocean_public_type

  type, public ::  ocean_state_type; private
     ! This type is private, and can therefore vary between different ocean models.
     ! All information entire ocean state may be contained here, although it is not
     ! necessary that this is implemented with all models.
     logical       :: is_ocean_pe = .false.       ! .true. on processors that run the ocean model.
  end type ocean_state_type

!-----------------------------------------------------------------------

   character(len=128) :: version = '$Id: ocean_model.F90,v 18.0.2.1 2010/09/07 18:44:36 pjp Exp $'
   character(len=128) :: tagname = '$Name: testing $'

contains

!#######################################################################

 subroutine update_ocean_model (Ice_boundary, Ocean_state, Ocean_sfc, &
       time_start_update, Ocean_coupling_time_step)

 type(ice_ocean_boundary_type), intent(in)    :: Ice_boundary
 type(ocean_state_type),        pointer       :: Ocean_state
 type(ocean_public_type),       intent(inout) :: Ocean_sfc
 type(time_type), intent(in)                  :: time_start_update
 type(time_type), intent(in)                  :: Ocean_coupling_time_step

 call error_mesg('ocean_model_mod', 'null ocean model should not be executed', FATAL )

 end subroutine update_ocean_model

!#######################################################################

 subroutine ocean_model_init (Ocean, Ocean_state, Time_init, Time)

 type(ocean_public_type), intent(inout) :: Ocean
 type(ocean_state_type),  pointer       :: Ocean_state
 type(time_type), intent(in) :: Time_init, Time

 real,    allocatable, dimension(:)     :: xgrid_area
 real,    allocatable, dimension(:,:)   :: geo_lonv, geo_latv, rmask, geo_lont, geo_latt, tmpx, tmpy, garea
 real,    allocatable, dimension(:,:,:) :: x_vert_t, y_vert_t
 integer, allocatable, dimension(:)     :: i1, j1, i2, j2
 integer                                :: i, j, ntiles, nfile_axo, nxgrid, n, m, grid_version, nlon, nlat
 integer                                :: isd, ied, jsd, jed
 integer                                :: layout(2), siz(4), nx(1), ny(1)
 character(len=256)                     :: grid_file = "INPUT/grid_spec.nc"
 character(len=256)                     :: ocean_mosaic, tile_file, ocn_mosaic_file, axo_file

 call write_version_number(version, tagname)


    if(field_exist(grid_file, 'geolon_t')) then
       grid_version = 0
       call field_size( grid_file, 'geolon_t', siz)
       nlon = siz(1)
       nlat = siz(2)
    else if(field_exist(grid_file, 'x_T')) then
       grid_version = 1
       call field_size( grid_file, 'x_T', siz)
       nlon = siz(1)
       nlat = siz(2)
    else if(field_exist(grid_file, 'ocn_mosaic_file') ) then ! read from mosaic file
       grid_version = 2
       call read_data(grid_file, "ocn_mosaic_file", ocean_mosaic)
       ocean_mosaic = "INPUT/"//trim(ocean_mosaic)
       ntiles = get_mosaic_ntiles(ocean_mosaic)
       if(ntiles .NE. 1) call error_mesg('ocean_model_init', ' ntiles should be 1 for ocean mosaic, contact developer', FATAL)
       call get_mosaic_grid_sizes( ocean_mosaic, nx, ny)
       nlon = nx(1)
       nlat = ny(1)
    else
       call error_mesg('ocean_model_init','x_T, geolon_t, ocn_mosaic_file does not exist in file '//trim(grid_file), FATAL )
    end if

    allocate (Ocean%Global%lon_bnd (nlon+1),     &
              Ocean%Global%lat_bnd (nlat+1),     &
              Ocean%Global%lon     (nlon, nlat), &
              Ocean%Global%lat     (nlon, nlat), &
              Ocean%Global%mask    (nlon, nlat))

    allocate (rmask(nlon,nlat), geo_lont(nlon,nlat), geo_latt(nlon,nlat), &
              geo_lonv(1:nlon+1,1:nlat+1), geo_latv(1:nlon+1,1:nlat+1) )

    layout = (/0,0/)
    call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes(), layout)
    call mpp_define_domains ( (/1,nlon,1,nlat/), layout, Ocean%Domain, name='NULL Ocean')

    select case (grid_version)
    case(0)
       call read_data(grid_file, "geolon_t",      geo_lont, no_domain=.TRUE. )
       call read_data(grid_file, "geolat_t",      geo_latt, no_domain=.TRUE. )
       call read_data(grid_file, "geolon_vert_t", geo_lonv, no_domain=.TRUE. )
       call read_data(grid_file, "geolat_vert_t", geo_latv, no_domain=.TRUE. )
       call read_data(grid_file, "wet",      rmask,     no_domain=.TRUE.)
    case(1)
       allocate (x_vert_t(nlon,nlat,4), y_vert_t(nlon,nlat,4) )
       call read_data(grid_file, "x_T", geo_lont, no_domain=.TRUE. )
       call read_data(grid_file, "y_T", geo_latt, no_domain=.TRUE. )
       call read_data(grid_file, "x_vert_T", x_vert_t, no_domain=.TRUE.)
       call read_data(grid_file, "y_vert_T", y_vert_t, no_domain=.TRUE. )
       geo_lonv(1:nlon,1:nlat) = x_vert_t(1:nlon,1:nlat,1)
       geo_lonv(nlon+1,1:nlat) = x_vert_t(nlon,1:nlat,2)
       geo_lonv(1:nlon,nlat+1) = x_vert_t(1:nlon,nlat,4)
       geo_lonv(nlon+1,nlat+1) = x_vert_t(nlon,nlat,3)
       geo_latv(1:nlon,1:nlat) = y_vert_t(1:nlon,1:nlat,1)
       geo_latv(nlon+1,1:nlat) = y_vert_t(nlon,1:nlat,2)
       geo_latv(1:nlon,nlat+1) = y_vert_t(1:nlon,nlat,4)
       geo_latv(nlon+1,nlat+1) = y_vert_t(nlon,nlat,3)
       deallocate(x_vert_t, y_vert_t)
       call read_data(grid_file, "wet",      rmask,     no_domain=.TRUE.)
    case(2)
       call get_mosaic_tile_grid(tile_file, ocean_mosaic, Ocean%Domain )
       allocate(tmpx(2*nlon+1, 2*nlat+1), tmpy(2*nlon+1, 2*nlat+1) )
       allocate(garea(nlon, nlat))
       call read_data(tile_file, "x", tmpx, no_domain=.TRUE.)
       call read_data(tile_file, "y", tmpy, no_domain=.TRUE.)
       do j = 1, nlat
          do i = 1, nlon
             geo_lont(i,j) = tmpx(i*2,j*2)
             geo_latt(i,j) = tmpy(i*2,j*2)
          end do
       end do
       do j = 1, nlat+1
          do i = 1, nlon+1
             geo_lonv(i,j) = tmpx(i*2-1,j*2-1)
             geo_latv(i,j) = tmpy(i*2-1,j*2-1)
          end do
       end do

       call calc_mosaic_grid_area(geo_lonv*pi/180., geo_latv*pi/180., garea )
       garea = garea/(4*PI*RADIUS*RADIUS)  ! scale the earth are to be 1
       call field_size(grid_file, "aXo_file", siz)
       nfile_axo = siz(2)
       rmask = 0.0
       do n = 1, nfile_axo
          call read_data(grid_file, "aXo_file", axo_file, level=n)
          axo_file = 'INPUT/'//trim(axo_file)
          nxgrid = get_mosaic_xgrid_size(axo_file)
          allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), xgrid_area(nxgrid))
          call get_mosaic_xgrid(aXo_file, i1, j1, i2, j2, xgrid_area)
          do m = 1, nxgrid
             i = i2(m); j = j2(m)
             rmask(i,j) = rmask(i,j) + xgrid_area(m)
          end do
          deallocate(i1, j1, i2, j2, xgrid_area)
       end do
       rmask = rmask/garea

       deallocate(tmpx, tmpy, garea)
    end select

   Ocean%Global%mask = .false.
   where (rmask > 0) Ocean%Global%mask = .true.
   Ocean%Global%lon_bnd = geo_lonv(:,1)*pi/180.
   Ocean%Global%lat_bnd = geo_latv(1,:)*pi/180.
   Ocean%Global%lon = geo_lont*pi/180.
   Ocean%Global%lat = geo_latt*pi/180.
   deallocate(geo_lonv, geo_latv, geo_lont, geo_latt, rmask)

   call mpp_get_compute_domain ( Ocean%Domain, isd, ied, jsd, jed )
   allocate ( Ocean%Data%lon_bnd (isd:ied+1),      &
              Ocean%Data%lat_bnd (jsd:jed+1),      &
              Ocean%Data%lon (isd:ied, jsd:jed),   &
              Ocean%Data%lat (isd:ied, jsd:jed),   &
              Ocean%Data%mask(isd:ied,jsd:jed))

   Ocean%Data%lon_bnd = Ocean%Global%lon_bnd(isd:ied+1)
   Ocean%Data%lat_bnd = Ocean%Global%lat_bnd(jsd:jed+1)
   Ocean%Data%lon     = Ocean%Global%lon(isd:ied,jsd:jed)
   Ocean%Data%lat     = Ocean%Global%lat(isd:ied,jsd:jed)
   Ocean%Data%mask    = .FALSE.

   allocate ( Ocean%t_surf (isd:ied,jsd:jed), &
              Ocean%u_surf (isd:ied,jsd:jed), &
              Ocean%v_surf (isd:ied,jsd:jed), &
              Ocean%frazil (isd:ied,jsd:jed), &
              Ocean%area   (isd:ied,jsd:jed), &
              Ocean%s_surf (isd:ied,jsd:jed), &
              Ocean%sea_lev(isd:ied,jsd:jed))

    Ocean%t_surf  = 280.
    Ocean%u_surf  = 0.0
    Ocean%v_surf  = 0.0
    Ocean%frazil  = 0.0
    Ocean%area    = 1.0
    Ocean%s_surf  = 0.0
    Ocean%sea_lev = 0.0

 end subroutine ocean_model_init

!#######################################################################

  subroutine ocean_model_end (Ocean, Ocean_state, Time_in)

  type(ocean_state_type),            pointer    :: Ocean_state
  type(time_type),                   intent(in) :: Time_in
  type(ocean_public_type), optional, intent(in) :: Ocean

     ! dummy routine

  end subroutine ocean_model_end

!#######################################################################
! <SUBROUTINE NAME="ocean_model_restart">
!
! <DESCRIPTION>
! dummy interface.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
  subroutine ocean_model_restart(Ocean_state, timestamp)
     type(ocean_state_type),    pointer     :: Ocean_state
     character(len=*), intent(in), optional :: timestamp

     ! dummy routine

  end subroutine ocean_model_restart
! </SUBROUTINE> NAME="ocean_model_restart"

!#######################################################################
  subroutine ocean_model_init_sfc(Ocean_state, Ocean)
  type(ocean_state_type),  pointer    :: Ocean_state    
  type(ocean_public_type), intent(in) :: Ocean

  

  end subroutine ocean_model_init_sfc
!#######################################################################
  subroutine ocean_model_flux_init(Ocean_state)
  type(ocean_state_type), pointer :: Ocean_state

 

  end subroutine ocean_model_flux_init
!#######################################################################
  subroutine ocean_stock_pe(Ocean_state, index, value, time_index)
  type(ocean_state_type), pointer    :: Ocean_state
  integer,               intent(in)  :: index
  real,                  intent(out) :: value
  integer, optional,     intent(in)  :: time_index

  value = 0.0

  end subroutine ocean_stock_pe
!#######################################################################

subroutine ocean_model_data2D_get(OS,Ocean, name, array2D,isc,jsc)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real, dimension(isc:,jsc:), intent(out):: array2D
  integer                   , intent(in) :: isc,jsc
  
  array2D(isc:,jsc:) = 0.0
  
end subroutine ocean_model_data2D_get
!#######################################################################

subroutine ocean_model_data1D_get(OS,Ocean, name, value)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real                      , intent(out):: value

  value = 0.0

end subroutine ocean_model_data1D_get
!#######################################################################

subroutine ice_ocn_bnd_type_chksum(id, timestep, Ice_ocean_boundary)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_ocean_boundary_type), intent(in) :: Ice_ocean_boundary
    return
end subroutine ice_ocn_bnd_type_chksum
!#######################################################################

subroutine ocean_public_type_chksum(id, timestep, Ocean)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ocean_public_type), intent(in) :: Ocean
    return
end subroutine ocean_public_type_chksum
!#######################################################################

end module ocean_model_mod
