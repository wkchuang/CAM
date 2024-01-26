module cb24cnn

!=====================================================================
!
! Purpose: Implement cb24cnn of the model state of U,V,T,Q, and/or PS
!          toward specified values from random analyses increments.
!
! Author: Wiliiam Chapman -- wchapman@ucar.edu
!
! Description:
!
!    This module assumes that the user has {U,V,T,Q,PS} values from analyses
!    which have been preprocessed onto the current model grid and adjusted
!    for differences in topography. It is also assumed that these resulting
!    values and are stored in individual files which are indexed with respect
!    to year, month, day, and second of the day. When the model is inbetween
!    the given begining and ending times, a relaxation forcing is added to
!    cb24cnn the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some beans analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the cb24cnn
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the cb24cnn module is set up to accomodate cb24cnn of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The cb24cnn of the model toward the analyses data is controlled by
!    the 'cb24cnn_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which cb24cnn is applied, the strength of the cb24cnn
!    tendencies, and its spatial distribution.
!

! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_curr_calday,get_step_size
  use, intrinsic :: ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_all, ieee_support_halting, ieee_overflow
  use phys_grid   ,   only:scatter_field_to_chunk, gather_chunk_to_field
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc, mstrid=>masterprocid, mpicom
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  use forpy_mod,                 only : module_py,list,ndarray,object,tuple
  use forpy_mod,                 only : err_print
  use forpy_mod,                 only : forpy_initialize,get_sys_path,import_py,print_py
  use forpy_mod,                 only : ndarray_create,tuple_create,call_py,cast
  use forpy_mod,                 only : forpy_finalize, dict, dict_create

! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: cb24cnn_Model!,cb24cnn_Model_Regress,cb24cnn_ON
  !public:: cb24cnn_readnl
  public:: cb24cnn_init
  public:: cb24cnn_timestep_init
  public:: cb24cnn_timestep_tend
  public:: cb24cnn_finalize


  !public:: regress_daml_timestep_tend
  !public:: regress_diurnal_daml_timestep_tend
  !public:: cb24cnn_diurnal_timestep_tend

  !private::cb24cnn_readin_weights_all
  !private::cb24cnn_readin_weights
  !private::cb24cnn_readin_weights_diurnal
  !private::cb24cnn_readin_normdict
  !private::cb24cnn_set_PSprofile
  !private::cb24cnn_set_profile


  !  cb24cnn Parameters
  !--------------------
  logical          :: cb24cnn_Model       =.true.
  logical          :: halting_mode(5)
  type(module_py)  :: pymodule



  integer :: &
    upwp_clubb_idx, &          ! Cloud fraction
    taux_idx,       &
    pblh_idx,       &
    ustarwec_idx,   &
    flntc_idx,      &
    fsntoa_idx


  ! cb24cnn State Arrays
  integer cb24cnn_nlon,cb24cnn_nlat,cb24cnn_ncol,cb24cnn_nlev
  real(r8),allocatable:: cb24cnn_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_PSstep(:,:)    !(pcols,begchunk:endchunk)



  !CNN input vars
  real(r8),allocatable:: cb24cnn_Model_V (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_T (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_U (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_Q (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_W (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_UPWP (:,:,:)  !(pcols,pverp,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_TAUX (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_TAUY (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_USTAR (:,:)   !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_TBOT (:,:) !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_PBLH (:,:) !(pcols,begchunk:endchunk)

  !derived type interfaces:

  !> Control structure for Python interface
  type, public :: python_interface ; private
    type(module_py) :: pymodule
    type(list) :: paths
  end type

contains

   !================================================================
  subroutine cb24cnn_init
    !
    ! cb24cnn_INIT: Allocate space and initialize cb24cnn values
    !===============================================================
    use ppgrid        ,only: pver,pverp,pcols,begchunk,endchunk
    use error_messages,only: alloc_err
    use dycore        ,only: dycore_is
    use dyn_grid      ,only: get_horiz_grid_dim_d
    use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
    use cam_history   ,only: addfld
    use shr_const_mod ,only: SHR_CONST_PI
    use filenames     ,only: interpret_filename_spec

    implicit none

    integer :: ierror
    type(list) :: my_list
    type(list) :: paths


    !local Values
    !-------------
    integer  istat,lchnk
    integer  hdim1_d,hdim2_d
    character(len=:), allocatable :: return_string



    ! Allocate Space for cb24 data arrays
    !-----------------------------------------

    !write(iulog,*) 'cb24cnn '
    ! Allocate Space for spatial dependence of
    ! cb24cnn Coefs and cb24cnn Forcing.
    !-------------------------------------------
    allocate(cb24cnn_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Ustep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Vstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Sstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Qstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_PSstep(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_PSstep',pcols*((endchunk-begchunk)+1)) !currently not used.



    allocate(cb24cnn_Model_U(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_U',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_V(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_V',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_T(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_T',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_Q',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_W(pcols,pverp,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_W',pcols*pverp*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_UPWP(pcols,pverp,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_Q',pcols*pverp*((endchunk-begchunk)+1))

    allocate(cb24cnn_Model_TAUX(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_TAUX',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24cnn_Model_TAUY(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_TAUY',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24cnn_Model_USTAR(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_USTAR',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24cnn_Model_TBOT(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_TBOT',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24cnn_Model_PBLH(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_PBLH',pcols*((endchunk-begchunk)+1)) !currently not used.


    ! Register output fields with the cam history module
     !-----------------------------------------------------
    call addfld( 'cb24cnn_U',(/ 'lev' /),'A','m/s/s'  ,'U DAMLining Tendency')
    call addfld( 'cb24cnn_V',(/ 'lev' /),'A','m/s/s'  ,'V DAMLining Tendency')
    call addfld( 'cb24cnn_T',(/ 'lev' /),'A','K/s'    ,'T DAMLining Tendency')
    call addfld( 'cb24cnn_Q',(/ 'lev' /),'A','kg/kg/s','Q DAMLining Tendency')


    ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     cb24cnn_nlon=hdim1_d
     cb24cnn_nlat=hdim2_d
     cb24cnn_ncol=hdim1_d*hdim2_d
     cb24cnn_nlev=pver


    ! Initialize the forpy module and import the python script "DAMLcnn.py"
     !-----------------------------------------------------
    if (masterproc) write(iulog,*) 'cb24cnn_init starting'
    if (.not. ieee_support_halting(ieee_overflow)) then
       call endrun("ieee_halting is not supported")
    endif
    call ieee_get_halting_mode(ieee_all, halting_mode)
    print *,__FILE__,__LINE__,halting_mode
    call ieee_set_halting_mode(ieee_all, .false.)

    ierror = forpy_initialize()
    ierror = get_sys_path(paths)
    ierror = paths%append(".") !stash it in the run directory; future NL

    ierror = import_py(pymodule,"DAMLcnn") !future NL
    call ieee_set_halting_mode(ieee_all, halting_mode)
    call paths%destroy

    ! Init forcing as zeros:
    do lchnk=begchunk,endchunk
      cb24cnn_Ustep(:pcols,:pver,lchnk)=0._r8
      cb24cnn_Vstep(:pcols,:pver,lchnk)=0._r8
      cb24cnn_Sstep(:pcols,:pver,lchnk)=0._r8
      cb24cnn_Qstep(:pcols,:pver,lchnk)=0._r8
      cb24cnn_PSstep(:pcols,lchnk)=0._r8
    end do

  end subroutine cb24cnn_init



  subroutine cb24cnn_timestep_init(phys_state,phys_tend,pbuf,cam_in,cam_out)

    !
    ! cb24cnn_TIMESTEP_INIT:
    !                 Check the current time and update Model/Nudging
    !                 arrays when necessary. Toggle the Nudging flag
    !                 when the time is withing the nudging window.
    !===============================================================
    use physconst    ,only: cpair
    use physics_types,only: physics_state,physics_ptend,physics_ptend_init, &
                            physics_state_copy
    use constituents ,only: cnst_get_ind
    use dycore       ,only: dycore_is
    use ppgrid       ,only: pver,pverp,pcols,begchunk,endchunk
    use filenames    ,only: interpret_filename_spec
    use ESMF
    use camsrfexch,     only: cam_in_t,cam_out_t
    use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              physics_buffer_desc


    ! Arguments
    !-------------
    type(physics_state), intent(in) :: phys_state
    type(cam_in_t),      intent(in)    :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    type(physics_ptend), intent(out):: phys_tend
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state) :: state1


    ! Local values
    !----------------
    integer Year,Month,Day,Sec
    integer YMD1,YMD2,YMD,YMD3,YMD4,YMD5,YMD6
    logical Update_Model,Update_Nudge,Sync_Error
    logical After_Beg   ,Before_End
    integer lchnk,ncol,indw,itim_old

    real(r8), pointer, dimension(:,:) :: upwp      ! upwp                               [no idea]
    real(r8), pointer, dimension(:) :: pblh     ! planetary boundary layer height       [m]
    real(r8), pointer, dimension(:) :: taux     ! Surface zonal momentum flux           [m]
    real(r8), pointer, dimension(:) :: tauy     ! Surface meridional momentum flux
    real(r8), pointer, dimension(:) :: ustarwec ! Surface meridional momentum flux

    ncol =phys_state%ncol

    itim_old = pbuf_old_tim_idx()

    upwp_clubb_idx     = pbuf_get_index('UPWP')         ! UPWP
    call pbuf_get_field(pbuf, upwp_clubb_idx,    upwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

    ustarwec_idx     = pbuf_get_index('ustarwec')         ! Ustar !16|33
    call pbuf_get_field(pbuf, ustarwec_idx,    ustarwec)

    pblh_idx     = pbuf_get_index('pblh')         ! PBL height
    call pbuf_get_field(pbuf, pblh_idx,    pblh)

    !cam_out%tbot
    !cam_in%wsx
    !cam_in%wsy



    ! Load values at Current into the Model arrays Predictors...
    !-----------------------------------------------
    ! look here: you need to update this:
    !/glade/work/wchapman/cesm/perfect_nudging/f.e21.FHIST.f09_f09_mg17_nudge_it_ERA5_1999/SourceMods/src.cam
    !nudging_calc_tend()

    call cnst_get_ind('Q',indw)

    cb24cnn_Model_V(:ncol,:pver,lchnk)=state1%u(:ncol,:pver)
    cb24cnn_Model_T(:ncol,:pver,lchnk)=state1%v(:ncol,:pver)
    cb24cnn_Model_U(:ncol,:pver,lchnk)=state1%t(:ncol,:pver)
    cb24cnn_Model_Q(:ncol,:pver,lchnk)=state1%q(:ncol,:pver,indw)
    cb24cnn_Model_W(:ncol,:pver,lchnk)=state1%omega(:ncol,:pver)
    cb24cnn_Model_UPWP(:ncol,:pverp,lchnk)=upwp(:ncol,:pverp)
    cb24cnn_Model_TAUX(:ncol,lchnk)=cam_in%wsx(:ncol)
    cb24cnn_Model_TAUY(:ncol,lchnk)=cam_in%wsy(:ncol)
    cb24cnn_Model_TBOT(:ncol,lchnk)=cam_out%tbot(:ncol)
    cb24cnn_Model_PBLH(:ncol,lchnk)=pblh(:ncol)
    cb24cnn_Model_USTAR(:ncol,lchnk)=ustarwec(:ncol)


  end subroutine cb24cnn_timestep_init



  subroutine cb24cnn_timestep_tend(phys_state,phys_tend,pbuf,cam_in,cam_out)

    use physics_types,only: physics_state,physics_ptend,physics_ptend_init,physics_state_copy
    use camsrfexch,     only: cam_in_t,cam_out_t
    use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              physics_buffer_desc


    !analagous to forpy_run_python

    integer :: ierror
    type(tuple) :: args
    type(dict) :: kwargs
    type(object) :: return_value

    ! Arguments
    !-------------
    type(physics_state), intent(in) :: phys_state
    type(cam_in_t),      intent(in)    :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    type(physics_ptend), intent(out):: phys_tend
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state) :: state1                ! Local copy of state variable

    !local values:
    !---------------
    integer lev
    integer nlon,nlat,plev,istat
    integer ncid,varid
    integer ilat,ilon,ilev
    real(r8) Xcnn(cb24cnn_nlon,cb24cnn_nlat,11)  !future NL
    real(r8) Vanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Uanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Tanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Qanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Wanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) UPWPanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) TAUXanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) TAUYanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) USTARanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) TBOTanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) PBLHanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) Lat_anal(cb24cnn_nlat)
    real(r8) Lon_anal(cb24cnn_nlon)
    real(r8) Xtrans(cb24cnn_nlon,cb24cnn_nlev,cb24cnn_nlat)
    integer  nn,Nindex
    real(r8) PI
    character(len=:), allocatable :: return_string

    if (masterproc) write(iulog,*) 'cb24cnn_timestep_tend starting....'

    PI    = 4.D0*DATAN(1.D0)
    !example of gather chunk
    !call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,state1%u,Uanal)


    !All your python stuff is here:
    call ieee_set_halting_mode(ieee_all, .false.)

    ierror = tuple_create(args, 3)
    ierror = args%setitem(0, 12)
    ierror = args%setitem(1, "Hi")
    ierror = args%setitem(2, .true.)

    ierror = dict_create(kwargs)
    ierror = kwargs%setitem("message", "Hello world!")

    ierror = call_py(return_value,pymodule,"DAMLcnn_run", args, kwargs)
    ierror = cast(return_string, return_value)
    if (masterproc) write(iulog,*)'cb24cnn_init ending ',return_string
    call ieee_set_halting_mode(ieee_all, halting_mode)

    call args%destroy
    call kwargs%destroy
    call return_value%destroy
    !All your python stuff is here:



  !write(iulog,*) 'done with that stuff, again.....'
   ! End Routine
   !------------
  end subroutine cb24cnn_timestep_tend ! cb24cnn_init


  subroutine cb24cnn_finalize()
    call pymodule%destroy
    call forpy_finalize
  end subroutine cb24cnn_finalize
  end module cb24cnn
