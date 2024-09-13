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
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL,r4=>SHR_KIND_R4
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

  use ftorch

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
  public:: cb24cnn_set_tend
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
  logical          :: skip_first          =.true.
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
  integer cb24cnn_nlon,cb24cnn_nlat,cb24cnn_ncol,cb24cnn_nlev,cb24cnn_nlevp
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
  !FTORCH vectors:
  type(torch_model) :: model_ftorch1 !ftorch
  type(torch_model) :: model_ftorch2 !ftorch
  type(torch_model) :: model_ftorch3 !ftorch
  type(torch_model) :: model_ftorch4 !ftorch
  type(torch_model) :: model_ftorch5 !ftorch
  type(torch_model) :: model_ftorch6 !ftorch
  type(torch_model) :: model_ftorch7 !ftorch
  type(torch_model) :: model_ftorch8 !ftorch
  type(torch_model) :: model_ftorch9 !ftorch
  type(torch_model) :: model_ftorch10 !ftorch
  type(torch_model) :: model_ftorch11 !ftorch
  type(torch_model) :: model_ftorch12 !ftorch
  type(torch_model) :: model_ftorch13 !ftorch
  type(torch_model) :: model_ftorch14 !ftorch
  type(torch_model) :: model_ftorch15 !ftorch
  type(torch_model) :: model_ftorch16 !ftorch
  type(torch_model) :: model_ftorch17 !ftorch
  type(torch_model) :: model_ftorch18 !ftorch
  type(torch_model) :: model_ftorch19 !ftorch
  type(torch_model) :: model_ftorch20 !ftorch
  type(torch_model) :: model_ftorch21 !ftorch
  type(torch_model) :: model_ftorch22 !ftorch
  type(torch_model) :: model_ftorch23 !ftorch

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
    use mpi

    implicit none

    integer :: ierror
    character(len=:), allocatable :: Torchpath1
    character(len=:), allocatable :: Torchpath2
    character(len=:), allocatable :: Torchpath3
    character(len=:), allocatable :: Torchpath4
    character(len=:), allocatable :: Torchpath5
    character(len=:), allocatable :: Torchpath6
    character(len=:), allocatable :: Torchpath7
    character(len=:), allocatable :: Torchpath8
    character(len=:), allocatable :: Torchpath9
    character(len=:), allocatable :: Torchpath10
    character(len=:), allocatable :: Torchpath11
    character(len=:), allocatable :: Torchpath12
    character(len=:), allocatable :: Torchpath13
    character(len=:), allocatable :: Torchpath14
    character(len=:), allocatable :: Torchpath15
    character(len=:), allocatable :: Torchpath16
    character(len=:), allocatable :: Torchpath17
    character(len=:), allocatable :: Torchpath18
    character(len=:), allocatable :: Torchpath19
    character(len=:), allocatable :: Torchpath20
    character(len=:), allocatable :: Torchpath21
    character(len=:), allocatable :: Torchpath22
    character(len=:), allocatable :: Torchpath23

    !local Values
    !-------------
    integer  istat,lchnk
    integer  hdim1_d,hdim2_d
    character(len=:), allocatable :: return_string

    ! MPI configuration
    integer :: rank, ierr, i
    
    ! Allocate Space for cb24 data arrays
    !-----------------------------------------

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
     cb24cnn_nlevp=pverp

    call mpi_comm_rank(mpi_comm_world, rank, ierr)

    ! Initialize the forpy module and import the python script "DAMLcnn.py"
     !-----------------------------------------------------
    if (masterproc) write(iulog,*) 'CHACHI cb24cnn_init starting'
    !Load the Pytorch Model
    if(masterproc) then
    Torchpath1 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_87.82123029232025mb_GPU.pt'
    call torch_model_load(model_ftorch1, Torchpath1 ,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath2 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_103.31712663173676mb_GPU.pt'
    call torch_model_load(model_ftorch2, Torchpath2 ,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath3 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_121.54724076390266mb_GPU.pt'
    call torch_model_load(model_ftorch3, Torchpath3,device_type=torch_kCUDA,device_index=1)
    if (masterproc) write(iulog,*) 'CHACHI Loaded 2'
    !Load the Pytorch Model
    Torchpath4 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_142.99403876066208mb_GPU.pt'
    call torch_model_load(model_ftorch4, Torchpath4,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath5 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_168.22507977485657mb_GPU.pt'
    call torch_model_load(model_ftorch5, Torchpath5,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath6 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_197.9080867022276mb_GPU.pt'
    call torch_model_load(model_ftorch6, Torchpath6,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath7 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_232.82861895859241mb_GPU.pt'
    call torch_model_load(model_ftorch7, Torchpath7,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath8 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_273.9108167588711mb_GPU.pt'
    call torch_model_load(model_ftorch8, Torchpath8,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath9 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_322.2419023513794mb_GPU.pt'
    call torch_model_load(model_ftorch9, Torchpath9,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath10 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_379.10090386867523mb_GPU.pt'
    call torch_model_load(model_ftorch10, Torchpath10,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath11 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_445.992574095726mb_GPU.pt'
    call torch_model_load(model_ftorch11, Torchpath11,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath12 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_524.6871747076511mb_GPU.pt'
    call torch_model_load(model_ftorch12, Torchpath12,device_type=torch_kCUDA,device_index=1)
    !Load the Pytorch Model
    Torchpath13 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_609.7786948084831mb_GPU.pt'
    call torch_model_load(model_ftorch13, Torchpath13,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath14 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_691.3894303143024mb_GPU.pt'
    call torch_model_load(model_ftorch14, Torchpath14,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath15 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_763.404481112957mb_GPU.pt'
    call torch_model_load(model_ftorch15, Torchpath15,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath16 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_820.8583686500788mb_GPU.pt'
    call torch_model_load(model_ftorch16, Torchpath16,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath17 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_859.5347665250301mb_GPU.pt'
    call torch_model_load(model_ftorch17, Torchpath17,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath18 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_887.0202489197254mb_GPU.pt'
    call torch_model_load(model_ftorch18, Torchpath18,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath19 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_912.644546944648mb_GPU.pt'
    call torch_model_load(model_ftorch19, Torchpath19,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath20 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_936.1983984708786mb_GPU.pt'
    call torch_model_load(model_ftorch20, Torchpath20,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath21 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_957.485479535535mb_GPU.pt'
    call torch_model_load(model_ftorch21, Torchpath21,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath22 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_976.325407391414mb_GPU.pt'
    call torch_model_load(model_ftorch22, Torchpath22,device_type=torch_kCUDA,device_index=2)
    !Load the Pytorch Model
    Torchpath23 = '/glade/work/wchapman/DA_ML/CESML_AI/Convert_To_Ftorch/UNET_lev_992.556095123291mb_GPU.pt'
    call torch_model_load(model_ftorch23, Torchpath23,device_type=torch_kCUDA,device_index=2)
    if (masterproc) write(iulog,*) 'CHACHI loaded models'
    end if
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
    use camsrfexch     ,only: cam_in_t,cam_out_t
    use physics_buffer ,only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              physics_buffer_desc
    use cam_history  ,only: outfld

    ! Arguments
    !-------------
    type(physics_state), intent(in) :: phys_state
    type(cam_in_t),      intent(in)    :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    type(physics_ptend), intent(out):: phys_tend
    type(physics_buffer_desc), pointer :: pbuf(:)
    
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

    !Grab Chunk and Columns
    !----------------------
    lchnk = phys_state%lchnk
    ncol  = phys_state%ncol

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

    cb24cnn_Model_V(:ncol,:pver,lchnk)=phys_state%v(:ncol,:pver)
    cb24cnn_Model_T(:ncol,:pver,lchnk)=phys_state%t(:ncol,:pver)
    cb24cnn_Model_U(:ncol,:pver,lchnk)=phys_state%u(:ncol,:pver)
    cb24cnn_Model_Q(:ncol,:pver,lchnk)=phys_state%q(:ncol,:pver,indw)
    cb24cnn_Model_W(:ncol,:pver,lchnk)=phys_state%omega(:ncol,:pver)
    cb24cnn_Model_UPWP(:ncol,:pverp,lchnk)=upwp(:ncol,:pverp)
    cb24cnn_Model_TAUX(:ncol,lchnk)=cam_in%wsx(:ncol)
    cb24cnn_Model_TAUY(:ncol,lchnk)=cam_in%wsy(:ncol)
    cb24cnn_Model_TBOT(:ncol,lchnk)=cam_out%tbot(:ncol)
    cb24cnn_Model_PBLH(:ncol,lchnk)=pblh(:ncol)
    cb24cnn_Model_USTAR(:ncol,lchnk)=ustarwec(:ncol)


  end subroutine cb24cnn_timestep_init



  subroutine cb24cnn_timestep_tend()

    use physics_types ,only: physics_state,physics_ptend,physics_ptend_init,physics_state_copy
    use physics_buffer,only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              physics_buffer_desc
    use time_manager  ,only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_curr_calday,get_step_size
    use ppgrid        ,only: pver,pverp,pcols,begchunk,endchunk
    use cam_history   ,only: outfld
    use netcdf
    
    ! Arguments
    !-------------
     !+++ Ftorch tensors
    type(torch_tensor), dimension(1) :: in_tensor, out_tensor !ftorch
    !--- Ftorch tensors
    
    !local values:
    !---------------
    integer lev
    integer c                               
    integer nlon,nlat,plev,istat,lchnk,indw
    integer ncid,varid
    integer ilat,ilon,ilev, jj
    integer londimid, latdimid, levdimid,vardimid,varoutdimid,tsoutdimid,rhvarid, rhvarid2, rhvarid3, rhvarid4, rhvarid5
    integer  Year,Month,Day,Sec,mdo
    real(r8) calday,Hour,windy
    real(r8) Xcnn(cb24cnn_nlon,cb24cnn_nlat,11)  !future NL
    real(r8) Vanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Uanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Tanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Qanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Wanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) UPWPanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlevp)
    real(r8) TAUXanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) TAUYanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) USTARanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) TBOTanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) PBLHanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) Lat_anal(cb24cnn_nlat)
    real(r8) Lon_anal(cb24cnn_nlon)
    real(r8) Xtrans(cb24cnn_nlon,cb24cnn_nlev,cb24cnn_nlat)
    real(r8) Xtransp(cb24cnn_nlon,cb24cnn_nlevp,cb24cnn_nlat)
    real(r8) Xtransf(1,cb24cnn_nlon,cb24cnn_nlat)
    real(r8) model_input_sub(1,13,cb24cnn_nlat,cb24cnn_nlon)
    real(r4) model_input_sub_single(1,13,cb24cnn_nlat,cb24cnn_nlon)
    real(r8) mean_1(13)  ,mean_2(13)  ,mean_3(13)  ,mean_4(13)  ,mean_5(13)  ,mean_6(13)  ,mean_7(13)  ,mean_8(13)  
    real(r8) mean_9(13)  ,mean_10(13)  ,mean_11(13)  ,mean_12(13)  ,mean_13(13)  ,mean_14(13)  ,mean_15(13)  ,mean_16(13)  
    real(r8) mean_17(13)  ,mean_18(13)  ,mean_19(13)  ,mean_20(13)  ,mean_21(13)  ,mean_22(13)  ,mean_23(13)  
    real(r8) std_1(13)  ,std_2(13)  ,std_3(13)  ,std_4(13)  ,std_5(13)  ,std_6(13)  ,std_7(13)  ,std_8(13)  
    real(r8) std_9(13)  ,std_10(13)  ,std_11(13)  ,std_12(13)  ,std_13(13)  ,std_14(13)  ,std_15(13)  ,std_16(13)  
    real(r8) std_17(13)  ,std_18(13)  ,std_19(13)  ,std_20(13)  ,std_21(13)  ,std_22(13)  ,std_23(13) 
    
    integer  nn,Nindex
    real(r8) PI

    !+++ Ftorch vectors
    integer, parameter :: in_dims = 4
    integer, parameter :: n = 13 !for softmax function
    integer :: in_shape(in_dims) = [1, 13, 192, 288]
    integer :: in_layout(in_dims) = [1, 2, 3, 4]
    integer, parameter :: out_dims = 4
    integer :: out_shape(out_dims) =  [1, 2, 192, 288]
    integer :: out_layout(out_dims) = [1, 2, 3, 4]
    integer, parameter :: n_inputs = 1
    real(r4), dimension(:, :, :, :), allocatable, target :: out_data !ftorch
    real(r8), dimension(:, :, :, :), allocatable, target :: out_data_double !ftorch
    real(r8), dimension(:, :, :, :), allocatable, target :: out_for !ftorch

    allocate(out_data(out_shape(1), out_shape(2),out_shape(3), out_shape(4)))
    allocate(out_data_double(1, out_shape(2), out_shape(3), out_shape(4)))
    allocate(out_for(2, 32, cb24cnn_nlat, cb24cnn_nlon))
    !--- Ftorch vectors

    !+++ Scaling Parameters
    mean_1 = [1.12751633e-02, 2.07608670e+02, 8.98119628e+00, &
                             2.56394071e-06, -2.41556448e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, 1.11329661e-05, &
                             2.77689311e+02, 5.95960542e+02, 8.85960407e-07, &
                             4.54997050e-07]
 
    std_1   = [7.10722233e+00, 1.24325691e+01, 1.30211619e+01, &
                             5.31865675e-07, 3.72748981e-02, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 5.44554598e-03, &
                             2.14271024e+01, 1.40864301e+03, 1.66069978e-04, &
                             1.78645431e-04]

    mean_2   = [1.58338996e-02, 2.08453720e+02, 1.03586555e+01, &
                             2.73094280e-06, 7.16691463e-06, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, 3.56053357e-05, &
                             2.77689311e+02, 5.95960542e+02, 4.91962505e-07, &
                             1.67221229e-07]
 
    std_2   = [7.56702052e+00, 1.21573871e+01, 1.37000098e+01, &
                             6.67411396e-07, 4.70150138e-02, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 1.16622374e-02, &
                             2.14271024e+01, 1.40864301e+03, 1.83074237e-04, &
                             2.01923169e-04]
                              
    mean_3   = [6.48410039e-03, 2.10219790e+02, 1.17343041e+01, &
                             3.32806817e-06, 3.93875646e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, 1.31653594e-04, &
                             2.77689311e+02, 5.95960542e+02, 1.71883758e-07, &
                             1.38802747e-07]
 
    std_3   = [8.24844626e+00, 1.10619296e+01, 1.47538158e+01, &
                             1.37038984e-06, 6.16590402e-02, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 3.05040308e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.00203343e-04, &
                             2.22724352e-04]

    mean_4   = [-1.60199771e-02, 2.12481925e+02, 1.28207477e+01, &
                             5.07180904e-06, 7.23543676e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, 5.43752353e-04, &
                             2.77689311e+02, 5.95960542e+02, -3.88493424e-08, &
                             1.58336756e-07]
 
    std_4   = [9.11539397e+00, 9.46257147e+00, 1.59010744e+01, &
                             3.66261665e-06, 8.24857796e-02, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 6.38417943e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.09930768e-04, &
                             2.32939072e-04]
                             
    mean_5   = [-4.48441400e-02, 2.14792942e+02, 1.35149381e+01, &
                             9.95051010e-06, 9.35198095e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, 9.92421562e-04, &
                             2.77689311e+02, 5.95960542e+02, 1.21403732e-07, &
                             9.38950870e-08]
 
    std_5   = [1.02188363e+01, 7.83061980e+00, 1.69471426e+01, &
                             9.69081709e-06, 1.10970587e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 9.17813517e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.18138940e-04, &
                             2.40748560e-04]

    mean_6   = [-4.73616676e-02, 2.17162661e+02, 1.37404746e+01, &
                             2.21207140e-05, 9.03779418e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, 7.62789906e-04, &
                             2.77689311e+02, 5.95960542e+02, 1.72015407e-07, &
                             2.05504195e-07]
 
    std_6   = [1.15526922e+01, 7.02165404e+00, 1.77971032e+01, &
                             2.43868228e-05, 1.49894167e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 9.34762906e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.29202518e-04, &
                             2.50572977e-04]

    mean_7   = [-4.88363673e-02, 2.20294966e+02, 1.33167948e+01, &
                             4.88314761e-05, 8.70787635e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -6.54157942e-04, &
                             2.77689311e+02, 5.95960542e+02, 3.16288631e-08, &
                             -7.29561245e-07]
 
    std_7   = [1.27127929e+01, 7.90688572e+00, 1.81418927e+01, &
                             5.90035836e-05, 2.00972495e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 6.54484445e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.33696642e-04, &
                             2.53048089e-04]

    mean_8   = [-4.74566175e-02, 2.24976132e+02, 1.22800173e+01, &
                             1.02081443e-04, 8.58252238e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -7.92602832e-04, &
                             2.77689311e+02, 5.95960542e+02, -2.50531979e-07, &
                             -1.37943113e-06]
 
    std_8   = [1.31771803e+01, 9.95660740e+00, 1.77858716e+01, &
                             1.32852573e-04, 2.60711512e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 5.79873934e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.42609347e-04, &
                             2.57763742e-04]


    mean_9   = [-3.02319180e-02, 2.31175244e+02, 1.08364086e+01, &
                             1.98063165e-04, 7.61023099e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -7.92602832e-04, &
                             2.77689311e+02, 5.95960542e+02, -1.25304287e-08, &
                             -1.63272735e-06]
 
    std_9   = [1.27194156e+01, 1.19972578e+01, 1.66960684e+01, &
                             2.71702714e-04, 3.21215358e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 5.79873934e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.47902509e-04, &
                             2.58753110e-04]
    
    mean_10   = [-1.10879880e-02, 2.38422402e+02, 9.23175704e+00, &
                             3.59904675e-04, 5.57535762e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -4.48603696e-04, &
                             2.77689311e+02, 5.95960542e+02, 2.76314659e-07, &
                             -1.60830636e-06]
 
    std_10   = [1.15996081e+01, 1.35364483e+01, 1.51364458e+01, &
                             5.07841002e-04, 3.72986997e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 4.22229813e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.39123353e-04, &
                             2.46867585e-04]

    mean_11   = [-3.86493706e-03, 2.46085920e+02, 7.60693507e+00, &
                             6.30048400e-04, 3.92355931e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -5.13951078e-04, &
                             2.77689311e+02, 5.95960542e+02, 7.65470920e-09, &
                             -9.79059107e-07]
 
    std_11   = [1.02094833e+01, 1.45268035e+01, 1.34309832e+01, &
                             8.79856883e-04, 4.07451398e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 4.02430729e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.21941827e-04, &
                             2.28224709e-04]
    
    mean_12   = [-2.23871343e-02, 2.53729682e+02, 6.08154861e+00, &
                             1.11928219e-03, 1.63689110e-05, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -7.47261182e-04, &
                             2.77689311e+02, 5.95960542e+02, -3.08401325e-07, &
                             -5.59135593e-07]
 
    std_12   = [8.83169701e+00, 1.50602461e+01, 1.17509268e+01, &
                             1.48128336e-03, 4.21177889e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 4.20755127e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.06057773e-04, &
                             2.12897158e-04]

    mean_13   = [-2.52487011e-02, 2.60423990e+02, 4.72581150e+00, &
                             1.71492938e-03, 1.38987625e-04, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -1.18231781e-03, &
                             2.77689311e+02, 5.95960542e+02, -2.86801728e-07, &
                             -3.36355304e-07]
 
    std_13   = [7.73537996e+00, 1.52837394e+01, 1.03674458e+01, &
                             1.99786834e-03, 4.22671587e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 6.09855000e-02, &
                             2.14271024e+01, 1.40864301e+03, 1.98750795e-04, &
                             2.05089449e-04]

    mean_14   = [-1.68129368e-02, 2.65773228e+02, 3.56617344e+00, &
                             2.40875539e-03, 2.10585971e-04, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -1.18231781e-03, &
                             2.77689311e+02, 5.95960542e+02, 1.04935775e-08, &
                             -4.51134672e-07]
 
    std_14   = [6.98206929e+00, 1.55853172e+01, 9.38453003e+00, &
                            2.56986950e-03, 4.18257061e-01, 3.84669990e-01, &
                            3.24386126e-01, 1.44692706e-01, 6.09855000e-02, &
                            2.14271024e+01, 1.40864301e+03, 2.02893813e-04, &
                             2.06726218e-04]

    mean_15   = [2.83461360e-03, 2.69604391e+02, 2.59273760e+00, &
                             3.18320962e-03, 1.02605344e-04, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -2.16956370e-03, &
                             2.77689311e+02, 5.95960542e+02, 2.97642283e-07, &
                             -3.27894633e-07]
 
    std_15   = [6.54419371e+00, 1.59749548e+01, 8.75288661e+00, &
                             3.07517437e-03, 4.09739112e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 8.75332358e-02, &
                             2.14271024e+01, 1.40864301e+03, 2.15455405e-04, &
                             2.18679274e-04]

    mean_16   = [2.75240649e-02, 2.71847169e+02, 1.83469388e+00, &
                             4.17634520e-03, -2.19614465e-04, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -4.39051210e-03, &
                             2.77689311e+02, 5.95960542e+02, 3.16367452e-07, &
                             1.44043399e-07]
 
    std_16   = [6.37842976e+00, 1.60728734e+01, 8.41841975e+00, &
                             3.69050395e-03, 3.96275145e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 1.15484776e-01, &
                             2.14271024e+01, 1.40864301e+03, 2.32427130e-04, &
                             2.34647741e-04]

    mean_17   = [4.81732752e-02, 2.73308828e+02, 1.32991844e+00, &
                             4.98709825e-03, -6.31851302e-04, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -9.80446146e-03, &
                             2.77689311e+02, 5.95960542e+02, 2.71214258e-07, &
                             3.30176503e-07]
 
    std_17   = [6.41962378e+00, 1.62584890e+01, 8.30943295e+00, &
                             4.25643811e-03, 3.79292000e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 1.50485547e-01, &
                             2.14271024e+01, 1.40864301e+03, 2.43919859e-04, &
                             2.47305117e-04]

    mean_18   = [6.45751488e-02, 2.74389077e+02, 9.85557550e-01, &
                             5.58614792e-03, -9.77205799e-04, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -1.09120129e-02, &
                             2.77689311e+02, 5.95960542e+02, 2.30710282e-07, &
                             4.92056744e-07]
 
    std_18   = [6.54457019e+00, 1.65295661e+01, 8.28523557e+00, &
                                 4.70803391e-03, 3.56891940e-01, 3.84669990e-01, &
                                 3.24386126e-01, 1.44692706e-01, 1.55942005e-01, &
                                 2.14271024e+01, 1.40864301e+03, 2.52736383e-04, &
                                 2.57875571e-04]
                             
    mean_19   = [7.77324687e-02, 2.75446864e+02, 6.79634757e-01, &
                             6.11801559e-03, -1.17622864e-03, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -1.09120129e-02, &
                             2.77689311e+02, 5.95960542e+02, 2.20294757e-07, &
                             8.08039070e-07]
 
    std_19   = [6.71011073e+00, 1.69478247e+01, 8.26973995e+00, &
                             5.13423431e-03, 3.29282373e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 1.55942005e-01, &
                             2.14271024e+01, 1.40864301e+03, 2.65767402e-04, &
                             2.72557290e-04]

    mean_20   = [9.47966007e-02, 2.76437403e+02, 4.09897939e-01, &
                             6.55689088e-03, -7.56357412e-04, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -1.09807339e-02, &
                             2.77689311e+02, 5.95960542e+02, 2.99895438e-07, &
                             1.36400135e-06]
 
    std_20   = [6.85373185e+00, 1.75246786e+01, 8.19790299e+00, &
                             5.52769798e-03, 3.01095492e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 1.61446598e-01, &
                             2.14271024e+01, 1.40864301e+03, 2.84940946e-04, &
                             2.93821028e-04]

    mean_21   = [1.27872690e-01, 2.77290153e+02, 1.88416093e-01, &
                             6.86851890e-03, 1.29395185e-03, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -8.92564820e-03, &
                             2.77689311e+02, 5.95960542e+02, 1.05891231e-06, &
                             1.91294436e-06]
 
    std_21   = [6.86795610e+00, 1.82960872e+01, 7.97747544e+00, &
                             5.81981324e-03, 2.83281366e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 2.04315161e-01, &
                             2.14271024e+01, 1.40864301e+03, 3.08772309e-04, &
                             3.19311561e-04]

    mean_22   = [1.79077334e-01, 2.77729799e+02, 2.79393881e-02, &
                             7.06290748e-03, 7.10115067e-03, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -8.92564820e-03, &
                             2.77689311e+02, 5.95960542e+02, 2.71649776e-06, &
                             1.48217739e-06]
 
    std_22   = [6.59001858e+00, 1.96449407e+01, 7.51399711e+00, &
                             6.00127275e-03, 3.00380582e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 2.04315161e-01, &
                             2.14271024e+01, 1.40864301e+03, 3.31403127e-04, &
                             3.41056316e-04]

    mean_23   = [2.06466766e-01, 2.77701554e+02, -6.54231915e-02, &
                             7.26872040e-03, 1.52488757e-02, -7.13681192e-03, &
                             -2.73477924e-03, 2.61954097e-01, -7.68694273e-03, &
                             2.77689311e+02, 5.95960542e+02, 1.43239934e-06, &
                             -5.07835787e-07]
 
    std_23   = [5.55351669e+00, 2.14303254e+01, 6.42386535e+00, &
                             6.19860987e-03, 3.39291081e-01, 3.84669990e-01, &
                             3.24386126e-01, 1.44692706e-01, 2.43212858e-01, &
                             2.14271024e+01, 1.40864301e+03, 3.49354670e-04, &
                             3.44564385e-04]

    
    if (masterproc) write(iulog,*) 'cb24cnn_timestep_tend starting....'

    PI    = 4.D0*DATAN(1.D0)
    !Grab Chunk and Columns
    !----------------------
    
    !---
    !V
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_V,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Vanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 

    !---
    !U
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_U,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Uanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 

    !---
    !T
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_T,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Tanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 
    
    !---
    !Q
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_Q,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Qanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 

    !---
    !OMEGA
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_W,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Wanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 
    
    !---
    !UPWP
    !---
    call gather_chunk_to_field(1,cb24cnn_nlevp,1,cb24cnn_nlon,cb24cnn_Model_UPWP,Xtransp)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlevp
        do ilon=1,cb24cnn_nlon
           UPWPanal(ilon,ilat,ilev)=Xtransp(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 

    !---
    !TAUX
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_TAUX,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           TAUXanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
        !if you want to see if you've done good kid. 
        !istat=nf90_create('CESM_dumpTAUX.nc',NF90_CLOBBER,ncid)
        !istat=nf90_def_dim(ncid, "lon",cb24cnn_nlon , londimid)
        !istat=nf90_def_dim(ncid, "lat",cb24cnn_nlat , latdimid)
        !istat=nf90_def_var(ncid, "TAUX", nf90_double, (/ londimid, latdimid /), rhvarid)
        !istat=nf90_enddef(ncid)
        !istat=nf90_put_var(ncid, rhvarid, TAUXanal) 
        !istat=nf90_close(ncid)
        !write(*,*) 'Finished TAUX dumping field.'
    end if

    !---
    !TAUY
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_TAUY,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           TAUYanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
    end if

    !---
    !USTAR
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_USTAR,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           USTARanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
        !if you want to see if you've done good kid. 
        !istat=nf90_create('CESM_dumpUSTAR.nc',NF90_CLOBBER,ncid)
        !istat=nf90_def_dim(ncid, "lon",cb24cnn_nlon , londimid)
        !istat=nf90_def_dim(ncid, "lat",cb24cnn_nlat , latdimid)
        !istat=nf90_def_var(ncid, "USTAR", nf90_double, (/ londimid, latdimid /), rhvarid2)
        !istat=nf90_enddef(ncid)
        !istat=nf90_put_var(ncid, rhvarid2, USTARanal) 
        !istat=nf90_close(ncid)
        !write(*,*) 'Finished USTAR dumping field.'
    end if

    !---
    !TBOT
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_TBOT,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           TBOTanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
        !!if you want to see if you've done good kid. 
        !istat=nf90_create('CESM_dumpTBOT.nc',NF90_CLOBBER,ncid)
        !istat=nf90_def_dim(ncid, "lon",cb24cnn_nlon , londimid)
        !istat=nf90_def_dim(ncid, "lat",cb24cnn_nlat , latdimid)
        !istat=nf90_def_var(ncid, "TBOT", nf90_double, (/ londimid, latdimid /), rhvarid4)
        !istat=nf90_enddef(ncid)
        !istat=nf90_put_var(ncid, rhvarid4, TBOTanal) 
        !istat=nf90_close(ncid)
        !write(*,*) 'Finished TBOT dumping field.'
    end if

    !---
    !PBLH
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_PBLH,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           PBLHanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
        !!if you want to see if you've done good kid. 
        !istat=nf90_create('CESM_dumpPBLH.nc',NF90_CLOBBER,ncid)
        !istat=nf90_def_dim(ncid, "lon",cb24cnn_nlon , londimid)
        !istat=nf90_def_dim(ncid, "lat",cb24cnn_nlat , latdimid)
        !istat=nf90_def_var(ncid, "PBLH", nf90_double, (/ londimid, latdimid /), rhvarid3)
        !istat=nf90_enddef(ncid)
        !istat=nf90_put_var(ncid, rhvarid3, PBLHanal) 
        !istat=nf90_close(ncid)
        !write(*,*) 'Finished USTAR dumping field.'
    end if

    
    
    if (masterproc) write(iulog,*) 'finished gather .....'

    !All your python stuff is here:

    if (masterproc) then

        calday = get_curr_calday()
        calday = 3.0_r8*SIN(2.0_r8*3.14159265359_r8*calday/365.0_r8)
        call get_curr_date(Year,Month,Day,Sec)
        Hour = Sec/3600.0_r8
        Hour = 3.0_r8*COS(2.0_r8*3.14159265359_r8*calday/24.0_r8)

        !+++++++++++++++++++++++++++++
        !set all the data for model 1
        !+++++++++++++++++++++++++++++
        !Vanal = lon,lat,lev
        !TAUXanal = lon,lat
        !mdo=10
        !call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
        !                       Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
        !                       UPWPanal, TBOTanal, PBLHanal, mean_1, std_1, calday, Hour, mdo)

        !!if you want to see if you've done good kid. 
        !istat=nf90_create('preddy.nc',NF90_CLOBBER,ncid)
        !istat=nf90_def_dim(ncid, "lon",cb24cnn_nlon , londimid)
        !istat=nf90_def_dim(ncid, "lat",cb24cnn_nlat , latdimid)
        !istat=nf90_def_dim(ncid, "uv" ,13            , varoutdimid)
        !istat=nf90_def_dim(ncid, "teto" ,1          , tsoutdimid)
        !istat=nf90_def_var(ncid, "Tender", nf90_double, (/tsoutdimid,varoutdimid,latdimid,londimid/), rhvarid4)
        !istat=nf90_enddef(ncid)
        !istat=nf90_put_var(ncid, rhvarid4, model_input_sub) 
        !istat=nf90_close(ncid)
        !write(*,*) 'Finished out_for dumping field.'

        !model_input_sub_single = real(model_input_sub, r4)
        !in_tensor(1) = torch_tensor_from_array(model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        !out_tensor(1) = torch_tensor_from_array(out_data, out_layout, torch_kCPU)
        !call torch_model_forward(model_ftorch1, in_tensor, out_tensor)
        !out_data_double = real(out_data, r8)
        !do ilat=1,cb24cnn_nlat
        !do ilon=1,cb24cnn_nlon
        !  windy = (windowingFactor(ilat))
        !  out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_1(12))+mean_1(12)
        !  out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_1(13))+mean_1(13)
        !end do
        !end do
        !call torch_delete(in_tensor)
        !call torch_delete(out_tensor)

        
        !+++++++++++++++++++++++++++++
        !set all the data for model 2
        !+++++++++++++++++++++++++++++
        !mdo=11
        !call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
        !                       Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
        !                       UPWPanal, TBOTanal, PBLHanal, mean_2, std_2, calday, Hour, mdo)

        !model_input_sub_single = real(model_input_sub, r4)
        !in_tensor(1) = torch_tensor_from_array(model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        !out_tensor(1) = torch_tensor_from_array(out_data, out_layout, torch_kCPU)
        !call torch_model_forward(model_ftorch2, in_tensor, out_tensor)
        !out_data_double = real(out_data, r8)
        !do ilat=1,cb24cnn_nlat
        !do ilon=1,cb24cnn_nlon
        !  windy = (windowingFactor(ilat))
        !  out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_2(12))+mean_2(12)
        !  out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_2(13))+mean_2(13)
        !end do
        !end do
        !call torch_delete(in_tensor)
        !call torch_delete(out_tensor)
        
        !+++++++++++++++++++++++++++++
        !set all the data for model 3
        !+++++++++++++++++++++++++++++
        mdo=12
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_3, std_3, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch3, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_3(12))+mean_3(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_3(13))+mean_3(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 4
        !+++++++++++++++++++++++++++++
        mdo=13
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_4, std_4, calday, Hour, mdo)

        

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch4, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_4(12))+mean_4(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_4(13))+mean_4(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 5
        !+++++++++++++++++++++++++++++
        mdo=14
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_5, std_5, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch5, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_5(12))+mean_5(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_5(13))+mean_5(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 6
        !+++++++++++++++++++++++++++++
        mdo=15
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_6, std_6, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch6, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_6(12))+mean_6(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_6(13))+mean_6(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)


         !+++++++++++++++++++++++++++++
        !set all the data for model 7
        !+++++++++++++++++++++++++++++
        mdo=16
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_7, std_7, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch7, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_7(12))+mean_7(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_7(13))+mean_7(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)


         !+++++++++++++++++++++++++++++
        !set all the data for model 8
        !+++++++++++++++++++++++++++++
        mdo=17
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_8, std_8, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch8, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_8(12))+mean_8(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_8(13))+mean_8(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)


         !+++++++++++++++++++++++++++++
        !set all the data for model 9
        !+++++++++++++++++++++++++++++
        mdo=18
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_9, std_9, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch9, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_9(12))+mean_9(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_9(13))+mean_9(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)


         !+++++++++++++++++++++++++++++
        !set all the data for model 10
        !+++++++++++++++++++++++++++++
        mdo=19
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_10, std_10, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch10, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_10(12))+mean_10(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_10(13))+mean_10(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)


         !+++++++++++++++++++++++++++++
        !set all the data for model 11
        !+++++++++++++++++++++++++++++
        mdo=20
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_11, std_11, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch11, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat)) 
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_11(12))+mean_11(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_11(13))+mean_11(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)


         !+++++++++++++++++++++++++++++
        !set all the data for model 12
        !+++++++++++++++++++++++++++++
        mdo=21
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_12, std_12, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=1)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch12, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_12(12))+mean_12(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_12(13))+mean_12(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 13
        !+++++++++++++++++++++++++++++
        mdo=22
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_13, std_13, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch13, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_13(12))+mean_13(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_13(13))+mean_13(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 14
        !+++++++++++++++++++++++++++++
        mdo=23
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_14, std_14, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch14, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_14(12))+mean_14(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_14(13))+mean_14(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 15
        !+++++++++++++++++++++++++++++
        mdo=24
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_15, std_15, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch15, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_15(12))+mean_15(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_15(13))+mean_15(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 16
        !+++++++++++++++++++++++++++++
        mdo=25
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_16, std_16, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch16, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_16(12))+mean_16(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_16(13))+mean_16(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 17
        !+++++++++++++++++++++++++++++
        mdo=26
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_17, std_17, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch17, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_17(12))+mean_17(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_17(13))+mean_17(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 18
        !+++++++++++++++++++++++++++++
        mdo=27
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_18, std_18, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch18, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_18(12))+mean_18(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_18(13))+mean_18(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 19
        !+++++++++++++++++++++++++++++
        mdo=28
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_19, std_19, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch19, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_19(12))+mean_19(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_19(13))+mean_19(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 20
        !+++++++++++++++++++++++++++++
        mdo=29
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_20, std_20, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch20, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_20(12))+mean_20(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_20(13))+mean_20(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 21
        !+++++++++++++++++++++++++++++
        mdo=30
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_21, std_21, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch21, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_21(12))+mean_21(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_21(13))+mean_21(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

         !+++++++++++++++++++++++++++++
        !set all the data for model 22
        !+++++++++++++++++++++++++++++
        mdo=31
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_22, std_22, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch22, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_22(12))+mean_22(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_22(13))+mean_22(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)

        !+++++++++++++++++++++++++++++
        !set all the data for model 23
        !+++++++++++++++++++++++++++++
        mdo=32
        call prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev, cb24cnn_nlevp, Vanal, Tanal, &
                               Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, &
                               UPWPanal, TBOTanal, PBLHanal, mean_23, std_23, calday, Hour, mdo)

        model_input_sub_single = real(model_input_sub, r4)
        call torch_tensor_from_array(in_tensor(1), model_input_sub_single, in_layout, torch_kCUDA,device_index=2)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_ftorch23, in_tensor, out_tensor)
        out_data_double = real(out_data, r8)
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
          windy = (windowingFactor(ilat))
          out_for(1,mdo,ilat,ilon) = (windy)*(out_data_double(1,1,ilat,ilon)*std_23(12))+mean_23(12)
          out_for(2,mdo,ilat,ilon) = (windy)*(out_data_double(1,2,ilat,ilon)*std_23(13))+mean_23(13)
        end do
        end do
        
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)
        if (masterproc) write(iulog,*)'cb24cnn_init ending '
        !here we need to do a scatter: 

    if(skip_first)then
    do ilat = 1, cb24cnn_nlat
        do ilon = 1, cb24cnn_nlon
            do ilev = 1, cb24cnn_nlev
                do jj = 1, 2
                    if (isnan(out_for(jj,ilev,ilat,ilon))) then
                        out_for(jj,ilev,ilat,ilon)=0.0_r8
                    end if
                    
                    if (ilev < 12) then
                        out_for(jj,ilev,ilat,ilon)=0.0_r8
                    end if
                    out_for(jj,ilev,ilat,ilon)=0.0_r8
                    
                end do
            end do
        end do
    end do
    end if 
    skip_first = .false.
    
    do ilat = 1, cb24cnn_nlat
        do ilon = 1, cb24cnn_nlon
            do ilev = 1, cb24cnn_nlev
                do jj = 1, 2
                    if (isnan(out_for(jj,ilev,ilat,ilon))) then
                        out_for(jj,ilev,ilat,ilon)=0.0_r8
                    end if
                    
                    if (ilev < 12) then
                        out_for(jj,ilev,ilat,ilon)=0.0_r8
                    end if

                    if (out_for(jj,ilev,ilat,ilon) > 0.1) then 
                        out_for(jj,ilev,ilat,ilon)=0.0_r8
                    end if 
                    
                    if (out_for(jj,ilev,ilat,ilon) < -0.1) then 
                        out_for(jj,ilev,ilat,ilon)=0.0_r8
                    end if 
                    
                end do
            end do
        end do
    end do
     
    !All your python stuff is here:
    end if
    
    if (masterproc) then
    do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
            do ilon=1,cb24cnn_nlon
               Xtrans(ilon,ilev,ilat)=out_for(1,ilev,ilat,ilon) !U
            end do
        end do
    end do
    end if
    call scatter_field_to_chunk(1,cb24cnn_nlev,1,cb24cnn_nlon,Xtrans,    &
                                cb24cnn_Ustep(1,1,begchunk))


    if (masterproc) then
    do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
            do ilon=1,cb24cnn_nlon
               Xtrans(ilon,ilev,ilat)=out_for(2,ilev,ilat,ilon) !V
            end do
        end do
    end do
    end if
    call scatter_field_to_chunk(1,cb24cnn_nlev,1,cb24cnn_nlon,Xtrans,    &
                                cb24cnn_Vstep(1,1,begchunk))

    if (masterproc) then
    do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
            do ilon=1,cb24cnn_nlon
               Xtrans(ilon,ilev,ilat)=0.0_r8 !T
            end do
        end do
    end do
    end if
    call scatter_field_to_chunk(1,cb24cnn_nlev,1,cb24cnn_nlon,Xtrans,    &
                                cb24cnn_Sstep(1,1,begchunk))

    
    if (masterproc) then
    do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
            do ilon=1,cb24cnn_nlon
               Xtrans(ilon,ilev,ilat)=0.0_r8 !Q
            end do
       end do
    end do
    end if
    call scatter_field_to_chunk(1,cb24cnn_nlev,1,cb24cnn_nlon,Xtrans,    &
                                cb24cnn_Qstep(1,1,begchunk))
    
    do c=begchunk,endchunk
        !update phys_tend
        call outfld('cb24cnn_U',cb24cnn_Ustep(:,:,c),pcols,c)
        call outfld('cb24cnn_V',cb24cnn_Vstep(:,:,c),pcols,c)
        call outfld('cb24cnn_T',cb24cnn_Sstep(:,:,c),pcols,c)
        call outfld('cb24cnn_Q',cb24cnn_Qstep(:,:,c),pcols,c)
    end do

  !write(iulog,*) 'done with that stuff, again.....'
   ! End Routine
   !------------
   
  end subroutine cb24cnn_timestep_tend 

  subroutine cb24cnn_set_tend(phys_state,phys_tend)
  !
   ! CB24cnn_set_teset:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the Nudge
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'cb24cnn',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   lchnk=phys_state%lchnk
   ncol =phys_state%ncol
   phys_tend%u(:ncol,:pver)     = cb24cnn_Ustep(:ncol,:pver,lchnk)*.35
   phys_tend%v(:ncol,:pver)     = cb24cnn_Vstep(:ncol,:pver,lchnk)*.35
   phys_tend%s(:ncol,:pver)     = cb24cnn_Sstep(:ncol,:pver,lchnk)*.35
   phys_tend%q(:ncol,:pver,indw)= cb24cnn_Qstep(:ncol,:pver,lchnk)*.35
  end subroutine cb24cnn_set_tend 

  subroutine prepareModelInput(model_input_sub, cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev,cb24cnn_nlevp , Vanal, Tanal, Uanal, Qanal, Wanal, TAUXanal, TAUYanal, USTARanal, UPWPanal, TBOTanal, PBLHanal, mean, std, calday, Hour, index)
    ! Assuming model_input_sub is declared in the calling scope or passed as an argument
    ! Declare the types and dimensions for the parameters
    integer, intent(in) :: cb24cnn_nlat, cb24cnn_nlon, cb24cnn_nlev,cb24cnn_nlevp, index
    real(r8), dimension(1,13,cb24cnn_nlat,cb24cnn_nlon), intent(inout) :: model_input_sub
    real(r8), dimension(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev), intent(in) :: Vanal, Tanal, Uanal, Qanal, Wanal
    real(r8), dimension(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlevp), intent(in) :: UPWPanal
    real(r8), dimension(cb24cnn_nlon,cb24cnn_nlat), intent(in) :: TAUXanal, TAUYanal, TBOTanal, PBLHanal, USTARanal
    
    real(r8), dimension(:), intent(in) :: mean, std
    real(r8), intent(in) :: calday, Hour
    integer :: ilat, ilon, j

    ! Populate model_input_sub
    do ilat = 1, cb24cnn_nlat
        do ilon = 1, cb24cnn_nlon                                
            model_input_sub(1,1,ilat,ilon) = (Vanal(ilon,ilat,index) - mean(1)) / std(1)
            model_input_sub(1,2,ilat,ilon) = (Tanal(ilon,ilat,index) - mean(2)) / std(2)
            model_input_sub(1,3,ilat,ilon) = (Uanal(ilon,ilat,index) - mean(3)) / std(3)
            model_input_sub(1,4,ilat,ilon) = (Qanal(ilon,ilat,index) - mean(4)) / std(4)
            model_input_sub(1,5,ilat,ilon) = (Wanal(ilon,ilat,index) - mean(5)) / std(5)
            model_input_sub(1,6,ilat,ilon) = (TAUXanal(ilon,ilat) - mean(6)) / std(6)
            model_input_sub(1,7,ilat,ilon) = (TAUYanal(ilon,ilat) - mean(7)) / std(7)
            model_input_sub(1,8,ilat,ilon) = (USTARanal(ilon,ilat) - mean(8)) / std(8)
            model_input_sub(1,9,ilat,ilon) = (UPWPanal(ilon,ilat,index) - mean(9)) / std(9)
            model_input_sub(1,10,ilat,ilon) = (TBOTanal(ilon,ilat) - mean(10)) / std(10)
            model_input_sub(1,11,ilat,ilon) = (PBLHanal(ilon,ilat) - mean(11)) / std(11)
            model_input_sub(1,12,ilat,ilon) = Hour
            model_input_sub(1,13,ilat,ilon) = calday
        end do
    end do

    do ilat = 1, cb24cnn_nlat
        do ilon = 1, cb24cnn_nlon
            ! Loop through each variable in the model_input_sub array
            do j = 1, 13
                ! Adjust each input based on mean and standard deviation
                ! Placeholder for specific computation per variable
                model_input_sub(1,j,ilat,ilon) = model_input_sub(1,j,ilat,ilon)
                ! Clamp values to the range [-15, 15]
                if (model_input_sub(1,j,ilat,ilon) > 15) then
                    model_input_sub(1,j,ilat,ilon) = 15
                else if (model_input_sub(1,j,ilat,ilon) < -15) then
                    model_input_sub(1,j,ilat,ilon) = -15
                endif
            end do
        end do
    end do


    ! Convert to single precision (if required here or can be moved to the main subroutine)
  end subroutine prepareModelInput


      ! Function to calculate the windowing factor based on index position
    function windowingFactor(index)
        implicit none
        integer, intent(in) :: index
        real(r8) :: windowingFactor
        real(r8) :: scaledIndex
    
        ! Check if index is in the lower ramp down zone (0 to 10)
        if (index <= 15) then
            ! Scale index from 0 at index 0 to -1 at index 10
            scaledIndex = (10.0 - index) / 10.0
            windowingFactor = 0.5 * (1.0 - tanh(scaledIndex))
    
        ! Check if index is in the upper ramp down zone (182 to 192)
        elseif (index >= 177 .and. index <= 192) then
            ! Scale index from 0 at index 192 to 1 at index 182
            scaledIndex = (index - 177.0) / 10.0
            windowingFactor = 0.5 * (1.0 - tanh(scaledIndex))
    
        else
            ! Keep as it is for indices between 11 and 181
            windowingFactor = 1.0
        endif
    
        return
    end function windowingFactor
  
  subroutine cb24cnn_finalize()
    ! Cleanup
    call torch_delete(model_ftorch1)
    call torch_delete(model_ftorch2)
    call torch_delete(model_ftorch3)
    call torch_delete(model_ftorch4)
    call torch_delete(model_ftorch5)
    call torch_delete(model_ftorch6)
    call torch_delete(model_ftorch7)
    call torch_delete(model_ftorch8)
    call torch_delete(model_ftorch9)
    call torch_delete(model_ftorch10)
    call torch_delete(model_ftorch11)
    call torch_delete(model_ftorch12)
    call torch_delete(model_ftorch13)
    call torch_delete(model_ftorch14)
    call torch_delete(model_ftorch15)
    call torch_delete(model_ftorch16)
    call torch_delete(model_ftorch17)
    call torch_delete(model_ftorch18)
    call torch_delete(model_ftorch19)
    call torch_delete(model_ftorch20)
    call torch_delete(model_ftorch21)
    call torch_delete(model_ftorch22)
    call torch_delete(model_ftorch23)
  end subroutine cb24cnn_finalize
  end module cb24cnn
