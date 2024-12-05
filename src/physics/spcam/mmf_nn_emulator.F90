module mmf_nn_emulator

  use constituents,    only: pcnst
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver
  use camsrfexch,      only: cam_out_t, cam_in_t
  use constituents,    only: cnst_get_ind
  use spmd_utils,      only: masterproc
  use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
  use physconst,       only: cpair
  use cam_logfile,     only: iulog
  use phys_grid,       only: get_rlat_p
  use physics_buffer,  only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
  
  !-------- torch fortran --------
  use ftorch
  use iso_fortran_env
  !-------------------------------

  implicit none
  save
  private
  ! Define variables for this module
  ! These are namelist variables
  ! If they are not specific in atm_in, then the following default values are used.
  integer :: inputlength
  integer :: outputlength
  logical :: input_rh = .false.
  character(len=256) :: cb_torch_model

  type(torch_model) :: model_pytorch  

  ! local
  logical :: cb_top_levels_zero_out = .true.
  integer :: cb_n_levels_zero = 12 ! top n levels to zero out
#ifdef MMF_NN_EMULATOR
  public init_neural_net, neural_net, mmf_nn_emulator_readnl
#else
  public mmf_nn_emulator_readnl
#endif /* MMF_NN_EMULATOR */

  contains

!===============================================================================
!===============================================================================
  
#ifdef MMF_NN_EMULATOR 
    subroutine init_neural_net()
      !!! Allocate space for neural network values

      implicit none

      !wkc Will have to return to this for CPU/GPU assessment 
      call torch_model_load(model_pytorch, trim(cb_torch_model))
    
      ! add diagnostic output files
      ! Called from cam_history
      ! call addfld ('TROP_IND',horiz_only,   'A', '1', 'lev index for tropopause')

    end subroutine init_neural_net
#endif /* MMF_NN_EMULATOR */

#ifdef MMF_NN_EMULATOR
    subroutine neural_net (ptend, state, state_aphys1, pbuf, cam_in, cam_out, coszrs, solin, ztodt, lchnk)
      ! note state is meant to have the "BP" state saved earlier.
      ! Purpose: Runs the MMF neural net (also called spcam)
      implicit none

      type(physics_ptend), intent(out)   :: ptend 
      type(physics_state), intent(in)    :: state
      type(physics_state), intent(in)    :: state_aphys1
      type(physics_buffer_desc), pointer :: pbuf(:)  
      type(cam_in_t),intent(in)          :: cam_in
      type(cam_out_t),     intent(inout) :: cam_out
      real(r8), intent(in)               :: coszrs(pcols)
      real(r8), intent(in)               :: solin(pcols)
      real(r8), intent(in)               :: ztodt
      integer, intent(in)                :: lchnk

      ! Local Variables
      real(r8), allocatable :: input(:, :)
      real(r8), allocatable :: output(:,:)
      integer, parameter    :: in_dims = 2
      integer, allocatable  :: in_shape(:)
      integer  :: in_layout(in_dims) = [1, 2]
      integer, allocatable  :: out_shape(:)
      integer  :: out_layout(in_dims) = [1, 2]
      
      integer   :: i, k, ncol, ixcldice, ixcldliq
      real (r8) :: s_bctend(pcols,pver), q_bctend(pcols,pver), qc_bctend(pcols,pver), qi_bctend(pcols,pver), qafter, safter, &
      u_bctend(pcols,pver), v_bctend(pcols,pver)
      logical   :: do_constraints
      logical   :: lq(pcnst)
      real(r8)  :: rh_loc

      ! NN variables
      ! real(real32) :: input_torch(inputlength, pcols)
      type(torch_tensor), dimension(1) :: in_tensor, out_tensor !
      real(real32), pointer :: input_torch(:, :) ! input_torch(inputlength, pcols)
      real(real32), pointer :: output_torch(:, :)

      real(r8), pointer, dimension(:)   :: lhflx, shflx, taux, tauy ! (/pcols/)
      real(r8), pointer, dimension(:,:) :: ozone, ch4, n2o ! (/pcols,pver/)

      real(r8) :: math_pi

      math_pi = 3.14159265358979323846_r8
      
      ncol = state%ncol
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)
      lq(:)        = .FALSE.
      lq(1)        = .TRUE. ! water vapor
      lq(ixcldliq) = .TRUE. ! cloud liquid
      lq(ixcldice) = .TRUE. ! cloud ice
      call physics_ptend_init(ptend, state%psetcols, 'neural-net', & ! Initialize local physics_ptend object
                              ls=.true., &
                              lq=lq,     &
                              lu=.true., &
                              lv=.true.  &
                            )
      
      do_constraints = .true.

      ! Initialize NN output variables 
      s_bctend(:,:)  = 0.
      q_bctend(:,:)  = 0.
      qc_bctend(:,:) = 0.
      qi_bctend(:,:) = 0.
      u_bctend(:,:)  = 0.
      v_bctend(:,:)  = 0.

      ! Associate pointers with pbuf fields
      call pbuf_get_field(pbuf, pbuf_get_index('LHFLX'), lhflx)
      call pbuf_get_field(pbuf, pbuf_get_index('SHFLX'), shflx)
      call pbuf_get_field(pbuf, pbuf_get_index('TAUX'),  taux)
      call pbuf_get_field(pbuf, pbuf_get_index('TAUY'),  tauy)
      call pbuf_get_field(pbuf, pbuf_get_index('ozone'), ozone)
      call pbuf_get_field(pbuf, pbuf_get_index('CH4'),   ch4)
      call pbuf_get_field(pbuf, pbuf_get_index('N2O'),   n2o)

      ! Set the input array. This is based on v4
      ! 1200 seconds is the timestep
      input(:ncol,0*pver+1:1*pver) = state%t(1:ncol,1:pver)          ! state_t
      input(:ncol,1*pver+1:2*pver) = state%q(1:ncol,1:pver,1)        ! state_q0001
      input(:ncol,2*pver+1:3*pver) = state%q(1:ncol,1:pver,ixcldliq)
      input(:ncol,3*pver+1:4*pver) = state%q(1:ncol,1:pver,ixcldice)
      input(:ncol,4*pver+1:5*pver) = state%u(1:ncol,1:pver)          ! state_u
      input(:ncol,5*pver+1:6*pver) = state%v(1:ncol,1:pver)          ! state_v
      ! state dynamics tendencies (large-scale forcings)
      input(:ncol,6*pver+1:7*pver) = (state%t(1:ncol,1:pver)-state_aphys1%t(1:ncol,1:pver))/1200. ! state_t_dyn
      ! state_q0_dyn
      input(:ncol,7*pver+1:8*pver) = (state%q(1:ncol,1:pver,1)-state_aphys1%q(1:ncol,1:pver,1) + state%q(1:ncol,1:pver,ixcldliq)-state_aphys1%q(1:ncol,1:pver,ixcldliq) + state%q(1:ncol,1:pver,ixcldice)-state_aphys1%q(1:ncol,1:pver,ixcldice))/1200.
      input(:ncol,8*pver+1:9*pver) = (state%u(1:ncol,1:pver)-state_aphys1%u(1:ncol,1:pver))/1200. ! state_u_dyn
      input(:ncol,9*pver+1:10*pver) = state%t_adv(2,1:ncol,1:pver)
      input(:ncol,10*pver+1:11*pver) = state%q_adv(2,1:ncol,1:pver,1) + state%q_adv(2,1:ncol,1:pver,ixcldliq) + state%q_adv(2,1:ncol,1:pver,ixcldice)
      input(:ncol,11*pver+1:12*pver) = state%u_adv(2,1:ncol,1:pver)
      ! previous state physics tendencies
      input(:ncol,12*pver+1:13*pver) = state%t_phy(1,1:ncol,1:pver)
      input(:ncol,13*pver+1:14*pver) = state%q_phy(1,1:ncol,1:pver,1)
      input(:ncol,14*pver+1:15*pver) = state%q_phy(1,1:ncol,1:pver,ixcldliq)
      input(:ncol,15*pver+1:16*pver) = state%q_phy(1,1:ncol,1:pver,ixcldice)
      input(:ncol,16*pver+1:17*pver) = state%u_phy(1,1:ncol,1:pver)
      ! 2-step in the past physics tendencies
      input(:ncol,17*pver+1:18*pver) = state%t_phy(2,1:ncol,1:pver)
      input(:ncol,18*pver+1:19*pver) = state%q_phy(2,1:ncol,1:pver,1)
      input(:ncol,19*pver+1:20*pver) = state%q_phy(2,1:ncol,1:pver,ixcldliq)
      input(:ncol,20*pver+1:21*pver) = state%q_phy(2,1:ncol,1:pver,ixcldice)
      input(:ncol,21*pver+1:22*pver) = state%u_phy(2,1:ncol,1:pver)
      !gas
      input(:ncol,22*pver+1:23*pver) = ozone(:ncol,1:pver)            ! pbuf_ozone
      input(:ncol,23*pver+1:24*pver) = ch4(:ncol,1:pver)             ! pbuf_CH4
      input(:ncol,24*pver+1:25*pver) = n2o(:ncol,1:pver)             ! pbuf_N2O
      ! 2d vars e.g., ps, solin
      input(:ncol,25*pver+1) = state%ps(1:ncol)                      ! state_ps
      input(:ncol,25*pver+2) = solin(1:ncol)                         ! pbuf_SOLIN
      input(:ncol,25*pver+3) = lhflx(1:ncol)                         ! pbuf_LHFLX
      input(:ncol,25*pver+4) = shflx(1:ncol)                         ! pbuf_SHFLX
      input(:ncol,25*pver+5) = taux(1:ncol)                          ! pbuf_TAUX
      input(:ncol,25*pver+6) = tauy(1:ncol)                          ! pbuf_TAUY
      input(:ncol,25*pver+7) = coszrs(1:ncol)                        ! pbuf_COSZRS
      input(:ncol,25*pver+8) = cam_in%ALDIF(:ncol)                   ! cam_in_ALDIF
      input(:ncol,25*pver+9) = cam_in%ALDIR(:ncol)                   ! cam_in_ALDIR
      input(:ncol,25*pver+10) = cam_in%ASDIF(:ncol)                  ! cam_in_ASDIF
      input(:ncol,25*pver+11) = cam_in%ASDIR(:ncol)                  ! cam_in_ASDIR
      input(:ncol,25*pver+12) = cam_in%LWUP(:ncol)                   ! cam_in_LWUP
      input(:ncol,25*pver+13) = cam_in%ICEFRAC(:ncol)                ! cam_in_ICEFRAC
      input(:ncol,25*pver+14) = cam_in%LANDFRAC(:ncol)               ! cam_in_LANDFRAC
      input(:ncol,25*pver+15) = cam_in%OCNFRAC(:ncol)                ! cam_in_OCNFRAC
      input(:ncol,25*pver+16) = cam_in%SNOWHICE(:ncol)               ! cam_in_SNOWHICE
      input(:ncol,25*pver+17) = cam_in%SNOWHLAND(:ncol)              ! cam_in_SNOWHLAND
      !5 placeholder for future input (which we plan to use some 2d vars at t-1, but in the end did not implement)
      input(:ncol,25*pver+18:25*pver+22) = 0._r8
      ! cos lat and sin lat
      do i = 1,ncol ! lat is get_lat_p(lchnk,i), 23/24 needs cos/sin
        input(i,25*pver+23) = cos(get_rlat_p(lchnk,i))
        input(i,25*pver+24) = sin(get_rlat_p(lchnk,i))
      end do
      input(:ncol,25*pver+25) = 0._r8               ! icol ! can be 1-384 in future but not implemented yet
      ! RH conversion
      if (input_rh) then ! relative humidity conversion for input
        do i = 1,ncol
          do k=1,pver
            ! Port of tom's RH =  Rv*p*qv/(R*esat(T))
            rh_loc = 461.*state%pmid(i,k)*state%q(i,k,1)/(287.*tom_esat(state%t(i,k))) ! note function tom_esat below refercing SAM's sat.F90
            input(i,1*pver+k) = rh_loc
          end do
        end do
      end if


      ! FTorch Inference. Start with switching the indices
      input_torch(:,:) = 0.
      do i=1,ncol
        do k=1,inputlength
          input_torch(k,i) = input(i,k)
        end do
      end do
      ! Set the input tensor 
      call torch_tensor_from_array(in_tensor(1), input_torch,  in_layout, torch_kCPU) ! Ftorch
      ! Set the out_tensor before moving forward one time step
      call torch_tensor_from_array(out_tensor(1), output_torch, out_layout, torch_kCPU) 
      if (masterproc) write(iulog,*) 'Calling Ftorch!' !ftorch
      print *, "Calling FTorch!"
      call torch_model_forward(model_pytorch, in_tensor, out_tensor) ! Ftorch
      ! End with switching the indices back
      do i=1, ncol
        do k=1,outputlength
          output(i,k) = output_torch(k,i)
        end do
      end do

      !!! do post processing some non-negative constraints + option to do stratosphere water denial

      ! Manually applying ReLU activation for positive-definite variables
      ! [TODO] for ensemble, ReLU should be moved before ens-averaging
      do i=1,ncol
        ! ReLU for the last 8 variables 
        do k=outputlength-7,outputlength
          output(i,k) = max(output(i,k), 0.)
        end do
        ! ! tiny flwds
        ! k=6*pver+2
        ! output(i,k) = max(output(i,k), tiny(output(i,k))) ! flwds
        !                                                   ! preventing flwds==0 error
  
        ! zero out surface solar fluxes when local time is at night
        if (coszrs(i) .le. 0.) then
          output(i,6*pver+1) = 0. ! netsw
          output(i,6*pver+5) = 0. ! sols
          output(i,6*pver+6) = 0. ! soll
          output(i,6*pver+7) = 0. ! solsd
          output(i,6*pver+8) = 0. ! solld
        endif
      end do

      ! ---------- 1. NN output to atmosphere forcing --------
      ! ['TBCTEND', 'QBCTEND','CLDLIQBCTEND','CLDICEBCTEND']
      s_bctend (:ncol,1:pver) = output(1:ncol,0*pver+1:1*pver)*cpair ! K/s --> J/kg/s (ptend expects that)
      q_bctend (:ncol,1:pver) = output(1:ncol,1*pver+1:2*pver)       ! kg/kg/s
      qc_bctend(:ncol,1:pver) = output(1:ncol,2*pver+1:3*pver)       ! kg/kg/s
      qi_bctend(:ncol,1:pver) = output(1:ncol,3*pver+1:4*pver)       ! kg/kg/s
      u_bctend (:ncol,1:pver) = output(1:ncol,4*pver+1:5*pver)       ! m/s/s
      v_bctend (:ncol,1:pver) = output(1:ncol,5*pver+1:6*pver)       ! m/s/s

      ! ! Remove any moisture activity and remove all clouds in the stratosphere:
      ! If true, zero out cloud (qc and qi) above tropopause in the NN output
      ! if (cb_strato_water_constraint) then
      !   do i=1,ncol
      !     call detect_tropopause(state%t(i,:),state%exner(i,:),state%zm(i,:),state%pmid(i,:),idx_trop(i))
      !     q_bctend (i,1:idx_trop(i)) = 0.
      !     qc_bctend(i,1:idx_trop(i)) = -state%q(i,1:idx_trop(i),ixcldliq)/ztodt
      !     qi_bctend(i,1:idx_trop(i)) = -state%q(i,1:idx_trop(i),ixcldice)/ztodt
      !   end do
      !   call outfld('TROP_IND', idx_trop(:ncol)*1._r8, ncol, state%lchnk)
      ! end if

      ! -- atmos positivity constraints ---- 
      if (do_constraints) then
        do i=1,ncol
          do k=1,pver
            ! energy positivity:
            safter = state%s(i,k) + s_bctend(i,k)*ztodt ! predicted DSE after NN tendency
            if (safter .lt. 0.) then ! can only happen when bctend < 0...
              s_bctend(i,k) = s_bctend(i,k) + abs(safter)/ztodt ! in which case reduce cooling rate
              write (iulog,*) 'HEY MMF_NN_EMULATOR made a negative absolute temperature, corrected but BEWARE!!!'
              write (iulog,*) '' ! [TODO] printout lat/lon and error magnitude
            endif
     
            ! vapor positivity:
            qafter = state%q(i,k,1) + q_bctend(i,k)*ztodt ! predicted vapor after NN tendency
            if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
              q_bctend(i,k) = q_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
              write (iulog,*) 'HEY MMF_NN_EMULATOR made a negative absolute q, corrected but BEWARE!!!'
            endif
     
            ! liquid positivity:
            qafter = state%q(i,k,ixcldliq) + qc_bctend(i,k)*ztodt ! predicted liquid after NN tendency
            if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
              qc_bctend(i,k) = qc_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
              !write (iulog,*) 'HEY MMF_NN_EMULATOR made a negative absolute qc, corrected but BEWARE!!!'
            endif

            ! ice positivity:
            qafter = state%q(i,k,ixcldice) + qi_bctend(i,k)*ztodt ! predicted ice after NN tendency
            if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
              qi_bctend(i,k) = qi_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
              !write (iulog,*) 'HEY MMF_NN_EMULATOR made a negative absolute qi, corrected but BEWARE!!!'            
            endif
          end do
        end do
      endif

      ! Wire to ptend:
      ptend%s(:ncol,:pver)          = s_bctend(:ncol,:pver)
      ptend%q(:ncol,:pver,1)        = q_bctend(:ncol,:pver)
      ptend%q(:ncol,:pver,ixcldliq) = qc_bctend(:ncol,:pver)
      ptend%q(:ncol,:pver,ixcldice) = qi_bctend(:ncol,:pver)
      ptend%u(:ncol,:pver)          = u_bctend(:ncol,:pver)
      ptend%v(:ncol,:pver)          = v_bctend(:ncol,:pver)
 
      ! for q{1,2,3}, u, v tendencies, top levels [1 to 12] are not coupled
      ! no std dev except ptend%q1[1:3], whose magnitudes are smaller that other levels by orders of magnitude
      if (cb_top_levels_zero_out) then
        ptend%q(:ncol,1:cb_n_levels_zero,1)        = 0. 
        ptend%q(:ncol,1:cb_n_levels_zero,ixcldliq) = 0. 
        ptend%q(:ncol,1:cb_n_levels_zero,ixcldice) = 0. 
        ptend%u(:ncol,1:cb_n_levels_zero)          = 0. 
        ptend%v(:ncol,1:cb_n_levels_zero)          = 0. 
      end if
    end subroutine neural_net
#endif /* MMF_NN_EMULATOR */

    real(r8) function tom_esat(T) 
      ! For consistency with the python port of Tom's RH-calculator, this is how it
      ! was done in the training environment (Caution: could be porting bugs here)
      implicit none
      real(r8) T
      real(r8), parameter :: T0 = 273.16
      real(r8), parameter :: T00 = 253.16
      real(r8), external :: esatw_crm, esati_crm ! register functions from crm source.
      real(r8) :: omtmp,omega
      omtmp = (T-T00)/(T0-T00)
      omega = max(0.,min(1.,omtmp))
      !tf.where(T>T0,eliq(T),tf.where(T<T00,eice(T),(omega*eliq(T)+(1-omega)*eice(T))))
      if (T .gt. T0) then
        tom_esat = tom_eliq(T)
      elseif (T .lt. T00) then
        tom_esat = tom_eice(T)
      else
        tom_esat = omega*tom_eliq(T) + (1.-omega)*tom_eice(T)
      endif
    end

    real(r8) function tom_eliq(T)
      implicit none
      real(r8) T
      real(r8), parameter :: T0 = 273.16
      real(r8), parameter :: cliq = -80. 
      real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
      real(r8) :: dt
      dt = max(cliq,T-T0)
      tom_eliq = 100.*(a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))))  
    end 


    real(r8) function tom_eice(T)
      implicit none
      real(r8) T
      real(r8), parameter :: T0 = 273.16
      real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
          6.11147274, 0.503160820, 0.188439774e-1, &
          0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
          0.385852041e-9, 0.146898966e-11, 0.252751365e-14/       
      real(r8) cice(6)
      real(r8) dt
      dt = T-T0
      cice(1) = 273.15
      cice(2) = 185.
      cice(3) = -100.
      cice(4) = 0.00763685
      cice(5) = 0.000151069
      cice(6) = 7.48215e-07
      if (T .gt. cice(1)) then
        tom_eice = tom_eliq(T)
      else if (T .le. cice(2)) then
        tom_eice = 100.*(cice(4) + max(cice(3),dt)*(cice(5)+max(cice(3),dt)*cice(6))) 
      else
        tom_eice = 100.*(a0 +dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))))
      end if
    end

    subroutine mmf_nn_emulator_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr, f
      character(len=*), parameter :: subname = 'mmf_nn_emulator_readnl'
      
      namelist /mmf_nn_emulator_nl/ inputlength, outputlength, input_rh, &
                          !  cb_partial_coupling, cb_partial_coupling_vars,&
                          !  cb_use_input_prectm1, &
                          !  cb_nn_var_combo, &
                           cb_torch_model 
                          !  cb_spinup_step, &
                          !  cb_do_ramp, cb_ramp_linear_steps, &
                          !  cb_ramp_option, cb_ramp_factor, cb_ramp_step_0steps, cb_ramp_step_1steps, &
                          !  cb_strato_water_constraint, dtheta_thred

      ! Initialize 'cb_partial_coupling_vars'
      ! do f = 1, pflds
      !   cb_partial_coupling_vars(f) = ' '
      ! end do

      ierr = 0
      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'mmf_nn_emulator_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, mmf_nn_emulator_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

      ! Broadcast namelist variables
      call mpibcast(inputlength,  1,                 mpiint,  0, mpicom)
      call mpibcast(outputlength, 1,                 mpiint,  0, mpicom)
      call mpibcast(input_rh,     1,                 mpilog,  0, mpicom)
      call mpibcast(cb_torch_model, len(cb_torch_model), mpichar, 0, mpicom)
   end subroutine mmf_nn_emulator_readnl

end module mmf_nn_emulator