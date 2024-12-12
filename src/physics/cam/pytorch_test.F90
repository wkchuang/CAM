module pytorch_test

	use shr_kind_mod, only: r8 => shr_kind_r8

  use ftorch, only : torch_model, torch_model_load, torch_model_forward, & 
	                   torch_tensor, torch_tensor_from_array, torch_kCPU, &
										 torch_delete
  use iso_fortran_env

  implicit none
  save

  character(len=256) :: cb_torch_model = "/glade/u/home/wchuang/spcam/notebooks/constant_model.pt"
  type(torch_model) :: model_pytorch

  public init_neural_net, neural_net
  contains

  subroutine init_neural_net()
    implicit none
    call torch_model_load(model_pytorch, trim(cb_torch_model))
  end subroutine init_neural_net

  subroutine neural_net()
    
		implicit none

		integer              :: in_layout(1) = [1]
		integer              :: out_layout(1) = [1]

		! NN variables
		type(torch_tensor), dimension(1) :: in_tensor, out_tensor
		real(real32) :: in_data(1) = 1.0
		real(real32), dimension(1) :: out_data

    call torch_tensor_from_array(in_tensor(1), in_data, in_layout, torch_kCPU) ! Ftorch
    call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU) 
    call torch_model_forward(model_pytorch, in_tensor, out_tensor)

  end subroutine neural_net

end module pytorch_test