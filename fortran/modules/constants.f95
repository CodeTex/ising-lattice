module constants

    implicit none
    ! static
    integer, parameter :: iter=selected_int_kind(8)                             ! kind for iterative variables
    integer, parameter :: p=selected_real_kind(16)                              ! float precision parameter
    integer, parameter :: char_max=80                                           ! max character length
    ! Program Interface Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    character(len=char_max) :: title='Monte Carlo Simulation of Ising Model'    ! title which is displayed in shell
    character(len=char_max), parameter :: dirname='data'                        ! path for data output
    logical, parameter :: default_lattice_array=.true.                          ! use shell prompt to enter lattice sizes
    logical, parameter :: output=.true.                                         ! generate output files
    logical, parameter :: verbose=.false.                                       ! console ouput with additional information
    ! Ising Model Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    character(len=char_max), parameter :: method='heat_up'                      ! 'metropolis' or 'heat_bath'
    integer, parameter :: T_min=1                                               ! min. temperature [J/k_B]
    integer, parameter :: T_max=4                                               ! max. temperature [J/k_B]
    integer, parameter :: T_steps=40
    integer, parameter :: sweeps=50000                                          ! N° of sweeps
    integer, parameter :: pre_sweeps=int(sweeps*0.15_p ,kind=iter)              ! 15% warm-up sweeps
    integer, parameter :: NN=4                                                  ! N° of neighbors                                                                   
    integer, dimension(NN), parameter :: Inn=(/ 1, -1, 0, 0/)                   ! Nearest neighbor array I
    integer, dimension(NN), parameter :: Jnn=(/ 0, 0, 1, -1/)                   ! Nearest neighbor array J
    real(kind=p), parameter :: JJ=1._p, HH=0._p                                 ! Coupling constant, magnetic field

    integer, dimension(6), parameter :: lattice_array=(/4, 8, 16, 32, 64, 128/)
    integer, parameter :: lattice_dim_max=128
    integer, dimension(lattice_dim_max, lattice_dim_max) :: lattice

end module constants