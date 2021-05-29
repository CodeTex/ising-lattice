program main

    use constants
    use io
    use shell_prompts

    ! use omp_lib                                                                                                         ! CHANGE DUE TO -fopenmp

    implicit none
    ! parameter
    integer, parameter :: unit_out=10
    integer, parameter :: T_diff=T_max-T_min
    ! dynamic
    character(len=char_max) :: output_file
    integer :: i,j,k
    integer(kind=iter) :: t1,t2
    integer :: threads
    integer :: lattice_dim                                                                                         ! CHANGE DUE TO -fopenmp
    ! integer :: lattice_dim, lattice_dim_max                                     ! lattice dimension
    ! integer, dimension(6) :: lattice_array                         ! lattice array (loop through)                       ! CHANGE DUE TO -fopenmp
    ! integer, dimension(128,128) :: lattice                             ! 2D Ising Lattice                               ! CHANGE DUE TO -fopenmp
    real(kind=p) :: time_start, time_end, time_diff
    real(kind=p) :: T, beta                                                     ! temperature, inverse temperature
    real(kind=p), dimension(5) :: M                                             ! magnetisation


    ! M(1) = M̅
    ! M(2) = <|M|>
    ! M(3) = <|M^2|>
    ! M(4) = <|M^4|>                1  <|M^4|>
    ! M(5) = Binder parameter = 1 - - ---------
    !                               3 <|M^2|>^2


    ! ########## PROGRAM START ############

    call program_prompt(title, '-')                                             ! shell_prompts

    ! threads = omp_get_max_threads()

    ! write(*, fmt='(a,i2,a,/)') 'Using up to ', threads, ' threads.'                                                     ! CHANGE DUE TO -fopenmp

    ! lattice_array = set_lattice_array(default_lattice_array)                  ! main                                  ! CHANGE DUE TO -fopenmp

    ! allocate max dim size to use same matrix for all dims
    ! lattice_dim_max = maxval(lattice_array, 1)                                                                        ! CHANGE DUE TO -fopenmp
    ! lattice_dim_max = 128
    ! allocate(lattice(lattice_dim_max,lattice_dim_max))                                                                ! CHANGE DUE TO -fopenmp
    ! lattice_array = (/ 4, 8, 16, 32, 64, 128/)

    write (*,fmt='(a,/)') 'Commencing sweeps...'

    do i = 1, size(lattice_array)

        lattice_dim = lattice_array(i)

        if (output) then
            output_file = create_path(dirname, lattice_dim)                     ! main
            call open_file_w(unit_out, output_file)                             ! io
        end if

        ! call cpu_time(time_start)
        call system_clock(t1)
        do j = 0, T_diff*100

            ! temperature is incremented by .1 every loop
            T = T_min + j*0.01_p
            beta = 1._p / T
            M = 0._p

            if (.not. verbose .and.  (mod(j, T_diff) .eq. 0)) then
                ! if one 10th share of workload is done
                call progress(                                              &   ! shell_prompts
                    int(j/T_diff, kind=iter),  &
                    lattice_dim                                             &
                )
            endif

            ! Initialize lattice with random magnetization (-1,1)
            call lattice_randomize(lattice, lattice_dim, lattice_dim_max)       ! main

            ! --- Start sweeps ---
            !$OMP PARALLEL DO PRIVATE(K)                                                                                ! CHANGE DUE TO -fopenmp
            do k = 1, pre_sweeps
                call sweep(                                                 &   ! main
                    lattice,lattice_dim,lattice_dim_max,                    &
                    method, beta                                            &
                )
            end do
            !$OMP END PARALLEL DO                                                                                       ! CHANGE DUE TO -fopenmp

            !$OMP PARALLEL DO PRIVATE(K)                                                                                ! CHANGE DUE TO -fopenmp
            do k = 1, sweeps
                call sweep(                                                 &   ! main
                    lattice,lattice_dim,lattice_dim_max,                    &
                    method, beta                                            &
                )

                M(1) = sum(lattice(:lattice_dim,:lattice_dim)) /            &
                    lattice_dim**2._p
                M(2) = M(2) + abs(M(1)) / sweeps
                M(3) = M(3) + abs(M(1))**2._p / sweeps
                M(4) = M(4) + abs(M(1))**4._p / sweeps
            end do
            !$OMP END PARALLEL DO                                                                                       ! CHANGE DUE TO -fopenmp

            M(5) = 1._p - M(4) / (3._p * M(3)**2._p)

            ! --- End sweeps ---

            if (output) then
                write (unit_out,*)                                          &
                    T, M(2:),                                               &
                    calcE(lattice,lattice_dim,lattice_dim_max)                  ! main
            end if

            if (verbose) then
                call visualize_magnetization(lattice, lattice_dim_max)          ! shell_prompts
                write (*, fmt='(a, f14.6)')                                 &
                    'Average energy:',                                      &
                    calcE(lattice,lattice_dim,lattice_dim_max)                  ! main
                write (*,*)   &
                    '<|M|> =', M(2),' <|M^2|> =', M(3),                     &
                    ' <|M^4|> =', M(4),' BK =', M(5)
            end if

        end do
        ! call cpu_time(time_end)
        call system_clock(t2)

        ! write(*,*) time_end-time_start

        if (.not. verbose) then
            time_diff = (t2 - t1) / 1000._p
            write (*, fmt='(a,i3,a,f6.3,a)') ' ✓ in ',                      &
                int(time_diff/60, kind=4), ' min ',                         &
                mod(time_diff,60._p), ' s'
        endif

        if (output) call close_file(unit_out)

    end do


    call program_prompt('Program ended successully', '-', len(trim(title)) + 2) ! shell_prompts


    return ! ########### PROGRAM END #############


    contains ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! function set_lattice_array(use_default) ! --------------------------------> FUNC set_lattice_array()
    !     ! set array of lattice dimensions through which to parse through

    !     integer(kind=iter) :: i
    !     integer :: lattice_n, default_lattice_n                                 ! lattice array length, default length
    !     integer, dimension(:), allocatable :: lattice_array_default             ! default lattice-dim-array
    !     integer, dimension(:), allocatable :: set_lattice_array                 ! lattice-dim-array
    !     logical, intent(in) :: use_default                                      ! use default or enter via input

    !     ! set default values
    !     default_lattice_n = 6                                                   ! makes sure allocate correctly!
    !     allocate(lattice_array_default(default_lattice_n))
    !     lattice_array_default = (/ 4, 8, 16, 32, 64, 128 /)                     ! ATTENTION: 64 ~50min, 128 ~5h

    !     if (use_default) then
    !         allocate(set_lattice_array(default_lattice_n))
    !         do i = 1, default_lattice_n
    !             set_lattice_array(i) = lattice_array_default(i)
    !         end do
    !     else
    !         call user_input(lattice_array, lattice_n)
    !     endif

    ! end function set_lattice_array


    subroutine lattice_randomize(lattice, dim, dim_max) ! --------------------> SUBR lattice_randomize()
        ! populate lattice with ranomized spin orientations (-1, 1)

        integer, intent(in) :: dim, dim_max                                     ! dim_max = allocated size, dim = actual used size
        integer, intent(inout) :: lattice(dim_max,dim_max)
        integer :: i,j
        real, save :: rand

        do i = 1, dim
            do j = 1, dim
                call random_number(rand)
                if (rand .gt. 0.5_p) then
                    lattice(i,j) = 1
                else
                    lattice(i,j) = -1
                endif
            enddo
        enddo

    end subroutine lattice_randomize


    function create_path(dirname, lattice_dim) ! -----------------------------> FUNC create_path()

        character(len=char_max), intent(in) :: dirname
        character(len=char_max) :: create_path, dim_formatted
        integer, intent(in) :: lattice_dim

        write(dim_formatted, fmt='(i3.3)') lattice_dim
        create_path = trim(adjustl(dirname)) // "/ising_lattice-"           &
                        // trim(dim_formatted) // ".dat"

    end function create_path


    subroutine sweep(lattice,dim,dim_max,method,beta) ! ----------------------> SUBR sweep()
        ! perform sweep routine

        character(len=*), intent(in) :: method
        integer, intent(in) :: dim, dim_max                                     ! dim_max = allocated size, dim = actual used size
        integer, intent(inout) :: lattice(dim_max,dim_max)
        integer :: i,j
        real, save :: rand
        real(kind=p), intent(in) :: beta
        real(kind=p) :: delta_E, r
        ! real(kind=p) :: E_n, E_ij

        do i=1,dim
            do j=1,dim

                ! E_n = calcE(lattice,dim,dim_max)
                ! E_ij = calc_dE(i,j,lattice,dim,dim_max)
                !
                ! delta_E = E_n - E_ij

                ! get energy delta
                delta_E = calc_dE(i,j,lattice,dim,dim_max)        ! main
                ! transition probability
                r = exp(-beta * delta_E)
                ! optional spin-flip
                call random_number(rand)
                if (rand .lt. acceptance_probability(method,r)) then
                    lattice(i,j) = -lattice(i,j)
                endif

            end do
        end do

    end subroutine sweep

    function acceptance_probability(method, r) ! -----------------------------> FUNC acceptance_probability()
        ! calculate acceptance probibilty based on passed method

        character(len=*) :: method
        real(kind=p) :: acceptance_probability, r

        select case (method)
            case ('metropolis')
                acceptance_probability =  min(1._p, r)
            case ('heat_bath')
                acceptance_probability = r / (1 + r)
            case default
                acceptance_probability =  min(1._p, r)
        end select

    end function acceptance_probability


    function calcE(lattice, dim, dim_max) ! ----------------------------------> FUNC calcE()
        ! calculate energy of current state

        integer, intent(in) :: dim, dim_max
        integer, intent(in) :: lattice(dim_max,dim_max)
        integer :: i,j,k
        integer :: Inew, Jnew
        real(kind=p) :: calcE, energy

        ! Determine the initial energy
        energy = 0.0

        do i = 1, dim
            do j = 1, dim

                ! Loop over nearest neighbors
                do k = 1, NN
                    Inew = i + Inn(k)
                    Jnew = j + Jnn(k)

                    ! periodic boundary conditions
                    if (Inew .le. 0) then
                            Inew = dim
                    else if(Inew .gt. dim) then
                        Inew = 1
                    endif

                    ! periodic boundary conditions
                    if (Jnew .le. 0) then
                        Jnew = dim
                    else if(Jnew .gt. dim) then
                        Jnew = 1
                    endif

                    ! Update the energy
                    energy = energy - JJ * lattice(i,j) * lattice(Inew,Jnew)
                    enddo
                ! Calculate the contribution from the field H
                energy = energy - 2.d0 * HH * lattice(i,j)

            enddo
        enddo

        ! Account for double counting
        energy = energy / 2._p
        calcE = energy

    end function calcE


    function calc_dE(i, j, lattice, dim, dim_max) ! --------------------------> FUNC calc_dE()
        ! calculate energy delta between previous and flipped state

        integer, intent(in) :: dim, dim_max
        integer, intent(in) :: lattice(dim_max,dim_max)
        integer :: i,j,k
        integer :: Inew, Jnew
        real(kind=p) :: calc_dE, energy_flipped, energy_prev

        energy_flipped = 0._p
        energy_prev = 0._p

        ! Loop over nearest neighbors
        do k = 1, NN
            Inew = i + Inn(k)
            Jnew = j + Jnn(k)

            ! periodic boundary conditions
            if (Inew .le. 0) then
                    Inew = dim
            else if (Inew .gt. dim) then
                Inew = 1
            endif

            ! periodic boundary conditions
            if (Jnew .le. 0) then
                Jnew = dim
            else if (Jnew .gt. dim) then
                Jnew = 1
            endif

            ! Calculate delta E, flip lattice(i,j) on site
            energy_prev =                                                   &
                energy_prev - JJ * lattice(i,j) * lattice(Inew,Jnew)
            energy_flipped =                                                &
                energy_flipped + JJ * lattice(i,j) * lattice(Inew,Jnew)
        enddo

        ! -> H part cancels for delta E calculation
        ! energy_prev = energy_prev - 2.d0 * HH * lattice(i,j)
        ! energy_flipped = energy_flipped - 2.d0 * HH * lattice(i,j)
        calc_dE = energy_flipped - energy_prev

    end function calc_dE

end program main
