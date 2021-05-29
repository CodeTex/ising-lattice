module shell_prompts

    use constants

    implicit none

    contains ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine program_prompt(text, div_sym, div_len) ! ----------------------> SUBR program_prompt()
        
        character(len=*), intent(in) :: text                                    ! text to be displayed
        character(len=1), intent(in) :: div_sym                                 ! symbol of the divider
        integer, optional :: div_len
        integer :: div_length
        integer :: text_shift

        if (present(div_len)) then
            div_length = div_len
        else
            div_length = len(trim(text)) + 2
        endif

        text_shift = (div_length - len(trim(text))) / 2 

        write(*, fmt='(/, a)')          repeat(div_sym, div_length)
        write(*, fmt='(a, a)')          repeat(' ', text_shift), trim(text)     ! align center horizontally
        write(*, fmt='(a, /)')          repeat(div_sym, div_length)

    end subroutine program_prompt


    subroutine progress(pos, lattice_dim) ! ----------------------------------> ! SUBR progress()

        character(len=28) :: bar="dddxddd -> ???% |          |"
        integer(kind=iter) :: i
        integer(kind=iter), intent(in) :: pos, lattice_dim

        write(unit=bar(1:3), fmt='(i3.3)') lattice_dim
        write(unit=bar(5:7), fmt='(i3.3)') lattice_dim
        write(unit=bar(12:14), fmt='(i3)') 10 * pos
        
        if (pos .gt. 0) then
            do i = 1, pos
                bar(17+i:17+i) = '*'
            enddo
        endif
        
        write(unit=*, fmt='(a1,a31,$)') char(13), bar

         if (pos .eq. 10) bar = "dddxddd -> ???% |          |"

    end subroutine progress


    subroutine user_input(array, length)

        integer(kind=iter) :: i
        integer, intent(inout) :: length
        integer, dimension(:), allocatable, intent(inout) :: array

        write (*, fmt='(a)', advance='no')                                  &   
            "Enter N° of lattices: "
        read *, length

        ! add do loop
        ! if ((length .lt. 1) .or. (length .gt. 10)) then
        !     write(*,*) "Lattice number must be between 1 and 10"
        !     stop
        ! endif

        allocate(array(length))
        
        write (*,*)
        do i = 1, length

            write (*, fmt='(a,i2.2,a)', advance='no')                       &
                " Lattice #", i, " dim:  "
            read (*,*) array(i)

            if (array(i) .ge. 1000) then
                write(*,fmt='(/,a,/)') 'ERR - max lattice dim: 999'
                stop
            endif

        end do
        write (*,*)

    end subroutine user_input


    subroutine visualize_magnetization(lattice, allocated_size) ! ------------> SUBR visualize_magnetization()
        ! display magnetization of lattice to shell with '->' and '<-' symbols
    
        integer, intent(in) :: allocated_size
        integer, intent(inout) :: lattice(allocated_size,allocated_size)
        integer(kind=iter) :: i,j
        integer(kind=4) :: lattice_dim

        lattice_dim = size(lattice, 1)

        write (*, fmt='(a, a, a)') '/', repeat('-', lattice_dim * 3), '\'

        do i = 1, lattice_dim
            write (*, fmt='(a)', advance='no') '|'
            do j = 1, lattice_dim
                if (lattice(i,j) .eq. 1) then
                    write (*, fmt='(a, 1x)', advance='no') '->'
                else
                    write (*, fmt='(a, 1x)', advance='no') '<-'
                endif
            enddo
            write (*, '(a)') '|'
        enddo

        write (*, fmt='(a, a, a)') '\', repeat('-', lattice_dim * 3), '/'

    end subroutine visualize_magnetization


    

end module shell_prompts