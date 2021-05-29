module io

    ! modules
    use constants

    implicit none
    ! static
    ! dynamic
    integer :: ios

    contains ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine open_file_w(unit, path) ! -------------------------------------> SUBR open_file_w
        
        character(len=*), intent(in) :: path                                    ! path to file
        integer, intent(in) :: unit
        
        open(                                                               &
            unit=unit, access='sequential', action='write',                 &
            file=path, iostat=ios, status='replace'                         &
        )

    end subroutine open_file_w


    subroutine close_file(unit) ! --------------------------------------------> SUBR close_file()

        integer, intent(in) :: unit

        close(unit=unit, iostat=ios, status='keep')

    end subroutine close_file


end module io