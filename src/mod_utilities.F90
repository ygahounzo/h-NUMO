!-------------------------------
!> @brief This utility module is based on the DART utility module
!-------------------------------
module mod_utilities   

    use mod_mpi_utilities

    implicit none
    private

    integer, parameter :: E_DBG  = -10, E_MSG = 0 , E_WARN  = 10, E_ERR  = 20

    integer, parameter :: TERMLEVEL = 11

    public :: find_namelist
    public :: E_DBG, E_MSG, E_WARN, E_ERR, &
        get_unit, error_handler,     &
        do_output
    
    !> Change string case
    public :: lowercase
    public :: uppercase

contains

    !------------------------------------------------------------------------------
    !>@brief This function gets the next available unit for i/o
    !------------------------------------------------------------------------------
    function get_unit() result(iunit)
        integer :: i, iunit
        logical :: available

        iunit = -1
        do i=10,100
            inquire(i,opened=available)
            if(.not.available) then
                iunit = i
                exit
            end if
        end do

        if(iunit == -1) call error_handler(E_ERR,'get_unit','No available units')

    end function


    !------------------------------------------------------------------------------
    !>@brief This subroutine is a basic error handler
    !------------------------------------------------------------------------------
    subroutine error_handler(level_in, routine, message)
        integer, intent(in), optional :: level_in
        character(len=*), intent(in) :: routine
        character(len=*), intent(in), optional :: message

        integer :: level

        if(present(level_in)) then
            level = level_in
        else
            level = E_ERR
        end if

        if(irank==irank0) then
            select case (level)
                case(E_MSG)
                    write(*,*) 'MESSAGE FROM:'
                case(E_WARN)
                    write(*,*) 'WARNING FROM:'
                case(E_ERR)
                    write(*,*) 'ERROR FROM:'
                case(E_DBG)
                    write(*,*) 'DEBUG FROM:'
            end select

            write(*,'(4x,A,1x,A)') 'routine:',trim(routine)
            if(present(message)) then
                write(*,'(4x,A,1x,A)') 'message:',trim(message)
            end if
        end if

        if(level >= TERMLEVEL) call exit_all(-99)

    end subroutine error_handler


    !------------------------------------------------------------------------------
    !>@brief Function returns true if output is writen to stdout
    !------------------------------------------------------------------------------
    function do_output()
        logical :: do_output

        if(irank==irank0) then
            do_output = .true.
        else
            do_output = .false.
        end if

    end function do_output

    !> lowercase
    !> ---------------
    !>@brief Wraper to change to case of a string to lower case
    !>  PARAMETERS
    !>   IN  str  Character string to change case of
    function lowercase(str) result(lstr)
        character(len=*), intent(in)        :: str
        character(len=len(str))             :: lstr
        logical                             :: TO_LOWER = .true.
        lstr = str
        call change_case(lstr,TO_LOWER)
    end function lowercase
  
    !> uppercase
    !> ---------------
    !>@brief Wraper to change to case of a string to upper case
    !> Sets a string to upper case
    !>  PARAMETERS
    !>   IN  str  Character string to change case of
    function uppercase(str) result(ustr)
        character(len=*), intent(in)       :: str
        character(len=len(str))            :: ustr
        logical                            :: TO_LOWER = .false.
        ustr = str
        call change_case(ustr,TO_LOWER)
    end function uppercase

    !> change_case
    !> ---------------
    !>@brief Changes the case of a string.
    !>  PARAMETERS
    !>   IN  string     Character string to change case of
    !>   IN  TO_LOWER   Change to lowercase or uppercase.
    subroutine change_case(string,TO_LOWER)
        character(len=*), intent(inout)  :: string
        logical, intent(in)              :: TO_LOWER

        integer                          :: idiff, i1, i
        integer                          :: ia, iz

        if(TO_LOWER) then
            idiff = ichar('a') - ichar('A')
            ia = ichar('A') ; iz = ichar('Z')
        else
            idiff = ichar('A') - ichar('a')
            ia = ichar('a') ; iz = ichar('z')
        end if

        do i=1,len(string)
            i1 = ichar(string(i:i))
            if( i1.ge.ia .and. i1.le.iz) then
                string(i:i) = char(ichar(string(i:i)) + idiff)
            end if
        end do
    end subroutine change_case

    !------------------------------------------------------------------------------
    !>@brief Subroutine determines if a namelist exists in a given file
    !> It it exists then open and return the unit number
    !------------------------------------------------------------------------------
    function find_namelist(nml_file_name, nml_name, iunit) result(nml_exists)
        character(len=*),  intent(in)  :: nml_file_name
        character(len=*),  intent(in)  :: nml_name
        integer,           intent(out) :: iunit
        logical :: nml_exists

        character(len = 169) :: nml_string, test_string, message
        character(len=*), parameter :: routine = 'find_namelist'
        logical :: is_opened
        integer :: io

        nml_exists = .false.

        if(file_exist(trim(nml_file_name))) then

            inquire(file=trim(nml_file_name), opened=is_opened, number=iunit)
   
            if(.not.is_opened) then
                iunit = get_unit()
                open(iunit, file = trim(nml_file_name), form='formatted', &
                    delim='apostrophe', action='read')
            end if

            test_string = '&' // trim(adjustl(nml_name))
            test_string=uppercase(test_string)

            !--------------------------------------------------------------------------
            ! Read each line until end of file is found
            !--------------------------------------------------------------------------
            do
                read(iunit, '(A)', iostat = io) nml_string
                if(io /= 0) then
                    write(message,'(A,1x,A)') trim(nml_name),'not found. Using default values.'
                    call error_handler(E_WARN, routine, message)
                    return
                else
                    if(trim(adjustl(uppercase(nml_string))) == trim(adjustl(test_string))) then
                        rewind(iunit)
                        nml_exists = .true.
                        return
                    endif
                endif
            end do
        else
            write(message,'(3(A,1x))') 'namelist file',trim(nml_file_name),'not found.'
            call error_handler(E_ERR, routine, message)
        endif

    end function find_namelist

    !------------------------------------------------------------------------------
    !>@brief This tests if a file exists
    !------------------------------------------------------------------------------
    function file_exist (file_name)
        character(len=*), intent(in) :: file_name
        logical  file_exist
        inquire (file=file_name(1:len_trim(file_name)), exist=file_exist)
    end function file_exist


    

end module mod_utilities
