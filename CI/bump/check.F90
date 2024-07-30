program code_verification
    implicit none

    integer :: nlayers, Nfield, nl, ifield
    integer :: index_layer, index_field
    character(len=256) :: output_file, ref_line, ref_field, ref_max, ref_min
    character(len=256) :: line, field_ref
    real :: field_max_ref, field_min_ref, field_max, field_min, error_field_max, error_field_min
    character(len=256) :: ref_filename, fin_filename
    integer :: ios
    character(len=256) :: ref_file_content(1000), fin_file_content(1000)
    integer :: ref_line_count, fin_line_count
    character(len=20) :: str
    real :: mass_loss_ref, mass_loss_fin

    ! Initialize variables
    nlayers = 2
    Nfield = 4
    output_file = 'output.log'
    ref_filename = 'mlswe_ref_FIN.txt'
    fin_filename = 'mlswe_FIN.txt'
    
    ! Read the content of the reference file
    call read_file(ref_filename, ref_file_content, ref_line_count)

    ! Read the content of the fin file
    call read_file(fin_filename, fin_file_content, fin_line_count)

    ! Open the output file for writing
    open(unit=10, file=output_file, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening output file"
        stop
    endif

    write(str, '(I2)') nlayers
    write(10, '(A)') ' '
    write(10, '(A)') 'Code verification for ' // trim(adjustl(str)) // ' layers model'
    write(10, '(A)') ' '

    do nl = 1, nlayers
        index_layer = 6 * (nl - 1) + 1

        ! Read the specific line from the reference file
        call get_line(ref_file_content, index_layer, ref_line)
        write(10, '(A)') trim(adjustl(ref_line))

        ! Extract data from the reference file
        !call parse_line_mass(ref_file_content(index_layer+1), mass_loss_ref)

        ! Print mass loss 
        call get_line(ref_file_content, index_layer+1, ref_line)
        write(10, '(A)') trim(adjustl(ref_line))

        ! Extract data from the fin file
        call parse_line_mass(fin_file_content(index_layer+1), mass_loss_fin)

        if(mass_loss_fin .gt. 1.0e-12) then 
            print*, 'Layer ', nl, 'mass_loss = ', mass_loss_fin, ' to large'

            stop
        endif 

        do ifield = 2, Nfield
            index_field = ifield + index_layer

            ! Extract data from the reference file
            call parse_line(ref_file_content(index_field), field_ref, field_max_ref, field_min_ref)

            ! Extract data from the fin file
            call parse_line(fin_file_content(index_field), field_ref, field_max, field_min)

            ! Calculate errors
            error_field_max = abs(field_max_ref - field_max)/abs(field_max_ref)
            error_field_min = abs(field_min_ref - field_min)/abs(field_min_ref)

            ! Write the results to the output file
            write(10, '(A3, A3, E10.5, E10.5)') field_ref, ' = ', error_field_max, error_field_min

        end do

        write(10, '(A)') ' '
    end do

    close(10)
end program code_verification

subroutine read_file(filename, content, line_count)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=256), intent(out) :: content(*)
    integer, intent(out) :: line_count
    integer :: ios, unit
    character(len=256) :: line

    unit = 15
    open(unit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file ", filename
        stop
    endif

    line_count = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        line_count = line_count + 1
        content(line_count) = line
    end do

    close(unit)
end subroutine read_file

subroutine get_line(content, line_number, line)
    implicit none
    character(len=256), intent(in) :: content(*)
    integer, intent(in) :: line_number
    character(len=256), intent(out) :: line

    line = content(line_number)
end subroutine get_line

subroutine parse_line(line, field, field_max, field_min)
    implicit none
    character(len=256), intent(in) :: line
    character(len=4), intent(out) :: field
    real, intent(out) :: field_max, field_min

    read(line, '(A23, 4X, E24.16, 7X, E24.16)') field, field_max, field_min

    !print*, field, field_max
end subroutine parse_line

subroutine parse_line_mass(line, field_value)
    implicit none
    character(len=256), intent(in) :: line
    real, intent(out) :: field_value
    character(len=4) :: field

    !print*, line

    read(line, '(A13,4X, E12.8)') field, field_value

    !print*, field
    !print*, field_value
end subroutine parse_line_mass
