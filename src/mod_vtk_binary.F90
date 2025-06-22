!--------------------------------------------------------------------------------------
!>@brief This module is used to write VTK files in binary format.
!>@details It extracted from the larger library LIBVTKIO by Stefano Zaghi, CNR, Italy
!> https://sites.google.com/site/stefanozaghi/lib_vtk_io 
!>
!> This module has been modified to adjust to NUMA; all the
!> unnecessary parts of the original library were erased.
!>
!>@author Simone Marras, May 2013
!>@todo I still need to clean it up a little (SM).
!>
!--------------------------------------------------------------------------------------
module mod_vtk_binary

  use mod_types, only: r8,r16

  implicit none
  private
  ! functions for VTK LEGACY
  public:: vtk_ini
  public:: vtk_geo
  public:: vtk_geo_UNST_R8
  public:: vtk_con
  public:: vtk_dat
  public:: vtk_var
  public:: vtk_var_scal_R8
  public:: vtk_var_vect_R8
  public:: vtk_var_vect2D_R8
  public:: vtk_end
  
  ! portable kind-precision
  public:: R16P, FR16P
  public:: R8P,  FR8P
  public:: R4P,  FR4P
  public:: R_P,  FR_P
  public:: I8P,  FI8P
  public:: I4P,  FI4P
  public:: I2P,  FI2P
  public:: I1P,  FI1P
  public:: I_P,  FI_P

  ! overloading of vtk_geo
  interface vtk_geo
     !  module procedure vtk_geo_unst_R8, & ! real(R8P) UNSTRUCTURED_GRID
     !                   vtk_geo_unst_R4, & ! real(R4P) UNSTRUCTURED_GRID
     !                   vtk_geo_STRP_R8, & ! real(R8P) STRUCTURED_POINTS
     !                   vtk_geo_STRP_R4, & ! real(R4P) STRUCTURED_POINTS
     !                   vtk_geo_STRG_R8, & ! real(R8P) STRUCTURED_GRID
     !                   vtk_geo_STRG_R4, & ! real(R4P) STRUCTURED_GRID
     !                   vtk_geo_RECT_R8, & ! real(R8P) RECTILINEAR_GRID
     !                   vtk_geo_RECT_R4    ! real(R4P) RECTILINEAR_GRID
     !
     !module procedure vtk_geo_unst_R4, & ! real(R4P) UNSTRUCTURED_GRID
     !     vtk_geo_STRP_R8, & ! real(R8P) STRUCTURED_POINTS
     !     vtk_geo_STRP_R4, & ! real(R4P) STRUCTURED_POINTS
     !     vtk_geo_STRG_R8, & ! real(R8P) STRUCTURED_GRID
     !     vtk_geo_STRG_R4, & ! real(R4P) STRUCTURED_GRID
     !     vtk_geo_RECT_R8, & ! real(R8P) RECTILINEAR_GRID
     !     vtk_geo_RECT_R4    ! real(R4P) RECTILINEAR_GRID

  endinterface vtk_geo

  ! overloading of vtk_var
  interface vtk_var
     !module procedure vtk_var_scal_R8, & ! real(R8P)    scalar
     !     vtk_var_scal_R4, & ! real(R4P)    scalar
     !     vtk_var_scal_I4, & ! integer(I4P) scalar
     !     vtk_var_VECT_R8, & ! real(R8P)    vectorial
     !     vtk_var_VECT_R4, & ! real(R4P)    vectorial
     !     vtk_var_VECT_I4, & ! integer(I4P) vectorial
     !     vtk_var_TEXT_R8, & ! real(R8P)    vectorial (texture)
     !     vtk_var_TEXT_R4    ! real(R4P)    vectorial (texture)vtk_var
  endinterface vtk_var
  

  !-----------------------------------------------------------------------------------------------
  !
  !Real precision definitions:
  !
  !>  integer, parameter:: R16P = selected_real_kind(33,4931) 
  !>  33  digits, range $[\pm 10^{-4931}  ,\pm 10^{+4931}   -1]$
  !>  integer, parameter:: R8P  = selected_real_kind(15,307)  
  !>  15  digits, range $[\pm 10^{-307}~~ ,\pm 10^{+307}~~  -1]$

  integer, parameter:: R16P = r16 !modified to comply with mod_types - MAK
  integer, parameter:: R8P  = r8 !modified to comply with mod_types - MAK
  integer, parameter:: R4P  = selected_real_kind(6,37) ! 6 digits, 
                                                      !  range $[\pm 10^{-37} ,\pm 10^{+37} -1]$
  integer, parameter:: R_P  = R8P                         ! default real precision
  
  !
  !Integer precision definitions:
  !
  integer, parameter:: I8P  = selected_int_kind(18)       ! range $[-2^{63} ,+2^{63}  -1]$
  integer, parameter:: I4P  = selected_int_kind(9)        ! range $[-2^{31} ,+2^{31}  -1]$
  integer, parameter:: I2P  = selected_int_kind(4)        ! range $[-2^{15} ,+2^{15}  -1]$
  integer, parameter:: I1P  = selected_int_kind(2)        ! range $[-2^{7}~~,+2^{7}~~ -1]$
  integer, parameter:: I_P  = I4P                         ! default integer precision
  
  !
  ! Real output formats:
  !
  character(10), parameter:: FR16P = '(E41.33E4)'         ! R16P  output format
  character(10), parameter:: FR8P  = '(E23.15E3)'         ! R8P   output format
  character(9),  parameter:: FR4P  = '(E14.6E2)'          ! R4P   output format
  character(10), parameter:: FR_P  = '(E23.15E3)'         ! R\_P  output format
 
  !
  ! Integer output formats:
  !
  character(5), parameter:: FI8P  = '(I21)'               ! I8P  output format
  character(5), parameter:: FI4P  = '(I12)'               ! I4P  output format
  character(4), parameter:: FI2P  = '(I7)'                ! I2P  output format
  character(4), parameter:: FI1P  = '(I5)'                ! I1P  output format
  character(5), parameter:: FI_P  = '(I12)'               ! I\_P output format
 
  !
  ! private variables:
  !
  integer(I4P), parameter:: maxlen       = 500         ! max number of characters os static string
  character(1), parameter:: end_rec      = char(10)    ! end-character for binary-record finalize
  integer(I4P), parameter:: f_out_ascii  = 0           ! ascii-output-format parameter identifier
  integer(I4P), parameter:: f_out_binary = 1           ! binary-output-format parameter identifier
  integer(I4P)::            f_out        = f_out_ascii ! current output-format 
                                                       ! (initialized to ascii format)
  character(len=maxlen)::   topology                   ! mesh topology
  integer(I4P)::            Unit_VTK                   ! internal logical unit
  integer(I4P)::            Unit_VTK_Append            ! internal logical unit for raw binary 
                                                       ! xml append file
  integer(I4P)::            N_Byte                     ! number of byte to be written/read
  real(R8P)::               tipo_R8                    ! prototype of R8P real
  real(R4P)::               tipo_R4                    ! prototype of R4P real
  integer(I8P)::            tipo_I8                    ! prototype of I8P integer
  integer(I4P)::            tipo_I4                    ! prototype of I4P integer
  integer(I2P)::            tipo_I2                    ! prototype of I2P integer
  integer(I1P)::            tipo_I1                    ! prototype of I1P integer
  integer(I4P)::            ioffset                    ! offset pointer
  integer(I4P)::            indent                     ! indent pointer

  !In the following chapters there is the API reference of all functions of \LIBVTKIO.
contains
  
  function GetUnit() result(Free_Unit)
    !----------------------------------------------------------------------------------------------
    !The GetUnit function is used to get a free logic unit. The users of \LIBVTKIO does not 
    !know which is the logical unit: \LIBVTKIO handels this information without boring the users. 
    !The logical unit used is safe-free: if the program calling \LIBVTKIO has others logical units 
    !used \LIBVTKIO will never use these units, but will choice one that is free.
    !---------------------------------------------------------------------------------------------

    implicit none
    !---------------------------------------------------------------------------------------------
    integer(I4P):: Free_Unit ! free logic unit
    integer(I4P):: n1        ! counter
    integer(I4P):: ios       ! inquiring flag
    logical(4)::   lopen     ! inquiring flag

    !The following is the code snippet of GetUnit function: the units 0, 5, 6, 9 and 
    !all non-free units are discarded.
    !
    Free_Unit = -1_I4P                                      ! initializing free logic unit
    n1=1_I4P                                                ! initializing counter
    do
       if ((n1/=5_I4P).AND.(n1/=6_I4P).AND.(n1/=9_I4P)) then
          inquire (unit=n1,opened=lopen,iostat=ios)           ! verify logic units
          if (ios==0_I4P) then
             if (.NOT.lopen) then
                Free_Unit = n1                                  ! assignment of free logic
                return
             endif
          endif
       endif
       n1=n1+1_I4P                                           ! updating counter
    enddo
    return
  
    !
    !GetUnit function is private and cannot be called outside \LIBVTKIO. If you are interested 
    !to use it change its scope to public.

  endfunction GetUnit

  function Upper_Case(string)
    !---------------------------------------------------------------------------------------------
    !The Upper\_Case function converts the lower case characters of a string to upper case one. 
    !\LIBVTKIO uses this function in
    !order to achieve case-insensitive: all character variables used within \LIBVTKIO functions 
    !are pre-processed by Uppper\_Case function before these variables are used. So the users 
    !can call \LIBVTKIO functions whitout pay attention of the
    !case of the kwywords passed to the functions: calling the function VTK\_INI with the string 
    !\code{E_IO = vtk_ini('Ascii',...)}
    !or with the string  \code{E_IO = vtk_ini('AscII',...)} is equivalent.
    !---------------------------------------------------------------------------------------------

    implicit none

    !--------------------------------------------------------------------------------------------
    character(len=*), intent(IN):: string     ! string to be converted
    character(len=len(string))::   Upper_Case ! converted string
    integer::                      n1         ! characters counter
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    !The following is the code snippet of Upper\_Case function.
    !
    !(\doc)codesnippet
    Upper_Case = string
    do n1=1,len(string)
       select case(ichar(string(n1:n1)))
       case(97:122)
          Upper_Case(n1:n1)=char(ichar(string(n1:n1))-32) ! Upper case conversion
       endselect
    enddo
    return
    !---------------------------------------------------------------------------------------------
  endfunction Upper_Case
  
  subroutine vtk_ini(output_format,filename,title,mesh_topology,time_value)
    !---------------------------------------------------------------------------------------------
    !The VTK\_INI subroutine is used for initializing file. This subroutine must 
    !be the first to be called.
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    character(*), intent(IN):: output_format ! output format: ASCII or BINARY
    character(*), intent(IN):: filename      ! name of file
    character(*), intent(IN):: title         ! title
    character(*), intent(IN):: mesh_topology ! mesh topology
    integer(I4P)            :: E_IO          ! Input/Output inquiring flag: $0$ if IO is done, 
                                             ! $> 0$ if IO is not done
    real(R8P),    intent(IN):: time_value    !Physical time of the simulation
    character(len=maxlen)   :: s_buffer      ! buffer string
    !---------------------------------------------------------------------------------------------
    topology = trim(mesh_topology)
    Unit_VTK=GetUnit()
    select case(trim(Upper_Case(output_format)))
    case('ASCII')
       print*,' VTK format: ', output_format,' filename: ',filename
       f_out = f_out_ascii
       open(unit     = Unit_VTK,       &
            file     = trim(filename), &
            form     = 'FORMATTED',    &
            access   = 'SEQUENTIAL',   &
            action   = 'WRITE',        &
#ifdef __INTEL_COMPILER
            buffered = 'YES',          &
#endif
            iostat   = E_IO)
       ! writing header of file
       write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'# vtk DataFile Version 3.0'
       write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)trim(title)
       write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)trim(Upper_Case(output_format))
       write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'DATASET '//trim(topology)
    case('BINARY')
       print*,' VTK format: ', output_format,' filename: ',filename
       f_out = f_out_binary
       open(unit       = Unit_VTK,       &
            file       = trim(filename), &
            form       = 'UNFORMATTED',  &
            access     = 'SEQUENTIAL',   &
            action     = 'WRITE',        &
            convert    = 'BIG_ENDIAN',   &
#ifdef __INTEL_COMPILER
            recordtype = 'STREAM',       &
            buffered   = 'YES',          &
#endif
            iostat     = E_IO)
       ! writing header of file
       write(unit=Unit_VTK,iostat=E_IO)'# vtk DataFile Version 3.0'//end_rec
       write(unit=Unit_VTK,iostat=E_IO)trim(title)//end_rec
       write(unit=Unit_VTK,iostat=E_IO)trim(Upper_Case(output_format))//end_rec
       write(unit=Unit_VTK,iostat=E_IO)'DATASET '//trim(topology)//end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_ini
  
  subroutine vtk_geo_STRP_R8(Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz)
    !---------------------------------------------------------------------------------------------
    ! Subroutine for saving mesh; topology = STRUCTURED\_POINTS (R8P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P), intent(IN):: Nx        ! number of nodes in x direction
    integer(I4P), intent(IN):: Ny        ! number of nodes in y direction
    integer(I4P), intent(IN):: Nz        ! number of nodes in z direction
    real(R8P),    intent(IN):: X0        ! x coordinate of origin
    real(R8P),    intent(IN):: Y0        ! y coordinate of origin
    real(R8P),    intent(IN):: Z0        ! z coordinate of origin
    real(R8P),    intent(IN):: Dx        ! space step in x direction
    real(R8P),    intent(IN):: Dy        ! space step in y direction
    real(R8P),    intent(IN):: Dz        ! space step in z direction
    integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, 
                                         ! $> 0$ if IO is not done
    character(len=maxlen)::    s_buffer  ! buffer string
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,fmt='(A,3'//FR8P//')', iostat=E_IO)'ORIGIN ',X0,Y0,Z0
       write(unit=Unit_VTK,fmt='(A,3'//FR8P//')', iostat=E_IO)'SPACING ',Dx,Dy,Dz
    case(f_out_binary)
       write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(s_buffer,     fmt='(A,3'//FR8P//')', iostat=E_IO)'ORIGIN ',X0,Y0,Z0
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(s_buffer,     fmt='(A,3'//FR8P//')', iostat=E_IO)'SPACING ',Dx,Dy,Dz
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_geo_STRP_R8

  subroutine vtk_geo_STRP_R4(Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz)
    !---------------------------------------------------------------------------------------------
    ! Subroutine for saving mesh; topology = STRUCTURED\_POINTS (R4P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P), intent(IN):: Nx        ! number of nodes in x direction
    integer(I4P), intent(IN):: Ny        ! number of nodes in y direction
    integer(I4P), intent(IN):: Nz        ! number of nodes in z direction
    real(R4P),    intent(IN):: X0        ! x coordinate of origin
    real(R4P),    intent(IN):: Y0        ! y coordinate of origin
    real(R4P),    intent(IN):: Z0        ! z coordinate of origin
    real(R4P),    intent(IN):: Dx        ! space step in x direction
    real(R4P),    intent(IN):: Dy        ! space step in y direction
    real(R4P),    intent(IN):: Dz        ! space step in z direction
    integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, 
                                         ! $> 0$ if IO is not done
    character(len=maxlen)::    s_buffer  ! buffer string
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,fmt='(A,3'//FR4P//')', iostat=E_IO)'ORIGIN ',X0,Y0,Z0
       write(unit=Unit_VTK,fmt='(A,3'//FR4P//')', iostat=E_IO)'SPACING ',Dx,Dy,Dz
    case(f_out_binary)
       write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(s_buffer,     fmt='(A,3'//FR4P//')', iostat=E_IO)'ORIGIN ',X0,Y0,Z0
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(s_buffer,     fmt='(A,3'//FR4P//')', iostat=E_IO)'SPACING ',Dx,Dy,Dz
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_geo_STRP_R4

  subroutine vtk_geo_STRG_R8(Nx,Ny,Nz,NN,X,Y,Z)
    !---------------------------------------------------------------------------------------------
    ! Subroutine for saving mesh; topology = STRUCTURED\_GRID (R8P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P), intent(IN):: Nx       ! number of nodes in x direction
    integer(I4P), intent(IN):: Ny       ! number of nodes in y direction
    integer(I4P), intent(IN):: Nz       ! number of nodes in z direction
    integer(I4P), intent(IN):: NN       ! number of all nodes
    real(R8P),    intent(IN):: X(1:NN)  ! x coordinates
    real(R8P),    intent(IN):: Y(1:NN)  ! y coordinates
    real(R8P),    intent(IN):: Z(1:NN)  ! z coordinates
    integer(I4P)::             E_IO     ! Input/Output inquiring flag: $0$ if IO is done, 
                                        ! $> 0$ if IO is not done
    character(len=maxlen)::    s_buffer ! buffer string
    integer(I4P)::             n1       ! counter
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
       write(unit=Unit_VTK,fmt='(3'//FR8P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    case(f_out_binary)
       write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
       write(unit=Unit_VTK,                       iostat=E_IO)end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_geo_STRG_R8

  subroutine vtk_geo_STRG_R4(Nx,Ny,Nz,NN,X,Y,Z)
    !---------------------------------------------------------------------------------------------
    ! Subroutine for saving mesh; topology = STRUCTURED\_GRID (R4P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P), intent(IN):: Nx       ! number of nodes in x direction
    integer(I4P), intent(IN):: Ny       ! number of nodes in y direction
    integer(I4P), intent(IN):: Nz       ! number of nodes in z direction
    integer(I4P), intent(IN):: NN       ! number of all nodes
    real(R4P),    intent(IN):: X(1:NN)  ! x coordinates
    real(R4P),    intent(IN):: Y(1:NN)  ! y coordinates
    real(R4P),    intent(IN):: Z(1:NN)  ! z coordinates
    integer(I4P)::             E_IO     ! Input/Output inquiring flag: $0$ if IO is done, 
                                        ! $> 0$ if IO is not done
    character(len=maxlen)::    s_buffer ! buffer string
    integer(I4P)::             n1       ! counter
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
       write(unit=Unit_VTK,fmt='(3'//FR4P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    case(f_out_binary)
       write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
       write(unit=Unit_VTK,                       iostat=E_IO)end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_geo_STRG_R4

  subroutine vtk_geo_RECT_R8(Nx,Ny,Nz,X,Y,Z)
    !---------------------------------------------------------------------------------------------
    ! Subroutine for saving mesh; topology = RECTILINEAR\_GRID (R8P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P), intent(IN):: Nx        ! number of nodes in x direction
    integer(I4P), intent(IN):: Ny        ! number of nodes in y direction
    integer(I4P), intent(IN):: Nz        ! number of nodes in z direction
    real(R8P),    intent(IN):: X(1:Nx)   ! x coordinates
    real(R8P),    intent(IN):: Y(1:Ny)   ! y coordinates
    real(R8P),    intent(IN):: Z(1:Nz)   ! z coordinates
    integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, 
                                         ! $> 0$ if IO is not done
    character(len=maxlen)::    s_buffer  ! buffer string
    integer(I4P)::             n1        ! counter
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' double'
       write(unit=Unit_VTK,fmt=FR8P,              iostat=E_IO)(X(n1),n1=1,Nx)
       write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' double'
       write(unit=Unit_VTK,fmt=FR8P,              iostat=E_IO)(Y(n1),n1=1,Ny)
       write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' double'
       write(unit=Unit_VTK,fmt=FR8P,              iostat=E_IO)(Z(n1),n1=1,Nz)
    case(f_out_binary)
       write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' double'
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),n1=1,Nx)
       write(unit=Unit_VTK,                       iostat=E_IO)end_rec
       write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' double'
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK,                       iostat=E_IO)(Y(n1),n1=1,Ny)
       write(unit=Unit_VTK,                       iostat=E_IO)end_rec
       write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' double'
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK,                       iostat=E_IO)(Z(n1),n1=1,Nz)
       write(unit=Unit_VTK,                       iostat=E_IO)end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_geo_RECT_R8
  
  subroutine vtk_geo_unst_R8(NN,X,Y,Z) 
    !---------------------------------------------------------------------------------------------
    ! Subroutine for saving mesh; topology = UNSTRUCTURED\_GRID (R8P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P), intent(IN):: NN        ! number of nodes
    real(R8P),    intent(IN):: X(1:NN)   ! x coordinates of all nodes
    real(R8P),    intent(IN):: Y(1:NN)   ! y coordinates of all nodes
    real(R8P),    intent(IN):: Z(1:NN)   ! z coordinates of all nodes
    integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, 
                                         ! $> 0$ if IO is not done
    character(len=maxlen)::    s_buffer  ! buffer string
    integer(I4P)::             n1        ! counter
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
       write(unit=Unit_VTK,fmt='(3'//FR8P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    case(f_out_binary)
       write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
       write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
       !write(*,*)(X(n1),Y(n1),Z(n1),n1=1,NN)
       write(unit=Unit_VTK,                       iostat=E_IO)end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_geo_unst_R8

  subroutine vtk_con(NC,nsize,ncon,connect,cell_type)
    !---------------------------------------------------------------------------------------------
    !Subroutine to write the connectivity matrix to the binary file
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P), intent(IN):: NC,ncon,nsize   ! number of cells
    integer(I4P) connect(ncon)                 ! mesh connectivity
    integer(I4P) cell_type(NC)                 ! VTK cell type
    integer(I4P)::             E_IO            ! Input/Output inquiring flag: $0$ if IO is done, 
                                               ! $> 0$ if IO is not done
    character(len=maxlen)::    s_buffer        ! buffer string
    integer(I4P)::             i
    !---------------------------------------------------------------------------------------------
    !ncon = size(connect,1)
    select case(f_out)
    case(f_out_ascii)
       write(unit=Unit_VTK,fmt='(A,2'//FI4P//')',iostat=E_IO)'CELLS ',NC,ncon
       write(unit=Unit_VTK,fmt=FI4P,             iostat=E_IO)connect
       write(unit=Unit_VTK,fmt='(A,'//FI4P//')', iostat=E_IO)'CELL_TYPES ',NC
       write(unit=Unit_VTK,fmt=FI4P,             iostat=E_IO)cell_type
    case(f_out_binary)
       write(s_buffer,     fmt='(A,2'//FI4P//')',iostat=E_IO)'CELLS ',NC,nsize
       write(unit=Unit_VTK,                      iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK,                       iostat=E_IO)(connect(i),i=1,ncon)
       !write(*,*)(connect(i),i=1,ncon)
       write(unit=Unit_VTK,                      iostat=E_IO)end_rec
       write(s_buffer,     fmt='(A,'//FI4P//')', iostat=E_IO)'CELL_TYPES ',NC
       write(unit=Unit_VTK,                      iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK,                      iostat=E_IO)(cell_type(i),i=1,NC)
       write(unit=Unit_VTK,                      iostat=E_IO)end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_con

  subroutine vtk_dat(NC_NN, var_location)
    !---------------------------------------------------------------------------------------------
    !@brief This subroutine is called before saving the data related to the 
    !geometric mesh. This subroutine initializes the data variables.
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes of field
    character(*), intent(IN):: var_location ! location of saving variables: cell for cell-centered,
                                            ! node for node-centered
    integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, 
                                            ! $> 0$ if IO is not done
    character(len=maxlen)::    s_buffer     ! buffer string
    !---------------------------------------------------------------------------------------------
    
    select case(f_out)
    case(f_out_ascii)
       select case(trim(Upper_Case(var_location)))
       case('CELL')
          write(unit=Unit_VTK,fmt='(A,'//FI4P//')',iostat=E_IO)'CELL_DATA ',NC_NN
       case('NODE')
          write(unit=Unit_VTK,fmt='(A,'//FI4P//')',iostat=E_IO)'POINT_DATA ',NC_NN
       endselect
    case(f_out_binary)
       select case(trim(Upper_Case(var_location)))
       case('CELL')
          write(s_buffer,     fmt='(A,'//FI4P//')',iostat=E_IO)'CELL_DATA ',NC_NN
          write(unit=Unit_VTK,                     iostat=E_IO)trim(s_buffer)//end_rec
       case('NODE')
          write(s_buffer,     fmt='(A,'//FI4P//')',iostat=E_IO)'POINT_DATA ',NC_NN
          write(unit=Unit_VTK,                     iostat=E_IO)trim(s_buffer)//end_rec
       endselect
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_dat
  
  subroutine vtk_var_scal_R8(NC_NN,variablename,variable)
    !---------------------------------------------------------------------------------------------
    ! Subroutine to write the field of a scalar variable (R8P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    integer(I4P)            :: n1
    integer(I4P), intent(IN):: NC_NN        ! number of nodes or cells
    character(*), intent(IN):: variablename      ! variable name
    character(len=maxlen)   :: s_buffer     ! buffer string
    real(R8P),    intent(IN):: variable(1:NC_NN) ! variable to be saved
    integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, 
                                            ! $> 0$ if IO is not done

    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'SCALARS '//trim(variablename)//' double 1'
       write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'LOOKUP_TABLE default'
       write(unit=Unit_VTK,fmt=FR8P, iostat=E_IO)variable
    case(f_out_binary)
       write(s_buffer,     fmt='(A,A,A)',iostat=E_IO)'SCALARS '//trim(variablename)//' double 1'
       write(unit=Unit_VTK              ,iostat=E_IO)trim(s_buffer)//end_rec
       write(s_buffer,     fmt='(A)'    ,iostat=E_IO)'LOOKUP_TABLE default'
       write(unit=Unit_VTK              ,iostat=E_IO)trim(s_buffer)//end_rec
       write(unit=Unit_VTK              ,iostat=E_IO)(variable(n1),n1=1,NC_NN)
       !write(*,*) 'variable: ', n1, (variable(n1),n1=1,NC_NN)
       write(unit=Unit_VTK              ,iostat=E_IO)end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_var_scal_R8
  
  subroutine vtk_var_vect_R8(vec_type,NC_NN,varname,varX,varY,varZ)
    !---------------------------------------------------------------------------------------------
    ! Subroutine to write the field of a vector variable (R8P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    character(*), intent(IN):: vec_type      ! vector type: vect = generic vector, 
                                             ! norm = normal vector
    integer(I4P), intent(IN):: NC_NN         ! number of nodes or cells
    character(*), intent(IN):: varname       ! variable name
    real(R8P),    intent(IN):: varX(1:NC_NN) ! x component of vector
    real(R8P),    intent(IN):: varY(1:NC_NN) ! y component of vector
    real(R8P),    intent(IN):: varZ(1:NC_NN) ! z component of vector
    integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, 
                                             ! $> 0$ if IO is not done
    integer(I8P)::             n1            ! counter
    !---------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       select case(Upper_Case(trim(vec_type)))
       case('VECT')
          write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'VECTORS '//trim(varname)//' double'
       case('NORM')
          write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'NORMALS '//trim(varname)//' double'
       endselect
       write(unit=Unit_VTK,fmt='(3'//FR8P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    case(f_out_binary)
       select case(Upper_Case(trim(vec_type)))
       case('VECT')
          write(unit=Unit_VTK,iostat=E_IO)'VECTORS '//trim(varname)//' double'//end_rec
       case('NORM')
          write(unit=Unit_VTK,iostat=E_IO)'NORMALS '//trim(varname)//' double'//end_rec
       endselect
       write(unit=Unit_VTK,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
       write(unit=Unit_VTK,iostat=E_IO)end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_var_vect_R8

  subroutine vtk_var_vect2D_R8(vec_type,NC_NN,varname,varX,varY)
    !---------------------------------------------------------------------------------------------
    ! Subroutine to write the field of a vector variable (R8P).
    !---------------------------------------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------------------------------------
    character(*), intent(IN):: vec_type      ! vector type: vect = generic vector , 
                                             ! norm = normal vector
    integer(I4P), intent(IN):: NC_NN         ! number of nodes or cells
    character(*), intent(IN):: varname       ! variable name
    real(R8P),    intent(IN):: varX(1:NC_NN) ! x component of vector
    real(R8P),    intent(IN):: varY(1:NC_NN) ! y component of vector
    integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, 
                                             ! $> 0$ if IO is not done
    integer(I8P)::             n1            ! counter
    !---------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------
    select case(f_out)
    case(f_out_ascii)
       select case(Upper_Case(trim(vec_type)))
       case('VECT')
          write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'VECTORS '//trim(varname)//' double'
       case('NORM')
          write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'NORMALS '//trim(varname)//' double'
       endselect
       write(unit=Unit_VTK,fmt='(2E41.33)',iostat=E_IO)(varX(n1),varY(n1),n1=1,NC_NN)
   !  case(f_out_binary)
   !     select case(Upper_Case(trim(vec_type)))
   !     case('VECT')
   !        write(unit=Unit_VTK,iostat=E_IO)'VECTORS '//trim(varname)//' double'//end_rec
   !     case('NORM')
   !        write(unit=Unit_VTK,iostat=E_IO)'NORMALS '//trim(varname)//' double'//end_rec
   !     endselect
   !     write(unit=Unit_VTK,iostat=E_IO)(varX(n1),  varY(n1),n1=1,NC_NN)
   !     write(unit=Unit_VTK,iostat=E_IO)end_rec
    endselect
    
    !---------------------------------------------------------------------------------------------
  end subroutine vtk_var_vect2D_R8
  
  subroutine vtk_end()
    !---------------------------------------------------------------------------------------------
    !@brief This subroutine is used to finalize the file opened and it has not inputs. 
    !The \LIBVTKIO manages the file unit without the user's action.
    !---------------------------------------------------------------------------------------------

    implicit none

    integer(I4P):: E_IO ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
    
    close(unit=Unit_VTK,iostat=E_IO)
    
  end subroutine vtk_end
  
end module mod_vtk_binary
