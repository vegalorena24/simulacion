
integer(8) function get_int_arg(posArg,defecto)
    !Funcion que coge el argumento i esimo y la convierte en un integer
    character(30) :: arg
    integer(8) :: defecto
    integer :: posArg
    call get_command_argument(posArg,arg)
    if(arg=='')then
        get_int_arg = defecto
    else
        read(arg,*) get_int_arg
    endif
end function
