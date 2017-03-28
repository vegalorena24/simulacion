module analisis2
contains


function T_compute(N,Velocitats) result(temperatura)
    !Funcio que calcula sa temperatura
    integer(4) :: N!Nombre particules!!
    real(8),dimension(N,3) :: Velocitats!Matriu de velocitats
    real(8) :: temperatura
    temperatura = sum(Velocitats**2.0d0)/(3*N)
end function

function P_compute(N,L,R,Forces,temperatura) result(Presion)
    !Calcula sa presio donada ses coordenades i ses forces
    integer(4) :: N,i
    real(8),dimension(N,3) :: Forces
    real(8),dimension(N,3) :: R
    real(8) :: Presion,L,temperatura
    Presion = 0.0d0
    
    !Calculem es treball i ho fiquem a P
    do i = 1,N
        Presion = Presion + dot_product(Forces(i,:),R(i,:))
    enddo
    !Acabem de calcular sa pressio
    Presion = (Presion/(N*3.0d0) + dble(N)*temperatura)/L**3.0d0
end function
end module
