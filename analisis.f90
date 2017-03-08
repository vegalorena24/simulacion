module analisis
implicit none

contains
function T_compute(N,V) result(temperatura)
    integer(8) :: N!Nombre particules
    real(8),dimension(N,3) :: Velocitats!Matriu de velocitats
    real(8) :: temperatura
    temperatura = sum(Velocitats**2.0d0)/(3*N)
end function

function P_compute(N,L,R,F) result(P)
    !Calcula sa presio donada ses coordenades i ses forces
    integer(8) :: N,i
    real(8),dimension(N,3) :: F
    real(8),dimension(N,3) :: R
    real(8) :: P,L
    P = 0.0d0
    
    !Calculem es treball i ho fiquem a P
    do i = 1,N
        P = P + dot_product(F(i,:),R(i,:))
    enddo
    !Acabem de calcular sa pressio
    P = (P/3.0d0 + dble(N)*T)/L**3.0d0
end function
end module
