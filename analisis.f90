module analisis
implicit none
include 'mpif.h'

contains

function P_compute_paralel(N,L,R,Forces,temperatura,rank,numproc,ierror,MASTER) result(Presion)
    !Calcula sa presio donada ses coordenades i ses forces
    integer(4) :: N,i
    real(4),dimension(N,3) :: Forces
    real(4),dimension(N,3) :: R
    real(4) :: Presion,L,Presion_parcial,temperatura!,presion
    !Declarem variables MPI
    integer :: comm,rank,numproc,ierror,MASTER
    !integer :: MASTER
    integer(4) :: interval_index,from_i,to_i
    Presion = 0.0
    
    !Paralelitzem calcul de traball per calcular P
    interval_index = N/numproc
    from_i = (rank)*interval_index+1
    to_i = (rank+1)*interval_index
    print*,from_i,to_i
    if (rank/=numproc-1)then
        Presion_parcial = 0
        print*,from_i,to_i
        do i=from_i,to_i
            Presion_parcial = Presion_parcial + dot_product(Forces(i,:),R(i,:))
        enddo
    else
        from_i = (rank)*interval_index+1
        Presion_parcial = 0
        do i=from_i,N
            Presion_parcial = Presion_parcial + dot_product(Forces(i,:),R(i,:))
        enddo
    endif
    !Sumo P parciales al MASTER
    call MPI_REDUCE(Presion_parcial,Presion,1,MPI_FLOAT,MPI_SUM,MASTER,MPI_COMM_WORLD,ierror)
    
    !Master acaba de calcular sa temperatura
    if(rank==MASTER)then
        Presion = (Presion/3 + dble(N)*temperatura)/L**3
    endif
end function


function T_compute_paralel(N,Velocitats,rank,numproc,ierror,MASTER) result(temperatura)
    !Funcio que calcula sa temperatura paralelament
    integer(4) :: N!Nombre particules!!
    real(4),dimension(N,3) :: Velocitats!Matriu de velocitats
    real(4) :: temperatura,temperatura_parcial
    !Declarem variables MPI
    integer :: rank,numproc,ierror,MASTER
    integer(4) :: interval_index
    
    !Paralelitzem suma de quadrats de Velocitats
    interval_index = N/numproc
    if (rank/=numproc-1)then
        temperatura_parcial = sum(Velocitats((rank)*interval_index+1:(rank+1)*interval_index,:)**2)
    else
        temperatura_parcial = sum(Velocitats((rank)*interval_index+1:,:)**2)
    endif
    call MPI_REDUCE(temperatura_parcial,temperatura,1,MPI_REAL,MPI_SUM,MASTER,MPI_COMM_WORLD,ierror)
    
    !Master acaba de calcular sa temperatura
    if(rank==MASTER)then
        temperatura = temperatura/(3*N)
    endif
end function


function T_compute(N,Velocitats) result(temperatura)
    !Funcio que calcula sa temperatura
    integer(8) :: N!Nombre particules!!
    real(8),dimension(N,3) :: Velocitats!Matriu de velocitats
    real(8) :: temperatura
    temperatura = sum(Velocitats**2.0d0)/(3*N)
end function

function P_compute(N,L,R,Forces,temperatura) result(Presion)
    !Calcula sa presio donada ses coordenades i ses forces
    integer(8) :: N,i
    real(8),dimension(N,3) :: Forces
    real(8),dimension(N,3) :: R
    real(8) :: Presion,L,temperatura
    integer :: comm,rank,numproc,ierror,MASTER
    Presion = 0.0d0
    
    !Calculem es treball i ho fiquem a P
    do i = 1,N
        Presion = Presion + dot_product(Forces(i,:),R(i,:))
    enddo
    !Acabem de calcular sa pressio
    Presion = (Presion/3.0d0 + dble(N)*temperatura)/L**3.0d0
end function
end module
