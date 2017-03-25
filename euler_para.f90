

! ======= EULER POSITIONS ==============================================================================
    subroutine EulerPositions(positions,velocities,forces,particles,dimnsion,BoxSize,mass,deltat)

    use paralelizar
    include 'mpif.h'

    real, dimension(:,:)  ::  positions, velocities, forces
    integer             ::  dimnsion, particles, i, request, Master = 0, iproc, partxproc
    real                ::  BoxSize, deltat, mass
    double precision    ::  start_time, lapso_time, end_time

    ! Update of the positions using Euler's algorithm
    start_time = MPI_Wtime()

    do i = ini(rank), fin(rank)
        positions(i,:) = positions(i,:) + ( velocities(i,:) + deltat * forces(i,:) / mass ) * deltat
    end do

    lapso_time = MPI_Wtime()

    ! Sending the calculations to the Master
    if (rank /= Master) then
        call MPI_ISEND(positions(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL,Master,1,MPI_COMM_WORLD,request,ierror)
    end if

    ! Wait for all workers
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Master receive from workers and merge all the array
    if ( rank == Master ) then

        ! Master receive and merge
        do iproc = 1, numproc-1
            call MPI_RECV(positions(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL,iproc,1,MPI_COMM_WORLD,request,ierror)
        end do

       ! Update now all the workers
        do iproc = 1, numproc-1
            call MPI_ISEND(positions(:,:), 3*particles, MPI_REAL, iproc,1,MPI_COMM_WORLD,request,ierror)
        end do

    end if

    ! Wait again for all the processors and synchronize
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! All the workers receive the updated positions from the master
    if (rank /= Master ) then
      call MPI_RECV(positions(:,:),3*particles,MPI_REAL,Master,1,MPI_COMM_WORLD,request,ierror)
    end if

    end_time = MPI_Wtime()

    ! Apply the PBC : Reput particles inside the box
    call Refold_Positions(positions,particles,dimnsion,BoxSize)

    end subroutine EulerPositions



! ==== EULER VELOCITIES =================================================================
    subroutine EulerVelocities(velocities,forces,particles,dimnsion,BoxSize,mass,deltat)

    use paralelizar
    include 'mpif.h'

    real,dimension(:,:)  ::  velocities, forces
    integer           ::  dimnsion, particles, i, request, Master = 0, iproc, partxproc
    real              ::  BoxSize, deltat, mass
    double precision  ::  start_time, lapso_time, end_time

    ! Update of the velocities using Euler's algorithm
    start_time = MPI_Wtime()

    do i = ini(rank), fin(rank)
     velocities(i,:) = velocities(i,:) + deltat * forces(i,:) / mass
    end do

    lapso_time = MPI_Wtime()

    ! Send the velocities from the workers to the Master
    if (rank /= Master) then
        call MPI_ISEND(velocities(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL,Master,1,MPI_COMM_WORLD,request,ierror)
    end if

    ! Wait until all workers have sent the velocities
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Master receive all the arrays of velocities and merge them
    if ( rank == Master ) then

        ! Master receive and merge
        do iproc = 1, numproc-1
            call MPI_RECV(velocities(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL,iproc,1,MPI_COMM_WORLD,request,ierror)
        end do

        ! Update the velocities to all the workers
        do iproc=1,numproc-1
            call MPI_ISEND(velocities(:,:), 3*particles, MPI_REAL, iproc,1,MPI_COMM_WORLD,request,ierror)
        end do

    end if

    ! Synchronize all processors
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Workers receive the updated velocities
    if (rank /= Master ) then
      call MPI_RECV(velocities(:,:),3*particles,MPI_REAL,Master,1,MPI_COMM_WORLD,request,ierror)
    end if

    end_time = MPI_Wtime()

    end subroutine EulerVelocities
