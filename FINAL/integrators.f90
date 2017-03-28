module integrators

use paralelizar

contains

! ======= EULER POSITIONS ==============================================================================
subroutine EulerPositions(pos,vel,forces,N,dimnsion,BoxSize,mass,deltat,ini,fin)

    integer, intent(in)                      ::  dimnsion, N
    integer, dimension(0:), intent(in)       ::  ini, fin
    real*8, intent(in)                       ::  BoxSize, deltat, mass
    real*8, dimension(:,:), intent(inout)    ::  pos
    real*8, dimension(:,:), intent(in)       ::  vel, forces
    double precision                         ::  start_time, lapso_time, end_time

    ! Update of the positions using Euler
    start_time=MPI_Wtime()
    do i=ini(rank),fin(rank)
        pos(i,:) = pos(i,:) + ( vel(i,:) + deltat * forces(i,:) / mass ) * deltat
    end do

    ! Send the positions from workers to Master
    if (rank /= MASTER) then
    call MPI_ISEND(pos(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),&
                         MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if

    ! Wait for all workers
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Master receives and merge all the calculations
    if ( rank == MASTER ) then
        do iproc=1,numproc-1
            call MPI_RECV(pos(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),&
                    MPI_REAL8,iproc,1,MPI_COMM_WORLD,stat,ierror)
        end do

        ! UPDATE
        do iproc=1,numproc-1
            call MPI_ISEND(pos(:,:), 3*N, MPI_REAL8, iproc,1,MPI_COMM_WORLD,request,ierror)
        end do
    end if

    ! Synchronize all processors
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! All workers receive the positions once updated and merged
    if (rank /= MASTER ) then
        call MPI_RECV(pos(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
    end if
    end_time=MPI_Wtime()

    !if (rank == MASTER) then
    !  print *,"Euler Positions Time : ",end_time-start_time," seconds"
    !endif

end subroutine EulerPositions


! ==== EULER VELOCITIES =================================================================
subroutine EulerVelocities(vel,forces,N,dimnsion,BoxSize,mass,deltat,ini,fin)

    integer, intent(in)                     ::  dimnsion, N
    integer                                 ::  j, i
    integer, dimension(0:), intent(in)      ::  ini, fin
    real*8, intent(in)                      ::  BoxSize, deltat, mass
    real*8, dimension(:,:), intent(inout)   ::  vel
    real*8, dimension(:,:), intent(in)      ::  forces
    double precision                        ::  start_time, lapso_time, end_time

    ! Update of the velocities using Euler
    start_time=MPI_Wtime()
    do i=ini(rank), fin(rank)
        vel(i,:) = vel(i,:) + deltat*forces(i,:)/mass
    end do

    ! Send velocities from the workers to the Main
    if (rank /= MASTER) then
        call MPI_ISEND(vel(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),&
                    MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if

    ! Waiting all workers
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Master receive and merge
    if ( rank == MASTER ) then
         do iproc=1,numproc-1
             call MPI_RECV(vel(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),&
                        MPI_REAL8,iproc,1,MPI_COMM_WORLD,stat,ierror)
         end do

        ! UPDATE
        do iproc=1,numproc-1
            call MPI_ISEND(vel(:,:), 3*N, MPI_REAL8, iproc,1,MPI_COMM_WORLD,request,ierror)
        end do
    end if

    ! Synchronize all processors
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! All workers receive the velocities once updated
    if (rank /= MASTER ) then
      call MPI_RECV(vel(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
    end if
    end_time=MPI_Wtime()

    !if (rank == MASTER) then
    !  print *,"Euler Velocities Time : ",end_time-start_time," seconds"
    !endif

end subroutine EulerVelocities


! ======================================================================================
subroutine Refold_Positions(pos,N,dimnsion,BoxSize,ini,fin)

    integer, intent(in)                    ::  dimnsion, N !N=Number of part.
    integer                                ::  i
    real*8, intent(in)                     ::  BoxSize
    integer, dimension(0:), intent(in)     ::  ini, fin
    real*8,dimension(:,:), intent(inout)   ::  pos
    double precision                       ::  start_time, lapso_time,end_time

    !start calculation
    start_time=MPI_Wtime()
    do i=ini(rank),fin(rank)
        pos(i,:) = pos(i,:) - BoxSize*nint(pos(i,:)/BoxSize) !periodic conditions
    end do

    !sending work of worker(i) to master
    if (rank /= MASTER) then
        call MPI_ISEND(pos(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),&
                    MPI_REAL8,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if

    !waiting all workers
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    !master receive and merge
    if ( rank == MASTER ) then
         do iproc=1,numproc-1
             call MPI_RECV(pos(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),&
                        MPI_REAL8,iproc,1,MPI_COMM_WORLD,stat,ierror)
         end do

        !UPDATE
        do iproc=1,numproc-1
            call MPI_ISEND(pos(:,:), 3*N, MPI_REAL8, iproc,1,MPI_COMM_WORLD,request,ierror)
        end do
    end if

    !Sincronize all processors
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !all workers receive the coord
    if (rank /= MASTER ) then
      call MPI_RECV(pos(:,:),3*N,MPI_REAL8,MASTER,1,MPI_COMM_WORLD,stat,ierror)
    end if
    end_time=MPI_Wtime()

    !if (rank == MASTER) then
    !  print *,"Refolf positions Time : ",end_time-start_time," seconds"
    !endif

end subroutine Refold_Positions

! ===================================================================================

end module
