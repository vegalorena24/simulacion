module paralelizar
integer:: ierror,rank,numproc,partxproc,request, MASTER
integer, dimension(:),allocatable :: ini, fin
end module paralelizar

program m
use paralelizar
implicit none
include 'mpif.h'
real:: deltat, BoxSize, mass,rc,epot, ekin
integer:: N,dimnsion,Nsteps,i,j,step
real, dimension(:,:), allocatable:: positions,accel,velocities
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
MASTER=0
!if ( rank == MASTER ) then
!datos de entrada de prueba
deltat=0.01
Nsteps=100
N=256
dimnsion=3
BoxSize=6.1984
mass=1.0
rc=2.5


allocate (positions(N,dimnsion))
allocate (accel(N,dimnsion))
allocate (velocities(N,dimnsion))

open (unit=10, File='coordenadas.dat')

do i=1,N
 read(10,*) positions(i,:)
end do
close (10)

velocities = 0.0

!MAIN

!call init()

allocate (ini(0:numproc),fin(0:numproc))
   !distribucion de particulas
    partxproc=nint(real(N)/real(numproc))
    do j=0,numproc-2
       ini(j)=j*partxproc+1
       fin(j)=(j+1)*partxproc
    end do
    ini(numproc-1)=(numproc-1)*partxproc+1
    fin(numproc-1)=N

open(unit=123,file='energy.dat',status='replace',action='write')

accel = 10.0

do step=1,Nsteps

 !call forces(positions,BoxSize,accel,rc,epot)

 call EulerPositions(positions,velocities,accel,N,dimnsion,BoxSize,mass,deltat)
 call EulerVelocities(velocities,accel,N,dimnsion,BoxSize,mass,deltat)

 call Refold_Positions(positions,N,dimnsion,BoxSize)


 !call sample

! ekin = sum(mass*norm2(velocities,2)/2.0)

 write(unit=123,fmt='(i10,3f20.10)') step, ekin+epot, ekin, epot

if ( rank == MASTER ) then
     print *, '************************************************'
     print *, 'step ', step
     print *, ' r = ', positions(100,:)
     print *, ' v = ', velocities(100,:)
     print *, ' F = ', accel(100,:)
endif

enddo

!end if

call MPI_FINALIZE(ierror)

contains

subroutine forces(positions,boxlength,accel,rc,epot)
real, dimension(:,:), intent(in)  :: positions
real, dimension(:,:), intent(out) :: accel
real, intent(in)                  :: boxlength, rc
real, intent(out)                 :: epot
integer                             :: is, js, natoms
real                             :: pot


      natoms = size(positions,1)

        accel = 0.0d0
        epot = 0.d0

!              atom-atom interactions

    do is = 1,natoms-1
           do js = is+1,natoms
              call lj(is,js,positions,boxlength,accel,rc,pot)
              epot = epot + pot
           end do
        end do

    end subroutine

    subroutine lj(is,js,positions,boxlength,accel,rc,pot)
    real, dimension(:,:), intent(in)      :: positions
    real, dimension(:,:), intent(inout)   :: accel
    real, intent(in)                      :: boxlength, rc
    integer, intent(in)                     :: is, js
    real, intent(out)                     :: pot
    real, dimension(size(positions,2))         :: rij
    real                                  :: rr2, rijl, rr, forcedist
    integer                                 :: l, dim

    dim = size(positions,2)

    rr2 = 0.d0
    pot = 0.d0
    do l = 1,dim
       rijl = positions(js,l) - positions(is,l)
       rij(l) = rijl - boxlength*nint(rijl/boxlength)
       rr2 = rr2 + rij(l)*rij(l)
    end do

    rr = sqrt(rr2)

    if (rr.lt.rc) then
       forcedist = 24*(2/rr**14-1.0/rr**8)
        pot = 4.0*(1.0/rr**12-1.0/rr**6)
        do l = 1,dim
            accel(is,l) = accel(is,l) - forcedist*rij(l)
            accel(js,l) = accel(js,l) + forcedist*rij(l)
        end do
    end if

    end subroutine

    ! ======= EULER POSITIONS ==============================================================================
        subroutine EulerPositions(positions,velocities,forces,particles,dimnsion,BoxSize,mass,deltat)

        use paralelizar
        include 'mpif.h'

        integer             ::  dimnsion, particles, i, request, Master = 0, iproc, partxproc
        real                ::  BoxSize, deltat, mass
        double precision    ::  start_time, lapso_time, end_time
        real, dimension(:,:)  ::  positions, velocities, forces

        ! Update of the positions using Euler's algorithm
        !start_time = MPI_Wtime()

        do i = ini(rank), fin(rank)
            positions(i,:) = positions(i,:) + ( velocities(i,:) + deltat * forces(i,:) / mass ) * deltat
        end do

        !lapso_time = MPI_Wtime()

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

        !end_time = MPI_Wtime()

        end subroutine EulerPositions

    ! ==== EULER VELOCITIES =================================================================
        subroutine EulerVelocities(velocities,forces,particles,dimnsion,BoxSize,mass,deltat)

        use paralelizar
        include 'mpif.h'

        integer           ::  dimnsion, particles, i, request, Master = 0, iproc, partxproc
        real              ::  BoxSize, deltat, mass
        double precision  ::  start_time, lapso_time, end_time
        real, dimension (:,:)  ::  velocities, forces

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

! ======================================================================================
    subroutine Refold_Positions(pos,N,dimnsion,BoxSize)

    use paralelizar
    include 'mpif.h'
    integer::dimnsion,N,i,request,MASTER=0,iproc,partxproc !N=Number of part.
    real:: BoxSize
    real,dimension(N,dimnsion):: pos !positions
    double precision:: start_time, lapso_time,end_time

    !start calculation
    start_time=MPI_Wtime()
    do i=ini(rank),fin(rank)
     pos(i,:) = pos(i,:) - BoxSize*nint(pos(i,:)/BoxSize) !periodic conditions
    end do
    lapso_time=MPI_Wtime()

    !sending work of worker(i) to master
    if (rank /= MASTER) then
    call MPI_ISEND(pos(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if

    !waiting all workers
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    !master receive and merge
    if ( rank == MASTER ) then
     do iproc=1,numproc-1
      call MPI_RECV(pos(ini(iproc):fin(iproc),:),3*(fin(iproc)-ini(iproc)+1),MPI_REAL,iproc,1,MPI_COMM_WORLD,request,ierror)
     end do

   !UPDATE
    do iproc=1,numproc-1
     call MPI_ISEND(pos(:,:), 3*N, MPI_REAL, iproc,1,MPI_COMM_WORLD,request,ierror)
    end do
    end if

    !Sincronize all processors
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !all workers receive the coord
    if (rank /= MASTER ) then
      call MPI_RECV(pos(:,:),3*N,MPI_REAL,MASTER,1,MPI_COMM_WORLD,request,ierror)
    end if
    end_time=MPI_Wtime()
    !print*, "rank:", rank, "time calculation:", lapso_time-start_time,"seconds"
    if (rank == MASTER) then
        open (unit=11, File='coord_paral.dat')
        do i=1,N
        write(11,*) pos(i,:)
        end do
        close(11)
    ! print*, "MASTER. UOPDATING INFO. TIME:", end_time-lapso_time,"seconds"
    end if

    end subroutine Refold_Positions
! ======================================================================================


end program m
