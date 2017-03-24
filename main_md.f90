module paralelizar
integer:: ierror,rank,numproc
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
!datos de entrada de prueba
deltat=0.0032
Nsteps=10000
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

do step=1,Nsteps

 call forces(positions,BoxSize,accel,rc,epot)

 call IntegrationEuler(positions,velocities,accel,deltat,N,mass,dimnsion,BoxSize)

 !call sample

! ekin = sum(mass*norm2(velocities,2)/2.0)

 write(unit=123,fmt='(i10,3f20.10)') step, ekin+epot, ekin, epot

enddo
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

        !      atom-atom interactions

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

    subroutine IntegrationEuler(positions, velocities, forces, deltat,N,mass, dim, length)
    real, dimension(N,3)  ::  positions, velocities, forces
    real, intent(in)     ::  deltat, mass, length
    integer                 :: k
    integer, intent(in)  ::N, dim
    !integer             ::  i, j, counter, ierror

    ! Define the communicator
    	!comm = MPI_COMM_WORLD
    ! Find out the number of processes
    	!call MPI_COMM_SIZE(comm, numproc, ierror)
    ! Obtain the taskid, the id of an individual process
    	!call MPI_COMM_RANK(comm, taskid, ierror)

    !counter = N / numproc

    ! Integrate positions
    do k = 1, N
    positions(k,1) = positions(k,1) + ( velocities(k,1) + forces(k,1)*deltat ) * deltat
    positions(k,2) = positions(k,2) + ( velocities(k,2) + forces(k,2)*deltat ) * deltat
    positions(k,3) = positions(k,3) + ( velocities(k,3) + forces(k,3)*deltat ) * deltat
    enddo

    ! Call PBC subroutine (Lorena)
    call Refold_Positions(positions, N, dim, length)

    ! Integrate velocities
    do k = 1, N
    velocities(k,1) = velocities(k,1) + ( forces(k,1) / mass ) * deltat
    velocities(k,2) = velocities(k,2) + ( forces(k,2) / mass ) * deltat
    velocities(k,3) = velocities(k,3) + ( forces(k,3) / mass ) * deltat
    enddo

    end subroutine IntegrationEuler

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
     pos(i,:)=pos(i,:)-BoxSize*nint(pos(i,:)/BoxSize) !periodic conditions
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
      call MPI_RECV(pos(ini(rank):fin(rank),:),3*(fin(rank)-ini(rank)+1),MPI_REAL,iproc,1,MPI_COMM_WORLD,request,ierror)
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
    print*, "rank:", rank, "time calculation:", lapso_time-start_time,"seconds"
    if (rank == MASTER) then
     print*, "MASTER. UOPDATING INFO. TIME:", end_time-lapso_time,"seconds"
    end if

    
    end subroutine Refold_Positions

end program m
