!********************************************************************************
!********************************** SUBROUTINE PBC ******************************
!********************************************************************************
!
!                                                          Lorena Vega Domínguez
!********************************************************************************
  subroutine Refold_Positions(pos,N,dimnsion,BoxSize)

    integer::dimnsion,N,i !N=Number of part.
    real*8:: BoxSize
    real*8,dimension(N,dimnsion):: pos !positions
 

    do i=1,N
     pos(i,:) = pos(i,:) - BoxSize*nint(pos(i,:)/BoxSize) !periodic conditions
    end do
     open (unit=11, File='coord_serie.dat')
    do i=1,N
     write(11,*) pos(i,:)
    end do
   close(11)

  end subroutine Refold_Positions