!********************************************************************************
!********************************** SUBROUTINE PBC ******************************
!********************************************************************************
!
!                                                          Lorena Vega Domínguez
!********************************************************************************



subroutine Refold_Positions(pos,N,dimnsion,BoxSize)
implicit none
integer::dimnsion,N,i !N=Number of part. 
real:: BoxSize
real,dimension(N,dimnsion):: pos !positions
do i=1,N
 pos(i,:)=pos(i,:)-BoxSize*nint(pos(i,:)/BoxSize) !get the part inside the box
end do
end subroutine Refold_Positions

