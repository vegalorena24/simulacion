  module xabi 
 contains 
 subroutine forces(positions,boxlength,accel,rc,epot)
real*8, dimension(:,:), intent(in)  :: positions
real*8, dimension(:,:), intent(out) :: accel
real*8, intent(in)                  :: boxlength, rc
real*8, intent(out)                 :: epot
integer                             :: is, js, natoms
real*8                             :: pot

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
    real*8, dimension(:,:), intent(in)      :: positions
    real*8, dimension(:,:), intent(inout)   :: accel
    real*8, intent(in)                      :: boxlength, rc
    integer, intent(in)                     :: is, js
    real*8, intent(out)                     :: pot
    real*8, dimension(size(positions,2))         :: rij
    real*8                                  :: rr2, rijl, rr, forcedist
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
    end module