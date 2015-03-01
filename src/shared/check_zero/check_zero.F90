MODULE check_zero_mod

use fms_mod, only :  error_mesg, FATAL, stdout, write_version_number, mpp_pe

implicit none
private

interface check_zero
   module procedure check_zero_3d, check_zero_2d, check_zero_1d, check_zero_0d
end interface check_zero

public check_zero

!-------------------------------------------------------------------------

CONTAINS

!#########################################################################

subroutine check_zero_3d (inarr, name)

real, dimension(:,:,:), intent(in)  :: inarr
character(len=*),       intent(in)  :: name

      integer :: i,j,k

!-----------------------------------------------------------------------
      do k=1,size(inarr,3)
        do j=1,size(inarr,2)
          do i=1,size(inarr,1)
            if (inarr(i,j,k) == 0.0) then
              write(stdout(),*) " ------------------------------------- "
              write(stdout(),*) " ZERO ERROR i,j,k=", i,j,k,"  pe=",mpp_pe()
              write(stdout(),*) " ZERO ERROR msg ", name 
              call error_mesg ('check_zero_3d', 'found zero', FATAL)
            end if
          end do
        end do
      end do
!-------------------------------------------------------------------------

end subroutine check_zero_3d

!########################################################################

subroutine check_zero_2d (inarr, name)

real, dimension(:,:),   intent(in)  :: inarr
character(len=*),       intent(in)  :: name

      integer :: i,j

!------------------------------------------------------------------------
      do j=1,size(inarr,2)
        do i=1,size(inarr,1)
          if (inarr(i,j) == 0.0) then
            write(stdout(),*) " ------------------------------------- "
            write(stdout(),*) " ZERO ERROR i1,i2=", i,j,"  pe=",mpp_pe()
            write(stdout(),*) " ZERO ERROR msg ", name 
            call error_mesg ( 'check_zero_2d', 'found zero', FATAL)
          end if
        end do
      end do
!------------------------------------------------------------------------

end subroutine check_zero_2d

!#########################################################################

subroutine check_zero_1d (inarr, name)

real, dimension(:),     intent(in)  :: inarr
character(len=*),       intent(in)  :: name

      integer :: i

!------------------------------------------------------------------------
      do i=1,size(inarr,1)
        if (inarr(i) == 0.0) then
          write(stdout(),*) " ------------------------------------- "
          write(stdout(),*) " ZERO ERROR i1=", i,"  pe=",mpp_pe()
          write(stdout(),*) " ZERO ERROR msg ", name 
          call error_mesg ( 'check_zero_1d', 'found zero', FATAL)
        end if
      end do
!------------------------------------------------------------------------

end subroutine check_zero_1d

!#########################################################################

subroutine check_zero_0d (inv, name)

real,             intent (in)    :: inv
character(len=*), intent(in)     :: name

!------------------------------------------------------------------------
      if (inv == 0.0) then
        write(stdout(),*) " ------------------------------------- "
        write(stdout(),*) " ZERO ERROR scalar ", name,"  pe=",mpp_pe()
        write(stdout(),*) " ZERO ERROR msg ", name 
        call error_mesg ( 'check_zero_0d', 'found zero', FATAL)
      end if
!-------------------------------------------------------------------------  

end subroutine check_zero_0d

!##########################################################################

END MODULE check_zero_mod
