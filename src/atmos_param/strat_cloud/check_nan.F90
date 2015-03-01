MODULE check_nan_mod

use fms_mod, only :  error_mesg, FATAL, stdout, write_version_number

implicit none
private

interface check_nan
   module procedure check_nan_3d, check_nan_2d, check_nan_1d, check_nan_0d
end interface check_nan

public check_nan, check_nan_init

!---------------version number---------------------------
Character(len=128) :: Version = '$Id: check_nan.F90,v 1.1.2.1 2011/03/02 08:18:07 Richard.Hemler Exp $'
Character(len=128) :: Tagname = '$Name: testing $'

!-------------------------------------------------------------------------

logical     :: module_is_initialized = .false.


CONTAINS


!#########################################################################

subroutine check_nan_init

      if (module_is_initialized) return

!------------------------------------------------------------------------
!    write version number to output file.
!------------------------------------------------------------------------
      call write_version_number (version, tagname)

!------------------------------------------------------------------------
!    declare module initialized.
!------------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------------

end subroutine check_nan_init


!#########################################################################

subroutine check_nan_3d (inarr, name)

real, dimension(:,:,:), intent(in)  :: inarr
character(len=*),       intent(in)  :: name

      integer :: i,j,k

!-----------------------------------------------------------------------
      do k=1,size(inarr,3)
        do j=1,size(inarr,2)
          do i=1,size(inarr,1)
            if (inarr(i,j,k) .ne. inarr(i,j,k)) then
              write(stdout(),*) " ------------------------------------- "
              write(stdout(),*) " NAN ERROR i,j,k ", i,j,k
              write(stdout(),*) " NAN ERROR msg ", name 
              call error_mesg ('check_nan_3d', 'found nan', FATAL)
            end if
          end do
        end do
      end do

!-------------------------------------------------------------------------

end subroutine check_nan_3d

!########################################################################

subroutine check_nan_2d (inarr, name)

real, dimension(:,:),   intent(in)  :: inarr
character(len=*),       intent(in)  :: name

      integer :: i,j

!------------------------------------------------------------------------
      do j=1,size(inarr,2)
        do i=1,size(inarr,1)
          if (inarr(i,j) .ne. inarr(i,j)) then
            write(stdout(),*) " ------------------------------------- "
            write(stdout(),*) " NAN ERROR i1,i2 ", i,j
            write(stdout(),*) " NAN ERROR msg ", name 
            call error_mesg ( 'check_nan_2d', 'found nan', FATAL)
          end if
        end do
      end do

!------------------------------------------------------------------------


end subroutine check_nan_2d



!#########################################################################

subroutine check_nan_1d (inarr, name)

real, dimension(:),     intent(in)  :: inarr
character(len=*),       intent(in)  :: name

      integer :: i

!------------------------------------------------------------------------
      do i=1,size(inarr,1)
        if (inarr(i) .ne. inarr(i)) then
          write(stdout(),*) " ------------------------------------- "
          write(stdout(),*) " NAN ERROR i1 ", i
          write(stdout(),*) " NAN ERROR msg ", name 
          call error_mesg ( 'check_nan_1d', 'found nan', FATAL)
        end if
      end do

!------------------------------------------------------------------------


end subroutine check_nan_1d



!#########################################################################

subroutine check_nan_0d (inv, name)

real,             intent (in)    :: inv
character(len=*), intent(in)     :: name


!------------------------------------------------------------------------
      if (inv .ne. inv) then
        write(stdout(),*) " ------------------------------------- "
        write(stdout(),*) " NAN ERROR msg ", name 
        call error_mesg ( 'check_nan_0d', 'found nan', FATAL)
      end if

!-------------------------------------------------------------------------  

end subroutine check_nan_0d



!##########################################################################


END MODULE check_nan_mod
