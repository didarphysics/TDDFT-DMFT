module convert !conver integer/real to string
  IMPLICIT NONE
contains
  function dirname_int(number)
    integer*8, intent(in) :: number
    character(len=8)  :: dirname_int

    ! Cast the (rounded) number to string using 6 digits and
    ! leading zeros
!    write (dirname, '(I6.6)')  nint(number)
    ! This is the same w/o leading zeros
    write (dirname_int, '(I6)')  number

    ! This is for one digit (no rounding)
    !write (dirname, '(F4.1)')  number
  end function

  function dirname_real(number)
    real*8, intent(in) :: number
    character(len=8)  :: dirname_real
    write (dirname_real, '(F5.1)')  number
  end function

end module
