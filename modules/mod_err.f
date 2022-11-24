module mod_err

!! Contains all error flags that are checked during VCMaker


!Descriptive string
      character(len=999)::err_msg

!Format error
      logical::flag_fhessian     ! OK
      logical::flag_fgradient    ! OK
      logical::flag_fxyz         ! OK

!Structure mismatch
      logical::flag_xyzn
      logical::flag_atunknown
      logical::flag_xyzlabel

!Inversion errors?
      logical::flag_invert
!Internal coordinaes errors
      logical::flag_internal

      !! Generic errors, in each subroutine ?

end module mod_err
