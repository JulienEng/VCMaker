!! This subroutine makes errors explicit.
!! You should like this subroutine.
        SUBROUTINE mk_error()
        use mod_err
        IMPLICIT NONE 


        write(*,*) trim(adjustl(err_msg))
        call endoftimes()
        STOP



END SUBROUTINE
