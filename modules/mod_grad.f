!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module mod_grad                                          !!!!!!
!!                                                          !!!!!! 
!! Variables requires if log_grad=.TRUE.                    !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::n_grad                                         !#of gradient files
character(len=666),dimension(:),allocatable::file_grad  ! Files containing the excited-states gradients
double precision,dimension(:,:),allocatable::grad       ! Excited-states gradients
double precision,dimension(:),allocatable::grad_tmp     ! Temporary array to call the subroutine
double precision,dimension(:),allocatable::kgrad_tmp    ! Temporary array to call the subroutine
double precision,dimension(:,:),allocatable::kappa_grad ! Intrastate coupling
double precision,dimension(:),allocatable::eshift_grad  ! dE for each DNC.
double precision::degrad                                ! Energy shift summed over all normal modes
end module mod_grad
