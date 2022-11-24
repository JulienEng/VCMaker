!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module mod_lambda                                         !!
!!                                                          !! 
!! Variables requires if log_lambda=.TRUE.                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::n_hess                                              !#of exc. state Hessian files
character(len=666),dimension(:),allocatable::file_eshess     ! Files containing the excited-states hessian matrices
double precision,dimension(:,:,:),allocatable::hess_es       ! Excited states Hessian matrices
double precision,dimension(:),allocatable::freq_es       ! Excited states Hessian matrices
double precision,dimension(:,:),allocatable::hess_es_tmp     ! Excited states Hessian matrices -- Temporary
double precision,dimension(:,:),allocatable::scratch_a       ! Those are matrices we discard. May be useful later! 
double precision,dimension(:,:),allocatable::scratch_b       ! Those are matrices we discard. May be useful later! 
integer::n_lambda                                            ! Number of interstate coupling to compute
integer,dimension(:,:),allocatable::lambda_pairs             ! Pairs of states to compute the interstate coupling
double precision,dimension(:,:),allocatable::lambda          ! Interstate coupling
double precision,dimension(:,:),allocatable::lambda_up          ! Interstate coupling

double precision,dimension(:,:),allocatable::lambda_down          ! Interstate coupling

double precision,dimension(:),allocatable::vib_nrj           ! Vertical Energies


end module mod_lambda
