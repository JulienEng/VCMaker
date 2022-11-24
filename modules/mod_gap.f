!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module mod_gap                                           !!!!!!
!!                                                          !!!!!! 
!! Variables required if log_gap=.TRUE.                     !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::n_gap                                          !#of gradient files
integer,dimension(:),allocatable::gap_grad              ! Index of the gradients
integer,dimension(:),allocatable::gap_disp              ! Index of the geomtries
double precision,dimension(:),allocatable::gap          ! GFI - size of n_gap
double precision,dimension(:,:),allocatable::rsma       ! RMSF - size of n_gi


end module mod_gap

