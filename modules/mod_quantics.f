!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module mod_quantics
!!
!! Contains the variables that required for 
!! any VCMaker runs
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


integer::n_qmode                                        ! Number of normal modes to be used
integer::n_qstate
integer,dimension(:),allocatable::qmode                 ! Index of the normal modes
integer,dimension(:),allocatable::qstate                ! Index of the electronic states
double precision,dimension(:),allocatable::qnrj         ! Vertical Energy of states
character(len=2),dimension(:),allocatable::qlabel       ! Role of the Normal mode


end module mod_quantics
