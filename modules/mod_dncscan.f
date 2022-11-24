!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module mod_dncscan                                        !!
!!                                                          !! 
!! Variables requires if log_dncscan=.TRUE.                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::n_dnc                                               !#of mode to scan along
integer,dimension(:),allocatable::dncscan_mode
         integer,dimension(:),allocatable::dncscan_ptleft
         integer,dimension(:),allocatable::dncscan_ptright
double precision,dimension(:),allocatable::dncscan_stepleft
double precision,dimension(:),allocatable::dncscan_stepright

!!! Grid specification
integer::n_dim                                            
integer,dimension(:),allocatable::grd_mode
         integer,dimension(:),allocatable::grd_ptleft
         integer,dimension(:),allocatable::grd_ptright
double precision,dimension(:),allocatable::grd_stepleft
double precision,dimension(:),allocatable::grd_stepright

!!! Gen specification
integer::n_gen
integer,dimension(:),allocatable::gen_mode
double precision,dimension(:),allocatable::gen_disto
         
end module mod_dncscan
