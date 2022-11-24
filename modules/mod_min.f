!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module mod_min
!!
!! Contains the variables that required for 
!! the generation of the "Theoretical" minimum
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer::min_n
    integer,dimension(:),allocatable::min_mode

    double precision,dimension(:,:),allocatable::min_outcoord

end module mod_min
