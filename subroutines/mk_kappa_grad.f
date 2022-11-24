!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mk_kappa_grad(&
                &kg_freq,&  ! Normal mode frequencies
                &kg_grad,&  ! Input Cartesian Gradient
                &kg_kappa,& ! Output DNC Gradient
                &kg_tdnc,&  ! Transformation matrice
                &kg_de,&    ! Total Energy Shift
                &kg_n)      ! Number of Atoms
use mod_constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer::l
integer,intent(IN)::kg_n
double precision,dimension(3*kg_n),intent(in)::kg_grad
double precision,dimension(3*kg_n),intent(in)::kg_freq
double precision,dimension(3*kg_n,3*kg_n),intent(in)::kg_tdnc

double precision,intent(inout)::kg_de
double precision,dimension(3*kg_n),intent(inout)::kg_kappa

   kg_kappa=MATMUL(kg_tdnc,kg_grad)
   do l=1,3*kg_n
   kg_kappa(l)=kg_kappa(l)*hartoeV
   if (l.gt.6) then
   kg_de=kg_de-(kg_kappa(l)/hartoeV)**2/(2*dsqrt(kg_freq(l)*me/amu))
   end if
   end do

end subroutine

