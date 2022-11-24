      subroutine mk_lambda(lb_hess_a,&
                          &lb_hess_b,&
                          &lb_nrj_a,&
                          &lb_nrj_b,&
                          &lb_tdnc,&
                          &lb_lambda,&
                          &lb_n)
      use mod_constants
      Implicit none
      integer::k                                                        !! Loop Integer
      integer,intent(in)::lb_n                                          !! Number of atom
      double precision,dimension(3*lb_n,3*lb_n),intent(in)::lb_hess_a   !! Cartesian Hessian of state a
      double precision,dimension(3*lb_n,3*lb_n),intent(in)::lb_hess_b   !! Cartesian Hessian of state b
      double precision,intent(in)::lb_nrj_a                             !! Energy of state a
      double precision,intent(in)::lb_nrj_b                             !! Energy of state b
      double precision,dimension(3*lb_n,3*lb_n),intent(in)::lb_tdnc     !! DNC matrice
      double precision,dimension(3*lb_n),intent(out)::lb_lambda         !! Lambda output
      double precision,dimension(3*lb_n,3*lb_n)::lb_mhess_a             !! DNC Hessian of state a
      double precision,dimension(3*lb_n,3*lb_n)::lb_mhess_b             !! DNC Hessian of state b


      !! Projection of the excited-states Hessian on to the ground state
      !! Normal modes. dE/(dxi)(dxj) -> dE/(dQk)(dQl)
      lb_mhess_a=MATMUL(MATMUL(lb_tdnc,lb_hess_a),TRANSPOSE(lb_tdnc))
      lb_mhess_b=MATMUL(MATMUL(lb_tdnc,lb_hess_b),TRANSPOSE(lb_tdnc))

      !! The first 6 lambda are set to zero. They correspond to the coupling
      !! Along the global rotations and translations. The difference of 
      !! hessian should be null along those dimensions (== flat potential)
      do k=1,6
         lb_lambda(k)=0.0d0
      end do

      !! A negative Lambda corresponds to an imaginary square root.
      !! A negative Lambda translate to two states "attracting" each other.
      do k=7,3*lb_n
         lb_lambda(k)=SIGN(0.5d0*dsqrt(ABS((lb_nrj_a-lb_nrj_b)*(lb_mhess_a(k,k)-lb_mhess_b(k,k)))), &
           & (lb_nrj_a-lb_nrj_b)*(lb_mhess_a(k,k)-lb_mhess_b(k,k)) )


      !! A quick check. State b needs to be higher in energy that state a.
      !! If not, lambda is the opposite.
         if (lb_nrj_a.gt.lb_nrj_b) then ! State A is the upper state
             lb_lambda(k)=-lb_lambda(k)*hartoev  ! Sign you be the opposite.
         else
             lb_lambda(k)=+lb_lambda(k)*hartoev
         end if
      end do

      end subroutine
