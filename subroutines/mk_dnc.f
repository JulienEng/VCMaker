      subroutine mk_dnc(dnc_dnc,dnc_tdnc,dnc_freq,dnc_evec,dnc_mass,dnc_n)
      use mod_constants
      implicit none
      integer::k,l                                                               !! Loop integers
      integer::dnc_n                                                             !! Number of atoms
      double precision,dimension(3*dnc_n,3*dnc_n),intent(IN):: dnc_evec          !! Eigenvectors of the Cartesian Hessian
      double precision,dimension(3*dnc_n),intent(IN)::         dnc_freq          !! Eigenvalues of the Cartesian Hessian (w^2)
      double precision,dimension(dnc_n),intent(IN)::           dnc_mass          !! Mass of atoms
      double precision,dimension(3*dnc_n,3*dnc_n),intent(OUT):: dnc_dnc          !! DNC matrice
      double precision,dimension(3*dnc_n,3*dnc_n),intent(OUT)::dnc_tdnc          !! Transpose DNC Matrice
      
      !! Loop over the columns
      do k=1,3*dnc_n 
      !! Loop over the rows
         do l=1,3*dnc_n
      !! Weighting the eigenvector of the hessian by sqrt(freq*mass)
            dnc_dnc(k,l)=dnc_evec(k,l)*dsqrt( dsqrt(dnc_freq(k))*&
            &dsqrt(me/amu)*dnc_mass(FLOOR((DBLE(l)+2)/3.0d0))*(amu/me))

      !! Weighting the eigenvector of the hessian by 1/sqrt(freq*mass)      
            dnc_tdnc(k,l)=dnc_evec(k,l)/dsqrt(dsqrt(dnc_freq(k))&
            &*dsqrt(me/amu)*dnc_mass(FLOOR((DBLE(l)+2)/3.0d0))*(amu/me))
        end do
      end do

      end subroutine

