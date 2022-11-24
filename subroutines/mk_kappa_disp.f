!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mk_kappa_disp(&
                      &kd_freq,&  ! Normal mode frequencies
                      &kd_geo,&   ! Input Cartesian Gradient
                     &kd_ref,&    ! Reference Cartesian Coord
                     &kd_kappa,&  ! Output DNC Gradient
                     &kd_dnc,&    ! Transformation matrice
                     &kd_de,&     ! Total Energy Shift
                     &kd_n)       ! Number of Atoms
      use mod_constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      integer::l
      integer,intent(IN)::kd_n
      double precision,dimension(3,kd_n),intent(in)::kd_geo
      double precision,dimension(3,kd_n),intent(in)::kd_ref
      double precision,dimension(3*kd_n),intent(in)::kd_freq
      double precision,dimension(3*kd_n,3*kd_n),intent(in)::kd_dnc

      double precision,dimension(3*kd_n)::kd_3n_geo
      double precision,dimension(3*kd_n)::kd_3n_ref
      double precision,intent(out)::kd_de
      double precision,dimension(3*kd_n),intent(out)::kd_kappa

      !! Convert the Cartesian coordinates array to a 3*kd_n vector
      call xyz_nto3n(kd_ref,kd_3n_ref,kd_n,1)   
      call xyz_nto3n(kd_geo,kd_3n_geo,kd_n,1) 

      !! Projection of the Cartesian geoemtry difference (in Bohr) on the DNCs.
      kd_kappa=MATMUL(kd_dnc,(kd_3n_geo(:)-kd_3n_ref(:))*angtoa0)
      !! Loop over the number of atoms
      do l=1,3*kd_n
         !! Conversion dQ to Kappa
         kd_kappa(l)=-kd_kappa(l)*dsqrt(kd_freq(l)*me/amu)*hartoeV
         !! if l >6: DNC is a vibration. We compute the associated shift in energy.
         if (l.gt.6) then
            kd_de=kd_de-(kd_kappa(l)/hartoeV)**2/(2*dsqrt(kd_freq(l)*me/amu)) ! In Hartree
         end if
      end do

end subroutine

