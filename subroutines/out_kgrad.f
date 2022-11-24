 !! This subroutine prints the results of the LVC intrastate coupling.
 
        subroutine out_kgrad(out_freq,out_kappa,out_de,out_redmass,out_n)
        use mod_constants
        implicit none
        integer::k
        integer,intent(IN)::out_n
        double precision,dimension(3*out_n),intent(in)::out_freq
        double precision,dimension(3*out_n),intent(in)::out_redmass
        double precision,dimension(3*out_n),intent(in)::out_kappa
        double precision,intent(in)::out_de
        write(*,*) 'Mode     Freq / cm-1   Red.Mass | Kappa / meV    dE / meV       dQ '
        write(*,*) '---------------------------------------------------------------------------'
         !! Loop over the vibrations (i.e. NOT the global translations and rotations)
        !! The Mode index printed is adjusted so that the first vibration is labelled 1.
        do k=7,3*out_n
        write(*,'(I6,2F12.2,A4,F11.3,2F12.3)') (k-6), 1/(2*pi*c)*dsqrt(out_freq(k)*hartoJ/(bohrtom**2*amu)),&
            &out_redmass(k),'   |', out_kappa(k)*1000, -(out_kappa(k)/hartoeV)**2/(2*dsqrt(out_freq(k)*me/amu))*1000,&
           &-(out_kappa(k)/hartoeV)/(dsqrt(out_freq(k)*me/amu))
        end do
        write(*,*) '---------------------------------------------------------------------------'
        write(*,'(A21,6X,F30.3,A3)') 'Sum of Energy shifts:', (out_de*hartoeV),' eV'
        write(*,'(A21,6X,F30.3,A3)') 'Stokes shift:', (2*out_de*hartoeV),' eV'
        write(*,*)




end subroutine



