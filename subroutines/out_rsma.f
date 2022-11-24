 !! This subroutine prints the results of the RSMF analysis.
 
        subroutine out_rsma(out_freq,out_rrsma,out_redmass,out_n)
        use mod_constants
        implicit none
        integer::k
        integer,intent(IN)::out_n
        double precision,dimension(3*out_n),intent(in)::out_freq
        double precision,dimension(3*out_n),intent(in)::out_rrsma
        double precision,dimension(3*out_n),intent(in)::out_redmass

        !write(*,*) out_n
        write(*,*) 'Mode     Freq / cm-1   Red.Mass | RSMA         '
        write(*,*) '-----------------------------------------------'
         !! Loop over the vibrations (i.e. NOT the global translations and rotations)
        !! The Mode index printed is adjusted so that the first vibration is labelled 1.
        do k=7,3*out_n
        write(*,'(I6,2F12.2,A4,F11.3)') (k-6), 1/(2*pi*c)*dsqrt(out_freq(k)*hartoJ/(bohrtom**2*amu)),&
            &out_redmass(k),'   |', out_rrsma(k) 
        end do
        write(*,*) '---------------------------------------------------------------------------'
        write(*,*)




end subroutine



