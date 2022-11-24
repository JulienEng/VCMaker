!! This subroutine prints the results of the LVC interstate coupling.

        subroutine out_lambda(out_freq,out_lmb,out_lmbp,out_lmbm,out_redmass,out_n)
        use mod_constants
        implicit none
        integer::k
        integer,intent(IN)::out_n
        double precision,dimension(3*out_n),intent(in)::out_freq
        double precision,dimension(3*out_n),intent(in)::out_redmass
        double precision,dimension(3*out_n),intent(in)::out_lmb
        double precision,dimension(3*out_n),intent(in)::out_lmbp
        double precision,dimension(3*out_n),intent(in)::out_lmbm


        write(*,*) 'Mode     Freq / cm-1   Red.Mass | Lambda / meV | Lambda+ / meV | Lambda- / meV    '
        write(*,*) '----------------------------------------------------------------------------------'
        !! Loop over the vibrations (i.e. NOT the global translations and rotations)
        !! The Mode index printed is adjusted so that the first vibration is labelled 1.
        do k=7,3*out_n
            write(*,'(I6,2F12.2,A4,F11.3,4X,2(F11.3,5X))') (k-6), 1/(2*pi*c)*dsqrt(out_freq(k)*hartoJ/(bohrtom**2*amu)),&
           &out_redmass(k),'   |', out_lmb(k)*1000 , out_lmbp(k)*1000 , out_lmbm(k)*1000
        end do
        write(*,*) '----------------------------------------------------------------------------------'
        write(*,*)
        end subroutine
