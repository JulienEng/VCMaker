    Subroutine mk_min(min_coord&      ! Initial Geometry
        &,min_label&      ! Atom labels
        &,min_tdnc&       ! DNC displacements
        &,min_grad&       ! Gradients
        &,min_freq&       ! Frequencies
        &,min_nat&
        &,min_out)          ! #Atom

        use mod_constants
        use mod_min
        implicit NONE
        integer::i
        integer,intent(in)::min_nat
        character(len=2),dimension(min_nat)::min_label
        double precision,dimension(3,min_nat),intent(inout)::min_coord
        double precision,dimension(3*min_nat,3*min_nat),intent(inout)::min_tdnc
        double precision,dimension(3*min_nat),intent(in)::min_grad
        double precision,dimension(3*min_nat),intent(in)::min_freq
        double precision,dimension(min_n)::min_dq
        character(len=*)::min_out

        logical::cleangen

        allocate(min_outcoord(3,min_nat))
        cleangen=.TRUE.

        min_outcoord=min_coord
 
        do i=1,min_n
            min_dq(i)=-min_grad(min_mode(i)+6)/(dsqrt(min_freq(min_mode(i)+6)*me/amu)*hartoeV)
        end do

        CALL mk_genxyz(min_coord&
            &,min_label&
            &,min_tdnc&
            &,min_mode&
            &,min_dq&
            &,min_n&
            &,min_nat&
            &,TRIM(ADJUSTL(min_out)))

            deallocate(min_outcoord)

    END SUBROUTINE
