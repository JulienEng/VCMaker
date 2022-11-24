        SUBROUTINE out_xyz(oxyz_coords,oxyz_labels,oxyz_n,oxyz_file)
        IMPLICIT NONE
        integer::i,j
        integer,intent(IN)::oxyz_n
        double precision,dimension(3,oxyz_n),intent(IN)::oxyz_coords
        character(len=2),dimension(oxyz_n),intent(IN)::oxyz_labels
        character(len=*),intent(IN)::oxyz_file
        integer::oxyz_unit

        if (oxyz_file=="screen") then
        oxyz_unit=6
        else
        oxyz_unit=68
        open(oxyz_unit,file=oxyz_file)
        end if

        write(oxyz_unit,'(i0)') oxyz_n
        write(oxyz_unit,*)
        do i=1,oxyz_n
        write(oxyz_unit,'(A2,3F18.10)') oxyz_labels(i),(oxyz_coords(j,i),j=1,3)
        end do

        if (oxyz_unit.eq.68) then
        close(oxyz_unit)
        end if
        END SUBROUTINE

