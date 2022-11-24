        subroutine out_sixcol(sixcol_data,sixcol_n)
        implicit none
        integer,intent(in)::sixcol_n
        double precision,dimension(3*sixcol_n,3*sixcol_n),intent(in) :: sixcol_data
        integer::k,rest,o,p,m

        !! Computing the number of blocks for the output formatting
        k=INT(3*sixcol_n/6.0d0) 
        !! Are there modes remaining ?
        rest=3*sixcol_n-k*6
        !! Loop over the number of blocks
        do o=1,k-1
        !! Write block header
        write(*,*) 'Mode: ',(p,p=(o-1)*6+1,(o-1)*6+6)
        !! Loop over the number of lines
        do m=1,sixcol_n
        write(*,'(A9,6F12.4)') '    | x  ',  (sixcol_data(p,(m-1)*3+1),p=(o)*6+1,(o)*6+6)
        write(*,'(I3,A6,6F12.4)') m,' | y  ',(sixcol_data(p,(m-1)*3+2),p=(o)*6+1,(o)*6+6)
        write(*,'(A9,6F12.4)') '    | z  ',  (sixcol_data(p,(m-1)*3+3),p=(o)*6+1,(o)*6+6)
        end do
        write(*,*)
        end do
        !! If there are any remaining modes to print!
        if (rest.ne.0) then
        !! Write block header
        write(*,*) 'Mode: ',(p,p=(k)*6+1,(k)*6+rest)
        !! Loop over the number of lines
        do m=1,sixcol_n
        write(*,'(A9,6F12.4)') '    | x  ',   (sixcol_data(p,(m-1)*3+1),p=(k)*6+1,k*6+rest)
        write(*,'(I3,A6,6F12.4)') m,' | y  ', (sixcol_data(p,(m-1)*3+2),p=(k)*6+1,k*6+rest)
        write(*,'(A9,6F12.4)') '    | z  ',   (sixcol_data(p,(m-1)*3+3),p=(k)*6+1,k*6+rest)
        end do
        end if

end subroutine
