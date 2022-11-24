        !! Conversion subroutine
        subroutine mk_convert(&
                             &val,&      
                             &unit&
                             &)

        use mod_constants
        implicit NONE
        double precision,intent(inout)::val   !! Value to convert
        character(len=*),intent(in)::unit     !! Initial energy unit

        CALL upper(unit)
        if (trim(adjustl(unit)).eq.'EV') then
            val=val                                      ! Nothing is happening. Why is this line here?
        else if (trim(adjustl(unit)).eq.'EH') then
            val=val*hartoeV                              ! [Eh] -> [eV]
        else if (trim(adjustl(unit)).eq.'CM-1') then 
            val=val*hjs*(hartoeV/hartoJ)*c               ! E= hc nu  | = h[J.s -> Eh.s -> eV.s] c[cm.s-1] nu[cm-1]
        else if (trim(adjustl(unit)).eq.'NM') then 
            val=(hjs*(hartoeV/hartoJ)*c*1E-2)/(val*1E-9) ! E= hc / lambda | = h[J.s -> Eh.s -> eV.s] c[cm.s-1 -> m.s-1] / lambda[nm -> m]
        else if (trim(adjustl(unit)).eq.'') then
            val=val                                      !! No unit detected, assuming eV !
        end if
        end subroutine