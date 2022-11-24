        subroutine get_input(filen)
        use mod_input
        use mod_constants
        use mod_grad
        use mod_disp
        use mod_gap
        use mod_lambda
        use mod_dncscan
        use mod_quantics
        use mod_min
        use mod_err
        !use mod_vsa
        implicit none
        integer::i
        character(*)::filen
        character(len=666)::reader,reader_a
        integer::iost
        character(len=60)::string



        !! Searching for the file that contains the reference geometry.
        open(67,file=filen)
        string="GEOMETRY"
        !! Initialisation of the string
        write(reader,*) 'Lemmy'
        do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                read(67,*,iostat=iost) reader
                !! Convert everything to an all UPPERCASE string
                call upper(reader)
                !! Check for the end-of-file
                if (iost.ne.0) then
                        !! Write the error
                        write(err_msg,*) 'Keyword "Geometry" not found.'
                        call mk_error() 
                end if
        end do
        !! Keyword "Geometry" found. We go back and read the whole line.
        BACKSPACE(67)
        !! Reading the filename.
        read(67,*) reader, file_gs
        close(67)

        !! Searching for the file that contains the Cartesian Hessian
        !! at the reference geometry.
        open(67,file=filen)
        string="HESSIAN"
        !! Initialisation of the string
        write(reader,*) 'Lemmy'
        do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                read(67,*,iostat=iost) reader
                !! Convert everything to an all UPPERCASE string
                call upper(reader)
                !! Check for the end-of-file
                if (iost.ne.0) then
                        !! Write the error
                        write(err_msg,*) 'Keyword "Hessian" not found.'
                        call mk_error() 
                end if
        end do
        !! Keyword "Hessian" found. We go back and read the whole line.
        BACKSPACE(67)
        !! Reading the filename.
        read(67,*) reader, hessian
        close(67)

        !! Retrieve job specifications
        !! Estimation of k from the Cartesian Gradient:
        string="lvc_grad"
        CALL get_logical(filen,log_grad,string)
        if (log_grad) then
                !! Generation of the theoretical minima (Based on the values of k)
                string="gen_min"
                call get_logical(filen,log_min,string)
        end if
        !! Estimation of k from the displaced geometry
        string="lvc_disp"
        CALL get_logical(filen,log_disp,string)
        if (log_grad.AND.log_disp) then
                string="lvc_gap"
                CALL get_logical(filen,log_gap,string)
        end if
        !! Estimation of the interstate vibronic coupling l
        string="lvc_lambda"
        CALL get_logical(filen,log_lambda,string)
        !! Generating a scan along one or many DNC.
        string="dnc_scan"
        CALL get_logical(filen,log_dncscan,string)
        if (log_dncscan) then
                !! Generating a 2D scan
                string="dnc_diag"
                CALL get_logical(filen,log_dncdiag,string)
        end if
        !! Generating a N-Dimension Grid of points
        string="dnc_grid"
        CALL get_logical(filen,log_dncgrid,string)
        !! Generate a single distorted geometry.
        string="gen_xyz"
        CALL get_logical(filen,log_genxyz,string)
        !! Generate an operator file compatible with Quantics.
        string='out_quantics'
        CALL get_logical(filen,log_quantics,string)

        !! Get For Gradient Kappas specifications
        if (log_grad) then
                open(67,file=filen)
                string="GRAD"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Estimation of K from the Cartesian gradient',&
                                                &' requested but no "Grad" Keyword has been found!'
                                call mk_error() 
                        end if
                end do
                BACKSPACE(67)
                read(67,*) reader, n_grad

                allocate(file_grad(n_grad))
                do i=1,n_grad
                        read(67,*) file_grad(i)
                end do
                close(67)
     
                if (log_min) then
                        open(67,file=filen)
                        string="MIN"
                        write(reader,*) 'Lemmy'
                        do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                                read(67,*,iostat=iost) reader
                                call upper(reader)
                                if (iost.ne.0) then
                                        write(err_msg,*) 'MIN block not found!'
                                        call mk_error()
                                end if  
                        end do
                        BACKSPACE(67)
                        read(67,*) reader, reader_a
                        call upper(reader_a)
                        if ( TRIM(ADJUSTL(reader_a)).eq."ALL") then
                                min_n=-1
                        else
                                read(reader_a,*) min_n
                                allocate(min_mode(min_n))
                                read(67,*) (min_mode(i),i=1,min_n)
                        end if
                        close(67)
                end if
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get For Displacement Kappas specifications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (log_disp) then
                open(67,file=filen)
                string="DISP"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Estimation of K from displaced geometries',&
                                                &' requested but no "Disp" Keyword has been found!'
                                call mk_error() 
                        end if
                end do
                BACKSPACE(67)
                read(67,*) reader, n_geo

                allocate(file_geo(n_geo))
                do i=1,n_geo
                        read(67,*) file_geo(i)
                end do
                close(67)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get For Specifications for the F-Analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (log_gap) then
                open(67,file=filen)
                string="GAP"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Floppiness analysis requested',&
                                                &' but no "Gap" Keyword has been found!'
                                call mk_error() 
                        end if
                end do
                BACKSPACE(67)
                read(67,*) reader, n_gap

                allocate(gap_disp(n_gap))
                allocate(gap_grad(n_gap))
                do i=1,n_gap
                        read(67,*,iostat=iost) gap_grad(i), gap_disp(i)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Error in the the GAP specifications'
                                call mk_error() 
                        end if
                end do
                close(67)       
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get For Lambda job specifications !!!! Think about all the possible way of mess that up!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (log_lambda) then 
                open(67,file=filen)
                string="LAMBDA"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Estimation of Lambda requested but no "Lambda" Keyword has been found!'
                                call mk_error() 
                        end if
                end do
                BACKSPACE(67)
                read(67,*,iostat=iost) reader, n_lambda ! n_lambda is the number of lambda to compute
                if (iost.ne.0) then
                        write(err_msg,*) 'Specified number of interstate vibronic coupling is not valid'
                        call mk_error() 
                end if
                allocate(lambda_pairs(n_lambda,2))
                do i=1,n_lambda
                        read(67,*) lambda_pairs(i,1),lambda_pairs(i,2) ! n_pairs define between which states
                end do
                REWIND 67
                string="ES_HESSIAN"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Estimation of Lambda requested but no "ES_Hessian" Keyword has been found!'
                                call mk_error() 
                        end if
                end do
                BACKSPACE(67)
                read(67,*,iostat=iost) reader, n_hess ! number of Hessian file
                if (iost.gt.0) then 
                        write(err_msg,*) 'Specified number of excited state hessian is not valid!'
                        call mk_error() 
                end if
                if (n_hess.lt.MAXVAL(lambda_pairs)) then
                        write(err_msg,*) 'The number of excited states considered is smaller than the highest states requested',&
                                        &' for the vibronic coupling calculations'
                        call mk_error() 
                end if
                allocate(file_eshess(n_hess))
                allocate(vib_nrj(n_hess)) 
                do i=1,n_hess
                        write(reader,*) ''

                        read(67,*,iostat=iost) file_eshess(i), vib_nrj(i)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Bad energy specified for file: ', trim(adjustl(file_eshess(i))),'!'
                                call mk_error() 
                        end if
                        BACKSPACE 67
                        read(67,*,iostat=iost) file_eshess(i), vib_nrj(i), reader
                        call upper(reader)    
                        call mk_convert(vib_nrj(i),reader)   ! Converting the unit to eV.
                end do
                close(67)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get For DNC Scan specifications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (log_dncscan) then

                open(67,file=filen)
                string="SCAN"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Scans along Normal Modes is requested but no "Scan" keyword has been found!'
                                call mk_error() 
                        end if
                end do
                BACKSPACE(67)
                read(67,*,iostat=iost) reader, reader_a 
                if (iost.ne.0) then
                        write(err_msg,*) 'Specified number of modes to scan along is invalid!'
                        call mk_error() 
                end if
                call upper(reader_a)
                if (reader_a.eq.'ALL') then
                        n_dnc=-1
                else
                        read(reader_a,*,iostat=iost) n_dnc

                        if (iost.ne.0) then
                                write(err_msg,*) 'Specified number of modes to scan along is invalid!'
                                call mk_error() 
                        end if
                end if

                if (n_dnc.ne.-1) then
                        allocate(dncscan_mode(n_dnc))
                        allocate(dncscan_ptleft(n_dnc))
                        allocate(dncscan_stepleft(n_dnc))
                        allocate(dncscan_ptright(n_dnc))
                        allocate(dncscan_stepright(n_dnc))
                        read(67,*,iostat=iost)      (dncscan_mode(i),i=1,n_dnc)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if
                        read(67,*,iostat=iost)    (dncscan_ptleft(i),i=1,n_dnc)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if

                        read(67,*,iostat=iost)  (dncscan_stepleft(i),i=1,n_dnc)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if
                        read(67,*,iostat=iost)   (dncscan_ptright(i),i=1,n_dnc)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if
                        read(67,*,iostat=iost) (dncscan_stepright(i),i=1,n_dnc)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if
                elseif (n_dnc.eq.-1) then
                        allocate(dncscan_ptleft(1))
                        allocate(dncscan_stepleft(1))
                        allocate(dncscan_ptright(1))
                        allocate(dncscan_stepright(1))
                        read(67,*,iostat=iost)    dncscan_ptleft(1)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if
                        read(67,*,iostat=iost)  dncscan_stepleft(1)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if
                        read(67,*,iostat=iost)   dncscan_ptright(1)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if
                        read(67,*,iostat=iost) dncscan_stepright(1)
                        if (iost.ne.0) then 
                                write(err_msg,*) 'Error in the DNC Scan specifications'
                                call mk_error()
                        end if
                end if
                close(67)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get For DNC Grid Specification
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (log_dncgrid) then
                open(67,file=filen)
                string="GRID"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Scans along a grid of Normal Modes',&
                                                &' is requested but no "Grid" keyword has been found!'
                                call mk_error() 
                        end if
                end do
                BACKSPACE 67
                read(67,*,iostat=iost) reader, n_dim 
                if (iost.ne.0) then
                        write(err_msg,*) 'Error in the dimension of the grid'
                        call mk_error()
                end if

                allocate(grd_mode(n_dim)) 
                allocate(grd_ptleft(n_dim)) 
                allocate(grd_ptright(n_dim)) 
                allocate(grd_stepleft(n_dim)) 
                allocate(grd_stepright(n_dim)) 

                read(67,*,iostat=iost)      (grd_mode(i),i=1,n_dim)
                if (iost.ne.0) then
                        write(err_msg,*) 'Error in the Grid specifications'
                        call mk_error()
                end if
                read(67,*,iostat=iost)    (grd_ptleft(i),i=1,n_dim)
                if (iost.ne.0) then
                        write(err_msg,*) 'Error in the Grid specifications'
                        call mk_error()
                end if
                read(67,*,iostat=iost)  (grd_stepleft(i),i=1,n_dim)
                if (iost.ne.0) then
                        write(err_msg,*) 'Error in the Grid specifications'
                        call mk_error()
                end if
                read(67,*,iostat=iost)   (grd_ptright(i),i=1,n_dim)
                if (iost.ne.0) then
                        write(err_msg,*) 'Error in the Grid specifications'
                        call mk_error()
                end if
                read(67,*,iostat=iost) (grd_stepright(i),i=1,n_dim)
                if (iost.ne.0) then
                        write(err_msg,*) 'Error in the Gris specifications'
                        call mk_error()
                end if

                close(67)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Check for output specification
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get For GEN_XYZ specification
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        if (log_genxyz) then
                open(67,file=filen)
                string="GEN"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Generation of a XYZ distorted structure requested',&
                                                &' by no "Gen" keyword has been found!'
                                call mk_error()
                        end if
                end do
                BACKSPACE 67
                read(67,*,iostat=iost) reader, n_gen
                allocate(gen_mode(n_gen))
                allocate(gen_disto(n_gen))
        
                do i=1,n_gen
                        read(67,*,iostat=iost) gen_mode(i)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Error in the Gen specifications'
                                call mk_error()
                        end if
                        BACKSPACE 67
                        read(67,*,iostat=iost) gen_mode(i),gen_disto(i)
                        if (iost.ne.0) then
                                write(err_msg,*) 'Error in the Gen specifications'
                                call mk_error()
                        end if
                end do
                close(67)
        end if

        if (log_quantics) then
                open(67,file=filen)
                string="QUANTICS_MODE"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)))
                        read(67,*) reader
                        call upper(reader)
                end do
                BACKSPACE(67)
                read(67,*) reader, reader_a 
                read(reader_a,*) n_qmode

                allocate(qmode(n_qmode))
                allocate(qlabel(n_qmode))

        ! Three scenarii! Let's check
                if ((log_grad.OR.log_disp).AND.(.NOT.log_lambda)) then !! is it XOR or  OR ?
                        ! We take the Kappas from log_grad by default.
                        do i=1,n_qmode
                                read(67,*) qmode(i)
                                write(qlabel(i),*) 'K'
                        end do
                else if (log_lambda.AND.(.NOT.(log_grad.OR.log_disp))) then
                        do i=1,n_qmode
                                read(67,*) qmode(i)
                                write(qlabel(i),*) 'L'
                        end do
                else if ((log_grad.OR.log_disp).AND.log_lambda) then
                        do i=1,n_qmode
                                read(67,*,iostat=iost) qmode(i), qlabel(i)
                        end do
                end if
                close(67)

                open(67,file=filen)
                string="QUANTICS_STATE"
                write(reader,*) 'Lemmy'
                do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)).and.iost.eq.0)
                        read(67,*,iostat=iost) reader
                        call upper(reader)
                end do
                BACKSPACE(67)
                read(67,*,iostat=iost) reader, reader_a 
                read(reader_a,*) n_qstate
                allocate(qstate(n_qstate))
                allocate(qnrj(n_qstate))

                do i=1,n_qstate
                        write(reader,*) ''
                        read(67,*,iostat=iost) qstate(i), qnrj(i), reader
                        call upper(reader)
                        call mk_convert(qnrj(i),reader)
                end do
                Close(67)
        end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine get_logical(filen,logic,string)     !!
        !! This subroutine look for 'string' in 'filen'!!
        !! if 'string' is found, then 'logic' is set   !!
        !! to true                                     !!
        use mod_input                                  !!
        use mod_constants                              !!
        implicit none                                  !!
        character(len=666)::reader                     !! 
        character(*)::string       !searched keyword   !!
        character(len=666)::str_value                  !! 
        character(*)::filen        !input file         !!                  
        integer::iost                                  !!
        logical::logic             !associated variable!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        open(67,file=filen)
        iost=0
        !! Resetting the reader
        write(reader,*) 'Lemmy'
        logic=.FALSE.
        call upper(string)
        do while (TRIM(ADJUSTL(reader)).ne.TRIM(ADJUSTL(string)).AND.iost.eq.0)
                read(67,*,iostat=iost) reader, str_value
                call upper(reader)
        end do
        call upper(str_value)
        if (TRIM(ADJUSTL(str_value)).eq.'TRUE') then
                logic=.TRUE.
        end if
        close(67)
end subroutine
