SUBROUTINE out_input()
         use mod_input
         use mod_constants
         use mod_grad
         use mod_gap
         use mod_disp
         use mod_dncscan
         use mod_min
         use mod_lambda
         IMPLICIT NONE
         integer::i,j,k,b,l,m ! Loop Integers
         double precision::scratch

Write(*,*) 'Job summary:'
Write(*,*) '-----------------------'
Write(*,*)
if (log_grad) then  
        write(*,'(t1,a,a,t22,a)') 'LVC_GRAD','................','True'
        Write(*,'(t5,a,i0)') 'Number of gradients: ',n_grad
        do i=1,n_grad
        write(*,'(t7,i0,a,t10,a)') i,'. ',TRIM(ADJUSTL(file_grad(i)))
        end do
        write(*,*)

        if (log_min) then
                write(*,'(t1,a,a,t22,a)') 'GEN_MIN','................','True'
                if (min_n.eq.-1) then
                        write(*,'(t7,a)') 'All modes are considered'
                        write(*,*)
                else
                        k=INT(min_n/6.0d0)
                        b=min_n-k*6 
                        if (k.gt.0) then
                                do l=1,k
                                        if (l.eq.1) then
                        write(*,'(a,6(i3,3X))') '            Modes considered:',     (min_mode((l-1)*6+m),m=1,6)
                                        else
                                                write(*,'(a,6(i3,3X))') '                               '&
                                                &,     (min_mode((l-1)*6+m),m=1,6)
                                        end if
                                end do
                        end if
                        if (b.gt.0) then
                                if (k.eq.0) then
                        write(*,'(a,6(i3,3X))') '            Modes considered:',     (min_mode(k*6+m),m=1,b)
                                else
                                        write(*,'(a,6(i3,3X))') '                             '&
                                        &,     (min_mode(k*6+m),m=1,b)
                                end if
                        end if
                        write(*,*)
                end if
        end if
else
        write(*,'(t1,a,a,t22,a)') 'LVC_GRAD','................','False'
end if
if (log_disp) then
        write(*,'(t1,a,a,t22,a)') 'LVC_DISP','................','True'
        Write(*,'(t5,a,i0)') 'Number of geometries: ',n_geo
        do i=1,n_geo
        write(*,'(t7,i0,a,t10,a)') i,'.',TRIM(ADJUSTL(file_geo(i)))
        end do
        write(*,*)
else
        write(*,'(t1,a,a,t22,a)') 'LVC_DISP','................','False'
end if
if (log_gap) then
        write(*,'(t1,a,a,t22,a)') 'LVC_GAP','..............','True'
        write(*,'(t5,a,i0)') 'Number of couples: ', n_gap
        do i=1,n_gap
        write(*,'(t7,i0,a,t10,a,a,a)') i,'.',TRIM(ADJUSTL(file_grad(gap_grad(i)))),&
                & ' and ',TRIM(ADJUSTL(file_geo(gap_disp(i))))
        end do
        write(*,*)
else
        write(*,'(t1,a,a,t22,a)') 'LVC_GAP','..............','False'
end if
if (log_lambda) then
        write(*,'(t1,a,a,t22,a)') 'LVC_LAMBDA','..............','True'
        write(*,'(t5,a,i0)') 'Number of couplings: ', n_lambda
        do i=1,n_lambda
        write(*,'(t7,i0,a,t10,a,a,a)') i,'.',TRIM(ADJUSTL(file_eshess(lambda_pairs(i,1)))),&
                & ' and ',TRIM(ADJUSTL(file_eshess(lambda_pairs(i,2))))
        end do
        write(*,*)
else
        write(*,'(t1,a,a,t22,a)') 'LVC_LAMBDA','..............','False'
end if
if (log_dncscan) then
        write(*,'(t1,a,a,t22,a)') 'DNC_SCAN','................','True'

        if (n_dnc.eq.-1) then

        scratch=DBLE(dncscan_ptleft(1))
        deallocate(dncscan_ptleft)
          allocate(dncscan_ptleft(3*n_at-6))
        do i=1,3*n_at-6
        dncscan_ptleft(i)=FLOOR(scratch)
        end do

        scratch=DBLE(dncscan_ptright(1))
        deallocate(dncscan_ptright)
          allocate(dncscan_ptright(3*n_at-6))
        do i=1,3*n_at-6
        dncscan_ptright(i)=FLOOR(scratch)
        end do

        scratch=dncscan_stepleft(1)
        deallocate(dncscan_stepleft)
          allocate(dncscan_stepleft(3*n_at-6))
        do i=1,3*n_at-6
        dncscan_stepleft(i)=scratch
        end do

        scratch=dncscan_stepright(1)
        deallocate(dncscan_stepright)
          allocate(dncscan_stepright(3*n_at-6))
        do i=1,3*n_at-6
        dncscan_stepright(i)=scratch
        end do

        n_dnc=3*n_at-6
        allocate(dncscan_mode(n_dnc))
        do i=1,3*n_at-6
        dncscan_mode(i)=i
        end do

        write(*,'(a)')      '     Number of modes to scan: All'
        write(*,'(a,i3)')   '    Number of negative steps:', dncscan_ptleft(1)
        write(*,'(a,i3)')   '    Number of positive steps:', dncscan_ptright(1)
        write(*,'(a,f6.2)') '          Negative increment:', dncscan_stepleft(1)
        write(*,'(a,f6.2)') '          Positive increment:', dncscan_stepright(1)
        

        elseif(n_dnc.gt.0) then
        write(*,'(a,i3)') '     Number of modes to scan:', n_dnc

        k=INT(n_dnc/6.0d0)
        b=n_dnc-k*6
        write(*,*)
        if (k.gt.0) then
                do l=1,k
        write(*,'(a,6(i3,3X))') '               Modes to scan:',     (dncscan_mode((l-1)*6+m),m=1,6)
        write(*,'(a,6(i3,3X))') '    Number of negative steps:',   (dncscan_ptleft((l-1)*6+m),m=1,6) 
        write(*,'(a,6(i3,3X))') '    Number of positive steps:',  (dncscan_ptright((l-1)*6+m),m=1,6)
        write(*,'(a,6f6.2)')    '          Negative increment:',    (dncscan_stepleft((l-1)*6+m),m=1,6)
        write(*,'(a,6f6.2)')    '          Positive increment:',   (dncscan_stepright((l-1)*6+m),m=1,6)
        write(*,*) ' '
                end do
        end if
        if (b.gt.0) then
        write(*,'(a,6(i3,3X))') '               Modes to scan:',     (dncscan_mode(k*6+m),m=1,b)
        write(*,'(a,6(i3,3X))') '    Number of negative steps:',   (dncscan_ptleft(k*6+m),m=1,b)
        write(*,'(a,6(i3,3X))') '    Number of positive steps:',  (dncscan_ptright(k*6+m),m=1,b)
        write(*,'(a,6f6.2)')    '          Negative increment:',    (dncscan_stepleft(k*6+m),m=1,b)
        write(*,'(a,6f6.2)')    '          Positive increment:',   (dncscan_stepright(k*6+m),m=1,b)
        write(*,*) ' '
        end if
 
        if (log_dncdiag) then
                write(*,'(t1,a,a,t22,a)') '     DNC_DIAG','................','TRUE'
        end if
        elseif(n_dnc.gt.0) then
                write(*,'(a)') 'No scan requested. Switching it off'
        end if
else
        write(*,'(t1,a,a,t22,a)') 'DNC_SCAN','................','False'
end if
if (log_dncgrid) then
        write(*,'(t1,a,a,t22,a)') 'DNC_GRID','................','TRUE'
        write(*,'(a,i3)') '     Number of modes to scan:', n_dim

        k=INT(n_dim/6.0d0)
        b=n_dim-k*6
        write(*,*)
        if (k.gt.0) then
                do l=1,k
        write(*,'(a,6(i3,3X))') '            Modes considered:',     (grd_mode((l-1)*6+m),m=1,6)
        write(*,'(a,6(i3,3X))') '    Number of negative steps:',   (grd_ptleft((l-1)*6+m),m=1,6) 
        write(*,'(a,6(i3,3X))') '    Number of positive steps:',  (grd_ptright((l-1)*6+m),m=1,6)
        write(*,'(a,6f6.2)')    '          Negative increment:',    (grd_stepleft((l-1)*6+m),m=1,6)
        write(*,'(a,6f6.2)')    '          Positive increment:',   (grd_stepright((l-1)*6+m),m=1,6)
        write(*,*) ' '
                end do
        end if
        if (b.gt.0) then
        write(*,'(a,6(i3,3X))') '            Modes considered:',     (grd_mode(k*6+m),m=1,b)
        write(*,'(a,6(i3,3X))') '    Number of negative steps:',   (grd_ptleft(k*6+m),m=1,b)
        write(*,'(a,6(i3,3X))') '    Number of positive steps:',  (grd_ptright(k*6+m),m=1,b)
        write(*,'(a,6f6.2)')    '          Negative increment:',    (grd_stepleft(k*6+m),m=1,b)
        write(*,'(a,6f6.2)')    '          Positive increment:',   (grd_stepright(k*6+m),m=1,b)
        write(*,*)
        end if
        
else
        write(*,'(t1,a,a,t22,a)') 'DNC_GRID','................','False'
end if
if (log_quantics) then
        write(*,'(t1,a,a,t22,a)') 'OUT_QUANTICS','..............','True'
else
        write(*,'(t1,a,a,t22,a)') 'OUT_QUANTICS','..............','False'
end if
write(*,*)
write(*,*) '================================='
Write(*,*) '-- Coordinates'
write(*,*) '================================='
write(*,*)
write(*,*) 'At        x           y           z         Mass / amu'
do i=1,n_at
write(*,'(A2,3F12.6,A5,F12.6)') label(i), (coord_gs(j,i),j=1,3),'  |  ', mass(i)!*(me/amu)
end do
write(*,*) '----------------------------------------------------------------------------'
write(*,'(A3,F11.6,2F12.6)') 'com', (com(j),j=1,3)
write(*,*)
write(*,*) 'Axis of Intertia --- '
write(*,'(A30)') '=============================='
write(*,'(A2,3F9.4)') 'x ', (inertia_axis(j,1),j=1,3)
write(*,'(A2,3F9.4)') 'y ', (inertia_axis(j,2),j=1,3)
write(*,'(A2,3F9.4)') 'z ', (inertia_axis(j,3),j=1,3)

        write(*,*)
        write(*,*)




END SUBROUTINE

