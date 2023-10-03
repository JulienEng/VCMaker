        SUBROUTINE header()
        IMPLICIT NONE
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone
        integer,dimension(8) :: values


        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)

        Write(*,*)
        Write(*,*)
        write(*,*) '  ██╗   ██╗ ██████╗███╗   ███╗'&
                &,' █████╗ ██╗  ██╗███████╗██████╗'
        write(*,*) '  ██║   ██║██╔════╝████╗ ████║'&
                &,'██╔══██╗██║ ██╔╝██╔════╝██╔══██╗'
        write(*,*) '  ██║   ██║██║     ██╔████╔██║'&
                &,'███████║█████╔╝ █████╗  ██████╔╝'
        write(*,*) '  ╚██╗ ██╔╝██║     ██║╚██╔╝██║'&
                &,'██╔══██║██╔═██╗ ██╔══╝  ██╔══██╗'
        write(*,*) '   ╚████╔╝ ╚██████╗██║ ╚═╝ ██║'&
                &,'██║  ██║██║  ██╗███████╗██║  ██║'
        write(*,*) '    ╚═══╝   ╚═════╝╚═╝     ╚═╝'&
                &,'╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝'
        write(*,*) '                                        v1.0          '
        Write(*,*)
        Write(*,*) '+=============================================================+'
        write(*,*) '|  J. Eng*                                                    |'             
        write(*,*) '|  * Chemistry, School of Natural and Environmental Sciences, |'
        write(*,*) '|    Newcastle University, NE1 7RU, Newcastle Upon Tyne, UK   |'
        write(*,*) '|                                                             |' 
        write(*,*) '| VCMaker website: https://vcmaker.glitch.me/                 |' 
        write(*,*) '|                                                             |' 
        write(*,*) '| Cite me as (bibtex entry at the end):                       |'
        write(*,*) '+=============================================================+'
        write(*,*)
        write(*,*) 'Started on the ', date(7:8),'/',date(5:6),'/',date(1:4),' at ', time(1:2),':',time(3:4),':',time(5:6),'.'
        write(*,*)


        END SUBROUTINE

        SUBROUTINE endoftimes()
        IMPLICIT NONE
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone
        integer,dimension(8) :: values

        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)

        write(*,*)
        write(*,*) 'Stopped on the ', date(7:8),'/',date(5:6),'/',date(1:4),' at ', time(1:2),':',time(3:4),':',time(5:6),'.'
        write(*,*) ' Thank you for using VCMaker!'
        END SUBROUTINE
