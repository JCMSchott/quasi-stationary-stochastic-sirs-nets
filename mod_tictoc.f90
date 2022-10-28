module mod_tictoc
implicit none

    type :: tictoc
        real*8 :: t_ini, t_fin, t_tot
        contains
            procedure :: reset => t_reset
            procedure :: start => t_reset
            procedure :: tic => t_tik
            procedure :: toc => t_tok
    end type

contains

    subroutine t_tik(this)
        class(tictoc) :: this
        call cpu_time(this%t_ini)
    end subroutine
    
    subroutine t_tok(this)
        class(tictoc) :: this
        call cpu_time(this%t_fin)
        this%t_tot = this%t_tot + (this%t_fin - this%t_ini )
    end subroutine

    subroutine t_reset(this)
        class(tictoc) :: this
        this%t_tot = 0d0
    end subroutine

end module
