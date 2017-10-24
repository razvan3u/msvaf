program prog
   use functii
   use kinds
   implicit none
   character(len=80) :: argument
   real(dp) :: opt_e0, opt_f, opt_gamma


   if (command_argument_count() > 1) then
       call get_command_argument(1,argument)
       call citesteExperiment(trim(argument))
       call get_command_argument(2,argument)
       call citesteCalcul(trim(argument))
   else 
       print *, "Eroare, am nevoie de 2 nume de fisier"
       print *, "        ca argument: nume fisier exp,"
       print *, "        urmat de ca nume fisier calc."
       stop
   end if

   ! Verific daca se doreste preconditionarea cu E0-0
    if (command_argument_count() > 2) then
        call get_command_argument(3,argument)
        read(argument,*) opt_e0
        call setParams(set_e0=opt_e0)
    end if  

   ! Verific daca se doreste preconditionarea cu factor de scala
    if (command_argument_count() > 3) then
        call get_command_argument(4,argument)
        read(argument,*) opt_f
        call setParams(set_f0=opt_f)
    end if  
   ! Verific daca se doreste preconditionarea cu gamma
    if (command_argument_count() > 4) then
        call get_command_argument(5,argument)
        read(argument,*) opt_gamma
        call setParams(set_gamma0=opt_gamma)
    end if 
   call initializare
   call fit(400)

end program
