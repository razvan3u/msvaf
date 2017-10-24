module functii

!------------------------------------------------------------------------------!
! Modul de functii si subrutine pentru calculul de regresie                    !
! neliniara !                                                                  !
!------------------------------------------------------------------------------!
! Copyright (c) 2015-2016, Razvan Cirdei <razvan.cirdei@chem.uaic.ro>          !
!                          Dan Maftei <dan.maftei@chem.uaic.ro>                !
!------------------------------------------------------------------------------!

    use kinds
    use jacobi

    implicit none
    integer :: nexp, nvib 

    real(dp) :: f, gama,  e0  
    real(dp) :: lambda, oldWSS, newWSS

    real(dp), parameter :: pi = 3.1415926_dp   
    real(dp), parameter :: ln2 = log(2.0_dp)   
    real(dp), parameter :: ln2pi = 0.4697186393498_dp
    real(dp), parameter :: h = 1.0d-8
    real(dp), dimension(:), allocatable :: energii, intensitati
    real(dp), dimension(:,:), allocatable :: J, JT, JTJ, dJTJ
    real(dp), dimension(:), allocatable :: eExp, yExp
    
    real(dp), parameter :: TOL = 1.0D-4
    
! Ordine coloane in Jacobian
!    f s gama e0 

contains
    
    subroutine citesteExperiment(fisier)
        implicit none
        character(len=*), intent(IN) :: fisier
        character(len=80) :: linie
        integer :: n,eroare
        real(dp), dimension(:), allocatable :: etmp, ytmp
        n = 0
        open(19,file=trim(fisier),action="read")
        do
            read(19,'(a)',iostat=eroare) linie
            if(eroare == 0) then 
                n = n+1
                if(allocated(eExp)) then
                    allocate(etmp(1:n))
                    allocate(ytmp(1:n))
                    read(linie,*) etmp(n),ytmp(n)
                    etmp(1:n-1) = eExp(1:n-1)
                    ytmp(1:n-1) = yExp(1:n-1)
                    deallocate(eExp)
                    deallocate(yExp)
                    allocate(eExp(1:n))
                    allocate(yExp(1:n))
                    eExp(1:n) = etmp(1:n)
                    yExp(1:n) = ytmp(1:n)
                    deallocate(etmp)
                    deallocate(ytmp)
                else
                    allocate(eExp(1:1))
                    allocate(yExp(1:1))
                    read(linie,*) eExp(n),yExp(n)
                 end if
            else    
                exit
            end if  
        end do
        close(19) 
        nexp = n
        print*,"-------------------Nr exp:" ,nexp
    end subroutine    

    subroutine citesteCalcul(fisier)
        implicit none
        character(len=*):: fisier
        character(len=80) :: linie
        integer :: n,eroare
        real(dp) :: x,y
        real(dp), dimension(:), allocatable :: etmp, ytmp
        n = 0
        open(19,file=trim(fisier),action="read")
        do
            read(19,fmt='(a)',iostat=eroare) linie
            if(eroare == 0) then 
            
                if(allocated(energii)) then
                    read(linie,*) x,y
                    if(abs(y)>1.0d-10) then
                      n = n+1
                      allocate(etmp(1:n))
                      allocate(ytmp(1:n))
                      etmp(n)= x
                      ytmp(n)= y
                      etmp(1:n-1) = energii(1:n-1)
                      ytmp(1:n-1) = intensitati(1:n-1)
                      deallocate(energii)
                      deallocate(intensitati)
                      allocate(energii(1:n))
                      allocate(intensitati(1:n))
                      energii(1:n) = etmp(1:n)
                      intensitati(1:n) = ytmp(1:n)
                      deallocate(etmp)
                      deallocate(ytmp)
                    end if
                else
                  read(linie,*) x,y
                  if(abs(y)>1.0d-10) then
                    n= n+1 
                      allocate(energii(1:1))
                      allocate(intensitati(1:1))
                      energii(n)= x
                      intensitati(n)= y                                                              
                    end if
                 end if
            else    
                exit
            end if   
        end do
        close(19) 
        nvib = n
        e0 = energii(1)
    print*, "-------------------Nr vib-calc:",nvib, e0    
end subroutine

    subroutine initializare(f0,gamma0)
        real(dp), intent(in), optional :: f0, gamma0
        allocate(J(nexp,3))
        allocate(JT(3,nexp))
        allocate(JTJ(3,3))
        allocate(dJTJ(3,3))
        
        f = 0.960d0
        gama = 0.01d0
        lambda =  1.0d+3
        oldWSS = -1.0d0
    
        if (present(f0)) f = f0
        if (present(gamma0)) gama = gamma0

    end subroutine


    real(dp) function G(x,E,I,f,gama,e0)
        implicit none
        real(dp), intent(IN) :: E, I
        real(dp), intent(in) :: x, f, gama, e0     
        real(dp) :: t1, t2, t3 
        t1 = ln2pi*I/gama
        t2 = x-(f*E+(1.0d0-f)*(e0))
        t3 = ln2*t2**2/gama**2        
        G = t1*exp(-t3)
     end function

    real(dp) function ycalc(x,f1,gama1,e01)
        real(dp), intent(in) :: x, f1,  gama1, e01
        integer :: i
        real(dp) :: y  
        y = 0.0d0
            do i=1,nvib
                 !if (intensitati(i) > TOL) y = y + G(x,energii(i),intensitati(i), f1,s1,gama1,e01 )
                 y = y + G(x,energii(i),intensitati(i), f1,gama1,e01 )
             end do
        ycalc = y
    end function ycalc

        !e^(x^2) -> e^(x^2) +2x
        
        ! (e^(x+4)^2 - e^(x-4)^2) / h
   ! real(dp) function dyds(x)
    !    real(dp), intent(in) :: x      
     ! dyds = (ycalc(x,f,s+h,gama,e0)-ycalc(x,f,s-h,gama,e0))/(2.0d0*h)
    !end function dyds
    
    real(dp) function dydf(x)
        real(dp), intent(in) :: x      
        dydf = (ycalc(x,f+h,gama,e0)-ycalc(x,f-h,gama,e0))/(2.0d0*h)      
    end function dydf
    
    real(dp) function dydg(x)
        real(dp), intent(in) :: x
        dydg = (ycalc(x,f,gama+h,e0)-ycalc(x,f,gama-h,e0))/(2.0d0*h)      
    end function dydg

     real(dp) function dyde(x)
        real(dp), intent(in) :: x
        dyde = (ycalc(x,f,gama,e0+h)-ycalc(x,f,gama,e0-h))/(2.0d0*h)      
    end function dyde
    
    
    subroutine jacobian()
      implicit none
      integer :: i
      
      !$OMP PARALLEL
      do i=1,nexp
        J(i,1) = dydf(eExp(i))
        J(i,2) = dydg(eExp(i))
        J(i,3) = dyde(eExp(i))
      end do
      !$OMP END PARALLEL
      
      JT = transpose(J)
      JTJ = matmul(JT,J)
    end subroutine
  
    !Calculeaza suma de patrate    
    real(dp) function sumsq(y_exp,y_calc)  
        real(dp), dimension(nexp), intent(IN) :: y_exp,y_calc       
        real(dp), dimension(nexp) :: dif

        dif = y_exp - y_calc
        sumsq=dot_product(dif,dif) 

    end function

   subroutine scrieDate(a_xexp, a_yexp, a_ycalc)
        implicit none
        real(dp), dimension(nExp), intent(IN) :: a_xexp, a_yexp, a_ycalc
        real(dp), dimension(nExp) :: a_dif
        integer :: i
     
        a_dif = a_ycalc - a_yexp     
     
        open(19, file="date.dat", action="write")
            do i=1, nExp
                write(19,"(4F16.8)"), a_xexp(i), a_yexp(i), a_ycalc(i), a_dif(i)
            end do
        close(19)
   end subroutine
   
    
    subroutine Fit(nIter, eps)
      implicit none
      ! Numar maxim de iteratii
      integer, intent(IN), optional :: nIter
      ! Toleranta fit
      real(dp), intent(IN), optional :: eps
      ! Variabile locale
      real(dp), dimension(3,3) :: a
      real(dp), dimension(3) :: dpar, JdY
      real(dp), dimension(:), allocatable :: dY
      integer :: piv(3), info, i, iter
      real(dp), dimension(:), allocatable :: yc, ycold
      real(dp) :: norm
      !call restore 
      allocate(dY(1:nexp))
      allocate(yc(1:nexp))
      allocate(ycold(1:nexp))
      iter=0
      dpar = 0.0d0
     print '(/,35x,a)', "-------------------------------------"
        print '(35x,a)',   "|        Starting new fit           |"
     print '(35x,a,/)', "-------------------------------------"
     
     print '(a,/)',"  Iter  |      f        |       s       |     gamma    &
   & |     E0-0      |    OLD_wss    |    NEW_wss   |    Lambda   " 

      do
         iter = iter+1
         !$OMP PARALLEL
         yc=0.0
         
         do i=1,nexp
            yc(i) = ycalc(eexp(i),f,gama,e0)
         end do    
         yc = yc / maxval(yc(1:nexp))

         !TODO
         call scrieDate(eExp,yExp,yc)
         ! END TODO

         !$OMP END PARALLEL
         dY =  yexp-yc
         ycold = yc      
     

         call jacobian 
         !deallocate(dy)        
         a = JTJ
         call Diag(JTJ,3)
         a = a + lambda*JTJ      

        
        JdY = matmul(JT, dY)
        ! Rezolv A * x = B, unde
        ! A = a
        ! B = jdy
        call dgesv(3,1,a,3,piv,jdy,3,info)
        if (info /= 0) then 
             print "(A,I0,A)", "Eroare ", info, " la rezolvarea sistemului de ecuatii!"
             return
        end if
          
          dpar = jdy
        !$OMP PARALLEL
        yc=0.0
        
         do i=1,nexp
             yc(i) = ycalc(eexp(i),f + dpar(1),gama + dpar(2),e0+ dpar(3))
         end do
         !yc = yc / maxval(yc(1:nexp))
             
         !$OMP END PARALLEL
        ! actualizare lambda ?
        newWSS = sumsq(yExp,yc)
        !s = maxval(yexp) / maxval(yc)

        if (iter > 1) then
       print '(2x,I2,$)', iter
            if (newWSS > oldWSS) then
               lambda = 10 * lambda
         !      print '(A, F20.4)', "Lambda will increase 10x -> ", lambda
              ! call scrieDate(eExp,yExp, ycold)
            else
             !actualizare parametru
               f = f + dpar(1)
               !s = s + dpar(2)
               gama = gama +dpar(2)
               e0 = e0 + dpar(3)
               lambda = lambda/10
      !         print '(A, F20.4)', "Lambda will decrease 10x -> ", lambda
               call scrieDate(eExp,yExp,yc)
            end if
        else
                cycle
          
        end if

        print "(3x,3es16.8,$)", f,gama,e0
        print '(2es16.4, $)', oldWSS, newWss
        print '(F15.4)', lambda

        oldWSS = newWSS    

        if (iter > niTer) exit

        if (present(eps)) then
            norm = sqrt(dot_product(dpar,dpar))
            if (norm < eps) then
                print '(/,35x,a,/)', "Toleranta atinsa. :)"
                exit
            end if
       end if

        
   ! call scrieDate(eExp,yExp,yc)

       end do
   end subroutine

 

end module functii 

