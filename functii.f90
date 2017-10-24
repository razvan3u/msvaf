!------------------------------------------------------------------------------!
! Modul de functii si subrutine pentru calculul de regresie                    !
! neliniara !                                                                  !
!------------------------------------------------------------------------------!
! Copyright (c) 2015-2016, Razvan Cirdei <razvan.cirdei@chem.uaic.ro>          !
!                          Dan Maftei    <dan.maftei@chem.uaic.ro>             !
!------------------------------------------------------------------------------!
module functii

     use kinds
     use jacobi

     implicit none
     integer  :: nexp, nvib 
     real(dp) :: f = 1.0_dp, gama=0.03_dp,  e0  
     real(dp) :: lambda, oldWSS, newWSS
     real(dp), parameter :: pi    = 3.1415926_dp   
     real(dp), parameter :: ln2   = log(2.0_dp)   
     real(dp), parameter :: ln2pi = 0.4697186393498_dp
     real(dp), parameter :: h     = 1.0d-10
     real(dp), parameter :: TOL = 1.0D-4
     real(dp), dimension(:),   allocatable :: energii, intensitati
     real(dp), dimension(:,:), allocatable :: J, JT, JTJ, dJTJ
     real(dp), dimension(:),   allocatable :: eExp, yExp

     ! Ordine coloane in Jacobian : f s gama e0 
     contains
          subroutine citesteExperiment(fisier)

               implicit none
               integer :: n, eroare
               character(len=* ), intent(IN) :: fisier
               character(len=80) :: linie
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
                                   read(linie,*) eExp(n), yExp(n)
                              end if
                         else    
                              exit
                    end if  
               end do
               close(19) 
               nexp = n
               ! Normalizare
               yExp = yExp/maxval(yExp(1:nexp))
               print*,"------Numar date experimentale:", nexp

          end subroutine    

          subroutine citesteCalcul(fisier)

               implicit none
               integer :: n, eroare, i
               character(len=* ) :: fisier
               character(len=80) :: linie
               real(dp) :: x, y, xmin, xmax
               real(dp), dimension(:), allocatable :: etmp, ytmp

               n = 0
               xmin = minval(eExp(1:nexp)) - 0.5
               xmin = maxval(eExp(1:nexp)) + 0.5
               open(19,file=trim(fisier),action="read")
                    do
                         read(19,fmt='(a)',iostat=eroare) linie
                         if(eroare == 0) then 
                              if(allocated(energii)) then
                                   read(linie,*) x,y
                                   !if((abs(y)>1.0d-10 ) .and. (x > xmin) .and. (x < xmax)) then
                                   if((abs(y)>1.0d-10 )) then
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
               ! Pastrez toate picurile vibronice in unitati de energie
               ! relative la energia de zero, pentru simplitate
               energii = energii - e0
               print*, "------Numar vibratii calculate:", nvib, e0   
 
          end subroutine

          subroutine setParams(set_e0, set_f0, set_gamma0)
               real(dp), intent(in), optional :: set_e0, set_f0, set_gamma0
               if (present(set_e0))     e0 = set_e0
               if (present(set_f0))     f = set_f0
               if (present(set_gamma0)) gama = set_gamma0
          end subroutine

          subroutine initializare()
          
               allocate(J   (nexp,3   ))
               allocate(JT  (3   ,nexp))
               allocate(JTJ (3   ,   3))
               allocate(dJTJ(3   ,   3))

          end subroutine

          real(dp) function G(x, E, I, f, gama, e0)
               ! Functie de profil Gauss:
               ! Argumente
               !    x - valoarea abscisei in care se evalueaza
               !    E - pozitia maximului, relativa la e0
               !    I - intensitatea 
               ! gama - latimea
               !   e0 - energia de zero

               implicit none
               real(dp) :: t1, t2, t3
               real(dp), intent(IN) :: E, I
               real(dp), intent(in) :: x, f, gama, e0      
           
!               if(abs(x-E) < 10 * gama) then
                    t1 = ln2pi * I / gama
                    t2 = x - (e0 + f*E)
                    t3 = ln2 * t2**2 / gama**2        
                    G = t1 * exp(-t3)
               
!               else
!                       G = 0.0
!               end if
          end function

          real(dp) function ycalc(x, f1, gama1, e01)
               
               implicit none
               integer  :: i
               real(dp) :: y  
               real(dp), intent(in) :: x, f1,  gama1, e01
          
               y = 0.0d0
               do i=1, nvib
                    !if (intensitati(i) > TOL) y = y + G(x,energii(i),intensitati(i), f1,s1,gama1,e01 )
                    y = y + G(x, energii(i), intensitati(i), f1, gama1, e01 )
               end do
               ycalc = y 

          end function ycalc


          subroutine jacobian()
               
               implicit none
               integer  :: i, jj
               real(dp) :: termexp, term2, sf, sg, se
          
               !$OMP PARALLEL

               do i=1,nexp
                    sg = 0.0
                    se = 0.0
                    sf = 0.0
                    do jj = 1, nvib
                         sf = sf + (G(eExp(i), energii(jj), intensitati(jj), f + h, gama,     e0     )- &
                                    G(eExp(i), energii(jj), intensitati(jj), f - h, gama,     e0    )) / (2*h)
                         sg = sg + (G(eExp(i), energii(jj), intensitati(jj), f,     gama + h, e0     )- &
                                    G(eExp(i), energii(jj), intensitati(jj), f,     gama - h, e0    )) / (2*h)         
                         se = se + (G(eExp(i), energii(jj), intensitati(jj), f,     gama,     e0 + h )- &
                                    G(eExp(i), energii(jj), intensitati(jj), f,     gama,     e0 - h)) / (2*h)     
                    end do
                    J(i,1) = sf
                    J(i,2) = sg
                    J(i,3) = se
              end do

               !$OMP END PARALLEL

               JT  = transpose(J)
               JTJ = matmul(JT,J)

          end subroutine

          !Calculeaza suma de patrate    
          real(dp) function sumsq(y_exp, y_calc)  

               implicit none
               real(dp), dimension(nexp), intent(IN) :: y_exp, y_calc       
               real(dp), dimension(nexp) :: dif

               dif   = y_exp - y_calc
               sumsq = dot_product(dif, dif) 

          end function

          subroutine scrieDate(a_xexp, a_yexp, a_ycalc)

               implicit none
               integer :: i
               real(dp), dimension(nExp), intent(IN) :: a_xexp, a_yexp, a_ycalc
               real(dp), dimension(nExp) :: a_dif
               
               a_dif = a_ycalc - a_yexp     
               open(19, file="date.dat", action="write")
                    do i=1, nExp
                         write(19,"(4F16.8)"), a_xexp(i), a_yexp(i), a_ycalc(i), a_dif(i)
                    end do
               close(19)
               open(20, file="stick.dat", action="write")
                    do i = 1, nVib
                         write(20,"(2f16.8)"), (f*energii(i) + (1.0d0-f)*(e0)), intensitati(i)
                    end do
               close(20)

          end subroutine

          subroutine Fit(nIter, eps)

               implicit none
               integer,  intent(IN), optional :: nIter !numar maxim de iteratii
               real(dp), intent(IN), optional :: eps !toleranta atinsa
               integer :: piv(3), info, i, iter
               real(dp), dimension(3,3) :: a
               real(dp), dimension(3)   :: dpar, JdY
               real(dp), dimension(:), allocatable :: dY
               real(dp), dimension(:), allocatable :: yc, ycold
               real(dp) :: norm
            
               allocate(dY(1:nexp))
               allocate(yc(1:nexp))
               allocate(ycold(1:nexp))
               iter=0
               dpar = 0.0d0
               print '(/,35x,a)', "-------------------------------------"
               print '(35x,a)',   "|        Starting new fit           |"
               print '(35x,a,/)', "-------------------------------------"
               print '(a,/)',"  Iter  |      f        |     gamma    &
               & |     E0-0      |    OLD_wss    |    NEW_wss   |    Lambda   " 

               lambda =  1.0d+2
               oldWSS = -1.0d0
               do
                    iter = iter + 1

                    !$OMP PARALLEL

                    yc=0.0
                    do i=1, nexp
                         yc(i) = ycalc(eexp(i), f, gama, e0)
                    end do    
                    ! Normalizare
                    yc = yc / maxval(yc(1:nexp))

                    call scrieDate(eExp, yExp, yc)
                    !$OMP END PARALLEL

                    dY    = yexp - yc
                    ycold = yc      
                    call jacobian   
                    a = JTJ
                    call Diag(JTJ,3)
                    a = a + lambda*JTJ      
                    JdY = matmul(JT, dY)
                    call dgesv(3, 1, a, 3, piv, jdy, 3, info)  ! Rezolv A * x = B, unde A = a B = jdy
                    if (info /= 0) then 
                         print "(A,I0,A)", "Eroare ", info, " la rezolvarea sistemului de ecuatii!"
                         return
                    end if
                    dpar =jdy

                    !$OMP PARALLEL

                    yc=0.0
                    do i=1, nexp
                         yc(i) = ycalc(eexp(i), f + dpar(1), gama + dpar(2), e0+ dpar(3))
                    end do
                    yc = yc / maxval(yc(1:nexp))

                    !$OMP END PARALLEL

                    newWSS = sumsq(yExp,yc)  ! actualizare lambda ?
                    !s = maxval(yexp) / maxval(yc)
                    if (iter > 1) then
                         if (newWSS > oldWSS) then
                              lambda = 10 * lambda
                              !print '(A, F20.4)', "Lambda will increase 10x -> ", lambda
                              !call scrieDate(eExp, yExp, ycold)
                         else  !actualizare parametru
                              f    = f    + dpar(1)
                             !s    = s    + dpar(2)
                              gama = gama + dpar(2)
                              e0   = e0   + dpar(3)
                              lambda = lambda/10
                              !print '(A, F20.4)', "Lambda will decrease 10x -> ", lambda
                              call scrieDate(eExp,yExp,yc)
                         end if
                    else
                         cycle
                    end if
                    print '(2x,I4,$)', iter
                    print '(6x,$)'
                    print "(3x,3F16.4,$)", f, gama, e0
                    print '(2es16.4, $)', oldWSS, newWss
                    print '(ES18.4)', lambda
                    oldWSS = newWSS    
                    if (iter > niTer) exit
                    if (present(eps)) then
                         norm = sqrt(dot_product(dpar,dpar))
                         if (norm < eps) then
                              print '(/,35x,a,/)', "Toleranta atinsa. :)"
                              !exit
                         end if
                    end if
               !call scrieDate(eExp,yExp,yc)
               end do

          end subroutine

end module functii 

