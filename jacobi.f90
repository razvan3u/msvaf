module jacobi

   use kinds

   implicit none
   contains
        subroutine Diag(A, N)
            implicit none

            integer, intent(IN) :: N
            real(dp), intent(INOUT), dimension(N,N)  :: A

            integer :: i,j,k,jup
            real(dp), dimension(N,N)  :: U
            real(dp)  :: diffr, aii, ajj, aod, asq, amax, tank, tden, s, c, singn, xj
            real(dp), parameter :: eps = 1.0d-10
        
          ! U - matrice unitate
            U = 0.0      
            do i=1,N
               U(i,i) = 1.0
            end do
                
          ! Necesar?
            do i=1,n
               A(i,i) = -A(i,i)
            end do
           
          do
            amax = 0.0
            
              do i=2,N
                 jup = i-1
                 do j=1,jup
                    aii = A(i,i)
                    ajj = A(j,j)
                    aod = A(i,j)
                    asq = aod*aod
                    if (asq > amax) amax = asq
                    if (asq > eps) then 
                        diffr = aii-ajj
                        if (abs(diffr) > 0.0) then
                            singn =  2.0
                        else 
                            singn = -2.0
                        end if
                        tden = diffr + sqrt(diffr*diffr + 4*asq)
                        tank = singn * aod / tden
                        c = 1 / (sqrt(1+tank*tank))
                        s = c * tank
                        do k=1,n
                           xj     = c * u(k,j) - s * u(k,i)
                           u(k,i) = s * u(k,j) + c * u(k,i)
                           u(k,j) = xj
                           if (k < j) then
                               xj     = c * A(j,k) - s * A(i,k)
                               A(i,k) = s * A(j,k) + c * A(i,k)
                               A(j,k) = xj
                           end if
                           if (k > j) then 
                               if (k < i) then
                                   xj     = c * A(k,j) - s * A(i,k)
                           A(i,k) = s * A(k,j) + c * A(i,k)
                           A(k,j) = xj
                               end if
                               if (k > i) then
                                   xj     = c * A(k,j) - s * A(k,i)
                           A(k,i) = s * A(k,j) + c * A(k,i)
                           A(k,j) = xj
                               end if
                           end if
                        end do 
                        A(i,i) = c*c*aii + s*s*ajj + 2*s*c*aod
                        A(j,j) = c*c*ajj + s*s*aii - 2*s*c*aod
                        A(i,j) = 0.0
                        A(j,i) = 0.0
                    end if
                 end do
              end do
              if (abs(amax) < eps) then 
                 A = U  
               exit 
            end if
            end do
               
        end subroutine Diag
end module

