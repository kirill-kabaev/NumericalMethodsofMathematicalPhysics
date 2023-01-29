program RungeKutta4
    implicit none
    integer, parameter :: N = 200
    real, parameter :: g1 = 0.5
    real :: Uk = 100.
    real :: Un = 110.
    real :: Xk = 12.
    real :: dd0 = 0.002
    real :: H0 = 3.8
    real :: H10 = 3.8
    real :: L = 1
    real :: Ue, C, x, dx
    real, dimension(2) :: F ! 1: d, 2: H1
    real, dimension(N) :: H
    integer :: i
    real, dimension(2) :: k1, k2, k3, k4
    
    F(1) = dd0*H0 ! d(t=0) = dd(t=0)*H
    F(2) = H10 ! H1(t=0)
    x = L
    dx = (Xk - L)/N
    H(1) = H0
    
    open(unit=1,status='unknown',file='dx.txt')
    open(unit=2,status='unknown',file='H1x.txt')
    open(unit=3,status='unknown',file='Hx.txt')

    
    
    do i = 1, N

        
        if (i>1) then
            if ((F(2) - 2)**2 - 3 <=0) then
                    F(2) = H10
                    H(i) = 1+1.12*(F(2)-2-sqrt((F(2)-2)*(F(2)-2)-3))**0.915
                else
                    H(i) = 1+1.12*(F(2)-2-sqrt((F(2)-2)*(F(2)-2)-3))**0.915
            endif
        endif
        
        if (x<=6*L) then 
            C=0.0299*(F(2)-3)**(-0.6169)  ! x<6L
            else if ((x>6*L) .and. (x<10*L)) then
                C=g1*(0.435*(H(i)-1)**0.907) + (1-g1)*(0.0299*(F(2)-3)**(-0.6169) )  ! 6L<x<10L
            else if (x>=10*L) then
                C=0.435*(H(i)-1)**0.907  ! x>10L
        endif
        
            
        k1 = dx*RK(x, F, H(i))
        k2 = dx*RK(x+0.5*dx, F+0.5*dx*k1, H(i))
        k3 = dx*RK(x+0.5*dx, F+0.5*dx*k2, H(i))
        k4 = dx*RK(x+dx, F+dx*k3, H(i))
        F = F + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)
        x = x + dx  
        
        
        WRITE(1,220) x, F(1)/H(i) ! dd(x)
        WRITE(2,220) x, F(2) ! H1(x)
        WRITE(3,220) x, H(i)
        
    end do
      
    print *, "Cx1 = ", 2*F(1)/(H(N)*L)
    print *, "Cx2 = ", 2*F(1)/(H(1)*L)*(0.99)**((2.918+5)*0.5)
    
    close(1)
    close(2)
    close(3)
    220 Format (F18.8, F18.8)
    !230 Format ('Progress: ', I3)
       
    pause
    contains
      
    function RK(x, F, H)
        real, dimension(2) :: RK
        real, intent(in) :: x
        ! real, dimension(2), intent(in) :: d ! dd=H/d
        real, dimension(2), intent(in) :: F
        real, intent(in) :: H
        real :: Ue
        Ue = Un+(Uk-Un)*x/Xk
        RK(1) = -(H+2)*(F(1)/H)*(Uk-Un)/(Xk*Ue)
        RK(2) = (H/F(1))*(2*C-F(2)*(-(H+1)*(F(1)/H)*(Uk-Un)))/(Xk*Ue)
    end function RK
       
end program RungeKutta4

