program integral
    implicit none
    include "omp_lib.h"
    integer :: N
    integer, parameter :: num_stream = 4
    real, parameter :: x0 = 0., x1 = 10.
    real :: Sum1, Sum2, t1, t2, t1p, t2p, dt
    integer :: i,k
    

    open(unit=1,status='UNKNOWN',file='numstreams2.txt')
    open(unit=2,status='UNKNOWN',file='numstreams1.txt')
    open(unit=3,status='UNKNOWN',file='N.txt')
    open(unit=4,status='UNKNOWN',file='res1.txt')
    
    N=1000000
    do k=0, 40
        N=N+k*100000
 
        t1p=omp_get_wtime() 
        call integralSimpsonparallel(F, x0, x1, N, sum1)
        t2p=omp_get_wtime() 
        
        t1=omp_get_wtime() 
        call integralSimpson(F, x0, x1, N, sum2)
        t2=omp_get_wtime()
        
        write(1, 1) t2p-t1p
        write(2, 1) t2-t1
        write(3, *) N
        write(4, *) Sum2-Sum1
        enddo
    1 Format(f18.8) 
    
100 Format (F18.8)
    
    pause
    contains
    
    !задаем произвольную функцию
    real function F(x)
        real :: x
        F = sin(x)*x
    end function F
    
    !процедура для интеграла
    subroutine integralSimpsonparallel(F, x0, x1, N, sum)
       real, external :: F
       real, intent(in) :: x0
       real, intent(in) :: x1
       integer, intent(in) :: N
       real, intent(inout) :: sum
       
       real, dimension(N+1) :: Ff
       
       dt = (x1-x0)/N 
       do i = 0, N
           Ff(i+1) = F(i*dt) !создаем массив значений функции с шагом dt
       enddo
       
    sum = 0.
    !$omp parallel reduction(+: sum) private(i) num_threads(num_stream)
    !$omp do
       do i = 0, int((N/2-1))
           sum = sum + dt/3*(Ff(2*i+1)+4*Ff(2*i+2)+Ff(2*i+3)) !суммируем по формуле
       enddo
    !$omp end parallel
    end subroutine integralSimpsonparallel
    
    subroutine integralSimpson(F, x0, x1, N, sum)
       real, external :: F
       real, intent(in) :: x0
       real, intent(in) :: x1
       integer, intent(in) :: N
       real, intent(inout) :: sum
       
       real, dimension(N+1) :: Ff
       
       dt = (x1-x0)/N
       do i = 0, N
           Ff(i+1) = F(i*dt)
       enddo
       
       sum = 0.
       do i = 0, int((N/2-1))
           sum = sum + dt/3*(Ff(2*i+1)+4*Ff(2*i+2)+Ff(2*i+3))
       enddo
    end subroutine integralSimpson
    
    end program integral

