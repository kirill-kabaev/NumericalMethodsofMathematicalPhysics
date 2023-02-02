program Jacobi
implicit none
include "omp_lib.h" 
integer :: i, j, k, p, l, nthreads, N
real, parameter :: eps=0.000001
integer, parameter :: num_stream=4
logical :: flag
real, dimension(10000) :: B, X
real, dimension(10000,10000) :: A
real :: Norm, ran
real :: sum1, sum2, t1, t2, t1par, t2par

open(unit=1,status='UNKNOWN',file='N.txt')
open(unit=2,status='UNKNOWN',file='resposl.txt')
open(unit=3,status='UNKNOWN',file='timeposl4.txt')
open(unit=4,status='UNKNOWN',file='timeparall4.txt')

!!
N = 5000
do k = 0, 50
    N=N+100
  
    call random_seed(N)
    do i=1, N
        B(i)= 1.
        do j=1, N 
            if (i==j) then
                CALL RANDOM_NUMBER(A(i, j))
                A(i, j)=A(i, j)*10+1
            else
                CALL RANDOM_NUMBER(A(i, j))
                A(i, j)=A(i, j)/(N-1)
            endif
        enddo
    enddo
    
    !проверяю норму матрицы, если <1 то сходимость есть
    
    print *, "N=", N
    call NormSyst(N)
    
    !t1=omp_get_wtime() 
    call LinearSystSolving(N, X)
    
    !print *, "X", X
    !t2=omp_get_wtime() 
    
    sum1=0.
    do i = 1, N
        sum1 = sum1 + X(i)
    end do

    !t1par=omp_get_wtime() 
    call LinearSystSolvingParal(N, X)
    !t2par=omp_get_wtime() 
    
    sum2=0.
    do i = 1, N
        sum2 = sum2 + X(i)
    end do

        write(4, 100) t2par-t1par
        write(3, 100) t2-t1
        write(2, 100) Sum2-Sum1
        write(1, *) N
!
end do

100 Format(f16.6) 

pause
contains

subroutine LinearSystSolving(N, X)
        integer, intent(in) :: N
        real, dimension(N), intent(inout) :: X 
        integer :: i, j, k
        real :: sum
        real :: max
        real, dimension(N) :: X0
        
        !инциализация нулевого решения X=X0
        do k = 1, N
            X(k)=B(k)/A(k,k)
            X0(k)=X(k)
        enddo
        p=0
        !бесконечный цикл
        !EXIT - точка выхода  
        !t1=omp_get_wtime() 
        do
            p=p+1
            !алгоритм в отчете
            do i = 1, N
                sum = 0.
                do j = 1, N
                    if (i /= j) then 
                        sum = sum + A(i,j)*X0(j)
                    endif
                enddo
                X(i)=(B(i) - sum)/A(i,i)
            enddo
          
            !print *, "X posl ", X
            !допустим максимум в разности между нулевым и первым приближением
            max = abs(X(1)-X0(1))
            !print *, "X posl", X
            !print *, "X0 posl", X0
            !ищем максимум в разности приближений X(i+1)-X(i)
            do i = 1, N
                if (abs(X(i)-X0(i))>max) then
                    max = abs(X(i)-X0(i))
                endif
                X0(i)=X(i)
            enddo
            !если этот макуимум меньше заданной точности, выходим из цикла
            !print *, "t2-t1 ", t2-t1
            print *, "itter posl ", p
            !print *, "max ", max
            print *, "max posl ", max
            if (max<Eps) EXIT
        enddo
        !t2=omp_get_wtime() 
end subroutine LinearSystSolving

subroutine LinearSystSolvingParal(N, X)
        integer, intent(in) :: N
        real, dimension(N), intent(inout) :: X 
        integer :: i, j, k
        real :: sum
        real :: max
        real, dimension(N) :: X0
        !Initial values X=X0
        do k = 1, N
            X(k)=B(k)/A(k,k)
            X0(k)=X(k)
        enddo
        p=0
        !t1par=omp_get_wtime() 
        do
            p=p+1
            !shared - переменные доступыне всем потокам
            !private - доступные одному потоку
            !$omp parallel shared(A, X0, B, X) private(i, j, sum) num_threads(num_stream)
            !$omp do
                do i = 1, N
                    sum = 0.
                    do j = 1, N
                        !print *, "i=", i, "j=", j
                        if (i /= j) then 
                            sum = sum + A(i,j)*X0(j)
                        endif
                    enddo
                    X(i)=(B(i) - sum)/A(i,i)
                enddo
            !$omp end parallel
            max = abs(X(1)-X0(1))
            !print *, "X par", X
            !print *, "X0 par", X0
            do i = 1, N
                if (abs(X(i)-X0(i))>max) then
                    max = abs(X(i)-X0(i))
                endif
                X0(i)=X(i)
            enddo
            print *, "itter par ", p
            !print *, "t2par-t1par ", t2par-t1par
            print *, "max par ", max
            if (max<Eps) EXIT
        enddo
        !t2par=omp_get_wtime()
end subroutine LinearSystSolvingParal

subroutine NormSyst(N)
            integer, intent(in) :: N
            real, dimension(N) :: BsumRow
            real :: Norm
        
            do i = 1, N
                BsumRow(i)=0.
                do j = 1, N
                    if (i /= j) then
                        BsumRow(i) = BsumRow(i) + (A(i,j)/A(i,i))
                    endif
                end do
            end do
            
            !Print *, "NormaS =", BsumRow
            Norm=MAXVAL(BsumRow)
            Print *, "Norma = ", Norm, " N= ", N
            !Print *, "BsumRow= ", BsumRow
end subroutine NormSyst

end 

