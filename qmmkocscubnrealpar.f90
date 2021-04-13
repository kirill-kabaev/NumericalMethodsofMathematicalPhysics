program MMK
    !include "mkl.fi"
    !include "omp_lib.h"
    
    include "omp_lib.h"
    
    integer, parameter :: num_stream = 4
    integer, parameter :: N = 2
    integer, parameter :: Nrealiz = 1
    complex, dimension(N, N) :: Hamiltonian, Ed!, Hh, Hc, Hi
    complex, dimension(2*N, 2*N) :: c, Hcu, phprod, enprod, bb
    complex, dimension(3,2*N) ::  Htt
    complex, dimension(2*N) :: k1, k2, k3, k4, fosc, foc, focn!, b
    real, dimension(N) :: a1, a2, a
    real, dimension(2) :: EigenValues, fcub
    complex, dimension(2,2) :: ff, H
    complex, dimension(4, 4) :: matrix
    complex, dimension(2,2) :: Sx,  Sz , II, Splus, Sminus, S1, Sed
    integer :: i,j,k,p, tn, l, u, z
    real :: s, ran, sumvec, ts, Pe, t, Am, temp1, temp2, Posc, P1, P2,  f
    real, parameter :: pi = 3.141592653589, dt = 0.1, ddt = dt/20., nt = 1.
    real, parameter :: wq = 6. !qubit frequency
    real, parameter :: w = 6. !frequency of the field exciting the qubit
    real, parameter :: Gf = 0.0 !phase relaxation
    real, parameter :: Ge = 0.0 !energy relaxation
    real, parameter :: Go= 0.0 !damping of the oscillator
    real, parameter :: wfo=1. !frequency of the field acting on the oscillator
    real, parameter :: f0=0.0!amplitude of the oscillator field
    !real, parameter :: f=0.9 !amplitude of the oscillator field 
    real, parameter :: Ams=0.5 !qubit field amplitude
    real, parameter :: Mu  = 0. !0.0001!non-linearity
    real, parameter :: lambda=0. !Interaction
    real, parameter :: wj=1. !1. !oscillator frequency
    real, parameter :: dwo=0. 
    real, parameter :: Wrab = ((wq-w)*(wq-w)+Ams*Ams) ** 0.5
    real, parameter :: period = 2*pi/Wrab
    integer, parameter :: tperiod = 30, Nstep = int(2 * pi /dt), tmax = Nstep*tperiod
    real, dimension(tmax) :: qp1, qp2, nr1, nr2, nr
    !real, dimension(Nstep*Tperiod*Nstep*Tperiod) :: mnr 
    real :: t1, t2
    
     !Sx = reshape((/(0.,0.), (1.,0.), (1., 0.), (0., 0.)/), shape(Sx))
    Sx(1,1) = (0.,0.)
    Sx(1,2) = (1.,0.)
    Sx(2,1) = (1.,0.)
    Sx(2,2) = (0.,0.)
    !Sz = reshape((/(1.,0.), (0.,0.), (0., 0.), (-1., 0.)/), shape(Sz))
    Sz(1,1) = (1.,0.)
    Sz(1,2) = (0.,0.)
    Sz(2,1) = (0.,0.)
    Sz(2,2) = (-1.,0.)
    !II = reshape((/(1.,0.), (0.,0.), (0., 0.), (0., 1.)/), shape(II))
    II(1,1) = (1.,0.)
    II(1,2) = (0.,0.)
    II(2,1) = (0.,0.)
    II(2,2) = (1.,0.)
    !Splus = reshape((/(0.,0.), (1.,0.), (0., 0.), (0., 0.)/), shape(Splus))
    Splus(1,1) = (0.,0.)
    Splus(1,2) = (1.,0.)
    Splus(2,1) = (0.,0.)
    Splus(2,2) = (0.,0.)
    !Sminus = reshape((/(0.,0.), (0.,0.), (1., 0.), (0., 0.)/), shape(Sminus))
    Sminus(1,1) = (0.,0.)
    Sminus(1,2) = (0.,0.)
    Sminus(2,1) = (1.,0.)
    Sminus(2,2) = (0.,0.)
    
    S1(1,1) = (1.,0.)
    S1(1,2) = (0.,0.)
    S1(2,1) = (0.,0.)
    S1(2,2) = (0.,0.)
    
    Sed(1,1) = (1.,0.)
    Sed(1,2) = (0.,0.)
    Sed(2,1) = (0.,0.)
    Sed(2,2) = (1.,0.)
   
    
    do i = 1, N
        do j = 1, N
            Ed(i,j) = (0., 0.)
        enddo
        Ed(i,i) = (1., 0.)
    enddo
    do j = 1, tmax
        nr(j) = 0.
        nr1(j) = 0.
        nr2(j) = 0.
        qp1(j) = 0.
    enddo
   
    call formedStartHamiltInter()
    t1 = omp_get_wtime()
    call QuantumMethodMonteCarlo(N)
    t2 = omp_get_wtime()
    print *, "Time: ", t2 - t1, " sec"
    



    pause
    contains

    subroutine QuantumMethodMonteCarlo(N)
        open(unit=1,status='unknown',file='qp1.txt')
        open(unit=2,status='unknown',file='qp2.txt')
        open(unit=3,status='unknown',file='nr1.txt')
        open(unit=4,status='unknown',file='nr2.txt')
        open(unit=5,status='unknown',file='nr.txt')
        open(unit=6,status='unknown',file='oscillRabbi.txt')
        

            do p = 1, N
                fosc(p)=(0.,0.) 
            end do
            fosc(1)=(1.,0.)
            
            !qubit wave functions
            call EigenVecVal(H0(), EigenValues, ff, 2) 
            do p = 1, 2
                fcub(p)=ff(p, 1) 
            end do
            
            do i = 1, N
                do j = 1, 2
                    foc(k) = ff(j,1)*fosc(i)
                    k=k+1
                enddo
            enddo
            sumvec=0.
    call KronProd(Ed, Sz , phprod)
    call KronProd(Ed, Sminus, enprod)
    call KronProd(Ed, S1 , bb)
        !print *, "bb_0 = ", bb
    u=0
    !$omp parallel reduction(+:nr, qp1, nr1, nr2) private(temp1, Posc, P1, P2, i, j, k, z, t, ts, k1, k2, k3, k4, G, Pe, focn) shared(ran ,s, u) num_threads(num_stream)
    !$omp do
        do k = 1, Nrealiz
            focn=foc
            call random_seed()
            do j = 1, tmax
                t=j*dt 
                temp1 = 0.
                Posc = 0.
                P1 = 0.
                P2 = 0.
                do i = 1, N
                    temp1 = temp1 + abs(CONJG(focn(2*i-1))*focn(2*i-1)) !population of the BASIC level of the qubit
                    Posc = Posc + (i-1)*(abs(CONJG(focn(2*i-1))*focn(2*i-1)) + abs(CONJG(focn(2*i))*focn(2*i)))
                    P1 = P1 + (i-1)*abs(CONJG(focn(2*i-1))*focn(2*i-1))
                    P2 = P2 + (i-1)*abs(CONJG(focn(2*i))*focn(2*i))
                    !print *, "P1=", P1, "P2=", P2
                enddo
                
                nr1(j) = nr1(j) + P1
                nr2(j) = nr2(j) + P2
                nr(j) = nr(j) + Posc
                !write (*,*) "nr1 = ", nr1(j), "nr2 = ", nr2(j)
                !write (*,*) "nr = ", nr(j)
                qp1(j) = qp1(j) + temp1  
                
                Pe = abs(DOT_PRODUCT(CONJG(focn), matmul(bb,focn)))
                G = Gf + Ge*Pe+Go*Posc
                call RANDOM_NUMBER(ran)
                if (ran < dt*G) then
                    call RANDOM_NUMBER(s)
                    if (s<Gf/G) then 
                        focn = matmul(phprod, focn)
                        !print*, "phase - ", j
                    else
                        if ((s > Gf/G).AND.(s <(Gf/G + Ge*Pe/G))) then
                            focn = matmul(enprod, focn)
                            !print*, "energy - ", j
                        else 
                            !print*, "osc - ", j
                            do l = 2, N
                                focn(2*l-3) = focn(2*l-1)/(l ** 0.5)
                                focn(2*l-2) = focn(2*l)/(l ** 0.5)
                            enddo
                        endif
                    endif         
                else
                    ts = t
                    do while (ts < (j+1)*dt)
                        k1 = ddt*RK(focn,ts)
                        k2 = ddt*RK(focn+0.5*k1,ts+0.5*ddt)
                        k3 = ddt*RK(focn+0.5*k2,ts+0.5*ddt)
                        k4 = ddt*RK(focn+k3,ts+ddt)
                        focn = focn + 1./6.*(k1 + 2.*k2 + 2.*k3 + k4)
                        ts = ts + ddt
                    end do
                
                endif
                !write (*,*) "focn(1) = ", focn(1), "focn2 = ", focn(2)
                !write (*,*) "Nrealiz = ", k, "RK(focn,ts) = ", RK(focn,ts)
                sumvec=0.
                do z = 1, 2*N
                    sumvec = sumvec + CONJG(focn(z))*focn(z)
                end do 
                focn = focn/sqrt(abs(sumvec))
            end do
            u=u+1
            write (*,*) "Nrealiz", u
        end do
    !$omp end parallel

        do i = 1, tmax
                WRITE(5,220) nr(i)/Nrealiz
                WRITE(3,220) nr1(i)/Nrealiz
                WRITE(4,220) nr2(i)/Nrealiz
                WRITE(1,220) qp1(i)/Nrealiz
                t=i*dt
                WRITE(6,220) Ams*Ams*cos(t*Wrab/4)*cos(t*Wrab/4)/(Wrab*Wrab)
        end do
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        220 Format (F18.8)
        230 Format ('Progress: ', I3)
    end 
    


    
    function Ht(t)
        complex, dimension(3,2*N) ::  Ht
        real :: t, Am
        Am = Ams
        !if (t<=period) then 
        !    Am = Ams
        !else 
        !    Am = 0.
        !endif
        Ht=Htt
        Ht(1,:) = Htt(1,:)*fun(t)
        Ht(2,:) = Htt(2,:)*Am*cos(w*t)/2.
    end function  Ht

    
    function RK(fs,t)
        complex, dimension(2*N) :: Rk
        real :: t
        complex, dimension(2*N) :: fs
        complex, dimension(2*N) :: Htfs
        !write (*,*) "Rk = ", Rk
        call chbmv('U', 2*N, 2, (0., -1.), Ht(t), 3, fs, 1, (0., 0.), Htfs, 1)
        Rk = Htfs
    end function RK
    

    

    function fun(t)
        real :: fun
        real :: t, f
        !if (t<=10*period) then 
        !    f = 0.
        !else
        !    f = f0
        !endif
        f = f0
        fun = f*cos(wfo*t)
    end function fun

    
    function H0() 
        complex, dimension(2,2) :: H0
        H0 = 1/2. * wq * Sz
    end function H0
    
    

    
    subroutine formedStartHamiltInter()
        complex, dimension(2*N,2*N) ::  Hinter
        complex, dimension(2*N,2*N) ::  Hqubit
        complex, dimension(2*N,2*N) ::  Hoscill
        complex, dimension(2*N,2*N) ::  Htmatrix
        integer i, j, m
        complex, dimension(N,N) :: Hh, Hc, Hi
        do i = 1, N
            do j = 1, N
                Hc(i,j) = (0.,0.)
                Hh(i,j) = (0.,0.)
                Hi(i,j) = (0.,0.)
                Hoscill(i,j) = (0.,0.)
            end do
        end do
        do i = 1, 2*N
            do j = 1, 2*N
                Hoscill(i,j) = (0.,0.)
            end do
        end do
        do i = 1, 2*N
            do j = 1, 2*N
                Hinter(i,j) = (0.,0.)
                Hqubit(i,j) = (0.,0.)
                Hi(i,j) = (0.,0.)
                Htmatrix(i,j) = (0.,0.)
            end do
        end do
        
        do i = 0, N-1
            Hc(i+1,i+1) = i*wj - Mu*i*i - (0.,1.)*i*Go/2.
        end do
        do i = 1, N-1
            Hh(i,i+1) = i ** 0.5
            Hh(i+1,i) = i ** 0.5
        enddo
        do i = 0, N-1
            Hi(i+1,i+1) = lambda*i*wj-Mu*i*i*lambda/4.
        enddo
        call KronProd(Hi, Sz , Hinter)
        call KronProd(Ed, H0() + Sx - (0.,1.)*(Gf/2.)*II - (0.,1.)*(Ge/2.)*S1 , Hqubit)
        call KronProd(Hh + Hc, II , Hoscill)
        Htmatrix = Hqubit + Hoscill + Hinter
        do j = 1, 2*N
            m=3-j
            do i = max(1, j - 2), j
                !print *, "m+i=", m + i, "j=", j, "****", "i=", i, "j=", j
                Htt(m + i, j) = Htmatrix(i, j)
            enddo
        enddo
    end
    
    subroutine KronProd(A,B,C)
       IMPLICIT NONE
       complex, dimension (:,:), intent(in)  :: A, B
       complex, dimension (size(A,1)*size(B,1),size(A,2)*size(B,2)), intent(inout) :: C
       integer :: i = 0, j = 0, k = 0, l = 0
       integer :: m = 0, nn = 0, p = 0, q = 0
       C=(0.,0.)
        do i = 1,size(A,1)
            do j = 1,size(A,2)
                nn=(i-1)*size(B,1) + 1
                m=nn+size(B,1) - 1
                p=(j-1)*size(B,2) + 1
                q=p+size(B,2) - 1
                C(nn:m,p:q) = A(i,j)*B
            enddo
        enddo 
    end 
    
    subroutine EigenVecVal(Matrix, EigenValues, EigenVectors, M)
        integer, intent(in) :: M
        complex, dimension(M,M), intent(inout) :: EigenVectors
        complex, dimension(M,M), intent(in) :: Matrix
        real, dimension(M), intent(inout) :: EigenValues
        integer, parameter :: LWMAX = 1000
        real, dimension(3*M-2) :: RWORK
        complex WORK(LWMAX)
        integer INFO, LWORK
        !documentation link:
        !software.intel.com/
        !cheev
        INFO=0
        LWORK = -1
        call cheev( 'V', 'U', M, Matrix, M, EigenValues, WORK, LWORK, RWORK, INFO)
        LWORK = MIN(LWMAX, INT(WORK(1)))
        call cheev( 'V', 'U', M, Matrix, M, EigenValues, WORK, LWORK, RWORK, INFO)
        EigenVectors = Matrix
        if( INFO.GT.0 ) then
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        end if
    end
    
    subroutine PRINT_Eigenvalues(EigenValues, N)
        integer, intent(in) :: N
        real, dimension(1, *), intent(in) :: EigenValues
        integer :: i, j
        write(*,*) 'Eigenvalues'
        do i = 1, 1
            write(*,2) (EigenValues(i,j), j = 1, N)
        end do
        2 FORMAT( 11(:,1X,F15.6) )
    end
    
    
    subroutine PRINT_Eigenvectors(EigenVectors, N)
        integer, intent(in) :: N
        complex, dimension(N,*), intent(in) ::Eigenvectors
        integer :: i, j
        write(*,*) 'Eigenvectors'
        do i = 1, N
            write(*,1) (EigenVectors(i, j), j = 1, N )
        end do
        1 FORMAT(11(:,1X,'(',F15.6,',',F15.6,')'))
    end

  

    end
