!*********Photoluminescence kinetics with Mn layer

program plkin
integer(4) i,it,ntot,N,i1s,ip,i2s
parameter(N=6,ntot=150000, nfin=1500, nfin2=50000)
real A,B,Cm,D,s1,s2,ttot,t,h,t1s,t1p,t2s,t2p,gammae,gammah,gammamn,Ne0,Nh0,Nmn, kh1,kh2,kh3,kh4
real C(N,ntot+1),k1(N),k2(N),k3(N),k4(N),Iplus(ntot+1),Iminus(ntot+1),P, IIm, IIp, Ff(nfin,nfin2), Nhsum(ntot+1), f1(ntot+1),f2(ntot+1)
real, parameter :: pi = 3.141592653589
!******************************************

!*********** Rates and other parameters************************************************


 !A=0.1		!*** laser field excitation rate for electrons and holes in 1/ps
 !B=0.0125		!*** recombination (or luminescence decay) rate for electrons and holes in 1/ps
 !Cm=1.9E-3         !*** hole polarization rate due to interaction with Mn layer
 !D=0.2          !*** Mn layer polarization rate due to interaction with holes      
 !s1=1.0         !*** first pulse laser field amplitude
 !s2=1.0         !*** second pulse laser field amplitude
 !ttot=1500.0    !*** total time of simulation in ps
 !h=ttot/ntot    !*** time step interval in ps
 !t1s=20.0      !*** first laser pulse start in ps
 !t1p=0.1        !*** first laser pulse duration in ps
 !t2s=520.0      !*** second laser pulse start in ps
 !t2p=0.1        !*** second laser pulse duration in ps
 !delayt=10.0      !*** camera delay time of writing the intensities
 !
 !gammae=1.923E-1 !*** electron spin relaxation rate in 1/ps
 !gammae=0.2/2. !*** electron spin relaxation rate in 1/ps
 !gammah=0.2     !*** hole spin relaxation rate in 1/ps
 !gammah=0.2/2.   !*** hole spin relaxation rate in 1/ps
 !
 !gammamn=2.5E-4 !*** Mn layer spin relaxation rate in 1/ps
 !Ne0=8.0E-3        !*** initial concentration of electrons in QW 
 !Nh0=0.0		!*** initial hole concentration in QW
 !Nmn=1.0        !*** total concentration of Mn spins
 !Nstrobe=ntot/ttot !** stroboscopic number of points for averaging
 !deltat=50.0    !*** camera window for convolution in ps 
 i1s=2000       !*** starting point for the first pulse
 i2s=52000      !*** starting point for the second pulse
 ip=10  
 !i1s= t1s*int(ntot/ttot) ! 2000       !*** starting point for the first pulse
 !i2s= t2s*int(ntot/ttot) ! 52000      !*** starting point for the second pulse
 !ip=t2p*int(ntot/ttot) ! 10          !*** each pulse duration in points


 A=0.2		!*** laser field excitation rate for electrons and holes in 1/ps
 B=0.0125	!*** recombination (or luminescence decay) rate for electrons and holes in 1/ps
 Cm=1.0E-3         !*** hole polarization rate due to interaction with Mn layer
 D=0.2          !*** Mn layer polarization rate due to interaction with holes      
 s1=1.0         !*** first pulse laser field amplitude
 s2=1.0         !*** second pulse laser field amplitude
 ttot=1500.0    !*** total time of simulation in ps
 h=ttot/ntot    !*** time step interval in ps
 t1s=20.0      !*** first laser pulse start in ps
 t1p=0.1        !*** first laser pulse duration in ps
 t2s=520.0      !*** second laser pulse start in ps
 t2p=0.1        !*** second laser pulse duration in ps
 gammae=1.923E-3 !*** electron spin relaxation rate in 1/ps
 gammah=0.2     !*** hole spin relaxation rate in 1/ps
 gammamn=2.0E-4 !*** Mn layer spin relaxation rate in 1/ps
 Ne0=9.0E-3        !*** initial concentration of electrons in QW 8.0E-3  
 Nh0=0.0		!*** initial hole concentration in QW
 Nmn=1.0        !*** total concentration of Mn spins
 !Nstrobe=1500 !** stroboscopic number of points for averaging
 deltat=50.0    !*** camera window for convolution in ps 
 
 
!***************************************************************************************
  
!open(unit=1,status='unknown',file='Ne1.DAT')
!open(unit=2,status='unknown',file='Ne2.DAT')
!open(unit=3,status='unknown',file='Nh1.DAT')
!open(unit=4,status='unknown',file='Nh2.DAT')

open(unit=7,status='unknown',file='Iminus2.DAT')
open(unit=8,status='unknown',file='Iplus2.DAT')
open(unit=9,status='unknown',file='Polariz2.DAT')
open(unit=10,status='unknown',file='Iminus1.DAT')
open(unit=11,status='unknown',file='Iplus1.DAT')
open(unit=12,status='unknown',file='Polariz1.DAT')
!open(unit=10,status='unknown',file='fminus.DAT')
!open(unit=11,status='unknown',file='fplus.DAT')

!********************************************************************************************
!*** C(1,it) is the spin-down electron concentration in the QW at the discrete moment of time it
!*** C(2,it) is the spin-up electron concentration in the QW
!*** C(3,it) is the spin-down hole concentration in the QW
!*** C(4,it) is the spin-up hole concentration in the QW
!*** C(5,it) is the spin-down Mn spin concentration in the delta layer 
!*** C(6,it) is the spin-up Mn spin concentration in the delta layer 



!************ Runge-Kutta 4th order method *****************************
!***** index i for k(i) is the number of kinetic equation, i=1...6 *******
  
  
call RungeKutta(Iminus,Iplus, .true.)
do it = 1, nfin
!        print *,"Tooooooooooooooooooooooooy", it
        call integralSvertka(it, Iminus, IIm)
        call integralSvertka(it, Iplus, IIp)
        WRITE(7,220) IIm
        WRITE(8,220) IIp
        IF (IIp+IIm) 4100,4000,4100
        4000 P=0.0
        GOTO 4200
        4100 P=(IIp-IIm)/(IIp+IIm)
        4200 CONTINUE
        !P=(Ip-Im)/(Ip+Im)
        WRITE(9,220) abs(P)
        print*, "Progress", int(100*it/nfin)
end do
!call RungeKutta(Iminus,Iplus, .true.)
!do it = 1, ntot
!    IF (Iplus(it)+Iminus(it)) 410,400,410
!    400 P=0.0
!    GOTO 420
!    410 P=(Iplus(it)-Iminus(it))/(Iplus(it)+Iminus(it))
!    420 CONTINUE
!   IF (mod(it, nfin) == 0) then
!        WRITE(7,220) Iminus(it)
!        WRITE(8,220) Iplus(it)
!        WRITE(9,220) abs(P)
!        print*, "Progress", int(100*it/ntot)   
!   endif
!end do
!call RungeKutta(Iminus,Iplus, .false.)
!do it = 1, ntot
!    IF (Iplus(it)+Iminus(it)) 4100,4000,4100
!    4000 P=0.0
!    GOTO 4200
!    4100 P=(Iplus(it)-Iminus(it))/(Iplus(it)+Iminus(it))
!    4200     CONTINUE
!   IF (mod(it, nfin) == 0) then
!        WRITE(10,220) Iminus(it)
!        WRITE(11,220) Iplus(it)
!        WRITE(12,220) abs(P)
!        print*, "Progress", int(100*it/ntot)   
!   endif
!end do

call RungeKutta(Iminus,Iplus, .false.)
do it = 1, nfin
        call integralSvertka(it, Iminus, IIm)
        call integralSvertka(it, Iplus, IIp)
        WRITE(10,220) IIm
        WRITE(11,220) IIp
        IF (IIp+IIm) 410,400,410
        400 P=0.0
        GOTO 420
        410 P=(IIp-IIm)/(IIp+IIm)
        420 CONTINUE
        WRITE(12,220) abs(P)
        print*, "Progress", int(100*it/nfin)
end do


close(7)
close(8)
close(9)
close(10)
close(11)
close(12)
220 Format (F18.8)
  
pause
contains

    subroutine integralSvertka(N, II, Sum)
        real, intent(out) :: Sum 
        real, dimension(ntot), intent(in) :: II
        integer, intent(in) :: N
        
        integer :: i
        real :: dt
        !real, dimension(N, N) :: Ff
        Sum = 0.
        !dt = h*ntot/Nstrobe!((t1 - t0)/N)*ntot/Nstrobe
        dt =h*ntot/nfin2
        !t=N*h*ntot/Nstrobe
        t=h*ntot/nfin*N
        !print *, "Sum strat", Sum
        DO i=1, INT(N*nfin2/nfin)
            Ff(N, i) = (1/sqrt(2.*pi*deltat*deltat)) * exp(-(t - i*dt)*(t - i*dt)/(2. * deltat * deltat))
            Sum = Sum + Ff(N,i)*II(INT(ntot/nfin2*i))
            !print *, "Sum", Sum
        end do
    end subroutine integralSvertka
        !if (MOD(int(N/Nstrobe), 3) == 0) then
        !
        !    do i = 1, int(N/Nstrobe)
        !        Ff(i) = (1/sqrt(2*pi*deltat*deltat)) * exp(-(t1 - i*dt)*(t1 - i*dt)/(2 * deltat * deltat))*II(i*ntot/Nstrobe)
        !    enddo
        !    Sum = 0.
        !    !do i = 0, int((Nstrobe-3)/2)
        !    do i = 0, int((N/Nstrobe-3)/2)
        !        Sum = Sum + dt/3 * (Ff(2*i+1) + 4*Ff(2*i + 2) + Ff(2*i + 3))
        !        print *, "Sum", Sum
        !    enddo
        !else
        !    if (int(N/Nstrobe)<3) then
        !        Sum = 0. 
        !    endif
        !endif
        !!print *, Sum

!*********** Cycle over time t, t is in ps **********************
!***************** "it" is the discrete time variable ***********

subroutine RungeKutta(Iplus, Iminus, flag)
    real, dimension(ntot+1), intent(out) :: Iplus
    real, dimension(ntot+1), intent(out) :: Iminus
        logical, intent(in) :: flag
        
    C(1,1)=Ne0/2.0     !*** initial concentration of spin-down electrons in the QW
    C(2,1)=Ne0/2.0     !*** initial concentration of spin-up electrons in the QW
    C(3,1)=0.0         !*** initial concentration of spin-down holes in the QW
    C(4,1)=0.0         !*** initial concentration of spin-up holes in the QW
    C(5,1)=Nmn/2.0     !*** initial concentration of spin-down Mn states in the delta layer
    C(6,1)=Nmn/2.0     !*** initial concentration of spin-up Mn states in the delta layer
    Nhsum(1)=C(3,1)+C(4,1)

        
    do it=1, ntot 
        t=h*it

        !!*** frames for sigma_minus and sigma_plus pulses, f1_2=1 only in window [ts...ts+tp] ********************
        !IF(t-t1s) 100,110,110
        !
        !100 f1(it)=0.0
        !GOTO 200
        !
        !110 IF(t-(t1s+t1p)) 120,120,130
        !!120 f1=EXP(-((t-(t1s+t1p/2.0))/(t1p/2.0))**2.0)
        ! 120 f1(it)=1.0
        !GOTO 200
        !
        !130 f1(it)=0.0
        !
        !
        !200 IF (t-t2s) 210,215,215
        !
        !210 f2(it)=0.0
        !GOTO 300
        !
        !215 IF(t-(t2s+t2p)) 230,230,240
        !230 f2(it)=1.0
        !GOTO 300
        !
        !240 f2(it)=0.0
        !
        !300     CONTINUE
        IF(it-i1s) 100,110,110

        100 f1(it)=0.0
        GOTO 200

        110 IF(it-(i1s+ip)) 120,120,130
        !120 f1=EXP(-((t-(t1s+t1p/2.0))/(t1p/2.0))**2.0)
         120 f1(it)=1.0
        GOTO 200

        130 f1(it)=0.0


        200 IF (it-i2s) 210,215,215

        210 f2(it)=0.0
        GOTO 300

        215 IF(it-(i2s+ip)) 230,230,240
        230 f2(it)=1.0
        GOTO 300

        240 f2(it)=0.0

        300 CONTINUE
        
    enddo
    
    do it = 1, ntot
    t=h*it

    if (flag == .true.) then
        open(unit=5,status='unknown',file='NMn1.DAT')
        open(unit=6,status='unknown',file='NMn2.DAT')
    !****** For the first pulse  with sigma_minus and the second pulse with sigma_plus

        sminus=f1(it)*s1
        sminus1=(f1(it)+f1(it+1))/2.0*s1
        sminus2=f1(it+1)*s1
        !**** sigma_plus pulse time dependence for three time points t, t+h/2, t+h **  
        splus=f2(it)*s2
        splus1=(f2(it)+f2(it+1))/2.0*s2
        splus2=f2(it+1)*s2

    else
    !****** For the first pulse  with sigma_minus and the second pulse with sigma_minus
        open(unit=5,status='unknown',file='NMn1.DAT')
        open(unit=6,status='unknown',file='NMn2.DAT')

        !splus1=0.0
        !splus2=0.0
        !
        !IF(t-(t1s+t1p)) 150,150,160
        !
        !150  sminus=f1(it)*s1
        !     sminus1=(f1(it)+f1(it+1))/2.0*s1
        !     sminus2=f1(it+1)*s1
        !GOTO 170
        !
        !160  sminus=f2(it)*s1
        !     sminus1=(f2(it)+f2(it+1))/2.0*s1
        !     sminus2=f2(it+1)*s1
        !
        !170      CONTINUE
        !    endif
        splus=0.0
        splus1=0.0
        splus2=0.0

        IF(it-(i1s+ip)) 150,150,160

        150  sminus=f1(it)*s1
             sminus1=(f1(it)+f1(it+1))/2.0*s1
             sminus2=f1(it+1)*s1
        GOTO 170

        160  sminus=f2(it)*s1
             sminus1=(f2(it)+f2(it+1))/2.0*s1
             sminus2=f2(it+1)*s1

        170      CONTINUE
    endif    

    !WRITE(10,220) sminus
    !WRITE(11,220) splus

    !**********************************************************************

    k1(1)=h*(A*splus-B*MIN(C(1,it),C(4,it))-gammae*(C(1,it)-C(2,it)))
    k1(2)=h*(A*sminus-B*MIN(C(2,it),C(3,it))-gammae*(C(2,it)-C(1,it)))
    k1(3)=h*(A*sminus+Cm*(C(5,it)-C(6,it))-B*MIN(C(2,it),C(3,it))-gammah*(C(3,it)-C(4,it)))
    k1(4)=h*(A*splus+Cm*(C(6,it)-C(5,it))-B*MIN(C(1,it),C(4,it))-gammah*(C(4,it)-C(3,it)))
    k1(5)=h*(D*(C(3,it)-C(4,it))-gammamn*(C(5,it)-C(6,it)))
    k1(6)=-k1(5)
    kh1=h*(A*(splus+sminus)-B*(MIN(C(1,it),C(4,it))+MIN(C(2,it),C(3,it))))

    k2(1)=h*(A*splus1-B*MIN(C(1,it)+k1(1)/2.0,C(4,it)+k1(4)/2.0)-gammae*(C(1,it)+k1(1)/2.0-(C(2,it)+k1(2)/2.0)))
    k2(2)=h*(A*sminus1-B*MIN(C(2,it)+k1(2)/2.0,C(3,it)+k1(3)/2.0)-gammae*(C(2,it)+k1(2)/2.0-(C(1,it)+k1(1)/2.0)))
    k2(3)=h*(A*sminus1+Cm*(C(5,it)+k1(5)/2.0-(C(6,it)+k1(6)/2.0))-B*MIN(C(2,it)+k1(2)/2.0,C(3,it)+k1(3)/2.0)-gammah*(C(3,it)+k1(3)/2.0-(C(4,it)+k1(4)/2.0)))
    k2(4)=h*(A*splus1+Cm*(C(6,it)+k1(6)/2.0-(C(5,it)+k1(5)/2.0))-B*MIN(C(1,it)+k1(1)/2.0,C(4,it)+k1(4)/2.0)-gammah*(C(4,it)+k1(4)/2.0-(C(3,it)+k1(3)/2.0)))
    k2(5)=h*(D*(C(3,it)+k1(3)/2.0-(C(4,it)+k1(4)/2.0))-gammamn*(C(5,it)+k1(5)/2.0-(C(6,it)+k1(6)/2.0)))
    k2(6)=-k2(5)
    kh2=h*(A*(splus1+sminus1)-B*(MIN(C(1,it)+k1(1)/2.0,C(4,it)+k1(4)/2.0)+MIN(C(2,it)+k1(2)/2.0,C(3,it)+k1(3)/2.0)))

    k3(1)=h*(A*splus1-B*MIN(C(1,it)+k2(1)/2.0,C(4,it)+k2(4)/2.0)-gammae*(C(1,it)+k2(1)/2.0-(C(2,it)+k2(2)/2.0)))
    k3(2)=h*(A*sminus1-B*MIN(C(2,it)+k2(2)/2.0,C(3,it)+k2(3)/2.0)-gammae*(C(2,it)+k2(2)/2.0-(C(1,it)+k2(1)/2.0)))
    k3(3)=h*(A*sminus1+Cm*(C(5,it)+k2(5)/2.0-(C(6,it)+k2(6)/2.0))-B*MIN(C(2,it)+k2(2)/2.0,C(3,it)+k2(3)/2.0)-gammah*(C(3,it)+k2(3)/2.0-(C(4,it)+k2(4)/2.0)))
    k3(4)=h*(A*splus1+Cm*(C(6,it)+k2(6)/2.0-(C(5,it)+k2(5)/2.0))-B*MIN(C(1,it)+k2(1)/2.0,C(4,it)+k2(4)/2.0)-gammah*(C(4,it)+k2(4)/2.0-(C(3,it)+k2(3)/2.0)))
    k3(5)=h*(D*(C(3,it)+k2(3)/2.0-(C(4,it)+k2(4)/2.0))-gammamn*(C(5,it)+k2(5)/2.0-(C(6,it)+k2(6)/2.0)))
    k3(6)=-k3(5)
    kh3=h*(A*(splus1+sminus1)-B*(MIN(C(1,it)+k2(1)/2.0,C(4,it)+k2(4)/2.0)+MIN(C(2,it)+k2(2)/2.0,C(3,it)+k2(3)/2.0)))

    k4(1)=h*(A*splus2-B*MIN(C(1,it)+k3(1),C(4,it)+k3(4))-gammae*(C(1,it)+k3(1)-(C(2,it)+k3(2))))
    k4(2)=h*(A*sminus2-B*MIN(C(2,it)+k3(2),C(3,it)+k3(3))-gammae*(C(2,it)+k3(2)-(C(1,it)+k3(1))))
    k4(3)=h*(A*sminus2+Cm*(C(5,it)+k3(5)-(C(6,it)+k3(6)))-B*MIN(C(2,it)+k3(2),C(3,it)+k3(3))-gammah*(C(3,it)+k3(3)-(C(4,it)+k3(4))))
    k4(4)=h*(A*splus2+Cm*(C(6,it)+k3(6)-(C(5,it)+k3(5)))-B*MIN(C(1,it)+k3(1),C(4,it)+k3(4))-gammah*(C(4,it)+k3(4)-(C(3,it)+k3(3))))
    k4(5)=h*(D*(C(3,it)+k3(3)-(C(4,it)+k3(4)))-gammamn*(C(5,it)+k3(5)-(C(6,it)+k3(6))))
    k4(6)=-k4(5)
    kh4=h*(A*(splus2+sminus2)-B*(MIN(C(1,it)+k3(1),C(4,it)+k3(4))+MIN(C(2,it)+k3(2),C(3,it)+k3(3))))
    

    DO i=1, N
    C(i,it+1)=C(i,it)+1.0/6.0*(k1(i)+2.0*k2(i)+2.0*k3(i)+k4(i))
    ENDDO
    Nhsum(it+1)=Nhsum(it)+1.0/6.0*(kh1+2.0*kh2+2.0*kh3+kh4)

    !******** Spin-resolved concentrations cannot exceed the full concentration chosen between 0 and 1 ********

    IF (ABS(C(1,it+1)).GT.1.0) THEN 
    C(1,it+1)=C(1,it)/ABS(C(1,it+1))
    END IF

    IF (ABS(C(2,it+1)).GT.1.0) THEN 
    C(2,it+1)=C(2,it)/ABS(C(2,it+1))
    END IF

    IF (ABS(C(3,it+1)).GT.1.0) THEN 
    C(3,it+1)=C(3,it)/ABS(C(3,it+1))
    END IF

    IF (ABS(C(4,it+1)).GT.1.0) THEN 
    C(4,it+1)=C(4,it)/ABS(C(4,it+1))
    END IF

    IF (ABS(C(5,it+1)).GT.1.0) THEN 
    C(5,it+1)=C(5,it)/ABS(C(5,it+1))
    END IF

    IF (ABS(C(6,it+1)).GT.1.0) THEN 
    C(6,it+1)=C(6,it)/ABS(C(6,it+1))
    END IF


    IF (C(1,it+1).LT.0.0) THEN 
    C(1,it+1)=0.0
    END IF

    IF (C(2,it+1).LT.0.0) THEN 
    C(2,it+1)=0.0
    END IF

    IF (C(3,it+1).LT.0.0) THEN 
    C(3,it+1)=0.0
    END IF

    IF (C(4,it+1).LT.0.0) THEN 
    C(4,it+1)=0.0
    END IF

    IF (C(5,it+1).LT.0.0) THEN 
    C(5,it+1)=0.0
    END IF

    IF (C(6,it+1).LT.0.0) THEN 
    C(6,it+1)=0.0
    END IF
    
    IF(C(3,it+1)+C(4,it+1).GT.Nhsum(it+1)) THEN
    C(3,it+1)=C(3,it+1)/(C(3,it+1)+C(4,it+1))*Nhsum(it+1)
    C(4,it+1)=C(4,it+1)/(C(3,it+1)+C(4,it+1))*Nhsum(it+1)
    END IF

    !******** Recombination intensities and polarization degree, IF for P to avoid division by 0 *****
    !if (C(1,it+1)<C(4,it+1)) then
    !    print *, "electron lum Iplus"
    !else
    !    print *, "dirki lum Iplus"
    !endif
    !if (C(2,it+1)<C(3,it+1)) then
    !    print *, "electron lum Iminus"
    !else
    !    print *, "dirki lum Iminus"
    !endif
    Iplus(it)=B*MIN(C(1,it+1),C(4,it+1))
    Iminus(it)=B*MIN(C(2,it+1),C(3,it+1))

    !IF (Iplus+Iminus) 410,400,410 
    !400 P=0.0
    !GOTO 420
    !410 P=(Iplus-Iminus)/(Iplus+Iminus)

    420 CONTINUE
    !********** Writing C(i,it+1) and polarization degree into files ************************************************
    !*** We use stroboscopic rate of writing data since the camera in experiments has the time resolution of 50 ps ***
    !***** We use the stroboscopic period of 15 ps **************************************************************

    IF ((1.0*it/100)-AINT(1.0*it/100)) 510,500,510

    !500 WRITE(1,220) C(1,it+1)
    !    WRITE(2,220) C(2,it+1)
    !    WRITE(3,220) C(3,it+1)
    !    WRITE(4,220) C(4,it+1)
    500 WRITE(5,220) C(5,it+1)
        WRITE(6,220) C(6,it+1)
        !WRITE(7,220) Iminus
        !WRITE(8,220) Iplus
        !WRITE(9,220) C(1,it+1) + C(2,it+1) - C(3,it+1) - C(4,it+1)

    510 CONTINUE

    WRITE(*,12) it, ntot

ENDDO
    !********* End of cycle over time **********************************
    !
    !close(1)
    !close(2)
    !close(3)
    !close(4)
    close(5)
    close(6)
    !close(7)
    !close(8)
    !close(9)
    !close(10)
    !close(11)
    
    220 Format (F18.8)
    12      Format (' Time step ', I6, ' of ', I6, '  ')
    end
    

!*********************************************************************
END
