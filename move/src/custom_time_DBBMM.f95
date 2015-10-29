!   
!   BBMM
!
!   This file contains routines to build a DLL which is callable from
!   R to fit the Brownian Bridge Movement Model.
!
!   Created by Ryan Nielson.
!
!   To compile with gfortran:
!
!   gfortran -shared -o bbmm.dll bbmm.f95 
!

! ---------------------------------------------------------------------------
subroutine custom_time_dBBMM(nLocs, gridSize, timeDiff, tTotal, X, Y, BMvar, LocationError, gridX, gridY, &
    timeStep, probability, startTime, endTime)
!
! Brownian Bridge Movement Model (does not include motion variance)
!
! Ryan Nielson
!

implicit none

! Input variables
integer :: nLocs, gridSize
double precision, dimension(nLocs) :: timeDiff
double precision :: tTotal, startTime, endTime
double precision, dimension(nLocs) :: X       
double precision, dimension(nLocs) :: Y       
double precision, dimension(nLocs) :: BMvar          
double precision, dimension(nLocs) :: LocationError  
double precision, dimension(gridSize) :: gridX
double precision, dimension(gridSize) :: gridY
double precision :: timeStep       

!character (len=180) MSG
! Local variables
double precision :: tm, alpha, muX, muY, sigma2, thisEnd
double precision, dimension(gridSize) :: int, theta, ZTZ
integer :: i               

! Output variables
double precision, dimension(gridSize) :: probability
int = 0.0
alpha = 0.0
muX = 0.0
muY = 0.0
sigma2 = 0.0

do i = 1, nLocs-1
    theta = 0.0                            
    ZTZ = 0.0
    if (i == 1) then
        tm = startTime 
    else 
        tm = 0.0
    end if
    if (i == (nLocs-1)) then
        thisEnd = endTime
    else 
        thisEnd = timeDiff(i)
    end if

    do while(tm <= thisEnd)
        call rchkusr()
        alpha = tm / timeDiff(i)
        muX = X(i) + alpha*(X(i+1) - X(i))
        muY = Y(i) + alpha*(Y(i+1) - Y(i))
        sigma2 = timeDiff(i)*alpha*(1-alpha)*BMvar(i) + &
                 ((1-alpha)**2)*(LocationError(i)**2) + &
                 (alpha**2)*(LocationError(i+1)**2)
        ZTZ = (gridX - muX)**2 + (gridY - muY)**2
        theta = (1/(2*3.14*sigma2))*exp(-ZTZ/(2*sigma2)) 
        int = int + theta
        tm = tm + timeStep
    end do
end do

!Scaling probabilities so they sum to 1.0
probability = int/tTotal
probability = probability/sum(probability)

end subroutine
