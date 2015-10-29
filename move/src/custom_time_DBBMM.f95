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
    timeStep, probability, nTimes, times, llocTimes)
!
! Brownian Bridge Movement Model (does not include motion variance)
!
! Ryan Nielson
!

implicit none

! Input variables
integer :: nLocs, gridSize, nTimes
double precision, dimension(nTimes) :: times
double precision, dimension(nLocs) :: timeDiff, llocTimes
double precision :: tTotal         
double precision, dimension(nLocs) :: X       
double precision, dimension(nLocs) :: Y       
double precision, dimension(nLocs) :: BMvar          
double precision, dimension(nLocs) :: LocationError  
double precision, dimension(gridSize) :: gridX
double precision, dimension(gridSize) :: gridY
double precision :: timeStep       

!character (len=180) MSG
! Local variables
double precision :: tm, this_end_time, alpha, muX, muY, sigma2
double precision, dimension(gridSize) :: int, theta, ZTZ
integer :: i, t, first_lloc, last_lloc               

! Output variables
double precision, dimension(nTimes-1, gridSize) :: probability
int = 0.0
alpha = 0.0
muX = 0.0
muY = 0.0
sigma2 = 0.0
first_lloc = 1
last_lloc = 2

do t = 1, nTimes-1
    int = 0.0 ! resetting

    do while(llocTimes(first_lloc+1)<times(t))
            first_lloc = first_lloc + 1
    end do
    do while(llocTimes(last_lloc)<times(t+1))
        last_lloc = last_lloc + 1
    end do

    do i = first_lloc, last_lloc-1

        theta = 0.0                            
        if(i == first_lloc) then
            tm = times(t) - llocTimes(first_lloc)
        else
            tm = 0
        end if
        if(i == (last_lloc-1)) then
            this_end_time = times(t+1) - llocTimes(last_lloc-1) 
        else
            this_end_time = timeDiff(i)
        end if
        ZTZ = 0.0
        do while(tm <= this_end_time)
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
            ! TO THINK ABOUT: SHOULDN'T WE CHECK IF THERE IS A CHUNK LEFT SMALLER THAN TIMESTEP?
        end do
    end do

    probability(t, :) = int/(times(t+1)-times(t))
    if(sum(probability(t, :)) > 0) then
        probability(t, :) = probability(t, :)/sum(probability(t, :))
    end if

    first_lloc = last_lloc
    last_lloc = first_lloc + 1

end do

end subroutine