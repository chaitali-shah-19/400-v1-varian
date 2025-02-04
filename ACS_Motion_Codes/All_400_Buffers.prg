#/ Controller version = 3.03
#/ Date = 7/31/2023 5:12 PM
#/ User remarks = All Buffers run on the 400
#0
! Auto-home script
int iAxis
iAxis=0

ENABLE iAxis
WAIT 1000

VEL(iAxis)=5
ACC(iAxis)=50
DEC(iAxis)=50
JERK(iAxis)=500

FDEF(iAxis).#LL=0									! Disable the axis left limit default response.
FDEF(iAxis).#RL=0									! Disable the axis right limit default response.
JOG (iAxis),+										! Move to the left limit switch.
TILL FAULT(iAxis).#RL									! Wait for the left limit switch activation.
HALT (iAxis)										! Stop motor.
TILL ^MST(iAxis).#MOVE									! Wait until motor stop moving.
JOG (iAxis),-										! Move to the encoder index.
TILL ^FAULT(iAxis).#RL									! Wait for the left limit release.
HALT (iAxis)										! Stop motor.
TILL ^MST(iAxis).#MOVE									! Wait until motor stop moving.
SET FPOS(iAxis)=0							! Set axis origin to the position of the index.
PTP/e (iAxis),-50										! Move to the origin.
FDEF(iAxis).#LL=1									! Enable the axis left limit default response
FDEF(iAxis).#RL=1									! Enable the axis right limit default response

STOP
#1
REAL POS_GROUND;
REAL WAIT_TIME;
REAL SHUTTLE_DELAY;

POS_GROUND = V0; 	!set "start" position
WAIT_TIME = V1; 	!time to wait
SHUTTLE_DELAY = V2;	!time for shuttling delay

TILL IN0.3= 1		!Wait for motion trigger

CALL G_TO_C_TO_G	!executes motion step 1 (detailed below)
STOP

G_TO_C_TO_G:
WAIT WAIT_TIME !WAIt for laser and MW irradiation
GO 0 !Start moving, position programmed by MATLAB
TILL MST(0).#INPOS !Wait till motion finishes
WAIT SHUTTLE_DELAY !Wait for spins to diffuse
OUT(0).2=1 ! Trigger NMR
WAIT 100
OUT(0).2=0
!WAIT 3000
!WAIT 60000 !For 30s exp
!WAIT 210000 !For 3 min exp
WAIT 90000 !For 70s exp
!WAIT 320000 !For 5 min exp
!WAIT 660000 !For 10 min exp
!WAIT 200000

VEL(0) = 5
ACC(0) = 50
DEC(0) = 50
JERK(0) = 500

PTP/w 0, POS_GROUND  ! Return to start position
GO 0

RET			!Returns to call statement, continues program execution.

#2
DISPCH=-2
ON ^IN (0).0		!Triggered on input 0, change number out of parentheses for different pin.   Ex. IN (0).1 for pin 1
	DISP "This was pressed."
	KILL ALL	!Kills motion on all axes
RET			!End of autoroutine
#3
REAL POS_GROUND;
REAL WAIT_TIME;

POS_GROUND = V0; 	!set "start" position
WAIT_TIME = V1; 	!time to wait

TILL IN0.3= 1		!Wait for motion trigger

CALL G_TO_C_TO_G	!executes motion step 1 (detailed below)
STOP

G_TO_C_TO_G:
WAIT WAIT_TIME !WAIt for laser and MW irradiation
GO 0 !Start moving, position programmed by MATLAB
TILL MST(0).#INPOS !Wait till motion finishes
OUT(0).2=1 ! Trigger NMR
WAIT 100
OUT(0).2=0
!WAIT 3000
!WAIT 10000
WAIT 320000

VEL(0) = 5
ACC(0) = 50
DEC(0) = 50
JERK(0) = 500

PTP/w 0, POS_GROUND  ! Return to start position
GO 0

RET			!Returns to call statement, continues program execution.

#4
INT POS_COIL;
INT POS_GROUND;
INT POS_HOME;
INT xx;
INT yy;
INT zz;
INT vel;

POS_HOME = -10;		!set start and end position
POS_COIL = -200;		!set coil position for NMR
POS_GROUND = -1550; 	!set nitrogen gun position

!TILL IN0.3= 1
xx = TIME
!WAIT 2000

!LOOP 2

CALL HOME_TO_COIL	!executes motion step 1 (detailed below)
!WAIT 2000
!PTP/w 0, POS_GROUND
!GO 0
!TILL MST(0).#INPOS
WAIT 2000
!END

TILL IN0.5 = 1		
CALL COIL_TO_GROUND	!executes motion step 2 

!TILL IN0.6 = 1	
!CALL GROUND_TO_COIL	!executes motion step 3 

TILL IN0.7 = 1		
!CALL COIL_TO_HOME	!executes motion step 4 


HOME_TO_COIL:
vel = 10
VEL(0) = vel			!Set velocity, acceleration, ?deceleration,? and jerk. (10, 100, 1000) 		
ACC(0) = 10* vel		!parameters were recommended by the company, but we will test 
DEC(0) = 10* vel			!what works best for us
JERK(0) = 100*vel	
PTP/w 0, POS_COIL
GO 0
TILL MST(0).#INPOS
yy = TIME
zz=(yy - xx)
!DISP "elapsed time: ", zz
OUT(0).2=1
WAIT 10
OUT(0).2=0

RET			!Returns to call statement, continues program execution.

COIL_TO_GROUND:		!analogous to above call.
VEL(0) = 2000
ACC(0) =30000
DEC(0) = 30000
JERK(0) = 200*vel
PTP/w 0, POS_GROUND
GO 0
TILL MST(0).#INPOS
OUT(0).2=1
WAIT 500
OUT(0).2=0
RET

GROUND_TO_COIL:
VEL(0) = 5
ACC(0) = 100
DEC(0) = 100
JERK(0) = 1000
PTP/w 0, POS_COIL
GO 0
TILL MST(0).#INPOS
OUT(0).2=1
WAIT 500
OUT(0).2=0
RET

COIL_TO_HOME:
VEL(0) = 10
ACC(0) = 100
DEC(0) = 100
JERK(0) = 1000
PTP/w 0, POS_HOME
GO 0
TILL MST(0).#INPOS
OUT(0).2=1
WAIT 500
OUT(0).2=0
RET
#5
REAL POS_GROUND;
REAL WAIT_TIME;
REAL DELTA_TIME;

POS_GROUND = V0; 	!set "start" position
WAIT_TIME = V1; 	!time to wait
DELTA_TIME = V8;

TILL IN0.3= 1		!Wait for motion trigger
CALL G_TO_C_TO_G	!executes motion step 1 (detailed below)

TILL IN0.3= 1		!Wait for motion trigger
WAIT DELTA_TIME
PTP/w 0, -617  ! Return to start position
CALL G_TO_C_TO_G	!executes motion step 1 (detailed below)

STOP

G_TO_C_TO_G:
WAIT WAIT_TIME !WAIt for laser and MW irradiation
GO 0 !Start moving, position programmed by MATLAB
TILL MST(0).#INPOS !Wait till motion finishes
OUT(0).2=1 ! Trigger NMR
WAIT 100
OUT(0).2=0
!WAIT 3000
!WAIT 10000
WAIT 10000

PTP/w 0, POS_GROUND  ! Return to start position
GO 0

RET			!Returns to call statement, continues program execution.
#6
REAL POS_GROUND;
REAL WAIT_TIME;
REAL POS_TIME;

POS_GROUND = V0; 	!set "start" position
WAIT_TIME = V1; 	!time to wait
POS_TIME = V2;		!time for relaxation

TILL IN0.3= 1		!Wait for motion trigger

CALL G_TO_C_TO_G	!executes motion step 1 (detailed below)
STOP

G_TO_C_TO_G:
WAIT WAIT_TIME !WAIt for laser and MW irradiation
GO 0 !Start moving, position programmed by MATLAB
TILL MST(0).#INPOS !Wait till motion finishes
WAIT POS_TIME !Wait for pos_time 
OUT(0).2=1 ! Trigger NMR
WAIT 100
OUT(0).2=0
WAIT 10000


PTP/w 0, POS_GROUND  ! Return to start position
GO 0

RET			!Returns to call statement, continues program execution.
#7
REAL POS_GROUND;
REAL WAIT_TIME;

POS_GROUND = V0; 	!set "start" position
WAIT_TIME = V1; 	!time to wait

!TILL IN0.3= 1		!Wait for motion trigger

CALL G_TO_C_TO_G	!executes motion step 1 (detailed below)
STOP

G_TO_C_TO_G:
WAIT WAIT_TIME !WAIt for thermal

OUT(0).2=1 ! Trigger NMR
WAIT 100
OUT(0).2=0
!WAIT 3000
WAIT 10000


RET			!Returns to call statement, continues program execution.

#8
REAL POS_GROUND;
REAL WAIT_TIME;
REAL POS_TIME;

POS_GROUND = V0; 	!set "start" position
WAIT_TIME = V1; 	!time to wait
POS_TIME = V2;		!time for relaxation

TILL IN0.3= 1		!Wait for motion trigger

CALL G_TO_C_TO_G	!executes motion step 1 (detailed below)
STOP

G_TO_C_TO_G:
WAIT WAIT_TIME !WAIt for laser and MW irradiation
GO 0 !Start moving, position programmed by MATLAB
TILL MST(0).#INPOS !Wait till motion finishes
OUT(0).2=1 ! Trigger NMR
WAIT 100
OUT(0).2=0

WAIT 5000 !Wait
PTP/w 0, -1300 !Go to NMR Coil position
GO 0 
TILL MST(0).#INPOS !Wait till motion finishes
WAIT 3000

VEL(0) = 5
ACC(0) = 50
DEC(0) = 50
JERK(0) = 500
PTP/w 0, POS_GROUND  ! Return to start position
GO 0

RET			!Returns to call statement, continues program execution.
#9
REAL POS_GROUND;
REAL WAIT_TIME;
REAL POS_TIME;

POS_GROUND = V0; 	!set "start" position
WAIT_TIME = V1; 	!time to wait
POS_TIME = V2;		!time for relaxation

TILL IN0.3= 1		!Wait for motion trigger

CALL G_TO_C_TO_G	!executes motion step 1 (detailed below)
STOP

G_TO_C_TO_G:
WAIT WAIT_TIME !WAIt for laser and MW irradiation
OUT(0).2=1 ! Trigger NMR
WAIT 100
OUT(0).2=0

RET			!Returns to call statement, continues program execution.
#A
!axisdef X=0,Y=1,Z=2,T=3,A=4,B=5,C=6,D=7
!axisdef x=0,y=1,z=2,t=3,a=4,b=5,c=6,d=7
global int I(100),I0,I1,I2,I3,I4,I5,I6,I7,I8,I9,I90,I91,I92,I93,I94,I95,I96,I97,I98,I99
global real V(100),V0,V1,V2,V3,V4,V5,V6,V7,V8,V9,V90,V91,V92,V93,V94,V95,V96,V97,V98,V99
