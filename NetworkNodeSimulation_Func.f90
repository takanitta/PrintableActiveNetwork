MODULE FUNC

USE PARAMETERS

IMPLICIT NONE

CONTAINS



!---Calculate Contraction Speed---------------------------------------

FUNCTION CalvContract(MovReverse,F)

!$acc routine seq

	USE PARAMETERS, ONLY : REAL_PRECISION, vContractMax, F_MotorStall


	INTEGER, VALUE, INTENT(IN) :: MovReverse
	REAL(KIND = REAL_PRECISION) :: CalvContract
	REAL(KIND = REAL_PRECISION), INTENT(IN) :: F
	


	IF (MovReverse == 1) THEN
		CalvContract = -vContractMax*(F_MotorStall + F)/F_MotorStall

		IF (F < -F_MotorStall) CalvContract = 0.0_REAL_PRECISION

		IF (CalvContract < -vContractMax) CalvContract = -vContractMax
	ELSE
		CalvContract = vContractMax*(F_MotorStall - F)/F_MotorStall

		IF (F > F_MotorStall) CalvContract = 0.0_REAL_PRECISION

		IF (CalvContract > vContractMax) CalvContract = vContractMax
	END IF

END FUNCTION CalvContract





!---Subroutine InitialDistributionCircular----------------------------
!Set Initial Distribution of Nodes within Circular Region
!---------------------------------------------------------------------

SUBROUTINE InitialDistributionCircular(X,Y)

	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode, Radius
	USE mtmod


	INTEGER :: Counter
	REAL(KIND = REAL_PRECISION) :: UR1, UR2, Xtemp, Ytemp
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode), INTENT(OUT) :: X,Y



	Counter = 0
	DO

		UR1 = grnd()
		UR2 = grnd()

		Xtemp = (UR1 - 0.5_REAL_PRECISION)*2.0_REAL_PRECISION*Radius
		Ytemp = (UR2 - 0.5_REAL_PRECISION)*2.0_REAL_PRECISION*Radius

		IF (Xtemp*Xtemp + Ytemp*Ytemp < Radius*Radius) THEN
			Counter = Counter + 1
			X(Counter) = Xtemp
			Y(Counter) = Ytemp
		END IF

		IF (Counter == NumNode) EXIT

	END DO

END SUBROUTINE InitialDistributionCircular





!---Subroutine InitialDistributionRectangular-------------------------
!Set Initial Distribution of Nodes within Rectangular Region
!---------------------------------------------------------------------

SUBROUTINE InitialDistributionRectangular(X,Y)

	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode, X_Size, Y_Size
	USE mtmod


	INTEGER :: I

	REAL(KIND = REAL_PRECISION) :: UR1, UR2
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode), INTENT(OUT) :: X,Y



	DO I=1,NumNode

		UR1 = grnd()
		UR2 = grnd()

		X(I) = (UR1 - 0.5_REAL_PRECISION) * X_Size
		Y(I) = (UR2 - 0.5_REAL_PRECISION) * Y_Size

	END DO

END SUBROUTINE InitialDistributionRectangular





!---Subroutine InitialDistributionRectangularGrid---------------------
!Set Initial Distribution of Nodes on Grid
!---------------------------------------------------------------------

SUBROUTINE InitialDistributionRectangularGrid(X,Y)

	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode, X_Size, Y_Size, GridSeparation
	USE mtmod


	INTEGER :: I, J, Nx, Ny, Counter
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode), INTENT(OUT) :: X,Y



	Nx = INT(X_Size/GridSeparation)
	Ny = INT(Y_Size/GridSeparation)


	Counter = 0
	DO I=1,Nx
		DO J=1,Ny
			Counter = Counter + 1
			X(Counter) = REAL(I-1)*GridSeparation
			Y(Counter) = REAL(J-1)*GridSeparation
		END DO
	END DO

	X = X - REAL(Nx-1)/2.0_REAL_PRECISION * GridSeparation
	Y = Y - REAL(Ny-1)/2.0_REAL_PRECISION * GridSeparation

END SUBROUTINE InitialDistributionRectangularGrid





!---Subroutine AssignSortNodeTypes------------------------------------
!Assign Node Types and Sort by Node Types
!---------------------------------------------------------------------

SUBROUTINE AssignNodeTypesSort(X,Y,NodeRadius,NodeType,NumMovableNode)



	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode
	USE mtmod


	INTEGER :: I, CounterMovable, CounterAnchor
	REAL(KIND = REAL_PRECISION) :: UR1, UR2
	INTEGER, INTENT(OUT) :: NumMovableNode
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode), INTENT(INOUT) :: X, Y, NodeRadius
	INTEGER, DIMENSION(NumNode), INTENT(OUT) :: NodeType
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode) :: Xtemp, Ytemp, NodeRadiusTemp




	DO I=1,NumNode
		Xtemp(I) = X(I)
		Ytemp(I) = Y(I)
		NodeRadiusTemp(I) = NodeRadius(I)
	END DO

	CounterMovable = 0
	CounterAnchor = 0


	DO I=1,NumNode
		IF ((Xtemp(I) >= -40.0_REAL_PRECISION) .AND. (Xtemp(I) <= 40.0_REAL_PRECISION)) THEN
			CounterMovable = CounterMovable + 1
			X(CounterMovable) = Xtemp(I)
			Y(CounterMovable) = Ytemp(I)
			NodeRadius(CounterMovable) = NodeRadiusTemp(I)
			NodeType(CounterMovable) = 0
		ELSE
			CounterAnchor = CounterAnchor + 1
			X(NumNode - CounterAnchor + 1) = Xtemp(I)
			Y(NumNode - CounterAnchor + 1) = Ytemp(I)
			NodeRadius(NumNode - CounterAnchor + 1) = NodeRadiusTemp(I)
			NodeType(NumNode - CounterAnchor + 1) = 1
		END IF

	END DO
	NumMovableNode = CounterMovable

END SUBROUTINE AssignNodeTypesSort











SUBROUTINE CalStats(BondState, RIJ, BondNatLength, NumMotorOnBond, ElasticE, NumContractingBond, NumContractedBond, NumBond)

	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode, k


	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode,NumNode), INTENT(IN) :: RIJ, BondNatLength, NumMotorOnBond
	INTEGER, DIMENSION(NumNode,NumNode), INTENT(IN) :: BondState
	INTEGER :: I, J
	INTEGER, INTENT(OUT) :: NumContractingBond, NumContractedBond, NumBond			! NumBoundNode, NumBond, BoundNodeIndexCounterJ
	REAL(KIND = REAL_PRECISION), INTENT(OUT) :: ElasticE



!---------------------------------------------------------------------



	NumContractingBond = 0
	NumContractedBond = 0
	NumBond = 0

	ElasticE = 0.0_REAL_PRECISION


	DO I=1,NumNode
		DO J=I+1,NumNode
			IF (BondState(I,J)==1) THEN

				ElasticE = ElasticE + 0.5_REAL_PRECISION * k * NumMotorOnBond(I,J) * (RIJ(I,J)-BondNatLength(I,J))*(RIJ(I,J)-BondNatLength(I,J))

				NumBond = NumBond + 1

				IF (BondNatLength(I,J) > BondNatLengthMin) THEN
					NumContractingBond = NumContractingBond + 1
				ELSE
					NumContractedBond = NumContractedBond + 1
				END IF

			END IF
		END DO
	END DO

END SUBROUTINE CalStats





!---Subroutine NormRNDVector-----------------------------------------
!Generate N-Dimensinal Normal Random Vector
!---------------------------------------------------------------------

SUBROUTINE NormRNDVector(N_Dim, NRV)



	USE PARAMETERS, ONLY : REAL_PRECISION
	USE mtmod, ONLY : grnd

	INTEGER, INTENT(IN) :: N_Dim
	REAL(KIND = REAL_PRECISION), DIMENSION(N_Dim), INTENT(OUT) :: NRV

	INTEGER :: I
	REAL(KIND = REAL_PRECISION) :: UR1, UR2, NR



	DO I=1, N_Dim

		UR1 = grnd()
		IF (UR1 == 0.0_REAL_PRECISION) UR1 = grnd()
		UR2 = grnd()

		NR = DSQRT(-2.0_REAL_PRECISION*DLOG(UR1)) * DCOS(2.0_REAL_PRECISION*pi*(UR2))

		NRV(I) = NR

	END DO

END SUBROUTINE NormRNDVector





!---Subroutine ReversePeriodic----------------------------------------			! Reverse
!Set logical variable for periodic reversal of contraction
!---------------------------------------------------------------------

FUNCTION ReversePeriodic(TS)

	USE PARAMETERS, ONLY : ReversalPeriod
	

	INTEGER, VALUE, INTENT(IN) :: TS
	INTEGER :: ReversePeriodic



	IF (MOD(TS/ReversalPeriod,2) == 0) THEN
		ReversePeriodic = 0
	ELSE
		ReversePeriodic = 1
	END IF

END FUNCTION ReversePeriodic





!---Subroutine ReverseOneTime-----------------------------------------			! Reverse
!Set logical variable for one-time reversal of contraction
!---------------------------------------------------------------------

FUNCTION ReverseOneTime(TS)

	USE PARAMETERS, ONLY : TSReverseON
	

	INTEGER, VALUE, INTENT(IN) :: TS
	INTEGER :: ReverseOneTime



	IF (TS >= TSReverseON) THEN
		ReverseOneTime = 1
	ELSE
		ReverseOneTime = 0
	END IF

END FUNCTION ReverseOneTime




!---Subroutine PeriodicReverse2Contract-------------------------------			! Reverse
!Set logical variable for periodic contraction reversal and then contraction
!---------------------------------------------------------------------

FUNCTION PeriodicReverse2Contract(TS)

	USE PARAMETERS, ONLY : ReversalPeriod, TSContractionON
	

	INTEGER, VALUE, INTENT(IN) :: TS
	INTEGER :: PeriodicReverse2Contract



	IF (MOD(TS/ReversalPeriod,2) == 0) THEN
		PeriodicReverse2Contract = 0
	ELSE
		PeriodicReverse2Contract = 1
	END IF

	IF (TS >= TSContractionON) PeriodicReverse2Contract = 0

END FUNCTION PeriodicReverse2Contract





END MODULE FUNC
