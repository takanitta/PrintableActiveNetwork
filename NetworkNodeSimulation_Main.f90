PROGRAM NodeSimulation


USE PARAMETERS
USE mtmod
USE OUTPUT
USE FUNC

USE openacc_curand			! GPU



IMPLICIT NONE

INTEGER	::	I, J, TS=0, OutputfileCounter, NumMovableNode, TSSeed, Seq, Offset, &
		NumContractingBond, NumContractedBond, NumBond, &
		MovReverse = 0		! 0 = Contraction; 1 = Expansion

INTEGER, DIMENSION(NumNode) :: NodeType

INTEGER, DIMENSION(NumNode,NumNode) :: BondState, BondStateNext

REAL(KIND = REAL_PRECISION) :: FSummation, ElasticE, FAnchorTotal, XLever

REAL(KIND = REAL_PRECISION), DIMENSION(NumNode) :: 	X, Y, X0, Y0, XNext, YNext, NodeRadius, ViscDragCoeff, &
							FxI, FyI, NormRandVector4NodeX, NormRandVector4NodeY

REAL(KIND = REAL_PRECISION), DIMENSION(NumNode,NumNode) :: 	RIJ, BondNatLength, BondNatLengthNext, &
								FIJTension, &					! FIJTension: Force per motor protein
								vContract, NumMotorOnBond

TYPE(curandStateMRG32k3a) :: h			! GPU





OPEN (UNIT = 21, FILE = "Stat.txt")




CALL sgrnd(seed)


OutputfileCounter = 0


!---Dummy Values for CPU----------------------------------------------

NodeRadius = 0.0_REAL_PRECISION
NodeType = 0
BondState = 0
RIJ = 0.0_REAL_PRECISION
BondNatLength = 0.0_REAL_PRECISION
FIJTension = 0.0_REAL_PRECISION
vContract = 0.0_REAL_PRECISION
NumMotorOnBond = 0.0_REAL_PRECISION
FxI = 0.0_REAL_PRECISION
FyI = 0.0_REAL_PRECISION
FAnchorTotal = 0.0_REAL_PRECISION
XLever = 0.0_REAL_PRECISION



!---Set Initial Distribution of Nodes---------------------------------

!CALL InitialDistributionCircular(X,Y)
CALL InitialDistributionRectangular(X,Y)
!CALL InitialDistributionRectangularGrid(X,Y)



!---Set Node Radius---------------------------------------------------

NodeRadius = 2.5_REAL_PRECISION



!---Assign Node Types and Sort by Node Types--------------------------

CALL AssignNodeTypesSort(X,Y,NodeRadius,NodeType,NumMovableNode)



!---Data to GPU-------------------------------------------------------

!$acc data &
!$acc copyin(X, Y, NodeRadius, NodeType, BondState, RIJ, BondNatLength, FIJTension, vContract, NumMotorOnBond, FxI, FyI, FAnchorTotal, XLever, MovReverse, NumMovableNode, h) &
!$acc create(ViscDragCoeff, FSummation, X0, Y0, XNext, YNext, NormRandVector4NodeX, NormRandVector4NodeY, BondNatLengthNext, BondStateNext)



!---Initilize Random Number Generator on GPU--------------------------

!$acc kernels
Seq = 0																		! GPU
Offset = 0																	! GPU
CALL curand_init(seed, Seq, Offset, h)										! GPU
!$acc end kernels



!---Set Viscous Drag Coefficient--------------------------------------

!$acc data present(ViscDragCoeff,NodeRadius)
!$acc kernels
DO I=1,NumNode
	ViscDragCoeff(I)=6.0_REAL_PRECISION*pi*VicosityWater*NodeRadius(I)
END DO
!$acc end kernels
!$acc end data



!---Save Initial Positions--------------------------------------------

!$acc data present(X, Y, X0, Y0)
!$acc kernels
DO I=1,NumNode
	X0(I) = X(I)
	Y0(I) = Y(I)
END DO
!$acc end kernels
!$acc end data



OutputfileCounter = 0



!---Calculate Distance between Nodes----------------------------------

!$acc data present(X, Y, RIJ)
!$acc kernels
	DO J=1,NumNode																	
		DO I=1,NumNode																
			RIJ(I,J) = SQRT((X(I)-X(J))*(X(I)-X(J)) + (Y(I)-Y(J))*(Y(I)-Y(J)))		
		END DO																		
	END DO																			
!$acc end kernels
!$acc end data



!---Check Contact-----------------------------------------------------

!$acc data present(RIJ,NodeRadius,BondState)
!$acc kernels
	DO J=1,NumNode														
			DO I=1,NumNode												
			IF (RIJ(I,J) <= NodeRadius(I) + NodeRadius(J)) THEN			
				BondState(I,J) = 1										
			ELSE														
				BondState(I,J) = 0										
			END IF														
		END DO															
	END DO																

	DO I=1,NumNode														
		BondState(I,I) = 0												
	END DO																
!$acc end kernels
!$acc end data



!---Set Initial Natural Length of Bond--------------------------------

!$acc data present(BondNatLength, RIJ)
!$acc kernels
	DO J=1,NumNode
		DO I=1,NumNode
			BondNatLength(I,J) = RIJ(I,J)
		END DO
	END DO
!$acc end kernels
!$acc end data



!---Calculate Number of Motors on Bond--------------------------------

!$acc data present(BondState,NumMotorOnBond,NodeRadius,BondNatLength,FIJTension)
!$acc kernels
	DO J=1,NumNode
	!$acc loop
		DO I=1,NumNode
			IF (BondState(I,J)==1) NumMotorOnBond(I,J) = (NodeRadius(I) + NodeRadius(J) - BondNatLength(I,J))*MotorDensity * &
				k_ON/(k_ON + k_OFF0*EXP(ABS(FIJTension(I,J))/F_MotorDetach))
		END DO
	END DO
!$acc end kernels
!$acc end data



!---Calculate Tension-------------------------------------------------

!$acc data present(BondState,RIJ,BondNatLength,FIJTension)
!$acc kernels
	DO J=1,NumNode														
		DO I=1,NumNode													
			IF (BondState(I,J)==0) THEN									
				FIJTension(I,J) = 0.0_REAL_PRECISION					
			ELSE														
				FIJTension(I,J) = k*(RIJ(I,J) - BondNatLength(I,J))		
			END IF														
		END DO															
	END DO																	
!$acc end kernels
!$acc end data



!---Calculate Contraction Speed---------------------------------------

!$acc data present(BondState,vContract,MovReverse,FIJTension)
!$acc kernels
	DO J=1,NumNode
		DO I=1,NumNode
			IF (BondState(I,J)==1) vContract(I,J) = CalvContract(MovReverse,FIJTension(I,J))
		END DO
	END DO
!$acc end kernels
!$acc end data



!---Calculate Summation of Forces for Each Nodes----------------------

!$acc data present(FSummation,BondState,NumMotorOnBond,FIJTension,X,Y,RIJ)
!$acc parallel
	!$acc loop private(FSummation,I)
	DO J=1,NumNode
		FSummation = 0.0_REAL_PRECISION
		!$acc loop reduction(+:FSummation)
		DO I=1,NumNode
			IF (BondState(I,J)==1) FSummation = FSummation + NumMotorOnBond(I,J)*FIJTension(I,J)*(X(I)-X(J))/RIJ(I,J)
		END DO
		FxI(J) = FSummation
	END DO

	!$acc loop private(FSummation,I)
	DO J=1,NumNode
		FSummation = 0.0_REAL_PRECISION
		!$acc loop reduction(+:FSummation)
		DO I=1,NumNode
			IF (BondState(I,J)==1) FSummation = FSummation + NumMotorOnBond(I,J)*FIJTension(I,J)*(Y(I)-Y(J))/RIJ(I,J)
		END DO
		FyI(J) = FSummation
	END DO
!$acc end parallel
!$acc end data



!-----Calculate Total Anchor Force-------------------------------------

!$acc data present(FAnchorTotal,X0,FxI)																		
!$acc kernels																								
FAnchorTotal = 0.0_REAL_PRECISION																			
!$acc loop reduction(+:FAnchorTotal)																		
	DO I=NumMovableNode+1,NumNode																			
		IF (X0(I) < 0.0_REAL_PRECISION) FAnchorTotal = FAnchorTotal + FxI(I)								
	END DO																									
!$acc end kernels																							
!$acc end data																								



!---Output Data-------------------------------------------------------

!$acc update self(X, Y, NodeRadius, NodeType, BondState, RIJ, BondNatLength, FIJTension, vContract, NumMotorOnBond, FxI, FyI, FAnchorTotal, XLever, MovReverse)

CALL OutPutNodeVTU(OutputfileCounter, X, Y, NodeRadius, NodeType)
CALL OutPutBondVTU3(OutputfileCounter, X, Y, BondState, RIJ, BondNatLength, FIJTension, NumMotorOnBond, vContract)

CALL OutPutNodeData(OutputfileCounter, X, Y, NodeRadius, NodeType, FxI, FyI)
CALL OutPutBondData2(OutputfileCounter, BondState, RIJ, BondNatLength, FIJTension, NumMotorOnBond, vContract)

CALL CalStats(BondState, RIJ, BondNatLength, NumMotorOnBond, ElasticE, NumContractingBond, NumContractedBond, NumBond)
WRITE(21,'(2I12,F24.8,3I12,2F16.8)') OutputfileCounter, MovReverse, ElasticE, NumContractingBond, NumContractedBond, NumBond, FAnchorTotal, XLever



!---Time Evolution----------------------------------------------------

DO TS = 1,TSEnd

!---Set Contraction Mode----------------------------------------------

!	MovReverse = 0									! No contraction reversal
!	MovReverse = ReversePeriodic(TS)				! Periodic contraction reversal
!	MovReverse = ReverseOneTime(TS)					! One-time contraction reversal
	MovReverse = PeriodicReverse2Contract(TS)		! Periodic contraction reversal and then contraction

	
!$acc update device(MovReverse)


!---Calculate X- and Y-Coordinates------------------------------------

!	CALL NormRNDVector(NumMovableNode, NormRandVector4NodeX)		! CPU
!	CALL NormRNDVector(NumMovableNode, NormRandVector4NodeY)		! CPU

	!$acc data present(h,NumMovableNode, NormRandVector4NodeX, NormRandVector4NodeY)		! GPU
	!$acc parallel																			! GPU
		!$acc loop seq																		! GPU
		DO I=1,NumMovableNode																! GPU
			NormRandVector4NodeX(I) = curand_normal(h)										! GPU
			NormRandVector4NodeY(I) = curand_normal(h)										! GPU
		END DO																				! GPU
	!$acc end parallel																		! GPU
	!$acc end data																			! GPU


	!$acc data present(XNext,YNext,X,Y,ViscDragCoeff,FxI,FyI,XLever,FAnchorTotal,NumMovableNode,NormRandVector4NodeX,NormRandVector4NodeY)
	!$acc kernels 
		DO I=1,NumMovableNode
			XNext(I) = X(I) + dt/ViscDragCoeff(I)*FxI(I) + SQRT(2.0_REAL_PRECISION*kBT*dt/ViscDragCoeff(I))*NormRandVector4NodeX(I)
			YNext(I) = Y(I) + dt/ViscDragCoeff(I)*FyI(I) + SQRT(2.0_REAL_PRECISION*kBT*dt/ViscDragCoeff(I))*NormRandVector4NodeY(I)
		END DO

!-----Fixed Anchor
		DO I=NumMovableNode+1,NumNode								! Fixed Anchor
			XNext(I) = X(I)											! Fixed Anchor
			YNext(I) = Y(I)											! Fixed Anchor
		END DO														! Fixed Anchor

!-----Flexible Anchor

!		XLever = XLever - kLever*dt*XLever/ViscDragCoeffLever + FAnchorTotal*dt/ViscDragCoeffLever														! Flexible Anchor

!																											! Flexible Anchor
!		DO I=NumMovableNode+1,NumNode																		! Flexible Anchor
!			IF (X0(I) < 0.0_REAL_PRECISION) THEN															! Flexible Anchor
!				XNext(I) = X0(I) + XLever																	! Flexible Anchor
!			ELSE																							! Flexible Anchor
!				XNext(I) = X0(I)																			! Flexible Anchor
!			END IF																							! Flexible Anchor
!		END DO																								! Flexible Anchor
!																											! Flexible Anchor
!		DO I=NumMovableNode+1,NumNode																		! Flexible Anchor
!			YNext(I) = Y0(I)																				! Flexible Anchor
!		END DO																								! Flexible Anchor

	!$acc end kernels
	!$acc end data



!---Calculate Natural Length of Bond----------------------------------

	!$acc data present(NumMovableNode,BondState,BondNatLength,FIJTension,vContract,MovReverse,BondNatLengthNext,XNext,YNext,RIJ)
	!$acc kernels
	DO J=1,NumNode
		DO I=1,NumNode
			IF (BondState(I,J)==0) THEN
				BondNatLengthNext(I,J) = BondNatLength(I,J)
			ELSE
				BondNatLengthNext(I,J) = BondNatLength(I,J) - dt*vContract(I,J)
			END IF
		END DO
	END DO


	DO J=NumMovableNode+1,NumNode
		DO I=NumMovableNode+1,NumNode
			BondNatLengthNext(I,J) = BondNatLength(I,J)
		END DO
	END DO


	DO J=1,NumNode
		DO I=1,NumNode
			IF ((BondNatLengthNext(I,J) < BondNatLength(I,J)) .AND. (BondNatLengthNext(I,J) <= BondNatLengthMin)) THEN
				BondNatLengthNext(I,J) = MIN(BondNatLength(I,J),BondNatLengthMin)
			END IF
		END DO
	END DO



!---Calculate Distance between Nodes----------------------------------

		DO J=1,NumNode																									
			DO I=1,NumNode																								
				RIJ(I,J) = SQRT((XNext(I)-XNext(J))*(XNext(I)-XNext(J)) + (YNext(I)-YNext(J))*(YNext(I)-YNext(J)))		
			END DO																										
		END DO																											
	!$acc end kernels
	!$acc end data



!---Find New Contact--------------------------------------------------

	!$acc data present(RIJ,NodeRadius,BondState,BondStateNext)
	!$acc kernels	
		DO J=1,NumNode		
			DO I=1,NumNode
				IF ((BondState(I,J) == 1) .OR. (RIJ(I,J) <= NodeRadius(I) + NodeRadius(J))) THEN
					BondStateNext(I,J) = 1
				ELSE
					BondStateNext(I,J) = 0
				END IF
			END DO
		END DO

		
		DO I=1,NumNode
			BondStateNext(I,I) = 0
		END DO
	!$acc end kernels
	!$acc end data



!---Set Natural Length of Bond of Newly Formed Bond-------------------

	!$acc data present(BondState,BondStateNext,BondNatLength,BondNatLengthNext,RIJ,NumMotorOnBond,NodeRadius)
	!$acc kernels
		DO J=1,NumNode
			DO I=1,NumNode
				IF ((BondState(I,J)==0) .AND. (BondStateNext(I,J)==1)) THEN
					BondNatLengthNext(I,J) = RIJ(I,J)
				END IF
			END DO
		END DO	


!---Check Detachment--------------------------------------------------

		DO J=1,NumNode
			DO I=1,NumNode
				IF ((BondState(I,J)==1) .AND. ((NumMotorOnBond(I,J) < CutoffNum) .OR. (NodeRadius(I) + NodeRadius(J) - BondNatLength(I,J) < 0.0_REAL_PRECISION))) THEN
					BondStateNext(I,J) = 0
				END IF
			END DO
		END DO
	!$acc end kernels
	!$acc end data



!---Update Variables--------------------------------------------------

	!$acc data present(X,Y,XNext,YNext,BondNatLength,BondNatLengthNext)
	!$acc kernels
		DO I=1,NumNode
			X(I) = XNext(I)
		END DO

		DO I=1,NumMovableNode
			Y(I) = YNext(I)
		END DO

		DO J=1,NumNode
			DO I=1,NumNode
				BondNatLength(I,J) = BondNatLengthNext(I,J)
			END DO
		END DO

		DO J=1,NumNode
			DO I=1,NumNode
				BondState(I,J) = BondStateNext(I,J)
			END DO
		END DO
	!$acc end kernels
	!$acc end data



!---Update Number of Motors on Bond-----------------------------------

	!$acc data present(BondState,NumMotorOnBond,NodeRadius,BondNatLength,FIJTension)
	!$acc kernels
		DO J=1,NumNode
			DO I=1,NumNode
				IF (BondState(I,J)==1) NumMotorOnBond(I,J) = (NodeRadius(I) + NodeRadius(J) - BondNatLength(I,J))*MotorDensity * &
					k_ON/(k_ON + k_OFF0*EXP(ABS(FIJTension(I,J))/F_MotorDetach))
			END DO
		END DO
	!$acc end kernels
	!$acc end data



!---Calculate Tension-------------------------------------------------

	!$acc data present(BondState,RIJ,BondNatLength,FIJTension)
	!$acc kernels
	DO J=1,NumNode														
		DO I=1,NumNode													
			IF (BondState(I,J)==0) THEN									
				FIJTension(I,J) = 0.0_REAL_PRECISION					
			ELSE														
				FIJTension(I,J) = k*(RIJ(I,J) - BondNatLength(I,J))		
			END IF														
		END DO															
	END DO																
	!$acc end kernels
	!$acc end data



!---Calculate Contraction Speed---------------------------------------

!$acc data present(BondState,vContract,MovReverse,FIJTension)
!$acc kernels
	DO J=1,NumNode
		DO I=1,NumNode
			IF (BondState(I,J)==1) vContract(I,J) = CalvContract(MovReverse,FIJTension(I,J))
		END DO
	END DO
!$acc end kernels
!$acc end data



!---Calculate Summation of Forces for Each Nodes----------------------

	!$acc data present(FSummation,BondState,NumMotorOnBond,FIJTension,X,Y,RIJ)
	!$acc parallel
		!$acc loop private(FSummation,I)
		DO J=1,NumNode
			FSummation = 0.0_REAL_PRECISION
			!$acc loop reduction(+:FSummation)
			DO I=1,NumNode
				IF (BondState(I,J)==1) FSummation = FSummation + NumMotorOnBond(I,J)*FIJTension(I,J)*(X(I)-X(J))/RIJ(I,J)
			END DO
			FxI(J) = FSummation
		END DO

		!$acc loop private(FSummation,I)
		DO J=1,NumNode
			FSummation = 0.0_REAL_PRECISION
			!$acc loop reduction(+:FSummation)
			DO I=1,NumNode
				IF (BondState(I,J)==1) FSummation = FSummation + NumMotorOnBond(I,J)*FIJTension(I,J)*(Y(I)-Y(J))/RIJ(I,J)
			END DO
			FyI(J) = FSummation
		END DO
	!$acc end parallel
	!$acc end data



!-----Calculate Total Anchor Force-------------------------------------

!$acc data present(FAnchorTotal,X0,FxI)																		
!$acc kernels																								
FAnchorTotal = 0.0_REAL_PRECISION																			
!$acc loop reduction(+:FAnchorTotal)																		
	DO I=NumMovableNode+1,NumNode																			
		IF (X0(I) < 0.0_REAL_PRECISION) FAnchorTotal = FAnchorTotal + FxI(I)								
	END DO																									
!$acc end kernels																							
!$acc end data																								



!---Output Data-------------------------------------------------------

!$acc update self(X, Y, NodeRadius, NodeType, BondState, RIJ, BondNatLength, FIJTension, vContract, NumMotorOnBond, FxI, FyI, FAnchorTotal, XLever, MovReverse) if (mod(TS,TSOutIncr)==0)

	IF (mod(TS,TSOutIncr)==0) THEN

		OutputfileCounter = OutputfileCounter + 1

		CALL OutPutNodeVTU(OutputfileCounter, X, Y, NodeRadius, NodeType)
		CALL OutPutBondVTU3(OutputfileCounter, X, Y, BondState, RIJ, BondNatLength, FIJTension, NumMotorOnBond, vContract)

		CALL OutPutNodeData(OutputfileCounter, X, Y, NodeRadius, NodeType, FxI, FyI)
		CALL OutPutBondData2(OutputfileCounter, BondState, RIJ, BondNatLength, FIJTension, NumMotorOnBond, vContract)

		CALL CalStats(BondState, RIJ, BondNatLength, NumMotorOnBond, ElasticE, NumContractingBond, NumContractedBond, NumBond)
		WRITE(21,'(2I12,F24.8,3I12,2F16.8)') OutputfileCounter, MovReverse, ElasticE, NumContractingBond, NumContractedBond, NumBond, FAnchorTotal, XLever

	END IF


END DO

!$acc end data


CLOSE(21)


END PROGRAM NodeSimulation
