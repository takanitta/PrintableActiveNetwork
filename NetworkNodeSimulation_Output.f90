MODULE OUTPUT

USE PARAMETERS

IMPLICIT NONE

CONTAINS





SUBROUTINE OutPutNodeData(OutputfileCounter, X, Y, NodeRadius, NodeType, FxI, FyI)

	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode


	INTEGER, INTENT(IN) :: OutputfileCounter
	INTEGER, DIMENSION(NumNode), INTENT(IN) :: NodeType
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode), INTENT(IN) :: X, Y, NodeRadius, FxI, FyI
	CHARACTER(LEN=50) :: FileName
	INTEGER :: Unit = 12, I


!---------------------------------------------------------------------

	WRITE(FileName,'(A,I8.8,A)') 'NodeData', OutputfileCounter, '.txt'
	OPEN(Unit,FILE=FileName)

	DO I=1, NumNode
		WRITE(Unit,'(2I8,3F20.8,I8,2F20.8)') OutputfileCounter, I, X(I), Y(I), NodeRadius(I), NodeType(I), FxI(I), FyI(I)
	END DO

	CLOSE(Unit)

END SUBROUTINE OutPutNodeData





SUBROUTINE OutPutBondData2(OutputfileCounter, BondState, RIJ, BondNatLength, FIJTension, NumMotorOnBond,vContract)

	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode


	INTEGER, INTENT(IN) :: OutputfileCounter
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode,NumNode), INTENT(IN) :: RIJ, BondNatLength, FIJTension, NumMotorOnBond,vContract
	INTEGER, DIMENSION(NumNode,NumNode), INTENT(IN) :: BondState
	CHARACTER(LEN=50) :: FileName
	INTEGER :: Unit = 12, I, J, NumBoundNode, NumBond, BoundNodeIndexCounterJ
	INTEGER, DIMENSION(NumNode) :: BoundNodeIndexI, BondNumI
	INTEGER, DIMENSION(NumNode,NumNode) :: BoundNodeIndex, TriUBondState, CompressedTriUBondState



!---------------------------------------------------------------------

	WRITE(FileName,'(A,I8.8,A)') 'BondData', OutputfileCounter, '.txt'
	OPEN(Unit,FILE=FileName)

	NumBoundNode = 0
	DO I=1, NumNode
		IF (SUM(BondState(I,:)) > 0) THEN
			NumBoundNode = NumBoundNode + 1
			BoundNodeIndexI(NumBoundNode) = I
		END IF
	END DO

	DO I=1,NumNode
		DO J=I,NumNode
			TriUBondState(I,J) = BondState(I,J)
		END DO
	END DO

	DO I=2,NumNode
		DO J=1,I-1
			TriUBondState(I,J) = 0
		END DO
	END DO

	NumBond = SUM(TriUBondState)


	CompressedTriUBondState = 0
	DO I=1, NumBoundNode
		BoundNodeIndexCounterJ = 0
		DO J=1, NumBoundNode
			IF (TriUBondState(BoundNodeIndexI(I),BoundNodeIndexI(J))==1) THEN
				BoundNodeIndexCounterJ = BoundNodeIndexCounterJ + 1
				BoundNodeIndex(I,BoundNodeIndexCounterJ) = BoundNodeIndexI(J)
				CompressedTriUBondState(I,J) = 1
			END IF
		END DO
		BondNumI(I) = BoundNodeIndexCounterJ
	END DO


	DO I=1,NumBoundNode
		DO J=1,BondNumI(I)
			WRITE(Unit,'(2I8,5F16.8)') BoundNodeIndexI(I),BoundNodeIndex(I,J),RIJ(BoundNodeIndexI(I),BoundNodeIndex(I,J)), BondNatLength(BoundNodeIndexI(I),BoundNodeIndex(I,J)), &
				FIJTension(BoundNodeIndexI(I),BoundNodeIndex(I,J)), NumMotorOnBond(BoundNodeIndexI(I),BoundNodeIndex(I,J)),vContract(BoundNodeIndexI(I),BoundNodeIndex(I,J))
		END DO
	END DO

	CLOSE(Unit)

END SUBROUTINE OutPutBondData2





SUBROUTINE OutPutNodeVTU(OutputfileCounter, X, Y, NodeRadius, NodeType)

	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode


	INTEGER, INTENT(IN) :: OutputfileCounter
	INTEGER, DIMENSION(NumNode), INTENT(IN) :: NodeType
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode), INTENT(IN) :: X, Y, NodeRadius
	CHARACTER(LEN=50) :: FileName
	INTEGER :: Unit = 11, I



!---------------------------------------------------------------------

	WRITE(FileName,'(A,I8.8,A)') 'Node', OutputfileCounter, '.vtu'
	OPEN(Unit,FILE=FileName)

	WRITE(Unit,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
	WRITE(Unit,'(A)') "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>"
	WRITE(Unit,'(A)') "<UnstructuredGrid>"
	WRITE(Unit,'(A,I,A,I,A)') "<Piece NumberOfCells=' ", NumNode, "' NumberOfPoints= '", NumNode, "' >"
	WRITE(Unit,'(A)') "<Points>"
	WRITE(Unit,'(A)') "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>"


	DO I=1, NumNode
		WRITE(Unit,'(3F16.8)') X(I), Y(I), 0.0
	END DO

	WRITE(Unit,'(A)') "</DataArray>"
	WRITE(Unit,'(A)') "</Points>"


	WRITE(Unit,'(A)') "<PointData>"

		WRITE(Unit,'(A)') "<DataArray NumberOfComponents='1' type='Float32' Name='NodeRadius' format='ascii'>"
		DO I=1, NumNode
			WRITE(Unit,'(F16.8)') NodeRadius(I)
		END DO
		WRITE(Unit,'(A)') "</DataArray>"


		WRITE(Unit,'(A)') "<DataArray NumberOfComponents='1' type='UInt8' Name='NodeType' format='ascii'>"
		DO I=1, NumNode
			WRITE(Unit,'(I)') NodeType(I)
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

	WRITE(Unit,'(A)') "</PointData>"


	WRITE(Unit,'(A)') "<Cells>"
		WRITE(Unit,'(A)') "<DataArray type='Int32' Name='connectivity' format='ascii'>"
		DO I=1, NumNode
			WRITE(Unit,'(I)') I-1
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

		WRITE(Unit,'(A)') "<DataArray type='Int32' Name='offsets' format='ascii'>"
		DO I=1, NumNode
			WRITE(Unit,'(I)') I
		END DO
		WRITE(Unit,'(A)') "</DataArray>"


		WRITE(Unit,'(A)') "<DataArray type='UInt8' Name='types' format='ascii'>"
		DO I=1, NumNode
			WRITE(Unit,'(I)') 1
		END DO
		WRITE(Unit,'(A)') "</DataArray>"
	WRITE(Unit,'(A)') "</Cells>"


	WRITE(Unit,'(A)') "</Piece>"
	WRITE(Unit,'(A)') "</UnstructuredGrid>"
	WRITE(Unit,'(A)') "</VTKFile>"


	CLOSE(Unit)

END SUBROUTINE OutPutNodeVTU





SUBROUTINE OutPutBondVTU3(OutputfileCounter, X, Y, BondState, RIJ, BondNatLength, FIJTension, NumMotorOnBond,vContract)

	USE PARAMETERS, ONLY : REAL_PRECISION, NumNode


	INTEGER, INTENT(IN) :: OutputfileCounter
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode), INTENT(IN) :: X, Y
	REAL(KIND = REAL_PRECISION), DIMENSION(NumNode,NumNode), INTENT(IN) :: RIJ, BondNatLength, FIJTension, NumMotorOnBond,vContract
	INTEGER, DIMENSION(NumNode,NumNode), INTENT(IN) :: BondState
	CHARACTER(LEN=50) :: FileName
	INTEGER :: Unit = 11, I, J, NumBoundNode, NumBond, BoundNodeIndexCounterJ
	INTEGER, DIMENSION(NumNode) :: BoundNodeIndexI, BondNumI
	INTEGER, DIMENSION(NumNode,NumNode) :: BoundNodeIndex, TriUBondState, CompressedTriUBondState



	NumBoundNode = 0
	DO I=1, NumNode
		IF (SUM(BondState(I,:)) > 0) THEN
			NumBoundNode = NumBoundNode + 1
			BoundNodeIndexI(NumBoundNode) = I
		END IF
	END DO


	DO I=1,NumNode
		DO J=I,NumNode
			TriUBondState(I,J) = BondState(I,J)
		END DO
	END DO
	DO I=2,NumNode
		DO J=1,I-1
			TriUBondState(I,J) = 0
		END DO
	END DO

	NumBond = SUM(TriUBondState)


	CompressedTriUBondState = 0
	DO I=1, NumBoundNode
		BoundNodeIndexCounterJ = 0
		DO J=1, NumBoundNode
			IF (TriUBondState(BoundNodeIndexI(I),BoundNodeIndexI(J))==1) THEN
				BoundNodeIndexCounterJ = BoundNodeIndexCounterJ + 1
				BoundNodeIndex(I,BoundNodeIndexCounterJ) = BoundNodeIndexI(J)
				CompressedTriUBondState(I,J) = 1
			END IF
		END DO
		BondNumI(I) = BoundNodeIndexCounterJ
	END DO



!---------------------------------------------------------------------

	WRITE(FileName,'(A,I8.8,A)') 'Bond', OutputfileCounter, '.vtu'
	OPEN(Unit,FILE=FileName)


	WRITE(Unit,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
	WRITE(Unit,'(A)') "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>"
	WRITE(Unit,'(A)') "<UnstructuredGrid>"
	WRITE(Unit,'(A,I,A,I,A)') "<Piece NumberOfCells=' ", NumBond, "' NumberOfPoints= '", NumBoundNode, "' >"
	WRITE(Unit,'(A)') "<Points>"
	WRITE(Unit,'(A)') "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>"
	DO I=1,NumBoundNode
		WRITE(Unit,'(3F16.8)') X(BoundNodeIndexI(I)), Y(BoundNodeIndexI(I)), 0.0
	END DO
	WRITE(Unit,'(A)') "</DataArray>"
	WRITE(Unit,'(A)') "</Points>"


	WRITE(Unit,'(A)') "<Cells>"
		WRITE(Unit,'(A)') "<DataArray type='Int32' Name='connectivity' format='ascii'>"
		DO I=1,NumBoundNode
			DO J=1,NumBoundNode
				IF (CompressedTriUBondState(I,J) == 1) WRITE(Unit,'(2I)') I-1, J-1
			END DO
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

		WRITE(Unit,'(A)') "<DataArray type='Int32' Name='offsets' format='ascii'>"
		DO I=1, NumBond
			WRITE(Unit,'(I)') I*2
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

		WRITE(Unit,'(A)') "<DataArray type='UInt8' Name='types' format='ascii'>"
		DO I=1, NumBond
			WRITE(Unit,'(I)') 3
		END DO
		WRITE(Unit,'(A)') "</DataArray>"
	WRITE(Unit,'(A)') "</Cells>"


	WRITE(Unit,'(A)') "<CellData>"
		WRITE(Unit,'(A)') "<DataArray NumberOfComponents='1' type='Float32' Name='NatLength' format='ascii'>"
		DO I=1,NumBoundNode
			DO J=1,BondNumI(I)
				WRITE(Unit,'(F16.8)') BondNatLength(BoundNodeIndexI(I),BoundNodeIndex(I,J))
			END DO
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

		WRITE(Unit,'(A)') "<DataArray NumberOfComponents='1' type='Float32' Name='RIJ' format='ascii'>"
		DO I=1, NumBoundNode
			DO J=1,BondNumI(I)
				WRITE(Unit,'(F16.8)') RIJ(BoundNodeIndexI(I),BoundNodeIndex(I,J))
			END DO
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

		WRITE(Unit,'(A)') "<DataArray NumberOfComponents='1' type='Float32' Name='TensionMotor' format='ascii'>"
		DO I=1, NumBoundNode
			DO J=1,BondNumI(I)
				WRITE(Unit,'(F16.8)') FIJTension(BoundNodeIndexI(I),BoundNodeIndex(I,J))
			END DO
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

		WRITE(Unit,'(A)') "<DataArray NumberOfComponents='1' type='Float32' Name='NumMotorOnBond' format='ascii'>"
		DO I=1, NumBoundNode
			DO J=1,BondNumI(I)
				WRITE(Unit,'(F16.8)') NumMotorOnBond(BoundNodeIndexI(I),BoundNodeIndex(I,J))
			END DO
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

		WRITE(Unit,'(A)') "<DataArray NumberOfComponents='1' type='Float32' Name='TensionBond' format='ascii'>"
		DO I=1, NumBoundNode
			DO J=1,BondNumI(I)
				WRITE(Unit,'(F16.8)') NumMotorOnBond(BoundNodeIndexI(I),BoundNodeIndex(I,J))*FIJTension(BoundNodeIndexI(I),BoundNodeIndex(I,J))
			END DO
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

		WRITE(Unit,'(A)') "<DataArray NumberOfComponents='1' type='Float32' Name='vContract' format='ascii'>"
		DO I=1, NumBoundNode
			DO J=1,BondNumI(I)
				WRITE(Unit,'(F16.8)') vContract(BoundNodeIndexI(I),BoundNodeIndex(I,J))
			END DO
		END DO
		WRITE(Unit,'(A)') "</DataArray>"

	WRITE(Unit,'(A)') "</CellData>"



	WRITE(Unit,'(A)') "</Piece>"
	WRITE(Unit,'(A)') "</UnstructuredGrid>"
	WRITE(Unit,'(A)') "</VTKFile>"

	CLOSE(Unit)

END SUBROUTINE OutPutBondVTU3





END MODULE OUTPUT
