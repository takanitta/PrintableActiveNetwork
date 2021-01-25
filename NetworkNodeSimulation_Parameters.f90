MODULE PARAMETERS

IMPLICIT NONE


! Unit: um, sec, pN



INTEGER, PARAMETER ::		INT_RANGE = SELECTED_INT_KIND(8), &
				REAL_PRECISION = SELECTED_REAL_KIND(8), &
				NumNode = 750, &									! Number of nodes
				TSEnd = 700E4, &									! Number of time steps
				ReversalPeriod = 50E4, &							! Duration of contraction or expansion in cyclic contraction/expansion
				TSReverseON = 50E4, &								! Time step of onset of expansion
				TSContractionON = 200E4, &							! Time step of onset of contraction after cyclic contraction/expansion
				TSOutIncr = 5E4, &									! Time step between output
				seed = 4397 ! 4397, 242, 6990, 52804, 2072544		! Random seed

REAL(KIND = REAL_PRECISION), PARAMETER :: 	pi = 3.14159265358979_REAL_PRECISION, &
				kBT = 0.0041_REAL_PRECISION, &						! Boltzmann constant multiplied with absolute temperature
				dt=2.0E-6_REAL_PRECISION, &							! Time between time step
				VicosityWater = 0.001_REAL_PRECISION, &				! Viscosity of water
				MotorDensity = 10.0_REAL_PRECISION, &				! Motor density (um^-1)
				k = 3.0_REAL_PRECISION, &							! Spring constant of kinesin
				vContractMax = 1.0_REAL_PRECISION, &				! Maximum speed of contraction
				k_ON = 5.0_REAL_PRECISION, &						! On-rate of kinesin binding to microtubule
				k_OFF0 = 0.3_REAL_PRECISION, &						! Off-rate of kinesin from microtubule without force
				CutoffNum = 0.5_REAL_PRECISION, &					! Cutoff value of separation between two binding nodes
				F_MotorStall = 5.0_REAL_PRECISION, &				! Kinesin stall force
				F_MotorDetach = 2.5_REAL_PRECISION, & 				! Kinesin unbinding force
				BondNatLengthMin = 2.5_REAL_PRECISION, &			! Minimum length of natural length of link between binding two nodes
				Radius = 15.0_REAL_PRECISION, &						! Initial radius of circular active network
				X_Size = 100.0_REAL_PRECISION, Y_Size = 30.0_REAL_PRECISION, &				! Size of rectangular active network
				kLever = 10.0_REAL_PRECISION, ViscDragCoeffLever = 0.5_REAL_PRECISION, &	! Spring constant and viscous drag coeeficient of flexible Lever
				GridSeparation = 2.0_REAL_PRECISION					! Initial separation between neighboring nodes of rectangular grid active network


END MODULE PARAMETERS
