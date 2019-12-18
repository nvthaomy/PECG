!======MODULE CODE FOR PAA0======

module m

    real(8), parameter :: kB = 0.0019858775
    real(8), parameter :: FPE = 0.0030114705
    real(8), parameter :: ECharge = 1.0
    real(8), parameter :: eV = 23.052347773
    integer, parameter :: Dim = 3
    integer, dimension(0:72-1), parameter :: BondOrdData = (/ 1, 2, 3, 4, 2, 1, 2, 3, 4, 3, 2, 1, 2, 3, 4, 4, 3, 2, 1, 2, 3, 4, &
      & 4, 3, 2, 1, 2, 3, 4, 4, 3, 2, 1, 2, 3, 4, 4, 3, 2, 1, 2, 3, 4, 4, 3, 2, 1, 2, 3, 4, 4, 3, 2, 1, 2, 3, 4, 4, 3, 2, 1, 2, &
      & 3, 4, 3, 2, 1, 2, 4, 3, 2, 1 /)
    integer, dimension(0:13-1), parameter :: BondOrdStart = (/ 0, 4, 9, 15, 22, 29, 36, 43, 50, 57, 63, 68, 72 /)
    integer, dimension(0:12-1), parameter :: BondOrdShift = (/ 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8 /)
    integer, parameter :: BondOrdLimit = 4
    integer, parameter :: NAID = 2
    integer, parameter :: NSID = 12
    integer, parameter :: NMID = 1
    integer, parameter :: NAtom = 180
    integer, parameter :: NMol = 15
    integer, dimension(0:1-1), parameter :: NDOFMID = (/ 36 /)
    integer, dimension(0:1-1), parameter :: AtomsPerMol = (/ 12 /)
    integer, parameter :: MaxAtomsPerMol = 12
    integer, parameter :: MaxRigidAtoms = 0
    logical, dimension(0:15-1), parameter :: MolIsRigid = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    real(8), dimension(0:1-1), parameter :: RBondLengthSq = (/ -1.0 /)
    integer, dimension(0:1-1, 0:2-1), parameter :: RBondInd = reshape( (/ -1,-1 /) , (/1, 2/) )
    integer, dimension(0:2-1), parameter :: RBondRange = (/ 0, 0 /)
    integer, parameter :: NTerm = 14
    integer, parameter :: NDParam = 47
    integer, parameter :: NDDParam = 273
    integer, parameter :: P0_ArgHistNBin = 10000
    integer, parameter :: P1_ArgHistNBin = 10000
    integer, parameter :: P2_ArgHistNBin = 10000
    integer, parameter :: P3_ArgHistNBin = 10000
    integer, parameter :: P4_ArgHistNBin = 10000
    integer, parameter :: P5_ArgHistNBin = 10000
    integer, parameter :: P6_ArgHistNBin = 10000
    integer, parameter :: P7_ArgHistNBin = 10000
    integer, parameter :: P8_ArgHistNBin = 10000
    integer, parameter :: P9_ArgHistNBin = 10000
    integer, parameter :: P10_ArgHistNBin = 10000
    integer, parameter :: P11_ArgHistNBin = 10000
    integer, parameter :: P12_ArgHistNBin = 10000
    integer, parameter :: P13_ArgHistNBin = 10000
    integer, parameter :: M0_NVal = 1
    integer, parameter :: M1_NVal = 1
    integer, parameter :: M2_NVal = 1
    integer, parameter :: M3_NVal = 1
    integer, parameter :: M4_NVal = 1
    integer, parameter :: M5_NVal = 1
    integer, parameter :: M6_NVal = 1
    integer, parameter :: M7_NVal = 1
    integer, parameter :: M8_NVal = 1
    integer, parameter :: M9_NVal = 1
    integer, parameter :: M10_NVal = 47
    integer, parameter :: M11_NVal = 273
    integer, parameter :: M12_NVal = 47
    integer, parameter :: M13_NVal = 273
    integer, parameter :: M14_NVal = 2209
    integer, parameter :: M15_NVal = 1
    integer, parameter :: M16_NVal = 1
    integer, parameter :: M17_NVal = 1
    integer, parameter :: NMeasure = 18
    integer, parameter :: VV_NH_N = 2
    integer, parameter :: NMCMoves = 4
    real(8), dimension(0:3-1) :: BoxL
    integer, dimension(0:2-1) :: AIDCount
    integer, dimension(0:16-1) :: MolRange
    integer, dimension(0:15-1) :: MolID
    integer :: NDOF
    real(8), dimension(0:1-1, 0:12-1, 0:3-1) :: COMPos
    real(8), dimension(0:180-1, 0:3-1) :: Pos
    real(8), dimension(0:180-1, 0:3-1) :: Vel
    real(8), dimension(0:180-1, 0:3-1) :: Force
    real(8), dimension(0:180-1) :: Mass
    real(8), dimension(0:180-1) :: iMass
    real(8), dimension(0:180-1) :: sqrtMass
    integer, dimension(0:180-1) :: MInd
    integer :: NActiveMol
    integer :: NInactiveMol
    integer, dimension(0:1-1) :: NActiveMID
    integer, dimension(0:1-1) :: NInactiveMID
    integer, dimension(0:180-1) :: AID
    integer, dimension(0:180-1) :: SID
    integer, dimension(0:180-1) :: MID
    integer, dimension(0:15-1) :: MolActive
    real(8) :: PEnergy
    real(8) :: KEnergy
    real(8) :: TEnergy
    real(8) :: Virial
    real(8) :: TempSet
    real(8) :: PresSet
    integer, dimension(0:1-1) :: MuSet
    integer :: TargetAtom
    integer :: TargetMol
    real(8), dimension(0:14-1) :: Terms
    real(8), dimension(0:14-1) :: Cut
    real(8), dimension(0:14-1) :: CutSq
    real(8) :: OldPEnergy
    real(8), dimension(0:14-1) :: OldTerms
    real(8) :: FluctE
    real(8) :: FluctE0
    real(8) :: FluctA0
    real(8) :: FluctBeta
    real(8) :: FluctTerm
    real(8), dimension(0:47-1) :: Param
    real(8), dimension(0:47-1) :: DUParam
    real(8), dimension(0:47-1) :: DWParam
    real(8), dimension(0:273-1) :: DDUParam
    real(8), dimension(0:273-1) :: DDWParam
    real(8) :: P0_ArgWeightSumHist
    real(8) :: P0_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P0_ArgMin
    real(8), dimension(0:1-1) :: P0_ArgMax
    real(8), dimension(0:1-1) :: P0_ArgCount
    real(8), dimension(0:1-1) :: P0_ArgSum
    real(8), dimension(0:1-1) :: P0_ArgSumSq
    real(8), dimension(0:1-1) :: P0_ArgReportMin
    real(8), dimension(0:1-1) :: P0_ArgReportMax
    real(8), dimension(0:1-1) :: P0_ArgReportDelta
    real(8), dimension(0:1-1) :: P0_ArgHistMin
    real(8), dimension(0:1-1) :: P0_ArgHistMax
    real(8), dimension(0:1-1) :: P0_ArgHistBinw
    real(8), dimension(0:1-1) :: P0_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P0_ArgHist
    real(8), dimension(0:13-1) :: P1_UShift
    real(8) :: P1_ArgWeightSumHist
    real(8) :: P1_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P1_ArgMin
    real(8), dimension(0:1-1) :: P1_ArgMax
    real(8), dimension(0:1-1) :: P1_ArgCount
    real(8), dimension(0:1-1) :: P1_ArgSum
    real(8), dimension(0:1-1) :: P1_ArgSumSq
    real(8), dimension(0:1-1) :: P1_ArgReportMin
    real(8), dimension(0:1-1) :: P1_ArgReportMax
    real(8), dimension(0:1-1) :: P1_ArgReportDelta
    real(8), dimension(0:1-1) :: P1_ArgHistMin
    real(8), dimension(0:1-1) :: P1_ArgHistMax
    real(8), dimension(0:1-1) :: P1_ArgHistBinw
    real(8), dimension(0:1-1) :: P1_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P1_ArgHist
    real(8), dimension(0:13-1) :: P2_UShift
    real(8) :: P2_ArgWeightSumHist
    real(8) :: P2_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P2_ArgMin
    real(8), dimension(0:1-1) :: P2_ArgMax
    real(8), dimension(0:1-1) :: P2_ArgCount
    real(8), dimension(0:1-1) :: P2_ArgSum
    real(8), dimension(0:1-1) :: P2_ArgSumSq
    real(8), dimension(0:1-1) :: P2_ArgReportMin
    real(8), dimension(0:1-1) :: P2_ArgReportMax
    real(8), dimension(0:1-1) :: P2_ArgReportDelta
    real(8), dimension(0:1-1) :: P2_ArgHistMin
    real(8), dimension(0:1-1) :: P2_ArgHistMax
    real(8), dimension(0:1-1) :: P2_ArgHistBinw
    real(8), dimension(0:1-1) :: P2_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P2_ArgHist
    real(8), dimension(0:13-1) :: P3_UShift
    real(8) :: P3_ArgWeightSumHist
    real(8) :: P3_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P3_ArgMin
    real(8), dimension(0:1-1) :: P3_ArgMax
    real(8), dimension(0:1-1) :: P3_ArgCount
    real(8), dimension(0:1-1) :: P3_ArgSum
    real(8), dimension(0:1-1) :: P3_ArgSumSq
    real(8), dimension(0:1-1) :: P3_ArgReportMin
    real(8), dimension(0:1-1) :: P3_ArgReportMax
    real(8), dimension(0:1-1) :: P3_ArgReportDelta
    real(8), dimension(0:1-1) :: P3_ArgHistMin
    real(8), dimension(0:1-1) :: P3_ArgHistMax
    real(8), dimension(0:1-1) :: P3_ArgHistBinw
    real(8), dimension(0:1-1) :: P3_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P3_ArgHist
    real(8), dimension(0:13-1) :: P4_UShift
    real(8) :: P4_ArgWeightSumHist
    real(8) :: P4_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P4_ArgMin
    real(8), dimension(0:1-1) :: P4_ArgMax
    real(8), dimension(0:1-1) :: P4_ArgCount
    real(8), dimension(0:1-1) :: P4_ArgSum
    real(8), dimension(0:1-1) :: P4_ArgSumSq
    real(8), dimension(0:1-1) :: P4_ArgReportMin
    real(8), dimension(0:1-1) :: P4_ArgReportMax
    real(8), dimension(0:1-1) :: P4_ArgReportDelta
    real(8), dimension(0:1-1) :: P4_ArgHistMin
    real(8), dimension(0:1-1) :: P4_ArgHistMax
    real(8), dimension(0:1-1) :: P4_ArgHistBinw
    real(8), dimension(0:1-1) :: P4_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P4_ArgHist
    real(8), dimension(0:13-1) :: P5_UShift
    real(8) :: P5_ArgWeightSumHist
    real(8) :: P5_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P5_ArgMin
    real(8), dimension(0:1-1) :: P5_ArgMax
    real(8), dimension(0:1-1) :: P5_ArgCount
    real(8), dimension(0:1-1) :: P5_ArgSum
    real(8), dimension(0:1-1) :: P5_ArgSumSq
    real(8), dimension(0:1-1) :: P5_ArgReportMin
    real(8), dimension(0:1-1) :: P5_ArgReportMax
    real(8), dimension(0:1-1) :: P5_ArgReportDelta
    real(8), dimension(0:1-1) :: P5_ArgHistMin
    real(8), dimension(0:1-1) :: P5_ArgHistMax
    real(8), dimension(0:1-1) :: P5_ArgHistBinw
    real(8), dimension(0:1-1) :: P5_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P5_ArgHist
    real(8), dimension(0:13-1) :: P6_UShift
    real(8) :: P6_ArgWeightSumHist
    real(8) :: P6_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P6_ArgMin
    real(8), dimension(0:1-1) :: P6_ArgMax
    real(8), dimension(0:1-1) :: P6_ArgCount
    real(8), dimension(0:1-1) :: P6_ArgSum
    real(8), dimension(0:1-1) :: P6_ArgSumSq
    real(8), dimension(0:1-1) :: P6_ArgReportMin
    real(8), dimension(0:1-1) :: P6_ArgReportMax
    real(8), dimension(0:1-1) :: P6_ArgReportDelta
    real(8), dimension(0:1-1) :: P6_ArgHistMin
    real(8), dimension(0:1-1) :: P6_ArgHistMax
    real(8), dimension(0:1-1) :: P6_ArgHistBinw
    real(8), dimension(0:1-1) :: P6_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P6_ArgHist
    integer :: P7_ExcludeBondOrd
    real(8), dimension(0:3-1) :: P7_EWBoxL
    real(8) :: P7_EWMinL
    integer :: P7_NKVec
    integer :: P7_KMax
    real(8) :: P7_AlphaL
    real(8) :: P7_Alpha
    real(8), dimension(0:297-1) :: P7_KSq
    real(8), dimension(0:297-1) :: P7_KUTerm
    real(8), dimension(0:3-1, 0:297-1) :: P7_KFTerm
    integer, dimension(0:297-1) :: P7_K_x
    integer, dimension(0:297-1) :: P7_K_y
    integer, dimension(0:297-1) :: P7_K_z
    real(8) :: P7_ShiftTerm
    real(8) :: P7_Term1
    real(8) :: P7_Term2
    complex*16, dimension(0:297-1) :: P7_Struc
    complex*16, dimension(0:297-1) :: P7_OldStruc
    real(8) :: P7_RPEnergy
    real(8) :: P7_RVirial
    real(8) :: P7_SPEnergy
    real(8) :: P7_OldRPEnergy
    real(8) :: P7_OldRVirial
    real(8) :: P7_OldSPEnergy
    complex*16, dimension(0:6-1, 0:180-1) :: P7_eik_x
    complex*16, dimension(0:11-1, 0:180-1) :: P7_eik_y
    complex*16, dimension(0:11-1, 0:180-1) :: P7_eik_z
    complex*16, dimension(0:6-1, 0:180-1) :: P7_oldeik_x
    complex*16, dimension(0:11-1, 0:180-1) :: P7_oldeik_y
    complex*16, dimension(0:11-1, 0:180-1) :: P7_oldeik_z
    real(8) :: P7_ArgWeightSumHist
    real(8) :: P7_ArgWeightSumStats
    real(8), dimension(0:3-1) :: P7_ArgMin
    real(8), dimension(0:3-1) :: P7_ArgMax
    real(8), dimension(0:3-1) :: P7_ArgCount
    real(8), dimension(0:3-1) :: P7_ArgSum
    real(8), dimension(0:3-1) :: P7_ArgSumSq
    real(8), dimension(0:3-1) :: P7_ArgReportMin
    real(8), dimension(0:3-1) :: P7_ArgReportMax
    real(8), dimension(0:3-1) :: P7_ArgReportDelta
    real(8), dimension(0:3-1) :: P7_ArgHistMin
    real(8), dimension(0:3-1) :: P7_ArgHistMax
    real(8), dimension(0:3-1) :: P7_ArgHistBinw
    real(8), dimension(0:3-1) :: P7_ArgHistiBinw
    real(8), dimension(0:3-1, 0:10000-1) :: P7_ArgHist
    real(8), dimension(0:3-1) :: P8_UShift
    real(8), dimension(0:3-1) :: P8_CoulShift
    real(8) :: P8_ArgWeightSumHist
    real(8) :: P8_ArgWeightSumStats
    real(8), dimension(0:3-1) :: P8_ArgMin
    real(8), dimension(0:3-1) :: P8_ArgMax
    real(8), dimension(0:3-1) :: P8_ArgCount
    real(8), dimension(0:3-1) :: P8_ArgSum
    real(8), dimension(0:3-1) :: P8_ArgSumSq
    real(8), dimension(0:3-1) :: P8_ArgReportMin
    real(8), dimension(0:3-1) :: P8_ArgReportMax
    real(8), dimension(0:3-1) :: P8_ArgReportDelta
    real(8), dimension(0:3-1) :: P8_ArgHistMin
    real(8), dimension(0:3-1) :: P8_ArgHistMax
    real(8), dimension(0:3-1) :: P8_ArgHistBinw
    real(8), dimension(0:3-1) :: P8_ArgHistiBinw
    real(8), dimension(0:3-1, 0:10000-1) :: P8_ArgHist
    real(8), dimension(0:3-1) :: P9_UShift
    real(8), dimension(0:3-1) :: P9_CoulShift
    real(8) :: P9_ArgWeightSumHist
    real(8) :: P9_ArgWeightSumStats
    real(8), dimension(0:3-1) :: P9_ArgMin
    real(8), dimension(0:3-1) :: P9_ArgMax
    real(8), dimension(0:3-1) :: P9_ArgCount
    real(8), dimension(0:3-1) :: P9_ArgSum
    real(8), dimension(0:3-1) :: P9_ArgSumSq
    real(8), dimension(0:3-1) :: P9_ArgReportMin
    real(8), dimension(0:3-1) :: P9_ArgReportMax
    real(8), dimension(0:3-1) :: P9_ArgReportDelta
    real(8), dimension(0:3-1) :: P9_ArgHistMin
    real(8), dimension(0:3-1) :: P9_ArgHistMax
    real(8), dimension(0:3-1) :: P9_ArgHistBinw
    real(8), dimension(0:3-1) :: P9_ArgHistiBinw
    real(8), dimension(0:3-1, 0:10000-1) :: P9_ArgHist
    real(8), dimension(0:3-1) :: P10_UShift
    real(8), dimension(0:3-1) :: P10_CoulShift
    real(8) :: P10_ArgWeightSumHist
    real(8) :: P10_ArgWeightSumStats
    real(8), dimension(0:3-1) :: P10_ArgMin
    real(8), dimension(0:3-1) :: P10_ArgMax
    real(8), dimension(0:3-1) :: P10_ArgCount
    real(8), dimension(0:3-1) :: P10_ArgSum
    real(8), dimension(0:3-1) :: P10_ArgSumSq
    real(8), dimension(0:3-1) :: P10_ArgReportMin
    real(8), dimension(0:3-1) :: P10_ArgReportMax
    real(8), dimension(0:3-1) :: P10_ArgReportDelta
    real(8), dimension(0:3-1) :: P10_ArgHistMin
    real(8), dimension(0:3-1) :: P10_ArgHistMax
    real(8), dimension(0:3-1) :: P10_ArgHistBinw
    real(8), dimension(0:3-1) :: P10_ArgHistiBinw
    real(8), dimension(0:3-1, 0:10000-1) :: P10_ArgHist
    real(8), dimension(0:3-1) :: P11_UShift
    real(8), dimension(0:3-1) :: P11_CoulShift
    real(8) :: P11_ArgWeightSumHist
    real(8) :: P11_ArgWeightSumStats
    real(8), dimension(0:3-1) :: P11_ArgMin
    real(8), dimension(0:3-1) :: P11_ArgMax
    real(8), dimension(0:3-1) :: P11_ArgCount
    real(8), dimension(0:3-1) :: P11_ArgSum
    real(8), dimension(0:3-1) :: P11_ArgSumSq
    real(8), dimension(0:3-1) :: P11_ArgReportMin
    real(8), dimension(0:3-1) :: P11_ArgReportMax
    real(8), dimension(0:3-1) :: P11_ArgReportDelta
    real(8), dimension(0:3-1) :: P11_ArgHistMin
    real(8), dimension(0:3-1) :: P11_ArgHistMax
    real(8), dimension(0:3-1) :: P11_ArgHistBinw
    real(8), dimension(0:3-1) :: P11_ArgHistiBinw
    real(8), dimension(0:3-1, 0:10000-1) :: P11_ArgHist
    real(8), dimension(0:3-1) :: P12_UShift
    real(8), dimension(0:3-1) :: P12_CoulShift
    real(8) :: P12_ArgWeightSumHist
    real(8) :: P12_ArgWeightSumStats
    real(8), dimension(0:3-1) :: P12_ArgMin
    real(8), dimension(0:3-1) :: P12_ArgMax
    real(8), dimension(0:3-1) :: P12_ArgCount
    real(8), dimension(0:3-1) :: P12_ArgSum
    real(8), dimension(0:3-1) :: P12_ArgSumSq
    real(8), dimension(0:3-1) :: P12_ArgReportMin
    real(8), dimension(0:3-1) :: P12_ArgReportMax
    real(8), dimension(0:3-1) :: P12_ArgReportDelta
    real(8), dimension(0:3-1) :: P12_ArgHistMin
    real(8), dimension(0:3-1) :: P12_ArgHistMax
    real(8), dimension(0:3-1) :: P12_ArgHistBinw
    real(8), dimension(0:3-1) :: P12_ArgHistiBinw
    real(8), dimension(0:3-1, 0:10000-1) :: P12_ArgHist
    real(8), dimension(0:3-1) :: P13_UShift
    real(8), dimension(0:3-1) :: P13_CoulShift
    real(8) :: P13_ArgWeightSumHist
    real(8) :: P13_ArgWeightSumStats
    real(8), dimension(0:3-1) :: P13_ArgMin
    real(8), dimension(0:3-1) :: P13_ArgMax
    real(8), dimension(0:3-1) :: P13_ArgCount
    real(8), dimension(0:3-1) :: P13_ArgSum
    real(8), dimension(0:3-1) :: P13_ArgSumSq
    real(8), dimension(0:3-1) :: P13_ArgReportMin
    real(8), dimension(0:3-1) :: P13_ArgReportMax
    real(8), dimension(0:3-1) :: P13_ArgReportDelta
    real(8), dimension(0:3-1) :: P13_ArgHistMin
    real(8), dimension(0:3-1) :: P13_ArgHistMax
    real(8), dimension(0:3-1) :: P13_ArgHistBinw
    real(8), dimension(0:3-1) :: P13_ArgHistiBinw
    real(8), dimension(0:3-1, 0:10000-1) :: P13_ArgHist
    integer :: M0_StepFreq
    integer :: M0_CycleFreq
    logical :: M0_Active
    real(8), dimension(0:1-1) :: M0_Val
    real(8), dimension(0:1-1) :: M0_ValSum
    real(8), dimension(0:1-1) :: M0_ValSumSq
    real(8) :: M0_Count
    integer :: M1_StepFreq
    integer :: M1_CycleFreq
    logical :: M1_Active
    real(8), dimension(0:1-1) :: M1_Val
    real(8), dimension(0:1-1) :: M1_ValSum
    real(8), dimension(0:1-1) :: M1_ValSumSq
    real(8) :: M1_Count
    integer :: M2_StepFreq
    integer :: M2_CycleFreq
    logical :: M2_Active
    real(8), dimension(0:1-1) :: M2_Val
    real(8), dimension(0:1-1) :: M2_ValSum
    real(8), dimension(0:1-1) :: M2_ValSumSq
    real(8) :: M2_Count
    integer :: M3_StepFreq
    integer :: M3_CycleFreq
    logical :: M3_Active
    real(8), dimension(0:1-1) :: M3_Val
    real(8), dimension(0:1-1) :: M3_ValSum
    real(8), dimension(0:1-1) :: M3_ValSumSq
    real(8) :: M3_Count
    integer :: M4_StepFreq
    integer :: M4_CycleFreq
    logical :: M4_Active
    real(8), dimension(0:1-1) :: M4_Val
    real(8), dimension(0:1-1) :: M4_ValSum
    real(8), dimension(0:1-1) :: M4_ValSumSq
    real(8) :: M4_Count
    integer :: M5_StepFreq
    integer :: M5_CycleFreq
    logical :: M5_Active
    real(8), dimension(0:1-1) :: M5_Val
    real(8), dimension(0:1-1) :: M5_ValSum
    real(8), dimension(0:1-1) :: M5_ValSumSq
    real(8) :: M5_Count
    integer :: M6_StepFreq
    integer :: M6_CycleFreq
    logical :: M6_Active
    real(8), dimension(0:1-1) :: M6_Val
    real(8), dimension(0:1-1) :: M6_ValSum
    real(8), dimension(0:1-1) :: M6_ValSumSq
    real(8) :: M6_Count
    integer :: M7_StepFreq
    integer :: M7_CycleFreq
    logical :: M7_Active
    real(8), dimension(0:1-1) :: M7_Val
    real(8), dimension(0:1-1) :: M7_ValSum
    real(8), dimension(0:1-1) :: M7_ValSumSq
    real(8) :: M7_Count
    integer :: M8_StepFreq
    integer :: M8_CycleFreq
    logical :: M8_Active
    real(8), dimension(0:1-1) :: M8_Val
    real(8), dimension(0:1-1) :: M8_ValSum
    real(8), dimension(0:1-1) :: M8_ValSumSq
    real(8) :: M8_Count
    integer :: M9_StepFreq
    integer :: M9_CycleFreq
    logical :: M9_Active
    real(8), dimension(0:1-1) :: M9_Val
    real(8), dimension(0:1-1) :: M9_ValSum
    real(8), dimension(0:1-1) :: M9_ValSumSq
    real(8) :: M9_Count
    integer :: M10_StepFreq
    integer :: M10_CycleFreq
    logical :: M10_Active
    real(8), dimension(0:47-1) :: M10_Val
    real(8), dimension(0:47-1) :: M10_ValSum
    real(8), dimension(0:47-1, 0:47-1) :: M10_ValSumSq
    real(8) :: M10_Count
    integer :: M11_StepFreq
    integer :: M11_CycleFreq
    logical :: M11_Active
    real(8), dimension(0:273-1) :: M11_Val
    real(8), dimension(0:273-1) :: M11_ValSum
    real(8), dimension(0:273-1) :: M11_ValSumSq
    real(8) :: M11_Count
    integer :: M12_StepFreq
    integer :: M12_CycleFreq
    logical :: M12_Active
    real(8), dimension(0:47-1) :: M12_Val
    real(8), dimension(0:47-1) :: M12_ValSum
    real(8), dimension(0:47-1) :: M12_ValSumSq
    real(8) :: M12_Count
    integer :: M13_StepFreq
    integer :: M13_CycleFreq
    logical :: M13_Active
    real(8), dimension(0:273-1) :: M13_Val
    real(8), dimension(0:273-1) :: M13_ValSum
    real(8), dimension(0:273-1) :: M13_ValSumSq
    real(8) :: M13_Count
    integer :: M14_StepFreq
    integer :: M14_CycleFreq
    logical :: M14_Active
    real(8), dimension(0:2209-1) :: M14_Val
    real(8), dimension(0:2209-1) :: M14_ValSum
    real(8), dimension(0:2209-1) :: M14_ValSumSq
    real(8) :: M14_Count
    integer :: M15_StepFreq
    integer :: M15_CycleFreq
    logical :: M15_Active
    real(8), dimension(0:1-1) :: M15_Val
    real(8), dimension(0:1-1) :: M15_ValSum
    real(8), dimension(0:1-1) :: M15_ValSumSq
    real(8) :: M15_Count
    integer :: M16_StepFreq
    integer :: M16_CycleFreq
    logical :: M16_Active
    real(8), dimension(0:1-1) :: M16_Val
    real(8), dimension(0:1-1) :: M16_ValSum
    real(8), dimension(0:1-1) :: M16_ValSumSq
    real(8) :: M16_Count
    integer :: M17_StepFreq
    integer :: M17_CycleFreq
    logical :: M17_Active
    real(8), dimension(0:1-1) :: M17_Val
    real(8), dimension(0:1-1) :: M17_ValSum
    real(8), dimension(0:1-1) :: M17_ValSumSq
    real(8) :: M17_Count
    real(8) :: VVQ_TimeStep
    integer :: VVQ_RattleMaxIterPerAtom
    real(8) :: VVQ_RattleTol1
    real(8) :: VVQ_RattleTol2
    real(8) :: VVQ_RattleTol3
    real(8) :: VV_TimeStep
    integer :: VV_RattleMaxIterPerAtom
    real(8) :: VV_RattleTol1
    real(8) :: VV_RattleTol2
    real(8) :: VV_RattleTol3
    integer :: VV_Thermostat
    integer :: VV_AndersenStep
    integer :: VV_AndersenStepFreq
    real(8) :: VV_AndersenCollisionFreq
    real(8), dimension(0:2-1) :: VV_Glogs
    real(8), dimension(0:2-1) :: VV_Vlogs
    real(8), dimension(0:2-1) :: VV_Xlogs
    real(8), dimension(0:2-1) :: VV_QMass
    real(8) :: VV_LangevinGamma
    integer :: VV_Barostat
    integer :: VV_BarostatStepFreq
    real(8) :: VV_BarostatDeltaV
    logical, dimension(0:3-1) :: VV_BarostatUseAxis
    logical :: VV_BarostatDoIsotropic
    real(8) :: VV_BarostatMaxVol
    real(8) :: VV_BarostatMinVol
    real(8) :: VV_BarostatNAtt
    real(8) :: VV_BarostatNAcc
    integer :: VV_RemoveCOMStep
    integer :: VV_RemoveCOMStepFreq
    real(8) :: VV_TEnergySum
    real(8) :: VV_TEnergySqSum
    real(8) :: MC0_NAtt
    real(8) :: MC0_NAcc
    real(8) :: MC0_NAtt2
    real(8) :: MC0_NAcc2
    real(8) :: MC0_Delta
    real(8) :: MC0_Delta2
    real(8) :: MC1_NAtt
    real(8) :: MC1_NAcc
    integer :: MC1_NInsertMols
    integer :: MC1_NDeleteMols
    integer, dimension(0:15-1) :: MC1_MolMoves
    real(8), dimension(0:1-1, 0:16-1) :: MC1_BoltzWeights
    integer, dimension(0:1-1) :: MC1_NAccMID
    integer, dimension(0:1-1) :: MC1_NAttMID
    real(8), dimension(0:16-1, 0:3-1) :: MC1_TM
    real(8) :: MC1_MFactor
    real(8) :: MC2_NAtt
    real(8) :: MC2_NAcc
    logical, dimension(0:3-1) :: MC2_UseAxis
    logical :: MC2_DoIsotropic
    real(8) :: MC2_Delta
    real(8) :: MC2_MaxVol
    real(8) :: MC2_MinVol
    real(8) :: MC3_NAtt
    real(8) :: MC3_NAcc
    integer :: MC3_MID1
    integer :: MC3_MID2
    real(8), dimension(0:16-1) :: MC3_BoltzWeights
    real(8), dimension(0:16-1, 0:3-1) :: MC3_TM
    real(8) :: MC3_MFactor
    real(8), dimension(0:4-1) :: MCMoveProbs
    real(8), dimension(0:180-1, 0:3-1) :: OldPos


contains

subroutine updateactive()
    implicit none
    integer :: i
    integer :: j
    integer :: Start
    integer :: Stop
    integer :: ThisMolID
    integer :: ThisAID

    NDOF = 0
    AIDCount = 0
    NActiveMol = 0
    NInactiveMol = 0
    NActiveMID = 0
    NInactiveMID = 0
    do i = 0, NMol-1
        Start = MolRange(i)
        Stop = MolRange(i+1) - 1
        if (MolActive(i) == 1) then
            NActiveMol = NActiveMol + 1
            ThisMolID = MolID(i)
            NActiveMID(ThisMolID) = NActiveMID(ThisMolID) + 1
            NDOF = NDOF + NDOFMID(ThisMolID)
            do j = Start, Stop
                ThisAID = AID(j)
                AIDCount(ThisAID) = AIDCount(ThisAID) + 1
            enddo
        elseif (MolActive(i) == -1) then
            NInactiveMol = NInactiveMol + 1
            ThisMolID = MolID(i)
            NInactiveMID(ThisMolID) = NInactiveMId(ThisMolID) + 1
        endif
    enddo

end subroutine


subroutine hidemol(MolInd)
    implicit none
    integer, intent(in) :: MolInd
    integer :: i
    integer :: ThisMolID
    integer :: ThisAID

    if (MolActive(MolInd) == 1) then
        ThisMolID = MID(MolInd)
        MolActive(MolInd) = -1
        NActiveMol = NActiveMol - 1
        NInactiveMol = NInactiveMol + 1
        NActiveMID(ThisMolID) = NActiveMID(ThisMolID) - 1
        NInactiveMID(ThisMolID) = NInactiveMID(ThisMolID) + 1
        NDOF = NDOF - NDOFMID(ThisMolID)
        do i = MolRange(MolInd), MolRange(MolInd + 1) - 1
            ThisAID = AID(i)
            AIDCount(ThisAID) = AIDCount(ThisAID) - 1
        enddo
    endif

end subroutine


subroutine showmol(MolInd)
    implicit none
    integer, intent(in) :: MolInd
    integer :: i
    integer :: ThisMolID
    integer :: ThisAID

    if (MolActive(MolInd) == -1) then
        ThisMolID = MID(MolInd)
        MolActive(MolInd) = 1
        NActiveMol = NActiveMol + 1
        NInactiveMol = NInactiveMol - 1
        NActiveMID(ThisMolID) = NActiveMID(ThisMolID) + 1
        NInactiveMID(ThisMolID) = NInactiveMID(ThisMolID) - 1
        NDOF = NDOF + NDOFMID(ThisMolID)
        do i = MolRange(MolInd), MolRange(MolInd+1) - 1
            ThisAID = AID(i)
            AIDCount(ThisAID) = AIDCount(ThisAID) + 1
        enddo
    endif

end subroutine


subroutine calcenergyforces(Mode, CalcForce, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    integer, intent(in) :: Mode
    logical, intent(in) :: CalcForce
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: k
    integer :: ForceDoMInd
    integer :: ForceDoAtom
    integer :: AIDi
    integer :: AIDj
    integer :: SIDi
    integer :: SIDj
    integer :: MIndi
    integer :: MIndj
    integer :: BondOrdIndShift
    integer :: BondOrdSIDjStart
    integer :: BondOrdSIDjStop
    integer :: BondOrdij
    integer :: istart
    integer :: istop
    integer :: jstart
    integer :: LoopMode
    logical :: SameMol
    logical :: Bonded
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8), parameter :: pi = 3.141592653589d0
    real(8), parameter :: sqrtpi = 1.772453850905d0
    real(8) :: Scale
    real(8), dimension(0:Dim-1) :: Forcei
    real(8) :: ThisU
    real(8) :: ThisW
    real(8) :: val1
    real(8) :: idist2
    real(8) :: idist6
    real(8) :: idist12
    real(8) :: val2
    real(8) :: Sig
    real(8) :: Eps
    real(8) :: val3
    real(8) :: val4
    real(8) :: val5
    real(8) :: val6
    real(8) :: val7
    real(8) :: Chargei
    real(8) :: Chargej
    real(8) :: idist
    real(8) :: erfcterm
    real(8) :: Temp1
    real(8) :: Temp2
    real(8) :: OldMinL
    real(8), dimension(0:Dim-1) :: OldBoxL
    real(8), dimension(0:Dim-1) :: Temp3
    real(8), dimension(0:Dim-1) :: Temp4
    real(8) :: EWScale
    integer :: KMaxSq
    integer :: KVecSq
    integer :: TotK
    complex*16 :: eikri
    complex*16 :: eikrj
    real(8) :: fac
    real(8) :: fac2
    real(8) :: iA
    real(8) :: tmp

    LoopMode = Mode
    if (Mode == 0) then
        !zero initial quantities
        PEnergy = 0.d0
        Terms = 0.d0
        if (CalcVirial) then
            Virial = 0.d0
        else
            Virial = 0.d0
        endif
        if (CalcForce) Force = 0.d0
        if (CalcDUParam .or. CalcFluct) then
            DUParam = 0.d0
            DDUParam = 0.d0
        endif
        if (CalcDWParam) then
            DWParam = 0.d0
            DDWParam = 0.d0
        endif
        if (CalcFluct) then
            FluctE0 = 0.
            FluctA0 = 0.
        endif
        Scale = 1.d0
    else
        if (CalcDUParam .or. CalcDWParam .or. CalcForce .or. CalcVirial .or. CalcFluct) then
            print *, "Can only use Calcvirial, CalcForce, CalcDUParam, CalcDWParam, CalcFluct for Mode=0."
            stop
        endif
        if (Mode > 0) then
            Scale = 1.d0
        elseif (Mode < 0) then
            Scale = -1.d0
        else
            print *, "Invalid value of Mode in calcenergyforces."
            stop
        endif
    endif

    !potential EW
    iBoxL = 1.d0 / BoxL
    if (any(BOXL /= P7_EWBoxL)) then
        OldMinL = P7_EWMinL
        OldBoxL = P7_EWBoxL
        P7_EWBoxL = BOXL
        P7_EWMinL = minval(P7_EWBoxL)
        P7_Alpha = P7_AlphaL/P7_EWMinL
        P7_Term1 = 2. * P7_Alpha / sqrtpi
        P7_Term2 = -P7_Alpha * P7_Alpha
        Temp1 = P7_EWMinL/OldMinL
        !update the truncation options, and the cutoff distance
        call erfc(P7_Alpha * Cut(7), erfcterm)
        P7_ShiftTerm = erfcterm / Cut(7)
        !check to see if we need to recalculate kvectors,
        !or can just scale them for an isotropic volume change
        Temp3 = P7_EWBoxL / OldBoxL
        if (abs(Temp3(0)-Temp3(1)) < 1.d-5 .and. abs(Temp3(1)-Temp3(2)) < 1.e-5 .and. all(OldBoxL > 0.)) then
            P7_KSq = P7_KSq / (Temp1 * Temp1)
            P7_KUTerm = P7_KUTerm * Temp1 * Temp1
            P7_KFTerm = P7_KFTerm / (Temp1 * Temp1)
        else
            TotK = 0
            P7_K_x = 0.
            P7_K_y = 0.
            P7_K_z = 0.
            P7_KSq = 0.
            P7_KUTerm = 0.
            P7_KFTerm = 0.
            KMaxSq = P7_KMax * P7_KMax
            do i = 0, P7_KMax
                do j = -P7_KMax, P7_KMax
                    do k = -P7_KMax, P7_KMax
                        KVecSq = i*i+j*j+k*k
                        if ( i==0 ) then
                            Temp1 = 1.
                        else
                            Temp1 = 2.
                        endif
                        if (KVecSq/=0 .and. KVecSq<=KMaxSq) then
                            P7_K_x(TotK) = i
                            P7_K_y(TotK) = j
                            P7_K_z(TotK) = k
                            !calculate k^2
                            P7_KSq(TotK) = 4.*Pi*Pi*sum( (real((/i, j, k/)) * iBoxL)**2 )
                            !calculate exp(-k^2/4a^2)/k^2
                            P7_KUTerm(TotK) = Temp1 * exp(-P7_KSq(TotK)/(4. * P7_Alpha * P7_Alpha))/P7_KSq(TotK)
                            !calculate (4pi/V)*kxyz*exp(-k^2/4a^2)/k^2 for forces
                            P7_KFTerm(:,TotK) = 8.*Pi*Pi * Param(34) * P7_KUTerm(TotK) * (product(iBoxL) * iBoxL)
                            P7_KFTerm(0,TotK) = P7_KFTerm(0,TotK) * real(P7_K_x(TotK))
                            P7_KFTerm(1,TotK) = P7_KFTerm(1,TotK) * real(P7_K_y(TotK))
                            P7_KFTerm(2,TotK) = P7_KFTerm(2,TotK) * real(P7_K_z(TotK))
                            TotK = TotK + 1
                        endif
                    enddo
                enddo
            enddo
            if (.not. TotK == P7_NKVec) then
                print *, "Total number of K vectors is unexpected in Ewald sum."
                stop
            endif
        endif
    endif
    if (Mode == 0) then
        !normal full atom loop
        istart = 0
        istop = NAtom - 1
        EWScale = 1.
    elseif (Mode == 1) then
        !single atom interactions for adding TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        EWScale = 1.
    elseif (Mode == -1) then
        !single atom interactions for deleting TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        EWScale = -1.
    elseif (Mode == 2) then
        !molecule interactions for adding TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        EWScale = 1.
    elseif (Mode == -2) then
        !molecule interactions for deleting TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        EWScale = -1.
    endif

    !calculate e^(ik*r) for kx = 0,1 and ky/kz = -1,0,1
    Temp4 = 2.*Pi * iBoxL
    do i = iStart, iStop
        !check if mol active
        if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
        Temp3 = -Pos(i,:) * Temp4
        P7_eik_y(P7_KMax-1,i) = dcmplx(cos(Temp3(1)), sin(Temp3(1)))
        P7_eik_z(P7_KMax-1,i) = dcmplx(cos(Temp3(2)), sin(Temp3(2)))
        P7_eik_x(0,i) = dcmplx( 1., 0. )
        P7_eik_y(P7_KMax,i) = dcmplx( 1., 0. )
        P7_eik_z(P7_KMax,i) = dcmplx( 1., 0. )
        Temp3 = Pos(i,:) * Temp4
        P7_eik_x(1,i) = dcmplx(cos(Temp3(0)), sin(Temp3(0)))
        P7_eik_y(P7_KMax+1,i) = dcmplx(cos(Temp3(1)), sin(Temp3(1)))
        P7_eik_z(P7_KMax+1,i) = dcmplx(cos(Temp3(2)), sin(Temp3(2)))
    enddo
    !calculate remaining e^(ik*r) by recursion
    do i = iStart, iStop
        if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
        do j = 2, P7_KMax
            !check if mol active
            P7_eik_x(j,i) = P7_eik_x(1,i) * P7_eik_x(j-1,i)
            P7_eik_y(P7_KMax+j,i) = P7_eik_y(P7_KMax+1,i) * P7_eik_y(P7_KMax+j-1,i)
            P7_eik_z(P7_KMax+j,i) = P7_eik_z(P7_KMax+1,i) * P7_eik_z(P7_KMax+j-1,i)
            P7_eik_y(P7_KMax-j,i) = dconjg(P7_eik_y(P7_KMax+j,i))
            P7_eik_z(P7_KMax-j,i) = dconjg(P7_eik_z(P7_KMax+j,i))
        enddo
    enddo
    !calculate the structure factor p(k) = sum( e^(ik*r) * q )
    if (Mode == 0) then
        P7_Struc = ( 0., 0. )
    endif
    do i = iStart, iStop
        !check if mol active
        if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
        Chargei = Param(AID(i))
        do j = 0, P7_NKVec - 1
            eikri = P7_eik_x(P7_K_x(j),i) * P7_eik_y(P7_K_y(j)+P7_KMax,i) * P7_eik_z(P7_K_z(j)+P7_KMax,i)
            P7_Struc(j) = P7_Struc(j) + eikri * Chargei * EWScale
        enddo
    enddo
    Temp2 = 0.5/(P7_Alpha*P7_Alpha)
    !val1 = reciprocal energy, val2 = recip virial, val3 = self energy
    val1 = 0.
    val2 = 0.
    do i = 0, P7_NKVec-1
        Temp1 = P7_KUTerm(i) * real(P7_Struc(i) * dconjg(P7_Struc(i)))
        val1 = val1 + Temp1
        val2 = val2 - Temp1 * ( 1. - Temp2 * P7_KSq(i) )
    enddo
    val1 = val1 * 2. * Pi * Param(34) * product(iBoxL)
    val2 = val2 * 2. * Pi * Param(34) * product(iBoxL)
    !potential self part
    val3 = 0.
    do i = iStart, iStop
        if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
        Chargei = Param(AID(i))
        val3 = val3 + Chargei**2
    enddo
    val3 = -P7_Alpha / sqrtpi * Param(34) * val3 * EWScale
    if (Mode == 0) then
        !a full update
        ThisU = val1 + val3
        ThisW = val2
        P7_RPEnergy = val1
        P7_RVirial = val2
        P7_SPEnergy = val3
    else
        !a partial update (ThisU, ThisW represent delta quantities)
        ThisU = (val1 - P7_RPEnergy) + val3
        ThisW = val2 - P7_RVirial
        P7_RPEnergy = val1
        P7_RVirial = val2
        P7_SPEnergy = P7_SPEnergy + val3
    endif
    PENERGY = PENERGY + ThisU
    Terms(7) = Terms(7) + ThisU
    if (CALCVIRIAL) VIRIAL = VIRIAL + ThisW
    if (CALCFORCE) then
        !forces reciprocal part
        do i = iStart, iStop
            if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
            Temp3 = 0.
            do j = 0, P7_NKVec-1
                eikri = P7_eik_x(P7_K_x(j), i) * P7_eik_y(P7_K_y(j)+P7_KMax, i) * P7_eik_z(P7_K_z(j)+P7_KMax, i)
                Temp3 = Temp3 - P7_KFTerm(:,j) * dimag(dconjg(eikri) * P7_Struc(j))
            enddo
            Temp3 = Temp3 * Param(AID(i))
            Force(i,:) = Force(i,:) + Temp3 * Scale
        enddo
    endif
    if (CALCDUPARAM .or. CALCDWPARAM) then
        val1 = 2. * Pi * product(iBoxL)
        val2 = Param(34) * val1
        val3 = -2. * P7_Alpha / sqrtpi
        Temp2 = 0.5/(P7_Alpha*P7_Alpha)
        DUParam(34) = DUParam(34) + ThisU / Param(34)
        if (CALCDWPARAM) then
            DWParam(34) = DWParam(34) + ThisW / Param(34)
        endif
        do i = iStart, iStop
            if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
            AIDI = AID(i)
            !reciprocal contribution
            do k = 0, P7_NKVec-1
                val5 = -(1. - Temp2 * P7_KSq(k))
                eikri = P7_eik_x(P7_K_x(k),i) * P7_eik_y(P7_K_y(k)+P7_KMax,i) * P7_eik_z(P7_K_z(k)+P7_KMax,i)
                eikri = dconjg(eikri)
                val4 = 2. * P7_KUTerm(k) * real(P7_Struc(k) * eikri)
                DUParam(AIDI) = DUParam(AIDI) + val2 * val4
                DDUParam(247 + AIDI + 0) = DDUParam(247 + AIDI + 0) + val1 * val4
                if (CALCDWPARAM) then
                    DWParam(AIDI) = DWParam(AIDI) + val2 * val4 * val5
                    DDWParam(247 + AIDI + 0) = DDWParam(247 + AIDI + 0) + val1 * val4 * val5
                endif
                do j = i, iStop
                    AIDJ = AID(j)
                    if (i == j) then
                        Temp1 = 2.
                    elseif (AIDI == AIDJ) then
                        Temp1 = 4.
                    else
                        Temp1 = 2.
                    endif
                    eikrj = P7_eik_x(P7_K_x(k),j) * P7_eik_y(P7_K_y(k)+P7_KMax,j) * P7_eik_z(P7_K_z(k)+P7_KMax,j)
                    val4 = Temp1 * val2 * P7_KUTerm(k) * real(eikrj * eikri)
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + val4
                    if (CALCDWPARAM) then
                        DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + val4 * val5
                    endif
                enddo
            enddo
            !self contribution
            DUParam(AIDI) = DUParam(AIDI) + Param(34) * val3 * Param(AIDI)
            DDUParam(AIDI + 2*AIDI) = DDUParam(AIDI + 2*AIDI) + Param(34) * val3
            DDUParam(247 + AIDI + 0) = DDUParam(247 + AIDI + 0) + val3 * Param(AIDI)
        enddo
    endif

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    ForceDoMInd = -1
    ForceDoAtom = -1
    if (LoopMode == 0) then
        !normal full atom loop
        istart = 0
        istop = NAtom - 1
    elseif (LoopMode == 1) then
        !single atom interactions for adding TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == -1) then
        !single atom interactions for deleting TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == 2) then
        !molecule interactions fo adding TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    elseif (LoopMode == -2) then
        !molecule interactions fo deleting TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    else
        print *, "Illegal value for LoopMode."
    endif

    !loop over i
    do i = istart, istop

        MIndi = MInd(i)

        if (MolActive(MIndi) < 0 .and. LoopMode == 0) cycle

        Posi = Pos(i,:)
        AIDi = AID(i)
        SIDi = SID(i)
        BondOrdSIDjStart = BondOrdShift(SIDi)
        BondOrdSIDjStop = BondOrdSIDjStart + BondOrdStart(SIDi+1) - BondOrdStart(SIDi)
        BondOrdIndShift = BondOrdStart(SIDi) - BondOrdShift(SIDi)

        !potential EW
        Chargei = Param(AIDI)

        if (AIDi==1) then
            !potential SmearCoulA_A
            Chargei = Param(AIDI)
        end if

        if (AIDi==1) then
            !potential SmearCoulA_A-
            Chargei = Param(AIDI)
        end if

        if (AIDi==1) then
            !potential SmearCoulA_Na+
            Chargei = Param(AIDI)
        end if

        if (AIDi==0) then
            !potential SmearCoulA-_A-
            Chargei = Param(AIDI)
        end if

        if (AIDi==0) then
            !potential SmearCoulA-_Na+
            Chargei = Param(AIDI)
        end if

        !potential SmearCoulNa+_Na+
        Chargei = Param(AIDI)

        if (LoopMode == 0) then
            jstart = i+1
        else
            jstart = 0
        endif

        !loop over j
        do j = jstart, NAtom - 1

            !check to see if same atom
            if (i==j) cycle

            MIndj = MInd(j)

            if (LoopMode == 2 .or. LoopMode == -2) then
                !!check to see if we need to skip because of double counting
                if (MIndj == MIndi .and. j < i) cycle
            endif

            if (MolActive(MIndj) < 0 .and. MIndj /= ForceDoMInd .and. j /= ForceDoAtom) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)
            SIDj = SID(j)
            SameMol = (MIndi == MIndj)
            BondOrdij = BondOrdLimit + 1
            if (SameMol) then
                if (SIDj >= BondOrdSIDjStart .and. SIDj < BondOrdSIDjStop) then
                    BondOrdij = BondOrdData(SIDj + BondOrdIndShift)
                endif
            endif
            Bonded = (SameMol .and. BondOrdij==2)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            if (dijsq < CutSq(0)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if ((((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) .and. Bonded) then
                    !potential BondA_A-
                    val1 = DIJ - Param(2)
                    THISU = Param(3) * val1*val1
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(0) = Terms(0) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = 2.d0 * Param(3) * val1*DIJ
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(3) = DUParam(3) + val1 * val1
                        DUParam(2) = DUParam(2) - 2.d0 * Param(3) * val1
                        DDUParam(4) = DDUParam(4) + 2.d0 * Param(3)
                        DDUParam(6) = DDUParam(6) - 2.d0 * val1
                    endif
                    if (CALCDWPARAM) then
                        DWParam(3) = DWParam(3) + 2.d0 * val1 * DIJ
                        DWParam(2) = DWParam(2) - 2.d0 * Param(3) * DIJ
                        DDWParam(6) = DDWParam(6) - 2.d0 * DIJ
                    endif
                end if
            end if

            if (dijsq < CutSq(1)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential LJGaussA_A
                    Sig = Param(5)
                    Eps = Param(4)
                    idist2 = Sig**2 / DIJSQ
                    idist6 = idist2 * idist2 * idist2
                    idist12 = idist6 * idist6
                    val1 = 4.d0 * (idist12 - idist6)
                    val2 = 24.d0 * idist6 - 48.d0 * idist12
                    val3 = (DIJ - Param(8))
                    val4 = val3*val3
                    val5 = exp(-Param(7) * val4)
                    THISU = Eps * val1 + Param(6) * val5 + P1_UShift(0)
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(1) = Terms(1) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = Eps * val2 - 2.d0 * Param(6) * Param(7) * DIJ * val3 * val5
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(4) = DUParam(4) + val1 + P1_UShift(1)
                        val6 = 1.d0 / Sig
                        DUParam(5) = DUParam(5) - Eps * val2 * val6 + P1_UShift(2)
                        DDUParam(14) = DDUParam(14) + Param(4) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P1_UShift(3)
                        DDUParam(9) = DDUParam(9) - val2 * val6 + P1_UShift(4)
                        DUParam(6) = DUParam(6) + val5 + P1_UShift(5)
                        val6 = -val4 * val5
                        DUParam(7) = DUParam(7) + Param(6) * val6 + P1_UShift(6)
                        DDUParam(25) = DDUParam(25) + val6 + P1_UShift(7)
                        val6 = 2.d0 * Param(7) * val3 * val5
                        DUParam(8) = DUParam(8) + Param(6) * val6 + P1_UShift(8)
                        DDUParam(30) = DDUParam(30) + val6 + P1_UShift(9)
                        DDUParam(26) = DDUParam(26) + Param(6) * val4*val4 * val5 + P1_UShift(10)
                        val6 = 2.d0*Param(6) * val5
                        DDUParam(31) = DDUParam(31) + val6 * val3 * (1.d0 - val4*Param(7)) + P1_UShift(11)
                        DDUParam(32) = DDUParam(32) + val6 * Param(7) * (2.d0 * Param(7) * val4 - 1.d0) + P1_UShift(12)
                    endif
                    if (CALCDWPARAM) then
                        DWParam(4) = DWParam(4) + val2
                        val6 = 1.d0 / Sig
                        val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                        DWParam(5) = DWParam(5) + Eps *  val7
                        DDWParam(14) = DDWParam(14) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                        DDWParam(9) = DDWParam(9) + val7
                        val7 = 2.d0 * DIJ * val3 * val5
                        DWParam(6) = DWParam(6) - val7 * Param(7)
                        val6 = val4 * Param(7)
                        DWParam(7) = DWParam(7) + Param(6) * val7 * (val6 - 1.d0)
                        DDWParam(26) = DDWParam(26) - Param(6) * val7 * val4 * (val6 - 2.d0)
                        DDWParam(25) = DDWParam(25) + val7 * (val6 - 1.d0)
                        val7 = 2.d0 * DIJ * val5
                        DWParam(8) = DWParam(8) - val7 * Param(6) * Param(7) * (2.d0 * Param(7) * val4 - 1.d0)
                        DDWParam(27) = DDWParam(27) + Param(6) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                        DDWParam(32) = DDWParam(32) - 2.d0 * val7 * Param(6) * Param(7)*Param(7) * val3 * (2.d0 * val6 - 3.d0)
                        DDWParam(22) = DDWParam(22) - val7 * Param(7) * (2.d0 * val6 - 1.d0)
                    endif
                end if
            end if

            if (dijsq < CutSq(2)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    !potential LJGaussA_A-
                    Sig = Param(10)
                    Eps = Param(9)
                    idist2 = Sig**2 / DIJSQ
                    idist6 = idist2 * idist2 * idist2
                    idist12 = idist6 * idist6
                    val1 = 4.d0 * (idist12 - idist6)
                    val2 = 24.d0 * idist6 - 48.d0 * idist12
                    val3 = (DIJ - Param(13))
                    val4 = val3*val3
                    val5 = exp(-Param(12) * val4)
                    THISU = Eps * val1 + Param(11) * val5 + P2_UShift(0)
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(2) = Terms(2) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = Eps * val2 - 2.d0 * Param(11) * Param(12) * DIJ * val3 * val5
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(9) = DUParam(9) + val1 + P2_UShift(1)
                        val6 = 1.d0 / Sig
                        DUParam(10) = DUParam(10) - Eps * val2 * val6 + P2_UShift(2)
                        DDUParam(39) = DDUParam(39) + Param(9) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P2_UShift(3)
                        DDUParam(34) = DDUParam(34) - val2 * val6 + P2_UShift(4)
                        DUParam(11) = DUParam(11) + val5 + P2_UShift(5)
                        val6 = -val4 * val5
                        DUParam(12) = DUParam(12) + Param(11) * val6 + P2_UShift(6)
                        DDUParam(50) = DDUParam(50) + val6 + P2_UShift(7)
                        val6 = 2.d0 * Param(12) * val3 * val5
                        DUParam(13) = DUParam(13) + Param(11) * val6 + P2_UShift(8)
                        DDUParam(55) = DDUParam(55) + val6 + P2_UShift(9)
                        DDUParam(51) = DDUParam(51) + Param(11) * val4*val4 * val5 + P2_UShift(10)
                        val6 = 2.d0*Param(11) * val5
                        DDUParam(56) = DDUParam(56) + val6 * val3 * (1.d0 - val4*Param(12)) + P2_UShift(11)
                        DDUParam(57) = DDUParam(57) + val6 * Param(12) * (2.d0 * Param(12) * val4 - 1.d0) + P2_UShift(12)
                    endif
                    if (CALCDWPARAM) then
                        DWParam(9) = DWParam(9) + val2
                        val6 = 1.d0 / Sig
                        val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                        DWParam(10) = DWParam(10) + Eps *  val7
                        DDWParam(39) = DDWParam(39) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                        DDWParam(34) = DDWParam(34) + val7
                        val7 = 2.d0 * DIJ * val3 * val5
                        DWParam(11) = DWParam(11) - val7 * Param(12)
                        val6 = val4 * Param(12)
                        DWParam(12) = DWParam(12) + Param(11) * val7 * (val6 - 1.d0)
                        DDWParam(51) = DDWParam(51) - Param(11) * val7 * val4 * (val6 - 2.d0)
                        DDWParam(50) = DDWParam(50) + val7 * (val6 - 1.d0)
                        val7 = 2.d0 * DIJ * val5
                        DWParam(13) = DWParam(13) - val7 * Param(11) * Param(12) * (2.d0 * Param(12) * val4 - 1.d0)
                        DDWParam(52) = DDWParam(52) + Param(11) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                        DDWParam(57) = DDWParam(57) - 2.d0 * val7 * Param(11) * Param(12)*Param(12) * val3 * (2.d0 * val6 - 3.d0)
                        DDWParam(47) = DDWParam(47) - val7 * Param(12) * (2.d0 * val6 - 1.d0)
                    endif
                end if
            end if

            if (dijsq < CutSq(3)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1)) .or. ((AIDj==1))) then
                    !potential LJGaussA_Na+
                    Sig = Param(15)
                    Eps = Param(14)
                    idist2 = Sig**2 / DIJSQ
                    idist6 = idist2 * idist2 * idist2
                    idist12 = idist6 * idist6
                    val1 = 4.d0 * (idist12 - idist6)
                    val2 = 24.d0 * idist6 - 48.d0 * idist12
                    val3 = (DIJ - Param(18))
                    val4 = val3*val3
                    val5 = exp(-Param(17) * val4)
                    THISU = Eps * val1 + Param(16) * val5 + P3_UShift(0)
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(3) = Terms(3) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = Eps * val2 - 2.d0 * Param(16) * Param(17) * DIJ * val3 * val5
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(14) = DUParam(14) + val1 + P3_UShift(1)
                        val6 = 1.d0 / Sig
                        DUParam(15) = DUParam(15) - Eps * val2 * val6 + P3_UShift(2)
                        DDUParam(64) = DDUParam(64) + Param(14) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P3_UShift(3)
                        DDUParam(59) = DDUParam(59) - val2 * val6 + P3_UShift(4)
                        DUParam(16) = DUParam(16) + val5 + P3_UShift(5)
                        val6 = -val4 * val5
                        DUParam(17) = DUParam(17) + Param(16) * val6 + P3_UShift(6)
                        DDUParam(75) = DDUParam(75) + val6 + P3_UShift(7)
                        val6 = 2.d0 * Param(17) * val3 * val5
                        DUParam(18) = DUParam(18) + Param(16) * val6 + P3_UShift(8)
                        DDUParam(80) = DDUParam(80) + val6 + P3_UShift(9)
                        DDUParam(76) = DDUParam(76) + Param(16) * val4*val4 * val5 + P3_UShift(10)
                        val6 = 2.d0*Param(16) * val5
                        DDUParam(81) = DDUParam(81) + val6 * val3 * (1.d0 - val4*Param(17)) + P3_UShift(11)
                        DDUParam(82) = DDUParam(82) + val6 * Param(17) * (2.d0 * Param(17) * val4 - 1.d0) + P3_UShift(12)
                    endif
                    if (CALCDWPARAM) then
                        DWParam(14) = DWParam(14) + val2
                        val6 = 1.d0 / Sig
                        val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                        DWParam(15) = DWParam(15) + Eps *  val7
                        DDWParam(64) = DDWParam(64) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                        DDWParam(59) = DDWParam(59) + val7
                        val7 = 2.d0 * DIJ * val3 * val5
                        DWParam(16) = DWParam(16) - val7 * Param(17)
                        val6 = val4 * Param(17)
                        DWParam(17) = DWParam(17) + Param(16) * val7 * (val6 - 1.d0)
                        DDWParam(76) = DDWParam(76) - Param(16) * val7 * val4 * (val6 - 2.d0)
                        DDWParam(75) = DDWParam(75) + val7 * (val6 - 1.d0)
                        val7 = 2.d0 * DIJ * val5
                        DWParam(18) = DWParam(18) - val7 * Param(16) * Param(17) * (2.d0 * Param(17) * val4 - 1.d0)
                        DDWParam(77) = DDWParam(77) + Param(16) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                        DDWParam(82) = DDWParam(82) - 2.d0 * val7 * Param(16) * Param(17)*Param(17) * val3 * (2.d0 * val6 - 3.d0)
                        DDWParam(72) = DDWParam(72) - val7 * Param(17) * (2.d0 * val6 - 1.d0)
                    endif
                end if
            end if

            if (dijsq < CutSq(4)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential LJGaussA-_A-
                    Sig = Param(20)
                    Eps = Param(19)
                    idist2 = Sig**2 / DIJSQ
                    idist6 = idist2 * idist2 * idist2
                    idist12 = idist6 * idist6
                    val1 = 4.d0 * (idist12 - idist6)
                    val2 = 24.d0 * idist6 - 48.d0 * idist12
                    val3 = (DIJ - Param(23))
                    val4 = val3*val3
                    val5 = exp(-Param(22) * val4)
                    THISU = Eps * val1 + Param(21) * val5 + P4_UShift(0)
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(4) = Terms(4) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = Eps * val2 - 2.d0 * Param(21) * Param(22) * DIJ * val3 * val5
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(19) = DUParam(19) + val1 + P4_UShift(1)
                        val6 = 1.d0 / Sig
                        DUParam(20) = DUParam(20) - Eps * val2 * val6 + P4_UShift(2)
                        DDUParam(89) = DDUParam(89) + Param(19) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P4_UShift(3)
                        DDUParam(84) = DDUParam(84) - val2 * val6 + P4_UShift(4)
                        DUParam(21) = DUParam(21) + val5 + P4_UShift(5)
                        val6 = -val4 * val5
                        DUParam(22) = DUParam(22) + Param(21) * val6 + P4_UShift(6)
                        DDUParam(100) = DDUParam(100) + val6 + P4_UShift(7)
                        val6 = 2.d0 * Param(22) * val3 * val5
                        DUParam(23) = DUParam(23) + Param(21) * val6 + P4_UShift(8)
                        DDUParam(105) = DDUParam(105) + val6 + P4_UShift(9)
                        DDUParam(101) = DDUParam(101) + Param(21) * val4*val4 * val5 + P4_UShift(10)
                        val6 = 2.d0*Param(21) * val5
                        DDUParam(106) = DDUParam(106) + val6 * val3 * (1.d0 - val4*Param(22)) + P4_UShift(11)
                        DDUParam(107) = DDUParam(107) + val6 * Param(22) * (2.d0 * Param(22) * val4 - 1.d0) + P4_UShift(12)
                    endif
                    if (CALCDWPARAM) then
                        DWParam(19) = DWParam(19) + val2
                        val6 = 1.d0 / Sig
                        val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                        DWParam(20) = DWParam(20) + Eps *  val7
                        DDWParam(89) = DDWParam(89) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                        DDWParam(84) = DDWParam(84) + val7
                        val7 = 2.d0 * DIJ * val3 * val5
                        DWParam(21) = DWParam(21) - val7 * Param(22)
                        val6 = val4 * Param(22)
                        DWParam(22) = DWParam(22) + Param(21) * val7 * (val6 - 1.d0)
                        DDWParam(101) = DDWParam(101) - Param(21) * val7 * val4 * (val6 - 2.d0)
                        DDWParam(100) = DDWParam(100) + val7 * (val6 - 1.d0)
                        val7 = 2.d0 * DIJ * val5
                        DWParam(23) = DWParam(23) - val7 * Param(21) * Param(22) * (2.d0 * Param(22) * val4 - 1.d0)
                        DDWParam(102) = DDWParam(102) + Param(21) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                        DDWParam(107) = DDWParam(107) - 2.d0 * val7 * Param(21) * Param(22)*Param(22) * val3 * (2.d0 * val6 - &
                          & 3.d0)
                        DDWParam(97) = DDWParam(97) - val7 * Param(22) * (2.d0 * val6 - 1.d0)
                    endif
                end if
            end if

            if (dijsq < CutSq(5)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0)) .or. ((AIDj==0))) then
                    !potential LJGaussA-_Na+
                    Sig = Param(25)
                    Eps = Param(24)
                    idist2 = Sig**2 / DIJSQ
                    idist6 = idist2 * idist2 * idist2
                    idist12 = idist6 * idist6
                    val1 = 4.d0 * (idist12 - idist6)
                    val2 = 24.d0 * idist6 - 48.d0 * idist12
                    val3 = (DIJ - Param(28))
                    val4 = val3*val3
                    val5 = exp(-Param(27) * val4)
                    THISU = Eps * val1 + Param(26) * val5 + P5_UShift(0)
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(5) = Terms(5) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = Eps * val2 - 2.d0 * Param(26) * Param(27) * DIJ * val3 * val5
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(24) = DUParam(24) + val1 + P5_UShift(1)
                        val6 = 1.d0 / Sig
                        DUParam(25) = DUParam(25) - Eps * val2 * val6 + P5_UShift(2)
                        DDUParam(114) = DDUParam(114) + Param(24) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P5_UShift(3)
                        DDUParam(109) = DDUParam(109) - val2 * val6 + P5_UShift(4)
                        DUParam(26) = DUParam(26) + val5 + P5_UShift(5)
                        val6 = -val4 * val5
                        DUParam(27) = DUParam(27) + Param(26) * val6 + P5_UShift(6)
                        DDUParam(125) = DDUParam(125) + val6 + P5_UShift(7)
                        val6 = 2.d0 * Param(27) * val3 * val5
                        DUParam(28) = DUParam(28) + Param(26) * val6 + P5_UShift(8)
                        DDUParam(130) = DDUParam(130) + val6 + P5_UShift(9)
                        DDUParam(126) = DDUParam(126) + Param(26) * val4*val4 * val5 + P5_UShift(10)
                        val6 = 2.d0*Param(26) * val5
                        DDUParam(131) = DDUParam(131) + val6 * val3 * (1.d0 - val4*Param(27)) + P5_UShift(11)
                        DDUParam(132) = DDUParam(132) + val6 * Param(27) * (2.d0 * Param(27) * val4 - 1.d0) + P5_UShift(12)
                    endif
                    if (CALCDWPARAM) then
                        DWParam(24) = DWParam(24) + val2
                        val6 = 1.d0 / Sig
                        val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                        DWParam(25) = DWParam(25) + Eps *  val7
                        DDWParam(114) = DDWParam(114) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                        DDWParam(109) = DDWParam(109) + val7
                        val7 = 2.d0 * DIJ * val3 * val5
                        DWParam(26) = DWParam(26) - val7 * Param(27)
                        val6 = val4 * Param(27)
                        DWParam(27) = DWParam(27) + Param(26) * val7 * (val6 - 1.d0)
                        DDWParam(126) = DDWParam(126) - Param(26) * val7 * val4 * (val6 - 2.d0)
                        DDWParam(125) = DDWParam(125) + val7 * (val6 - 1.d0)
                        val7 = 2.d0 * DIJ * val5
                        DWParam(28) = DWParam(28) - val7 * Param(26) * Param(27) * (2.d0 * Param(27) * val4 - 1.d0)
                        DDWParam(127) = DDWParam(127) + Param(26) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                        DDWParam(132) = DDWParam(132) - 2.d0 * val7 * Param(26) * Param(27)*Param(27) * val3 * (2.d0 * val6 - &
                          & 3.d0)
                        DDWParam(122) = DDWParam(122) - val7 * Param(27) * (2.d0 * val6 - 1.d0)
                    endif
                end if
            end if

            if (dijsq < CutSq(6)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential LJGaussNa+_Na+
                Sig = Param(30)
                Eps = Param(29)
                idist2 = Sig**2 / DIJSQ
                idist6 = idist2 * idist2 * idist2
                idist12 = idist6 * idist6
                val1 = 4.d0 * (idist12 - idist6)
                val2 = 24.d0 * idist6 - 48.d0 * idist12
                val3 = (DIJ - Param(33))
                val4 = val3*val3
                val5 = exp(-Param(32) * val4)
                THISU = Eps * val1 + Param(31) * val5 + P6_UShift(0)
                PEnergy = PEnergy + ThisU * Scale
                Terms(6) = Terms(6) + ThisU * Scale
                if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                    THISW = Eps * val2 - 2.d0 * Param(31) * Param(32) * DIJ * val3 * val5
                    Virial = Virial + ThisW * Scale
                endif
                if (CALCFORCE) then
                    FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                    FORCEI = FORCEI * SCALE
                    FORCE(I,:) = FORCE(I,:) + FORCEI
                    FORCE(J,:) = FORCE(J,:) - FORCEI
                endif
                if (CALCDUPARAM .or. CALCFLUCT) then
                    DUParam(29) = DUParam(29) + val1 + P6_UShift(1)
                    val6 = 1.d0 / Sig
                    DUParam(30) = DUParam(30) - Eps * val2 * val6 + P6_UShift(2)
                    DDUParam(139) = DDUParam(139) + Param(29) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P6_UShift(3)
                    DDUParam(134) = DDUParam(134) - val2 * val6 + P6_UShift(4)
                    DUParam(31) = DUParam(31) + val5 + P6_UShift(5)
                    val6 = -val4 * val5
                    DUParam(32) = DUParam(32) + Param(31) * val6 + P6_UShift(6)
                    DDUParam(150) = DDUParam(150) + val6 + P6_UShift(7)
                    val6 = 2.d0 * Param(32) * val3 * val5
                    DUParam(33) = DUParam(33) + Param(31) * val6 + P6_UShift(8)
                    DDUParam(155) = DDUParam(155) + val6 + P6_UShift(9)
                    DDUParam(151) = DDUParam(151) + Param(31) * val4*val4 * val5 + P6_UShift(10)
                    val6 = 2.d0*Param(31) * val5
                    DDUParam(156) = DDUParam(156) + val6 * val3 * (1.d0 - val4*Param(32)) + P6_UShift(11)
                    DDUParam(157) = DDUParam(157) + val6 * Param(32) * (2.d0 * Param(32) * val4 - 1.d0) + P6_UShift(12)
                endif
                if (CALCDWPARAM) then
                    DWParam(29) = DWParam(29) + val2
                    val6 = 1.d0 / Sig
                    val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                    DWParam(30) = DWParam(30) + Eps *  val7
                    DDWParam(139) = DDWParam(139) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                    DDWParam(134) = DDWParam(134) + val7
                    val7 = 2.d0 * DIJ * val3 * val5
                    DWParam(31) = DWParam(31) - val7 * Param(32)
                    val6 = val4 * Param(32)
                    DWParam(32) = DWParam(32) + Param(31) * val7 * (val6 - 1.d0)
                    DDWParam(151) = DDWParam(151) - Param(31) * val7 * val4 * (val6 - 2.d0)
                    DDWParam(150) = DDWParam(150) + val7 * (val6 - 1.d0)
                    val7 = 2.d0 * DIJ * val5
                    DWParam(33) = DWParam(33) - val7 * Param(31) * Param(32) * (2.d0 * Param(32) * val4 - 1.d0)
                    DDWParam(152) = DDWParam(152) + Param(31) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                    DDWParam(157) = DDWParam(157) - 2.d0 * val7 * Param(31) * Param(32)*Param(32) * val3 * (2.d0 * val6 - 3.d0)
                    DDWParam(147) = DDWParam(147) - val7 * Param(32) * (2.d0 * val6 - 1.d0)
                endif
            end if

            if (dijsq < CutSq(7)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential EW
                !pair interactions for potential <EW>
                Chargej = Param(AIDJ)
                call erfc(P7_Alpha * DIJ, erfcterm)
                if (BondOrdij > 0 .and. BondOrdij <= P7_ExcludeBondOrd) then
                    !exclude this interaction so subtract off 1/r
                    erfcterm = erfcterm - 1.0
                endif
                idist = 1.d0 / DIJ
                erfcterm = erfcterm * idist
                val1 = Param(34) * Chargei * Chargej
                val2 = erfcterm - P7_ShiftTerm
                THISU = val1 * val2
                PEnergy = PEnergy + ThisU * Scale
                Terms(7) = Terms(7) + ThisU * Scale
                if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                    val3 = -(erfcterm + P7_Term1 * exp(P7_Term2 * DIJSQ))
                    THISW = val1 * val3
                    Virial = Virial + ThisW * Scale
                    if (CALCFORCE) then
                        FORCEI = (val1 * val3 * idist*idist) * RIJ
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                endif
                if (CALCDUPARAM) then
                    DUParam(34) = DUParam(34) + Chargei * Chargej * val2
                    if (AIDI==AIDJ) then
                        DUParam(AIDI) = DUParam(AIDI) + 2.d0 * Param(34) * Chargej * val2
                        DDUParam(247 + AIDI + 0) = DDUParam(247 + AIDI + 0) + 2.d0 * Chargej * val2
                        DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + 2.d0 * Param(34) * val2
                    else
                        DUParam(AIDI) = DUParam(AIDI) + Param(34) * Chargej * val2
                        DUParam(AIDJ) = DUParam(AIDJ) + Param(34) * Chargei * val2
                        DDUParam(247 + AIDI + 0) = DDUParam(247 + AIDI + 0) + Chargej * val2
                        DDUParam(247 + AIDJ + 0) = DDUParam(247 + AIDJ + 0) + Chargei * val2
                        DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + Param(34) * val2
                    endif
                endif
                if (CALCDWPARAM) then
                    DWParam(34) = DWParam(34) + Chargei * Chargej * val3
                    if (AIDI==AIDJ) then
                        DWParam(AIDI) = DWParam(AIDI) + 2.d0 * Param(34) * Chargej * val3
                        DDWParam(247 + AIDI + 0) = DDWParam(247 + AIDI + 0) + 2.d0 * Chargej * val3
                        DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + 2.d0 * Param(34) * val3
                    else
                        DWParam(AIDI) = DWParam(AIDI) + Param(34) * Chargej * val3
                        DWParam(AIDJ) = DWParam(AIDJ) + Param(34) * Chargei * val3
                        DDWParam(247 + AIDI + 0) = DDWParam(247 + AIDI + 0) + Chargej * val3
                        DDWParam(247 + AIDJ + 0) = DDWParam(247 + AIDJ + 0) + Chargei * val3
                        DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + Param(34) * val3
                    endif
                endif
            end if

            if (dijsq < CutSq(8)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential SmearCoulA_A
                    !pair interactions for potential <SmearCoulA_A>
                    Chargej = Param(AIDJ)
                    idist = 1.d0 / DIJ
                    fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(36)
                    fac2 = fac*fac
                    val1 = Chargei * Chargej
                    call erf(fac*DIJ, tmp)
                    val2 = tmp * idist
                    val3 = Param(35) * val1
                    THISU = val3 * (val2 - idist + P8_UShift(0) + P8_CoulShift(0))
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(8) = Terms(8) + ThisU * Scale
                    iA = 1.d0 / Param(36)
                    val5 = exp(-fac2*DIJSQ)*iA*iA
                    val6 = exp(-fac2*DIJSQ)/Param(36)
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = val3 * ( val6 - val2 + idist)
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCDUPARAM) then
                        DUParam(35) = DUParam(35) + val1 * (val2 - idist + P8_UShift(0) + P8_CoulShift(0))
                        DUParam(36) = DUParam(36) + val3*(-val5 + P8_UShift(1)  )
                        DDUParam(162) = DDUParam(162) + val3 * ( val6*2.d0/Param(36)**2.d0 - &
                          & val6*2.d0/Param(36)**2.d0*fac2*DIJSQ  + P8_UShift(2)  )
                        DDUParam(160) = DDUParam(160) + val1*(-val5 + P8_UShift(1) )
                        if (AIDI==AIDJ) then
                            DUParam(AIDI) = DUParam(AIDI) + 2.d0 * Param(35) * Chargej * (val2 + P8_UShift(0) - idist + &
                              & P8_CoulShift(0))
                            DDUParam(251 + AIDI + 0) = DDUParam(251 + AIDI + 0) + 2.d0 * Param(35) * Chargej * (-val5 + &
                              & P8_UShift(1))
                            DDUParam(249 + AIDI + 0) = DDUParam(249 + AIDI + 0) + 2.d0 * Chargej * (val2 + P8_UShift(0) - idist &
                              & + P8_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + 2.d0 * Param(35) * (val2 + P8_UShift(0) - idist &
                              & + P8_CoulShift(0))
                        else
                            DUParam(AIDI) = DUParam(AIDI) + Param(35) * Chargej * (val2 + P8_UShift(0) -idist + P8_CoulShift(0))
                            DUParam(AIDJ) = DUParam(AIDJ) + Param(35) * Chargei * (val2 + P8_UShift(0) -idist + P8_CoulShift(0))
                            DDUParam(251 + AIDI + 0) = DDUParam(251 + AIDI + 0) + Param(35) * Chargej * (-val5 + P8_UShift(1))
                            DDUParam(251 + AIDJ + 0) = DDUParam(251 + AIDJ + 0) + Param(35) * Chargei * (-val5 + P8_UShift(1))
                            DDUParam(249 + AIDI + 0) = DDUParam(249 + AIDI + 0) + Chargej * (val2 + P8_UShift(0) - idist + &
                              & P8_CoulShift(0))
                            DDUParam(249 + AIDJ + 0) = DDUParam(249 + AIDJ + 0) + Chargei * (val2 + P8_UShift(0) - idist + &
                              & P8_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + Param(35) * (val2 + P8_UShift(0) - idist + &
                              & P8_CoulShift(0))
                        endif
                    endif
                    if (CALCDWPARAM) then

                        DWParam(35) = DWParam(35) + val1 * (val6 - val2 + idist)
                        val4 = 2.d0 * val5 * fac2 * DIJSQ
                        DWParam(36) = DWParam(36) + val3*val4

                        DDWParam(162) = DDWParam(162) + 2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4
                        DDWParam(160) = DDWParam(160) + val1 * val4

                        if (AIDI==AIDJ) then
                            DWParam(AIDI) = DWParam(AIDI) + 2.d0 * Chargej * Param(35) * (val6 - val2 + idist)
                            DDWParam(251 + AIDI + 0) = DDWParam(251 + AIDI + 0) + 2.d0 * Chargej * Param(35) * val4
                            DDWParam(249 + AIDI + 0) = DDWParam(249 + AIDI + 0) + 2.d0 * Chargej * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + 2.d0 * Param(35) * (val6 - val2 + idist)
                        else
                            DWParam(AIDI) = DWParam(AIDI) + Param(35) * Chargej * (val6 - val2 + idist)
                            DWParam(AIDJ) = DWParam(AIDJ) + Param(35) * Chargei * (val6 - val2 + idist)

                            DDWParam(251 + AIDI + 0) = DDWParam(251 + AIDI + 0) + Param(35) * val4 * Chargej
                            DDWParam(251 + AIDJ + 0) = DDWParam(251 + AIDJ + 0) + Param(35) * val4 * Chargei
                            DDWParam(249 + AIDI + 0) = DDWParam(249 + AIDI + 0) + Chargej * (val6 - val2 + idist)
                            DDWParam(249 + AIDJ + 0) = DDWParam(249 + AIDJ + 0) + Chargei * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + Param(35) * (val6 - val2 + idist)
                        endif
                    endif
                    if (CALCFORCE) then
                        FORCEI = RIJ * THISW / DIJSQ
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                end if
            end if

            if (dijsq < CutSq(9)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    !potential SmearCoulA_A-
                    !pair interactions for potential <SmearCoulA_A->
                    Chargej = Param(AIDJ)
                    idist = 1.d0 / DIJ
                    fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(38)
                    fac2 = fac*fac
                    val1 = Chargei * Chargej
                    call erf(fac*DIJ, tmp)
                    val2 = tmp * idist
                    val3 = Param(37) * val1
                    THISU = val3 * (val2 - idist + P9_UShift(0) + P9_CoulShift(0))
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(9) = Terms(9) + ThisU * Scale
                    iA = 1.d0 / Param(38)
                    val5 = exp(-fac2*DIJSQ)*iA*iA
                    val6 = exp(-fac2*DIJSQ)/Param(38)
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = val3 * ( val6 - val2 + idist)
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCDUPARAM) then
                        DUParam(37) = DUParam(37) + val1 * (val2 - idist + P9_UShift(0) + P9_CoulShift(0))
                        DUParam(38) = DUParam(38) + val3*(-val5 + P9_UShift(1)  )
                        DDUParam(166) = DDUParam(166) + val3 * ( val6*2.d0/Param(38)**2.d0 - &
                          & val6*2.d0/Param(38)**2.d0*fac2*DIJSQ  + P9_UShift(2)  )
                        DDUParam(164) = DDUParam(164) + val1*(-val5 + P9_UShift(1) )
                        if (AIDI==AIDJ) then
                            DUParam(AIDI) = DUParam(AIDI) + 2.d0 * Param(37) * Chargej * (val2 + P9_UShift(0) - idist + &
                              & P9_CoulShift(0))
                            DDUParam(255 + AIDI + 0) = DDUParam(255 + AIDI + 0) + 2.d0 * Param(37) * Chargej * (-val5 + &
                              & P9_UShift(1))
                            DDUParam(253 + AIDI + 0) = DDUParam(253 + AIDI + 0) + 2.d0 * Chargej * (val2 + P9_UShift(0) - idist &
                              & + P9_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + 2.d0 * Param(37) * (val2 + P9_UShift(0) - idist &
                              & + P9_CoulShift(0))
                        else
                            DUParam(AIDI) = DUParam(AIDI) + Param(37) * Chargej * (val2 + P9_UShift(0) -idist + P9_CoulShift(0))
                            DUParam(AIDJ) = DUParam(AIDJ) + Param(37) * Chargei * (val2 + P9_UShift(0) -idist + P9_CoulShift(0))
                            DDUParam(255 + AIDI + 0) = DDUParam(255 + AIDI + 0) + Param(37) * Chargej * (-val5 + P9_UShift(1))
                            DDUParam(255 + AIDJ + 0) = DDUParam(255 + AIDJ + 0) + Param(37) * Chargei * (-val5 + P9_UShift(1))
                            DDUParam(253 + AIDI + 0) = DDUParam(253 + AIDI + 0) + Chargej * (val2 + P9_UShift(0) - idist + &
                              & P9_CoulShift(0))
                            DDUParam(253 + AIDJ + 0) = DDUParam(253 + AIDJ + 0) + Chargei * (val2 + P9_UShift(0) - idist + &
                              & P9_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + Param(37) * (val2 + P9_UShift(0) - idist + &
                              & P9_CoulShift(0))
                        endif
                    endif
                    if (CALCDWPARAM) then

                        DWParam(37) = DWParam(37) + val1 * (val6 - val2 + idist)
                        val4 = 2.d0 * val5 * fac2 * DIJSQ
                        DWParam(38) = DWParam(38) + val3*val4

                        DDWParam(166) = DDWParam(166) + 2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4
                        DDWParam(164) = DDWParam(164) + val1 * val4

                        if (AIDI==AIDJ) then
                            DWParam(AIDI) = DWParam(AIDI) + 2.d0 * Chargej * Param(37) * (val6 - val2 + idist)
                            DDWParam(255 + AIDI + 0) = DDWParam(255 + AIDI + 0) + 2.d0 * Chargej * Param(37) * val4
                            DDWParam(253 + AIDI + 0) = DDWParam(253 + AIDI + 0) + 2.d0 * Chargej * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + 2.d0 * Param(37) * (val6 - val2 + idist)
                        else
                            DWParam(AIDI) = DWParam(AIDI) + Param(37) * Chargej * (val6 - val2 + idist)
                            DWParam(AIDJ) = DWParam(AIDJ) + Param(37) * Chargei * (val6 - val2 + idist)

                            DDWParam(255 + AIDI + 0) = DDWParam(255 + AIDI + 0) + Param(37) * val4 * Chargej
                            DDWParam(255 + AIDJ + 0) = DDWParam(255 + AIDJ + 0) + Param(37) * val4 * Chargei
                            DDWParam(253 + AIDI + 0) = DDWParam(253 + AIDI + 0) + Chargej * (val6 - val2 + idist)
                            DDWParam(253 + AIDJ + 0) = DDWParam(253 + AIDJ + 0) + Chargei * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + Param(37) * (val6 - val2 + idist)
                        endif
                    endif
                    if (CALCFORCE) then
                        FORCEI = RIJ * THISW / DIJSQ
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                end if
            end if

            if (dijsq < CutSq(10)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1)) .or. ((AIDj==1))) then
                    !potential SmearCoulA_Na+
                    !pair interactions for potential <SmearCoulA_Na+>
                    Chargej = Param(AIDJ)
                    idist = 1.d0 / DIJ
                    fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(40)
                    fac2 = fac*fac
                    val1 = Chargei * Chargej
                    call erf(fac*DIJ, tmp)
                    val2 = tmp * idist
                    val3 = Param(39) * val1
                    THISU = val3 * (val2 - idist + P10_UShift(0) + P10_CoulShift(0))
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(10) = Terms(10) + ThisU * Scale
                    iA = 1.d0 / Param(40)
                    val5 = exp(-fac2*DIJSQ)*iA*iA
                    val6 = exp(-fac2*DIJSQ)/Param(40)
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = val3 * ( val6 - val2 + idist)
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCDUPARAM) then
                        DUParam(39) = DUParam(39) + val1 * (val2 - idist + P10_UShift(0) + P10_CoulShift(0))
                        DUParam(40) = DUParam(40) + val3*(-val5 + P10_UShift(1)  )
                        DDUParam(170) = DDUParam(170) + val3 * ( val6*2.d0/Param(40)**2.d0 - &
                          & val6*2.d0/Param(40)**2.d0*fac2*DIJSQ  + P10_UShift(2)  )
                        DDUParam(168) = DDUParam(168) + val1*(-val5 + P10_UShift(1) )
                        if (AIDI==AIDJ) then
                            DUParam(AIDI) = DUParam(AIDI) + 2.d0 * Param(39) * Chargej * (val2 + P10_UShift(0) - idist + &
                              & P10_CoulShift(0))
                            DDUParam(259 + AIDI + 0) = DDUParam(259 + AIDI + 0) + 2.d0 * Param(39) * Chargej * (-val5 + &
                              & P10_UShift(1))
                            DDUParam(257 + AIDI + 0) = DDUParam(257 + AIDI + 0) + 2.d0 * Chargej * (val2 + P10_UShift(0) - &
                              & idist + P10_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + 2.d0 * Param(39) * (val2 + P10_UShift(0) - &
                              & idist + P10_CoulShift(0))
                        else
                            DUParam(AIDI) = DUParam(AIDI) + Param(39) * Chargej * (val2 + P10_UShift(0) -idist + &
                              & P10_CoulShift(0))
                            DUParam(AIDJ) = DUParam(AIDJ) + Param(39) * Chargei * (val2 + P10_UShift(0) -idist + &
                              & P10_CoulShift(0))
                            DDUParam(259 + AIDI + 0) = DDUParam(259 + AIDI + 0) + Param(39) * Chargej * (-val5 + P10_UShift(1))
                            DDUParam(259 + AIDJ + 0) = DDUParam(259 + AIDJ + 0) + Param(39) * Chargei * (-val5 + P10_UShift(1))
                            DDUParam(257 + AIDI + 0) = DDUParam(257 + AIDI + 0) + Chargej * (val2 + P10_UShift(0) - idist + &
                              & P10_CoulShift(0))
                            DDUParam(257 + AIDJ + 0) = DDUParam(257 + AIDJ + 0) + Chargei * (val2 + P10_UShift(0) - idist + &
                              & P10_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + Param(39) * (val2 + P10_UShift(0) - idist + &
                              & P10_CoulShift(0))
                        endif
                    endif
                    if (CALCDWPARAM) then

                        DWParam(39) = DWParam(39) + val1 * (val6 - val2 + idist)
                        val4 = 2.d0 * val5 * fac2 * DIJSQ
                        DWParam(40) = DWParam(40) + val3*val4

                        DDWParam(170) = DDWParam(170) + 2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4
                        DDWParam(168) = DDWParam(168) + val1 * val4

                        if (AIDI==AIDJ) then
                            DWParam(AIDI) = DWParam(AIDI) + 2.d0 * Chargej * Param(39) * (val6 - val2 + idist)
                            DDWParam(259 + AIDI + 0) = DDWParam(259 + AIDI + 0) + 2.d0 * Chargej * Param(39) * val4
                            DDWParam(257 + AIDI + 0) = DDWParam(257 + AIDI + 0) + 2.d0 * Chargej * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + 2.d0 * Param(39) * (val6 - val2 + idist)
                        else
                            DWParam(AIDI) = DWParam(AIDI) + Param(39) * Chargej * (val6 - val2 + idist)
                            DWParam(AIDJ) = DWParam(AIDJ) + Param(39) * Chargei * (val6 - val2 + idist)

                            DDWParam(259 + AIDI + 0) = DDWParam(259 + AIDI + 0) + Param(39) * val4 * Chargej
                            DDWParam(259 + AIDJ + 0) = DDWParam(259 + AIDJ + 0) + Param(39) * val4 * Chargei
                            DDWParam(257 + AIDI + 0) = DDWParam(257 + AIDI + 0) + Chargej * (val6 - val2 + idist)
                            DDWParam(257 + AIDJ + 0) = DDWParam(257 + AIDJ + 0) + Chargei * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + Param(39) * (val6 - val2 + idist)
                        endif
                    endif
                    if (CALCFORCE) then
                        FORCEI = RIJ * THISW / DIJSQ
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                end if
            end if

            if (dijsq < CutSq(11)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential SmearCoulA-_A-
                    !pair interactions for potential <SmearCoulA-_A->
                    Chargej = Param(AIDJ)
                    idist = 1.d0 / DIJ
                    fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(42)
                    fac2 = fac*fac
                    val1 = Chargei * Chargej
                    call erf(fac*DIJ, tmp)
                    val2 = tmp * idist
                    val3 = Param(41) * val1
                    THISU = val3 * (val2 - idist + P11_UShift(0) + P11_CoulShift(0))
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(11) = Terms(11) + ThisU * Scale
                    iA = 1.d0 / Param(42)
                    val5 = exp(-fac2*DIJSQ)*iA*iA
                    val6 = exp(-fac2*DIJSQ)/Param(42)
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = val3 * ( val6 - val2 + idist)
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCDUPARAM) then
                        DUParam(41) = DUParam(41) + val1 * (val2 - idist + P11_UShift(0) + P11_CoulShift(0))
                        DUParam(42) = DUParam(42) + val3*(-val5 + P11_UShift(1)  )
                        DDUParam(174) = DDUParam(174) + val3 * ( val6*2.d0/Param(42)**2.d0 - &
                          & val6*2.d0/Param(42)**2.d0*fac2*DIJSQ  + P11_UShift(2)  )
                        DDUParam(172) = DDUParam(172) + val1*(-val5 + P11_UShift(1) )
                        if (AIDI==AIDJ) then
                            DUParam(AIDI) = DUParam(AIDI) + 2.d0 * Param(41) * Chargej * (val2 + P11_UShift(0) - idist + &
                              & P11_CoulShift(0))
                            DDUParam(263 + AIDI + 0) = DDUParam(263 + AIDI + 0) + 2.d0 * Param(41) * Chargej * (-val5 + &
                              & P11_UShift(1))
                            DDUParam(261 + AIDI + 0) = DDUParam(261 + AIDI + 0) + 2.d0 * Chargej * (val2 + P11_UShift(0) - &
                              & idist + P11_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + 2.d0 * Param(41) * (val2 + P11_UShift(0) - &
                              & idist + P11_CoulShift(0))
                        else
                            DUParam(AIDI) = DUParam(AIDI) + Param(41) * Chargej * (val2 + P11_UShift(0) -idist + &
                              & P11_CoulShift(0))
                            DUParam(AIDJ) = DUParam(AIDJ) + Param(41) * Chargei * (val2 + P11_UShift(0) -idist + &
                              & P11_CoulShift(0))
                            DDUParam(263 + AIDI + 0) = DDUParam(263 + AIDI + 0) + Param(41) * Chargej * (-val5 + P11_UShift(1))
                            DDUParam(263 + AIDJ + 0) = DDUParam(263 + AIDJ + 0) + Param(41) * Chargei * (-val5 + P11_UShift(1))
                            DDUParam(261 + AIDI + 0) = DDUParam(261 + AIDI + 0) + Chargej * (val2 + P11_UShift(0) - idist + &
                              & P11_CoulShift(0))
                            DDUParam(261 + AIDJ + 0) = DDUParam(261 + AIDJ + 0) + Chargei * (val2 + P11_UShift(0) - idist + &
                              & P11_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + Param(41) * (val2 + P11_UShift(0) - idist + &
                              & P11_CoulShift(0))
                        endif
                    endif
                    if (CALCDWPARAM) then

                        DWParam(41) = DWParam(41) + val1 * (val6 - val2 + idist)
                        val4 = 2.d0 * val5 * fac2 * DIJSQ
                        DWParam(42) = DWParam(42) + val3*val4

                        DDWParam(174) = DDWParam(174) + 2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4
                        DDWParam(172) = DDWParam(172) + val1 * val4

                        if (AIDI==AIDJ) then
                            DWParam(AIDI) = DWParam(AIDI) + 2.d0 * Chargej * Param(41) * (val6 - val2 + idist)
                            DDWParam(263 + AIDI + 0) = DDWParam(263 + AIDI + 0) + 2.d0 * Chargej * Param(41) * val4
                            DDWParam(261 + AIDI + 0) = DDWParam(261 + AIDI + 0) + 2.d0 * Chargej * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + 2.d0 * Param(41) * (val6 - val2 + idist)
                        else
                            DWParam(AIDI) = DWParam(AIDI) + Param(41) * Chargej * (val6 - val2 + idist)
                            DWParam(AIDJ) = DWParam(AIDJ) + Param(41) * Chargei * (val6 - val2 + idist)

                            DDWParam(263 + AIDI + 0) = DDWParam(263 + AIDI + 0) + Param(41) * val4 * Chargej
                            DDWParam(263 + AIDJ + 0) = DDWParam(263 + AIDJ + 0) + Param(41) * val4 * Chargei
                            DDWParam(261 + AIDI + 0) = DDWParam(261 + AIDI + 0) + Chargej * (val6 - val2 + idist)
                            DDWParam(261 + AIDJ + 0) = DDWParam(261 + AIDJ + 0) + Chargei * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + Param(41) * (val6 - val2 + idist)
                        endif
                    endif
                    if (CALCFORCE) then
                        FORCEI = RIJ * THISW / DIJSQ
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                end if
            end if

            if (dijsq < CutSq(12)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0)) .or. ((AIDj==0))) then
                    !potential SmearCoulA-_Na+
                    !pair interactions for potential <SmearCoulA-_Na+>
                    Chargej = Param(AIDJ)
                    idist = 1.d0 / DIJ
                    fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(44)
                    fac2 = fac*fac
                    val1 = Chargei * Chargej
                    call erf(fac*DIJ, tmp)
                    val2 = tmp * idist
                    val3 = Param(43) * val1
                    THISU = val3 * (val2 - idist + P12_UShift(0) + P12_CoulShift(0))
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(12) = Terms(12) + ThisU * Scale
                    iA = 1.d0 / Param(44)
                    val5 = exp(-fac2*DIJSQ)*iA*iA
                    val6 = exp(-fac2*DIJSQ)/Param(44)
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = val3 * ( val6 - val2 + idist)
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCDUPARAM) then
                        DUParam(43) = DUParam(43) + val1 * (val2 - idist + P12_UShift(0) + P12_CoulShift(0))
                        DUParam(44) = DUParam(44) + val3*(-val5 + P12_UShift(1)  )
                        DDUParam(178) = DDUParam(178) + val3 * ( val6*2.d0/Param(44)**2.d0 - &
                          & val6*2.d0/Param(44)**2.d0*fac2*DIJSQ  + P12_UShift(2)  )
                        DDUParam(176) = DDUParam(176) + val1*(-val5 + P12_UShift(1) )
                        if (AIDI==AIDJ) then
                            DUParam(AIDI) = DUParam(AIDI) + 2.d0 * Param(43) * Chargej * (val2 + P12_UShift(0) - idist + &
                              & P12_CoulShift(0))
                            DDUParam(267 + AIDI + 0) = DDUParam(267 + AIDI + 0) + 2.d0 * Param(43) * Chargej * (-val5 + &
                              & P12_UShift(1))
                            DDUParam(265 + AIDI + 0) = DDUParam(265 + AIDI + 0) + 2.d0 * Chargej * (val2 + P12_UShift(0) - &
                              & idist + P12_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + 2.d0 * Param(43) * (val2 + P12_UShift(0) - &
                              & idist + P12_CoulShift(0))
                        else
                            DUParam(AIDI) = DUParam(AIDI) + Param(43) * Chargej * (val2 + P12_UShift(0) -idist + &
                              & P12_CoulShift(0))
                            DUParam(AIDJ) = DUParam(AIDJ) + Param(43) * Chargei * (val2 + P12_UShift(0) -idist + &
                              & P12_CoulShift(0))
                            DDUParam(267 + AIDI + 0) = DDUParam(267 + AIDI + 0) + Param(43) * Chargej * (-val5 + P12_UShift(1))
                            DDUParam(267 + AIDJ + 0) = DDUParam(267 + AIDJ + 0) + Param(43) * Chargei * (-val5 + P12_UShift(1))
                            DDUParam(265 + AIDI + 0) = DDUParam(265 + AIDI + 0) + Chargej * (val2 + P12_UShift(0) - idist + &
                              & P12_CoulShift(0))
                            DDUParam(265 + AIDJ + 0) = DDUParam(265 + AIDJ + 0) + Chargei * (val2 + P12_UShift(0) - idist + &
                              & P12_CoulShift(0))
                            DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + Param(43) * (val2 + P12_UShift(0) - idist + &
                              & P12_CoulShift(0))
                        endif
                    endif
                    if (CALCDWPARAM) then

                        DWParam(43) = DWParam(43) + val1 * (val6 - val2 + idist)
                        val4 = 2.d0 * val5 * fac2 * DIJSQ
                        DWParam(44) = DWParam(44) + val3*val4

                        DDWParam(178) = DDWParam(178) + 2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4
                        DDWParam(176) = DDWParam(176) + val1 * val4

                        if (AIDI==AIDJ) then
                            DWParam(AIDI) = DWParam(AIDI) + 2.d0 * Chargej * Param(43) * (val6 - val2 + idist)
                            DDWParam(267 + AIDI + 0) = DDWParam(267 + AIDI + 0) + 2.d0 * Chargej * Param(43) * val4
                            DDWParam(265 + AIDI + 0) = DDWParam(265 + AIDI + 0) + 2.d0 * Chargej * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + 2.d0 * Param(43) * (val6 - val2 + idist)
                        else
                            DWParam(AIDI) = DWParam(AIDI) + Param(43) * Chargej * (val6 - val2 + idist)
                            DWParam(AIDJ) = DWParam(AIDJ) + Param(43) * Chargei * (val6 - val2 + idist)

                            DDWParam(267 + AIDI + 0) = DDWParam(267 + AIDI + 0) + Param(43) * val4 * Chargej
                            DDWParam(267 + AIDJ + 0) = DDWParam(267 + AIDJ + 0) + Param(43) * val4 * Chargei
                            DDWParam(265 + AIDI + 0) = DDWParam(265 + AIDI + 0) + Chargej * (val6 - val2 + idist)
                            DDWParam(265 + AIDJ + 0) = DDWParam(265 + AIDJ + 0) + Chargei * (val6 - val2 + idist)
                            DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + Param(43) * (val6 - val2 + idist)
                        endif
                    endif
                    if (CALCFORCE) then
                        FORCEI = RIJ * THISW / DIJSQ
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                end if
            end if

            if (dijsq < CutSq(13)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential SmearCoulNa+_Na+
                !pair interactions for potential <SmearCoulNa+_Na+>
                Chargej = Param(AIDJ)
                idist = 1.d0 / DIJ
                fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(46)
                fac2 = fac*fac
                val1 = Chargei * Chargej
                call erf(fac*DIJ, tmp)
                val2 = tmp * idist
                val3 = Param(45) * val1
                THISU = val3 * (val2 - idist + P13_UShift(0) + P13_CoulShift(0))
                PEnergy = PEnergy + ThisU * Scale
                Terms(13) = Terms(13) + ThisU * Scale
                iA = 1.d0 / Param(46)
                val5 = exp(-fac2*DIJSQ)*iA*iA
                val6 = exp(-fac2*DIJSQ)/Param(46)
                if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                    THISW = val3 * ( val6 - val2 + idist)
                    Virial = Virial + ThisW * Scale
                endif
                if (CALCDUPARAM) then
                    DUParam(45) = DUParam(45) + val1 * (val2 - idist + P13_UShift(0) + P13_CoulShift(0))
                    DUParam(46) = DUParam(46) + val3*(-val5 + P13_UShift(1)  )
                    DDUParam(182) = DDUParam(182) + val3 * ( val6*2.d0/Param(46)**2.d0 - val6*2.d0/Param(46)**2.d0*fac2*DIJSQ   &
                      & + P13_UShift(2)  )
                    DDUParam(180) = DDUParam(180) + val1*(-val5 + P13_UShift(1) )
                    if (AIDI==AIDJ) then
                        DUParam(AIDI) = DUParam(AIDI) + 2.d0 * Param(45) * Chargej * (val2 + P13_UShift(0) - idist + &
                          & P13_CoulShift(0))
                        DDUParam(271 + AIDI + 0) = DDUParam(271 + AIDI + 0) + 2.d0 * Param(45) * Chargej * (-val5 + &
                          & P13_UShift(1))
                        DDUParam(269 + AIDI + 0) = DDUParam(269 + AIDI + 0) + 2.d0 * Chargej * (val2 + P13_UShift(0) - idist + &
                          & P13_CoulShift(0))
                        DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + 2.d0 * Param(45) * (val2 + P13_UShift(0) - idist + &
                          & P13_CoulShift(0))
                    else
                        DUParam(AIDI) = DUParam(AIDI) + Param(45) * Chargej * (val2 + P13_UShift(0) -idist + P13_CoulShift(0))
                        DUParam(AIDJ) = DUParam(AIDJ) + Param(45) * Chargei * (val2 + P13_UShift(0) -idist + P13_CoulShift(0))
                        DDUParam(271 + AIDI + 0) = DDUParam(271 + AIDI + 0) + Param(45) * Chargej * (-val5 + P13_UShift(1))
                        DDUParam(271 + AIDJ + 0) = DDUParam(271 + AIDJ + 0) + Param(45) * Chargei * (-val5 + P13_UShift(1))
                        DDUParam(269 + AIDI + 0) = DDUParam(269 + AIDI + 0) + Chargej * (val2 + P13_UShift(0) - idist + &
                          & P13_CoulShift(0))
                        DDUParam(269 + AIDJ + 0) = DDUParam(269 + AIDJ + 0) + Chargei * (val2 + P13_UShift(0) - idist + &
                          & P13_CoulShift(0))
                        DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + Param(45) * (val2 + P13_UShift(0) - idist + &
                          & P13_CoulShift(0))
                    endif
                endif
                if (CALCDWPARAM) then

                    DWParam(45) = DWParam(45) + val1 * (val6 - val2 + idist)
                    val4 = 2.d0 * val5 * fac2 * DIJSQ
                    DWParam(46) = DWParam(46) + val3*val4

                    DDWParam(182) = DDWParam(182) + 2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4
                    DDWParam(180) = DDWParam(180) + val1 * val4

                    if (AIDI==AIDJ) then
                        DWParam(AIDI) = DWParam(AIDI) + 2.d0 * Chargej * Param(45) * (val6 - val2 + idist)
                        DDWParam(271 + AIDI + 0) = DDWParam(271 + AIDI + 0) + 2.d0 * Chargej * Param(45) * val4
                        DDWParam(269 + AIDI + 0) = DDWParam(269 + AIDI + 0) + 2.d0 * Chargej * (val6 - val2 + idist)
                        DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + 2.d0 * Param(45) * (val6 - val2 + idist)
                    else
                        DWParam(AIDI) = DWParam(AIDI) + Param(45) * Chargej * (val6 - val2 + idist)
                        DWParam(AIDJ) = DWParam(AIDJ) + Param(45) * Chargei * (val6 - val2 + idist)

                        DDWParam(271 + AIDI + 0) = DDWParam(271 + AIDI + 0) + Param(45) * val4 * Chargej
                        DDWParam(271 + AIDJ + 0) = DDWParam(271 + AIDJ + 0) + Param(45) * val4 * Chargei
                        DDWParam(269 + AIDI + 0) = DDWParam(269 + AIDI + 0) + Chargej * (val6 - val2 + idist)
                        DDWParam(269 + AIDJ + 0) = DDWParam(269 + AIDJ + 0) + Chargei * (val6 - val2 + idist)
                        DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + Param(45) * (val6 - val2 + idist)
                    endif
                endif
                if (CALCFORCE) then
                    FORCEI = RIJ * THISW / DIJSQ
                    FORCE(I,:) = FORCE(I,:) + FORCEI
                    FORCE(J,:) = FORCE(J,:) - FORCEI
                endif
            end if

        !end of loop j
        enddo

    !end of loop i
    enddo

end subroutine


subroutine saveenergystate(Mode)
    implicit none
    integer, intent(in) :: Mode
    integer :: istart
    integer :: istop

    if (Mode == 0) then
        !save for all atoms
        istart = 0
        istop = NAtom - 1
    elseif (Mode == 1) then
        !save for single atom TargetAtom
        istart = TargetAtom
        istop = TargetAtom
    elseif (Mode == 2) then
        !save for single molecule TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
    else
        print *, "Illegal value for Mode in SaveEnergyState."
    endif
    OldPEnergy = PEnergy
    OldTerms = Terms

    !potential EW
    P7_OldRPEnergy = P7_RPEnergy
    P7_OldRVirial = P7_RVirial
    P7_OldSPEnergy = P7_SPEnergy
    P7_OldStruc(istart:istop) = P7_Struc(istart:istop)
    P7_oldeik_x(:,istart:istop) = P7_eik_x(:,istart:istop)
    P7_oldeik_y(:,istart:istop) = P7_eik_y(:,istart:istop)
    P7_oldeik_z(:,istart:istop) = P7_eik_z(:,istart:istop)

end subroutine


subroutine revertenergystate(Mode)
    implicit none
    integer, intent(in) :: Mode
    integer :: istart
    integer :: istop

    if (Mode == 0) then
        !revert for all atoms
        istart = 0
        istop = NAtom - 1
    elseif (Mode == 1) then
        !revert for single atom TargetAtom
        istart = TargetAtom
        istop = TargetAtom
    elseif (Mode == 2) then
        !revert for single molecule TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
    else
        print *, "Illegal value for Mode in RevertEnergyState."
    endif
    PEnergy = OldPEnergy
    Terms = OldTerms

    !potential EW
    P7_RPEnergy = P7_OldRPEnergy
    P7_RVirial = P7_OldRVirial
    P7_SPEnergy = P7_OldSPEnergy
    P7_Struc(istart:istop) = P7_OldStruc(istart:istop)
    P7_eik_x(:,istart:istop) = P7_oldeik_x(:,istart:istop)
    P7_eik_y(:,istart:istop) = P7_oldeik_y(:,istart:istop)
    P7_eik_z(:,istart:istop) = P7_oldeik_z(:,istart:istop)

end subroutine


subroutine calcargstats(Weight)
    implicit none
    real(8), intent(in) :: Weight
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: ForceDoMInd
    integer :: ForceDoAtom
    integer :: AIDi
    integer :: AIDj
    integer :: AIDij
    integer :: SIDi
    integer :: SIDj
    integer :: MIndi
    integer :: MIndj
    integer :: BondOrdIndShift
    integer :: BondOrdSIDjStart
    integer :: BondOrdSIDjStop
    integer :: BondOrdij
    integer :: istart
    integer :: istop
    integer :: jstart
    integer :: LoopMode
    logical :: SameMol
    logical :: Bonded
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8) :: Scale
    integer :: ArgType
    real(8) :: ArgVal

    Scale = 1.d0
    LoopMode = 0

    !potential BondA_A-
    P0_ArgWeightSumStats = P0_ArgWeightSumStats + Weight
    !potential LJGaussA_A
    P1_ArgWeightSumStats = P1_ArgWeightSumStats + Weight
    !potential LJGaussA_A-
    P2_ArgWeightSumStats = P2_ArgWeightSumStats + Weight
    !potential LJGaussA_Na+
    P3_ArgWeightSumStats = P3_ArgWeightSumStats + Weight
    !potential LJGaussA-_A-
    P4_ArgWeightSumStats = P4_ArgWeightSumStats + Weight
    !potential LJGaussA-_Na+
    P5_ArgWeightSumStats = P5_ArgWeightSumStats + Weight
    !potential LJGaussNa+_Na+
    P6_ArgWeightSumStats = P6_ArgWeightSumStats + Weight
    !potential EW
    P7_ArgWeightSumStats = P7_ArgWeightSumStats + Weight
    !potential SmearCoulA_A
    P8_ArgWeightSumStats = P8_ArgWeightSumStats + Weight
    !potential SmearCoulA_A-
    P9_ArgWeightSumStats = P9_ArgWeightSumStats + Weight
    !potential SmearCoulA_Na+
    P10_ArgWeightSumStats = P10_ArgWeightSumStats + Weight
    !potential SmearCoulA-_A-
    P11_ArgWeightSumStats = P11_ArgWeightSumStats + Weight
    !potential SmearCoulA-_Na+
    P12_ArgWeightSumStats = P12_ArgWeightSumStats + Weight
    !potential SmearCoulNa+_Na+
    P13_ArgWeightSumStats = P13_ArgWeightSumStats + Weight

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    ForceDoMInd = -1
    ForceDoAtom = -1
    if (LoopMode == 0) then
        !normal full atom loop
        istart = 0
        istop = NAtom - 1
    elseif (LoopMode == 1) then
        !single atom interactions for adding TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == -1) then
        !single atom interactions for deleting TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == 2) then
        !molecule interactions fo adding TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    elseif (LoopMode == -2) then
        !molecule interactions fo deleting TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    else
        print *, "Illegal value for LoopMode."
    endif

    !loop over i
    do i = istart, istop

        MIndi = MInd(i)

        if (MolActive(MIndi) < 0 .and. LoopMode == 0) cycle

        Posi = Pos(i,:)
        AIDi = AID(i)
        SIDi = SID(i)
        BondOrdSIDjStart = BondOrdShift(SIDi)
        BondOrdSIDjStop = BondOrdSIDjStart + BondOrdStart(SIDi+1) - BondOrdStart(SIDi)
        BondOrdIndShift = BondOrdStart(SIDi) - BondOrdShift(SIDi)

        if (LoopMode == 0) then
            jstart = i+1
        else
            jstart = 0
        endif

        !loop over j
        do j = jstart, NAtom - 1

            !check to see if same atom
            if (i==j) cycle

            MIndj = MInd(j)

            if (LoopMode == 2 .or. LoopMode == -2) then
                !!check to see if we need to skip because of double counting
                if (MIndj == MIndi .and. j < i) cycle
            endif

            if (MolActive(MIndj) < 0 .and. MIndj /= ForceDoMInd .and. j /= ForceDoAtom) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)
            SIDj = SID(j)
            if (AIDi > AIDj) then
                AIDij = AIDi * (AIDi + 1) / 2 + AIDj
            else
                AIDij = AIDj * (AIDj + 1) / 2 + AIDi
            endif
            SameMol = (MIndi == MIndj)
            BondOrdij = BondOrdLimit + 1
            if (SameMol) then
                if (SIDj >= BondOrdSIDjStart .and. SIDj < BondOrdSIDjStop) then
                    BondOrdij = BondOrdData(SIDj + BondOrdIndShift)
                endif
            endif
            Bonded = (SameMol .and. BondOrdij==2)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            if (dijsq < CutSq(0)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if ((((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) .and. Bonded) then
                    !potential BondA_A-
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(0)) cycle
                    if (Weight > 0.d0) then
                        P0_ArgMin(ArgType) = min(P0_ArgMin(ArgType), ArgVal)
                        P0_ArgMax(ArgType) = max(P0_ArgMax(ArgType), ArgVal)
                        P0_ArgCount(ArgType) = P0_ArgCount(ArgType) + Weight
                        P0_ArgSum(ArgType) = P0_ArgSum(ArgType) + ArgVal * Weight
                        P0_ArgSumSq(ArgType) = P0_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(1)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential LJGaussA_A
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(1)) cycle
                    if (Weight > 0.d0) then
                        P1_ArgMin(ArgType) = min(P1_ArgMin(ArgType), ArgVal)
                        P1_ArgMax(ArgType) = max(P1_ArgMax(ArgType), ArgVal)
                        P1_ArgCount(ArgType) = P1_ArgCount(ArgType) + Weight
                        P1_ArgSum(ArgType) = P1_ArgSum(ArgType) + ArgVal * Weight
                        P1_ArgSumSq(ArgType) = P1_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(2)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    !potential LJGaussA_A-
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(2)) cycle
                    if (Weight > 0.d0) then
                        P2_ArgMin(ArgType) = min(P2_ArgMin(ArgType), ArgVal)
                        P2_ArgMax(ArgType) = max(P2_ArgMax(ArgType), ArgVal)
                        P2_ArgCount(ArgType) = P2_ArgCount(ArgType) + Weight
                        P2_ArgSum(ArgType) = P2_ArgSum(ArgType) + ArgVal * Weight
                        P2_ArgSumSq(ArgType) = P2_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(3)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1)) .or. ((AIDj==1))) then
                    !potential LJGaussA_Na+
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(3)) cycle
                    if (Weight > 0.d0) then
                        P3_ArgMin(ArgType) = min(P3_ArgMin(ArgType), ArgVal)
                        P3_ArgMax(ArgType) = max(P3_ArgMax(ArgType), ArgVal)
                        P3_ArgCount(ArgType) = P3_ArgCount(ArgType) + Weight
                        P3_ArgSum(ArgType) = P3_ArgSum(ArgType) + ArgVal * Weight
                        P3_ArgSumSq(ArgType) = P3_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(4)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential LJGaussA-_A-
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(4)) cycle
                    if (Weight > 0.d0) then
                        P4_ArgMin(ArgType) = min(P4_ArgMin(ArgType), ArgVal)
                        P4_ArgMax(ArgType) = max(P4_ArgMax(ArgType), ArgVal)
                        P4_ArgCount(ArgType) = P4_ArgCount(ArgType) + Weight
                        P4_ArgSum(ArgType) = P4_ArgSum(ArgType) + ArgVal * Weight
                        P4_ArgSumSq(ArgType) = P4_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(5)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0)) .or. ((AIDj==0))) then
                    !potential LJGaussA-_Na+
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(5)) cycle
                    if (Weight > 0.d0) then
                        P5_ArgMin(ArgType) = min(P5_ArgMin(ArgType), ArgVal)
                        P5_ArgMax(ArgType) = max(P5_ArgMax(ArgType), ArgVal)
                        P5_ArgCount(ArgType) = P5_ArgCount(ArgType) + Weight
                        P5_ArgSum(ArgType) = P5_ArgSum(ArgType) + ArgVal * Weight
                        P5_ArgSumSq(ArgType) = P5_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(6)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential LJGaussNa+_Na+
                ARGVAL = DIJ
                ARGTYPE = 0
                if (ArgVal*ArgVal > CutSq(6)) cycle
                if (Weight > 0.d0) then
                    P6_ArgMin(ArgType) = min(P6_ArgMin(ArgType), ArgVal)
                    P6_ArgMax(ArgType) = max(P6_ArgMax(ArgType), ArgVal)
                    P6_ArgCount(ArgType) = P6_ArgCount(ArgType) + Weight
                    P6_ArgSum(ArgType) = P6_ArgSum(ArgType) + ArgVal * Weight
                    P6_ArgSumSq(ArgType) = P6_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                endif
            end if

            if (dijsq < CutSq(7)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential EW
                ARGVAL = DIJ
                ARGTYPE = AIDIJ
                if (Weight > 0.d0) then
                    P7_ArgMin(ArgType) = min(P7_ArgMin(ArgType), ArgVal)
                    P7_ArgMax(ArgType) = max(P7_ArgMax(ArgType), ArgVal)
                    P7_ArgCount(ArgType) = P7_ArgCount(ArgType) + Weight
                    P7_ArgSum(ArgType) = P7_ArgSum(ArgType) + ArgVal * Weight
                    P7_ArgSumSq(ArgType) = P7_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                endif
            end if

            if (dijsq < CutSq(8)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential SmearCoulA_A
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (Weight > 0.d0) then
                        P8_ArgMin(ArgType) = min(P8_ArgMin(ArgType), ArgVal)
                        P8_ArgMax(ArgType) = max(P8_ArgMax(ArgType), ArgVal)
                        P8_ArgCount(ArgType) = P8_ArgCount(ArgType) + Weight
                        P8_ArgSum(ArgType) = P8_ArgSum(ArgType) + ArgVal * Weight
                        P8_ArgSumSq(ArgType) = P8_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(9)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    !potential SmearCoulA_A-
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (Weight > 0.d0) then
                        P9_ArgMin(ArgType) = min(P9_ArgMin(ArgType), ArgVal)
                        P9_ArgMax(ArgType) = max(P9_ArgMax(ArgType), ArgVal)
                        P9_ArgCount(ArgType) = P9_ArgCount(ArgType) + Weight
                        P9_ArgSum(ArgType) = P9_ArgSum(ArgType) + ArgVal * Weight
                        P9_ArgSumSq(ArgType) = P9_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(10)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1)) .or. ((AIDj==1))) then
                    !potential SmearCoulA_Na+
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (Weight > 0.d0) then
                        P10_ArgMin(ArgType) = min(P10_ArgMin(ArgType), ArgVal)
                        P10_ArgMax(ArgType) = max(P10_ArgMax(ArgType), ArgVal)
                        P10_ArgCount(ArgType) = P10_ArgCount(ArgType) + Weight
                        P10_ArgSum(ArgType) = P10_ArgSum(ArgType) + ArgVal * Weight
                        P10_ArgSumSq(ArgType) = P10_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(11)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential SmearCoulA-_A-
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (Weight > 0.d0) then
                        P11_ArgMin(ArgType) = min(P11_ArgMin(ArgType), ArgVal)
                        P11_ArgMax(ArgType) = max(P11_ArgMax(ArgType), ArgVal)
                        P11_ArgCount(ArgType) = P11_ArgCount(ArgType) + Weight
                        P11_ArgSum(ArgType) = P11_ArgSum(ArgType) + ArgVal * Weight
                        P11_ArgSumSq(ArgType) = P11_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(12)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0)) .or. ((AIDj==0))) then
                    !potential SmearCoulA-_Na+
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (Weight > 0.d0) then
                        P12_ArgMin(ArgType) = min(P12_ArgMin(ArgType), ArgVal)
                        P12_ArgMax(ArgType) = max(P12_ArgMax(ArgType), ArgVal)
                        P12_ArgCount(ArgType) = P12_ArgCount(ArgType) + Weight
                        P12_ArgSum(ArgType) = P12_ArgSum(ArgType) + ArgVal * Weight
                        P12_ArgSumSq(ArgType) = P12_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(13)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential SmearCoulNa+_Na+
                ARGVAL = DIJ
                ARGTYPE = AIDIJ
                if (Weight > 0.d0) then
                    P13_ArgMin(ArgType) = min(P13_ArgMin(ArgType), ArgVal)
                    P13_ArgMax(ArgType) = max(P13_ArgMax(ArgType), ArgVal)
                    P13_ArgCount(ArgType) = P13_ArgCount(ArgType) + Weight
                    P13_ArgSum(ArgType) = P13_ArgSum(ArgType) + ArgVal * Weight
                    P13_ArgSumSq(ArgType) = P13_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                endif
            end if

        !end of loop j
        enddo

    !end of loop i
    enddo

end subroutine


subroutine calcarghist(Weight)
    implicit none
    real(8), intent(in) :: Weight
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: ForceDoMInd
    integer :: ForceDoAtom
    integer :: AIDi
    integer :: AIDj
    integer :: AIDij
    integer :: SIDi
    integer :: SIDj
    integer :: MIndi
    integer :: MIndj
    integer :: BondOrdIndShift
    integer :: BondOrdSIDjStart
    integer :: BondOrdSIDjStop
    integer :: BondOrdij
    integer :: istart
    integer :: istop
    integer :: jstart
    integer :: LoopMode
    logical :: SameMol
    logical :: Bonded
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8) :: Scale
    integer :: ArgType
    real(8) :: ArgVal
    integer :: BinInd

    Scale = 1.d0
    LoopMode = 0

    !potential BondA_A-
    P0_ArgWeightSumHist = P0_ArgWeightSumHist + Weight
    !potential LJGaussA_A
    P1_ArgWeightSumHist = P1_ArgWeightSumHist + Weight
    !potential LJGaussA_A-
    P2_ArgWeightSumHist = P2_ArgWeightSumHist + Weight
    !potential LJGaussA_Na+
    P3_ArgWeightSumHist = P3_ArgWeightSumHist + Weight
    !potential LJGaussA-_A-
    P4_ArgWeightSumHist = P4_ArgWeightSumHist + Weight
    !potential LJGaussA-_Na+
    P5_ArgWeightSumHist = P5_ArgWeightSumHist + Weight
    !potential LJGaussNa+_Na+
    P6_ArgWeightSumHist = P6_ArgWeightSumHist + Weight
    !potential EW
    P7_ArgWeightSumHist = P7_ArgWeightSumHist + Weight
    !potential SmearCoulA_A
    P8_ArgWeightSumHist = P8_ArgWeightSumHist + Weight
    !potential SmearCoulA_A-
    P9_ArgWeightSumHist = P9_ArgWeightSumHist + Weight
    !potential SmearCoulA_Na+
    P10_ArgWeightSumHist = P10_ArgWeightSumHist + Weight
    !potential SmearCoulA-_A-
    P11_ArgWeightSumHist = P11_ArgWeightSumHist + Weight
    !potential SmearCoulA-_Na+
    P12_ArgWeightSumHist = P12_ArgWeightSumHist + Weight
    !potential SmearCoulNa+_Na+
    P13_ArgWeightSumHist = P13_ArgWeightSumHist + Weight

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    ForceDoMInd = -1
    ForceDoAtom = -1
    if (LoopMode == 0) then
        !normal full atom loop
        istart = 0
        istop = NAtom - 1
    elseif (LoopMode == 1) then
        !single atom interactions for adding TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == -1) then
        !single atom interactions for deleting TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == 2) then
        !molecule interactions fo adding TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    elseif (LoopMode == -2) then
        !molecule interactions fo deleting TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    else
        print *, "Illegal value for LoopMode."
    endif

    !loop over i
    do i = istart, istop

        MIndi = MInd(i)

        if (MolActive(MIndi) < 0 .and. LoopMode == 0) cycle

        Posi = Pos(i,:)
        AIDi = AID(i)
        SIDi = SID(i)
        BondOrdSIDjStart = BondOrdShift(SIDi)
        BondOrdSIDjStop = BondOrdSIDjStart + BondOrdStart(SIDi+1) - BondOrdStart(SIDi)
        BondOrdIndShift = BondOrdStart(SIDi) - BondOrdShift(SIDi)

        if (LoopMode == 0) then
            jstart = i+1
        else
            jstart = 0
        endif

        !loop over j
        do j = jstart, NAtom - 1

            !check to see if same atom
            if (i==j) cycle

            MIndj = MInd(j)

            if (LoopMode == 2 .or. LoopMode == -2) then
                !!check to see if we need to skip because of double counting
                if (MIndj == MIndi .and. j < i) cycle
            endif

            if (MolActive(MIndj) < 0 .and. MIndj /= ForceDoMInd .and. j /= ForceDoAtom) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)
            SIDj = SID(j)
            if (AIDi > AIDj) then
                AIDij = AIDi * (AIDi + 1) / 2 + AIDj
            else
                AIDij = AIDj * (AIDj + 1) / 2 + AIDi
            endif
            SameMol = (MIndi == MIndj)
            BondOrdij = BondOrdLimit + 1
            if (SameMol) then
                if (SIDj >= BondOrdSIDjStart .and. SIDj < BondOrdSIDjStop) then
                    BondOrdij = BondOrdData(SIDj + BondOrdIndShift)
                endif
            endif
            Bonded = (SameMol .and. BondOrdij==2)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            if (dijsq < CutSq(0)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if ((((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) .and. Bonded) then
                    !potential BondA_A-
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(0)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P0_ArgHistMin(ArgType)) * P0_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P0_ArgHist(ArgType,BinInd) = P0_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(1)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential LJGaussA_A
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(1)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P1_ArgHistMin(ArgType)) * P1_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P1_ArgHist(ArgType,BinInd) = P1_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(2)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    !potential LJGaussA_A-
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(2)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P2_ArgHistMin(ArgType)) * P2_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P2_ArgHist(ArgType,BinInd) = P2_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(3)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1)) .or. ((AIDj==1))) then
                    !potential LJGaussA_Na+
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(3)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P3_ArgHistMin(ArgType)) * P3_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P3_ArgHist(ArgType,BinInd) = P3_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(4)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential LJGaussA-_A-
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(4)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P4_ArgHistMin(ArgType)) * P4_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P4_ArgHist(ArgType,BinInd) = P4_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(5)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0)) .or. ((AIDj==0))) then
                    !potential LJGaussA-_Na+
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(5)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P5_ArgHistMin(ArgType)) * P5_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P5_ArgHist(ArgType,BinInd) = P5_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(6)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential LJGaussNa+_Na+
                ARGVAL = DIJ
                ARGTYPE = 0
                if (ArgVal*ArgVal > CutSq(6)) cycle
                if (10000 > 0) then
                    BinInd = int((ArgVal - P6_ArgHistMin(ArgType)) * P6_ArgHistiBinw(ArgType))
                    if (BinInd >= 0 .and. BinInd < 10000) P6_ArgHist(ArgType,BinInd) = P6_ArgHist(ArgType,BinInd) + Weight
                endif
            end if

            if (dijsq < CutSq(7)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential EW
                ARGVAL = DIJ
                ARGTYPE = AIDIJ
                if (10000 > 0) then
                    BinInd = int((ArgVal - P7_ArgHistMin(ArgType)) * P7_ArgHistiBinw(ArgType))
                    if (BinInd >= 0 .and. BinInd < 10000) P7_ArgHist(ArgType,BinInd) = P7_ArgHist(ArgType,BinInd) + Weight
                endif
            end if

            if (dijsq < CutSq(8)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential SmearCoulA_A
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P8_ArgHistMin(ArgType)) * P8_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P8_ArgHist(ArgType,BinInd) = P8_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(9)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    !potential SmearCoulA_A-
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P9_ArgHistMin(ArgType)) * P9_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P9_ArgHist(ArgType,BinInd) = P9_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(10)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1)) .or. ((AIDj==1))) then
                    !potential SmearCoulA_Na+
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P10_ArgHistMin(ArgType)) * P10_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P10_ArgHist(ArgType,BinInd) = P10_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(11)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential SmearCoulA-_A-
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P11_ArgHistMin(ArgType)) * P11_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P11_ArgHist(ArgType,BinInd) = P11_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(12)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0)) .or. ((AIDj==0))) then
                    !potential SmearCoulA-_Na+
                    ARGVAL = DIJ
                    ARGTYPE = AIDIJ
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P12_ArgHistMin(ArgType)) * P12_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P12_ArgHist(ArgType,BinInd) = P12_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(13)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                !potential SmearCoulNa+_Na+
                ARGVAL = DIJ
                ARGTYPE = AIDIJ
                if (10000 > 0) then
                    BinInd = int((ArgVal - P13_ArgHistMin(ArgType)) * P13_ArgHistiBinw(ArgType))
                    if (BinInd >= 0 .and. BinInd < 10000) P13_ArgHist(ArgType,BinInd) = P13_ArgHist(ArgType,BinInd) + Weight
                endif
            end if

        !end of loop j
        enddo

    !end of loop i
    enddo

end subroutine


subroutine calcargeval(CalcForce, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    logical, intent(in) :: CalcForce
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    integer :: i
    integer :: AIDi
    integer :: AIDj
    integer :: AIDij
    real(8) :: dijsq
    real(8) :: dij
    real(8) :: Scale
    real(8) :: ThisU
    real(8) :: ThisW
    integer :: ArgType
    real(8) :: ArgVal
    real(8) :: ThisHist
    real(8) :: val1
    real(8) :: idist2
    real(8) :: idist6
    real(8) :: idist12
    real(8) :: val2
    real(8) :: Sig
    real(8) :: Eps
    real(8) :: val3
    real(8) :: val4
    real(8) :: val5
    real(8) :: val6
    real(8) :: val7
    real(8) :: Chargei
    real(8) :: Chargej
    real(8) :: idist
    real(8) :: fac
    real(8) :: fac2
    real(8) :: iA
    real(8) :: tmp

    !compute initial quantities
    PEnergy = 0.d0
    Virial = 0.d0
    if (CalcDUParam) then
        DUParam = 0.d0
        DDUParam = 0.d0
    endif
    if (CalcDWParam) then
        DWParam = 0.d0
        DDWParam = 0.d0
    endif
    Terms = 0.d0
    Scale = 1.d0

    !potential BondA_A-
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P0_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P0_ArgHistMin(ArgType) + P0_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            val1 = DIJ - Param(2)
            THISU = Param(3) * val1*val1
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(0) = Terms(0) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = 2.d0 * Param(3) * val1*DIJ
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(3) = DUParam(3) + (ThisHist) * (val1 * val1)
                DUParam(2) = DUParam(2) + (ThisHist) * (-2.d0 * Param(3) * val1)
                DDUParam(4) = DDUParam(4) + (ThisHist) * (2.d0 * Param(3))
                DDUParam(6) = DDUParam(6) + (ThisHist) * (-2.d0 * val1)
            endif
            if (CALCDWPARAM) then
                DWParam(3) = DWParam(3) + (ThisHist) * (2.d0 * val1 * DIJ)
                DWParam(2) = DWParam(2) + (ThisHist) * (-2.d0 * Param(3) * DIJ)
                DDWParam(6) = DDWParam(6) + (ThisHist) * (-2.d0 * DIJ)
            endif
        enddo
    enddo

    !potential LJGaussA_A
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P1_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P1_ArgHistMin(ArgType) + P1_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(5)
            Eps = Param(4)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(8))
            val4 = val3*val3
            val5 = exp(-Param(7) * val4)
            THISU = Eps * val1 + Param(6) * val5 + P1_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(1) = Terms(1) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(6) * Param(7) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(4) = DUParam(4) + (ThisHist) * (val1 + P1_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(5) = DUParam(5) + (ThisHist) * (-Eps * val2 * val6 + P1_UShift(2))
                DDUParam(14) = DDUParam(14) + (ThisHist) * (Param(4) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P1_UShift(3))
                DDUParam(9) = DDUParam(9) + (ThisHist) * (-val2 * val6 + P1_UShift(4))
                DUParam(6) = DUParam(6) + (ThisHist) * (val5 + P1_UShift(5))
                val6 = -val4 * val5
                DUParam(7) = DUParam(7) + (ThisHist) * (Param(6) * val6 + P1_UShift(6))
                DDUParam(25) = DDUParam(25) + (ThisHist) * (val6 + P1_UShift(7))
                val6 = 2.d0 * Param(7) * val3 * val5
                DUParam(8) = DUParam(8) + (ThisHist) * (Param(6) * val6 + P1_UShift(8))
                DDUParam(30) = DDUParam(30) + (ThisHist) * (val6 + P1_UShift(9))
                DDUParam(26) = DDUParam(26) + (ThisHist) * (Param(6) * val4*val4 * val5 + P1_UShift(10))
                val6 = 2.d0*Param(6) * val5
                DDUParam(31) = DDUParam(31) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(7)) + P1_UShift(11))
                DDUParam(32) = DDUParam(32) + (ThisHist) * (val6 * Param(7) * (2.d0 * Param(7) * val4 - 1.d0) + P1_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(4) = DWParam(4) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(5) = DWParam(5) + (ThisHist) * (Eps *  val7)
                DDWParam(14) = DDWParam(14) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(9) = DDWParam(9) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(6) = DWParam(6) + (ThisHist) * (-val7 * Param(7))
                val6 = val4 * Param(7)
                DWParam(7) = DWParam(7) + (ThisHist) * (Param(6) * val7 * (val6 - 1.d0))
                DDWParam(26) = DDWParam(26) + (ThisHist) * (-Param(6) * val7 * val4 * (val6 - 2.d0))
                DDWParam(25) = DDWParam(25) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(8) = DWParam(8) + (ThisHist) * (-val7 * Param(6) * Param(7) * (2.d0 * Param(7) * val4 - 1.d0))
                DDWParam(27) = DDWParam(27) + (ThisHist) * (Param(6) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(32) = DDWParam(32) + (ThisHist) * (-2.d0 * val7 * Param(6) * Param(7)*Param(7) * val3 * (2.d0 * val6 - &
                  & 3.d0))
                DDWParam(22) = DDWParam(22) + (ThisHist) * (-val7 * Param(7) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential LJGaussA_A-
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P2_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P2_ArgHistMin(ArgType) + P2_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(10)
            Eps = Param(9)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(13))
            val4 = val3*val3
            val5 = exp(-Param(12) * val4)
            THISU = Eps * val1 + Param(11) * val5 + P2_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(2) = Terms(2) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(11) * Param(12) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(9) = DUParam(9) + (ThisHist) * (val1 + P2_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(10) = DUParam(10) + (ThisHist) * (-Eps * val2 * val6 + P2_UShift(2))
                DDUParam(39) = DDUParam(39) + (ThisHist) * (Param(9) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P2_UShift(3))
                DDUParam(34) = DDUParam(34) + (ThisHist) * (-val2 * val6 + P2_UShift(4))
                DUParam(11) = DUParam(11) + (ThisHist) * (val5 + P2_UShift(5))
                val6 = -val4 * val5
                DUParam(12) = DUParam(12) + (ThisHist) * (Param(11) * val6 + P2_UShift(6))
                DDUParam(50) = DDUParam(50) + (ThisHist) * (val6 + P2_UShift(7))
                val6 = 2.d0 * Param(12) * val3 * val5
                DUParam(13) = DUParam(13) + (ThisHist) * (Param(11) * val6 + P2_UShift(8))
                DDUParam(55) = DDUParam(55) + (ThisHist) * (val6 + P2_UShift(9))
                DDUParam(51) = DDUParam(51) + (ThisHist) * (Param(11) * val4*val4 * val5 + P2_UShift(10))
                val6 = 2.d0*Param(11) * val5
                DDUParam(56) = DDUParam(56) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(12)) + P2_UShift(11))
                DDUParam(57) = DDUParam(57) + (ThisHist) * (val6 * Param(12) * (2.d0 * Param(12) * val4 - 1.d0) + P2_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(9) = DWParam(9) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(10) = DWParam(10) + (ThisHist) * (Eps *  val7)
                DDWParam(39) = DDWParam(39) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(34) = DDWParam(34) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(11) = DWParam(11) + (ThisHist) * (-val7 * Param(12))
                val6 = val4 * Param(12)
                DWParam(12) = DWParam(12) + (ThisHist) * (Param(11) * val7 * (val6 - 1.d0))
                DDWParam(51) = DDWParam(51) + (ThisHist) * (-Param(11) * val7 * val4 * (val6 - 2.d0))
                DDWParam(50) = DDWParam(50) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(13) = DWParam(13) + (ThisHist) * (-val7 * Param(11) * Param(12) * (2.d0 * Param(12) * val4 - 1.d0))
                DDWParam(52) = DDWParam(52) + (ThisHist) * (Param(11) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(57) = DDWParam(57) + (ThisHist) * (-2.d0 * val7 * Param(11) * Param(12)*Param(12) * val3 * (2.d0 * &
                  & val6 - 3.d0))
                DDWParam(47) = DDWParam(47) + (ThisHist) * (-val7 * Param(12) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential LJGaussA_Na+
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P3_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P3_ArgHistMin(ArgType) + P3_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(15)
            Eps = Param(14)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(18))
            val4 = val3*val3
            val5 = exp(-Param(17) * val4)
            THISU = Eps * val1 + Param(16) * val5 + P3_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(3) = Terms(3) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(16) * Param(17) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(14) = DUParam(14) + (ThisHist) * (val1 + P3_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(15) = DUParam(15) + (ThisHist) * (-Eps * val2 * val6 + P3_UShift(2))
                DDUParam(64) = DDUParam(64) + (ThisHist) * (Param(14) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P3_UShift(3))
                DDUParam(59) = DDUParam(59) + (ThisHist) * (-val2 * val6 + P3_UShift(4))
                DUParam(16) = DUParam(16) + (ThisHist) * (val5 + P3_UShift(5))
                val6 = -val4 * val5
                DUParam(17) = DUParam(17) + (ThisHist) * (Param(16) * val6 + P3_UShift(6))
                DDUParam(75) = DDUParam(75) + (ThisHist) * (val6 + P3_UShift(7))
                val6 = 2.d0 * Param(17) * val3 * val5
                DUParam(18) = DUParam(18) + (ThisHist) * (Param(16) * val6 + P3_UShift(8))
                DDUParam(80) = DDUParam(80) + (ThisHist) * (val6 + P3_UShift(9))
                DDUParam(76) = DDUParam(76) + (ThisHist) * (Param(16) * val4*val4 * val5 + P3_UShift(10))
                val6 = 2.d0*Param(16) * val5
                DDUParam(81) = DDUParam(81) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(17)) + P3_UShift(11))
                DDUParam(82) = DDUParam(82) + (ThisHist) * (val6 * Param(17) * (2.d0 * Param(17) * val4 - 1.d0) + P3_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(14) = DWParam(14) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(15) = DWParam(15) + (ThisHist) * (Eps *  val7)
                DDWParam(64) = DDWParam(64) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(59) = DDWParam(59) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(16) = DWParam(16) + (ThisHist) * (-val7 * Param(17))
                val6 = val4 * Param(17)
                DWParam(17) = DWParam(17) + (ThisHist) * (Param(16) * val7 * (val6 - 1.d0))
                DDWParam(76) = DDWParam(76) + (ThisHist) * (-Param(16) * val7 * val4 * (val6 - 2.d0))
                DDWParam(75) = DDWParam(75) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(18) = DWParam(18) + (ThisHist) * (-val7 * Param(16) * Param(17) * (2.d0 * Param(17) * val4 - 1.d0))
                DDWParam(77) = DDWParam(77) + (ThisHist) * (Param(16) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(82) = DDWParam(82) + (ThisHist) * (-2.d0 * val7 * Param(16) * Param(17)*Param(17) * val3 * (2.d0 * &
                  & val6 - 3.d0))
                DDWParam(72) = DDWParam(72) + (ThisHist) * (-val7 * Param(17) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential LJGaussA-_A-
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P4_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P4_ArgHistMin(ArgType) + P4_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(20)
            Eps = Param(19)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(23))
            val4 = val3*val3
            val5 = exp(-Param(22) * val4)
            THISU = Eps * val1 + Param(21) * val5 + P4_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(4) = Terms(4) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(21) * Param(22) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(19) = DUParam(19) + (ThisHist) * (val1 + P4_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(20) = DUParam(20) + (ThisHist) * (-Eps * val2 * val6 + P4_UShift(2))
                DDUParam(89) = DDUParam(89) + (ThisHist) * (Param(19) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P4_UShift(3))
                DDUParam(84) = DDUParam(84) + (ThisHist) * (-val2 * val6 + P4_UShift(4))
                DUParam(21) = DUParam(21) + (ThisHist) * (val5 + P4_UShift(5))
                val6 = -val4 * val5
                DUParam(22) = DUParam(22) + (ThisHist) * (Param(21) * val6 + P4_UShift(6))
                DDUParam(100) = DDUParam(100) + (ThisHist) * (val6 + P4_UShift(7))
                val6 = 2.d0 * Param(22) * val3 * val5
                DUParam(23) = DUParam(23) + (ThisHist) * (Param(21) * val6 + P4_UShift(8))
                DDUParam(105) = DDUParam(105) + (ThisHist) * (val6 + P4_UShift(9))
                DDUParam(101) = DDUParam(101) + (ThisHist) * (Param(21) * val4*val4 * val5 + P4_UShift(10))
                val6 = 2.d0*Param(21) * val5
                DDUParam(106) = DDUParam(106) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(22)) + P4_UShift(11))
                DDUParam(107) = DDUParam(107) + (ThisHist) * (val6 * Param(22) * (2.d0 * Param(22) * val4 - 1.d0) + &
                  & P4_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(19) = DWParam(19) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(20) = DWParam(20) + (ThisHist) * (Eps *  val7)
                DDWParam(89) = DDWParam(89) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(84) = DDWParam(84) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(21) = DWParam(21) + (ThisHist) * (-val7 * Param(22))
                val6 = val4 * Param(22)
                DWParam(22) = DWParam(22) + (ThisHist) * (Param(21) * val7 * (val6 - 1.d0))
                DDWParam(101) = DDWParam(101) + (ThisHist) * (-Param(21) * val7 * val4 * (val6 - 2.d0))
                DDWParam(100) = DDWParam(100) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(23) = DWParam(23) + (ThisHist) * (-val7 * Param(21) * Param(22) * (2.d0 * Param(22) * val4 - 1.d0))
                DDWParam(102) = DDWParam(102) + (ThisHist) * (Param(21) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(107) = DDWParam(107) + (ThisHist) * (-2.d0 * val7 * Param(21) * Param(22)*Param(22) * val3 * (2.d0 * &
                  & val6 - 3.d0))
                DDWParam(97) = DDWParam(97) + (ThisHist) * (-val7 * Param(22) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential LJGaussA-_Na+
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P5_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P5_ArgHistMin(ArgType) + P5_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(25)
            Eps = Param(24)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(28))
            val4 = val3*val3
            val5 = exp(-Param(27) * val4)
            THISU = Eps * val1 + Param(26) * val5 + P5_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(5) = Terms(5) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(26) * Param(27) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(24) = DUParam(24) + (ThisHist) * (val1 + P5_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(25) = DUParam(25) + (ThisHist) * (-Eps * val2 * val6 + P5_UShift(2))
                DDUParam(114) = DDUParam(114) + (ThisHist) * (Param(24) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P5_UShift(3))
                DDUParam(109) = DDUParam(109) + (ThisHist) * (-val2 * val6 + P5_UShift(4))
                DUParam(26) = DUParam(26) + (ThisHist) * (val5 + P5_UShift(5))
                val6 = -val4 * val5
                DUParam(27) = DUParam(27) + (ThisHist) * (Param(26) * val6 + P5_UShift(6))
                DDUParam(125) = DDUParam(125) + (ThisHist) * (val6 + P5_UShift(7))
                val6 = 2.d0 * Param(27) * val3 * val5
                DUParam(28) = DUParam(28) + (ThisHist) * (Param(26) * val6 + P5_UShift(8))
                DDUParam(130) = DDUParam(130) + (ThisHist) * (val6 + P5_UShift(9))
                DDUParam(126) = DDUParam(126) + (ThisHist) * (Param(26) * val4*val4 * val5 + P5_UShift(10))
                val6 = 2.d0*Param(26) * val5
                DDUParam(131) = DDUParam(131) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(27)) + P5_UShift(11))
                DDUParam(132) = DDUParam(132) + (ThisHist) * (val6 * Param(27) * (2.d0 * Param(27) * val4 - 1.d0) + &
                  & P5_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(24) = DWParam(24) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(25) = DWParam(25) + (ThisHist) * (Eps *  val7)
                DDWParam(114) = DDWParam(114) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(109) = DDWParam(109) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(26) = DWParam(26) + (ThisHist) * (-val7 * Param(27))
                val6 = val4 * Param(27)
                DWParam(27) = DWParam(27) + (ThisHist) * (Param(26) * val7 * (val6 - 1.d0))
                DDWParam(126) = DDWParam(126) + (ThisHist) * (-Param(26) * val7 * val4 * (val6 - 2.d0))
                DDWParam(125) = DDWParam(125) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(28) = DWParam(28) + (ThisHist) * (-val7 * Param(26) * Param(27) * (2.d0 * Param(27) * val4 - 1.d0))
                DDWParam(127) = DDWParam(127) + (ThisHist) * (Param(26) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(132) = DDWParam(132) + (ThisHist) * (-2.d0 * val7 * Param(26) * Param(27)*Param(27) * val3 * (2.d0 * &
                  & val6 - 3.d0))
                DDWParam(122) = DDWParam(122) + (ThisHist) * (-val7 * Param(27) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential LJGaussNa+_Na+
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P6_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P6_ArgHistMin(ArgType) + P6_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(30)
            Eps = Param(29)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(33))
            val4 = val3*val3
            val5 = exp(-Param(32) * val4)
            THISU = Eps * val1 + Param(31) * val5 + P6_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(6) = Terms(6) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(31) * Param(32) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(29) = DUParam(29) + (ThisHist) * (val1 + P6_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(30) = DUParam(30) + (ThisHist) * (-Eps * val2 * val6 + P6_UShift(2))
                DDUParam(139) = DDUParam(139) + (ThisHist) * (Param(29) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P6_UShift(3))
                DDUParam(134) = DDUParam(134) + (ThisHist) * (-val2 * val6 + P6_UShift(4))
                DUParam(31) = DUParam(31) + (ThisHist) * (val5 + P6_UShift(5))
                val6 = -val4 * val5
                DUParam(32) = DUParam(32) + (ThisHist) * (Param(31) * val6 + P6_UShift(6))
                DDUParam(150) = DDUParam(150) + (ThisHist) * (val6 + P6_UShift(7))
                val6 = 2.d0 * Param(32) * val3 * val5
                DUParam(33) = DUParam(33) + (ThisHist) * (Param(31) * val6 + P6_UShift(8))
                DDUParam(155) = DDUParam(155) + (ThisHist) * (val6 + P6_UShift(9))
                DDUParam(151) = DDUParam(151) + (ThisHist) * (Param(31) * val4*val4 * val5 + P6_UShift(10))
                val6 = 2.d0*Param(31) * val5
                DDUParam(156) = DDUParam(156) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(32)) + P6_UShift(11))
                DDUParam(157) = DDUParam(157) + (ThisHist) * (val6 * Param(32) * (2.d0 * Param(32) * val4 - 1.d0) + &
                  & P6_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(29) = DWParam(29) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(30) = DWParam(30) + (ThisHist) * (Eps *  val7)
                DDWParam(139) = DDWParam(139) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(134) = DDWParam(134) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(31) = DWParam(31) + (ThisHist) * (-val7 * Param(32))
                val6 = val4 * Param(32)
                DWParam(32) = DWParam(32) + (ThisHist) * (Param(31) * val7 * (val6 - 1.d0))
                DDWParam(151) = DDWParam(151) + (ThisHist) * (-Param(31) * val7 * val4 * (val6 - 2.d0))
                DDWParam(150) = DDWParam(150) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(33) = DWParam(33) + (ThisHist) * (-val7 * Param(31) * Param(32) * (2.d0 * Param(32) * val4 - 1.d0))
                DDWParam(152) = DDWParam(152) + (ThisHist) * (Param(31) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(157) = DDWParam(157) + (ThisHist) * (-2.d0 * val7 * Param(31) * Param(32)*Param(32) * val3 * (2.d0 * &
                  & val6 - 3.d0))
                DDWParam(147) = DDWParam(147) + (ThisHist) * (-val7 * Param(32) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential EW
    do ArgType = 0, 3 - 1
        AIDIJ = ArgType
        call GetijFromPairInd(AIDIJ, AIDi, AIDj)
        Chargei = Param(AIDI)
        Chargej = Param(AIDJ)
        do i = 0, 10000 - 1
            ThisHist = P7_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P7_ArgHistMin(ArgType) + P7_ArgHistBinw(ArgType) * (0.5d0 + i)
            print *, "Evaluation of energy from pair distance distributions does not work for Ewald sum."
            stop
        enddo
    enddo

    !potential SmearCoulA_A
    do ArgType = 0, 3 - 1
        AIDIJ = ArgType
        call GetijFromPairInd(AIDIJ, AIDi, AIDj)
        Chargei = Param(AIDI)
        Chargej = Param(AIDJ)
        do i = 0, 10000 - 1
            ThisHist = P8_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P8_ArgHistMin(ArgType) + P8_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            idist = 1.d0 / DIJ
            fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(36)
            fac2 = fac*fac
            val1 = Chargei * Chargej
            call erf(fac*DIJ, tmp)
            val2 = tmp * idist
            val3 = Param(35) * val1
            THISU = val3 * (val2 - idist + P8_UShift(0) + P8_CoulShift(0))
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(8) = Terms(8) + ThisU * ThisHist
            iA = 1.d0 / Param(36)
            val5 = exp(-fac2*DIJSQ)*iA*iA
            val6 = exp(-fac2*DIJSQ)/Param(36)
            if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                THISW = val3 * ( val6 - val2 + idist)
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM) then
                DUParam(35) = DUParam(35) + (ThisHist) * (val1 * (val2 - idist + P8_UShift(0) + P8_CoulShift(0)))
                DUParam(36) = DUParam(36) + (ThisHist) * (val3*(-val5 + P8_UShift(1)  ))
                DDUParam(162) = DDUParam(162) + (ThisHist) * (val3 * ( val6*2.d0/Param(36)**2.d0 - &
                  & val6*2.d0/Param(36)**2.d0*fac2*DIJSQ  + P8_UShift(2)  ))
                DDUParam(160) = DDUParam(160) + (ThisHist) * (val1*(-val5 + P8_UShift(1) ))
                if (AIDI==AIDJ) then
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (2.d0 * Param(35) * Chargej * (val2 + P8_UShift(0) - idist + &
                      & P8_CoulShift(0)))
                    DDUParam(251 + AIDI + 0) = DDUParam(251 + AIDI + 0) + (ThisHist) * (2.d0 * Param(35) * Chargej * (-val5 + &
                      & P8_UShift(1)))
                    DDUParam(249 + AIDI + 0) = DDUParam(249 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val2 + P8_UShift(0) - &
                      & idist + P8_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(35) * (val2 + P8_UShift(0) - &
                      & idist + P8_CoulShift(0)))
                else
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (Param(35) * Chargej * (val2 + P8_UShift(0) -idist + &
                      & P8_CoulShift(0)))
                    DUParam(AIDJ) = DUParam(AIDJ) + (ThisHist) * (Param(35) * Chargei * (val2 + P8_UShift(0) -idist + &
                      & P8_CoulShift(0)))
                    DDUParam(251 + AIDI + 0) = DDUParam(251 + AIDI + 0) + (ThisHist) * (Param(35) * Chargej * (-val5 + &
                      & P8_UShift(1)))
                    DDUParam(251 + AIDJ + 0) = DDUParam(251 + AIDJ + 0) + (ThisHist) * (Param(35) * Chargei * (-val5 + &
                      & P8_UShift(1)))
                    DDUParam(249 + AIDI + 0) = DDUParam(249 + AIDI + 0) + (ThisHist) * (Chargej * (val2 + P8_UShift(0) - idist  &
                      & + P8_CoulShift(0)))
                    DDUParam(249 + AIDJ + 0) = DDUParam(249 + AIDJ + 0) + (ThisHist) * (Chargei * (val2 + P8_UShift(0) - idist  &
                      & + P8_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(35) * (val2 + P8_UShift(0) - idist  &
                      & + P8_CoulShift(0)))
                endif
            endif
            if (CALCDWPARAM) then

                DWParam(35) = DWParam(35) + (ThisHist) * (val1 * (val6 - val2 + idist))
                val4 = 2.d0 * val5 * fac2 * DIJSQ
                DWParam(36) = DWParam(36) + (ThisHist) * (val3*val4)

                DDWParam(162) = DDWParam(162) + (ThisHist) * (2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4)
                DDWParam(160) = DDWParam(160) + (ThisHist) * (val1 * val4)

                if (AIDI==AIDJ) then
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (2.d0 * Chargej * Param(35) * (val6 - val2 + idist))
                    DDWParam(251 + AIDI + 0) = DDWParam(251 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * Param(35) * val4)
                    DDWParam(249 + AIDI + 0) = DDWParam(249 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(35) * (val6 - val2 + idist))
                else
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (Param(35) * Chargej * (val6 - val2 + idist))
                    DWParam(AIDJ) = DWParam(AIDJ) + (ThisHist) * (Param(35) * Chargei * (val6 - val2 + idist))

                    DDWParam(251 + AIDI + 0) = DDWParam(251 + AIDI + 0) + (ThisHist) * (Param(35) * val4 * Chargej)
                    DDWParam(251 + AIDJ + 0) = DDWParam(251 + AIDJ + 0) + (ThisHist) * (Param(35) * val4 * Chargei)
                    DDWParam(249 + AIDI + 0) = DDWParam(249 + AIDI + 0) + (ThisHist) * (Chargej * (val6 - val2 + idist))
                    DDWParam(249 + AIDJ + 0) = DDWParam(249 + AIDJ + 0) + (ThisHist) * (Chargei * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(35) * (val6 - val2 + idist))
                endif
            endif
        enddo
    enddo

    !potential SmearCoulA_A-
    do ArgType = 0, 3 - 1
        AIDIJ = ArgType
        call GetijFromPairInd(AIDIJ, AIDi, AIDj)
        Chargei = Param(AIDI)
        Chargej = Param(AIDJ)
        do i = 0, 10000 - 1
            ThisHist = P9_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P9_ArgHistMin(ArgType) + P9_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            idist = 1.d0 / DIJ
            fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(38)
            fac2 = fac*fac
            val1 = Chargei * Chargej
            call erf(fac*DIJ, tmp)
            val2 = tmp * idist
            val3 = Param(37) * val1
            THISU = val3 * (val2 - idist + P9_UShift(0) + P9_CoulShift(0))
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(9) = Terms(9) + ThisU * ThisHist
            iA = 1.d0 / Param(38)
            val5 = exp(-fac2*DIJSQ)*iA*iA
            val6 = exp(-fac2*DIJSQ)/Param(38)
            if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                THISW = val3 * ( val6 - val2 + idist)
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM) then
                DUParam(37) = DUParam(37) + (ThisHist) * (val1 * (val2 - idist + P9_UShift(0) + P9_CoulShift(0)))
                DUParam(38) = DUParam(38) + (ThisHist) * (val3*(-val5 + P9_UShift(1)  ))
                DDUParam(166) = DDUParam(166) + (ThisHist) * (val3 * ( val6*2.d0/Param(38)**2.d0 - &
                  & val6*2.d0/Param(38)**2.d0*fac2*DIJSQ  + P9_UShift(2)  ))
                DDUParam(164) = DDUParam(164) + (ThisHist) * (val1*(-val5 + P9_UShift(1) ))
                if (AIDI==AIDJ) then
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (2.d0 * Param(37) * Chargej * (val2 + P9_UShift(0) - idist + &
                      & P9_CoulShift(0)))
                    DDUParam(255 + AIDI + 0) = DDUParam(255 + AIDI + 0) + (ThisHist) * (2.d0 * Param(37) * Chargej * (-val5 + &
                      & P9_UShift(1)))
                    DDUParam(253 + AIDI + 0) = DDUParam(253 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val2 + P9_UShift(0) - &
                      & idist + P9_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(37) * (val2 + P9_UShift(0) - &
                      & idist + P9_CoulShift(0)))
                else
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (Param(37) * Chargej * (val2 + P9_UShift(0) -idist + &
                      & P9_CoulShift(0)))
                    DUParam(AIDJ) = DUParam(AIDJ) + (ThisHist) * (Param(37) * Chargei * (val2 + P9_UShift(0) -idist + &
                      & P9_CoulShift(0)))
                    DDUParam(255 + AIDI + 0) = DDUParam(255 + AIDI + 0) + (ThisHist) * (Param(37) * Chargej * (-val5 + &
                      & P9_UShift(1)))
                    DDUParam(255 + AIDJ + 0) = DDUParam(255 + AIDJ + 0) + (ThisHist) * (Param(37) * Chargei * (-val5 + &
                      & P9_UShift(1)))
                    DDUParam(253 + AIDI + 0) = DDUParam(253 + AIDI + 0) + (ThisHist) * (Chargej * (val2 + P9_UShift(0) - idist  &
                      & + P9_CoulShift(0)))
                    DDUParam(253 + AIDJ + 0) = DDUParam(253 + AIDJ + 0) + (ThisHist) * (Chargei * (val2 + P9_UShift(0) - idist  &
                      & + P9_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(37) * (val2 + P9_UShift(0) - idist  &
                      & + P9_CoulShift(0)))
                endif
            endif
            if (CALCDWPARAM) then

                DWParam(37) = DWParam(37) + (ThisHist) * (val1 * (val6 - val2 + idist))
                val4 = 2.d0 * val5 * fac2 * DIJSQ
                DWParam(38) = DWParam(38) + (ThisHist) * (val3*val4)

                DDWParam(166) = DDWParam(166) + (ThisHist) * (2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4)
                DDWParam(164) = DDWParam(164) + (ThisHist) * (val1 * val4)

                if (AIDI==AIDJ) then
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (2.d0 * Chargej * Param(37) * (val6 - val2 + idist))
                    DDWParam(255 + AIDI + 0) = DDWParam(255 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * Param(37) * val4)
                    DDWParam(253 + AIDI + 0) = DDWParam(253 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(37) * (val6 - val2 + idist))
                else
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (Param(37) * Chargej * (val6 - val2 + idist))
                    DWParam(AIDJ) = DWParam(AIDJ) + (ThisHist) * (Param(37) * Chargei * (val6 - val2 + idist))

                    DDWParam(255 + AIDI + 0) = DDWParam(255 + AIDI + 0) + (ThisHist) * (Param(37) * val4 * Chargej)
                    DDWParam(255 + AIDJ + 0) = DDWParam(255 + AIDJ + 0) + (ThisHist) * (Param(37) * val4 * Chargei)
                    DDWParam(253 + AIDI + 0) = DDWParam(253 + AIDI + 0) + (ThisHist) * (Chargej * (val6 - val2 + idist))
                    DDWParam(253 + AIDJ + 0) = DDWParam(253 + AIDJ + 0) + (ThisHist) * (Chargei * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(37) * (val6 - val2 + idist))
                endif
            endif
        enddo
    enddo

    !potential SmearCoulA_Na+
    do ArgType = 0, 3 - 1
        AIDIJ = ArgType
        call GetijFromPairInd(AIDIJ, AIDi, AIDj)
        Chargei = Param(AIDI)
        Chargej = Param(AIDJ)
        do i = 0, 10000 - 1
            ThisHist = P10_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P10_ArgHistMin(ArgType) + P10_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            idist = 1.d0 / DIJ
            fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(40)
            fac2 = fac*fac
            val1 = Chargei * Chargej
            call erf(fac*DIJ, tmp)
            val2 = tmp * idist
            val3 = Param(39) * val1
            THISU = val3 * (val2 - idist + P10_UShift(0) + P10_CoulShift(0))
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(10) = Terms(10) + ThisU * ThisHist
            iA = 1.d0 / Param(40)
            val5 = exp(-fac2*DIJSQ)*iA*iA
            val6 = exp(-fac2*DIJSQ)/Param(40)
            if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                THISW = val3 * ( val6 - val2 + idist)
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM) then
                DUParam(39) = DUParam(39) + (ThisHist) * (val1 * (val2 - idist + P10_UShift(0) + P10_CoulShift(0)))
                DUParam(40) = DUParam(40) + (ThisHist) * (val3*(-val5 + P10_UShift(1)  ))
                DDUParam(170) = DDUParam(170) + (ThisHist) * (val3 * ( val6*2.d0/Param(40)**2.d0 - &
                  & val6*2.d0/Param(40)**2.d0*fac2*DIJSQ  + P10_UShift(2)  ))
                DDUParam(168) = DDUParam(168) + (ThisHist) * (val1*(-val5 + P10_UShift(1) ))
                if (AIDI==AIDJ) then
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (2.d0 * Param(39) * Chargej * (val2 + P10_UShift(0) - idist + &
                      & P10_CoulShift(0)))
                    DDUParam(259 + AIDI + 0) = DDUParam(259 + AIDI + 0) + (ThisHist) * (2.d0 * Param(39) * Chargej * (-val5 + &
                      & P10_UShift(1)))
                    DDUParam(257 + AIDI + 0) = DDUParam(257 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val2 + P10_UShift(0)  &
                      & - idist + P10_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(39) * (val2 + P10_UShift(0)  &
                      & - idist + P10_CoulShift(0)))
                else
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (Param(39) * Chargej * (val2 + P10_UShift(0) -idist + &
                      & P10_CoulShift(0)))
                    DUParam(AIDJ) = DUParam(AIDJ) + (ThisHist) * (Param(39) * Chargei * (val2 + P10_UShift(0) -idist + &
                      & P10_CoulShift(0)))
                    DDUParam(259 + AIDI + 0) = DDUParam(259 + AIDI + 0) + (ThisHist) * (Param(39) * Chargej * (-val5 + &
                      & P10_UShift(1)))
                    DDUParam(259 + AIDJ + 0) = DDUParam(259 + AIDJ + 0) + (ThisHist) * (Param(39) * Chargei * (-val5 + &
                      & P10_UShift(1)))
                    DDUParam(257 + AIDI + 0) = DDUParam(257 + AIDI + 0) + (ThisHist) * (Chargej * (val2 + P10_UShift(0) - idist &
                      & + P10_CoulShift(0)))
                    DDUParam(257 + AIDJ + 0) = DDUParam(257 + AIDJ + 0) + (ThisHist) * (Chargei * (val2 + P10_UShift(0) - idist &
                      & + P10_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(39) * (val2 + P10_UShift(0) - idist &
                      & + P10_CoulShift(0)))
                endif
            endif
            if (CALCDWPARAM) then

                DWParam(39) = DWParam(39) + (ThisHist) * (val1 * (val6 - val2 + idist))
                val4 = 2.d0 * val5 * fac2 * DIJSQ
                DWParam(40) = DWParam(40) + (ThisHist) * (val3*val4)

                DDWParam(170) = DDWParam(170) + (ThisHist) * (2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4)
                DDWParam(168) = DDWParam(168) + (ThisHist) * (val1 * val4)

                if (AIDI==AIDJ) then
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (2.d0 * Chargej * Param(39) * (val6 - val2 + idist))
                    DDWParam(259 + AIDI + 0) = DDWParam(259 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * Param(39) * val4)
                    DDWParam(257 + AIDI + 0) = DDWParam(257 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(39) * (val6 - val2 + idist))
                else
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (Param(39) * Chargej * (val6 - val2 + idist))
                    DWParam(AIDJ) = DWParam(AIDJ) + (ThisHist) * (Param(39) * Chargei * (val6 - val2 + idist))

                    DDWParam(259 + AIDI + 0) = DDWParam(259 + AIDI + 0) + (ThisHist) * (Param(39) * val4 * Chargej)
                    DDWParam(259 + AIDJ + 0) = DDWParam(259 + AIDJ + 0) + (ThisHist) * (Param(39) * val4 * Chargei)
                    DDWParam(257 + AIDI + 0) = DDWParam(257 + AIDI + 0) + (ThisHist) * (Chargej * (val6 - val2 + idist))
                    DDWParam(257 + AIDJ + 0) = DDWParam(257 + AIDJ + 0) + (ThisHist) * (Chargei * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(39) * (val6 - val2 + idist))
                endif
            endif
        enddo
    enddo

    !potential SmearCoulA-_A-
    do ArgType = 0, 3 - 1
        AIDIJ = ArgType
        call GetijFromPairInd(AIDIJ, AIDi, AIDj)
        Chargei = Param(AIDI)
        Chargej = Param(AIDJ)
        do i = 0, 10000 - 1
            ThisHist = P11_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P11_ArgHistMin(ArgType) + P11_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            idist = 1.d0 / DIJ
            fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(42)
            fac2 = fac*fac
            val1 = Chargei * Chargej
            call erf(fac*DIJ, tmp)
            val2 = tmp * idist
            val3 = Param(41) * val1
            THISU = val3 * (val2 - idist + P11_UShift(0) + P11_CoulShift(0))
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(11) = Terms(11) + ThisU * ThisHist
            iA = 1.d0 / Param(42)
            val5 = exp(-fac2*DIJSQ)*iA*iA
            val6 = exp(-fac2*DIJSQ)/Param(42)
            if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                THISW = val3 * ( val6 - val2 + idist)
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM) then
                DUParam(41) = DUParam(41) + (ThisHist) * (val1 * (val2 - idist + P11_UShift(0) + P11_CoulShift(0)))
                DUParam(42) = DUParam(42) + (ThisHist) * (val3*(-val5 + P11_UShift(1)  ))
                DDUParam(174) = DDUParam(174) + (ThisHist) * (val3 * ( val6*2.d0/Param(42)**2.d0 - &
                  & val6*2.d0/Param(42)**2.d0*fac2*DIJSQ  + P11_UShift(2)  ))
                DDUParam(172) = DDUParam(172) + (ThisHist) * (val1*(-val5 + P11_UShift(1) ))
                if (AIDI==AIDJ) then
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (2.d0 * Param(41) * Chargej * (val2 + P11_UShift(0) - idist + &
                      & P11_CoulShift(0)))
                    DDUParam(263 + AIDI + 0) = DDUParam(263 + AIDI + 0) + (ThisHist) * (2.d0 * Param(41) * Chargej * (-val5 + &
                      & P11_UShift(1)))
                    DDUParam(261 + AIDI + 0) = DDUParam(261 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val2 + P11_UShift(0)  &
                      & - idist + P11_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(41) * (val2 + P11_UShift(0)  &
                      & - idist + P11_CoulShift(0)))
                else
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (Param(41) * Chargej * (val2 + P11_UShift(0) -idist + &
                      & P11_CoulShift(0)))
                    DUParam(AIDJ) = DUParam(AIDJ) + (ThisHist) * (Param(41) * Chargei * (val2 + P11_UShift(0) -idist + &
                      & P11_CoulShift(0)))
                    DDUParam(263 + AIDI + 0) = DDUParam(263 + AIDI + 0) + (ThisHist) * (Param(41) * Chargej * (-val5 + &
                      & P11_UShift(1)))
                    DDUParam(263 + AIDJ + 0) = DDUParam(263 + AIDJ + 0) + (ThisHist) * (Param(41) * Chargei * (-val5 + &
                      & P11_UShift(1)))
                    DDUParam(261 + AIDI + 0) = DDUParam(261 + AIDI + 0) + (ThisHist) * (Chargej * (val2 + P11_UShift(0) - idist &
                      & + P11_CoulShift(0)))
                    DDUParam(261 + AIDJ + 0) = DDUParam(261 + AIDJ + 0) + (ThisHist) * (Chargei * (val2 + P11_UShift(0) - idist &
                      & + P11_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(41) * (val2 + P11_UShift(0) - idist &
                      & + P11_CoulShift(0)))
                endif
            endif
            if (CALCDWPARAM) then

                DWParam(41) = DWParam(41) + (ThisHist) * (val1 * (val6 - val2 + idist))
                val4 = 2.d0 * val5 * fac2 * DIJSQ
                DWParam(42) = DWParam(42) + (ThisHist) * (val3*val4)

                DDWParam(174) = DDWParam(174) + (ThisHist) * (2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4)
                DDWParam(172) = DDWParam(172) + (ThisHist) * (val1 * val4)

                if (AIDI==AIDJ) then
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (2.d0 * Chargej * Param(41) * (val6 - val2 + idist))
                    DDWParam(263 + AIDI + 0) = DDWParam(263 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * Param(41) * val4)
                    DDWParam(261 + AIDI + 0) = DDWParam(261 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(41) * (val6 - val2 + idist))
                else
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (Param(41) * Chargej * (val6 - val2 + idist))
                    DWParam(AIDJ) = DWParam(AIDJ) + (ThisHist) * (Param(41) * Chargei * (val6 - val2 + idist))

                    DDWParam(263 + AIDI + 0) = DDWParam(263 + AIDI + 0) + (ThisHist) * (Param(41) * val4 * Chargej)
                    DDWParam(263 + AIDJ + 0) = DDWParam(263 + AIDJ + 0) + (ThisHist) * (Param(41) * val4 * Chargei)
                    DDWParam(261 + AIDI + 0) = DDWParam(261 + AIDI + 0) + (ThisHist) * (Chargej * (val6 - val2 + idist))
                    DDWParam(261 + AIDJ + 0) = DDWParam(261 + AIDJ + 0) + (ThisHist) * (Chargei * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(41) * (val6 - val2 + idist))
                endif
            endif
        enddo
    enddo

    !potential SmearCoulA-_Na+
    do ArgType = 0, 3 - 1
        AIDIJ = ArgType
        call GetijFromPairInd(AIDIJ, AIDi, AIDj)
        Chargei = Param(AIDI)
        Chargej = Param(AIDJ)
        do i = 0, 10000 - 1
            ThisHist = P12_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P12_ArgHistMin(ArgType) + P12_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            idist = 1.d0 / DIJ
            fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(44)
            fac2 = fac*fac
            val1 = Chargei * Chargej
            call erf(fac*DIJ, tmp)
            val2 = tmp * idist
            val3 = Param(43) * val1
            THISU = val3 * (val2 - idist + P12_UShift(0) + P12_CoulShift(0))
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(12) = Terms(12) + ThisU * ThisHist
            iA = 1.d0 / Param(44)
            val5 = exp(-fac2*DIJSQ)*iA*iA
            val6 = exp(-fac2*DIJSQ)/Param(44)
            if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                THISW = val3 * ( val6 - val2 + idist)
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM) then
                DUParam(43) = DUParam(43) + (ThisHist) * (val1 * (val2 - idist + P12_UShift(0) + P12_CoulShift(0)))
                DUParam(44) = DUParam(44) + (ThisHist) * (val3*(-val5 + P12_UShift(1)  ))
                DDUParam(178) = DDUParam(178) + (ThisHist) * (val3 * ( val6*2.d0/Param(44)**2.d0 - &
                  & val6*2.d0/Param(44)**2.d0*fac2*DIJSQ  + P12_UShift(2)  ))
                DDUParam(176) = DDUParam(176) + (ThisHist) * (val1*(-val5 + P12_UShift(1) ))
                if (AIDI==AIDJ) then
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (2.d0 * Param(43) * Chargej * (val2 + P12_UShift(0) - idist + &
                      & P12_CoulShift(0)))
                    DDUParam(267 + AIDI + 0) = DDUParam(267 + AIDI + 0) + (ThisHist) * (2.d0 * Param(43) * Chargej * (-val5 + &
                      & P12_UShift(1)))
                    DDUParam(265 + AIDI + 0) = DDUParam(265 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val2 + P12_UShift(0)  &
                      & - idist + P12_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(43) * (val2 + P12_UShift(0)  &
                      & - idist + P12_CoulShift(0)))
                else
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (Param(43) * Chargej * (val2 + P12_UShift(0) -idist + &
                      & P12_CoulShift(0)))
                    DUParam(AIDJ) = DUParam(AIDJ) + (ThisHist) * (Param(43) * Chargei * (val2 + P12_UShift(0) -idist + &
                      & P12_CoulShift(0)))
                    DDUParam(267 + AIDI + 0) = DDUParam(267 + AIDI + 0) + (ThisHist) * (Param(43) * Chargej * (-val5 + &
                      & P12_UShift(1)))
                    DDUParam(267 + AIDJ + 0) = DDUParam(267 + AIDJ + 0) + (ThisHist) * (Param(43) * Chargei * (-val5 + &
                      & P12_UShift(1)))
                    DDUParam(265 + AIDI + 0) = DDUParam(265 + AIDI + 0) + (ThisHist) * (Chargej * (val2 + P12_UShift(0) - idist &
                      & + P12_CoulShift(0)))
                    DDUParam(265 + AIDJ + 0) = DDUParam(265 + AIDJ + 0) + (ThisHist) * (Chargei * (val2 + P12_UShift(0) - idist &
                      & + P12_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(43) * (val2 + P12_UShift(0) - idist &
                      & + P12_CoulShift(0)))
                endif
            endif
            if (CALCDWPARAM) then

                DWParam(43) = DWParam(43) + (ThisHist) * (val1 * (val6 - val2 + idist))
                val4 = 2.d0 * val5 * fac2 * DIJSQ
                DWParam(44) = DWParam(44) + (ThisHist) * (val3*val4)

                DDWParam(178) = DDWParam(178) + (ThisHist) * (2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4)
                DDWParam(176) = DDWParam(176) + (ThisHist) * (val1 * val4)

                if (AIDI==AIDJ) then
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (2.d0 * Chargej * Param(43) * (val6 - val2 + idist))
                    DDWParam(267 + AIDI + 0) = DDWParam(267 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * Param(43) * val4)
                    DDWParam(265 + AIDI + 0) = DDWParam(265 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(43) * (val6 - val2 + idist))
                else
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (Param(43) * Chargej * (val6 - val2 + idist))
                    DWParam(AIDJ) = DWParam(AIDJ) + (ThisHist) * (Param(43) * Chargei * (val6 - val2 + idist))

                    DDWParam(267 + AIDI + 0) = DDWParam(267 + AIDI + 0) + (ThisHist) * (Param(43) * val4 * Chargej)
                    DDWParam(267 + AIDJ + 0) = DDWParam(267 + AIDJ + 0) + (ThisHist) * (Param(43) * val4 * Chargei)
                    DDWParam(265 + AIDI + 0) = DDWParam(265 + AIDI + 0) + (ThisHist) * (Chargej * (val6 - val2 + idist))
                    DDWParam(265 + AIDJ + 0) = DDWParam(265 + AIDJ + 0) + (ThisHist) * (Chargei * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(43) * (val6 - val2 + idist))
                endif
            endif
        enddo
    enddo

    !potential SmearCoulNa+_Na+
    do ArgType = 0, 3 - 1
        AIDIJ = ArgType
        call GetijFromPairInd(AIDIJ, AIDi, AIDj)
        Chargei = Param(AIDI)
        Chargej = Param(AIDJ)
        do i = 0, 10000 - 1
            ThisHist = P13_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P13_ArgHistMin(ArgType) + P13_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            idist = 1.d0 / DIJ
            fac  = sqrt(4.d0*atan(1.0_8))/2.d0/Param(46)
            fac2 = fac*fac
            val1 = Chargei * Chargej
            call erf(fac*DIJ, tmp)
            val2 = tmp * idist
            val3 = Param(45) * val1
            THISU = val3 * (val2 - idist + P13_UShift(0) + P13_CoulShift(0))
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(13) = Terms(13) + ThisU * ThisHist
            iA = 1.d0 / Param(46)
            val5 = exp(-fac2*DIJSQ)*iA*iA
            val6 = exp(-fac2*DIJSQ)/Param(46)
            if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                THISW = val3 * ( val6 - val2 + idist)
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM) then
                DUParam(45) = DUParam(45) + (ThisHist) * (val1 * (val2 - idist + P13_UShift(0) + P13_CoulShift(0)))
                DUParam(46) = DUParam(46) + (ThisHist) * (val3*(-val5 + P13_UShift(1)  ))
                DDUParam(182) = DDUParam(182) + (ThisHist) * (val3 * ( val6*2.d0/Param(46)**2.d0 - &
                  & val6*2.d0/Param(46)**2.d0*fac2*DIJSQ  + P13_UShift(2)  ))
                DDUParam(180) = DDUParam(180) + (ThisHist) * (val1*(-val5 + P13_UShift(1) ))
                if (AIDI==AIDJ) then
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (2.d0 * Param(45) * Chargej * (val2 + P13_UShift(0) - idist + &
                      & P13_CoulShift(0)))
                    DDUParam(271 + AIDI + 0) = DDUParam(271 + AIDI + 0) + (ThisHist) * (2.d0 * Param(45) * Chargej * (-val5 + &
                      & P13_UShift(1)))
                    DDUParam(269 + AIDI + 0) = DDUParam(269 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val2 + P13_UShift(0)  &
                      & - idist + P13_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(45) * (val2 + P13_UShift(0)  &
                      & - idist + P13_CoulShift(0)))
                else
                    DUParam(AIDI) = DUParam(AIDI) + (ThisHist) * (Param(45) * Chargej * (val2 + P13_UShift(0) -idist + &
                      & P13_CoulShift(0)))
                    DUParam(AIDJ) = DUParam(AIDJ) + (ThisHist) * (Param(45) * Chargei * (val2 + P13_UShift(0) -idist + &
                      & P13_CoulShift(0)))
                    DDUParam(271 + AIDI + 0) = DDUParam(271 + AIDI + 0) + (ThisHist) * (Param(45) * Chargej * (-val5 + &
                      & P13_UShift(1)))
                    DDUParam(271 + AIDJ + 0) = DDUParam(271 + AIDJ + 0) + (ThisHist) * (Param(45) * Chargei * (-val5 + &
                      & P13_UShift(1)))
                    DDUParam(269 + AIDI + 0) = DDUParam(269 + AIDI + 0) + (ThisHist) * (Chargej * (val2 + P13_UShift(0) - idist &
                      & + P13_CoulShift(0)))
                    DDUParam(269 + AIDJ + 0) = DDUParam(269 + AIDJ + 0) + (ThisHist) * (Chargei * (val2 + P13_UShift(0) - idist &
                      & + P13_CoulShift(0)))
                    DDUParam(AIDI + 2*AIDJ) = DDUParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(45) * (val2 + P13_UShift(0) - idist &
                      & + P13_CoulShift(0)))
                endif
            endif
            if (CALCDWPARAM) then

                DWParam(45) = DWParam(45) + (ThisHist) * (val1 * (val6 - val2 + idist))
                val4 = 2.d0 * val5 * fac2 * DIJSQ
                DWParam(46) = DWParam(46) + (ThisHist) * (val3*val4)

                DDWParam(182) = DDWParam(182) + (ThisHist) * (2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4)
                DDWParam(180) = DDWParam(180) + (ThisHist) * (val1 * val4)

                if (AIDI==AIDJ) then
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (2.d0 * Chargej * Param(45) * (val6 - val2 + idist))
                    DDWParam(271 + AIDI + 0) = DDWParam(271 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * Param(45) * val4)
                    DDWParam(269 + AIDI + 0) = DDWParam(269 + AIDI + 0) + (ThisHist) * (2.d0 * Chargej * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (2.d0 * Param(45) * (val6 - val2 + idist))
                else
                    DWParam(AIDI) = DWParam(AIDI) + (ThisHist) * (Param(45) * Chargej * (val6 - val2 + idist))
                    DWParam(AIDJ) = DWParam(AIDJ) + (ThisHist) * (Param(45) * Chargei * (val6 - val2 + idist))

                    DDWParam(271 + AIDI + 0) = DDWParam(271 + AIDI + 0) + (ThisHist) * (Param(45) * val4 * Chargej)
                    DDWParam(271 + AIDJ + 0) = DDWParam(271 + AIDJ + 0) + (ThisHist) * (Param(45) * val4 * Chargei)
                    DDWParam(269 + AIDI + 0) = DDWParam(269 + AIDI + 0) + (ThisHist) * (Chargej * (val6 - val2 + idist))
                    DDWParam(269 + AIDJ + 0) = DDWParam(269 + AIDJ + 0) + (ThisHist) * (Chargei * (val6 - val2 + idist))
                    DDWParam(AIDI + 2*AIDJ) = DDWParam(AIDI + 2*AIDJ) + (ThisHist) * (Param(45) * (val6 - val2 + idist))
                endif
            endif
        enddo
    enddo

end subroutine


subroutine calcmeasures(MeasureAll, StepNum, CycleNum, Weight, ConservesMomentum)
    implicit none
    logical, intent(in) :: MeasureAll
    integer, intent(in) :: StepNum
    integer, intent(in) :: CycleNum
    real(8), intent(in) :: Weight
    logical, intent(in) :: ConservesMomentum
    integer :: i
    integer :: j
    logical, dimension(0:NMeasure-1) :: UseMeasure
    integer :: v1
    integer :: v2
    real(8) :: Val0

    !measure KEnergy
    if (M0_Active) then
        if (MeasureAll) then
            UseMeasure(0) = .true.
            M0_Val = 0.d0
        elseif (M0_StepFreq > 0 .and. mod(StepNum, M0_StepFreq)==0) then
            UseMeasure(0) = .true.
            M0_Val = 0.d0
        elseif (M0_CycleFreq > 0 .and. mod(CycleNum, M0_CycleFreq)==0) then
            UseMeasure(0) = .true.
            M0_Val = 0.d0
        else
            UseMeasure(0) = .false.
        endif
    else
        UseMeasure(0) = .false.
    endif
    !measure PEnergy
    if (M1_Active) then
        if (MeasureAll) then
            UseMeasure(1) = .true.
            M1_Val = 0.d0
        elseif (M1_StepFreq > 0 .and. mod(StepNum, M1_StepFreq)==0) then
            UseMeasure(1) = .true.
            M1_Val = 0.d0
        elseif (M1_CycleFreq > 0 .and. mod(CycleNum, M1_CycleFreq)==0) then
            UseMeasure(1) = .true.
            M1_Val = 0.d0
        else
            UseMeasure(1) = .false.
        endif
    else
        UseMeasure(1) = .false.
    endif
    !measure TEnergy
    if (M2_Active) then
        if (MeasureAll) then
            UseMeasure(2) = .true.
            M2_Val = 0.d0
        elseif (M2_StepFreq > 0 .and. mod(StepNum, M2_StepFreq)==0) then
            UseMeasure(2) = .true.
            M2_Val = 0.d0
        elseif (M2_CycleFreq > 0 .and. mod(CycleNum, M2_CycleFreq)==0) then
            UseMeasure(2) = .true.
            M2_Val = 0.d0
        else
            UseMeasure(2) = .false.
        endif
    else
        UseMeasure(2) = .false.
    endif
    !measure KTemp
    if (M3_Active) then
        if (MeasureAll) then
            UseMeasure(3) = .true.
            M3_Val = 0.d0
        elseif (M3_StepFreq > 0 .and. mod(StepNum, M3_StepFreq)==0) then
            UseMeasure(3) = .true.
            M3_Val = 0.d0
        elseif (M3_CycleFreq > 0 .and. mod(CycleNum, M3_CycleFreq)==0) then
            UseMeasure(3) = .true.
            M3_Val = 0.d0
        else
            UseMeasure(3) = .false.
        endif
    else
        UseMeasure(3) = .false.
    endif
    !measure Vol
    if (M4_Active) then
        if (MeasureAll) then
            UseMeasure(4) = .true.
            M4_Val = 0.d0
        elseif (M4_StepFreq > 0 .and. mod(StepNum, M4_StepFreq)==0) then
            UseMeasure(4) = .true.
            M4_Val = 0.d0
        elseif (M4_CycleFreq > 0 .and. mod(CycleNum, M4_CycleFreq)==0) then
            UseMeasure(4) = .true.
            M4_Val = 0.d0
        else
            UseMeasure(4) = .false.
        endif
    else
        UseMeasure(4) = .false.
    endif
    !measure Pressure
    if (M5_Active) then
        if (MeasureAll) then
            UseMeasure(5) = .true.
            M5_Val = 0.d0
        elseif (M5_StepFreq > 0 .and. mod(StepNum, M5_StepFreq)==0) then
            UseMeasure(5) = .true.
            M5_Val = 0.d0
        elseif (M5_CycleFreq > 0 .and. mod(CycleNum, M5_CycleFreq)==0) then
            UseMeasure(5) = .true.
            M5_Val = 0.d0
        else
            UseMeasure(5) = .false.
        endif
    else
        UseMeasure(5) = .false.
    endif
    !measure Virial
    if (M6_Active) then
        if (MeasureAll) then
            UseMeasure(6) = .true.
            M6_Val = 0.d0
        elseif (M6_StepFreq > 0 .and. mod(StepNum, M6_StepFreq)==0) then
            UseMeasure(6) = .true.
            M6_Val = 0.d0
        elseif (M6_CycleFreq > 0 .and. mod(CycleNum, M6_CycleFreq)==0) then
            UseMeasure(6) = .true.
            M6_Val = 0.d0
        else
            UseMeasure(6) = .false.
        endif
    else
        UseMeasure(6) = .false.
    endif
    !measure N_PAA
    if (M7_Active) then
        if (MeasureAll) then
            UseMeasure(7) = .true.
            M7_Val = 0.d0
        elseif (M7_StepFreq > 0 .and. mod(StepNum, M7_StepFreq)==0) then
            UseMeasure(7) = .true.
            M7_Val = 0.d0
        elseif (M7_CycleFreq > 0 .and. mod(CycleNum, M7_CycleFreq)==0) then
            UseMeasure(7) = .true.
            M7_Val = 0.d0
        else
            UseMeasure(7) = .false.
        endif
    else
        UseMeasure(7) = .false.
    endif
    !measure r_PAA
    if (M8_Active) then
        if (MeasureAll) then
            UseMeasure(8) = .true.
            M8_Val = 0.d0
        elseif (M8_StepFreq > 0 .and. mod(StepNum, M8_StepFreq)==0) then
            UseMeasure(8) = .true.
            M8_Val = 0.d0
        elseif (M8_CycleFreq > 0 .and. mod(CycleNum, M8_CycleFreq)==0) then
            UseMeasure(8) = .true.
            M8_Val = 0.d0
        else
            UseMeasure(8) = .false.
        endif
    else
        UseMeasure(8) = .false.
    endif
    !measure x_PAA
    if (M9_Active) then
        if (MeasureAll) then
            UseMeasure(9) = .true.
            M9_Val = 0.d0
        elseif (M9_StepFreq > 0 .and. mod(StepNum, M9_StepFreq)==0) then
            UseMeasure(9) = .true.
            M9_Val = 0.d0
        elseif (M9_CycleFreq > 0 .and. mod(CycleNum, M9_CycleFreq)==0) then
            UseMeasure(9) = .true.
            M9_Val = 0.d0
        else
            UseMeasure(9) = .false.
        endif
    else
        UseMeasure(9) = .false.
    endif
    !measure DUParam
    if (M10_Active) then
        if (MeasureAll) then
            UseMeasure(10) = .true.
            M10_Val = 0.d0
        elseif (M10_StepFreq > 0 .and. mod(StepNum, M10_StepFreq)==0) then
            UseMeasure(10) = .true.
            M10_Val = 0.d0
        elseif (M10_CycleFreq > 0 .and. mod(CycleNum, M10_CycleFreq)==0) then
            UseMeasure(10) = .true.
            M10_Val = 0.d0
        else
            UseMeasure(10) = .false.
        endif
    else
        UseMeasure(10) = .false.
    endif
    !measure DDUParam
    if (M11_Active) then
        if (MeasureAll) then
            UseMeasure(11) = .true.
            M11_Val = 0.d0
        elseif (M11_StepFreq > 0 .and. mod(StepNum, M11_StepFreq)==0) then
            UseMeasure(11) = .true.
            M11_Val = 0.d0
        elseif (M11_CycleFreq > 0 .and. mod(CycleNum, M11_CycleFreq)==0) then
            UseMeasure(11) = .true.
            M11_Val = 0.d0
        else
            UseMeasure(11) = .false.
        endif
    else
        UseMeasure(11) = .false.
    endif
    !measure DWParam
    if (M12_Active) then
        if (MeasureAll) then
            UseMeasure(12) = .true.
            M12_Val = 0.d0
        elseif (M12_StepFreq > 0 .and. mod(StepNum, M12_StepFreq)==0) then
            UseMeasure(12) = .true.
            M12_Val = 0.d0
        elseif (M12_CycleFreq > 0 .and. mod(CycleNum, M12_CycleFreq)==0) then
            UseMeasure(12) = .true.
            M12_Val = 0.d0
        else
            UseMeasure(12) = .false.
        endif
    else
        UseMeasure(12) = .false.
    endif
    !measure DDWParam
    if (M13_Active) then
        if (MeasureAll) then
            UseMeasure(13) = .true.
            M13_Val = 0.d0
        elseif (M13_StepFreq > 0 .and. mod(StepNum, M13_StepFreq)==0) then
            UseMeasure(13) = .true.
            M13_Val = 0.d0
        elseif (M13_CycleFreq > 0 .and. mod(CycleNum, M13_CycleFreq)==0) then
            UseMeasure(13) = .true.
            M13_Val = 0.d0
        else
            UseMeasure(13) = .false.
        endif
    else
        UseMeasure(13) = .false.
    endif
    !measure DUParamDWParam
    if (M14_Active) then
        if (MeasureAll) then
            UseMeasure(14) = .true.
            M14_Val = 0.d0
        elseif (M14_StepFreq > 0 .and. mod(StepNum, M14_StepFreq)==0) then
            UseMeasure(14) = .true.
            M14_Val = 0.d0
        elseif (M14_CycleFreq > 0 .and. mod(CycleNum, M14_CycleFreq)==0) then
            UseMeasure(14) = .true.
            M14_Val = 0.d0
        else
            UseMeasure(14) = .false.
        endif
    else
        UseMeasure(14) = .false.
    endif
    !measure FluctTerm
    if (M15_Active) then
        if (MeasureAll) then
            UseMeasure(15) = .true.
            M15_Val = 0.d0
        elseif (M15_StepFreq > 0 .and. mod(StepNum, M15_StepFreq)==0) then
            UseMeasure(15) = .true.
            M15_Val = 0.d0
        elseif (M15_CycleFreq > 0 .and. mod(CycleNum, M15_CycleFreq)==0) then
            UseMeasure(15) = .true.
            M15_Val = 0.d0
        else
            UseMeasure(15) = .false.
        endif
    else
        UseMeasure(15) = .false.
    endif
    !measure FluctE0
    if (M16_Active) then
        if (MeasureAll) then
            UseMeasure(16) = .true.
            M16_Val = 0.d0
        elseif (M16_StepFreq > 0 .and. mod(StepNum, M16_StepFreq)==0) then
            UseMeasure(16) = .true.
            M16_Val = 0.d0
        elseif (M16_CycleFreq > 0 .and. mod(CycleNum, M16_CycleFreq)==0) then
            UseMeasure(16) = .true.
            M16_Val = 0.d0
        else
            UseMeasure(16) = .false.
        endif
    else
        UseMeasure(16) = .false.
    endif
    !measure FluctA0
    if (M17_Active) then
        if (MeasureAll) then
            UseMeasure(17) = .true.
            M17_Val = 0.d0
        elseif (M17_StepFreq > 0 .and. mod(StepNum, M17_StepFreq)==0) then
            UseMeasure(17) = .true.
            M17_Val = 0.d0
        elseif (M17_CycleFreq > 0 .and. mod(CycleNum, M17_CycleFreq)==0) then
            UseMeasure(17) = .true.
            M17_Val = 0.d0
        else
            UseMeasure(17) = .false.
        endif
    else
        UseMeasure(17) = .false.
    endif

    if (UseMeasure(0)) then
        !measure KEnergy
        M0_Val = KEnergy
        Val0 = M0_Val(0)
        M0_ValSum(0) = M0_ValSum(0) + Val0*Weight
        M0_ValSumSq(0) = M0_ValSumSq(0) + Val0*Val0*Weight
        M0_Count = M0_Count + Weight
    end if

    if (UseMeasure(1)) then
        !measure PEnergy
        M1_Val = PEnergy
        Val0 = M1_Val(0)
        M1_ValSum(0) = M1_ValSum(0) + Val0*Weight
        M1_ValSumSq(0) = M1_ValSumSq(0) + Val0*Val0*Weight
        M1_Count = M1_Count + Weight
    end if

    if (UseMeasure(2)) then
        !measure TEnergy
        M2_Val = TEnergy
        Val0 = M2_Val(0)
        M2_ValSum(0) = M2_ValSum(0) + Val0*Weight
        M2_ValSumSq(0) = M2_ValSumSq(0) + Val0*Val0*Weight
        M2_Count = M2_Count + Weight
    end if

    if (UseMeasure(3)) then
        !measure KTemp
        if (merge(NDOF-Dim, NDOF, ConservesMomentum) <= 0) then
            M3_Val = 0.
        else
            M3_Val = 2.d0 * KEnergy / (kB * merge(NDOF-Dim, NDOF, ConservesMomentum))
        endif
        Val0 = M3_Val(0)
        M3_ValSum(0) = M3_ValSum(0) + Val0*Weight
        M3_ValSumSq(0) = M3_ValSumSq(0) + Val0*Val0*Weight
        M3_Count = M3_Count + Weight
    end if

    if (UseMeasure(4)) then
        !measure Vol
        M4_Val = product(BoxL)
        Val0 = M4_Val(0)
        M4_ValSum(0) = M4_ValSum(0) + Val0*Weight
        M4_ValSumSq(0) = M4_ValSumSq(0) + Val0*Val0*Weight
        M4_Count = M4_Count + Weight
    end if

    if (UseMeasure(5)) then
        !measure Pressure
        M5_Val = (kB*NDOF*TempSet - Virial) / (Dim * product(BoxL))
        Val0 = M5_Val(0)
        M5_ValSum(0) = M5_ValSum(0) + Val0*Weight
        M5_ValSumSq(0) = M5_ValSumSq(0) + Val0*Val0*Weight
        M5_Count = M5_Count + Weight
    end if

    if (UseMeasure(6)) then
        !measure Virial
        M6_Val = Virial
        Val0 = M6_Val(0)
        M6_ValSum(0) = M6_ValSum(0) + Val0*Weight
        M6_ValSumSq(0) = M6_ValSumSq(0) + Val0*Val0*Weight
        M6_Count = M6_Count + Weight
    end if

    if (UseMeasure(7)) then
        !measure N_PAA
        M7_Val = float(NActiveMID(0))
        Val0 = M7_Val(0)
        M7_ValSum(0) = M7_ValSum(0) + Val0*Weight
        M7_ValSumSq(0) = M7_ValSumSq(0) + Val0*Val0*Weight
        M7_Count = M7_Count + Weight
    end if

    if (UseMeasure(8)) then
        !measure r_PAA
        M8_Val = float(NActiveMID(0)) / product(BoxL)
        Val0 = M8_Val(0)
        M8_ValSum(0) = M8_ValSum(0) + Val0*Weight
        M8_ValSumSq(0) = M8_ValSumSq(0) + Val0*Val0*Weight
        M8_Count = M8_Count + Weight
    end if

    if (UseMeasure(9)) then
        !measure x_PAA
        M9_Val = float(NActiveMID(0)) / sum(NActiveMID)
        Val0 = M9_Val(0)
        M9_ValSum(0) = M9_ValSum(0) + Val0*Weight
        M9_ValSumSq(0) = M9_ValSumSq(0) + Val0*Val0*Weight
        M9_Count = M9_Count + Weight
    end if

    if (UseMeasure(10)) then
        !measure DUParam
        M10_Val = DUParam
        M10_ValSum = M10_ValSum + (M10_Val)*Weight
        do v1 = 0, 47 - 1
            do v2 = 0, 47 - 1
                M10_ValSumSq(v1,v2) = M10_ValSumSq(v1,v2) + M10_Val(v1)*M10_Val(v2)*Weight
            enddo
        enddo
        M10_Count = M10_Count + Weight
    end if

    if (UseMeasure(11)) then
        !measure DDUParam
        M11_Val = DDUParam
        M11_ValSum = M11_ValSum + (M11_Val)*Weight
        M11_ValSumSq = M11_ValSumSq + (M11_Val)*(M11_Val)*Weight
        M11_Count = M11_Count + Weight
    end if

    if (UseMeasure(12)) then
        !measure DWParam
        M12_Val = DWParam
        M12_ValSum = M12_ValSum + (M12_Val)*Weight
        M12_ValSumSq = M12_ValSumSq + (M12_Val)*(M12_Val)*Weight
        M12_Count = M12_Count + Weight
    end if

    if (UseMeasure(13)) then
        !measure DDWParam
        M13_Val = DDWParam
        M13_ValSum = M13_ValSum + (M13_Val)*Weight
        M13_ValSumSq = M13_ValSumSq + (M13_Val)*(M13_Val)*Weight
        M13_Count = M13_Count + Weight
    end if

    if (UseMeasure(14)) then
        !measure DUParamDWParam
        do i = 0, NDParam - 1
            do j = 0, NDParam - 1
                M14_Val(i*NDParam + j) = DUParam(i) * DWParam(j)
            enddo
        enddo
        M14_ValSum = M14_ValSum + (M14_Val)*Weight
        M14_ValSumSq = M14_ValSumSq + (M14_Val)*(M14_Val)*Weight
        M14_Count = M14_Count + Weight
    end if

    if (UseMeasure(15)) then
        !measure FluctTerm
        M15_Val = FluctTerm
        Val0 = M15_Val(0)
        M15_ValSum(0) = M15_ValSum(0) + Val0*Weight
        M15_ValSumSq(0) = M15_ValSumSq(0) + Val0*Val0*Weight
        M15_Count = M15_Count + Weight
    end if

    if (UseMeasure(16)) then
        !measure FluctE0
        M16_Val = FluctE0
        Val0 = M16_Val(0)
        M16_ValSum(0) = M16_ValSum(0) + Val0*Weight
        M16_ValSumSq(0) = M16_ValSumSq(0) + Val0*Val0*Weight
        M16_Count = M16_Count + Weight
    end if

    if (UseMeasure(17)) then
        !measure FluctA0
        M17_Val = FluctA0
        Val0 = M17_Val(0)
        M17_ValSum(0) = M17_ValSum(0) + Val0*Weight
        M17_ValSumSq(0) = M17_ValSumSq(0) + Val0*Val0*Weight
        M17_Count = M17_Count + Weight
    end if

end subroutine


subroutine vvquench(NSteps, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    integer, intent(in) :: NSteps
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    real(8) :: dtsq2
    real(8) :: dt2
    real(8) :: idt
    real(8), dimension(0:Dim-1) :: Accel
    integer :: i
    integer :: Step
    logical :: ThisCalcVirial
    logical :: ThisCalcDUParam
    logical :: ThisCalcDWParam

    !velocity verlet quench 1
    dtsq2 = VVQ_TimeStep*VVQ_TimeStep*0.5
    dt2 = 0.5*VVQ_TimeStep
    idt = 1./VVQ_TimeStep

    do Step = 0, NSteps-1
        KEnergy = 0.d0
        TEnergy = KEnergy + PEnergy
        ThisCalcVirial= (CalcVirial .and. Step == NSteps - 1)
        ThisCalcDUParam = (CalcDUParam .and. Step == NSteps - 1)
        ThisCalcDWParam = (CalcDWParam .and. Step == NSteps - 1)

        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            Accel = Force(i,:) * iMass(i)
            Vel(i,:) = 0.d0
            Pos(i,:) = Pos(i,:) + dtsq2*Accel
        enddo

        !no rattle for this system

        call calcenergyforces(0, .true., ThisCalcVirial, ThisCalcDUParam, ThisCalcDWParam, CalcFluct)
    enddo

end subroutine


subroutine vvupdatekenergy()
    implicit none



end subroutine


subroutine vvintegrate(NSteps, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    integer, intent(in) :: NSteps
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    real(8) :: dtsq2
    real(8) :: dt2
    real(8) :: idt
    real(8), dimension(0:Dim-1) :: Accel
    integer :: i
    integer :: j
    integer :: istart
    integer :: istop
    integer :: ind
    integer :: Step
    integer :: m
    real(8) :: NH_wdti0
    real(8) :: NH_wdti1
    real(8) :: NH_wdti2
    real(8) :: NH_wdti3
    real(8) :: NH_kT
    real(8) :: NH_NkT
    real(8) :: NH_scale
    real(8) :: AA
    real(8) :: NH_akin
    integer :: inos
    real(8) :: rn
    real(8), dimension(0:Dim-1) :: ranvec
    real(8) :: dtfreq
    real(8) :: sqrtkt
    real(8) :: langevin1
    real(8) :: langevin2
    real(8), parameter :: randomvelclip = -1.d0
    logical :: ThisCalcVirial
    logical :: ThisCalcDUParam
    logical :: ThisCalcDWParam
    integer :: NActiveAxes
    real(8) :: Beta
    real(8), dimension(0:Dim-1) :: DeltaPos
    real(8), dimension(0:Dim-1) :: CentroidPos
    real(8) :: DeltaVol
    real(8), dimension(0:Dim-1) :: ScaleFactor
    real(8) :: OldVol
    real(8) :: NewVol
    real(8) :: OldE
    real(8) :: NewE
    real(8) :: r
    real(8) :: lnP
    real(8), dimension(0:Dim-1) :: OldBoxL
    logical :: Acc
    logical, dimension(0:Dim-1) :: AxisMask
    real(8) :: Scale

    dtsq2 = VV_TimeStep*VV_TimeStep*0.5
    dt2 = 0.5*VV_TimeStep
    idt = 1./VV_TimeStep

    do Step = 0, NSteps-1

        ThisCalcVirial = (CalcVirial .and. Step == NSteps - 1)
        ThisCalcDUParam = (CalcDUParam .and. Step == NSteps - 1)
        ThisCalcDWParam = (CalcDWParam .and. Step == NSteps - 1)

        if (VV_Barostat == 1 .and. mod(Step, VV_BarostatStepFreq)==0) then
            Beta = 1. / (TempSet * kB)
            !update the attempt
            VV_BarostatNAtt = VV_BarostatNAtt + 1.

            NActiveAxes = count(VV_BarostatUseAxis)
            OldBoxL = BoxL
            OldVol = product(BoxL)
            OldE = PEnergy
            call saveenergystate(0)

            !check if we need to scale independently
            if (.not. VV_BarostatDoIsotropic) then
                !find a random active axis
                call ran2int(NActiveAxes, i)
                j = -1
                do ind = 0, Dim - 1
                    if (VV_BarostatUseAxis(ind)) j = j + 1
                    if (j == i) exit
                enddo
                AxisMask = .false.
                AxisMask(ind) = .true.
            else
                AxisMask = VV_BarostatUseAxis
            endif

            !choose a random volume change
            call ran2(r)
            DeltaVol = VV_BarostatDeltaV * (2.d0 * r - 1.d0)
            NewVol = OldVol + DeltaVol

            if (NewVol > 0. .and. NewVol >= VV_BarostatMinVol .and. NewVol <= VV_BarostatMaxVol) then

                ScaleFactor = merge((NewVol / OldVol)**(1.d0 / dble(count(AxisMask))), 1.d0, AxisMask)
                BoxL = BoxL * ScaleFactor
                OldPos = Pos

                !now scale the molecule centers of mass
                do m = 0, NMol - 1
                    !skip frozen and inactive mols
                    if (MolActive(m) < 1) cycle

                    !find the current centroid
                    istart = MolRange(m)
                    istop = MolRange(m+1) - 1
                    CentroidPos = sum(Pos(istart:istop,:), dim=1) / AtomsPerMol(MolID(m))

                    !find displacement
                    DeltaPos = CentroidPos * (ScaleFactor - 1.d0)

                    !update atom positions
                    do i = istart, istop
                        Pos(i,:) = Pos(i,:) + DeltaPos
                    enddo
                enddo

                !update energy
                call calcenergyforces(0, .false., .false., .false., .false., CalcFluct)
                NewE = PEnergy

                lnP = Beta * (OldE - NewE) - Beta * PresSet * DeltaVol
                lnP = lnP + NActiveMol * log(NewVol / OldVol)
                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif

                if (Acc) then
                    VV_BarostatNAcc = VV_BarostatNAcc + 1.
                else
                    BoxL = OldBoxL
                    Pos = OldPos
                    call revertenergystate(0)
                endif
            endif
        endif

        !Andersen thermostats
        if (VV_Thermostat == 1) then
            !do an andersen massive collision update
            if (mod(VV_AndersenStep, VV_AndersenStepFreq)==0) then
                VV_AndersenStep = 0
                sqrtkt = sqrt(kB * TempSet)
                do i = 0, NAtom-1
                    if (.not. MolActive(MInd(i))==1) cycle
                    call ran2normarray(Dim, ranvec)
                    !clip the limits on the random variate for stability
                    if (randomvelclip > 0.d0) then
                        ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                    endif
                    ranvec = ranvec * sqrtkt / sqrtMass(i)
                    Vel(i,:) = ranvec
                enddo
            endif
            VV_AndersenStep = VV_AndersenStep + 1
        endif

        if (VV_Thermostat == 2) then
            !do an andersen particle collision update
            dtfreq = VV_TimeStep * VV_AndersenCollisionFreq
            sqrtkt = sqrt(kB * TempSet)
            do m = 0, NMol-1
                if (.not. MolActive(m)==1) cycle
                call ran2(rn)
                if (rn < dtfreq) then
                    do i = MolRange(m), MolRange(m+1)-1
                        call ran2normarray(Dim, ranvec)
                        !clip the limits on the random variate for stability
                        if (randomvelclip > 0.d0) then
                            ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                        endif
                        ranvec = ranvec * sqrtkt / sqrtMass(i)
                        Vel(i,:) = ranvec
                    enddo
                endif
            enddo
        endif

        !remove the center of mass
        if (VV_RemoveCOMStepFreq > 0) then
            if (mod(VV_RemoveCOMStep, VV_RemoveCOMStepFreq)==0) then
                VV_RemoveCOMStep = 0
                ranvec = 0.d0
                do i = 0, NAtom-1
                    if (.not. MolActive(MInd(i))==1) cycle
                    ranvec = ranvec + Mass(i) * Vel(i,:)
                enddo
                ranvec = ranvec / NAtom
                do i = 0, NAtom-1
                    if (.not. MolActive(MInd(i))==1) cycle
                    Vel(i,:) = Vel(i,:) - ranvec * iMass(i)
                enddo
            endif
            VV_RemoveCOMStep = VV_RemoveCOMStep + 1
        endif

        !NOSE-HOOVER ROUTINES
        if (VV_Thermostat == 3) then
            !update kinetic energy
            KEnergy = 0.d0
            do i = 0, NAtom-1
                KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
            enddo
            KEnergy = KEnergy * 0.5d0
            TEnergy = KEnergy + PEnergy
            !set frequently used variables
            NH_wdti0 = VV_TimeStep
            NH_wdti1 = VV_TimeStep * 0.5d0
            NH_wdti2 = VV_TimeStep * 0.25d0
            NH_wdti3 = VV_TimeStep * 0.125
            NH_kT = TempSet * kB
            NH_NkT = NH_kT * dble(NDOF - Dim)
            NH_scale = 1.D0
            !get kinetic energy
            NH_akin = 2.d0 * KEnergy
            !update the forces
            VV_Glogs(0) = (NH_akin - NH_NkT) / VV_QMass(0)
            !update the VV_Thermostat velocities
            VV_Vlogs(VV_NH_N-1) = VV_Vlogs(VV_NH_N-1) + VV_Glogs(VV_NH_N-1) * NH_wdti2
            do inos = 1, VV_NH_N - 1
                AA = exp( -NH_wdti3 * VV_Vlogs(VV_NH_N-inos) )
                VV_Vlogs(VV_NH_N-inos-1) = VV_Vlogs(VV_NH_N-inos-1)*AA*AA + NH_wdti2*VV_Glogs(VV_NH_N-inos-1)*AA
            enddo
            !update the particle velocities
            AA = exp( -NH_wdti1 * VV_Vlogs(0) )
            NH_scale = NH_scale * AA
            !update the forces
            VV_Glogs(0) = (NH_scale*NH_scale*NH_akin - NH_NkT) / VV_QMass(0)
            !update the VV_Thermostat positions
            do inos = 0, VV_NH_N - 1
                VV_Xlogs(inos) = VV_Xlogs(inos) + VV_Vlogs(inos) * NH_wdti1
            enddo
            !update the VV_Thermostat velocities
            do inos = 1, VV_NH_N-1
                AA = exp( -NH_wdti3 * VV_Vlogs(inos) )
                VV_Vlogs(inos-1) = VV_Vlogs(inos-1)*AA*AA + NH_wdti2*VV_Glogs(inos-1)*AA
                VV_Glogs(inos) = (VV_QMass(inos-1)*VV_Vlogs(inos-1)*VV_Vlogs(inos-1)-NH_kT) / VV_QMass(inos)
            enddo
            VV_Vlogs(VV_NH_N-1) = VV_Vlogs(VV_NH_N-1) + NH_wdti2*VV_Glogs(VV_NH_N-1)
            !update the particle velocities
            if (NH_scale > 0.) then
                Vel = Vel * NH_scale
            endif
        endif

        !velocity verlet integration part 1
        if (VV_Thermostat == 4) then
            !langevin VV_Thermostat
            !based on Bussi and Parrinello, Physical Review E 75, 056707, 2007
            langevin1 = exp(-0.5d0 * VV_LangevinGamma * VV_TimeStep)
            langevin2 = sqrt((1.d0 - langevin1**2) * (kB * TempSet))
            do i = 0, NAtom-1
                if (.not. MolActive(MInd(i))==1) cycle
                Accel = Force(i,:) * iMass(i)
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
                Pos(i,:) = Pos(i,:) + VV_TimeStep*Vel(i,:) + dtsq2*Accel
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        else
            !normal constant energy dynamics
            do i = 0, NAtom-1
                if (.not. MolActive(MInd(i))==1) cycle
                Accel = Force(i,:) * iMass(i)
                Pos(i,:) = Pos(i,:) + VV_TimeStep*Vel(i,:) + dtsq2*Accel
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        endif

        !no rattle for this system

        call calcenergyforces(0, .true., ThisCalcVirial, ThisCalcDUParam, ThisCalcDWParam, CalcFluct)

        !velocity verlet integration part 2

        if (VV_Thermostat == 4) then
            !langevin VV_Thermostat
            langevin1 = exp(-0.5d0 * VV_LangevinGamma * VV_TimeStep)
            langevin2 = sqrt((1.d0 - langevin1**2) * (kB * TempSet))
            do i = 0, NAtom-1
                if (.not. MolActive(MInd(i))==1) cycle
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                Accel = Force(i,:) * iMass(i)
                Vel(i,:) = Vel(i,:) + dt2*Accel
                Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
            enddo
        else
            !normal constant energy dynamics
            do i = 0, NAtom-1
                if (.not. MolActive(MInd(i))==1) cycle
                Accel = Force(i,:) * iMass(i)
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        endif

        !no rattle for this system

        !NOSE-HOOVER ROUTINES
        if (VV_Thermostat == 3) then
            !update kinetic energy
            KEnergy = 0.d0
            do i = 0, NAtom-1
                KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
            enddo
            KEnergy = KEnergy * 0.5d0
            TEnergy = KEnergy + PEnergy
            !set frequently used variables
            NH_wdti0 = VV_TimeStep
            NH_wdti1 = VV_TimeStep * 0.5d0
            NH_wdti2 = VV_TimeStep * 0.25d0
            NH_wdti3 = VV_TimeStep * 0.125
            NH_kT = TempSet * kB
            NH_NkT = NH_kT * dble(NDOF - Dim)
            NH_scale = 1.D0
            !get kinetic energy
            NH_akin = 2.d0 * KEnergy
            !update the forces
            VV_Glogs(0) = (NH_akin - NH_NkT) / VV_QMass(0)
            !update the VV_Thermostat velocities
            VV_Vlogs(VV_NH_N-1) = VV_Vlogs(VV_NH_N-1) + VV_Glogs(VV_NH_N-1) * NH_wdti2
            do inos = 1, VV_NH_N - 1
                AA = exp( -NH_wdti3 * VV_Vlogs(VV_NH_N-inos) )
                VV_Vlogs(VV_NH_N-inos-1) = VV_Vlogs(VV_NH_N-inos-1)*AA*AA + NH_wdti2*VV_Glogs(VV_NH_N-inos-1)*AA
            enddo
            !update the particle velocities
            AA = exp( -NH_wdti1 * VV_Vlogs(0) )
            NH_scale = NH_scale * AA
            !update the forces
            VV_Glogs(0) = (NH_scale*NH_scale*NH_akin - NH_NkT) / VV_QMass(0)
            !update the VV_Thermostat positions
            do inos = 0, VV_NH_N - 1
                VV_Xlogs(inos) = VV_Xlogs(inos) + VV_Vlogs(inos) * NH_wdti1
            enddo
            !update the VV_Thermostat velocities
            do inos = 1, VV_NH_N-1
                AA = exp( -NH_wdti3 * VV_Vlogs(inos) )
                VV_Vlogs(inos-1) = VV_Vlogs(inos-1)*AA*AA + NH_wdti2*VV_Glogs(inos-1)*AA
                VV_Glogs(inos) = (VV_QMass(inos-1)*VV_Vlogs(inos-1)*VV_Vlogs(inos-1)-NH_kT) / VV_QMass(inos)
            enddo
            VV_Vlogs(VV_NH_N-1) = VV_Vlogs(VV_NH_N-1) + NH_wdti2*VV_Glogs(VV_NH_N-1)
            !update the particle velocities
            if (NH_scale > 0.) then
                Vel = Vel * NH_scale
            endif
        endif

        !update kinetic energy
        KEnergy = 0.d0
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
        enddo
        KEnergy = KEnergy * 0.5d0
        TEnergy = KEnergy + PEnergy

        !update running sums for total energy
        VV_TEnergySum = VV_TEnergySum + TEnergy
        VV_TEnergySqSum = VV_TEnergySqSum + TEnergy*TEnergy

    enddo

end subroutine


subroutine montecarlocycle(NSteps, CalcForce, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    integer, intent(in) :: NSteps
    logical, intent(in) :: CalcForce
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    integer :: StepNum
    real(8) :: mcmoverannum
    real(8) :: Beta
    real(8), dimension(0:Dim-1) :: iBoxL
    integer :: i
    integer :: istart
    integer :: istop
    integer :: m
    integer :: ActiveInd
    integer :: ind
    integer :: ThisAtomsPerMol
    real(8), dimension(0:Dim-1) :: OldAtomPos
    real(8), dimension(0:Dim-1) :: dPos
    real(8), dimension(0:Dim-1) :: CentroidPos
    real(8) :: OldE
    real(8) :: NewE
    real(8) :: r
    real(8) :: lnP
    real(8), dimension(0:Dim-1,0:Dim-1) :: RotMat
    logical :: Acc
    logical :: Rotate
    integer :: ThisMID
    real(8), dimension(0:Dim-1) :: NewPos
    real(8) :: Pacc0
    real(8) :: ProbIns
    real(8) :: ProbDel
    logical :: DoInsert
    integer :: j
    integer :: NActiveAxes
    real(8), dimension(0:Dim-1) :: DeltaPos
    real(8) :: DeltaVol
    real(8), dimension(0:Dim-1) :: ScaleFactor
    real(8) :: OldVol
    real(8) :: NewVol
    real(8), dimension(0:Dim-1) :: OldBoxL
    logical, dimension(0:Dim-1) :: AxisMask
    integer :: DeleteMID
    integer :: DeleteMol
    integer :: DeleteAtomsPerMol
    integer :: InsertMol
    integer :: InsertMID
    integer :: InsertAtomsPerMol
    integer :: Delta1
    integer :: N1
    real(8), dimension(0:Dim-1) :: DeleteCentroidPos
    real(8), dimension(0:Dim-1) :: InsertCentroidPos
    real(8) :: Scale

    Beta = 1.d0 / (kB * TempSet)
    iBoxL = 1.d0 / BoxL
    KEnergy = 0.d0
    TEnergy = 0.d0

    do StepNum = 0, NSteps - 1

        !draw a random number to decide which move to do
        call ran2(mcmoverannum)

        if (mcmoverannum <= MCMoveProbs(0)) then

            !check if active
            if (NActiveMol == 0) then
                MC0_NAtt = MC0_NAtt + 1
                return
            endif

            !pick a random molecule that's active
            call ran2int(NActiveMol, ind)
            ActiveInd = -1
            TargetMol = -1
            do while (ActiveInd < ind)
                TargetMol = TargetMol + 1
                if (MolActive(TargetMol) == 1) ActiveInd = ActiveInd + 1
            enddo
            ThisAtomsPerMol = AtomsPerMol(MolID(TargetMol))

            !check if this is a rigid molecule or not
            if (MolIsRigid(TargetMol)) then

                !rigid molecule; find start and stop atoms
                istart = MolRange(TargetMol)
                istop = MolRange(TargetMol+1) - 1

                !save the old information
                OldE = PEnergy
                OldPos(istart:istop,:) = Pos(istart:istop,:)
                call saveenergystate(2)

                !calculate energy without molecule
                call calcenergyforces(-2, .false., .false., .false., .false., CalcFluct)

                !decide whether to move or rotate
                call ran2(r)
                Rotate = (r > 0.5d0)

                if (Rotate) then

                    MC0_NAtt2 = MC0_NAtt2 + 1
                    !rotate about a random axis
                    call RandomRotMatxyz(MC0_Delta2, RotMat)
                    CentroidPos = sum(Pos(istart:istop,:), dim=1) / ThisAtomsPerMol
                    do i = istart, istop
                        Pos(i,:) = matmul(RotMat, Pos(i,:) - CentroidPos) + CentroidPos
                    enddo

                else

                    MC0_NAtt = MC0_NAtt + 1
                    !translation
                    call ran2array(Dim, dPos)
                    dPos = MC0_Delta * (2.d0 * dPos - 1.d0)
                    do i = istart, istop
                        Pos(i,:) = Pos(i,:) + dPos
                    enddo

                endif

                !calculate energy with molecule
                call calcenergyforces(+2, .false., .false., .false., .false., CalcFluct)
                NewE = PEnergy

                lnP = Beta * (OldE - NewE)
                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif
                if (Acc) then
                    if (Rotate) then
                        MC0_NAcc2 = MC0_NAcc2 + 1
                    else
                        MC0_NAcc = MC0_NAcc + 1
                    endif
                else
                    !revert to old state
                    Pos(istart:istop,:) = OldPos(istart:istop,:)
                    call revertenergystate(2)
                endif

            else

                !this is a non-rigid molecule; pick a random atom to displace
                call ran2int(ThisAtomsPerMol, i)
                TargetAtom = i + MolRange(TargetMol)

                !save the old information
                OldE = PEnergy
                OldAtomPos = Pos(TargetAtom,:)
                call saveenergystate(1)

                !calculate energy without atom
                call calcenergyforces(-1, .false., .false., .false., .false., CalcFluct)

                MC0_NAtt = MC0_NAtt + 1

                !random displacement
                call ran2array(Dim, dPos)
                dPos = MC0_Delta * (2.d0 * dPos - 1.d0)
                Pos(TargetAtom,:) = Pos(TargetAtom,:) + dPos

                !calculate energy with atom
                call calcenergyforces(+1, .false., .false., .false., .false., CalcFluct)
                NewE = PEnergy

                lnP = Beta * (OldE - NewE)
                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif
                if (Acc) then
                    MC0_NAcc = MC0_NAcc + 1
                else
                    !revert to old state
                    Pos(TargetAtom,:) = OldAtomPos
                    call revertenergystate(1)
                endif

            endif

        elseif (mcmoverannum <= MCMoveProbs(1)) then

            !choose whether to insert or delete
            call ran2(r)
            DoInsert = (r > 0.5d0)

            !update the attempt
            MC1_NAtt = MC1_NAtt + 1

            !do an insertion
            if (DoInsert) then

                !make sure there are available molecules to insert
                if (MC1_NInsertMols > 0) then

                    call ran2int(MC1_NInsertMols, ind)
                    TargetMol = MC1_MolMoves(ind)
                    ThisMID = MolID(TargetMol)
                    ThisAtomsPerMol = AtomsPerMol(ThisMID)
                    istart = MolRange(TargetMol)
                    istop = MolRange(TargetMol+1) - 1

                    MC1_NAttMID(ThisMID) = MC1_NAttMID(ThisMID) + 1
                    OldE = PEnergy
                    call saveenergystate(2)

                    !pick a random location
                    call ran2array(Dim, NewPos)
                    NewPos = (NewPos - 0.5d0) * BoxL
                    if (ThisAtomsPerMol > 1) then
                        !pick a random rotation
                        CentroidPos = sum(Pos(istart:istop,:), dim=1) / ThisAtomsPerMol
                        call RandomRotMat3D(3.1415926535897931D0, RotMat)
                        do i = istart, istop
                            Pos(i,:) = matmul(RotMat, Pos(i,:) - CentroidPos) + NewPos
                        enddo
                    else
                        Pos(istart,:) = NewPos
                    endif

                    !update energy for adding this molecule
                    call calcenergyforces(+2, .false., .false., .false., .false., CalcFluct)
                    NewE = PEnergy

                    lnP = Beta * (OldE - NewE) + Beta * MuSet(ThisMID)
                    ProbIns = dble(NInactiveMID(ThisMID)) / dble(MC1_NInsertMols) / product(BoxL)
                    ProbDel = 1.d0 / dble(MC1_NDeleteMols + 1)
                    lnP = lnP + log(ProbDel / ProbIns)
                    Pacc0 = exp(min(0.d0, lnP))
                    lnP = lnP + MC1_BoltzWeights(ThisMID, NActiveMID(ThisMID)+1) - MC1_BoltzWeights(ThisMID, NActiveMID(ThisMID))

                    !update transition matrix
                    MC1_TM(NActiveMol, 1) = MC1_TM(NActiveMol, 1) + (1.d0 - Pacc0)
                    MC1_TM(NActiveMol, 2) = MC1_TM(NActiveMol, 2) + Pacc0

                    if (lnP >= 0) then
                        Acc = .true.
                    else
                        call ran2(r)
                        Acc = (exp(lnP) > r)
                    endif

                    if (Acc) then
                        MC1_NAcc = MC1_NAcc + 1
                        NActiveMol = NActiveMol + 1
                        NActiveMID(ThisMID) = NActiveMID(ThisMID) + 1
                        NInactiveMol = NInactiveMol - 1
                        NInactiveMID(ThisMID) = NInactiveMID(ThisMID) - 1
                        MolActive(TargetMol) = 1
                        MC1_NAccMID(ThisMID) = MC1_NAccMID(ThisMID) + 1
                        !swap the indices of the inserted atom the last possible inserted
                        MC1_MolMoves(ind) = MC1_MolMoves(MC1_NInsertMols - 1)
                        MC1_MolMoves(MC1_NInsertMols - 1) = TargetMol
                        MC1_NInsertMols = MC1_NInsertMols - 1
                        MC1_NDeleteMols = MC1_NDeleteMols + 1
                    else
                        MolActive(TargetMol) = -1
                        call revertenergystate(2)
                    endif

                else

                    !update transition matrix
                    MC1_TM(NActiveMol, 1) = MC1_TM(NActiveMol, 1) + 1.d0

                endif

            !do deletion
            else

                !make sure there are available molecules to delete
                if (MC1_NDeleteMols > 0) then

                    call ran2int(MC1_NDeleteMols, ind)
                    ind = ind + MC1_NInsertMols
                    TargetMol = MC1_MolMoves(ind)
                    ThisMID = MolID(TargetMol)
                    istart = MolRange(TargetMol)
                    istop = MolRange(TargetMol+1) - 1

                    MC1_NAttMID(ThisMID) = MC1_NAttMID(ThisMID) + 1
                    OldE = PEnergy
                    call saveenergystate(2)

                    !update energy for deleting this molecule
                    call calcenergyforces(-2, .false., .false., .false., .false., CalcFluct)
                    NewE = PEnergy

                    lnP = Beta * (OldE - NewE) - Beta * MuSet(ThisMID)
                    ProbDel = 1.d0 / dble(MC1_NDeleteMols)
                    ProbIns = dble(NInactiveMID(ThisMID) + 1) / dble(MC1_NInsertMols+1) / product(BoxL)
                    lnP = lnP + log(ProbIns / ProbDel)
                    Pacc0 = exp(min(0.d0, lnP))
                    lnP = lnP + MC1_BoltzWeights(ThisMID, NActiveMID(ThisMID)-1) - MC1_BoltzWeights(ThisMID, NActiveMID(ThisMID))

                    !update transition matrix
                    MC1_TM(NActiveMol, 1) = MC1_TM(NActiveMol, 1) + (1.d0 - Pacc0)
                    MC1_TM(NActiveMol, 0) = MC1_TM(NActiveMol, 0) + Pacc0

                    if (lnP >= 0) then
                        Acc = .true.
                    else
                        call ran2(r)
                        Acc = (exp(lnP) > r)
                    endif

                    if (Acc) then
                        MC1_NAcc = MC1_NAcc + 1
                        NActiveMol = NActiveMol - 1
                        NActiveMID(ThisMID) = NActiveMID(ThisMID) - 1
                        NInactiveMol = NInactiveMol + 1
                        NInactiveMID(ThisMID) = NInactiveMID(ThisMID) + 1
                        MolActive(TargetMol) = -1
                        MC1_NAccMID(ThisMID) = MC1_NAccMID(ThisMID) + 1
                        !swap the indices of the deleted atom and the first possible deleted
                        MC1_MolMoves(ind) = MC1_MolMoves(MC1_NInsertMols)
                        MC1_MolMoves(MC1_NInsertMols) = TargetMol
                        MC1_NInsertMols = MC1_NInsertMols + 1
                        MC1_NDeleteMols = MC1_NDeleteMols - 1
                    else
                        call revertenergystate(2)
                    endif

                else

                    !update transition matrix
                    MC1_TM(NActiveMol, 1) = MC1_TM(NActiveMol, 1) + 1.d0

                endif

            endif

            !update weights
            do i = 0, NMID-1
                MC1_BoltzWeights(i, NActiveMID(i)) = MC1_BoltzWeights(i, NActiveMID(i)) + MC1_MFactor
            enddo

        elseif (mcmoverannum <= MCMoveProbs(2)) then

            !update the attempt
            MC2_NAtt = MC2_NAtt + 1

            NActiveAxes = count(MC2_UseAxis)
            OldBoxL = BoxL
            OldVol = product(BoxL)
            OldE = PEnergy
            call saveenergystate(0)

            !check if we need to scale independently
            if (.not. MC2_DoIsotropic) then
                !find a random active axis
                call ran2int(NActiveAxes, i)
                j = -1
                do ind = 0, Dim - 1
                    if (MC2_UseAxis(ind)) j = j + 1
                    if (j == i) exit
                enddo
                AxisMask = .false.
                AxisMask(ind) = .true.
            else
                AxisMask = MC2_UseAxis
            endif

            !choose a random volume change
            call ran2(r)
            DeltaVol = MC2_Delta * (2.d0 * r - 1.d0)
            NewVol = OldVol + DeltaVol

            if (NewVol > 0. .and. NewVol >= MC2_MinVol .and. NewVol <= MC2_MaxVol) then

                ScaleFactor = merge((NewVol / OldVol)**(1.d0 / dble(count(AxisMask))), 1.d0, AxisMask)
                BoxL = BoxL * ScaleFactor
                OldPos = Pos

                !now scale the molecule centers of mass
                do m = 0, NMol - 1
                    !skip frozen and inactive mols
                    if (MolActive(m) < 1) cycle

                    !find the current centroid
                    istart = MolRange(m)
                    istop = MolRange(m+1) - 1
                    CentroidPos = sum(Pos(istart:istop,:), dim=1) / AtomsPerMol(MolID(m))

                    !find displacement
                    DeltaPos = CentroidPos * (ScaleFactor - 1.d0)

                    !update atom positions
                    do i = istart, istop
                        Pos(i,:) = Pos(i,:) + DeltaPos
                    enddo
                enddo

                !update energy
                call calcenergyforces(0, .false., .false., .false., .false., CalcFluct)
                NewE = PEnergy

                lnP = Beta * (OldE - NewE) - Beta * PresSet * DeltaVol
                lnP = lnP + NActiveMol * log(NewVol / OldVol)
                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif

                if (Acc) then
                    MC2_NAcc = MC2_NAcc + 1
                else
                    BoxL = OldBoxL
                    Pos = OldPos
                    call revertenergystate(0)
                endif
            endif

        elseif (mcmoverannum <= MCMoveProbs(3)) then

            !update the attempt
            MC3_NAtt = MC3_NAtt + 1

            !pick a random molecule to mutate that's active and type MC3_MID1 or MC3_MID2
            call ran2int(NActiveMID(MC3_MID1) + NActiveMID(MC3_MID2), ind)
            ActiveInd = -1
            DeleteMol = -1
            do while (ActiveInd < ind)
                DeleteMol = DeleteMol + 1
                DeleteMID = MolID(DeleteMol)
                if (MolActive(DeleteMol) == 1 .and. (DeleteMID == MC3_MID1 .or. DeleteMID == MC3_MID2)) then
                    ActiveInd = ActiveInd + 1
                endif
            enddo
            DeleteAtomsPerMol = AtomsPerMol(DeleteMID)

            !find the molecule type to replace with
            if (DeleteMID == MC3_MID1 .and. NInactiveMID(MC3_MID2) > 0) then
                InsertMID = MC3_MID2
                Delta1 = -1
            elseif (DeleteMID == MC3_MID2 .and. NInactiveMID(MC3_MID1) > 0) then
                InsertMID = MC3_MID1
                Delta1 = 1
            else
                InsertMID = -1
            endif

            !order parameter for weights and MC3_TM
            N1 = NActiveMID(MC3_MID1)

            if (InsertMID >= 0) then

                !find a replacement molecule
                call ran2int(NInactiveMID(InsertMID), ind)
                ActiveInd = -1
                InsertMol = -1
                do while (ActiveInd < ind)
                    InsertMol = InsertMol + 1
                    if (MolActive(InsertMol) == -1 .and. MolID(InsertMol) == InsertMID) then
                        ActiveInd = ActiveInd + 1
                    endif
                enddo
                InsertAtomsPerMol = AtomsPerMol(InsertMID)

                !save current energy
                OldE = PEnergy
                call saveenergystate(0)

                !get the center of mass of the original molecule
                istart = MolRange(DeleteMol)
                istop = MolRange(DeleteMol+1) - 1
                DeleteCentroidPos = sum(Pos(istart:istop,:), dim=1) / DeleteAtomsPerMol

                !update energy for deleting the original
                TargetMol = DeleteMol
                MolActive(DeleteMol) = -1
                call calcenergyforces(-2, .false., .false., .false., .false., CalcFluct)

                !now add the new molecule at the same location
                istart = MolRange(InsertMol)
                istop = MolRange(InsertMol+1) - 1
                InsertCentroidPos = sum(Pos(istart:istop,:), dim=1) / InsertAtomsPerMol

                if (InsertAtomsPerMol > 1) then
                    !pick a random rotation
                    call RandomRotMat3D(3.1415926535897931D0, RotMat)
                    do i = istart, istop
                        Pos(i,:) = matmul(RotMat, Pos(i,:) - InsertCentroidPos) + DeleteCentroidPos
                    enddo
                else
                    Pos(istart,:) = DeleteCentroidPos
                endif

                !update energy for adding this molecule
                TargetMol = InsertMol
                MolActive(InsertMol) = 1
                call calcenergyforces(+2, .false., .false., .false., .false., CalcFluct)

                !get the new energy
                NewE = PEnergy

                lnP = Beta * (OldE - NewE) + Beta * (MuSet(InsertMID) - MuSet(DeleteMID))
                Pacc0 = exp(min(0.d0, lnP))
                lnP = lnP + MC3_BoltzWeights(N1 + Delta1) - MC3_BoltzWeights(N1)

                !update transition matrix
                MC3_TM(N1, 1) = MC3_TM(N1, 1) + (1.d0 - Pacc0)
                MC3_TM(N1, 1+Delta1) = MC3_TM(N1, 1+Delta1) + Pacc0

                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif

                if (Acc) then
                    MC3_NAcc = MC3_NAcc + 1
                    NActiveMID(MC3_MID1) = NActiveMID(MC3_MID1) + Delta1
                    NActiveMID(MC3_MID2) = NActiveMID(MC3_MID2) - Delta1
                    NInactiveMID(MC3_MID1) = NInactiveMID(MC3_MID1) - Delta1
                    NInactiveMID(MC3_MID2) = NInactiveMID(MC3_MID2) + Delta1
                    N1 = N1 + Delta1
                    MC3_BoltzWeights(N1) = MC3_BoltzWeights(N1) + MC3_MFactor
                else
                    MolActive(InsertMol) = -1
                    MolActive(DeleteMol) = 1
                    MC3_BoltzWeights(N1) = MC3_BoltzWeights(N1) + MC3_MFactor
                    call revertenergystate(0)
                endif

            else

                !update transition matrix
                MC3_TM(N1, 1) = MC3_TM(N1, 1) + 1.d0

            endif

        end if

    enddo

    call calcenergyforces(0, CalcForce, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)

end subroutine





subroutine ran2(r)
    implicit none
    real(8), intent(out) :: r
    integer, parameter :: NTAB = 32
    integer, parameter :: IM1=2147483563, IM2=2147483399, IMM1=2147483562
    integer, parameter :: IA1=40014, IA2=40692, IQ1=53668, IQ2=52774
    integer, parameter :: IR1=12211, IR2=3791, NDIV=1+IMM1/NTAB  
    real(8), parameter :: AM=1.d0/IM1, EPS=1.2d-7, RNMX=1.d0-EPS
    integer, dimension(3+32) :: state 
    COMMON /ran2data/ state
    integer :: j, k

    if (state(1).le.0) then
        state(1)=max(-state(1),1)
        state(2)=state(1)
        do j=NTAB+8,1,-1
             k=state(1)/IQ1
             state(1)=IA1*(state(1)-k*IQ1)-k*IR1
             if (state(1).lt.0) state(1)=state(1)+IM1
             if (j.le.NTAB) state(3+j)=state(1)
        enddo
        state(3)=state(4)
    endif
    k=state(1)/IQ1
    state(1)=IA1*(state(1)-k*IQ1)-k*IR1
    if (state(1).lt.0) state(1)=state(1)+IM1
    k=state(2)/IQ2
    state(2)=IA2*(state(2)-k*IQ2)-k*IR2
    if (state(2).lt.0) state(2)=state(2)+IM2
    j=1+state(3)/NDIV
    state(3)=state(3+j)-state(2)
    state(3+j)=state(1)
    if(state(3).lt.1)state(3)=state(3)+IMM1
    r=min(AM*state(3),RNMX)
end subroutine

subroutine ran2seed(seedval)
    implicit none
    integer, intent(in) :: seedval
    integer, dimension(3+32) :: state 
    COMMON /ran2data/ state
    state(1) = abs(seedval)
end subroutine

subroutine ran2array(n, rarray)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: rarray
    integer :: i
    real(8) :: r
    do i = 1, n
        call ran2(r)
        rarray(i) = r
    enddo
end subroutine

subroutine ran2int(n, i)
    !returns an integer on [0,n)
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: i
    real(8) :: r
    call ran2(r)
    i = int(real(n) * r)
    i = min(n-1, max(0, i))
end subroutine

subroutine ran2norm(r)
    implicit none
    real(8), intent(out) :: r
    real(8) :: r1, r2, rsq
    real(8), save :: rsaved = 0.1
    logical, save :: hassaved = .false.
    if (hassaved) then
        r = rsaved
        hassaved = .false.
    else
        rsq = 2.d0
        do while (rsq == 0.d0 .or. rsq >= 1.d0)
            call ran2(r1)
            call ran2(r2)
            r1 = 2.d0 * r1 - 1.d0
            r2 = 2.d0 * r2 - 1.d0
            rsq = r1*r1 + r2*r2
        enddo
        rsq = sqrt(-2.d0 * log(rsq) / rsq)
        r = r1 * rsq
        rsaved = r2 * rsq
        hassaved = .true.
    endif
end subroutine

subroutine ran2normarray(n, rarray)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: rarray
    integer :: i
    real(8) :: r
    do i = 1, n
        call ran2norm(r)
        rarray(i) = r
    enddo
end subroutine


subroutine erfc(x, value)
    real(8), intent(in) :: x
    real(8), intent(out) :: value
    !Returns the complementary error function erfc(x) with fractional 
    !error everywhere less than 1:2^10-7
    real(8) :: t, z
    z = abs(x)
    t = 1./(1.+0.5*z)
    value = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(.37409196 + &
          & t*(.09678418 + t*(-.18628806 + t*(.27886807 + t*(-1.13520398 + &
          & t*(1.48851587 + t*(-.82215223 + t*.17087277)))))))))
    if (x < 0.) value = 2. - value
end subroutine erfc

subroutine erf(x, value)
    real(8), intent(in) :: x
    real(8), intent(out) :: value
    !Returns the error function erf(x) with fractional 
    !error everywhere less than 1:2^10-7
    real(8) :: t, z
    z = abs(x)
    t = 1./(1.+0.5*z)
    value = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(.37409196 + &
          & t*(.09678418 + t*(-.18628806 + t*(.27886807 + t*(-1.13520398 + &
          & t*(1.48851587 + t*(-.82215223 + t*.17087277)))))))))
    if (x < 0.) value = 2. - value
    value = 1. - value
end subroutine erf


integer function GetPairIndFromij(i, j)
    integer, intent(in) :: i, j
    if (i > j) then
        GetPairIndFromij = i * (i + 1) / 2 + j
    else
        GetPairIndFromij = j * (j + 1) / 2 + i
    endif      
end function

subroutine GetijFromPairInd(ind, i, j)
    integer, intent(in) :: ind
    integer, intent(out) :: i, j
    real(8) :: r
    r = dble(ind) + 0.001d0
    r = sqrt(1.d0 + 8.d0 * r) * 0.5d0 - 0.5d0
    i = int(r)
    j = ind - i*(i+1)/2
end subroutine

subroutine RandomRotMat3D(MaxAng, RotMat)
    real(8), intent(in) :: MaxAng
    real(8), dimension(3,3) :: RotMat
    real(8) :: theta, sphi, cphi
    real(8) :: Ang, q0, q1, q2, q3
    real(8), dimension(3) :: Vec
    real(8), parameter :: pi = 3.1415926535897931D0
    !get a random rotation angle
    call ran2(Ang)
    Ang = (2.*Ang - 1) * MaxAng
    !get a random vector about which to rotate
    call ran2(cphi)
    cphi = 2.*cphi-1.
    cphi = max(-1.,min(1.,cphi))
    sphi = sqrt(1.-cphi*cphi)
    call ran2(theta)
    theta = 2.d0 * pi * theta
    Vec = (/ cos(theta)*sphi, sin(theta)*sphi, cphi /) 
    !make intermediate variables
    q0 = cos(0.5*Ang)
    Vec = sin(0.5*Ang) * Vec
    q1 = Vec(1)
    q2 = Vec(2)
    q3 = Vec(3)
    !assemble the rotation matrix
    RotMat(1,1) = q0*q0 + q1*q1 - q2*q2 - q3*q3
    RotMat(1,2) = 2.*(q1*q2 + q0*q3)
    RotMat(1,3) = 2.*(q1*q3 - q0*q2)
    RotMat(2,1) = 2.*(q1*q2 - q0*q3)
    RotMat(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3
    RotMat(2,3) = 2.*(q2*q3 + q0*q1)
    RotMat(3,1) = 2.*(q1*q3 + q0*q2)
    RotMat(3,2) = 2.*(q2*q3 - q0*q1)
    RotMat(3,3) = q0*q0 - q1*q1 - q2*q2 + q3*q3
end subroutine

subroutine RandomRotMatxyz(MaxAng, RotMat)
    real(8), intent(in) :: MaxAng
    real(8), dimension(3,3) :: RotMat
    real(8) :: theta, costheta, sintheta, axis
    !get a random rotation angle
    call ran2(theta)
    theta = (2.d0 * theta - 1.d0) * MaxAng
    sintheta = sin(theta)
    costheta = cos(theta)
    !get a random axis about which to rotate
    call ran2(axis)
    RotMat = 0.d0
    if (axis < 0.33333333333d0) then
        RotMat(1,1) = costheta
        RotMat(1,2) = -sintheta
        RotMat(2,1) = sintheta
        RotMat(2,2) = costheta
        RotMat(3,3) = 1.d0
    elseif (axis < 0.66666666667d0) then
        RotMat(2,2) = costheta
        RotMat(2,3) = -sintheta
        RotMat(3,2) = sintheta
        RotMat(3,3) = costheta
        RotMat(1,1) = 1.d0    
    else
        RotMat(1,1) = costheta
        RotMat(1,3) = -sintheta
        RotMat(3,1) = sintheta
        RotMat(3,3) = costheta
        RotMat(2,2) = 1.d0
    endif
end subroutine



function modulehash()
    character(len=40) :: modulehash
    modulehash = '5dd0f3866a418e4ebd43f6843b486c19882be909'
end function


end module

