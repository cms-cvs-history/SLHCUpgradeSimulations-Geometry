# Phase 1 R34F16 minbias pileup files
# E34 cm-2s-1
import FWCore.ParameterSet.Config as cms

# this is the configuration to model pileup in the design LHC (10**34)
from SimGeneral.MixingModule.mixObjects_cfi import *
mix = cms.EDProducer("MixingModule",
    LabelPlayback = cms.string(''),
    maxBunch = cms.int32(3),
    minBunch = cms.int32(-5), ## in terms of 25 ns

    bunchspace = cms.int32(25), ## nsec
    mixProdStep1 = cms.bool(False),
    mixProdStep2 = cms.bool(False),

    playback = cms.untracked.bool(False),
    useCurrentProcessOnly = cms.bool(False),
                   
    input = cms.SecSource("PoolSource",
    nbPileupEvents = cms.PSet(
            sigmaInel = cms.double(80.0),
            Lumi = cms.double(10.0)
        ),
        type = cms.string('poisson'),
    sequential = cms.untracked.bool(False),
        fileNames = cms.untracked.vstring('file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_1_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_2_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_3_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_4_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_5_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_6_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_7_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_8_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_9_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_10_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_11_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_12_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_13_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_14_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_15_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_16_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_17_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_18_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_19_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_20_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_21_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_22_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_23_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_24_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_25_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_26_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_27_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_28_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_29_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_30_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_31_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_32_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_33_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_34_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_35_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_36_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_37_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_38_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_39_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_40_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_41_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_42_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_43_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_44_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_45_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_46_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_47_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_48_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_49_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_50_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_51_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_52_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_53_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_54_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_55_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_56_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_57_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_58_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_59_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_60_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_61_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_62_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_63_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_64_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_65_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_66_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_67_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_68_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_69_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_70_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_71_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_72_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_73_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_74_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_75_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_76_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_77_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_78_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_79_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_80_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_81_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_82_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_83_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_84_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_85_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_86_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_87_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_88_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_89_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_90_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_91_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_92_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_93_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_94_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_95_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_96_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_97_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_98_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_99_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_100_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_101_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_102_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_103_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_104_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_105_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_106_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_107_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_108_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_109_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_110_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_111_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_112_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_113_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_114_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_115_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_116_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_117_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_118_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_119_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_120_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_121_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_122_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_123_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_124_1.root',
        'file:/uscms_data/d2/lpcsimu/slhc/phase1/R34F16/minbias4pu/MinBias_cfi_py_GEN_SIM_125_1.root'
    )
    ),
    mixObjects = cms.PSet(
        mixCH = cms.PSet(
            mixCaloHits
        ),
        mixTracks = cms.PSet(
            mixSimTracks
        ),
        mixVertices = cms.PSet(
            mixSimVertices
        ),
        mixSH = cms.PSet(
            mixSimHits
        ),
        mixHepMC = cms.PSet(
            mixHepMCProducts
        )
    )
)
