#!/bin/bash
#
set -x
WOMTOOL_PATH="${WOMTOOL_PATH:-$(command -v womtool)}"
if [[ -z "${WOMTOOL_PATH}" ]]; then
    echo "ERROR: Set WOMTOOL_PATH or install womtool on PATH." >&2
    exit 1
fi

if [[ "${WOMTOOL_PATH}" == *.jar ]]; then
    WOMTOOL=(java -jar "${WOMTOOL_PATH}")
else
    WOMTOOL=("${WOMTOOL_PATH}")
fi

"${WOMTOOL[@]}" validate -l SV_Integration_UltralongAnnotate.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_BndAnnotate.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_BndGetTrainingIntervals.wdl
"${WOMTOOL[@]}" validate -l UltralongAnnotate.wdl
"${WOMTOOL[@]}" validate -l UltralongRecordsInTrack.wdl
"${WOMTOOL[@]}" validate -l InsRemap.wdl
"${WOMTOOL[@]}" validate -l UltralongGetTrainingIntervalsSvim.wdl
"${WOMTOOL[@]}" validate -l SvimAsm.wdl
"${WOMTOOL[@]}" validate -l UltralongGetTrainingIntervals.wdl
"${WOMTOOL[@]}" validate -l UltralongCanonizeDipcall.wdl
"${WOMTOOL[@]}" validate -l UltralongScore.wdl
"${WOMTOOL[@]}" validate -l UltralongMerge.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_UltralongAnalysis.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_RegenotypingAnalysis.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_PlotPrme.wdl
"${WOMTOOL[@]}" validate -l InvestigateMaleSamples2.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_PlotHwe.wdl
"${WOMTOOL[@]}" validate -l InvestigateMaleSamples.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage15.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage14.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage13.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage12.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage11.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage9.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage8.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage9_families.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage9_trios.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage7.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage6.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage5.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage3.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_BuildTrainingResource.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage1.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_PlotHwe_SNVs.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage10.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage8_Prime.wdl
"${WOMTOOL[@]}" validate -l CollapseHapsVcf.wdl
"${WOMTOOL[@]}" validate -l Trgt2Kanpig.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage2.wdl
"${WOMTOOL[@]}" validate -l RegenotypeKanpigHapsVcf.wdl
"${WOMTOOL[@]}" validate -l BuildKanpigHapVcf.wdl
"${WOMTOOL[@]}" validate -l RegenotypeShapeit4.wdl
"${WOMTOOL[@]}" validate -l PhabRegenotypedCohort.wdl
"${WOMTOOL[@]}" validate -l SV_Integration_Workpackage5_Bnd.wdl
"${WOMTOOL[@]}" validate -l BenchCohortSamples_PersonalizedCohortVcf.wdl
"${WOMTOOL[@]}" validate -l HGSVC3Dipcall2BAMs.wdl
"${WOMTOOL[@]}" validate -l GetCompositeSvs.wdl
"${WOMTOOL[@]}" validate -l PersonalizedCohortVcf.wdl
"${WOMTOOL[@]}" validate -l Bam2Fastq.wdl
"${WOMTOOL[@]}" validate -l GetNCalls2.wdl
"${WOMTOOL[@]}" validate -l GetNCalls.wdl
"${WOMTOOL[@]}" validate -l BenchCohortTriosSquish5.wdl
"${WOMTOOL[@]}" validate -l CountTrSubstratifications.wdl
"${WOMTOOL[@]}" validate -l BenchCohortTriosSquish4.wdl
"${WOMTOOL[@]}" validate -l SubsetToAncestry.wdl
"${WOMTOOL[@]}" validate -l FilterTruvariIntersample2.wdl
"${WOMTOOL[@]}" validate -l BenchCohortSamples_windowed.wdl
"${WOMTOOL[@]}" validate -l FilterTruvariIntersample.wdl
"${WOMTOOL[@]}" validate -l BenchCohortTriosSquish3.wdl
"${WOMTOOL[@]}" validate -l MapCCSPhase2.wdl
"${WOMTOOL[@]}" validate -l MapCCSPhase2Prime.wdl
"${WOMTOOL[@]}" validate -l MapR10Phase2.wdl
"${WOMTOOL[@]}" validate -l SubsampleAlignedBam.wdl
"${WOMTOOL[@]}" validate -l Workpackage9Squish.wdl
"${WOMTOOL[@]}" validate -l DeNovoByRegion.wdl
"${WOMTOOL[@]}" validate -l PhabTrios.wdl
"${WOMTOOL[@]}" validate -l TestKanpigIntersample.wdl
"${WOMTOOL[@]}" validate -l TestBetaBinomial3.wdl
"${WOMTOOL[@]}" validate -l TestBetaBinomial_Merge.wdl
"${WOMTOOL[@]}" validate -l BenchCohortSamplesBetaBinomial.wdl
"${WOMTOOL[@]}" validate -l TestBetaBinomial2.wdl
"${WOMTOOL[@]}" validate -l TestBetaBinomial.wdl
"${WOMTOOL[@]}" validate -l Workpackage9Palt.wdl
"${WOMTOOL[@]}" validate -l GetGenotypingPriors.wdl
"${WOMTOOL[@]}" validate -l FixUnsupportedGts.wdl
"${WOMTOOL[@]}" validate -l AnalyzeGtAdMatrix.wdl
"${WOMTOOL[@]}" validate -l GetGtAdMatrix.wdl
"${WOMTOOL[@]}" validate -l Workpackage9Subsets.wdl
"${WOMTOOL[@]}" validate -l GetMapqDistribution.wdl
"${WOMTOOL[@]}" validate -l BenchCohortTriosSquish2.wdl
"${WOMTOOL[@]}" validate -l BenchCohortTriosSquish.wdl
"${WOMTOOL[@]}" validate -l PlotHweFocusedTruvariCollapse.wdl
"${WOMTOOL[@]}" validate -l BenchCohortTrios.wdl
"${WOMTOOL[@]}" validate -l PlotHweFocusedAc.wdl
"${WOMTOOL[@]}" validate -l BenchCohortSamples.wdl
"${WOMTOOL[@]}" validate -l PlotHweFocused.wdl
"${WOMTOOL[@]}" validate -l MapR10Phase2ScatteredLrhq.wdl
"${WOMTOOL[@]}" validate -l MapR10Phase2Scattered.wdl
"${WOMTOOL[@]}" validate -l AddReadGroup.wdl
"${WOMTOOL[@]}" validate -l DownloadAssembly.wdl
"${WOMTOOL[@]}" validate -l SubsampleSimple.wdl
"${WOMTOOL[@]}" validate -l ReadLengthDistribution.wdl
"${WOMTOOL[@]}" validate -l GetLongCalls.wdl
"${WOMTOOL[@]}" validate -l CheckDeNovo.wdl
"${WOMTOOL[@]}" validate -l BenchHprcSamples.wdl
"${WOMTOOL[@]}" validate -l GetPresentCalls.wdl
"${WOMTOOL[@]}" validate -l InterCenterMerge.wdl
"${WOMTOOL[@]}" validate -l CheckMendelian.wdl
"${WOMTOOL[@]}" validate -l PlotHwe.wdl
"${WOMTOOL[@]}" validate -l QcPlots.wdl
"${WOMTOOL[@]}" validate -l Workpackage13.wdl
"${WOMTOOL[@]}" validate -l Workpackage12.wdl
"${WOMTOOL[@]}" validate -l Workpackage11.wdl
"${WOMTOOL[@]}" validate -l Workpackage10.wdl
"${WOMTOOL[@]}" validate -l Workpackage9.wdl
"${WOMTOOL[@]}" validate -l Workpackage8.wdl
"${WOMTOOL[@]}" validate -l Workpackage7.wdl
"${WOMTOOL[@]}" validate -l Workpackage6.wdl
"${WOMTOOL[@]}" validate -l Workpackage5.wdl
"${WOMTOOL[@]}" validate -l Workpackage4.wdl
"${WOMTOOL[@]}" validate -l Workpackage3.wdl
"${WOMTOOL[@]}" validate -l Workpackage2.wdl
"${WOMTOOL[@]}" validate -l Workpackage1.wdl
"${WOMTOOL[@]}" validate -l InterCenterBench.wdl
"${WOMTOOL[@]}" validate -l PasteGTs.wdl
"${WOMTOOL[@]}" validate -l KanpigMerged.wdl
"${WOMTOOL[@]}" validate -l RemoveSamples.wdl
"${WOMTOOL[@]}" validate -l TruvariIntersamplePhase2.wdl
"${WOMTOOL[@]}" validate -l Split.wdl
"${WOMTOOL[@]}" validate -l HGSVC3ExtractHapsFromAssemblies.wdl
"${WOMTOOL[@]}" validate -l DipcallPhase2.wdl
"${WOMTOOL[@]}" validate -l FilterIntrasampleDevPhase2.wdl
"${WOMTOOL[@]}" validate -l Kanpig.wdl
"${WOMTOOL[@]}" validate -l TruvariIntrasample.wdl
"${WOMTOOL[@]}" validate -l Resolve.wdl
"${WOMTOOL[@]}" validate -l PAV2SVs.wdl
