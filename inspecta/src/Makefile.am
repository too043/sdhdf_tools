LIBS = -lcfitsio -lhdf5 -lerfa -lcalceph

bin_PROGRAMS = sdhdf_describe \
  sdhdf_checkFile \
  sdhdf_simulate \
  sdhdf_extractBand \
  sdhdf_extractDump \
  sdhdf_quickdump \
  sdhdf_onoff \
  sdhdf_fluxDensity \
  sdhdf_sum \
  sdhdf_modify \
  sdhdf_calibrate \
  sdhdf_calibrateOld \
  sdhdf_join \
  sdhdf_joinBands \
  sdhdf_joinFiles \
  sdhdf_identify \
  sdhdf_plotNoiseCal \
  sdhdf_plotPulsar \
  sdhdf_plotPulsarSearch \
  sdhdf_tsys \
  sdhdf_convertFromFITS \
  sdhdf_convertPulsarSearchMode \
  sdhdf_convertTo \
  sdhdf_convertToPSRFITS \
  sdhdf_convertFromPsr \
  sdhdf_listLines

COMMON_SRCS = sdhdfProc_maths.c \
  sdhdfProc_fileManipulation.c \
  sdhdfProc_metadata.c \
  sdhdfProc_calibration.c \
  sdhdfProc_ephemeris.c \
  T2toolkit.c \
  TKfit.c

sdhdf_describe_SOURCES = sdhdf_describe.c \
  $(COMMON_SRCS)

sdhdf_checkFile_SOURCES = sdhdf_checkFile.c \
  $(COMMON_SRCS)

sdhdf_simulate_SOURCES = sdhdf_simulate.c \
  $(COMMON_SRCS)

sdhdf_extractBand_SOURCES = sdhdf_extractBand.c \
  $(COMMON_SRCS)

sdhdf_extractDump_SOURCES = sdhdf_extractDump.c \
  $(COMMON_SRCS)

sdhdf_quickdump_SOURCES = sdhdf_quickdump.c \
  $(COMMON_SRCS)

sdhdf_onoff_SOURCES = sdhdf_onoff.c \
  $(COMMON_SRCS)

sdhdf_fluxDensity_SOURCES = sdhdf_fluxDensity.c \
  $(COMMON_SRCS)

sdhdf_sum_SOURCES = sdhdf_sum.c \
  $(COMMON_SRCS)

sdhdf_modify_SOURCES = sdhdf_modify.c \
  $(COMMON_SRCS)

sdhdf_calibrate_SOURCES = sdhdf_calibrate.c \
  sdhdf_opacity.c \
  $(COMMON_SRCS)

sdhdf_calibrateOld_SOURCES = sdhdf_calibrateOld.c \
  $(COMMON_SRCS)

sdhdf_join_SOURCES = sdhdf_join.c \
  $(COMMON_SRCS)

sdhdf_joinFiles_SOURCES = sdhdf_joinFiles.c \
  $(COMMON_SRCS)

sdhdf_joinBands_SOURCES = sdhdf_joinBands.c \
  $(COMMON_SRCS)

sdhdf_identify_SOURCES = sdhdf_identify.c \
  $(COMMON_SRCS)

sdhdf_convertFromFITS_SOURCES = sdhdf_convertFromFITS.c \
  $(COMMON_SRCS)

sdhdf_convertPulsarSearchMode_SOURCES = sdhdf_convertPulsarSearchMode.c \
  $(COMMON_SRCS)

sdhdf_convertTo_SOURCES = sdhdf_convertTo.c \
  $(COMMON_SRCS)

sdhdf_convertToPSRFITS_SOURCES = sdhdf_convertToPSRFITS.c \
  $(COMMON_SRCS)

sdhdf_convertFromPsr_SOURCES = sdhdf_convertFromPsr.c \
  $(COMMON_SRCS)

sdhdf_plotNoiseCal_SOURCES = sdhdf_plotNoiseCal.c \
  $(COMMON_SRCS)

sdhdf_plotPulsar_SOURCES = sdhdf_plotPulsar.c \
  $(COMMON_SRCS)

sdhdf_plotPulsarSearch_SOURCES = sdhdf_plotPulsarSearch.c \
  $(COMMON_SRCS)

sdhdf_tsys_SOURCES = sdhdf_tsys.c \
  $(COMMON_SRCS)

sdhdf_listLines_SOURCES = sdhdf_listLines.c \
  $(COMMON_SRCS)

LDADD = $(LIBS)


if HAVE_PGPLOT

bin_PROGRAMS += sdhdf_primaryFluxCal \
  sdhdf_gainCurve \
  sdhdf_applyCal \
  sdhdf_calProc \
  sdhdf_quicklook \
  sdhdf_processBWS \
  sdhdf_fitBaseline \
  sdhdf_statistics \
  sdhdf_pointing \
  sdhdf_autoFlag \
  sdhdf_rfiTime \
  sdhdf_waterfall \
  sdhdf_flag \
  sdhdf_plotScan \
  sdhdf_plotWide \
  sdhdf_plotMultiSpec \
  sdhdf_map \
  sdhdf_weightSpectra \
  sdhdf_plotSpectrum \
  sdhdf_model 

sdhdf_primaryFluxCal_SOURCES = sdhdf_primaryFluxCal.c \
  $(COMMON_SRCS)

sdhdf_gainCurve_SOURCES = sdhdf_gainCurve.c \
  $(COMMON_SRCS)

sdhdf_applyCal_SOURCES = sdhdf_applyCal.c \
  $(COMMON_SRCS)

sdhdf_calProc_SOURCES = sdhdf_calProc.c \
  $(COMMON_SRCS)

sdhdf_pointing_SOURCES = sdhdf_pointing.c \
  $(COMMON_SRCS)

sdhdf_autoFlag_SOURCES = sdhdf_autoFlag.c \
  $(COMMON_SRCS)

sdhdf_rfiTime_SOURCES = sdhdf_rfiTime.c \
  $(COMMON_SRCS)

sdhdf_waterfall_SOURCES = sdhdf_waterfall.c \
  $(COMMON_SRCS)

sdhdf_flag_SOURCES = sdhdf_flag.c \
  $(COMMON_SRCS)

sdhdf_map_SOURCES = sdhdf_map.c \
  $(COMMON_SRCS)

sdhdf_weightSpectra_SOURCES = sdhdf_weightSpectra.c \
  $(COMMON_SRCS)

sdhdf_plotScan_SOURCES = sdhdf_plotScan.c \
  $(COMMON_SRCS)

sdhdf_quicklook_SOURCES = sdhdf_quicklook.c \
  $(COMMON_SRCS)

sdhdf_processBWS_SOURCES = sdhdf_processBWS.c \
  $(COMMON_SRCS)

sdhdf_fitBaseline_SOURCES = sdhdf_fitBaseline.c \
  $(COMMON_SRCS)

sdhdf_statistics_SOURCES = sdhdf_statistics.c \
  $(COMMON_SRCS)

sdhdf_plotWide_SOURCES = sdhdf_plotWide.c \
  $(COMMON_SRCS)

sdhdf_plotMultiSpec_SOURCES = sdhdf_plotMultiSpec.c \
  $(COMMON_SRCS)

sdhdf_plotSpectrum_SOURCES = sdhdf_plotSpectrum.c \
  $(COMMON_SRCS)

sdhdf_model_SOURCES = sdhdf_model.c \
  $(COMMON_SRCS)

AM_CFLAGS = @PGPLOT_CFLAGS@
LDADD += @PGPLOT_LIBS@

endif
