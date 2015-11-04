SUBDIRS = 

if build_monitors
SUBDIRS	+=	\
	  BicoMon BitTest blrms burstMon callineMon Cumulus DEnvCorr dewarMon \
	  DuoTone dvTest endtimes eqMon GainMon glitchMon HistCompr IRIG-B \
	  InspiralMon kleineWelle LIGOLwMon LightMon LineMonitor LockLoss \
	  MultiVolt NoiseFloorMonitor OmegaMon Omega_c PCalMon \
	  PSLmon PhotonCal PlaneMon PulsarMon RangeMonitor SegGener \
	  SenseMonitor ShapeMon SixtyHertzMon SpectrumArchiver \
	  SpectrumFold Station StochMon StrainbandsMon \
	  suspensionMon TimeMon TrigDsply TrigSpec WaveMon

#
# --> difficult
SUBDIRS +=	RayleighMonitor

else
if build_omega
SUBDIRS +=	\
	Omega_c
endif
endif
#
#  The following monitors are still available, but are not built by default.
OBSOLETE= hMon hoftMon LscMonitor NdasWriter PDnMon SatScan Slice2 StrainWatch
