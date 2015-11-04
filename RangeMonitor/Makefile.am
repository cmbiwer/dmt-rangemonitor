## Process this file with automake to produce Makefile.in
## Author: Junwei Cao

MAINTAINERCLEANFILES 	= Makefile.in

bin_PROGRAMS		= RangeMonitor
RangeMonitor_SOURCES	= \
			RangeMonitor.cc \
			Integrand.cc

noinst_HEADERS		= \
			Integrand.hh \
			Integrate.hh \
			Range.hh \
			RangeSummary.hh \
			RangeMonitor.hh

AM_CPPFLAGS		= \
			-I$(top_builddir)/config \
			-I$(top_srcdir)/Base/time \
			-I$(top_srcdir)/Base/complex \
			-I$(top_srcdir)/Base/misc \
			-I$(top_srcdir)/Containers \
			-I$(top_srcdir)/Containers/TSeries \
			-I$(top_srcdir)/Containers/DVector \
			-I$(top_srcdir)/Containers/CWVec \
			-I$(top_srcdir)/Containers/FilterIO \
			-I$(top_srcdir)/Containers/FSpectrum \
			-I$(top_srcdir)/Containers/Histogram1 \
			-I$(top_srcdir)/SignalProcessing/DecimateBy2 \
			-I$(top_srcdir)/SignalProcessing/FilterBase \
			-I$(top_srcdir)/SignalProcessing/Window \
			-I$(top_srcdir)/SignalProcessing/IIRFilter \
			-I$(top_srcdir)/SignalProcessing/FIRFilter \
			-I$(top_srcdir)/SignalProcessing/FilterDesign \
			-I$(top_srcdir)/SignalProcessing/FDFilter \
			-I$(top_srcdir)/IO \
			-I$(top_srcdir)/IO/parseline \
			-I$(top_srcdir)/Triggers/trig \
			-I$(top_srcdir)/Triggers/Par \
			-I$(top_srcdir)/Services/TrigClient \
			-I$(top_srcdir)/Services/MonServer \
			-I$(top_srcdir)/Services/AppServer \
			-I$(top_srcdir)/Services/NameClient \
			-I$(top_srcdir)/Services/AppClient \
			-I$(top_srcdir)/Services/dacc \
			-I$(top_srcdir)/Services/trend \
			-I$(top_srcdir)/DMT/DatEnv \
			-I$(top_srcdir)/DMT/osc3 \
			-I$(top_srcdir)/DMT/EasyCalibrate \
			$(framecppcflags)

LDADD                  	= \
                        ../../Base/libgdsbase.la \
			../../Containers/libgdscntr.la \
			../../Triggers/trig/libgdstrig.la \
			../../SignalProcessing/libdmtsigp.la \
			../../IO/frameutils/libframeutil.la \
			../../IO/parseline/libparsl.la \
			../../IO/html/libhtml.la \
			../../IO/lxr/liblxr.la \
			../../IO/xsil/libxsil.la \
			../../Services/dacc/libframeio.la \
			../../Services/libmonitor.la \
			../../Services/libtclient.la \
			../../DMT/DatEnv/libdmtenv.a \
			../../DMT/osc3/libosc.la \
			../../DMT/EasyCalibrate/libezcalib.la \
			$(metaiolibs) $(framecpplibs)
if online
LDADD			+= \
			../../IO/lmsg/liblmsg.la \
			../../IO/lsmp/liblsmp.la
endif
AM_CXXFLAGS		= $(WARNINGCXXFLAGS) $(opt)
if offline
DEFS                    += -DDMTOFFLINE
endif
