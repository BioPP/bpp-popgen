H_FILES = \
AlleleInfo.h\
AlleleInfoContainerExceptions.h\
AnalyzedLoci.h\
Coord.h\
Date.h\
Individual.h\
Locality.h\
LocusInfo.h\
LocusInfoContainerExceptions.h\
MonolocusGenotype.h

CPP_FILES = \
AlleleInfo.cpp\
AlleleInfoContainerExceptions.cpp\
AnalyzedLoci.cpp\
Date.cpp\
Individual.cpp\
LocusInfo.cpp\
LocusInfoContainerExceptions.cpp\
MonolocusGenotype.cpp

lib_LTLIBRARIES = libpop.la
libpop_la_SOURCES = $(H_FILES) $(CPP_FILES)
libpop_la_LDFLAGS = -version-info 1:0:0

library_includedir=$(includedir)/Pop
library_include_HEADERS = $(H_FILES)
