H_FILES = \
AlleleInfoContainerExceptions.h\
AlleleInfo.h\
AnalyzedLoci.h\
Coord.h\
CoordsTools.h\
DataSet.h\
Date.h\
DiploidMonolocusGenotype.h\
Genotype.h\
Group.h\
HaplodiploidMonolocusGenotype.h\
HaploidMonolocusGenotype.h\
Individual.h\
Locality.h\
LocusInfoContainerExceptions.h\
LocusInfo.h\
MonolocusGenotype.h

CPP_FILES = \
AlleleInfoContainerExceptions.cpp\
AlleleInfo.cpp\
AnalyzedLoci.cpp\
DataSet.cpp\
Date.cpp\
DiploidMonolocusGenotype.cpp\
Genotype.cpp\
Group.cpp\
HaplodiploidMonolocusGenotype.cpp\
HaploidMonolocusGenotype.cpp\
Individual.cpp\
LocusInfoContainerExceptions.cpp\
LocusInfo.cpp\
MonolocusGenotype.cpp

lib_LTLIBRARIES = libpop.la
libpop_la_SOURCES = $(H_FILES) $(CPP_FILES)
libpop_la_LDFLAGS = -version-info 1:0:0

library_includedir=$(includedir)/Pop
library_include_HEADERS = $(H_FILES)
