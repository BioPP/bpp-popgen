H_FILES = \
AbstractIDataSet.h\
AbstractODataSet.h\
AlleleInfo.h\
AnalyzedLoci.h\
AnalyzedSequences.h\
BasicAlleleInfo.h\
BiAlleleMonolocusGenotype.h\
Coord.h\
CoordsTools.h\
DataSet.h\
DataSetTools.h\
Date.h\
GeneralExceptions.h\
Genepop.h\
Genetix.h\
Group.h\
Individual.h\
IDataSet.h\
IODataSet.h\
Locality.h\
LocusInfo.h\
MonoAlleleMonolocusGenotype.h\
MonolocusGenotype.h\
MultilocusGenotype.h\
MultilocusGenotypeStatistics.h\
MultiSeqIndividual.h\
ODataSet.h\
PolymorphismMultiGContainer.h\
PolymorphismMultiGContainerTools.h\
PolymorphismSequenceContainer.h\
PolymorphismSequenceContainerTools.h\
PopgenlibIO.h\
SequenceStatistics.h

CPP_FILES = \
AbstractIDataSet.cpp\
AbstractODataSet.cpp\
AlleleInfo.cpp\
AnalyzedLoci.cpp\
AnalyzedSequences.cpp\
BasicAlleleInfo.cpp\
BiAlleleMonolocusGenotype.cpp\
DataSet.cpp\
DataSetTools.cpp\
Date.cpp\
GeneralExceptions.cpp\
Genepop.cpp\
Genetix.cpp\
Group.cpp\
Individual.cpp\
IDataSet.cpp\
IODataSet.cpp\
ODataSet.cpp\
LocusInfo.cpp\
MonoAlleleMonolocusGenotype.cpp\
MonolocusGenotype.cpp\
MultilocusGenotype.cpp\
MultilocusGenotypeStatistics.cpp\
MultiSeqIndividual.cpp\
PolymorphismMultiGContainer.cpp\
PolymorphismMultiGContainerTools.cpp\
PolymorphismSequenceContainer.cpp\
PolymorphismSequenceContainerTools.cpp\
PopgenlibIO.cpp\
SequenceStatistics.cpp

lib_LTLIBRARIES = libpop.la
libpop_la_SOURCES = $(H_FILES) $(CPP_FILES)
libpop_la_LDFLAGS = -version-info 0:1:0

library_includedir=$(includedir)/Pop
library_include_HEADERS = $(H_FILES)
