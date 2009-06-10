#Copyright or © or Copr. CNRS, (November 16, 2004-2009)
#
#This software is a computer program whose purpose is to provide classes
#for population genetics analysis.
#
#This software is governed by the CeCILL  license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL
#license as circulated by CEA, CNRS and INRIA at the following URL
#"http://www.cecill.info". 
#
#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 
#
#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 
#
#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL license and that you accept its terms.

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
DarwinDon.h\
DarwinVarSingle.h\
DataSet.h\
DataSetTools.h\
Date.h\
GeneralExceptions.h\
GeneMapperCsvExport.h\
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
MonolocusGenotypeTools.h\
MultiAlleleMonolocusGenotype.h\
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
AnalyzedLoci.cpp\
AnalyzedSequences.cpp\
BasicAlleleInfo.cpp\
BiAlleleMonolocusGenotype.cpp\
DarwinDon.cpp\
DarwinVarSingle.cpp\
DataSet.cpp\
DataSetTools.cpp\
Date.cpp\
GeneralExceptions.cpp\
GeneMapperCsvExport.cpp\
Genepop.cpp\
Genetix.cpp\
Group.cpp\
Individual.cpp\
LocusInfo.cpp\
MonoAlleleMonolocusGenotype.cpp\
MonolocusGenotypeTools.cpp\
MultiAlleleMonolocusGenotype.cpp\
MultilocusGenotype.cpp\
MultilocusGenotypeStatistics.cpp\
MultiSeqIndividual.cpp\
PolymorphismMultiGContainer.cpp\
PolymorphismMultiGContainerTools.cpp\
PolymorphismSequenceContainer.cpp\
PolymorphismSequenceContainerTools.cpp\
PopgenlibIO.cpp\
SequenceStatistics.cpp

lib_LTLIBRARIES = libbpp-popgen.la
libbpp_popgen_la_SOURCES = $(H_FILES) $(CPP_FILES)
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current: The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age: The difference between the newest and oldest interfaces that this library implements.
# In other words, the library implements all the interface numbers in the range from number current - age to current.
libbpp_popgen_la_LDFLAGS = -version-info 5:0:2

library_includedir=$(includedir)/Pop
library_include_HEADERS = $(H_FILES)
