# CMake script for Bio++ PopGen
# Authors:
#   Sylvain Gaillard
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 21/08/2009

# File list
set (CPP_FILES
  Bpp/PopGen/BasicAlleleInfo.cpp
  Bpp/PopGen/BiAlleleMonolocusGenotype.cpp
  Bpp/PopGen/DataSet/AnalyzedLoci.cpp
  Bpp/PopGen/DataSet/AnalyzedSequences.cpp
  Bpp/PopGen/DataSet/DataSet.cpp
  Bpp/PopGen/DataSet/DataSetTools.cpp
  Bpp/PopGen/DataSet/Date.cpp
  Bpp/PopGen/DataSet/Group.cpp
  Bpp/PopGen/DataSet/Individual.cpp
  Bpp/PopGen/DataSet/Io/AbstractIDataSet.cpp
  Bpp/PopGen/DataSet/Io/AbstractODataSet.cpp
  Bpp/PopGen/DataSet/Io/Darwin/DarwinDon.cpp
  Bpp/PopGen/DataSet/Io/Darwin/DarwinVarSingle.cpp
  Bpp/PopGen/DataSet/Io/GeneMapper/GeneMapperCsvExport.cpp
  Bpp/PopGen/DataSet/Io/Genepop/Genepop.cpp
  Bpp/PopGen/DataSet/Io/Genetix/Genetix.cpp
  Bpp/PopGen/DataSet/Io/PopgenlibIO.cpp
  Bpp/PopGen/DataSet/MultiSeqIndividual.cpp
  Bpp/PopGen/GeneralExceptions.cpp
  Bpp/PopGen/LocusInfo.cpp
  Bpp/PopGen/MonoAlleleMonolocusGenotype.cpp
  Bpp/PopGen/MonolocusGenotypeTools.cpp
  Bpp/PopGen/MultiAlleleMonolocusGenotype.cpp
  Bpp/PopGen/MultilocusGenotype.cpp
  Bpp/PopGen/MultilocusGenotypeStatistics.cpp
  Bpp/PopGen/PolymorphismMultiGContainer.cpp
  Bpp/PopGen/PolymorphismMultiGContainerTools.cpp
  Bpp/PopGen/PolymorphismSequenceContainer.cpp
  Bpp/PopGen/PolymorphismSequenceContainerTools.cpp
  Bpp/PopGen/SequenceStatistics.cpp
  )

# Build the static lib
add_library (${PROJECT_NAME}-static STATIC ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-static PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )
set_target_properties (${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries (${PROJECT_NAME}-static ${BPP_LIBS_STATIC})

# Build the shared lib
add_library (${PROJECT_NAME}-shared SHARED ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-shared PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )
set_target_properties (${PROJECT_NAME}-shared
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}
  MACOSX_RPATH 1
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
  )
target_link_libraries (${PROJECT_NAME}-shared ${BPP_LIBS_SHARED})

# Install libs and headers
install (
  TARGETS ${PROJECT_NAME}-static ${PROJECT_NAME}-shared
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX}
  )
install (
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Bpp
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
  )
# Install cmake file representing targets
install (EXPORT ${PROJECT_NAME}-targets DESTINATION ${cmake-package-location})
