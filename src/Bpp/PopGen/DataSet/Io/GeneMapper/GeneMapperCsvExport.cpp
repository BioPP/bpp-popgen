//
// File:    GeneMapperCsvExport.cpp
// Author:  Sylvain Gaillard
// Created: April 2, 2008
//

/*
   Copyright or © or Copr. Bio++ Development Team, (April 2, 2008)

   This software is a computer program whose purpose is to provide classes
   for population genetics analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "GeneMapperCsvExport.h"

using namespace bpp;
using namespace std;

const std::string GeneMapperCsvExport::SAMPLE_FILE_H = "Sample File";
const std::string GeneMapperCsvExport::SAMPLE_NAME_H = "Sample Name";
const std::string GeneMapperCsvExport::PANEL_H = "Panel";
const std::string GeneMapperCsvExport::MARKER_H = "Marker";
const std::string GeneMapperCsvExport::DYE_H = "Dye";
const std::string GeneMapperCsvExport::ALLELE_H = "Allele ";
const std::string GeneMapperCsvExport::SIZE_H = "Size ";
const std::string GeneMapperCsvExport::HEIGHT_H = "Height ";
const std::string GeneMapperCsvExport::PEAK_AREA_H = "Peak Area ";
const std::string GeneMapperCsvExport::DAC_H = "DAC";
const std::string GeneMapperCsvExport::AN_H = "AN";

// GeneMapperCsvExport::GeneMapperCsvExport(bool ia) : IndependentAlleles_(ia) {}

GeneMapperCsvExport::~GeneMapperCsvExport() {}

void GeneMapperCsvExport::read(std::istream& is, DataSet& data_set)
{
  if (!is)
    throw IOException("GeneMapperCsvExport::read: fail to open stream.");

  /*
   * Feed a DataTable with the data
   */
  auto dtp = DataTable::read(is, "\t", true, -1);
  DataTable& dt = *dtp;

  /*
   * Fixe the individuals' name if there is duplicate in the file
   */
  vector<string> ind_names;
  vector<string> markers;
  try
  {
    ind_names = dt.getColumn(SAMPLE_NAME_H);
    markers = dt.getColumn(MARKER_H);
  }
  catch (Exception& e)
  {
    throw e;
  }
  map<string, int> indname_marker;
  for (size_t i = 0; i < dt.getNumberOfRows(); i++)
  {
    string test_lab = dt(i, SAMPLE_NAME_H) + dt(i, MARKER_H);
    if (indname_marker.find(test_lab) != indname_marker.end())
    {
      string new_lab = dt(i, SAMPLE_NAME_H) + "_" + TextTools::toString(indname_marker[test_lab] + 1);
      dt (i, SAMPLE_NAME_H) = new_lab;
    }
    indname_marker[test_lab]++;
  }
  ind_names = dt.getColumn(SAMPLE_NAME_H);

  map<string, size_t> ind_count = VectorTools::countValues(ind_names);
  ind_names = VectorTools::unique(ind_names);
  markers = VectorTools::unique(markers);
  size_t loc_nbr = markers.size();

  /*
   * Loci number
   */
  data_set.initAnalyzedLoci(loc_nbr);

  /*
   * Group of individuals
   */
  data_set.addEmptyGroup(0);
  for (unsigned int i = 0; i < ind_names.size(); i++)
  {
    Individual ind(ind_names[i]);
    data_set.addIndividualToGroup(data_set.getGroupPosition(0), ind);
  }

  /*
   * Loci data
   */
  AnalyzedLoci al(markers.size());
  vector<string> col_names = dt.getColumnNames();

  // Finds columns containing allele data
  vector<size_t> alleles_cols;
  for (size_t i = 0; i < col_names.size(); i++)
  {
    if (TextTools::startsWith(col_names[i], ALLELE_H))
      alleles_cols.push_back(i);
  }
  // Set LocusInfo
  vector<vector<size_t> > alleles_pos;
  for (size_t i = 0; i < markers.size(); i++)
  {
    al.setLocusInfo(i, LocusInfo(markers[i], LocusInfo::UNKNOWN));
  }
  std::map< std::string, std::set< std::string > > markerAlleles;
  for (size_t i = 0; i < dt.getNumberOfRows(); ++i)
  {
    for (size_t j = 0; j < alleles_cols.size(); ++j)
    {
      if (dt(i, alleles_cols[j]) != "")
      {
        markerAlleles[dt(i, MARKER_H)].insert(dt(i, alleles_cols[j]));
      }
    }
  }
  for (std::map< std::string, std::set< std::string > >::iterator itm = markerAlleles.begin(); itm != markerAlleles.end(); itm++)
  {
    std::set< std::string >& s = itm->second;
    for (std::set< std::string >::iterator its = s.begin(); its != s.end(); its++)
    {
      al.addAlleleInfoByLocusName(itm->first, BasicAlleleInfo(*its));
    }
  }
  data_set.setAnalyzedLoci(al);

  /*
   * Individuals informations
   */
  size_t ind_col_index = VectorTools::which(dt.getColumnNames(), SAMPLE_NAME_H);
  size_t mark_col_index = VectorTools::which(dt.getColumnNames(), MARKER_H);
  for (size_t i = 0; i < dt.getNumberOfRows(); i++)
  {
    vector<size_t> alleles;
    for (size_t j = 0; j < alleles_cols.size(); j++)
    {
      if (!TextTools::isEmpty(dt(i, alleles_cols[j])))
      {
        unsigned int num = (data_set.getLocusInfoByName(dt(i, mark_col_index))).getAlleleInfoKey(dt(i, alleles_cols[j]));
        alleles.push_back(num);
      }
    }
    alleles = VectorTools::unique(alleles);
    MultiAlleleMonolocusGenotype ma(alleles);
    if (!data_set.getIndividualByIdFromGroup(0, dt(i, ind_col_index))->hasGenotype())
      data_set.initIndividualGenotypeInGroup(0, data_set.getIndividualPositionInGroup(0, dt(i, ind_col_index)));
    if (alleles.size())
      data_set.setIndividualMonolocusGenotypeInGroup(0, data_set.getIndividualPositionInGroup(0, dt(i, ind_col_index)), data_set.getAnalyzedLoci()->getLocusInfoPosition(dt(i, mark_col_index)), ma);
  }
}

void GeneMapperCsvExport::read(const std::string& path, DataSet& data_set)
{
  AbstractIDataSet::read(path, data_set);
}

DataSet* GeneMapperCsvExport::read(std::istream& is)
{
  return AbstractIDataSet::read(is);
}

DataSet* GeneMapperCsvExport::read(const std::string& path)
{
  return AbstractIDataSet::read(path);
}

// --- GeneMapperCsvExport::Record ---
GeneMapperCsvExport::Record::Record(const std::string& row) : sampleFile_(),
  sampleName_(),
  panel_(),
  markerName_(),
  dye_(),
  alleles_(),
  dac_(),
  an_(0.)
{
  StringTokenizer st(row, "\t", true, false);
  /*
     if (st.numberOfRemainingTokens() != 7 + 4 * alleleNumber) {
     throw Exception("GeneMapperCsvExport::Record::Record: bad number of allele");
     }
   */
  size_t itemNum = st.numberOfRemainingTokens();
  size_t alleleNum = (itemNum - 7) / 4;
  sampleFile_ = st.getToken(0);
  sampleName_ = st.getToken(1);
  panel_ = st.getToken(2);
  markerName_ = st.getToken(3);
  dye_ = st.getToken(4);
  dac_ = st.getToken(itemNum - 2);
  an_ = TextTools::toDouble(st.getToken(itemNum - 1));
  for (unsigned int i = 0; i < alleleNum; ++i)
  {
    GeneMapperCsvExport::Allele al(
      st.getToken(5 + i),
      TextTools::toDouble(st.getToken(5 + alleleNum + i)),
      TextTools::to<unsigned int>(st.getToken(5 + (2 * alleleNum) + i)),
      TextTools::toDouble(st.getToken(5 + (3 * alleleNum) + i))
      );
    alleles_.push_back(al);
  }
}
