//
// File:    GeneMapperCsvExport.h
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

#ifndef _BPP_POPGEN_GENEMAPPERCSVEXPORT_H_
#define _BPP_POPGEN_GENEMAPPERCSVEXPORT_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/DataTable.h>

// From local Pop
#include "AbstractIDataSet.h"
#include "BasicAlleleInfo.h"
#include "MultiAlleleMonolocusGenotype.h"

namespace bpp
{
/**
 * @brief The GeneMapperCsvExport input format for popgenlib.
 *
 * This input format takes a csv file exported from GeneMapper® (Applied Biosystems).
 *
 * @author Sylvain Gaillard
 */
class GeneMapperCsvExport : public AbstractIDataSet
{
public:
  static const std::string SAMPLE_FILE_H;
  static const std::string SAMPLE_NAME_H;
  static const std::string PANEL_H;
  static const std::string MARKER_H;
  static const std::string DYE_H;
  static const std::string ALLELE_H;
  static const std::string SIZE_H;
  static const std::string HEIGHT_H;
  static const std::string PEAK_AREA_H;
  static const std::string DAC_H;
  static const std::string AN_H;

private:
  //bool IndependentAlleles_; //jdutheilon 19/09/14: this does not seem to be used anywhere!

public:
  // Constructor and destructor
  //GeneMapperCsvExport(bool ia = false);
  GeneMapperCsvExport() {}
  ~GeneMapperCsvExport();

  // public:
  /**
   * @brief Set if allels are concidered as independent markers.
   *
   */
  // SetAllelsAsIndependent(bool flag);

public:
  /**
   * @name The IDataSet interface.
   * @{
   */
  void read(std::istream& is, DataSet& data_set) throw (Exception);
  void read(const std::string& path, DataSet& data_set) throw (Exception);
  DataSet* read(std::istream& is) throw (Exception);
  DataSet* read(const std::string& path) throw (Exception);
  /**
   * @}
   */

  /**
   * @name The IOFormat interface
   * @{
   */
  virtual const std::string getFormatName() const
  {
    return "GeneMapper® cvs export";
  }
  virtual const std::string getFormatDescription() const
  {
    return "GeneMapper® is a flexible genotyping software package that provides DNA sizing and quality allele calls for all Applied Biosystems electrophoresis-based genotyping systems.";
  }
  /**
   * @}
   */

  /**
   * @brief Store data for one allele
   */
  class Allele
  {
private:
    std::string name_;
    double size_;
    unsigned int height_;
    double peakArea_;

public:
    Allele(const std::string& name, double size, unsigned int height, double peakArea) : name_(name),
      size_(size),
      height_(height),
      peakArea_(peakArea) {}

    const std::string& getName() const
    {
      return name_;
    }
    const double& getSize() const
    {
      return size_;
    }
    const unsigned int& getHeight() const
    {
      return height_;
    }
    const double& getPeakArea() const
    {
      return peakArea_;
    }
  };

  /**
   * @brief Store one line of the GeneMapper file
   */
  class Record
  {
private:
    std::string sampleFile_;
    std::string sampleName_;
    std::string panel_;
    std::string markerName_;
    std::string dye_;
    std::vector< GeneMapperCsvExport::Allele > alleles_;
    std::string dac_;
    double an_;

public:
    /**
     * @brief Constructor
     *
     * @param row One row of the file as a std::string
     */
    Record(const std::string& row);

    const std::string& getSampleFileName() const
    {
      return sampleFile_;
    }
    const std::string& getSampleName() const
    {
      return sampleName_;
    }
    const std::string& getPanel() const
    {
      return panel_;
    }
    const std::string& getMarkerName() const
    {
      return markerName_;
    }
    const std::string& getDye() const
    {
      return dye_;
    }
    const size_t getNumberOfAllele() const
    {
      return alleles_.size();
    }
    const GeneMapperCsvExport::Allele& getAllele(size_t allelePos) const
    {
      return alleles_[allelePos];
    }
  };
};
} // end of namespace bpp;

#endif // _BPP_POPGEN_GENEMAPPERCSVEXPORT_H_

