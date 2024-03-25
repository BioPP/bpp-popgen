// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _BPP_POPGEN_GENEMAPPERCSVEXPORT_H_
#define _BPP_POPGEN_GENEMAPPERCSVEXPORT_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/DataTable.h>

// From local Pop
#include "../AbstractIDataSet.h"
#include "../../../BasicAlleleInfo.h"
#include "../../../MultiAlleleMonolocusGenotype.h"

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
  // bool IndependentAlleles_; //jdutheilon 19/09/14: this does not seem to be used anywhere!

public:
  // Constructor and destructor
  // GeneMapperCsvExport(bool ia = false);
  GeneMapperCsvExport() {}
  ~GeneMapperCsvExport();

  // public:
  /**
   * @brief Set if allels are considered as independent markers.
   *
   */
  // SetAllelsAsIndependent(bool flag);

public:
  /**
   * @name The IDataSet interface.
   * @{
   */
  void read(std::istream& is, DataSet& data_set);
  void read(const std::string& path, DataSet& data_set);
  DataSet* read(std::istream& is);
  DataSet* read(const std::string& path);
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
