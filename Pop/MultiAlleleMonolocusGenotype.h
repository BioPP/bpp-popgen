//
// File MultiAlleleMonolocusGenotype.h
// Author : Sylvain Gaillard <sylvain.gaillard@angers.inra.fr>
// Last modification : Wednesday March 5 2008
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

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
// Secured inclusion of header's file
#ifndef _MULTIALLELEMONOLOCUSGENOTYPE_H_
#define _MULTIALLELEMONOLOCUSGENOTYPE_H_

// From STL
#include <vector>

using namespace std;

// From Utils
#include <Utils/Exceptions.h>

//From local
#include "MonolocusGenotype.h"

namespace bpp
{

  /**
   * @brief The MultiAlleleMonolocusGenotype class.
   *
   * This class is intended to handle monolocus genotype with many alleles
   * like polyploid loci or loci obtained by trace file without cutoff on
   * peaks or other filter.
   *
   * @author Sylvain Gaillard
   */
  class MultiAlleleMonolocusGenotype:
    public MonolocusGenotype
  {
    public: // Constructors and destructor

      /**
       * @brief Build a monolocus genotype containing many alleles.
       */
      MultiAlleleMonolocusGenotype(vector<unsigned int> allele_index);

      /**
       * @brief Copy constructor.
       */
      MultiAlleleMonolocusGenotype(const MultiAlleleMonolocusGenotype & mmg);

      /**
       * @brief Destroy the MultiAlleleMonolocusGenotype.
       */
      ~MultiAlleleMonolocusGenotype();

    public: // Other methodes
      /**
       * @brief The affectation operator.
       */
      MultiAlleleMonolocusGenotype & operator= (const MultiAlleleMonolocusGenotype & mmg);

      /**
       * @brief The == operator.
       */
      bool operator== (const MultiAlleleMonolocusGenotype & mmg) const;

      /**
       * @brief Test the homozygozity of the locus (i.e. all allele are identical).
       */
      bool isHomozygous() const;

      /**
       * @name The MonolocusGenotype interface:
       *
       * @{
       */
      vector<unsigned int> getAlleleIndex() const;
      /** @} */

      /**
       * @name The Clonable interface:
       *
       * @{
       */
      Clonable * clone() const;
      /** @} **/
    protected:
      vector<unsigned int> _allele_index;
  };

} //end of namespace bpp;

#endif // _MULTIALLELEMONOLOCUSGENOTYPE_H_

