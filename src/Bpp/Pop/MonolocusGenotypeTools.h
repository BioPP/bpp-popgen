//
// File MonolocusGenotypeTools.h
// Author : Sylvain Gaillard <sylvain.gaillard@angers.inra.fr>
// Last modification : April 4, 2008
//

/*
   Copyright or © or Copr. CNRS, (April 4, 2008)

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
#ifndef _MonolocusGenotypeTools_h_
#define _MonolocusGenotypeTools_h_

// From STL
#include <vector>
#include <memory>

#include <Bpp/Exceptions.h>

//From Pop
#include "MonolocusGenotype.h"

namespace bpp
{

  /**
   * @brief The MonolocusGenotypeTools static class.
   *
   * This class provides tools for MonolocusGenotype manipulation or creation.
   *
   * @author Sylvain Gaillard
   */
  class MonolocusGenotypeTools
  {
    public:
    /**
     * @brief Build a proper MonolocusGenotype accordig to the number of alleles.
     *
     * Return a MonolocusGenotype build according to the number of allels.
     * If one allele key, send a MonoAlleleMonolocusGenotype,
     * if two allele keys, send a BiAlleleMonolocusGenotype,
     * if more allele keys, send a MultiAlleleMonolocusGenotype.
     *
     * @param allele_keys A vector containing thes allele keys to put in the MonolocusGenotype.
     * @return A MonolocusGenotype according to the number of alleles
     */
    static std::auto_ptr<MonolocusGenotype> buildMonolocusGenotypeByAlleleKey(const std::vector<unsigned int> allele_keys) throw (Exception);
  };

} //end of namespace bpp;

#endif // _MonolocusGenotypeTools_h_

