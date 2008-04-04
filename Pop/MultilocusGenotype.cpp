//
// File MultilocusGenotype.cpp
// Author : Sylvain Gaillard <sylvain.gaillard@angers.inra.fr>
// Last modification : April 4, 2008
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

#include "MultilocusGenotype.h"

using namespace bpp;

//** Class constructor: *******************************************************/

MultilocusGenotype::MultilocusGenotype(unsigned int loci_number) throw (BadIntegerException)
{
  if (loci_number < 1)
    throw BadIntegerException("MultilocusGenotype::MultilocusGenotype: loci_number must be > 0.", loci_number);

  // Set the _loci size to the right number of loci
  _loci.resize(loci_number);

  // Set all the _loci pointers to NULL
  for (unsigned int i=0 ; i<loci_number ; i++)
    _loci[i] = NULL;
}

MultilocusGenotype::MultilocusGenotype(const MultilocusGenotype & genotype)
{
  for(unsigned int i = 0; i < genotype.size(); i++)
  {
    if (! genotype.isMonolocusGenotypeMissing(i))
      _loci.push_back(dynamic_cast<MonolocusGenotype *>((genotype.getMonolocusGenotype(i))->clone()));
    else
      _loci.push_back(NULL);
  }
}

//** Class destructor: *******************************************************/

MultilocusGenotype::~MultilocusGenotype()
{
  for(unsigned int i = 0; i < _loci.size(); i++)
    delete _loci[i];
  _loci.clear();
}

//** Other methodes: *********************************************************/

void MultilocusGenotype::setMonolocusGenotype(unsigned int locus_position,
    const MonolocusGenotype & monogen) throw (IndexOutOfBoundsException)
{
  if (locus_position < _loci.size())
    _loci[locus_position] = dynamic_cast<MonolocusGenotype *>(monogen.clone());
  else
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locus_position out of bounds.",
        locus_position, 0, _loci.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleKey(unsigned int locus_position,
    const vector<unsigned int> allele_keys) throw (Exception)
{
  if (allele_keys.size() < 1)
    throw Exception("MultilocusGenotype::setMonolocusGenotypeByAlleleKey: no key in allele_keys.");

  if (locus_position < _loci.size())
  {
    setMonolocusGenotype(locus_position, * MonolocusGenotypeTools::buildMonolocusGenotypeByAlleleKey(allele_keys));
  }
  else
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locus_position out of bounds.",
        locus_position, 0, _loci.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleId(unsigned int locus_position,
    const vector<string> allele_id, const LocusInfo & locus_info) throw (Exception) {
  vector<unsigned int> allele_keys;
  for (unsigned int i = 0 ; i < allele_id.size() ; i++) {
    try {
      allele_keys.push_back(locus_info.getAlleleInfoKey(allele_id[i]));
    }
    catch (AlleleNotFoundException & anfe) {
      throw AlleleNotFoundException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
    }
  }
  try {
    setMonolocusGenotypeByAlleleKey(locus_position, allele_keys);
  }
  catch (IndexOutOfBoundsException & ioobe) {
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

void MultilocusGenotype::setMonolocusGenotypeAsMissing(unsigned int locus_position) throw (IndexOutOfBoundsException)
{
  if (locus_position >= _loci.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeAsMissing: locus_position out of bounds.", locus_position, 0, _loci.size());
  if (_loci[locus_position] != NULL)
    delete _loci[locus_position];
  _loci[locus_position] = NULL;
}

bool MultilocusGenotype::isMonolocusGenotypeMissing(unsigned int locus_position) const throw (IndexOutOfBoundsException)
{
  if (locus_position >= _loci.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::isMonolocusGenotypeMissing: locus_position out of bounds.", locus_position, 0, _loci.size());
  return _loci[locus_position] == NULL;
}

const MonolocusGenotype * MultilocusGenotype::getMonolocusGenotype(unsigned int locus_position) const throw (IndexOutOfBoundsException)
{
  if (locus_position >= _loci.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::getMonolocusGenotype: locus_position out of bounds", locus_position, 0, _loci.size());
  return _loci[locus_position];
}

unsigned int MultilocusGenotype::size() const
{
  return _loci.size();
}

unsigned int MultilocusGenotype::countNonMissingLoci() const
{
  unsigned int count = 0;
  for (unsigned int i = 0 ; i < _loci.size() ; i++)
    if (_loci[i] != NULL) count++;
  return count;
}

unsigned int MultilocusGenotype::countHomozygousLoci() const
{
  unsigned int count = 0;
  for (unsigned int i = 0 ; i < _loci.size() ; i++) {
    try {
      if (dynamic_cast<BiAlleleMonolocusGenotype *>(_loci[i])->isHomozygous()) count++;
    }
    catch (...) {}
  }
  return count;
}

unsigned int MultilocusGenotype::countHeterozygousLoci() const
{
  unsigned int count = 0;
  for(unsigned int i = 0; i < _loci.size(); i++)
  {
    try {
      if (!(dynamic_cast<BiAlleleMonolocusGenotype *>(_loci[i])->isHomozygous())) count++;
    }
    catch (...) {}
  }
  return count;
}

