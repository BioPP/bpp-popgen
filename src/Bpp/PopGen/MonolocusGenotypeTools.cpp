// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From Pop
#include "MonoAlleleMonolocusGenotype.h"
#include "BiAlleleMonolocusGenotype.h"
#include "MultiAlleleMonolocusGenotype.h"

#include "MonolocusGenotypeTools.h"

using namespace bpp;
using namespace std;

std::unique_ptr<MonolocusGenotypeInterface> MonolocusGenotypeTools::buildMonolocusGenotypeByAlleleKey(const std::vector<size_t> alleleKeys)
{
  if (alleleKeys.size() < 1)
    throw Exception("MonolocusGenotypeTools::buildMonolocusGenotypeByAlleleKey: no key in allele_keys.");

  if (alleleKeys.size() == 1)
    return unique_ptr<MonolocusGenotypeInterface>(new MonoAlleleMonolocusGenotype(alleleKeys));
  if (alleleKeys.size() == 2)
    return unique_ptr<MonolocusGenotypeInterface>(new BiAlleleMonolocusGenotype(alleleKeys));
  // for all other cases (allele_keys.size() > 2)
  return unique_ptr<MonolocusGenotypeInterface>(new MultiAlleleMonolocusGenotype(alleleKeys));
}

