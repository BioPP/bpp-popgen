// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "BasicAlleleInfo.h"

using namespace bpp;

// ** Class constructor: *******************************************************/

BasicAlleleInfo::BasicAlleleInfo(const std::string& id) : id_(id) {}

BasicAlleleInfo::BasicAlleleInfo(const BasicAlleleInfo& allele) : id_(allele.getId()) {}

// ** Class destructor: *******************************************************/

BasicAlleleInfo::~BasicAlleleInfo() {}

// ** Other methodes: *********************************************************/

BasicAlleleInfo& BasicAlleleInfo::operator=(const BasicAlleleInfo& allele)
{
  id_ = allele.getId();
  return *this;
}

bool BasicAlleleInfo::operator==(const BasicAlleleInfo& allele) const
{
  return id_ == allele.getId();
}

bool BasicAlleleInfo::operator!=(const BasicAlleleInfo& allele) const
{
  return !(id_ == allele.getId());
}

void BasicAlleleInfo::setId(const std::string& allele_id)
{
  id_ = allele_id;
}

const std::string& BasicAlleleInfo::getId() const
{
  return id_;
}
