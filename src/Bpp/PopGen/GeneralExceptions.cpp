//
// File GeneralExceptions.cpp
// Author : Sylvain Gaillard
// Last modification: Thursday July 29 2004
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

#include "GeneralExceptions.h"

#include <Bpp/Text/TextTools.h>

using namespace bpp;
using namespace std;

// ** BadIdentifierException **************************************************/

BadIdentifierException::BadIdentifierException(const char* text,
                                               const size_t id) : Exception("BadIdentifierException: " +
                                                                                  string(text) + "(" + TextTools::toString(id) + ")"),
  id_(TextTools::toString(id)) {}

BadIdentifierException::BadIdentifierException(const std::string& text,
                                               const size_t id) : Exception("BadIdentifierException: " +
                                                                                  text + "(" + TextTools::toString(id) + ")"),
  id_(TextTools::toString(id)) {}

BadIdentifierException::BadIdentifierException(const char* text,
                                               const std::string& id) : Exception("BadIdentifierException: " + string(text) +
                                                                                  "(" + id + ")"),
  id_(id) {}

BadIdentifierException::BadIdentifierException(const std::string& text,
                                               const std::string& id) : Exception("BadIdentifierException: " + text +
                                                                                  "(" + id + ")"),
  id_(id) {}

BadIdentifierException::~BadIdentifierException() throw () {}

const std::string BadIdentifierException::getIdentifier() const
{
  return id_;
}

// ** LocusNotFoundException **************************************************/

LocusNotFoundException::LocusNotFoundException(const char* text,
                                               const size_t id) : BadIdentifierException("LocusNotFoundException: " +
                                                                                               string(text) + "(" + TextTools::toString(id) + ")",
                                                                                               id) {}

LocusNotFoundException::LocusNotFoundException(const std::string& text,
                                               const size_t id) : BadIdentifierException("LocusNotFoundException: " +
                                                                                               text + "(" + TextTools::toString(id) + ")",
                                                                                               id) {}

LocusNotFoundException::LocusNotFoundException(const char* text,
                                               const std::string& id) : BadIdentifierException("LocusNotFoundException: " + string(text) +
                                                                                               "(" + id + ")",
                                                                                               id) {}

LocusNotFoundException::LocusNotFoundException(const std::string& text,
                                               const std::string& id) : BadIdentifierException("LocusNotFoundException: " + text +
                                                                                               "(" + id + ")",
                                                                                               id) {}

LocusNotFoundException::~LocusNotFoundException() throw () {}

const std::string LocusNotFoundException::getIdentifier() const
{
  return BadIdentifierException::getIdentifier();
}

// ** AlleleNotFoundException **************************************************/

AlleleNotFoundException::AlleleNotFoundException(const char* text,
                                                 const size_t id) : BadIdentifierException("AlleleNotFoundException: " +
                                                                                                 string(text) + "(" + TextTools::toString(id) + ")",
                                                                                                 id) {}

AlleleNotFoundException::AlleleNotFoundException(const std::string& text,
                                                 const size_t id) : BadIdentifierException("AlleleNotFoundException: " +
                                                                                                 text + "(" + TextTools::toString(id) + ")",
                                                                                                 id) {}

AlleleNotFoundException::AlleleNotFoundException(const char* text,
                                                 const std::string& id) : BadIdentifierException("AlleleNotFoundException: " + string(text) +
                                                                                                 "(" + id + ")",
                                                                                                 id) {}

AlleleNotFoundException::AlleleNotFoundException(const std::string& text,
                                                 const std::string& id) : BadIdentifierException("AlleleNotFoundException: " + text +
                                                                                                 "(" + id + ")",
                                                                                                 id) {}

AlleleNotFoundException::~AlleleNotFoundException() throw () {}

const std::string AlleleNotFoundException::getIdentifier() const
{
  return BadIdentifierException::getIdentifier();
}

// ** LocalityNotFoundException **************************************************/

LocalityNotFoundException::LocalityNotFoundException(const char* text,
                                                     const size_t id) : BadIdentifierException("LocalityNotFoundException: " +
                                                                                                     string(text) + "(" + TextTools::toString(id) + ")",
                                                                                                     id) {}

LocalityNotFoundException::LocalityNotFoundException(const std::string& text,
                                                     const size_t id) : BadIdentifierException("LocalityNotFoundException: " +
                                                                                                     text + "(" + TextTools::toString(id) + ")",
                                                                                                     id) {}

LocalityNotFoundException::LocalityNotFoundException(const char* text,
                                                     const std::string& id) : BadIdentifierException("LocalityNotFoundException: " + string(text) +
                                                                                                     "(" + id + ")",
                                                                                                     id) {}

LocalityNotFoundException::LocalityNotFoundException(const std::string& text,
                                                     const std::string& id) : BadIdentifierException("LocalityNotFoundException: " + text +
                                                                                                     "(" + id + ")",
                                                                                                     id) {}

LocalityNotFoundException::~LocalityNotFoundException() throw () {}

const std::string LocalityNotFoundException::getIdentifier() const
{
  return BadIdentifierException::getIdentifier();
}

// ** IndividualNotFoundException **************************************************/

IndividualNotFoundException::IndividualNotFoundException(const char* text,
                                                         const size_t id) : BadIdentifierException("IndividualNotFoundException: " +
                                                                                                         string(text) + "(" + TextTools::toString(id) + ")",
                                                                                                         id) {}

IndividualNotFoundException::IndividualNotFoundException(const std::string& text,
                                                         const size_t id) : BadIdentifierException("IndividualNotFoundException: " +
                                                                                                         text + "(" + TextTools::toString(id) + ")",
                                                                                                         id) {}

IndividualNotFoundException::IndividualNotFoundException(const char* text,
                                                         const std::string& id) : BadIdentifierException("IndividualNotFoundException: " + string(text) +
                                                                                                         "(" + id + ")",
                                                                                                         id) {}

IndividualNotFoundException::IndividualNotFoundException(const std::string& text,
                                                         const std::string& id) : BadIdentifierException("IndividualNotFoundException: " + text +
                                                                                                         "(" + id + ")",
                                                                                                         id) {}

IndividualNotFoundException::~IndividualNotFoundException() throw () {}

const std::string IndividualNotFoundException::getIdentifier() const
{
  return BadIdentifierException::getIdentifier();
}

// ** GroupNotFoundException **************************************************/

GroupNotFoundException::GroupNotFoundException(const char* text,
                                               const size_t id) : BadIdentifierException("GroupNotFoundException: " +
                                                                                               string(text) + "(" + TextTools::toString(id) + ")",
                                                                                               id) {}

GroupNotFoundException::GroupNotFoundException(const std::string& text,
                                               const size_t id) : BadIdentifierException("GroupNotFoundException: " +
                                                                                               text + "(" + TextTools::toString(id) + ")",
                                                                                               id) {}

GroupNotFoundException::GroupNotFoundException(const char* text,
                                               const std::string& id) : BadIdentifierException("GroupNotFoundException: " + string(text) +
                                                                                               "(" + id + ")",
                                                                                               id) {}

GroupNotFoundException::GroupNotFoundException(const std::string& text,
                                               const std::string& id) : BadIdentifierException("GroupNotFoundException: " + text +
                                                                                               "(" + id + ")",
                                                                                               id) {}

GroupNotFoundException::~GroupNotFoundException() throw () {}

const std::string GroupNotFoundException::getIdentifier() const
{
  return BadIdentifierException::getIdentifier();
}

