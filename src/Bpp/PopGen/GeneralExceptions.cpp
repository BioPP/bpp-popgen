// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
