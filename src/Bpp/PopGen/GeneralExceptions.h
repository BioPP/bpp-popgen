// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _GENERALEXCEPTIONS_H_
#define _GENERALEXCEPTIONS_H_

// From STL
#include <string>

#include <Bpp/Exceptions.h>

namespace bpp
{
// ****************************************************************************
//
/**
 * @brief The BadIdentifierException class.
 *
 * This exception is used when an identifier is not found.
 * The identifier can be either a string or an integer but its
 * value is stored as a string.
 *
 * @author Sylvain Gaillard
 */
class BadIdentifierException :
  public Exception
{
public:
  // Class constructor
  /**
   * @brief Build the exception with a numerical identifier.
   */
  BadIdentifierException(const char* text, size_t id);
  /**
   * @brief Build the exception with a numerical identifier.
   */
  BadIdentifierException(const std::string& text, size_t id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  BadIdentifierException(const char* text, const std::string& id);
  /**
   * @brief Build the exception with a textual identifier.
   */
  BadIdentifierException(const std::string& text, const std::string& id);

  // Class destructor
  ~BadIdentifierException() throw ();

public:
  /**
   * @brief Return the value of the identifier as a string.
   */
  virtual const std::string getIdentifier() const;

protected:
  const std::string id_;
};

// *****************************************************************************

/**
 * @brief The LocusNotFoundException class.
 */
class LocusNotFoundException :
  public BadIdentifierException
{
public:
  // Class constructor
  /**
   * @brief Build the exception with a numerical identifier.
   */
  LocusNotFoundException(const char* text, size_t id);

  /**
   * @brief Build the exception with a numerical identifier.
   */
  LocusNotFoundException(const std::string& text, size_t id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  LocusNotFoundException(const char* text, const std::string& id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  LocusNotFoundException(const std::string& text, const std::string& id);

  // Class destructor
  ~LocusNotFoundException() throw ();

public:
  /**
   * @brief Return the value of the identifier as a string.
   */
  virtual const std::string getIdentifier() const;
};

// *****************************************************************************

/**
 * @brief The AlleleNotFoundException class.
 */
class AlleleNotFoundException :
  public BadIdentifierException
{
public:
  // Class constructor
  /**
   * @brief Build the exception with a numerical identifier.
   */
  AlleleNotFoundException(const char* text, size_t id);

  /**
   * @brief Build the exception with a numerical identifier.
   */
  AlleleNotFoundException(const std::string& text, size_t id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  AlleleNotFoundException(const char* text, const std::string& id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  AlleleNotFoundException(const std::string& text, const std::string& id);

  // Class destructor
  ~AlleleNotFoundException() throw ();

public:
  /**
   * @brief Return the value of the identifier as a string.
   */
  virtual const std::string getIdentifier() const;
};

// *****************************************************************************

/**
 * @brief The LocalityNotFoundException class.
 */
class LocalityNotFoundException :
  public BadIdentifierException
{
public:
  // Class constructor
  /**
   * @brief Build the exception with a numerical identifier.
   */
  LocalityNotFoundException(const char* text, size_t id);

  /**
   * @brief Build the exception with a numerical identifier.
   */
  LocalityNotFoundException(const std::string& text, size_t id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  LocalityNotFoundException(const char* text, const std::string& id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  LocalityNotFoundException(const std::string& text, const std::string& id);

  // Class destructor
  ~LocalityNotFoundException() throw ();

public:
  /**
   * @brief Return the value of the identifier as a string.
   */
  virtual const std::string getIdentifier() const;
};

// *****************************************************************************

/**
 * @brief The IndividualNotFoundException class.
 */
class IndividualNotFoundException :
  public BadIdentifierException
{
public:
  // Class constructor
  /**
   * @brief Build the exception with a numerical identifier.
   */
  IndividualNotFoundException(const char* text, size_t id);

  /**
   * @brief Build the exception with a numerical identifier.
   */
  IndividualNotFoundException(const std::string& text, size_t id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  IndividualNotFoundException(const char* text, const std::string& id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  IndividualNotFoundException(const std::string& text, const std::string& id);

  // Class destructor
  ~IndividualNotFoundException() throw ();

public:
  /**
   * @brief Return the value of the identifier as a string.
   */
  virtual const std::string getIdentifier() const;
};

// *****************************************************************************

/**
 * @brief The GroupNotFoundException class.
 */
class GroupNotFoundException :
  public BadIdentifierException
{
public:
  // Class constructor
  /**
   * @brief Build the exception with a numerical identifier.
   */
  GroupNotFoundException(const char* text, size_t id);

  /**
   * @brief Build the exception with a numerical identifier.
   */
  GroupNotFoundException(const std::string& text, size_t id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  GroupNotFoundException(const char* text, const std::string& id);

  /**
   * @brief Build the exception with a textual identifier.
   */
  GroupNotFoundException(const std::string& text, const std::string& id);

  // Class destructor
  ~GroupNotFoundException() throw ();

public:
  /**
   * @brief Return the value of the identifier as a string.
   */
  virtual const std::string getIdentifier() const;
};
} // end of namespace bpp;

#endif// _GENERALEXCEPTIONS_H_
