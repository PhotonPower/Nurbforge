#pragma once
#include <stdexcept>
#include <string>

namespace nf {

/// Thrown for invalid geometric input (bad knot vectors, degenerate axes, ...).
class GeometryError : public std::runtime_error {
 public:
  explicit GeometryError(const std::string& what) : std::runtime_error(what) {}
};

/// Thrown for features that are planned but not yet implemented.
class NotImplementedError : public std::logic_error {
 public:
  explicit NotImplementedError(const std::string& what) : std::logic_error(what) {}
};

}  // namespace nf
