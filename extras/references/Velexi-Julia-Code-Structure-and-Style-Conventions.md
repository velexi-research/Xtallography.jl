Velexi Julia Code Structure and Style Conventions
=================================================

--------------------------------------------------------------------------------------------
## Code Organization

* __Modules__. Define modules within a file with the same name as the module (with a ".jl"
  suffix) that includes the source files that define types and methods for the module.

* __Types__. Use a separate file to define each type and the methods that are "tightly
  associated" with it (e.g., the method's only non-primitive argument is the type, the
  method modifies the data contained in the type).

  For types with a small number of methods, it makes sense to group multiple types and
  their methods into a single file.

* __General Methods__. Organized general methods (i.e., methods not tightly associated
  with a single type or group of types) in a logical manner into one or more files that
  are not "overly large".

--------------------------------------------------------------------------------------------
## Docstrings

### Format

* Function signature

  * Include argument and return types in signatures (using `::` operator)

  * Follow the function signature with a blank line.

* Description of what the function does.

* (Optional) Notes

* (Optional) Examples

* (Optional) Implementation

### Conventions

* Always include a docstring for (1) empty generic functions used to define interfaces and
  (2) methods with no argument type specialization.

* For methods that specialize generic functions/methods, only include a separate docstring
  if the specialized method behaves differently than described in the docstring of the
  generic function/method.

  * _Note_: This convention also applies to methods that specialize functions defined in
    external modules (e.g., Base).

--------------------------------------------------------------------------------------------
## Abstract Types

* Include a docstring for the abstract type that provides a description of the type's
  purpose.

* Define the interface for the type.

  * __Attributes__: functions that return information about an instance of the type

  * __Operations/Computations__: functions that perform operations, computations, etc.
    that "dependent" on the value of an instance of the type

* Use empty generic functions to define the signatures (i.e., arguments and return values)
  for all interface functions. Use the `::` operator to indicate argument and
  return value types (when appropriate).

--------------------------------------------------------------------------------------------
