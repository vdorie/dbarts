#ifndef MISC_TYPES_H
#define MISC_TYPES_H

/* misc/types.h for a standalone build, which cannot run configure: written
   from src/include/misc/types.h.in with configure's one value. */

#include <stdint.h> // type should be replaced with a uintXX_t

#ifndef XINT_TYPE
/* The integer type used to map the floating point x to integer types,
   which limits the number of cut points allowable. configure fixes this
   at uint16_t with no build-time override (--with-xint-size errors out). */
#define XINT_TYPE uint16_t
#endif

typedef XINT_TYPE misc_xint_t;

#endif // MISC_TYPES_H
