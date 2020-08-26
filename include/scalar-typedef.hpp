#ifndef SCALAR_TYPEDEF_H_
#define SCALAR_TYPEDEF_H_

#if defined(USE_BASIC_TYPES) && !defined(USE_DCO_TYPES)
#elif !defined(USE_BASIC_TYPES) && defined(USE_DCO_TYPES)

#include "unsupported/Eigen/dco"

typedef dco::gt1s<double> DCO_GT1S_MODE;
typedef DCO_GT1S_MODE::type gt1s_scalar;

typedef dco::ga1s<double> DCO_GA1S_MODE;
typedef DCO_GA1S_MODE::type ga1s_scalar;

typedef dco::ga1s<gt1s_scalar> DCO_GT2S_GA1S_MODE;
typedef DCO_GT2S_GA1S_MODE::type gt2s_ga1s_scalar;

#ifndef DCO_MAX_ALLOCATION
# error "DCO_MAX_ALLOCATION is not defined"
#endif

#ifndef DCO_MEM_RATIO
# error "DCO_MEM_RATIO is not defined"
#endif

#else
# error "Exactly one of USE_BASIC_TYPES or USE_DCO_TYPES needs to be defined"
#endif

#endif
