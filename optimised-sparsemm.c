#include "utils.h"

void basic_sparsemm(const COO, const COO, COO *);

void optimised_sparsemm(const COO a, const COO b, COO *c)
{
    return basic_sparsemm(a, b, c);
}
