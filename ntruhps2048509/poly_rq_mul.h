#ifndef POLY_RQ_MUL_H
#define POLY_RQ_MUL_H

#include "params.h"

typedef struct{
  uint16_t coeffs[POLY_N];
} poly;

void poly_Rq_mul(poly *r, const poly *a, const poly *b);

#endif

