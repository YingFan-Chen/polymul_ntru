
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "params.h"
#include "poly_rq_mul.h"

#include "hal.h"

char outs[128];
char *out_ptr;
unsigned long long oldcount, newcount;

int main(void){

    hal_setup(CLOCK_BENCHMARK);

    hal_send_str("\n============ IGNORE OUTPUT BEFORE THIS LINE ============\n");

    poly res, a, b;

    hal_send_str("Toom-4, Karatsuba start");

    oldcount = hal_get_time();
    poly_Rq_mul(&res, &a, &b);
    newcount = hal_get_time();

    sprintf(outs, "poly_Rq_mul: %llu cycles", newcount - oldcount);
    hal_send_str(outs);



}
