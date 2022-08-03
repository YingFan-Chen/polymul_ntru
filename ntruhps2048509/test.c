
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

static void polymul_schoolbook(uint16_t *res, uint16_t *a, uint16_t *b){

    uint16_t buf[2 * POLY_N];

    for(size_t i = 0; i < 2 * POLY_N; i++){
        buf[i] = 0;
    }

    for(size_t i = 0; i < NTRU_N; i++){
        for(size_t j = 0; j < NTRU_N; j++){
            buf[i + j] += a[i] * b[j];
        }
    }

    for(size_t i = 2 * POLY_N - 1; i >= NTRU_N; i--){
        buf[i - NTRU_N] += buf[i];
    }

    for(size_t i = 0; i < NTRU_N; i++){
        res[i] = buf[i];
    }

}

int main(void){

    hal_setup(CLOCK_FAST);

    hal_send_str("\n============ IGNORE OUTPUT BEFORE THIS LINE ============\n");

    uint16_t poly_res[POLY_N], poly_a[POLY_N], poly_b[POLY_N];
    
    poly res, a, b;

    for(size_t i = 0; i < NTRU_N; i++){
        poly_a[i] = rand() % NTRU_Q;
        poly_b[i] = rand() % 3;
        if(poly_b[i] == 2){
            poly_b[i] = 2047;
        }
        a.coeffs[i] = poly_a[i];
        b.coeffs[i] = poly_b[i];
    }
    for(size_t i = NTRU_N; i < POLY_N; i++){
        poly_a[i] = poly_b[i] = a.coeffs[i] = b.coeffs[i] = 0;
    }

    hal_send_str("schoolbook start");
    polymul_schoolbook(poly_res, poly_a, poly_b);

    hal_send_str("Toom-4, Karatsuba start");
    poly_Rq_mul(&res, &a, &b);

    for(size_t i = 0; i < NTRU_N; i++){
        if((res.coeffs[i] & 2047) != (poly_res[i] & 2047)){
            sprintf(outs, "%4zu: %8d, %8d\n", i, (res.coeffs[i] & 2047), (poly_res[i] & 2047));
            hal_send_str(outs);
        }
    }

    hal_send_str("test finished");

}
