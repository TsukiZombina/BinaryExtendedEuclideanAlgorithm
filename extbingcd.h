//https://gmplib.org/manual/Concept-Index.html#Concept-Index
//https://gmplib.org/manual/Integer-Functions.html
	//#Integer-Functions
//mpz_gcdext(result, s, t, x, y);
#include <gmp.h>

typedef struct{
  mpz_t x, y, gcd;
} gcd_t;

// // @brief: This function initializes u, v, s, t, and r
// // with values 
void initialize_parameters(mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t);

void free_parameters(mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t);

// @brief: This function computes the greatest common divisor of two 
// arbitrary-precision numbers as well the coefficients of Bezout's identity
// The result put into a three element struct called gcd_result. This struct
// is make up of three arbitrary-precision numbers named x, y, and gcd.
void extended_binary_gcd(gcd_t*, const mpz_t, const mpz_t);

void save_result(gcd_t*, mpz_t, mpz_t, mpz_t, mpz_t);

void free_result(gcd_t*);

void least_common_multiple(mpz_t, const mpz_t, const mpz_t);

//(gcd != 1) ? null : (x%b + b) % b
void mod_inverse(mpz_t, const mpz_t, const mpz_t);