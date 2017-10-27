//Please, pass values x and y from the command line as this:
//     ./gcd 11 13
// where 11 is a and 13 is b
#include <stdio.h>

#include "extbingcd.h"

int main(int argc, char* argv[]) {
  mpz_t a, b, lcm, inv;
  gcd_t result;
  mpz_init(lcm);
  mpz_init(inv);
  mpz_init_set_str(a, argv[1], 10);
  mpz_init_set_str(b, argv[2], 10);
  extended_binary_gcd(&result, a, b);
  least_common_multiple(lcm, a, b);
  mod_inverse(inv, a, b);
  gmp_printf("gcd: %Zd\n", result.gcd);
  gmp_printf("x: %Zd\n", result.x);
  gmp_printf("y: %Zd\n", result.y);
  gmp_printf("lcm: %Zd\n", lcm);
  gmp_printf("inverse: %Zd\n", inv);
  free_result(&result);
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(lcm);
  mpz_clear(inv
  );
  return 0;
}