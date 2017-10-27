#include "extbingcd.h"

void initialize_parameters(
  mpz_t a, mpz_t b, mpz_t alpha, mpz_t beta, 
  mpz_t u, mpz_t v, mpz_t s, mpz_t t, mpz_t r
){
  mpz_init(a);        //reserves memory from a
  mpz_init(b);
  mpz_init(alpha);
  mpz_init(beta);
  mpz_init(u);
  mpz_init(v);
  mpz_init(s);
  mpz_init(t);
  mpz_init(r);
  mpz_set_si(u, 1);   //u = 1
  mpz_set_si(v, 0);   //v = 0
  mpz_set_si(s, 0);   //s = 0
  mpz_set_si(t, 1);   //t = 1
  mpz_set_si(r, 0);   //r = 0
}

void save_result(gcd_t* result, mpz_t a, mpz_t r, mpz_t s, mpz_t t){
  mpz_init(result->gcd);
  mpz_init(result->x);
  mpz_init(result->y);
  mpz_mul_2exp(result->gcd, a, mpz_get_ui(r)); // a << r
  mpz_set(result->x, s);  //s
  mpz_set(result->y, t);  //t
}

void free_result(gcd_t* result){
  //frees memory
  mpz_clear(result->gcd);
  mpz_clear(result->x);
  mpz_clear(result->y);
}

void free_parameters(
  mpz_t a, mpz_t b, mpz_t alpha, mpz_t beta,
  mpz_t u, mpz_t v, mpz_t s, mpz_t t, mpz_t r
){
  //frees memory
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(alpha);
  mpz_clear(beta);
  mpz_clear(u);
  mpz_clear(v);
  mpz_clear(s);
  mpz_clear(t);
  mpz_clear(r);
}

void extended_binary_gcd(gcd_t* result, const mpz_t A, const mpz_t B){
  mpz_t a, b, alpha, beta, u, v, s, t, r;
  initialize_parameters(a, b, alpha, beta, u, v, s, t, r);
  mpz_set(a, A); //moves A to a
  mpz_set(b, B); //moves B to b
  while (mpz_even_p(a) && mpz_even_p(b)){
    mpz_tdiv_q_2exp(a, a, 1); // a = a >> 1
    mpz_tdiv_q_2exp(b, b, 1); // b = b >> 1
    mpz_add_ui(r, r, 1);      // r++
  }
  mpz_set(alpha, a); //alpha = a
  mpz_set(beta, b);  //beta  = b
  while (mpz_even_p(a)){ //a is even
    mpz_tdiv_q_2exp(a, a, 1); // a = a >> 1
    if (mpz_even_p(u) && mpz_even_p(v)){
      mpz_tdiv_q_2exp(u, u, 1); //u = u >> 1
      mpz_tdiv_q_2exp(v, v, 1); //v = v >> 1
    }
    else{
      mpz_add(u, u, beta);      //u = u + beta
      mpz_tdiv_q_2exp(u, u, 1); //u = u >> 1
      mpz_sub(v, v, alpha);     //v = v - alpha
      mpz_tdiv_q_2exp(v, v, 1); //v = v >> 1
    }
  }
  while(mpz_cmp(a, b)){   //a != b
    if (mpz_even_p(b)){   //b is even
      mpz_tdiv_q_2exp(b, b, 1); //b = b >> 1
      if (mpz_even_p(s) && mpz_even_p(t)){
        mpz_tdiv_q_2exp(s, s, 1); //s = s >> 1
        mpz_tdiv_q_2exp(t, t, 1); //t = t >> 1
      }
      else{
        mpz_add(s, s, beta);      //s = s + beta
        mpz_tdiv_q_2exp(s, s, 1); //s = s >> 1
        mpz_sub(t, t, alpha);     //t = t - alpha
        mpz_tdiv_q_2exp(t, t, 1); //t = t >> 1
      }
    }
    else if (mpz_cmp(b, a) < 0){ // b < a
      mpz_set(a, b); //a = b
      mpz_set(b, a); //b = a
      mpz_set(u, s); //u = s
      mpz_set(v, t); //v = t
      mpz_set(s, u); //s = u
      mpz_set(t, v); //t = v
    }
    else{
      mpz_sub(b, b, a); //b = b - a
      mpz_sub(s, s, u); //s = s - u
      mpz_sub(t, t, v); //t = t - v
    }
  }
  save_result(result, a, r, s, t);
  free_parameters(a, b, alpha, beta, u, v, s, t, r);
}

void least_common_multiple(mpz_t result, const mpz_t a, const mpz_t b){
  gcd_t gcd;
  extended_binary_gcd(&gcd, a, b);
  mpz_mul(result, a, b);
  mpz_div(result, result, gcd.gcd);
  free_result(&gcd);
}

//(gcd != 1) ? null : (x%b + b) % b
void mod_inverse(mpz_t result, const mpz_t a, const mpz_t b){
  gcd_t gcd;
  extended_binary_gcd(&gcd, a, b);
  if(!mpz_cmp_ui(gcd.gcd, 1)){
    mpz_mod(result, gcd.x, b);
    mpz_add(result, result, b);
    mpz_mod(result, result, b);
  }
  else
    gmp_printf("Inverse does not exit\n");
  free_result(&gcd);
}