/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                integrer.c                                  *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
 *                                                                            *
 * This library is free software; you can redistribute it and/or modify it    *
 * under the terms of the GNU Lesser General Public License as published by   *
 * the Free Software Foundation; either version 2.1 of the License, or (at    *
 * your option) any later version.                                            *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this library; if not, write to the Free Software Foundation,    *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA         *
 *                                                                            *
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "pip.h"

/*  The routines in this file are used to build a Gomory cut from
    a non-integral row of the problem tableau                             */

extern long int cross_product, limit;
extern int verbose;
extern int deepest_cut;
extern FILE * dump;
char compose[256];


// From osl_int
long long int piplib_llgcd(long long int const a, long long int const b) {
  return (b ? piplib_llgcd(b, a % b) : a);
}
long long int piplib_llgcd_llabs(long long int const a, long long int const b) {
  return llabs(piplib_llgcd(a, b));
}
size_t piplib_lllog2(long long int x) {
  size_t n = 0;

  x = llabs(x);

  while (x) { x >>= 1; ++n; }

  return ((n == 0) ? 1 : n);
}
long long int piplib_llmod(long long int const a, long long int const b) {
  long long mod = a % b;
   if (mod < 0) { mod += llabs(b); }
  return mod;
}
long long int piplib_ll_floor_div_q(long long int const a,
                                    long long int const b) {
  long long int q = a / b;
  if (q < 0) { if (a % b != 0) --q; }
  else if (q == 0) {
    if ((a > 0 && b < 0) || (a < 0 && b > 0)) { --q; }
  }
}
long long int piplib_ll_floor_div_r(long long int const a,
                                    long long int const b) {
  long long int q = piplib_ll_floor_div_q(a, b);
  return (a - q * b);
}

/* this routine is useless at present                                     */

int non_borne(tp, nvar, D, bigparm)
Tableau *tp;
int nvar, bigparm;
piplib_int_t D;
{int i, ff;
 for(i = 0; i<nvar; i++)
     {ff = Flag(tp, i);
      if(bigparm > 0)
	 {if(ff & Unit)return(Pip_True);
          if(Index(tp, i, bigparm) != D) return(Pip_True);
	 }
      }
 return(Pip_False);
}

/* This routine solve for z in the equation z.y = x (mod delta), provided
   y and delta are mutually prime. Remember that for multiple precision
   operation, the responsibility of creating and destroying <<z>> is the 
   caller's.                                                                */

void bezout(piplib_int_t x, piplib_int_t y, piplib_int_t delta, piplib_int_t *z){
  piplib_int_t a, b, c, d, e, f, u, v, q, r;
#if defined(PIPLIB_INT_GMP)
  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
  mpz_init(e); mpz_init(f); mpz_init(u); mpz_init(v);
  mpz_init(q); mpz_init(r);
  mpz_set_ui(a, 1); mpz_set_ui(b, 0); mpz_set_ui(c, 0);
  mpz_set_ui(d, 1); mpz_set(u, y); mpz_set(v, delta);
#else
  a = 1; b = 0; c = 0; d = 1;
  u = y; v = delta;
#endif
  for(;;){
#if defined(PIPLIB_INT_GMP)
    mpz_fdiv_qr(q, r, u, v);
    if(mpz_cmp_ui(r, 0) == 0) break;
    mpz_set(u, v);
    mpz_set(v, r);
    mpz_mul(e, q, c);
    mpz_sub(e, a, e);
    mpz_mul(f, q, d);
    mpz_sub(f, b, f);
    mpz_set(a, c);
    mpz_set(b, d);
    mpz_set(c, e);
    mpz_set(d, f);
#else
    q = u / v;
    piplib_int_mod(r, u, v);
    if(r == 0) break;
    u = v;
    v = r;
    e = a - q*c;
    f = b - q*d;
    a = c;
    b = d;
    c = e;
    d = f;
#endif
  }
#if defined(PIPLIB_INT_GMP)
  if(mpz_cmp_ui(v, 1) != 0)
    mpz_set_ui(*z, 0);
  else {
    mpz_mul(a, c, x);
    mpz_mod(*z, a, delta);
  }
  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);
  mpz_clear(e); mpz_clear(f); mpz_clear(u); mpz_clear(v);
  mpz_clear(q); mpz_clear(r);
  
#else
  if (v != 1) { *z = 0; } /* y and delta are not mutually prime */
  else { piplib_int_mod(*z, c * x, delta); }
#endif
}

Tableau *expanser();

/* cut: constant parameters denominator */
static int add_parm(Tableau **pcontext, int nr, int *pnparm, int *pni, int *pnc,
		    piplib_int_t *cut)
{
    int nparm = *pnparm;
    int i, j, k;
    piplib_int_t x;

    piplib_int_init(x);

/*        Build the definition of the new parameter into the solution :
      p_{nparm} = -(sum_{j=0}^{nparm-1} c_{nvar + 1 + j} p_j 
                     + c_{nvar})/D                             (3)
         The minus sign is there to compensate the one in (1)     */

    sol_new(nparm);
    sol_div();
    sol_forme(nparm+1);
    for (j = 0; j < nparm; j++) {
	piplib_int_oppose(x, cut[1+j]);
        sol_val(x, UN);
    }
    piplib_int_oppose(x, cut[0]);
    sol_val(x, UN);
    sol_val(cut[1+nparm], UN);		    /* The divisor                */

    if (nr+2 > (*pcontext)->height || nparm+1+1 > (*pcontext)->width) {
	int dcw, dch;
	dcw = piplib_int_size_in_base_2(cut[1+nparm]);
	dch = 2 * dcw + *pni;
	*pcontext = expanser(*pcontext, 0, nr, nparm+1, 0, dch, dcw);
    }

/* Since a new parameter is to be added, the constant term has to be moved
   right and a zero has to be inserted in all rows of the old context    */

    for (k = 0; k < nr; k++) {
	piplib_int_assign(Index(*pcontext, k, nparm+1), Index(*pcontext, k, nparm));
	piplib_int_set_si(Index(*pcontext, k, nparm), 0);
    }

/* The value of the new parameter is specified by applying the definition of
   Euclidean division to (3) :

 0<= - sum_{j=0}^{nparm-1} c_{nvar+1+j} p_j - c_{nvar} - D * p_{nparm} < D (4)

   This formula gives two inequalities which are stored in the context    */

    for (j = 0; j < nparm; j++) {
	piplib_int_oppose(Index(*pcontext, nr, j), cut[1+j]);
	piplib_int_assign(Index(*pcontext, nr+1, j), cut[1+j]);
    }
    piplib_int_oppose(Index(*pcontext, nr, nparm), cut[1+nparm]);
    piplib_int_assign(Index(*pcontext, nr+1, nparm), cut[1+nparm]);
    piplib_int_assign(x, cut[0]);
    piplib_int_oppose(Index(*pcontext, nr, nparm+1), x);
    piplib_int_decrement(x, x);
    piplib_int_add(Index(*pcontext, nr+1, nparm+1), x, cut[1+nparm]);

    Flag(*pcontext, nr) = Unknown;
    Flag(*pcontext, nr+1) = Unknown;
    piplib_int_set_si(Denom(*pcontext, nr), 1);
    piplib_int_set_si(Denom(*pcontext, nr+1), 1);
    (*pnparm)++;
    (*pnc) += 2;
    if (verbose > 0) {
	fprintf(dump, "enlarged context %d x %d\n", *pnparm, *pnc);
	fflush(dump);
    }

    piplib_int_clear(x);
}

static int has_cut(Tableau *context, int nr, int nparm, int p, piplib_int_t *cut)
{
    int row, col;

    for (row = 0; row < nr; ++row) {
	if (piplib_int_ne(Index(context, row, p), cut[1+nparm]))
	    continue;
	if (piplib_int_ne(Index(context, row, nparm), cut[0]))
	    continue;
	for (col = p+1; col < nparm; ++col)
	    if (piplib_int_zero(Index(context, row, col)) == 0)
		break;
	if (col < nparm)
	    continue;
	for (col = 0; col < p; ++col)
	    if (piplib_int_ne(Index(context, row, col), cut[1+col]))
		break;
	if (col < p)
	    continue;
	return 1;
    }
    return 0;
}

/* cut: constant parameters denominator */
static int find_parm(Tableau *context, int nr, int nparm, piplib_int_t *cut)
{
    int p;
    int col;
    int found;

    if (piplib_int_zero(cut[1+nparm-1]) == 0)
	return -1;

    piplib_int_add(cut[0], cut[0], cut[1+nparm]);
    piplib_int_decrement(cut[0], cut[0]);
    for (p = nparm-1; p >= 0; --p) {
	if (piplib_int_zero(cut[1+p]) == 0)
	    break;
	if (!has_cut(context, nr, nparm, p, cut))
	    continue;
	piplib_int_increment(cut[0], cut[0]);
	piplib_int_sub(cut[0], cut[0], cut[1+nparm]);
	for (col = 0; col < 1+nparm+1; ++col)
	    piplib_int_oppose(cut[col], cut[col]);
	found = has_cut(context, nr, nparm, p, cut);
	for (col = 0; col < 1+nparm+1; ++col)
	    piplib_int_oppose(cut[col], cut[col]);
	if (found)
	    return p;
	piplib_int_add(cut[0], cut[0], cut[1+nparm]);
	piplib_int_decrement(cut[0], cut[0]);
    }
    piplib_int_increment(cut[0], cut[0]);
    piplib_int_sub(cut[0], cut[0], cut[1+nparm]);
    return -1;
}

/* integrer(.....) add a cut to the problem tableau, or return 0 when an
   integral solution has been found, or -1 when no integral solution
   exists.

   Since integrer may add rows and columns to the problem tableau, its
   arguments are pointers rather than values. If a cut is constructed,
   ni increases by 1. If the cut is parametric, nparm increases by 1 and
   nc increases by 2.
									 */

int integrer(Tableau **ptp, Tableau **pcontext, 
	     int *pnvar, int *pnparm, int *pni, int *pnc, int bigparm)
{int ncol = *pnvar+*pnparm+1;
 int nligne = *pnvar + *pni;
 int nparm = *pnparm;
 int nvar = *pnvar;
 int ni = *pni;
 int nc = *pnc;
 piplib_int_t coupure[MAXCOL];
 int i, j, k, ff;
 piplib_int_t x, d;
 int ok_var, ok_const, ok_parm;
 piplib_int_t D;
    int parm;

 piplib_int_t t, delta, tau, lambda;

    if (ncol >= MAXCOL) {
	fprintf(stderr, "Too many variables: %d\n", ncol);
	exit(3);
    }
 
 #if defined(PIPLIB_INT_GMP)
 for(i=0; i<=ncol; i++)
   mpz_init(coupure[i]);

 mpz_init(x); mpz_init(d); mpz_init(D);
 mpz_init(t); mpz_init(delta); mpz_init(tau); mpz_init(lambda);
 #endif


/* search for a non-integral row */
 for(i = 0; i<nvar; i++) {
      #if defined(PIPLIB_INT_GMP)
      mpz_set(D, Denom(*ptp, i));
      if(mpz_cmp_ui(D, 1) == 0) continue;
      #else
      D = Denom(*ptp,i);
      if(D == 1) continue;
      #endif
/*                          If the common denominator of the row is 1
                            the row is integral                         */
      ff = Flag(*ptp, i);
      if(ff & Unit)continue;
/*                          If the row is a Unit, it is integral        */

/*                          Here a portential candidate has been found.
                            Build the cut by reducing each coefficient
                            modulo D, the common denominator            */
      ok_var = Pip_False;
      for(j = 0; j<nvar; j++) {
         #if defined(PIPLIB_INT_GMP)
         mpz_fdiv_r(x, Index(*ptp, i, j), D);
         mpz_set(coupure[j], x);
         #else
         piplib_int_mod(coupure[j], Index(*ptp, i, j), D);
         x = coupure[j];
         #endif
	    if (piplib_int_pos(x))
		ok_var = Pip_True;
          }
/*                          Done for the coefficient of the variables.  */

      #if defined(PIPLIB_INT_GMP)
      mpz_neg(x, Index(*ptp, i, nvar));
      mpz_fdiv_r(x, x, D);
      mpz_neg(x, x);
      mpz_set(coupure[nvar], x);
      ok_const = mpz_cmp_ui(x, 0);
      #else
      piplib_int_mod(coupure[nvar], -Index(*ptp, i, nvar), D);
	  coupure[nvar] = - coupure[nvar];
      x = coupure[nvar];
      ok_const = (x != 0);
      #endif
/*                          This is the constant term                   */
      ok_parm = Pip_False;
      for(j = nvar+1; j<ncol; j++) {
	 /* We assume that the big parameter is divisible by any number. */
	 if (j == bigparm) {
	    piplib_int_set_si(coupure[j], 0);
	    continue;
	 }
	 piplib_int_oppose(x, Index(*ptp, i, j));
	 piplib_int_mod(x, x, D);
	 piplib_int_oppose(coupure[j], x);
	 if (piplib_int_zero(coupure[j]) == 0)
	    ok_parm = Pip_True;
      }
/*                          These are the parametric terms              */

      #if defined(PIPLIB_INT_GMP)
      mpz_set(coupure[ncol], D);
      #else
      coupure[ncol] = D;    /* Just in case                             */
      #endif

/* The question now is whether the cut is valid. The answer is given
by the following decision table:

ok_var   ok_parm   ok_const

  F        F         F       (a) continue, integral row
  F        F         T       (b) return -1, no solution
  F        T         F       
                             (c) if the <<constant>> part is not divisible
                             by D then bottom else ....
  F        T         T
  T        F         F       (a) continue, integral row
  T        F         T       (d) constant cut
  T        T         F
                             (e) parametric cut
  T        T         T

                                                                case (a)  */

      if(!ok_parm && !ok_const) continue;
      if(!ok_parm)
          if(ok_var) {                                   /*     case (d)  */
              if(nligne >= (*ptp)->height) {
		  int d, dth, dtw;
                  #if defined(PIPLIB_INT_GMP)
	          d = mpz_sizeinbase(D, 2);
                  #else
                  d = llog(D);
                  #endif
                  dth = d;
		  *ptp = expanser(*ptp, nvar, ni, ncol, 0, dth, 0);
                  }
	      /* Find the deepest cut*/
	      if(deepest_cut){
#if defined(PIPLIB_INT_GMP)
	      mpz_neg(t, coupure[nvar]);
              mpz_gcd(delta, t, D);
	      mpz_divexact(tau, t, delta);
	      mpz_divexact(d, D, delta);
              mpz_sub_ui(t, d, 1);
              bezout(t, tau, d, &lambda);
	      mpz_gcd(t, lambda, D);
              while(mpz_cmp_ui(t, 1) != 0){
		mpz_add(lambda, lambda, d);
		mpz_gcd(t, lambda, D);
	      }
	      for(j=0; j<nvar; j++){
		mpz_mul(t, lambda, coupure[j]);
		mpz_fdiv_r(coupure[j], t, D);
	      }
	      mpz_mul(t, coupure[nvar], lambda);
	      mpz_mod(t, t, D);
	      mpz_sub(t, D, t);
	      mpz_neg(coupure[nvar], t);
#else
	      t = -coupure[nvar];
	      delta = piplib_llgcd_llabs(t,D);
	      tau = t/delta;
	      d = D/delta;
	      bezout(d-1, tau, d, &lambda);
	      while(piplib_llgcd_llabs(lambda, D) != 1) { lambda += d; }
	      for(j=0; j<nvar; j++) {
	        piplib_int_mod(coupure[j], lambda*coupure[j], D);
		  }
	      piplib_int_mod(coupure[nvar], -lambda*coupure[nvar], D);
	      coupure[nvar] = - coupure[nvar];
#endif
	      }
                         /* The cut has a negative <<constant>> part      */
              Flag(*ptp, nligne) = Minus; 
              #if defined(PIPLIB_INT_GMP)
              mpz_set(Denom(*ptp, nligne), D);
              #else
              Denom(*ptp, nligne) = D;
              #endif
                         /* Insert the cut */
	      for(j = 0; j<ncol; j++)
                  #if defined(PIPLIB_INT_GMP)
	          mpz_set(Index(*ptp, nligne, j), coupure[j]);
                  #else
                  Index(*ptp, nligne, j) = coupure[j];
                  #endif
                      /* A new row has been added to the problem tableau. */
	      (*pni)++;
              if(verbose > 0) {
		fprintf(dump, "just cut ");
                if(deepest_cut){
		  fprintf(dump, "Bezout multiplier ");
#if defined(PIPLIB_INT_GMP)
		  mpz_out_str(dump, 10, lambda);
#else
		  fprintf(dump, piplib_int_format, lambda);
#endif
		}
                fprintf(dump, "\n");
		k=0;
                for(i=0; i<nvar; i++){
                  if(Flag(*ptp, i) & Unit){
#if defined(PIPLIB_INT_GMP)
		    fprintf(dump, "0 ");
#else
		    sprintf(compose+k, "0 ");
#endif
		    k += 2;
		  }
		  else {
#if defined(PIPLIB_INT_GMP)
		    k += mpz_out_str(dump, 10, Index(*ptp, i, nvar));
		    fprintf(dump, "/");
		    k++;
		    k += mpz_out_str(dump, 10, Denom(*ptp, i));
		    fprintf(dump, " ");
		    k++;
		    if(k > 60){
		      putc('\n', dump);
		      k = 0;
		    }
#else
		    sprintf(compose+k, piplib_int_format, Index(*ptp, i, nvar));
		    k = strlen(compose);
		    sprintf(compose+k, "/");
		    k++;
		    sprintf(compose+k, piplib_int_format, Denom(*ptp, i));
		    k = strlen(compose);
		    sprintf(compose+k, " ");
		    k++;
		    if(k>60)  {
		      fputs(compose, dump);
		      putc('\n', dump);
		      k=0;
		    }
#endif
		  }
		}
		fputs(compose, dump);
		putc('\n', dump);
	      }
	      if(verbose > 2) tab_display(*ptp, dump);
	      goto clear;
              }
          else {                                         /*   case (b)    */
            nligne = -1; 
            goto clear;
          }
/* In case (e), one has to introduce a new parameter and
   introduce its defining inequalities into the context.
   
   Let the cut be    sum_{j=0}^{nvar-1} c_j x_j + c_{nvar} +             (2)
                     sum_{j=0}^{nparm-1} c_{nvar + 1 + j} p_j >= 0.       */
           
	parm = find_parm(*pcontext, nc, nparm, coupure+nvar);
	if (parm == -1) {
	    add_parm(pcontext, nc, pnparm, pni, pnc, coupure+nvar);
	    parm = nparm;
	}

	assert(ok_var);
          if(nligne >= (*ptp)->height || ncol >= (*ptp)->width) {
              int d, dth, dtw;
             #if defined(PIPLIB_INT_GMP)
             d = mpz_sizeinbase(D, 2);
             #else
             d = llog(D);
             #endif
              dth = d + ni;
	      dtw = d;
	      *ptp = expanser(*ptp, nvar, ni, ncol, 0, dth, dtw);
              }
                         /* Zeroing out the new column seems to be useless
			    since <<expanser>> does it anyway            */
                            
			 /* The cut has a negative <<constant>> part    */
	  Flag(*ptp, nligne) = Minus;
          #if defined(PIPLIB_INT_GMP)
          mpz_set(Denom(*ptp, nligne), D);
          #else
	  Denom(*ptp, nligne) = D;
          #endif
              	 /* Insert the cut */
	for (j = 0; j < ncol; j++)
              #if defined(PIPLIB_INT_GMP)
              mpz_set(Index(*ptp, nligne, j), coupure[j]);
              #else
	      Index(*ptp, nligne, j) = coupure[j];
              #endif
	piplib_int_add(Index(*ptp, nligne, nvar+1+parm),
		    Index(*ptp, nligne, nvar+1+parm), coupure[ncol]);
		 /* A new row has been added to the problem tableau.    */
	  (*pni)++;
          goto clear;
      }
 /* The solution is integral.                              */
    nligne = 0;
clear: 
   for(i=0; i <= ncol; i++)
	piplib_int_clear(coupure[i]);
    piplib_int_clear(x); piplib_int_clear(d); piplib_int_clear(D);
    piplib_int_clear(t); piplib_int_clear(tau); piplib_int_clear(lambda); piplib_int_clear(delta);
    return nligne;
}

