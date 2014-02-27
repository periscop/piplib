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

/* this routine is useless at present                                     */
int non_borne(tp, nvar, D, bigparm)
Tableau *tp;
int nvar, bigparm;
osl_int_t D;
{int i, ff;
 for(i = 0; i<nvar; i++)
     {ff = Flag(tp, i);
      if(bigparm > 0)
	 {if(ff & Unit)return(Pip_True);
          if(osl_int_ne(PIPLIB_INT_PRECISION, Index(tp, i, bigparm), D)) return(Pip_True);
	 }
      }
 return(Pip_False);
}

/* This routine solve for z in the equation z.y = x (mod delta), provided
   y and delta are mutually prime. Remember that for multiple precision
   operation, the responsibility of creating and destroying <<z>> is the 
   caller's.                                                                */

void bezout(osl_int_t x, osl_int_t y, osl_int_t delta, osl_int_t *z){
  osl_int_t a, b, c, d, e, f, u, v, q, r;
  osl_int_init(PIPLIB_INT_PRECISION, &a);
  osl_int_init(PIPLIB_INT_PRECISION, &b);
  osl_int_init(PIPLIB_INT_PRECISION, &c);
  osl_int_init(PIPLIB_INT_PRECISION, &d);
  osl_int_init(PIPLIB_INT_PRECISION, &e);
  osl_int_init(PIPLIB_INT_PRECISION, &f);
  osl_int_init(PIPLIB_INT_PRECISION, &u);
  osl_int_init(PIPLIB_INT_PRECISION, &v);
  osl_int_init(PIPLIB_INT_PRECISION, &q);
  osl_int_init(PIPLIB_INT_PRECISION, &r);
  osl_int_set_si(PIPLIB_INT_PRECISION, &a, 1);
  osl_int_set_si(PIPLIB_INT_PRECISION, &b, 0);
  osl_int_set_si(PIPLIB_INT_PRECISION, &c, 0);
  osl_int_set_si(PIPLIB_INT_PRECISION, &d, 1);
  osl_int_assign(PIPLIB_INT_PRECISION, &u, y);
  osl_int_assign(PIPLIB_INT_PRECISION, &v, delta);
  for(;;){
    osl_int_floor_div_q_r(PIPLIB_INT_PRECISION, &q, &r, u, v);
    if (osl_int_zero(PIPLIB_INT_PRECISION, r)) break;
    osl_int_assign(PIPLIB_INT_PRECISION, &u, v);
    osl_int_assign(PIPLIB_INT_PRECISION, &v, r);
    osl_int_mul(PIPLIB_INT_PRECISION, &e, q, c);
    osl_int_sub(PIPLIB_INT_PRECISION, &e, a, e);
    osl_int_mul(PIPLIB_INT_PRECISION, &f, q, d);
    osl_int_sub(PIPLIB_INT_PRECISION, &f, b, f);
    osl_int_assign(PIPLIB_INT_PRECISION, &a, c);
    osl_int_assign(PIPLIB_INT_PRECISION, &b, d);
    osl_int_assign(PIPLIB_INT_PRECISION, &c, e);
    osl_int_assign(PIPLIB_INT_PRECISION, &d, f);
  }
  if (! osl_int_one(PIPLIB_INT_PRECISION, v))
    osl_int_set_si(PIPLIB_INT_PRECISION, z, 0);
  else {
    osl_int_mul(PIPLIB_INT_PRECISION, &a, c, x);
    osl_int_mod(PIPLIB_INT_PRECISION, z, a, delta);
  }
  osl_int_clear(PIPLIB_INT_PRECISION, &a);
  osl_int_clear(PIPLIB_INT_PRECISION, &b);
  osl_int_clear(PIPLIB_INT_PRECISION, &c);
  osl_int_clear(PIPLIB_INT_PRECISION, &d);
  osl_int_clear(PIPLIB_INT_PRECISION, &e);
  osl_int_clear(PIPLIB_INT_PRECISION, &f);
  osl_int_clear(PIPLIB_INT_PRECISION, &u);
  osl_int_clear(PIPLIB_INT_PRECISION, &v);
  osl_int_clear(PIPLIB_INT_PRECISION, &q);
  osl_int_clear(PIPLIB_INT_PRECISION, &r);
}

Tableau *expanser();

/* cut: constant parameters denominator */
static void add_parm(Tableau **pcontext, int nr, int *pnparm, int *pni, int *pnc,
		    osl_int_t *cut)
{
    int nparm = *pnparm;
    int j, k;
    osl_int_t x;

    osl_int_init(PIPLIB_INT_PRECISION, &x);

/*        Build the definition of the new parameter into the solution :
      p_{nparm} = -(sum_{j=0}^{nparm-1} c_{nvar + 1 + j} p_j 
                     + c_{nvar})/D                             (3)
         The minus sign is there to compensate the one in (1)     */

    sol_new(nparm);
    sol_div();
    sol_forme(nparm+1);
    for (j = 0; j < nparm; j++) {
	osl_int_oppose(PIPLIB_INT_PRECISION, &x, cut[1+j]);
        sol_val(x, UN);
    }
    osl_int_oppose(PIPLIB_INT_PRECISION, &x, cut[0]);
    sol_val(x, UN);
    sol_val(cut[1+nparm], UN);		    /* The divisor                */

    if (nr+2 > (*pcontext)->height || nparm+1+1 > (*pcontext)->width) {
	int dcw, dch;
	dcw = osl_int_size_in_base_2(PIPLIB_INT_PRECISION, cut[1+nparm]);
	dch = 2 * dcw + *pni;
	*pcontext = expanser(*pcontext, 0, nr, nparm+1, 0, dch, dcw);
    }

/* Since a new parameter is to be added, the constant term has to be moved
   right and a zero has to be inserted in all rows of the old context    */

    for (k = 0; k < nr; k++) {
	osl_int_assign(PIPLIB_INT_PRECISION, &Index(*pcontext, k, nparm+1), Index(*pcontext, k, nparm));
	osl_int_set_si(PIPLIB_INT_PRECISION, &Index(*pcontext, k, nparm), 0);
    }

/* The value of the new parameter is specified by applying the definition of
   Euclidean division to (3) :

 0<= - sum_{j=0}^{nparm-1} c_{nvar+1+j} p_j - c_{nvar} - D * p_{nparm} < D (4)

   This formula gives two inequalities which are stored in the context    */

    for (j = 0; j < nparm; j++) {
	osl_int_oppose(PIPLIB_INT_PRECISION, &Index(*pcontext, nr, j), cut[1+j]);
	osl_int_assign(PIPLIB_INT_PRECISION, &Index(*pcontext, nr+1, j), cut[1+j]);
    }
    osl_int_oppose(PIPLIB_INT_PRECISION, &Index(*pcontext, nr, nparm), cut[1+nparm]);
    osl_int_assign(PIPLIB_INT_PRECISION, &Index(*pcontext, nr+1, nparm), cut[1+nparm]);
    osl_int_assign(PIPLIB_INT_PRECISION, &x, cut[0]);
    osl_int_oppose(PIPLIB_INT_PRECISION, &Index(*pcontext, nr, nparm+1), x);
    osl_int_decrement(PIPLIB_INT_PRECISION, &x, x);
    osl_int_add(PIPLIB_INT_PRECISION, &Index(*pcontext, nr+1, nparm+1), x, cut[1+nparm]);

    Flag(*pcontext, nr) = Unknown;
    Flag(*pcontext, nr+1) = Unknown;
    osl_int_set_si(PIPLIB_INT_PRECISION, &Denom(*pcontext, nr), 1);
    osl_int_set_si(PIPLIB_INT_PRECISION, &Denom(*pcontext, nr+1), 1);
    (*pnparm)++;
    (*pnc) += 2;
    if (verbose > 0) {
	fprintf(dump, "enlarged context %d x %d\n", *pnparm, *pnc);
	fflush(dump);
    }

    osl_int_clear(PIPLIB_INT_PRECISION, &x);
}

static int has_cut(Tableau *context, int nr, int nparm, int p, osl_int_t *cut)
{
    int row, col;

    for (row = 0; row < nr; ++row) {
	if (osl_int_ne(PIPLIB_INT_PRECISION, Index(context, row, p), cut[1+nparm]))
	    continue;
	if (osl_int_ne(PIPLIB_INT_PRECISION, Index(context, row, nparm), cut[0]))
	    continue;
	for (col = p+1; col < nparm; ++col)
	    if (! osl_int_zero(PIPLIB_INT_PRECISION, Index(context, row, col)))
		break;
	if (col < nparm)
	    continue;
	for (col = 0; col < p; ++col)
	    if (osl_int_ne(PIPLIB_INT_PRECISION, Index(context, row, col), cut[1+col]))
		break;
	if (col < p)
	    continue;
	return 1;
    }
    return 0;
}

/* cut: constant parameters denominator */
static int find_parm(Tableau *context, int nr, int nparm, osl_int_t *cut)
{
    int p;
    int col;
    int found;

    if (! osl_int_zero(PIPLIB_INT_PRECISION, cut[1+nparm-1]))
	return -1;

    osl_int_add(PIPLIB_INT_PRECISION, &cut[0], cut[0], cut[1+nparm]);
    osl_int_decrement(PIPLIB_INT_PRECISION, &cut[0], cut[0]);
    for (p = nparm-1; p >= 0; --p) {
	if (! osl_int_zero(PIPLIB_INT_PRECISION, cut[1+p]))
	    break;
	if (!has_cut(context, nr, nparm, p, cut))
	    continue;
	osl_int_increment(PIPLIB_INT_PRECISION, &cut[0], cut[0]);
	osl_int_sub(PIPLIB_INT_PRECISION, &cut[0], cut[0], cut[1+nparm]);
	for (col = 0; col < 1+nparm+1; ++col)
	    osl_int_oppose(PIPLIB_INT_PRECISION, &cut[col], cut[col]);
	found = has_cut(context, nr, nparm, p, cut);
	for (col = 0; col < 1+nparm+1; ++col)
	    osl_int_oppose(PIPLIB_INT_PRECISION, &cut[col], cut[col]);
	if (found)
	    return p;
	osl_int_add(PIPLIB_INT_PRECISION, &cut[0], cut[0], cut[1+nparm]);
	osl_int_decrement(PIPLIB_INT_PRECISION, &cut[0], cut[0]);
    }
    osl_int_increment(PIPLIB_INT_PRECISION, &cut[0], cut[0]);
    osl_int_sub(PIPLIB_INT_PRECISION, &cut[0], cut[0], cut[1+nparm]);
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
 osl_int_t coupure[MAXCOL];
 int i, j, k, ff;
 osl_int_t x, d;
 int ok_var, ok_const, ok_parm;
 osl_int_t D;
    int parm;

 osl_int_t t, delta, tau, lambda;

    if (ncol >= MAXCOL) {
	fprintf(stderr, "Too many variables: %d\n", ncol);
	exit(3);
    }

 for(i=0; i<=ncol; i++)
   osl_int_init(PIPLIB_INT_PRECISION, &coupure[i]);

 osl_int_init(PIPLIB_INT_PRECISION, &x);
 osl_int_init(PIPLIB_INT_PRECISION, &d);
 osl_int_init(PIPLIB_INT_PRECISION, &D);
 osl_int_init(PIPLIB_INT_PRECISION, &t);
 osl_int_init(PIPLIB_INT_PRECISION, &delta);
 osl_int_init(PIPLIB_INT_PRECISION, &tau);
 osl_int_init(PIPLIB_INT_PRECISION, &lambda);


/* search for a non-integral row */
 for(i = 0; i<nvar; i++) {
      osl_int_assign(PIPLIB_INT_PRECISION, &D, Denom(*ptp, i));
      if (osl_int_one(PIPLIB_INT_PRECISION, D)) continue;
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
         osl_int_floor_div_r(PIPLIB_INT_PRECISION, &x, Index(*ptp, i, j), D);
         osl_int_assign(PIPLIB_INT_PRECISION, &coupure[j], x);
	    if (osl_int_pos(PIPLIB_INT_PRECISION, x))
		ok_var = Pip_True;
          }
/*                          Done for the coefficient of the variables.  */

      osl_int_oppose(PIPLIB_INT_PRECISION, &x, Index(*ptp, i, nvar));
      osl_int_floor_div_r(PIPLIB_INT_PRECISION, &x, x, D);
      osl_int_oppose(PIPLIB_INT_PRECISION, &x, x);
      osl_int_assign(PIPLIB_INT_PRECISION, &coupure[nvar], x);
      ok_const = ! osl_int_zero(PIPLIB_INT_PRECISION, x);
/*                          This is the constant term                   */
      ok_parm = Pip_False;
      for(j = nvar+1; j<ncol; j++) {
	 /* We assume that the big parameter is divisible by any number. */
	 if (j == bigparm) {
	    osl_int_set_si(PIPLIB_INT_PRECISION, &coupure[j], 0);
	    continue;
	 }
	 osl_int_oppose(PIPLIB_INT_PRECISION, &x, Index(*ptp, i, j));
	 osl_int_floor_div_r(PIPLIB_INT_PRECISION, &x, x, D);
	 osl_int_oppose(PIPLIB_INT_PRECISION, &coupure[j], x);
	 if (! osl_int_zero(PIPLIB_INT_PRECISION, coupure[j]))
	    ok_parm = Pip_True;
      }
/*                          These are the parametric terms              */

      osl_int_assign(PIPLIB_INT_PRECISION, &coupure[ncol], D);

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
      if(!ok_parm) {
          if(ok_var) {                                   /*     case (d)  */
              if(nligne >= (*ptp)->height) {
		  int d, dth;
	          d = osl_int_size_in_base_2(PIPLIB_INT_PRECISION, D);
                  dth = d;
		  *ptp = expanser(*ptp, nvar, ni, ncol, 0, dth, 0);
                  }
	      /* Find the deepest cut*/
	      if(deepest_cut){
	      osl_int_oppose(PIPLIB_INT_PRECISION, &t, coupure[nvar]);
              osl_int_gcd(PIPLIB_INT_PRECISION, &delta, t, D);
	      osl_int_div_exact(PIPLIB_INT_PRECISION, &tau, t, delta);
	      osl_int_div_exact(PIPLIB_INT_PRECISION, &d, D, delta);
              osl_int_add_si(PIPLIB_INT_PRECISION, &t, d, -1);
              bezout(t, tau, d, &lambda);
	      osl_int_gcd(PIPLIB_INT_PRECISION, &t, lambda, D);
              while(! osl_int_one(PIPLIB_INT_PRECISION, t)){
		osl_int_add(PIPLIB_INT_PRECISION, &lambda, lambda, d);
		osl_int_gcd(PIPLIB_INT_PRECISION, &t, lambda, D);
	      }
	      for(j=0; j<nvar; j++){
		osl_int_mul(PIPLIB_INT_PRECISION, &t, lambda, coupure[j]);
		osl_int_floor_div_r(PIPLIB_INT_PRECISION, &coupure[j], t, D);
	      }
	      osl_int_mul(PIPLIB_INT_PRECISION, &t, coupure[nvar], lambda);
	      osl_int_mod(PIPLIB_INT_PRECISION, &t, t, D);
	      osl_int_sub(PIPLIB_INT_PRECISION, &t, D, t);
	      osl_int_oppose(PIPLIB_INT_PRECISION, &coupure[nvar], t);
	      }
                         /* The cut has a negative <<constant>> part      */
              Flag(*ptp, nligne) = Minus; 
              osl_int_assign(PIPLIB_INT_PRECISION, &Denom(*ptp, nligne), D);
                         /* Insert the cut */
	      for(j = 0; j<ncol; j++)
	          osl_int_assign(PIPLIB_INT_PRECISION, &Index(*ptp, nligne, j), coupure[j]);
                      /* A new row has been added to the problem tableau. */
	      (*pni)++;
              if(verbose > 0) {
		fprintf(dump, "just cut ");
                if(deepest_cut){
		  fprintf(dump, "Bezout multiplier ");
		  osl_int_print(dump, PIPLIB_INT_PRECISION, lambda);
		}
                fprintf(dump, "\n");
		k=0;
                for(i=0; i<nvar; i++){
                  if(Flag(*ptp, i) & Unit){
		    fprintf(dump, "0 ");
		    k += 2;
		  }
		  else {
		    osl_int_print(dump, PIPLIB_INT_PRECISION, Index(*ptp, i, nvar));
		    fprintf(dump, "/");
		    osl_int_print(dump, PIPLIB_INT_PRECISION, Denom(*ptp, i));
		    fprintf(dump, " ");
			k += osl_int_size_in_base_10(PIPLIB_INT_PRECISION, Index(*ptp, i, nvar)) + 1 +
			     osl_int_size_in_base_10(PIPLIB_INT_PRECISION, Denom(*ptp, i)) + 1;
		    if(k > 60){
		      putc('\n', dump);
		      k = 0;
		    }
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
             d = osl_int_size_in_base_2(PIPLIB_INT_PRECISION, D);
              dth = d + ni;
	      dtw = d;
	      *ptp = expanser(*ptp, nvar, ni, ncol, 0, dth, dtw);
              }
                         /* Zeroing out the new column seems to be useless
			    since <<expanser>> does it anyway            */
                            
			 /* The cut has a negative <<constant>> part    */
	  Flag(*ptp, nligne) = Minus;
          osl_int_assign(PIPLIB_INT_PRECISION, &Denom(*ptp, nligne), D);
              	 /* Insert the cut */
	for (j = 0; j < ncol; j++)
              osl_int_assign(PIPLIB_INT_PRECISION, &Index(*ptp, nligne, j), coupure[j]);
	osl_int_add(PIPLIB_INT_PRECISION, &Index(*ptp, nligne, nvar+1+parm),
		    Index(*ptp, nligne, nvar+1+parm), coupure[ncol]);
		 /* A new row has been added to the problem tableau.    */
	  (*pni)++;
          goto clear;
      }
 /* The solution is integral.                              */
    nligne = 0;
clear: 
   for(i=0; i <= ncol; i++)
	osl_int_clear(PIPLIB_INT_PRECISION, &coupure[i]);
    osl_int_clear(PIPLIB_INT_PRECISION, &x); osl_int_clear(PIPLIB_INT_PRECISION, &d); osl_int_clear(PIPLIB_INT_PRECISION, &D);
    osl_int_clear(PIPLIB_INT_PRECISION, &t); osl_int_clear(PIPLIB_INT_PRECISION, &tau); osl_int_clear(PIPLIB_INT_PRECISION, &lambda); osl_int_clear(PIPLIB_INT_PRECISION, &delta);
    return nligne;
}
