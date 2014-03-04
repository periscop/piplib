/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 traiter.c                                  *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988-2005                                        *
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

#include <osl/macros.h>
#include <osl/int.h>

#include "pip.h"


extern long int cross_product, limit;
extern int verbose;
extern FILE* dump;
extern int profondeur;
extern int compa_count;


int chercher(Tableau* p, int masque, int n)
{
  int i;
  for (i = 0; i < n; i++)
    if (p->row[i].flags & masque) break;
  return i;
}

/* il est convenu que traiter ne doit modifier ni le tableau, ni le contexte;
   le tableau peut grandir en cas de coupure (+1 en hauteur et +1 en largeur
   si nparm != 0) et en cas de partage (+1 en hauteur)(seulement si nparm != 0).
   le contexte peut grandir en cas de coupure (+2 en hauteur et +1 en largeur)
   (seulement si nparm !=0) et en cas de partage (+1 en hauteur)(nparm !=0).
   On estime le nombre de coupures a llog(D) et le nombre de partages a
   ni. */
Tableau* expanser(
  Tableau* tp, int virt, int reel, int ncol, int off, int dh, int dw
)
{
  int i, j, ff;
  osl_int_t* pq;
  osl_int_t* pp;
  osl_int_t* qq;
  Tableau* rp;
  if(tp == NULL) return(NULL);
  rp = tab_alloc(reel+dh, ncol+dw, virt);

  osl_int_assign(PIPLIB_INT_PRECISION, &rp->determinant, tp->determinant);
  pq = (osl_int_t *) & (rp->row[virt+reel+dh]);
  for(i = off; i<virt + reel; i++)
  {
    ff = Flag(rp, i) = Flag(tp, i-off);
    osl_int_assign(PIPLIB_INT_PRECISION, &Denom(rp, i), Denom(tp, i-off));
    if (ff & Unit) rp->row[i].objet.unit = tp->row[i-off].objet.unit;
    else {
      rp->row[i].objet.val = pq;
      pq +=(ncol + dw);
      pp = tp->row[i-off].objet.val;
      qq = rp->row[i].objet.val;
      for (j = 0; j<ncol; j++)
        osl_int_assign(PIPLIB_INT_PRECISION, qq++, *pp++);
      }
  }

  return(rp);
}

/* Check for "obvious" signs of the parametric constant terms
 * of the inequalities.  As soon as a negative sign is found
 * we return from this function and handle this constraint
 * in the calling function.  The signs of the other constraints
 * are then mostly irrelevant.
 * If any of the negative signs is due to the "big parameter",
 * then we want to use this constraint first.
 * We therefore check for signs determined by the coefficient
 * of the big parameter first. */
int exam_coef(Tableau *tp, int nvar, int ncol, int bigparm)
{
  int i, j ;
  int ff, fff;
  osl_int_t* p;

  if (bigparm >= 0)
    for (i = 0; i<tp->height; i++) {
      if (Flag(tp, i) != Unknown) continue;

      if (osl_int_neg(PIPLIB_INT_PRECISION, Index(tp,i, bigparm))) {
        Flag(tp, i) = Minus;
      }
      else if (osl_int_pos(PIPLIB_INT_PRECISION, Index(tp,i, bigparm))) {
        Flag(tp, i) = Plus;
      }
    }

  for(i = 0; i<tp->height; i++) {
    ff = Flag(tp,i);
    if(ff == 0) break;
    if(ff == Unknown) {
      ff = Zero;
      p = &(tp->row[i].objet.val[nvar+1]);
      for(j = nvar+1; j<ncol; j++) {

        if (osl_int_pos(PIPLIB_INT_PRECISION, *p)) fff = Plus;
        else if (osl_int_neg(PIPLIB_INT_PRECISION, *p)) fff = Minus;
        else fff = Zero;

        ++p;

        if (fff != Zero && fff != ff) {
          if(ff == Zero) ff = fff;
          else { ff = Unknown; break; }
        }
      }
      /* bug de'tecte' par [paf], 16/2/93 !
         Si tous les coefficients des parame`tres sont ne'gatifs
         et si le terme constant est nul, le signe est inconnu!!
         On traite donc spe'cialement le terme constant. */
      if (osl_int_pos(PIPLIB_INT_PRECISION, Index(tp, i, nvar))) fff = Plus;
      else if (osl_int_neg(PIPLIB_INT_PRECISION, Index(tp, i, nvar))) fff = Minus;
      else fff = Zero;

      /* ici on a le signe du terme constant */
      switch(ff) {
        /* le signe est inconnu si les coefficients sont positifs et
           le terme constant ne'gatif */
        case Plus: if(fff == Minus) ff = Unknown; break;
        /* si les coefficients sont tous nuls, le signe est celui
           du terme constant */
        case Zero: ff = fff; break;
        /* le signe est inconnu si les coefficients sont ne'gatifs,
           sauf si le terme constant est egalement negatif. */
        case Minus: if(fff != Minus) ff = Unknown; break;
        /* enfin, il n'y a rien a` dire si le signe des coefficients est inconnu */
      }
      Flag(tp, i) = ff;
      if(ff == Minus) return i;
    }
  }

  return i;
}

void compa_test(Tableau *tp, Tableau *context,
		int ni, int nvar, int nparm, int nc)
{
 int i, j;
 int ff;
 int cPlus, cMinus, isCritic;
 Tableau *tPlus, *tMinus;
 int p;
 struct high_water_mark q;

 if(nparm == 0) return;
 if(nparm >= MAXPARM) {
     fprintf(stderr, "Too much parameters : %d\n", nparm);
     exit(1);
     }
 q = tab_hwm();

 for(i = 0; i<ni + nvar; i++)
     {ff = Flag(tp,i);
      if(ff & (Critic | Unknown))
	  {isCritic = Pip_True;
	   for(j = 0; j<nvar; j++)
		 if(osl_int_pos(PIPLIB_INT_PRECISION, Index(tp, i, j)))
		 {isCritic = Pip_False;
		  break;
		 }
           compa_count++;
	   tPlus = expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tPlus, nparm+nc) = Unknown;
	   for (j = 0; j < nparm; j++)
	       osl_int_assign(PIPLIB_INT_PRECISION, &Index(tPlus, nparm+nc, j), Index(tp, i, j+nvar+1));
	   osl_int_assign(PIPLIB_INT_PRECISION, &Index(tPlus, nparm+nc, nparm), Index(tp, i, nvar));
	   if (!isCritic)
	       osl_int_decrement(PIPLIB_INT_PRECISION, &Index(tPlus, nparm+nc, nparm),
				    Index(tPlus, nparm+nc, nparm));
	   osl_int_assign(PIPLIB_INT_PRECISION, &Denom(tPlus, nparm+nc), UN);
	   
	   p = sol_hwm();
	   traiter(tPlus, NULL, nparm, 0, nc+1, 0, -1, TRAITER_INT);
	   cPlus = is_not_Nil(p);
	   if(verbose>0){
	     fprintf(dump, "\nThe positive case has been found ");
	     fprintf(dump, cPlus? "possible\n": "impossible\n");
	     fflush(dump);
	   }

	   sol_reset(p);
	   tMinus = expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tMinus, nparm+nc) = Unknown;
	   for (j = 0; j < nparm; j++)
	       osl_int_oppose(PIPLIB_INT_PRECISION, &Index(tMinus, nparm+nc, j), Index(tp, i, j+nvar+1));
	   osl_int_oppose(PIPLIB_INT_PRECISION, &Index(tMinus, nparm+nc, nparm), Index(tp, i, nvar));
	   osl_int_decrement(PIPLIB_INT_PRECISION, &Index(tMinus, nparm+nc, nparm),
				Index(tMinus, nparm+nc, nparm));
	   osl_int_assign(PIPLIB_INT_PRECISION, &Denom(tMinus, nparm+nc), UN);
	   traiter(tMinus, NULL, nparm, 0, nc+1, 0, -1, TRAITER_INT);
	   cMinus = is_not_Nil(p);
	   if(verbose>0){
	     fprintf(dump, "\nThe negative case has been found ");
	     fprintf(dump, cMinus? "possible\n": "impossible\n");
	     fflush(dump);
	   }

	   sol_reset(p);
	   if (cPlus && cMinus) {
	       Flag(tp,i) = isCritic ? Critic : Unknown;
	     }
	   else if (cMinus)
	      {Flag(tp,i) = Minus;
	       break;
	      }
	   else {
	     Flag(tp,i) = cPlus ? Plus : Zero;
	   }
	  }
     }
 tab_reset(q);
 
 return;
}

osl_int_t *valeur(Tableau *tp, int i, int j)
{
 if(Flag(tp, i) & Unit)
     return(tp->row[i].objet.unit == j ? &Denom(tp,i) : &ZERO);
 else return(&Index(tp, i, j));
}

void solution(Tableau *tp, int nvar, int nparm)
{int i, j;
 int ncol = nvar + nparm + 1;

 sol_list(nvar);
 for(i = 0; i<nvar; i++)
     {sol_forme(nparm+1);
      for(j = nvar+1; j<ncol; j++)
	 sol_val(*valeur(tp, i, j), Denom(tp,i));
      sol_val(*valeur(tp, i, nvar), Denom(tp,i));
     }
}

static void solution_dual(Tableau *tp, int nvar, int *pos)
{
    int i;

    sol_list(tp->height - nvar);
    for (i = 0; i < tp->height - nvar; ++i) {
	sol_forme(1);
	if (Flag(tp, pos[i]) & Unit)
	    sol_val(*valeur(tp, 0, tp->row[pos[i]].objet.unit), Denom(tp, 0));
	else
	    sol_val(ZERO, UN);
    }
}

int choisir_piv(Tableau *tp, int pivi, int nvar, int nligne)
{
  int j, k;
  osl_int_t pivot, foo, x, y;
  int pivj = -1;

  osl_int_init(PIPLIB_INT_PRECISION, &pivot);
  osl_int_init(PIPLIB_INT_PRECISION, &foo);
  osl_int_init(PIPLIB_INT_PRECISION, &x);
  osl_int_init(PIPLIB_INT_PRECISION, &y);

  for(j = 0; j<nvar; j++) {
    osl_int_assign(PIPLIB_INT_PRECISION, &foo, Index(tp, pivi, j));
    if(!osl_int_pos(PIPLIB_INT_PRECISION, foo)) continue;
    if(pivj < 0) {
      pivj = j;
      osl_int_assign(PIPLIB_INT_PRECISION, &pivot, foo);
      continue;
    }
    for(k = 0; k<nligne; k++) {
      osl_int_mul(PIPLIB_INT_PRECISION, &x, pivot, *valeur(tp, k, j)); 
      osl_int_mul(PIPLIB_INT_PRECISION, &y, *valeur(tp, k, pivj), foo);
      osl_int_sub(PIPLIB_INT_PRECISION, &x, x, y);
      cross_product++;

       if (! osl_int_zero(PIPLIB_INT_PRECISION, x)) break;
    }
    if(osl_int_neg(PIPLIB_INT_PRECISION, x)) {
      pivj = j;
      osl_int_assign(PIPLIB_INT_PRECISION, &pivot, foo);
    }
  }

  osl_int_clear(PIPLIB_INT_PRECISION, &pivot);
  osl_int_clear(PIPLIB_INT_PRECISION, &foo);
  osl_int_clear(PIPLIB_INT_PRECISION, &x);
  osl_int_clear(PIPLIB_INT_PRECISION, &y);

  return pivj;
}


int pivoter(Tableau* tp, int pivi, int nvar, int nparm, int ni)
{
 int pivj;
 int ncol = nvar + nparm + 1;
 int nligne = nvar + ni;
 int i, j, k;
 osl_int_t x, y, d, gcd, u, dpiv;
 int ff, fff;
 osl_int_t pivot, foo, z;
 osl_int_t ppivot, dppiv;
 osl_int_t new[MAXCOL], *p, *q;
 osl_int_t lpiv;

 if(ncol >= MAXCOL) {
   fprintf(stdout, "Too much variables\n");
   exit(1);
 }
 if(0 > pivi || pivi >= nligne || Flag(tp, pivi) == Unit) {
   fprintf(stdout, "Syserr : pivoter : wrong pivot row\n");
   exit(1);
 }

 pivj = choisir_piv(tp, pivi, nvar, nligne);
 if(pivj < 0) return(-1);
 if(pivj >= nvar) {
   fprintf(stdout, "Syserr : pivoter : wrong pivot\n");
   exit(1);
 }

 osl_int_init(PIPLIB_INT_PRECISION, &x);
 osl_int_init(PIPLIB_INT_PRECISION, &y);
 osl_int_init(PIPLIB_INT_PRECISION, &d); 
 osl_int_init(PIPLIB_INT_PRECISION, &gcd);
 osl_int_init(PIPLIB_INT_PRECISION, &u);
 osl_int_init(PIPLIB_INT_PRECISION, &dpiv);
 osl_int_init(PIPLIB_INT_PRECISION, &lpiv);
 osl_int_init(PIPLIB_INT_PRECISION, &pivot);
 osl_int_init(PIPLIB_INT_PRECISION, &foo);
 osl_int_init(PIPLIB_INT_PRECISION, &z);
 osl_int_init(PIPLIB_INT_PRECISION, &ppivot);
 osl_int_init(PIPLIB_INT_PRECISION, &dppiv);

 for(i=0; i<ncol; i++)
   osl_int_init(PIPLIB_INT_PRECISION, &new[i]);

 osl_int_assign(PIPLIB_INT_PRECISION, &pivot, Index(tp, pivi, pivj));
 osl_int_assign(PIPLIB_INT_PRECISION, &dpiv, Denom(tp, pivi));
 osl_int_gcd(PIPLIB_INT_PRECISION, &d, pivot, dpiv);
 osl_int_div_exact(PIPLIB_INT_PRECISION, &ppivot, pivot, d);
 osl_int_div_exact(PIPLIB_INT_PRECISION, &dppiv, dpiv, d);
 
 if(verbose>1){
   fprintf(dump, "Pivot ");
   osl_int_print(dump, PIPLIB_INT_PRECISION, ppivot);
   putc('/', dump);
   osl_int_print(dump, PIPLIB_INT_PRECISION, dppiv);
   putc('\n', dump);
   fprintf(dump, "%d x %d\n", pivi, pivj);
 }

 osl_int_floor_div_q_r(PIPLIB_INT_PRECISION, &x, &y, tp->determinant, dppiv);

 if (!osl_int_zero(PIPLIB_INT_PRECISION, y)) {
   fprintf(stderr, "Integer overflow\n");
   if(verbose>0) fflush(dump);
   exit(1);
 }
 
 osl_int_mul(PIPLIB_INT_PRECISION, &tp->determinant, x, ppivot);

 if(verbose>1){
   fprintf(dump, "determinant ");
   osl_int_print(dump, PIPLIB_INT_PRECISION, tp->determinant);
   fprintf(dump, "\n");
 }

 
 for(j = 0; j<ncol; j++)
   if(j==pivj)
     osl_int_assign(PIPLIB_INT_PRECISION, &new[j], dpiv);
   else 
     osl_int_oppose(PIPLIB_INT_PRECISION, &new[j], Index(tp, pivi, j));

 for(k = 0; k<nligne; k++){
   if(Flag(tp,k) & Unit)continue;
   if(k == pivi)continue;
   osl_int_assign(PIPLIB_INT_PRECISION, &foo, Index(tp, k, pivj));
   osl_int_gcd(PIPLIB_INT_PRECISION, &d, pivot, foo);
   osl_int_div_exact(PIPLIB_INT_PRECISION, &lpiv, pivot, d);
   osl_int_div_exact(PIPLIB_INT_PRECISION, &foo, foo, d);
   osl_int_assign(PIPLIB_INT_PRECISION, &d, Denom(tp,k));
   osl_int_mul(PIPLIB_INT_PRECISION, &gcd, lpiv, d);
   osl_int_assign(PIPLIB_INT_PRECISION, &Denom(tp, k), gcd);
   p = tp->row[k].objet.val;
   q = tp->row[pivi].objet.val;
   for(j = 0; j<ncol; j++){
     if(j == pivj)
       osl_int_mul(PIPLIB_INT_PRECISION, &z, dpiv, foo);
     else {
       osl_int_mul(PIPLIB_INT_PRECISION, &z, *p, lpiv);
       osl_int_mul(PIPLIB_INT_PRECISION, &y, *q, foo);
       osl_int_sub(PIPLIB_INT_PRECISION, &z, z, y);
     }
     q++;
     cross_product++;
     osl_int_assign(PIPLIB_INT_PRECISION, p, z);
     p++;
     if (! osl_int_one(PIPLIB_INT_PRECISION, gcd))
       osl_int_gcd(PIPLIB_INT_PRECISION, &gcd, gcd, z);
   }
   if (! osl_int_one(PIPLIB_INT_PRECISION, gcd)) {
     p = tp->row[k].objet.val;
     for(j = 0; j<ncol; j++){
       osl_int_div_exact(PIPLIB_INT_PRECISION, p, *p, gcd);
       p++;
     }
   }
   osl_int_div_exact(PIPLIB_INT_PRECISION, &Denom(tp,k), Denom(tp,k), gcd);
 }
 p = tp->row[pivi].objet.val;
 for(k = 0; k<nligne; k++)
   if((Flag(tp, k) & Unit) && tp->row[k].objet.unit == pivj) break;
 Flag(tp, k) = Plus;
 tp->row[k].objet.val = p;
 for(j = 0; j<ncol; j++)
   osl_int_assign(PIPLIB_INT_PRECISION, p++, new[j]);

 osl_int_assign(PIPLIB_INT_PRECISION, &Denom(tp, k), pivot);
 Flag(tp, pivi) = Unit | Zero;
 osl_int_assign(PIPLIB_INT_PRECISION, &Denom(tp, pivi), UN);
 tp->row[pivi].objet.unit = pivj;

 for(k = 0; k<nligne; k++){
   ff = Flag(tp, k);
   if(ff & Unit) continue;

   if (osl_int_pos(PIPLIB_INT_PRECISION, Index(tp, k, pivj))) fff = Plus;
   else if (osl_int_neg(PIPLIB_INT_PRECISION, Index(tp, k, pivj))) fff = Minus;
   else fff = Zero;

   if (fff != Zero && fff != ff) {
     if(ff == Zero) ff = (fff == Minus ? Unknown : fff);
     else ff = Unknown;
   }
   Flag(tp, k) = ff;
 }

 if(verbose>2){
   fprintf(dump, "just pivoted\n");
   tab_display(tp, dump);
 }

osl_int_clear(PIPLIB_INT_PRECISION, &x);
osl_int_clear(PIPLIB_INT_PRECISION, &y);
osl_int_clear(PIPLIB_INT_PRECISION, &d);
osl_int_clear(PIPLIB_INT_PRECISION, &gcd);
osl_int_clear(PIPLIB_INT_PRECISION, &u);
osl_int_clear(PIPLIB_INT_PRECISION, &dpiv);
osl_int_clear(PIPLIB_INT_PRECISION, &lpiv);
osl_int_clear(PIPLIB_INT_PRECISION, &pivot);
osl_int_clear(PIPLIB_INT_PRECISION, &foo);
osl_int_clear(PIPLIB_INT_PRECISION, &z);
osl_int_clear(PIPLIB_INT_PRECISION, &ppivot);
osl_int_clear(PIPLIB_INT_PRECISION, &dppiv);

 for(i=0; i<ncol; i++)
  osl_int_clear(PIPLIB_INT_PRECISION, &new[i]);

 return 0;
}

/*
 * Sort the rows in increasing order of the largest coefficient
 * and (if TRAITER_DUAL is set) return the new position of the
 * original constraints.
 */
static int* tab_sort_rows(Tableau* tp, int nvar, int nligne, int flags)
{
    int i, j;
    int pivi;
    double s, t, d, smax = 0;
    struct L temp;
    int *pos = NULL, *ineq = NULL;

    if (flags & TRAITER_DUAL) {
	ineq = malloc(tp->height * sizeof(int));
	pos = malloc((tp->height-nvar) * sizeof(int));
	if (!ineq || !pos) {
	    fprintf(stderr, "Memory Overflow.\n") ;
	    exit(1) ;
	}
    }

    for (i = nvar; i < nligne; i++) {
	if (Flag(tp,i) & Unit)
	    continue;
	s = 0;
	d = osl_int_get_si(PIPLIB_INT_PRECISION, Denom(tp, i));
	for (j = 0; j < nvar; j++) {
	    t = osl_int_get_si(PIPLIB_INT_PRECISION, Index(tp,i,j))/d;
	    s = OSL_max(s, abs(t));
	}
	tp->row[i].size = s;
	smax = OSL_max(s, smax);
	if (flags & TRAITER_DUAL)
	    ineq[i] = i-nvar;
    }

    for (i = nvar; i < nligne; i++) {
	if (Flag(tp,i) & Unit)
	    continue;
	s = smax;
	pivi = i;
	for (j = i; j < nligne; j++) {
	    if (Flag(tp,j) & Unit)
		continue;
	    if (tp->row[j].size < s) {
		s = tp->row[j].size;
		pivi = j;
	    }
	}
	if (pivi != i) {
	    temp = tp->row[pivi];
	    tp->row[pivi] = tp->row[i];
	    tp->row[i] = temp;
	    if (flags & TRAITER_DUAL) {
		j = ineq[i];
		ineq[i] = ineq[pivi];
		ineq[pivi] = j;
	     }
	}
    }

    if (flags & TRAITER_DUAL) {
	for (i = nvar; i < nligne; i++)
	    pos[ineq[i]] = i;
	free(ineq);
    }

    return pos;
}

/* dans cette version, "traiter" modifie ineq;
   par contre le contexte est immediatement recopie' */
void traiter(
  Tableau* tp, Tableau *ctxt, int nvar, int nparm, int ni, int nc,
  int bigparm, int flags
)
{
 int j;
 int pivi, nligne, ncol;
 struct high_water_mark x;
 Tableau *context;
 int dch, dcw;
 int *pos;

 dcw = osl_int_size_in_base_2(PIPLIB_INT_PRECISION, tp->determinant);
 dch = 2 * dcw + 1;
 x = tab_hwm();
 nligne = nvar+ni;

 context = expanser(ctxt, 0, nc, nparm+1, 0, dch, dcw);

 pos = tab_sort_rows(tp, nvar, nligne, flags);

 for(;;) {
   if(verbose>2){
     fprintf(dump, "debut for\n");
     tab_display(tp, dump);
     fflush(dump);
   }
   nligne = nvar+ni; ncol = nvar+nparm+1;
   if(nligne > tp->height || ncol > tp->width) {
     fprintf(stdout, "Syserr : traiter : tableau too small\n");
     exit(1);
   }
   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;	       /* There is a negative row   */
   
   pivi = exam_coef(tp, nvar, ncol, bigparm);

   if(verbose>2){
     fprintf(dump, "coefs examined\n");
     tab_display(tp, dump);
     fflush(dump);
   }

   if(pivi < nligne) goto pirouette;
   /* There is a row whose coefficients are negative */
   compa_test(tp, context, ni, nvar, nparm, nc);
   if(verbose>2){
     fprintf(dump, "compatibility tested\n");
     tab_display(tp, dump);
     fflush(dump);
   }

   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;
   /* The compatibility test has found a negative row */
   pivi = chercher(tp, Critic, nligne);
   if(pivi >= nligne)pivi = chercher(tp, Unknown, nligne);
   /* Here, the problem tree splits        */
   if(pivi < nligne) {
     Tableau * ntp;
     osl_int_t com_dem;
     struct high_water_mark q;
     if(nc >= context->height) {
       dcw = osl_int_size_in_base_2(PIPLIB_INT_PRECISION, context->determinant);
       dch = 2 * dcw + 1;
       context = expanser(context, 0, nc, nparm+1, 0, dch, dcw);
     }
     if(nparm >= MAXPARM) {
       fprintf(stdout, "Too much parameters : %d\n", nparm);
       exit(2);
     }
     q = tab_hwm();
     if(verbose>1)
       fprintf(stdout,"profondeur %d %p\n", profondeur, q.top);
     ntp = expanser(tp, nvar, ni, ncol, 0, 0, 0);
     fflush(stdout);
     sol_if();
     sol_forme(nparm+1);
     osl_int_init_set_si(PIPLIB_INT_PRECISION, &com_dem, 0);
     for (j = 0; j < nparm; j++)
       osl_int_gcd(PIPLIB_INT_PRECISION, &com_dem, com_dem, Index(tp, pivi, j + nvar +1));
     if (!(flags & TRAITER_INT))
	 osl_int_gcd(PIPLIB_INT_PRECISION, &com_dem, com_dem, Index(tp, pivi, nvar));
     for (j = 0; j < nparm; j++) {
       osl_int_div_exact(PIPLIB_INT_PRECISION, &Index(context, nc, j), Index(tp, pivi, j + nvar + 1), com_dem);
       sol_val(Index(context, nc, j), UN);
     }
     if (!(flags & TRAITER_INT))
	 osl_int_div_exact(PIPLIB_INT_PRECISION, &Index(context, nc, nparm), Index(tp, pivi, nvar), com_dem);
     else
	 osl_int_floor_div_q(PIPLIB_INT_PRECISION, &Index(context, nc, nparm), Index(tp, pivi, nvar), com_dem);
     sol_val(Index(context, nc, nparm), UN);
     osl_int_clear(PIPLIB_INT_PRECISION, &com_dem);
     Flag(context, nc) = Unknown;
     osl_int_set_si(PIPLIB_INT_PRECISION, &Denom(context, nc), 1);
     Flag(ntp, pivi) = Plus;
     profondeur++;
     fflush(stdout);
     if(verbose > 0) fflush(dump);
     traiter(ntp, context, nvar, nparm, ni, nc+1, bigparm, flags);
     profondeur--;
     tab_reset(q);
     if(verbose>1)
       fprintf(stdout, "descente %d %p\n", profondeur, tab_hwm().top);
     for(j = 0; j<nparm; j++)
       osl_int_oppose(PIPLIB_INT_PRECISION, &Index(context, nc, j), Index(context, nc, j));
     osl_int_increment(PIPLIB_INT_PRECISION, &Index(context, nc, nparm), Index(context, nc, nparm));
     osl_int_oppose(PIPLIB_INT_PRECISION, &Index(context, nc, nparm), Index(context, nc, nparm));
     Flag(tp, pivi) = Minus;
     osl_int_assign(PIPLIB_INT_PRECISION, &Denom(context, nc), UN);
     nc++;
     goto pirouette;
   }
/* Here, all rows are positive. Do we need an integral solution?      */
   if (!(flags & TRAITER_INT)) {
     solution(tp, nvar, nparm);
     if (flags & TRAITER_DUAL)
	solution_dual(tp, nvar, /*nparm,*/ pos);
     break;
   }
/* Yes we do! */
   pivi = integrer(&tp, &context, &nvar, &nparm, &ni, &nc, bigparm);
   if(pivi > 0) goto pirouette;
		    /* A cut has been inserted and is always negative */
/* Here, either there is an integral solution, */
   if(pivi == 0) solution(tp, nvar, nparm);
/* or no solution exists */
   else sol_nil();
   break;

/* Here, a negative row has been found. The call to <<pivoter>> executes
      a pivoting step                                                 */

pirouette :
     if (pivoter(tp, pivi, nvar, nparm, ni) < 0) {
       sol_nil();
       break;
     }
 }
/* Danger : a premature return would induce memory leaks   */
 tab_reset(x);
 free(pos);
 return;
}
