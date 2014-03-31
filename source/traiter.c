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

#include "pip.h"


/*extern long int cross_product;*/
extern int PIPLIB_NAME(verbose);
extern FILE *PIPLIB_NAME(dump);
/*extern int compa_count;*/

int chercher(PIPLIB_NAME(Tableau) *p, int masque, int n)
{int i;
 for(i = 0; i<n; i++)
     if(p->row[i].flags & masque) break;
 return(i);
}

/* il est convenu que traiter ne doit modifier ni le tableau, ni le contexte;
   le tableau peut grandir en cas de coupure (+1 en hauteur et +1 en largeur
   si nparm != 0) et en cas de partage (+1 en hauteur)(seulement si nparm != 0).
   le contexte peut grandir en cas de coupure (+2 en hauteur et +1 en largeur)
   (seulement si nparm !=0) et en cas de partage (+1 en hauteur)(nparm !=0).
   On estime le nombre de coupures a llog(D) et le nombre de partages a
   ni.
*/

PIPLIB_NAME(Tableau) *PIPLIB_NAME(expanser)(PIPLIB_NAME(Tableau) *tp, int virt, int reel, int ncol, 
                               int off, int dh, int dw)
{
 int i, j, ff;
 PIPLIB_NAME(piplib_int_t) *pq;
 PIPLIB_NAME(piplib_int_t) *pp, *qq;
 PIPLIB_NAME(Tableau) *rp;
 if(tp == NULL) return(NULL);
 rp = PIPLIB_NAME(tab_alloc)(reel+dh, ncol+dw, virt);

 #if defined(PIPLIB_ONE_DETERMINANT)
 piplib_int_assign(rp->determinant, tp->determinant);
 #else
 rp->l_determinant = tp->l_determinant;
 for(i=0; i<tp->l_determinant; i++)
     rp->determinant[i] = tp->determinant[i];
 #endif
 pq = (PIPLIB_NAME(piplib_int_t) *) & (rp->row[virt+reel+dh]);
 for(i = off; i<virt + reel; i++)
     {ff = Flag(rp, i) = Flag(tp, i-off);
      piplib_int_assign(Denom(rp, i), Denom(tp, i-off));
      if(ff & Unit) rp->row[i].objet.unit = tp->row[i-off].objet.unit;
      else {
	  rp->row[i].objet.val = pq;
	  pq +=(ncol + dw);
	  pp = tp->row[i-off].objet.val;
	  qq = rp->row[i].objet.val;
	  for(j = 0; j<ncol; j++)
	     piplib_int_assign(*qq++, *pp++);
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
 * of the big parameter first.
 */
int PIPLIB_NAME(exam_coef)(PIPLIB_NAME(Tableau) *tp, int nvar, int ncol, int bigparm)
{int i, j ;
 int ff, fff;
 PIPLIB_NAME(piplib_int_t) *p;
 
 if (bigparm >= 0)
    for (i = 0; i<tp->height; i++) {
	if (Flag(tp, i) != Unknown)
	    continue;
	if (piplib_int_neg(Index(tp,i, bigparm))) {
	    Flag(tp, i) = Minus;
	    return i;
	} else if (piplib_int_pos(Index(tp,i, bigparm)))
	    Flag(tp, i) = Plus;
    }

 for(i = 0; i<tp->height; i++)
     {ff = Flag(tp,i);
      if(ff == 0) break;
      if(ff == Unknown) {
	   ff = Zero;
	   p = &(tp->row[i].objet.val[nvar+1]);
	   for(j = nvar+1; j<ncol; j++) {
		if (piplib_int_neg(*p)) fff = Minus;
		else if (piplib_int_pos(*p)) fff = Plus;
		else fff = Zero;
		p++;
		if(fff != Zero && fff != ff) {
		    if(ff == Zero) { ff = fff; }
		    else { ff = Unknown; break; }
		}
	       }
/* bug de'tecte' par [paf], 16/2/93 !
   Si tous les coefficients des parame`tres sont ne'gatifs
   et si le terme constant est nul, le signe est inconnu!!
   On traite donc spe'cialement le terme constant. */
	   if(piplib_int_neg(Index(tp, i, nvar))) fff = Minus;
	   else if(piplib_int_pos(Index(tp, i, nvar))) fff = Plus;
	   else fff = Zero;
/* ici on a le signe du terme constant */
	   switch(ff){
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
	   if(ff == Minus) return(i);
	  }
      }
 return(i);
}

void PIPLIB_NAME(compa_test)(PIPLIB_NAME(Tableau) *tp, PIPLIB_NAME(Tableau) *context,
		int ni, int nvar, int nparm, int nc)
{
 int i, j;
 int ff;
 int cPlus, cMinus, isCritic;
 PIPLIB_NAME(Tableau) *tPlus, *tMinus;
 int p;
 struct PIPLIB_NAME(high_water_mark) q;

 if(nparm == 0) return;
 if(nparm >= MAXPARM) {
     fprintf(stderr, "Too much parameters : %d\n", nparm);
     exit(1);
     }
 q = PIPLIB_NAME(tab_hwm)();

 for(i = 0; i<ni + nvar; i++)
     {ff = Flag(tp,i);
      if(ff & (Critic | Unknown))
	  {isCritic = Pip_True;
	   for(j = 0; j<nvar; j++) {
	         if(piplib_int_pos(Index(tp, i, j)))
		 {isCritic = Pip_False;
		  break;
		 }
	   }
           /*compa_count++;*/
	   tPlus = PIPLIB_NAME(expanser)(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tPlus, nparm+nc) = Unknown;
	   for (j = 0; j < nparm; j++)
	       piplib_int_assign(Index(tPlus, nparm+nc, j), Index(tp, i, j+nvar+1));
	   piplib_int_assign(Index(tPlus, nparm+nc, nparm), Index(tp, i, nvar));
	   if (!isCritic)
	       piplib_int_decrement(Index(tPlus, nparm+nc, nparm),
				    Index(tPlus, nparm+nc, nparm));
	   piplib_int_set_si(Denom(tPlus, nparm+nc), 1);
	   
	   p = PIPLIB_NAME(sol_hwm)();
	   PIPLIB_NAME(traiter)(tPlus, NULL, nparm, 0, nc+1, 0, -1, TRAITER_INT);
	   cPlus = PIPLIB_NAME(is_not_Nil)(p);
	   if(PIPLIB_NAME(verbose)>0){
	     fprintf(PIPLIB_NAME(dump), "\nThe positive case has been found ");
	     fprintf(PIPLIB_NAME(dump), cPlus? "possible\n": "impossible\n");
	     fflush(PIPLIB_NAME(dump));
	   }

	   PIPLIB_NAME(sol_reset)(p);
	   tMinus = PIPLIB_NAME(expanser)(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tMinus, nparm+nc) = Unknown;
	   for (j = 0; j < nparm; j++)
	       piplib_int_oppose(Index(tMinus, nparm+nc, j), Index(tp, i, j+nvar+1));
	   piplib_int_oppose(Index(tMinus, nparm+nc, nparm), Index(tp, i, nvar));
	   piplib_int_decrement(Index(tMinus, nparm+nc, nparm),
				Index(tMinus, nparm+nc, nparm));
	   piplib_int_set_si(Denom(tMinus, nparm+nc), 1);
	   PIPLIB_NAME(traiter)(tMinus, NULL, nparm, 0, nc+1, 0, -1, TRAITER_INT);
	   cMinus = PIPLIB_NAME(is_not_Nil)(p);
	   if(PIPLIB_NAME(verbose)>0){
	     fprintf(PIPLIB_NAME(dump), "\nThe negative case has been found ");
	     fprintf(PIPLIB_NAME(dump), cMinus? "possible\n": "impossible\n");
	     fflush(PIPLIB_NAME(dump));
	   }

	   PIPLIB_NAME(sol_reset)(p);
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
 PIPLIB_NAME(tab_reset)(q);
 
 return;
}

PIPLIB_NAME(piplib_int_t) *PIPLIB_NAME(valeur)(PIPLIB_NAME(Tableau) *tp, int i, int j, PIPLIB_NAME(piplib_int_t)* zero)
{
 if(Flag(tp, i) & Unit)
     return(tp->row[i].objet.unit == j ? &Denom(tp,i) : zero);
 else return(&Index(tp, i, j));
}

void solution(PIPLIB_NAME(Tableau) *tp, int nvar, int nparm)
{int i, j;
 int ncol = nvar + nparm + 1;

 PIPLIB_NAME(piplib_int_t) zero;
 piplib_int_init_set_si(zero, 0);

 PIPLIB_NAME(sol_list)(nvar);
 for(i = 0; i<nvar; i++)
     {PIPLIB_NAME(sol_forme)(nparm+1);
      for(j = nvar+1; j<ncol; j++)
	 PIPLIB_NAME(sol_val)(*PIPLIB_NAME(valeur)(tp, i, j, &zero), Denom(tp,i));
      PIPLIB_NAME(sol_val)(*PIPLIB_NAME(valeur)(tp, i, nvar, &zero), Denom(tp,i));
     }

 piplib_int_clear(zero);
}

static void PIPLIB_NAME(solution_dual)(PIPLIB_NAME(Tableau) *tp, int nvar/*, int nparm*/, int *pos)
{
    int i;

    PIPLIB_NAME(piplib_int_t) zero;
    piplib_int_init_set_si(zero, 0);

    PIPLIB_NAME(sol_list)(tp->height - nvar);
    for (i = 0; i < tp->height - nvar; ++i) {
	PIPLIB_NAME(sol_forme)(1);
	if (Flag(tp, pos[i]) & Unit)
	    PIPLIB_NAME(sol_val)(*PIPLIB_NAME(valeur)(tp, 0, tp->row[pos[i]].objet.unit, &zero), Denom(tp, 0));
	else
	    PIPLIB_NAME(sol_val_zero_one)();
    }

    piplib_int_clear(zero);
}

int PIPLIB_NAME(choisir_piv)(PIPLIB_NAME(Tableau) *tp, int pivi, int nvar, int nligne)
{
 int j, k;
 PIPLIB_NAME(piplib_int_t) pivot, foo, x, y;
 int pivj = -1;

 PIPLIB_NAME(piplib_int_t) zero;
 piplib_int_init_set_si(zero, 0);

 piplib_int_init(pivot);
 piplib_int_init(foo);
 piplib_int_init(x);
 piplib_int_init(y);
 
 for(j = 0; j<nvar; j++) {
    piplib_int_assign(foo, Index(tp, pivi, j));
    if(piplib_int_pos(foo) == 0) continue;
    if(pivj < 0)
	{pivj = j;
         piplib_int_assign(pivot, foo);
	 continue;
	}
    for(k = 0; k<nligne; k++)
        {piplib_int_mul(x, pivot, *PIPLIB_NAME(valeur)(tp, k, j, &zero)); 
         piplib_int_mul(y, *PIPLIB_NAME(valeur)(tp, k, pivj, &zero), foo);
         piplib_int_sub(x, x, y);
         /*cross_product++;*/
         if(piplib_int_zero(x) == 0) break;
	}
    if(piplib_int_neg(x))
        {pivj = j;
         piplib_int_assign(pivot, foo);
        }
 }
 
 piplib_int_clear(pivot);
 piplib_int_clear(foo);
 piplib_int_clear(x);
 piplib_int_clear(y);

 piplib_int_clear(zero);

 return(pivj);
}


int PIPLIB_NAME(pivoter)(PIPLIB_NAME(Tableau) *tp, int pivi, int nvar, int nparm, int ni)

{int pivj;
 int ncol = nvar + nparm + 1;
 int nligne = nvar + ni;
 int i, j, k;
 #if defined(PIPLIB_ONE_DETERMINANT)
  PIPLIB_NAME(piplib_int_t) x;
 #endif
 PIPLIB_NAME(piplib_int_t) y, d, gcd, dpiv;
 int ff, fff;
 PIPLIB_NAME(piplib_int_t) pivot, foo, z;
 PIPLIB_NAME(piplib_int_t) ppivot, dppiv;
 PIPLIB_NAME(piplib_int_t) new[MAXCOL], *p, *q;
 PIPLIB_NAME(piplib_int_t) lpiv;

 if(ncol >= MAXCOL) {
   fprintf(stdout, "Too much variables\n");
   exit(1);
 }
 if(0 > pivi || pivi >= nligne || Flag(tp, pivi) == Unit) {
   fprintf(stdout, "Syserr : PIPLIB_NAME(pivoter) : wrong pivot row\n");
   exit(1);
 }

 pivj = PIPLIB_NAME(choisir_piv)(tp, pivi, nvar, nligne);
 if(pivj < 0) return(-1);
 if(pivj >= nvar) {
   fprintf(stdout, "Syserr : PIPLIB_NAME(pivoter) : wrong pivot\n");
   exit(1);
 }

 #if defined(PIPLIB_ONE_DETERMINANT)
  piplib_int_init(x);
 #endif
 piplib_int_init(y); piplib_int_init(d); 
 piplib_int_init(gcd); piplib_int_init(dpiv);
 piplib_int_init(lpiv); piplib_int_init(pivot); piplib_int_init(foo);
 piplib_int_init(z); piplib_int_init(ppivot); piplib_int_init(dppiv);

 for(i=0; i<ncol; i++)
   piplib_int_init(new[i]);

 piplib_int_assign(pivot, Index(tp, pivi, pivj));
 piplib_int_assign(dpiv, Denom(tp, pivi));
 piplib_int_gcd(d, pivot, dpiv);
 piplib_int_div_exact(ppivot, pivot, d);
 piplib_int_div_exact(dppiv, dpiv, d);
 
 if(PIPLIB_NAME(verbose)>1){
   fprintf(PIPLIB_NAME(dump), "Pivot ");
   piplib_int_print(PIPLIB_NAME(dump), ppivot);
   putc('/', PIPLIB_NAME(dump));
   piplib_int_print(PIPLIB_NAME(dump), dppiv);
   putc('\n', PIPLIB_NAME(dump));
   fprintf(PIPLIB_NAME(dump), "%d x %d\n", pivi, pivj);
 }

 #if defined(PIPLIB_ONE_DETERMINANT)
 piplib_int_floor_div_q_r(x, y, tp->determinant, dppiv); 
 #else
 for(i=0; i< tp->l_determinant; i++){
     piplib_int_gcd(d, tp->determinant[i], dppiv);
     tp->determinant[i] /= d;
     dppiv /= d;
     }
 #endif

 #if defined(PIPLIB_ONE_DETERMINANT)
 if (piplib_int_zero(y) == 0) {
 #else
 if(dppiv != 1) {
 #endif
   fprintf(stderr, "Integer overflow\n");
   if(PIPLIB_NAME(verbose)>0) fflush(PIPLIB_NAME(dump));
   exit(1);
 }
 
 #if defined(PIPLIB_ONE_DETERMINANT)
 piplib_int_mul(tp->determinant, x, ppivot);
 #else
 for(i=0; i<tp->l_determinant; i++)
     if(PIPLIB_NAME(piplib_lllog2)(tp->determinant[i]) + PIPLIB_NAME(piplib_lllog2)(ppivot) < 8*sizeof(PIPLIB_NAME(piplib_int_t))){
	 tp->determinant[i] *= ppivot;
	 break;
	 }
 if(i >= tp->l_determinant){
     tp->l_determinant++;
     if(tp->l_determinant >= MAX_DETERMINANT){
	 fprintf(stderr, "Integer overflow : %d\n", tp->l_determinant);
	 exit(1);
	 }
     tp->determinant[i] = ppivot;
     }
 #endif

 if(PIPLIB_NAME(verbose)>1){
   fprintf(PIPLIB_NAME(dump), "determinant ");
   #if defined(PIPLIB_ONE_DETERMINANT)
   piplib_int_print(PIPLIB_NAME(dump), tp->determinant);
   #else
   for(i=0; i<tp->l_determinant; i++)
	fprintf(PIPLIB_NAME(dump), piplib_int_format, tp->determinant[i]);
   #endif
   fprintf(PIPLIB_NAME(dump), "\n");
 }

 
 for(j = 0; j<ncol; j++)
   if(j==pivj)
     piplib_int_assign(new[j], dpiv);
   else 
     piplib_int_oppose(new[j], Index(tp, pivi, j));

 for(k = 0; k<nligne; k++){
   if(Flag(tp,k) & Unit)continue;
   if(k == pivi)continue;
   piplib_int_assign(foo, Index(tp, k, pivj));
   piplib_int_gcd(d, pivot, foo);
   piplib_int_div_exact(lpiv, pivot, d);
   piplib_int_div_exact(foo, foo, d);
   piplib_int_assign(d, Denom(tp,k));
   piplib_int_mul(gcd, lpiv, d);
   piplib_int_assign(Denom(tp, k), gcd);
   p = tp->row[k].objet.val;
   q = tp->row[pivi].objet.val;
   for(j = 0; j<ncol; j++){
     if(j == pivj)
       piplib_int_mul(z, dpiv, foo);
     else {
       piplib_int_mul(z, *p, lpiv);
       piplib_int_mul(y, *q, foo);
       piplib_int_sub(z, z, y);
     }
     q++;
     /*cross_product++;*/
     piplib_int_assign(*p, z);
     p++;
     if (piplib_int_one(gcd) == 0)
       piplib_int_gcd(gcd, gcd, z);
   }
   if(piplib_int_one(gcd) == 0){
     p = tp->row[k].objet.val;
     for(j = 0; j<ncol; j++){
       piplib_int_div_exact(*p, *p, gcd);
       p++;
     }
   }
   piplib_int_div_exact(Denom(tp,k), Denom(tp,k), gcd);
 }
 p = tp->row[pivi].objet.val;
 for(k = 0; k<nligne; k++)
   if((Flag(tp, k) & Unit) && tp->row[k].objet.unit == pivj) break;
 Flag(tp, k) = Plus;
 tp->row[k].objet.val = p;
 for(j = 0; j<ncol; j++)
   piplib_int_assign(*p++, new[j]);

 piplib_int_assign(Denom(tp, k), pivot);
 Flag(tp, pivi) = Unit | Zero;
 piplib_int_set_si(Denom(tp, pivi), 1);
 tp->row[pivi].objet.unit = pivj;

 for(k = 0; k<nligne; k++){
   ff = Flag(tp, k);
   if(ff & Unit) continue;
   if (piplib_int_neg(Index(tp, k, pivj))) fff = Minus;
   else if(piplib_int_zero(Index(tp, k, pivj))) fff = Zero;
   else fff = Plus;
   if (fff != Zero && fff != ff)  {
     if(ff == Zero) { ff = (fff == Minus ? Unknown : fff); }
     else { ff = Unknown; }
   }
   Flag(tp, k) = ff;
 }

 if(PIPLIB_NAME(verbose)>2){
   fprintf(PIPLIB_NAME(dump), "just pivoted\n");
   PIPLIB_NAME(tab_display)(tp, PIPLIB_NAME(dump));
 }

 #if defined(PIPLIB_ONE_DETERMINANT)
  piplib_int_clear(x);
 #endif
 piplib_int_clear(y); piplib_int_clear(d); piplib_int_clear(gcd);
 piplib_int_clear(dpiv); piplib_int_clear(lpiv);
 piplib_int_clear(pivot); piplib_int_clear(foo); piplib_int_clear(z);
 piplib_int_clear(ppivot); piplib_int_clear(dppiv);

 for(i=0; i<ncol; i++)
   piplib_int_clear(new[i]);

 return(0);
}

/*
 * Sort the rows in increasing order of the largest coefficient
 * and (if TRAITER_DUAL is set) return the new position of the
 * original constraints.
 */
static int *PIPLIB_NAME(tab_sort_rows)(PIPLIB_NAME(Tableau) *tp, int nvar, int nligne, int flags)
{
    #define piplib_max(x,y) ((x) > (y)? (x) : (y))
	
    int i, j;
    int pivi;
    double s, t, d, smax = 0;
    struct PIPLIB_NAME(L) temp;
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
	d = piplib_int_get_d(Denom(tp, i));
	for (j = 0; j < nvar; j++) {
	    t = piplib_int_get_d(Index(tp,i,j))/d;
	    s = piplib_max(s, abs(t));
	}
	tp->row[i].size = s;
	smax = piplib_max(s, smax);
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

/* dans cette version, "traiter" modifie ineq; par contre
   le contexte est immediatement recopie' */

void PIPLIB_NAME(traiter)(PIPLIB_NAME(Tableau) *tp, PIPLIB_NAME(Tableau) *ctxt, int nvar, int nparm, int ni, int nc,
	     int bigparm, int flags)
{
 static int profondeur = 1;
 int j;
 int pivi, nligne, ncol;
 struct PIPLIB_NAME(high_water_mark) x;
 PIPLIB_NAME(Tableau) *context;
 int dch, dcw;
 int *pos;
 #if !defined(PIPLIB_ONE_DETERMINANT)
 int i;
 #endif

 #if defined(PIPLIB_ONE_DETERMINANT)
 dcw = piplib_int_size_in_base_2(tp->determinant);
 #else
 dcw = 0;
 for(i=0; i<tp->l_determinant; i++)
   dcw += PIPLIB_NAME(piplib_lllog2)(tp->determinant[i]);
 #endif
 dch = 2 * dcw + 1;
 x = PIPLIB_NAME(tab_hwm)();
 nligne = nvar+ni;

 context = PIPLIB_NAME(expanser)(ctxt, 0, nc, nparm+1, 0, dch, dcw);

 pos = PIPLIB_NAME(tab_sort_rows)(tp, nvar, nligne, flags);

 for(;;) {
   if(PIPLIB_NAME(verbose)>2){
     fprintf(PIPLIB_NAME(dump), "debut for\n");
     PIPLIB_NAME(tab_display)(tp, PIPLIB_NAME(dump));
     fflush(PIPLIB_NAME(dump));
   }
   nligne = nvar+ni; ncol = nvar+nparm+1;
   if(nligne > tp->height || ncol > tp->width) {
     fprintf(stdout, "Syserr : PIPLIB_NAME(traiter) : tableau too small\n");
     exit(1);
   }
   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;	       /* There is a negative row   */
   
   pivi = PIPLIB_NAME(exam_coef)(tp, nvar, ncol, bigparm);

   if(PIPLIB_NAME(verbose)>2){
     fprintf(PIPLIB_NAME(dump), "coefs examined\n");
     PIPLIB_NAME(tab_display)(tp, PIPLIB_NAME(dump));
     fflush(PIPLIB_NAME(dump));
   }

   if(pivi < nligne) goto pirouette;
   /* There is a row whose coefficients are negative */
   PIPLIB_NAME(compa_test)(tp, context, ni, nvar, nparm, nc);
   if(PIPLIB_NAME(verbose)>2){
     fprintf(PIPLIB_NAME(dump), "compatibility tested\n");
     PIPLIB_NAME(tab_display)(tp, PIPLIB_NAME(dump));
     fflush(PIPLIB_NAME(dump));
   }

   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;
   /* The compatibility test has found a negative row */
   pivi = chercher(tp, Critic, nligne);
   if(pivi >= nligne)pivi = chercher(tp, Unknown, nligne);
   /* Here, the problem tree splits        */
   if(pivi < nligne) {
     PIPLIB_NAME(Tableau) * ntp;
     PIPLIB_NAME(piplib_int_t) com_dem;
     struct PIPLIB_NAME(high_water_mark) q;
     if(nc >= context->height) {
       #if defined(PIPLIB_ONE_DETERMINANT)
       dcw = piplib_int_size_in_base_2(context->determinant);
       #else
       dcw = 0;
       for(i=0; i<tp->l_determinant; i++)
       dcw += PIPLIB_NAME(piplib_lllog2)(tp->determinant[i]);
       #endif
       dch = 2 * dcw + 1;
       context = PIPLIB_NAME(expanser)(context, 0, nc, nparm+1, 0, dch, dcw);
     }
     if(nparm >= MAXPARM) {
       fprintf(stdout, "Too much parameters : %d\n", nparm);
       exit(2);
     }
     q = PIPLIB_NAME(tab_hwm)();
     if(PIPLIB_NAME(verbose)>1)
       fprintf(stdout,"profondeur %d %p\n", profondeur, q.top);
     ntp = PIPLIB_NAME(expanser)(tp, nvar, ni, ncol, 0, 0, 0);
     fflush(stdout);
     PIPLIB_NAME(sol_if)();
     PIPLIB_NAME(sol_forme)(nparm+1);
     piplib_int_init(com_dem);
     for (j = 0; j < nparm; j++)
       piplib_int_gcd(com_dem, com_dem, Index(tp, pivi, j + nvar +1));
     if (!(flags & TRAITER_INT))
	 piplib_int_gcd(com_dem, com_dem, Index(tp, pivi, nvar));
     for (j = 0; j < nparm; j++) {
       piplib_int_div_exact(Index(context, nc, j), Index(tp, pivi, j + nvar + 1), com_dem);
       PIPLIB_NAME(sol_val_one)(Index(context, nc, j));
     }
     if (!(flags & TRAITER_INT))
	 piplib_int_div_exact(Index(context, nc, nparm), Index(tp, pivi, nvar), com_dem);
     else
	 piplib_int_floor_div_q(Index(context, nc, nparm), Index(tp, pivi, nvar), com_dem);
     PIPLIB_NAME(sol_val_one)(Index(context, nc, nparm));
     piplib_int_clear(com_dem);
     Flag(context, nc) = Unknown;
     piplib_int_set_si(Denom(context, nc), 1);
     Flag(ntp, pivi) = Plus;
     profondeur++;
     fflush(stdout);
     if(PIPLIB_NAME(verbose) > 0) fflush(PIPLIB_NAME(dump));
     PIPLIB_NAME(traiter)(ntp, context, nvar, nparm, ni, nc+1, bigparm, flags);
     profondeur--;
     PIPLIB_NAME(tab_reset)(q);
     if(PIPLIB_NAME(verbose)>1)
       fprintf(stdout, "descente %d %p\n", profondeur, PIPLIB_NAME(tab_hwm)().top);
     for(j = 0; j<nparm; j++)
       piplib_int_oppose(Index(context, nc, j), Index(context, nc, j));
     piplib_int_increment(Index(context, nc, nparm), Index(context, nc, nparm));
     piplib_int_oppose(Index(context, nc, nparm), Index(context, nc, nparm));
     Flag(tp, pivi) = Minus;
     piplib_int_set_si(Denom(context, nc), 1);
     nc++;
     goto pirouette;
   }
/* Here, all rows are positive. Do we need an integral solution?      */
   if (!(flags & TRAITER_INT)) {
     solution(tp, nvar, nparm);
     if (flags & TRAITER_DUAL)
	PIPLIB_NAME(solution_dual)(tp, nvar/*, nparm*/, pos);
     break;
   }
/* Yes we do! */
   pivi = PIPLIB_NAME(integrer)(&tp, &context, &nvar, &nparm, &ni, &nc, bigparm);
   if(pivi > 0) goto pirouette;
		    /* A cut has been inserted and is always negative */
/* Here, either there is an integral solution, */
   if(pivi == 0) solution(tp, nvar, nparm);
/* or no solution exists */
   else PIPLIB_NAME(sol_nil)();
   break;

/* Here, a negative row has been found. The call to <<pivoter>> executes
      a pivoting step                                                 */

pirouette :
     if (PIPLIB_NAME(pivoter)(tp, pivi, nvar, nparm, ni) < 0) {
       PIPLIB_NAME(sol_nil)();
       break;
     }
 }
/* Danger : a premature return would induce memory leaks   */
 PIPLIB_NAME(tab_reset)(x);
 free(pos);
 return;
}
