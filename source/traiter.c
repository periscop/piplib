/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 traiter.h                                  *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.							              *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <piplib/piplib.h>

#define max(x,y) ((x) > (y)? (x) : (y))

extern long int cross_product, limit;
extern int verbose;
extern FILE *dump;
extern int profondeur;
extern int compa_count;

int llog(Entier x)
{int n = 0;
/* x must be positive, you dummy */
 if(x<0) x=-x;
 while(x) x >>= 1, n++;
 return(n);
}

int chercher(Tableau *p, int masque, int n)
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

Tableau *expanser(Tableau *tp, int virt, int reel, int ncol, 
                               int off, int dh, int dw)
{
 int i, j, ff;
 char *q; Entier *pq;
 Entier *pp, *qq;
 Tableau *rp;
 if(tp == NULL) return(NULL);
 rp = tab_alloc(reel+dh, ncol+dw, virt);
 rp->l_determinant = tp->l_determinant;
 for(i=0; i<tp->l_determinant; i++)
     rp->determinant[i] = tp->determinant[i];
 pq = (Entier *) & (rp->row[virt+reel+dh]);
 for(i = off; i<virt + reel; i++)
     {ff = Flag(rp, i) = Flag(tp, i-off);
      Denom(rp, i) = Denom(tp, i-off);
      if(ff & Unit) rp->row[i].objet.unit = tp->row[i-off].objet.unit;
      else {
	  rp->row[i].objet.val = pq;
	  pq +=(ncol + dw);
	  pp = tp->row[i-off].objet.val;
	  qq = rp->row[i].objet.val;
	  for(j = 0; j<ncol; j++) *qq++ = *pp++;
	  }
      }
 return(rp);
}

int exam_coef(Tableau *tp, int nvar, int ncol, int bigparm)
{int i, j;
 int ff, fff;
 Entier x, *p;
 for(i = 0; i<tp->height; i++)
     {ff = Flag(tp,i);
      if(ff == 0) break;
      if(ff == Unknown)
	  {if(bigparm >= 0 && (x = Index(tp,i, bigparm)))
	      {if(x<0) {Flag(tp, i) = Minus;
			return(i);
		       }
	       else    Flag(tp, i) = Plus;
	       continue;
	      }
	   ff = Zero;
	   p = &(tp->row[i].objet.val[nvar+1]);
	   for(j = nvar+1; j<ncol; j++)
	       {x = *p++;
		if(x<0) fff = Minus;
		else if (x>0) fff = Plus;
		else fff = Zero;
		if(fff != Zero && fff != ff)
		    if(ff == Zero) ff = fff;
		    else {ff = Unknown;
			  break;
			 }
	       }
/* bug de'tecte' par [paf], 16/2/93 !
   Si tous les coefficients des parame`tres sont ne'gatifs
   et si le terme constant est nul, le signe est inconnu!!
   On traite donc spe'cialement le terme constant. */
	   x = Index(tp, i, nvar);
	   if(x<0) fff = Minus;
	   else if(x>0) fff = Plus;
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
   sauf si le terme constant est ‚galement n‚gatif. */
	   case Minus: if(fff != Minus) ff = Unknown; break;
/* enfin, il n'y a rien a` dire si le signe des coefficients est inconnu */
	   }
	   Flag(tp, i) = ff;
	   if(ff == Minus) return(i);
	  }
      }
 return(i);
}

void compa_test(Tableau *tp, Tableau *context,
		int ni, int nvar, int nparm, int nc)
{
 Entier discr[MAXPARM];
 int i, j;
 int ff;
 int cPlus, cMinus, isCritic;
 int verbold;
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
	  {isCritic = True;
	   for(j = 0; j<nvar; j++) if(Index(tp, i, j) > 0)
		 {isCritic = False;
		  break;
		 }
           compa_count++;
	   for(j = 0; j < nparm; j++) discr[j] = Index(tp, i, j+nvar+1);
	   discr[nparm] = Index(tp, i, nvar)- (isCritic ? 0 : 1);
           /* NdCed : Attention au contexte == NULL ! */
	   tPlus = expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tPlus, nparm+nc) = Unknown;
	   for(j = 0; j<=nparm; j++)Index(tPlus, nparm+nc, j) = discr[j];
	   Denom(tPlus, nparm+nc) = UN;
	   
	   p = sol_hwm();
	   traiter(tPlus, NULL, True, UN, nparm, 0, nc+1, 0, -1);
	   cPlus = is_not_Nil(p);
	   if(verbose>0){
	     fprintf(dump, "\nThe positive case has been found ");
	     fprintf(dump, cPlus? "possible\n": "impossible\n");
	     fflush(dump);
	   }

	   sol_reset(p);
	   for(j = 0; j<nparm+1; j++) discr[j] = -discr[j];
	   discr[nparm] = discr[nparm] - (isCritic ? 1 : 2);
	   tMinus = expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tMinus, nparm+nc) = Unknown;
	   for(j = 0; j<= nparm; j++)Index(tMinus, nparm+nc, j) = discr[j];
	   Denom(tMinus, nparm+nc) = UN;
	   traiter(tMinus, NULL, True, UN, nparm, 0, nc+1, 0, -1);
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

Entier valeur(Tableau *tp, int i, int j, Entier D)
{
 if(Flag(tp, i) & Unit)
     return(tp->row[i].objet.unit == j ? Denom(tp,i) : 0);
 else return(Index(tp, i, j));
}

void solution(Tableau *tp, int nvar, int nparm, Entier D)
{int i, j;
 int ncol = nvar + nparm + 1;
 sol_list(nvar);
 for(i = 0; i<nvar; i++)
     {sol_forme(nparm+1);
      for(j = nvar+1; j<ncol; j++)
	 sol_val(valeur(tp, i, j, D), Denom(tp,i));
      sol_val(valeur(tp, i, nvar, D), Denom(tp,i));
     }
}

int choisir_piv(Tableau *tp, int pivi, int nvar, int nligne, Entier D)
{
 int j, k;
 Entier pivot, foo, x;
 int pivj = -1;
 for(j = 0; j<nvar; j++)
     {if((foo = Index(tp, pivi, j)) <= 0) continue;
      if(pivj < 0)
	  {pivj = j;
	   pivot = foo;
	   continue;
	  }
      for(k = 0; k<nligne; k++)
	  {x = pivot * valeur(tp, k, j, D) - valeur(tp, k, pivj, D) * foo;
	   cross_product++;
	   if(x) break;
	  }
      if(x < 0)
	  {pivj = j;
	   pivot = foo;
	  }
     }
 return(pivj);
}

Entier pivoter(Tableau *tp, int pivi, Entier D, int nvar,
	       int nparm, int ni, int iq)
{int pivj;
 int ncol = nvar + nparm + 1;
 int nligne = nvar + ni;
 int i, j, k;
 Entier x, d, gcd, u, dpiv;
 int ff, fff;
 Entier pivot, foo, z;
 Entier ppivot, dppiv;
 Entier new[MAXCOL], *p, *q;
 char format_format[32];

 sprintf(format_format, "\nPivot %s/%s\n", FORMAT, FORMAT);
 if(ncol >= MAXCOL) {
     fprintf(stderr, "Too much variables\n");
     exit(1);
     }
 if(0 > pivi || pivi >= nligne || Flag(tp, pivi) == Unit) {
     fprintf(stderr, "Syserr : pivoter : wrong pivot row\n");
     exit(1);
     }
 pivj = choisir_piv(tp, pivi, nvar, nligne, D);
 if(pivj < 0) return(-1);
 if(pivj >= nvar) {
     fprintf(stderr, "Syserr : pivoter : wrong pivot\n");
     exit(1);
     }
 pivot = Index(tp, pivi, pivj);
 dpiv = Denom(tp, pivi);
 d = pgcd(pivot, dpiv);
 ppivot = pivot/d;
 dppiv = dpiv/d;
 if(verbose>0){
     fprintf(dump, format_format, ppivot, dppiv);
     fprintf(dump, "%d x %d\n", pivi, pivj);
     }
 for(i=0; i< tp->l_determinant; i++){
     d=pgcd(tp->determinant[i], dppiv);
     tp->determinant[i] /= d;
     dppiv /= d;
     }
 if(dppiv != 1) {
     fprintf(stderr, "Integer overflow : %d\n", D);
/*   tab_display(tp, stdout); */
     if(verbose>0) fflush(dump);
     exit(1);
     }
 for(i=0; i<tp->l_determinant; i++)
     if(llog(tp->determinant[i]) + llog(ppivot) < 8*sizeof(Entier)){
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
 if(verbose>0){
    fprintf(dump, "determinant ");
    for(i=0; i<tp->l_determinant; i++)
	fprintf(dump, FORMAT, tp->determinant[i]);
    fprintf(dump, "\n");
    }
 for(j = 0; j<ncol; j++)
      new[j] = (j == pivj ? dpiv : -Index(tp, pivi, j));
 for(k = 0; k<nligne; k++)
     {Entier lpiv;
      if(Flag(tp,k) & Unit)continue;
      if(k == pivi)continue;
      foo = Index(tp, k, pivj);
      d = pgcd(pivot, foo);
      lpiv = pivot/d;
      foo /= d;
      d = Denom(tp,k);
      gcd = lpiv * d;
      Denom(tp, k) = gcd;
      p = tp->row[k].objet.val;
      q = tp->row[pivi].objet.val;
      for(j = 0; j<ncol; j++)
	   {if(j == pivj)
		z = dpiv * foo;
	    else
		z = (*p) * lpiv - (*q) * foo;
	    q++;
	    cross_product++;
	    *p++ = z;
	    if(gcd != 1) gcd = pgcd(gcd, z);
	   }
       if(gcd != 1) {
	   p = tp->row[k].objet.val;
	   for(j = 0; j<ncol; j++)
	       *p++ /= gcd;
	   }
       Denom(tp,k) = Denom(tp,k)/gcd;
      }
 p = tp->row[pivi].objet.val;
 for(k = 0; k<nligne; k++)
     if((Flag(tp, k) & Unit) && tp->row[k].objet.unit == pivj) break;
 Flag(tp, k) = Plus;
 tp->row[k].objet.val = p;
 for(j = 0; j<ncol; j++)
      *p++ = new[j];
 Denom(tp, k) = pivot; 
 Flag(tp, pivi) = Unit | Zero;
 Denom(tp, pivi) = UN;
 tp->row[pivi].objet.unit = pivj;

 for(k = 0; k<nligne; k++)
    {ff = Flag(tp, k);
     if(ff & Unit) continue;
     x = Index(tp, k, pivj);
     if(x < 0) fff = Minus;
     else if(x == 0) fff = Zero;
     else fff = Plus;
     if(fff != Zero && fff != ff)
	 if(ff == Zero) ff = (fff == Minus ? Unknown : fff);
	 else ff = Unknown;
     Flag(tp, k) = ff;
    }
/*
 if(verbose>0){
     for(k=0; k<nligne; k++){
	 if(Flag(tp,k) & Unit)
	     fprintf(dump, "1 ");
	 else
	     fprintf(dump, FORMAT, Denom(tp,k));
	 if(k%10 == 9) putc('\n', dump);
	 else putc(' ', dump);
	 }
     putc('\n', dump);
     }

 if(verbose>0){
     for(i=0; i<nligne; i++){
         if(Flag(tp, i) & Unit)
               fprintf(dump, "0 ");
         else {
             fprintf(dump, FORMAT, Index(tp,i, nvar));
             putc('/', dump);
             fprintf(dump, FORMAT, Denom(tp,i));
             }
         if((i + 1) % 10 == 0) putc('\n', dump);
         else putc(' ', dump);
         }
     putc('\n', dump);
     }
*/
 if(verbose>0){
   fprintf(dump, "just pivoted\n");
   tab_display(tp, dump);
 }

 return(D);
}

/* dans cette version, "traiter" modifie ineq; par contre
   le contexte est immediatement recopie' */

static Entier discr[MAXPARM];

Entier traiter(tp, ctxt, iq, D, nvar, nparm, ni, nc, bigparm)
Tableau *tp, *ctxt;
int iq, nvar, nparm, ni, nc, bigparm;
Entier D;
{
 int j;
 int pivi, nligne, ncol;
 struct high_water_mark x;
 Tableau *context;
 int dch, dcw;
 float s, t, d, smax;
 int i;
 struct L temp;
 
 dcw = llog(D);
 dch = 2 * dcw + 1;
 x = tab_hwm();
 nligne = nvar+ni;

 context = expanser(ctxt, 0, nc, nparm+1, 0, dch, dcw);
/*
 sort the rows in increasing order of the largest coefficient
*/
 smax = 0.;

 for(i=nvar; i<nligne; i++){
     if(Flag(tp,i) & Unit) continue;
     s = 0.;
     d = (float) Denom(tp,i);     
     for(j=0; j<nvar; j++){
         t = Index(tp,i,j)/d;
         s = max(s, abs(t));
         }
     tp->row[i].size = s;
     smax = max(s, smax);
     }
    
 for(i=nvar; i<nligne; i++){
      if(Flag(tp,i) & Unit) continue;
      s = smax;
      pivi = i;
      for(j=i; j<nligne; j++){
          if(Flag(tp,j) & Unit) continue;
          if(tp->row[j].size < s){
              s = tp->row[i].size;
              pivi = j;
              }
          }
       if(pivi != i) {
           temp = tp->row[pivi];
           tp->row[pivi] = tp->row[i];
           tp->row[i]=temp;
           }
       }   
 
 if(verbose>0){
   fprintf(dump, "just sorted\n");
   tab_display(tp, dump);
 }
      
 for(;;) {
     nligne = nvar+ni; ncol = nvar+nparm+1;
     if(nligne > tp->height || ncol > tp->width) {
	 fprintf(stderr, "Syserr : traiter : tableau too small\n");
	 exit(1);
	 }
     pivi = chercher(tp, Minus, nligne);
     if(pivi < nligne) goto pirouette;	       /* There is a negative row   */

     pivi = exam_coef(tp, nvar, ncol, bigparm);

     if(verbose > 0){
       fprintf(dump, "coefs examined\n");
       tab_display(tp, dump);
     }

     if(pivi < nligne) goto pirouette;
			/* There is a row whose coefficients are negative */
     compa_test(tp, context, ni, nvar, nparm, nc);
     if(verbose>0){
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
	 Entier com_dem;
	 struct high_water_mark q;
	 if(nc >= context->height) {
	     dcw = llog(D);
	     dch = 2 * dcw + 1;
 	     context = expanser(context, 0, nc, nparm+1, 0, dch, dcw);
	     }
	 if(nparm >= MAXPARM) {
	     fprintf(stderr, "Too much parameters : %d\n", nparm);
	     exit(2);
	     }
	 q = tab_hwm();
	 if(verbose>0)
	     fprintf(stderr,"profondeur %d %lx\n", profondeur, q.top);
	 ntp = expanser(tp, nvar, ni, ncol, 0, 0, 0);
	 sol_if();
	 sol_forme(nparm+1);
	 com_dem = 0;
	 for(j = 0; j<nparm; j++) {
	     discr[j] = Index(tp, pivi, j + nvar +1);
	     com_dem = pgcd(com_dem, discr[j]);
	     }
	 discr[nparm] = Index(tp, pivi, nvar);
	 com_dem = pgcd(com_dem, discr[nparm]);
	 for(j = 0; j<=nparm; j++) {
	     discr[j] /= com_dem;
	     Index(context, nc, j) = discr[j];
	     sol_val(discr[j], UN);
	     }
	 Flag(context, nc) = Unknown;
	 Denom(context, nc) = UN;
	 Flag(ntp, pivi) = Plus;
	 profondeur++;
        traiter(ntp, context, iq, D, nvar, nparm, ni, nc+1, bigparm);
 	 profondeur--;
	 tab_reset(q);
	 if(verbose>0)
	     fprintf(stderr, "descente %d %lx\n", profondeur, tab_hwm().top);
	 for(j = 0; j<nparm; j++)
	     Index(context, nc, j) = - Index(context, nc, j);
	 Index(context, nc, nparm) = - Index(context, nc, nparm) -1;
	 Flag(tp, pivi) = Minus;
	 Denom(context, nc) = UN;
	 nc++;
	 goto pirouette;
	 }
/* Here, all rows are positive. Do we need an integral solution?      */
    if(!iq) {
	solution(tp, nvar, nparm, D);
	break;
	}
/* Yes we do! */
    pivi = integrer(&tp, &context, D, &nvar, &nparm, &ni, &nc);
    if(pivi > 0) goto pirouette;
		    /* A cut has been inserted and is always negative */
/* Here, either there is an integral solution, */
    if(pivi == 0) solution(tp, nvar, nparm, D);
/* or no solution exists */
    else sol_nil();
    break;

/* Here, a negative row has been found. The call to <<pivoter>> executes
      a pivoting step                                                 */

pirouette :
    if((D = pivoter(tp, pivi, D, nvar, nparm, ni, iq)) < 0) {
	sol_nil();
	break;
	}
    }
/* Danger : a premature return would induce memory leaks   */
 tab_reset(x);
 return D;
}







