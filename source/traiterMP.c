/********************************************************/
/* Part of MultiPrecision PIP port by Zbigniew Chamski  */
/* and Paul Feautrier.                                  */
/* Based on straight PIP E.1 version by Paul Feautrier  */
/* <Paul.Feautrier@inria.fr>                            */
/*                                                      */
/* and a previous port (C) Zbigniew CHAMSKI, 1993.      */
/* <Zbigniew.Chamski@philips.com>                       */
/*                                                      */
/* Copying subject to the terms and conditions of the   */
/* GNU General Public License.                          */
/*                                                      */
/* Send questions, bug reports and fixes to:            */
/* <Paul.Feautrier@inria.fr>                            */
/********************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <piplib/piplib.h>

#define max(x,y) ((x) > (y)? (x) : (y))

extern long int cross_product, limit;
extern int verbose;
extern FILE *dump;
extern int profondeur;
extern float clock;
extern int compa_count;

/* The function llog is replaced by mpz_sizeinbase. */

int chercher(Tableau *p, int masque, int n)
{int i;
 for(i = 0; i<n; i++)
     if(p->row[i].flags & masque) break;
 return(i);
}

/* il est convenu que traiter ne doit modifier ni le tableau, ni le contexte;
   le tableau peut grandir en cas de coupure (+1 en hauteur et +1 en largeur
   si nparm != 0) et en cas de partage (+1 en hauteur)
   (seulement si nparm != 0).
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

 mpz_set(rp->determinant, tp->determinant);
 pq = (Entier *) & (rp->row[virt+reel+dh]);
 for(i = off; i<virt + reel; i++)
     {ff = Flag(rp, i) = Flag(tp, i-off);
      mpz_set(Denom(rp, i), Denom(tp, i-off));
      if(ff & Unit) rp->row[i].objet.unit = tp->row[i-off].objet.unit;
      else {
	  rp->row[i].objet.val = pq;
	  pq +=(ncol + dw);
	  pp = tp->row[i-off].objet.val;
	  qq = rp->row[i].objet.val;
	  for(j = 0; j<ncol; j++) mpz_set(*qq++, *pp++);
	  }
      }
 return(rp);
}

int exam_coef(Tableau *tp, int nvar, int ncol, int bigparm)
{int i, j;
 int ff, fff;
 int x;
 Entier *p;
 for(i = 0; i<tp->height; i++)
     {ff = Flag(tp,i);
      if(ff == 0) break;
      if(ff == Unknown)
	  {if(bigparm >= 0){
	    x = mpz_sgn(Index(tp,i, bigparm));
	    if(x<0){
	      Flag(tp, i) = Minus;
	      return(i);
	    }
	    else if(x>0)
	      Flag(tp, i) = Plus;
	    continue;
	  }
	  ff = Zero;
	  p = &(tp->row[i].objet.val[nvar+1]);
	  for(j = nvar+1; j<ncol; j++){
	    x = mpz_sgn(*p); p++;
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
	  x = mpz_sgn(Index(tp, i, nvar));
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
     fprintf(stdout, "Too much parameters : %d\n", nparm);
     exit(1);
     }
 q = tab_hwm();
 for(i=0; i<=nparm; i++)
   mpz_init(discr[i]);
 
 for(i = 0; i<ni + nvar; i++)
     {ff = Flag(tp,i);
      if(ff & (Critic | Unknown))
	  {isCritic = True;
	   for(j = 0; j<nvar; j++) if(mpz_sgn(Index(tp, i, j)) > 0)
		 {isCritic = False;
		  break;
		 }
           compa_count++;
	   for(j = 0; j<nparm; j++)mpz_set(discr[j], Index(tp, i, j+nvar+1));
	   mpz_set(discr[nparm], Index(tp, i, nvar));
           mpz_sub_ui(discr[nparm], discr[nparm], (isCritic ? 0 : 1));
	   tPlus = expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tPlus, nparm+nc) = Unknown;
	   for(j = 0;j<=nparm; j++)mpz_set(Index(tPlus, nparm+nc, j),discr[j]);
	   mpz_set(Denom(tPlus, nparm+nc), UN);
	   p = sol_hwm();
	   traiter(tPlus, NULL, True, nparm, 0, nc+1, 0, -1);
	   cPlus = is_not_Nil(p);
	   if(verbose>0){
	     fprintf(dump, "\nThe positive case has been found ");
	     fprintf(dump, cPlus? "possible\n": "impossible\n");
	     fflush(dump);
	   }

	   sol_reset(p);
	   for(j = 0; j<nparm+1; j++) mpz_neg(discr[j], discr[j]);
	   mpz_sub_ui(discr[nparm], discr[nparm], (isCritic ? 1 : 2));
	   tMinus = expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tMinus, nparm+nc) = Unknown;
	   for(j=0;j<=nparm; j++)mpz_set(Index(tMinus, nparm+nc, j),discr[j]);
	   mpz_set(Denom(tMinus, nparm+nc), UN);
	   traiter(tMinus, NULL, True, nparm, 0, nc+1, 0, -1);
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

 for(i=0; i<=nparm; i++)
   mpz_clear(discr[i]);

 return;
}

Entier *valeur(Tableau *tp, int i, int j)
{
 if(Flag(tp, i) & Unit)
     return(tp->row[i].objet.unit == j ? &Denom(tp,i) : &ZERO);
 else return(&Index(tp, i, j));
}

void solution(Tableau *tp, int nvar, int nparm)
{int i, j;
 int ncol = nvar + nparm + 1;
 Entier *denom;

 sol_list(nvar);
 for(i = 0; i<nvar; i++)
     {sol_forme(nparm+1);
      for(j = nvar+1; j<ncol; j++)
	 sol_val(*valeur(tp, i, j), Denom(tp,i));
      sol_val(*valeur(tp, i, nvar), Denom(tp,i));
     }
}

int choisir_piv(Tableau *tp, int pivi, int nvar, int nligne)
{
 int j, k;
 Entier pivot, foo, x, y;
 int sgn_x;
 int pivj = -1;

 mpz_init(pivot); mpz_init(foo); mpz_init(x); mpz_init(y);

 for(j = 0; j<nvar; j++){
   mpz_set(foo, Index(tp, pivi, j));
   if(mpz_sgn(foo) <= 0) continue;
   if(pivj < 0){
     pivj = j;
     mpz_set(pivot, foo);
     continue;
   }
   for(k = 0; k<nligne; k++){
     mpz_mul(x, pivot, *valeur(tp, k, j)); 
     mpz_mul(y, *valeur(tp, k, pivj), foo);
     mpz_sub(x, x, y);
     cross_product++;
     sgn_x = mpz_sgn(x);
     if(sgn_x) break;
   }
   if(sgn_x < 0){
     pivj = j;
     mpz_set(pivot, foo);
   }
 }
 mpz_clear(pivot); mpz_clear(foo); mpz_clear(x); mpz_clear(y);
 return(pivj);
}

int pivoter(Tableau *tp, int pivi, int nvar,
	       int nparm, int ni, int iq)
{int pivj;
 int ncol = nvar + nparm + 1;
 int nligne = nvar + ni;
 int i, j, k;
 Entier x, y, d, gcd, u, dpiv;
 int ff, fff;
 Entier pivot, foo, z;
 Entier ppivot, dppiv;
 Entier new[MAXCOL], *p, *q;
 Entier lpiv;
 int sgn_x;


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

 mpz_init(x); mpz_init(y); mpz_init(d); 
 mpz_init(gcd); mpz_init(u); mpz_init(dpiv);
 mpz_init(lpiv); mpz_init(pivot); mpz_init(foo);
 mpz_init(z); mpz_init(ppivot); mpz_init(dppiv);

 for(i=0; i<ncol; i++)
   mpz_init(new[i]);

 mpz_set(pivot, Index(tp, pivi, pivj));
 mpz_set(dpiv, Denom(tp, pivi));
 mpz_gcd(d, pivot, dpiv);
 mpz_divexact(ppivot, pivot, d);
 mpz_divexact(dppiv, dpiv, d);
 if(verbose>0){
   fprintf(dump, "Pivot ");
   mpz_out_str(dump, 10, ppivot);
   putc('/', dump);
   mpz_out_str(dump, 10, dppiv);
   putc('\n', dump);
   fprintf(dump, "%d x %d\n", pivi, pivj);
 }
 mpz_fdiv_qr(x, y, tp->determinant, dppiv); 

 if(mpz_sgn(y) != 0){ 
   fprintf(stdout, "Computation error\n");
   if(verbose>0) fflush(dump);
   exit(1);
 }
 mpz_mul(tp->determinant, x, ppivot);

 if(verbose>0){
   fprintf(dump, "determinant ");
   mpz_out_str(dump, 10, tp->determinant);
   fprintf(dump, "\n");
 }

 
 for(j = 0; j<ncol; j++)
   if(j==pivj)
     mpz_set(new[j], dpiv);
   else 
     mpz_neg(new[j], Index(tp, pivi, j));
 for(k = 0; k<nligne; k++){
   if(Flag(tp,k) & Unit)continue;
   if(k == pivi)continue;
   /* foo = Index(tp, k, pivj) */
   mpz_set(foo, Index(tp, k, pivj));
   /* d = gcd(pivot, foo); */
   mpz_gcd(d, pivot, foo);
   /* lpiv = pivot/d; */
   mpz_divexact(lpiv, pivot, d);
   /* foo /= d; */
   mpz_divexact(foo, foo, d);
   /* d = Denom(tp,k); */
   mpz_set(d, Denom(tp,k));
   mpz_mul(gcd, lpiv, d);
   mpz_set(Denom(tp, k), gcd);
   p = tp->row[k].objet.val;
   q = tp->row[pivi].objet.val;
   for(j = 0; j<ncol; j++){
     if(j == pivj)
       mpz_mul(z, dpiv, foo);
     else {
       mpz_mul(z, *p, lpiv);
       mpz_mul(y, *q, foo);
       mpz_sub(z, z, y);
     }
     q++;
     cross_product++;
     mpz_set(*p, z);
     p++;
     if(mpz_cmp_ui(gcd, 1) != 0)
       mpz_gcd(gcd, gcd, z);
   }
   if(mpz_cmp_ui(gcd, 1) != 0){
     p = tp->row[k].objet.val;
     for(j = 0; j<ncol; j++){
       mpz_divexact(*p, *p, gcd);
       p++;
     }
   }
   mpz_divexact(Denom(tp,k), Denom(tp,k), gcd);
 }
 p = tp->row[pivi].objet.val;
 for(k = 0; k<nligne; k++)
   if((Flag(tp, k) & Unit) && tp->row[k].objet.unit == pivj) break;
 Flag(tp, k) = Plus;
 tp->row[k].objet.val = p;
 for(j = 0; j<ncol; j++)
   /* *p++ = new[j] */
   mpz_set(*p++, new[j]);
 mpz_set(Denom(tp, k), pivot);
 Flag(tp, pivi) = Unit | Zero;
 mpz_set(Denom(tp, pivi), UN);
 tp->row[pivi].objet.unit = pivj;

 for(k = 0; k<nligne; k++){
   ff = Flag(tp, k);
   if(ff & Unit) continue;
   sgn_x = mpz_sgn(Index(tp, k, pivj));
   if(sgn_x < 0) fff = Minus;
   else if(sgn_x == 0) fff = Zero;
   else fff = Plus;
   if(fff != Zero && fff != ff)
     if(ff == Zero) ff = (fff == Minus ? Unknown : fff);
     else ff = Unknown;
   Flag(tp, k) = ff;
 }

 if(verbose>0){
   fprintf(dump, "just pivoted\n");
   tab_display(tp, dump);
 }

 mpz_clear(x); mpz_clear(y); mpz_clear(d); mpz_clear(gcd);
 mpz_clear(u); mpz_clear(dpiv); mpz_clear(lpiv);
 mpz_clear(pivot); mpz_clear(foo); mpz_clear(z);
 mpz_clear(ppivot); mpz_clear(dppiv);

 for(i=0; i<ncol; i++)
   mpz_clear(new[i]);

 return(0);
}

/* dans cette version, "traiter" modifie ineq; par contre
   le contexte est immediatement recopie' */

void traiter(tp, ctxt, iq, nvar, nparm, ni, nc, bigparm)
Tableau *tp, *ctxt;
int iq, nvar, nparm, ni, nc, bigparm;

{
 int j;
 int pivi, nligne, ncol;
 struct high_water_mark x;
 Tableau *context;
 int dch, dcw;
 double s, t, d, smax;
 int i;
 int sgn_x;
 struct L temp;
 Entier discr[MAXPARM];
 for(i=0; i<=MAXPARM; i++)
   mpz_init(discr[i]);

 dcw = mpz_sizeinbase(tp->determinant, 2);
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
   d = mpz_get_d(Denom(tp,i));     
   for(j=0; j<nvar; j++){
     t = mpz_get_d(Index(tp,i,j))/d;
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
   fflush(dump);
 }
   
 for(;;) {
   nligne = nvar+ni; ncol = nvar+nparm+1;
   if(nligne > tp->height || ncol > tp->width) {
     fprintf(stdout, "Syserr : traiter : tableau too small\n");
     exit(1);
   }
   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;	       /* There is a negative row   */
   
   pivi = exam_coef(tp, nvar, ncol, bigparm);

   if(verbose>0){
     fprintf(dump, "coefs examined\n");
     tab_display(tp, dump);
     fflush(dump);
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
       dcw = mpz_sizeinbase(context->determinant,2);
       dch = 2 * dcw + 1;
       context = expanser(context, 0, nc, nparm+1, 0, dch, dcw);
     }
     if(nparm >= MAXPARM) {
       fprintf(stdout, "Too much parameters : %d\n", nparm);
       exit(2);
     }
     q = tab_hwm();
     if(verbose>0)
       fprintf(stdout,"profondeur %d %lx\n", profondeur, q.top);
     ntp = expanser(tp, nvar, ni, ncol, 0, 0, 0);
     fflush(stdout);
     sol_if();
     sol_forme(nparm+1);
     mpz_init_set_ui(com_dem, 0);
     for(j = 0; j<nparm; j++) {
       mpz_set(discr[j], Index(tp, pivi, j + nvar +1));
       mpz_gcd(com_dem, com_dem, discr[j]);
     }
     mpz_set(discr[nparm], Index(tp, pivi, nvar));
     mpz_gcd(com_dem, com_dem, discr[nparm]);
     for(j = 0; j<=nparm; j++) {
       mpz_divexact(discr[j], discr[j], com_dem);
       mpz_set(Index(context, nc, j), discr[j]);
       sol_val(discr[j], UN);
     }
     mpz_clear(com_dem);
     Flag(context, nc) = Unknown;
     mpz_set(Denom(context, nc), UN);
     Flag(ntp, pivi) = Plus;
     profondeur++;
     fflush(stdout);
     if(verbose > 0) fflush(dump);
     traiter(ntp, context, iq, nvar, nparm, ni, nc+1, bigparm);
     profondeur--;
     tab_reset(q);
     if(verbose>0)
       fprintf(stdout, "descente %d %lx\n", profondeur, tab_hwm().top);
     for(j = 0; j<nparm; j++)
       mpz_neg(Index(context, nc, j), Index(context, nc, j));
     mpz_add_ui(Index(context, nc, nparm), Index(context, nc, nparm), 1);
     mpz_neg(Index(context, nc, nparm), Index(context, nc, nparm));
     Flag(tp, pivi) = Minus;
     mpz_set(Denom(context, nc), UN);
     nc++;
     goto pirouette;
   }
/* Here, all rows are positive. Do we need an integral solution?      */
   if(!iq) {
     solution(tp, nvar, nparm);
     break;
   }
/* Yes we do! */
   pivi = integrer(&tp, &context, &nvar, &nparm, &ni, &nc);
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
     if(pivoter(tp, pivi, nvar, nparm, ni, iq) < 0) {
       sol_nil();
       break;
     }
 }
/* Danger : a premature return would induce memory leaks   */
 tab_reset(x);
 for(i=0; i<MAXPARM; i++)
   mpz_clear(discr[i]);
 return;
}







