/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                   tab.h                                    *
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
 * Written by Paul Feautrier and Cedric Bastoul                               *
 *                                                                            *
 *****************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "pip.h"

#define TAB_CHUNK 4096*sizeof(PIPLIB_NAME(piplib_int_t))

static char *PIPLIB_NAME(tab_free), *PIPLIB_NAME(tab_top);
static struct PIPLIB_NAME(A) *PIPLIB_NAME(tab_base);

/*extern long int cross_product;*/
static int PIPLIB_NAME(chunk_count);

int PIPLIB_NAME(dgetc)(FILE *);

extern FILE * PIPLIB_NAME(dump);

#define sizeof_struct_A ((sizeof(struct PIPLIB_NAME(A)) % sizeof(PIPLIB_NAME(piplib_int_t))) ?		    \
			 (sizeof(struct PIPLIB_NAME(A)) + sizeof(PIPLIB_NAME(piplib_int_t))		    \
				- (sizeof(struct PIPLIB_NAME(A)) % sizeof(PIPLIB_NAME(piplib_int_t)))) :    \
			  sizeof(struct PIPLIB_NAME(A)))

void PIPLIB_NAME(tab_init)(void)
{
 PIPLIB_NAME(tab_free) = malloc(sizeof_struct_A);
 if(PIPLIB_NAME(tab_free) == NULL)
     {fprintf(stderr, "Your computer doesn't have enough memory\n");
      exit(1);
     }
 PIPLIB_NAME(tab_top) = PIPLIB_NAME(tab_free) + sizeof_struct_A;
 PIPLIB_NAME(tab_base) = (struct PIPLIB_NAME(A) *)PIPLIB_NAME(tab_free);
 PIPLIB_NAME(tab_free) += sizeof_struct_A;
 PIPLIB_NAME(tab_base)->precedent = NULL;
 PIPLIB_NAME(tab_base)->bout = PIPLIB_NAME(tab_top);
 PIPLIB_NAME(tab_base)->free = PIPLIB_NAME(tab_free);
 PIPLIB_NAME(chunk_count) = 1;
}
 
 
void PIPLIB_NAME(tab_close)(void)
{
  if (PIPLIB_NAME(tab_base)) free(PIPLIB_NAME(tab_base));
}


struct PIPLIB_NAME(high_water_mark) PIPLIB_NAME(tab_hwm)(void)
{struct PIPLIB_NAME(high_water_mark) p;
 p.chunk = PIPLIB_NAME(chunk_count);
 p.top = PIPLIB_NAME(tab_free);
 return p;
}


#if defined(PIPLIB_ONE_DETERMINANT)
/* the clear_tab routine clears the GMP objects which may be referenced
   in the given PIPLIB_NAME(Tableau).
*/
void PIPLIB_NAME(tab_clear)(PIPLIB_NAME(Tableau) *tp)
{
  int i, j;
  /* clear the determinant */
  piplib_int_clear(tp->determinant);

  for(i=0; i<tp->height; i++){
    /* clear the denominator */
    piplib_int_clear(Denom(tp, i));
    if((Flag(tp, i) & Unit) == 0)
      for(j=0; j<tp->width; j++)
        piplib_int_clear(Index(tp,i,j));
  }
}
#endif

void PIPLIB_NAME(tab_reset)(struct PIPLIB_NAME(high_water_mark) by_the_mark) {
 struct PIPLIB_NAME(A) *g;
 #if defined(PIPLIB_ONE_DETERMINANT)
  char *p;
 #endif
 while(PIPLIB_NAME(chunk_count) > by_the_mark.chunk)
     {
      g = PIPLIB_NAME(tab_base)->precedent;
      
      #if defined(PIPLIB_ONE_DETERMINANT)
      /* Before actually freeing the memory, one has to clear the
       * included PIPLIB_NAME(Tableau)x. If this is not done, the GMP objects
       * referenced in the PIPLIB_NAME(Tableau)x will be orphaned.
       */

      /* Enumerate the included tableaux. */
      p = (char *)PIPLIB_NAME(tab_base) + sizeof_struct_A;
      while(p < PIPLIB_NAME(tab_base)->free){
        PIPLIB_NAME(Tableau) *pt;
        pt = (PIPLIB_NAME(Tableau) *) p;
	PIPLIB_NAME(tab_clear)(pt);
        p += pt->taille;
      } 
      #endif
      
      free(PIPLIB_NAME(tab_base));
      PIPLIB_NAME(tab_base) = g;
      PIPLIB_NAME(tab_top) = PIPLIB_NAME(tab_base)->bout;
      PIPLIB_NAME(chunk_count)--;
     }
 if(PIPLIB_NAME(chunk_count) > 0) {
     #if defined(PIPLIB_ONE_DETERMINANT)
     /* Do not forget to clear the tables in the current chunk above the
        high water mark */
     p = (char *)by_the_mark.top;
     while(p < PIPLIB_NAME(tab_base)->free) {
        PIPLIB_NAME(Tableau) *pt;
        pt = (PIPLIB_NAME(Tableau) *) p;
        PIPLIB_NAME(tab_clear)(pt);
        p += pt->taille;
        } 
     #endif   
     PIPLIB_NAME(tab_free) = by_the_mark.top;
     PIPLIB_NAME(tab_base)->free = PIPLIB_NAME(tab_free);
     }
 else {
     fprintf(stderr,
             "Syserr: PIPLIB_NAME(tab_reset) : error in memory allocation\n");
     exit(1);
     }
}

PIPLIB_NAME(Tableau) * PIPLIB_NAME(tab_alloc)(int h, int w, int n)

/* h : le nombre de ligne reelles;
   n : le nombre de lignes virtuelles
*/
{
 char *p; PIPLIB_NAME(Tableau) *tp;
 PIPLIB_NAME(piplib_int_t) *q;
 unsigned long taille;
 int i, j;
 taille = sizeof(PIPLIB_NAME(Tableau))
          + (h+n-1) * sizeof (struct PIPLIB_NAME(L))
	  + h * w * sizeof (PIPLIB_NAME(piplib_int_t));
 if(PIPLIB_NAME(tab_free) + taille >= PIPLIB_NAME(tab_top))
     {struct PIPLIB_NAME(A) * g;
      unsigned long d;
      d = taille + sizeof_struct_A;
      if(d < TAB_CHUNK) d = TAB_CHUNK;
      PIPLIB_NAME(tab_free) = malloc(d);
      if(PIPLIB_NAME(tab_free) == NULL)
	  {printf("Memory overflow\n");
	   exit(23);
	  }
      PIPLIB_NAME(chunk_count)++;
      g = (struct PIPLIB_NAME(A) *)PIPLIB_NAME(tab_free);
      g->precedent = PIPLIB_NAME(tab_base);
      PIPLIB_NAME(tab_top) = PIPLIB_NAME(tab_free) + d;
      PIPLIB_NAME(tab_free) += sizeof_struct_A;
      PIPLIB_NAME(tab_base) = g;
      g->bout = PIPLIB_NAME(tab_top);
     }
 p = PIPLIB_NAME(tab_free);
 PIPLIB_NAME(tab_free) += taille;
 PIPLIB_NAME(tab_base)->free = PIPLIB_NAME(tab_free);
 tp = (PIPLIB_NAME(Tableau) *)p;
 q = (PIPLIB_NAME(piplib_int_t) *)(p +  sizeof(PIPLIB_NAME(Tableau))
     + (h+n-1) * sizeof (struct PIPLIB_NAME(L)));
 #if defined(PIPLIB_ONE_DETERMINANT)
 piplib_int_init_set_si(tp->determinant,1);
 #else
 tp->determinant[0] = (PIPLIB_NAME(piplib_int_t)) 1;
 tp->l_determinant = 1;
 #endif
 for(i = 0; i<n; i++){
   tp->row[i].flags = Unit;
   tp->row[i].objet.unit = i;
   piplib_int_init_set_si(Denom(tp, i), 1);
 }
 for(i = n; i < (h+n); i++){
   tp->row[i].flags = 0;
   tp->row[i].objet.val = q;
   tp->row[i].size = 0;
   for(j = 0; j < w; j++)
   piplib_int_init_set_si(*q++, 0); /* loop body. */
   piplib_int_init_set_si(Denom(tp, i), 0);
 }
 tp->height = h + n; tp->width = w;
 #if defined(PIPLIB_ONE_DETERMINANT)
 tp->taille = taille ;
 #endif
 
 return(tp);
}

PIPLIB_NAME(Tableau) * PIPLIB_NAME(tab_get)(foo, h, w, n)
FILE * foo;
int h, w, n;
{
 PIPLIB_NAME(Tableau) *p;
 int i, j, c;
 PIPLIB_NAME(piplib_int_t) x;
 piplib_int_init(x);
 
 p = PIPLIB_NAME(tab_alloc)(h, w, n);
 while((c = PIPLIB_NAME(dgetc)(foo)) != EOF)
      if(c == '(')break;
 for(i = n; i<h+n; i++)
     {p->row[i].flags = Unknown;
      piplib_int_set_si(Denom(p, i), 1);
      while((c = PIPLIB_NAME(dgetc)(foo)) != EOF)if(c == '[')break;
      for(j = 0; j<w; j++){
        if(PIPLIB_NAME(dscanf)(foo, &x) < 0) return NULL;
        else piplib_int_assign(p->row[i].objet.val[j], x);
        }
      } 
      while((c = PIPLIB_NAME(dgetc)(foo)) != EOF)if(c == ']')break;
 
 piplib_int_clear(x);
     
 return(p);
}


/* Fonction tab_Matrix2Tableau :
 * Cette fonction effectue la conversion du format de matrice de la polylib
 * vers le format de traitement de Pip. matrix est la matrice a convertir.
 * Nineq est le nombre d'inequations necessaires (dans le format de la
 * polylib, le premier element d'une ligne indique si l'equation decrite
 * est une inequation ou une egalite. Pip ne gere que les inequations. On
 * compte donc le nombre d'inequations total pour reserver la place
 * necessaire, et on scinde toute egalite p(x)=0 en p(x)>=0 et -p(x)>=0).
 * Nv est le nombre de variables dans la premiere serie de variables (c'est
 * a dire que si les premiers coefficients dans les lignes de la matrice
 * sont ceux des inconnues, Nv est le nombre d'inconnues, resp. parametres).
 * n est le nombre de lignes 'virtuelles' contenues dans la matrice (c'est
 * a dire en fait le nombre d'inconnues). Si Shift vaut 0, on va rechercher
 * le minimum lexicographique non-negatif, sinon on recherche le maximum 
 * (Shift = 1) ou bien le minimum tout court (Shift = -1). La fonction
 * met alors en place le bignum s'il n'y est pas deja et prepare les
 * contraintes au calcul du maximum lexicographique.
 *
 * This function is called both for both the context (only parameters)
 * and the actual domain (variables + parameters).
 * Let Np be the number of parameters and Nn the number of variables.
 *
 * For the context, the columns in matrix are
 *		1 Np 1
 * while the result has
 *		Np Bg Urs_parms 1
 * Nv = Np + Bg; n = -1
 * 
 * For the domain, matrix has
 *		1 Nn Np 1
 * while the result has
 *		Nn 1 Np Bg Urs_parms
 * Nv = Nn; n >= 0
 *
 * 27 juillet 2001 : Premiere version, Ced.
 * 30 juillet 2001 : Nombreuses modifications. Le calcul du nombre total
 *                   d'inequations (Nineq) se fait a present a l'exterieur.
 *  3 octobre 2001 : Pas mal d'ameliorations.
 * 18 octobre 2003 : Mise en place de la possibilite de calculer le
 *                   maximum lexicographique (parties 'if (Max)').
 */
PIPLIB_NAME(Tableau) * PIPLIB_NAME(tab_Matrix2Tableau)(matrix, Nineq, Nv, n,
                                                       Shift, Bg, Urs_parms)
PIPLIB_NAME(PipMatrix) * matrix ;
int Nineq, Nv, n, Shift, Bg, Urs_parms;
{ PIPLIB_NAME(Tableau) * p ;
  unsigned i, k, current, new, nb_columns, decal=0, bignum_is_new;
  int j;
  unsigned cst;
  int inequality, ctx;
  PIPLIB_NAME(piplib_int_t) bignum;

  /* Are we dealing with the context? */
  ctx = n == -1;
  if (ctx)
    n = 0;
  piplib_int_init(bignum);
  nb_columns = matrix->NbColumns - 1 ;
  
  /* S'il faut un BigNum et qu'il n'existe pas, on lui reserve sa place. */
  bignum_is_new = Shift
  && (Bg+ctx > 0)
  && ((unsigned int)(Bg+ctx) > (matrix->NbColumns - 2));
  
  if (bignum_is_new)
    nb_columns++;
  if (ctx) {
    Shift = 0;
    cst = Nv + Urs_parms;
  } else
    cst = Nv;

  p = PIPLIB_NAME(tab_alloc)(Nineq,nb_columns+Urs_parms,n) ;
    
  /* La variable decal sert a prendre en compte les lignes supplementaires
   * issues des egalites.
   */
  for (i = 0; i < matrix->NbRows; i++) {
    current = i + n + decal;
    Flag(p,current) = Unknown ;
    piplib_int_set_si(Denom(p,current), 1);
    if (Shift)
      piplib_int_set_si(bignum, 0);
    /* Pour passer l'indicateur d'egalite/inegalite. */
    inequality = (piplib_int_zero(matrix->p[i][0]) == 0);
         
    /* Dans le format de la polylib, l'element constant est place en
     * dernier. Dans le format de Pip, il se trouve apres la premiere
     * serie de variables (inconnues ou parametres). On remet donc les
     * choses dans l'ordre de Pip. Ici pour p(x) >= 0.
     */
    for (j=0;j<Nv;j++) {
      if (bignum_is_new && j == Bg)
	continue;
      if (Shift)
	piplib_int_add(bignum, bignum, matrix->p[i][1+j]);
      if (Shift > 0)
	piplib_int_oppose(p->row[current].objet.val[j], matrix->p[i][1+j]);
      else
	piplib_int_assign(p->row[current].objet.val[j], matrix->p[i][1+j]);
    }
    for (k=j=Nv+1;(unsigned int)j<nb_columns;j++) {
	if (bignum_is_new && j == Bg)
	  continue;
	piplib_int_assign(p->row[current].objet.val[j], matrix->p[i][k]);
	k++;
    }
    for (j=0; j < Urs_parms; ++j) {
	int pos_n = nb_columns - ctx + j;
	int pos = pos_n - Urs_parms;
	if (pos <= Bg)
	    --pos;
	piplib_int_oppose(p->row[current].objet.val[pos_n],
		     p->row[current].objet.val[pos]);
    }
    piplib_int_assign(p->row[current].objet.val[cst], 
		 matrix->p[i][matrix->NbColumns-1]);
    if (Shift) {
      if (Shift < 0)
	piplib_int_oppose(bignum, bignum);

      if (bignum_is_new)
	piplib_int_assign(p->row[current].objet.val[Bg], bignum);
      else
	piplib_int_add(p->row[current].objet.val[Bg], 
		    p->row[current].objet.val[Bg], bignum);
    }
    
    /* Et ici lors de l'ajout de -p(x) >= 0 quand on traite une egalite. */
    if (!inequality) {
      decal ++ ;
      new = current + 1 ;
      Flag(p,new)= Unknown ;
      piplib_int_set_si(Denom(p,new), 1);
      
      for (j=0;(unsigned int)j<nb_columns+Urs_parms;j++)
	piplib_int_oppose(p->row[new].objet.val[j], p->row[current].objet.val[j]);
    }
  }
  piplib_int_clear(bignum);

  return(p);
}


int PIPLIB_NAME(tab_simplify)(PIPLIB_NAME(Tableau) *tp, int cst)
{
    int i, j;
    PIPLIB_NAME(piplib_int_t) gcd;

    piplib_int_init(gcd);
    for (i = 0; i < tp->height; ++i) {
	if (Flag(tp, i) & Unit)
	    continue;
	piplib_int_set_si(gcd, 0);
	for (j = 0; j < tp->width; ++j) {
	    if (j == cst)
		continue;
	    piplib_int_gcd(gcd, gcd, Index(tp, i, j));
	    if (piplib_int_one(gcd))
		break;
	}
	if (piplib_int_zero(gcd))
	    continue;
	if (piplib_int_one(gcd))
	    continue;
	for (j = 0; j < tp->width; ++j) {
	    if (j == cst)
		piplib_int_floor_div_q(Index(tp, i, j), Index(tp, i, j), gcd);
	    else
		piplib_int_div_exact(Index(tp, i, j), Index(tp, i, j), gcd);
	}
    }
    piplib_int_clear(gcd);

    return 0;
}


void PIPLIB_NAME(tab_display)(p, foo)
FILE *foo;
PIPLIB_NAME(Tableau) *p;
{
 char const * const Attr[] = {"Unit", "+", "-", "0", "*", "?"};

 int i, j, ff, fff, n;
 PIPLIB_NAME(piplib_int_t) d;
 piplib_int_init(d);

 fprintf(foo, "cross_product (%ld) /[%d * %d]\n", 0L/*cross_product*/,
              p->height, p->width);
 for(i = 0; i<p->height; i++){
   fff = ff = p->row[i].flags;
   /* if(fff ==0) continue; */
   piplib_int_assign(d, Denom(p, i));
   n = 0;
   while(fff){
     if(fff & 1) fprintf(foo, "%s ",Attr[n]);
     n++; fff >>= 1;
   }
   fprintf(foo, "%f #[", p->row[i].size);
   if(ff & Unit)
     for(j = 0; j<p->width; j++)
       fprintf(foo, " /%d/",(j == p->row[i].objet.unit)? 1: 0);
   else
     for(j = 0; j<p->width; j++){
       piplib_int_print(foo, Index(p, i, j));
       putc(' ', foo);
     }
   fprintf(foo, "]/");
   piplib_int_print(foo, d);
   putc('\n', foo);
 }
 piplib_int_clear(d);
}
