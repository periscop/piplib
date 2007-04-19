/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                   tab.h                                    *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.                                                                   *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Written by Paul Feautrier and Cedric Bastoul                               *
 *                                                                            *
 *****************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "pip.h"

#define TAB_CHUNK 4096*sizeof(Entier)

static char *tab_free, *tab_top;
static struct A *tab_base;

extern int allocation;
extern long int cross_product, limit;
static int chunk_count;

int dgetc(FILE *);
#if defined(LINEAR_VALUE_IS_MP)
int dscanf(FILE *, Entier);
#else
int dscanf(FILE *, Entier *);
#endif

extern FILE * dump;

#define sizeof_struct_A ((sizeof(struct A) % sizeof(Entier)) ?		    \
			 (sizeof(struct A) + sizeof(Entier)		    \
				- (sizeof(struct A) % sizeof(Entier))) :    \
			  sizeof(struct A))

void tab_init(void)
{
 tab_free = malloc(sizeof_struct_A);
 if(tab_free == NULL)
     {fprintf(stderr, "Your computer doesn't have enough memory\n");
      exit(1);
     }
 allocation = 1;
 tab_top = tab_free + sizeof_struct_A;
 tab_base = (struct A *)tab_free;
 tab_free += sizeof_struct_A;
 tab_base->precedent = NULL;
 tab_base->bout = tab_top;
 tab_base->free = tab_free;
 chunk_count = 1;
}
 
 
void tab_close(void)
{
  if (tab_base) free(tab_base);
}


struct high_water_mark tab_hwm(void)
{struct high_water_mark p;
 p.chunk = chunk_count;
 p.top = tab_free;
 return p;
}


#if defined(LINEAR_VALUE_IS_MP)
/* the clear_tab routine clears the GMP objects which may be referenced
   in the given Tableau.
*/
void tab_clear(Tableau *tp)
{
  int i, j;
  /* clear the determinant */
  mpz_clear(tp->determinant);

  for(i=0; i<tp->height; i++){
    /* clear the denominator */
    mpz_clear(Denom(tp, i));
    if((Flag(tp, i) & Unit) == 0)
      for(j=0; j<tp->width; j++)
        mpz_clear(Index(tp,i,j));
  }
}
#endif

void tab_reset(struct high_water_mark by_the_mark)

{struct A *g;
 char *p;
 while(chunk_count > by_the_mark.chunk)
     {
      g = tab_base->precedent;
      
      #if defined(LINEAR_VALUE_IS_MP)
      /* Before actually freeing the memory, one has to clear the
       * included Tableaux. If this is not done, the GMP objects
       * referenced in the Tableaux will be orphaned.
       */

      /* Enumerate the included tableaux. */
      p = (char *)tab_base + sizeof_struct_A;
      while(p < tab_base->free){
        Tableau *pt;
        pt = (Tableau *) p;
	tab_clear(pt);
        p += pt->taille;
      } 
      #endif
      
      free(tab_base);
      tab_base = g;
      tab_top = tab_base->bout;
      chunk_count--;
     }
 if(chunk_count > 0) {
     #if defined(LINEAR_VALUE_IS_MP)
     /* Do not forget to clear the tables in the current chunk above the
        high water mark */
     p = (char *)by_the_mark.top;
     while(p < tab_base->free) {
        Tableau *pt;
        pt = (Tableau *) p;
        tab_clear(pt);
        p += pt->taille;
        } 
     #endif   
     tab_free = by_the_mark.top;
     tab_base->free = tab_free;
     }
 else {
     fprintf(stderr, "Syserr: tab_reset : error in memory allocation\n");
     exit(1);
     }
}

Tableau * tab_alloc(int h, int w, int n)

/* h : le nombre de ligne reelles;
   n : le nombre de lignes virtuelles
*/
{
 char *p; Tableau *tp;
 Entier *q;
 unsigned long taille;
 int i, j;
 taille = sizeof(struct T) + (h+n-1) * sizeof (struct L)
	  + h * w * sizeof (Entier);
 if(tab_free + taille >= tab_top)
     {struct A * g;
      unsigned long d;
      d = taille + sizeof_struct_A;
      if(d < TAB_CHUNK) d = TAB_CHUNK;
      tab_free = malloc(d);
      if(tab_free == NULL)
	  {printf("Memory overflow\n");
	   exit(23);
	  }
      chunk_count++;
      g = (struct A *)tab_free;
      g->precedent = tab_base;
      tab_top = tab_free + d;
      tab_free += sizeof_struct_A;
      tab_base = g;
      g->bout = tab_top;
     }
 p = tab_free;
 tab_free += taille;
 tab_base->free = tab_free;
 tp = (Tableau *)p;
 q = (Entier *)(p +  sizeof(struct T) + (h+n-1) * sizeof (struct L));
 #if defined(LINEAR_VALUE_IS_MP)
 mpz_init_set_ui(tp->determinant,1);
 #else
 tp->determinant[0] = (Entier) 1;
 tp->l_determinant = 1;
 #endif
 for(i = 0; i<n ; i++){
   tp->row[i].flags = Unit;
   tp->row[i].objet.unit = i;
   #if defined(LINEAR_VALUE_IS_MP)
   mpz_init_set_ui(Denom(tp, i), 1);
   #else
   Denom(tp, i) = UN ;
   #endif
 }
 for(i = n; i < (h+n); i++){
   tp->row[i].flags = 0;
   tp->row[i].objet.val = q;
   for(j = 0; j < w; j++)
   #if defined(LINEAR_VALUE_IS_MP)
   mpz_init_set_ui(*q++, 0); /* loop body. */
   mpz_init_set_ui(Denom(tp, i), 0);
   #else
   *q++ = 0;                 /* loop body. */
   Denom(tp, i) = ZERO ;
   #endif
 }
 tp->height = h + n; tp->width = w;
 #if defined(LINEAR_VALUE_IS_MP)
 tp->taille = taille ;
 #endif
 
 return(tp);
}

Tableau * tab_get(foo, h, w, n)
FILE * foo;
int h, w, n;
{
 Tableau *p;
 int i, j, c;
 Entier x;
 #if defined(LINEAR_VALUE_IS_MP)
 mpz_init(x);
 #endif
 
 p = tab_alloc(h, w, n);
 while((c = dgetc(foo)) != EOF)
      if(c == '(')break;
 for(i = n; i<h+n; i++)
     {p->row[i].flags = Unknown;
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_set_ui(Denom(p, i), 1);
      #else
      Denom(p, i) = UN;
      #endif
      while((c = dgetc(foo)) != EOF)if(c == '[')break;
      for(j = 0; j<w; j++){
        #if defined(LINEAR_VALUE_IS_MP)
	if(dscanf(foo, x) < 0)
          return NULL;
        else
	  mpz_set(p->row[i].objet.val[j], x);
        #else
	if(dscanf(foo, &x) < 0)
          return NULL;
        else
	  p->row[i].objet.val[j] = x;
        #endif
        }
      } 
      while((c = dgetc(foo)) != EOF)if(c == ']')break;
 
 #if defined(LINEAR_VALUE_IS_MP)
 mpz_clear(x);
 #endif
     
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
 * 27 juillet 2001 : Premiere version, Ced.
 * 30 juillet 2001 : Nombreuses modifications. Le calcul du nombre total
 *                   d'inequations (Nineq) se fait a present a l'exterieur.
 *  3 octobre 2001 : Pas mal d'ameliorations.
 * 18 octobre 2003 : Mise en place de la possibilite de calculer le
 *                   maximum lexicographique (parties 'if (Max)').
 */
Tableau * tab_Matrix2Tableau(matrix, Nineq, Nv, n, Shift, Bg, Urs_parms)
PipMatrix * matrix ;
int Nineq, Nv, n, Shift, Bg, Urs_parms;
{ Tableau * p ;
  unsigned i, j, k, current, new, nb_columns, decal=0, bignum_is_new ;
  int inequality;
  Entier bignum;
  
  entier_init(bignum);
  nb_columns = matrix->NbColumns - 1 ;
  /* S'il faut un BigNum et qu'il n'existe pas, on lui reserve sa place. */
  bignum_is_new = Shift && (Bg > (matrix->NbColumns - 2));
  if (bignum_is_new)
    nb_columns++;
  /* Ce sont juste des parametres. */
  if (Bg <= Nv)
    Shift = 0;

  p = tab_alloc(Nineq,nb_columns+Urs_parms,n) ;
    
  /* La variable decal sert a prendre en compte les lignes supplementaires
   * issues des egalites.
   */
  for (i = 0; i < matrix->NbRows; i++) {
    current = i + n + decal;
    Flag(p,current) = Unknown ;
    entier_set_si(Denom(p,current), 1);
    if (Shift)
      entier_set_si(bignum, 0);
    /* Pour passer l'indicateur d'egalite/inegalite. */
    inequality = entier_notzero_p(matrix->p[i][0]);
         
    /* Dans le format de la polylib, l'element constant est place en
     * dernier. Dans le format de Pip, il se trouve apres la premiere
     * serie de variables (inconnues ou parametres). On remet donc les
     * choses dans l'ordre de Pip. Ici pour p(x) >= 0.
     */
    for (j=0;j<Nv;j++) {
      if (bignum_is_new && 1+j == Bg)
	continue;
      if (Shift)
	entier_addto(bignum, bignum, matrix->p[i][1+j]);
      if (Shift > 0)
	entier_oppose(p->row[current].objet.val[j], matrix->p[i][1+j]);
      else
	entier_assign(p->row[current].objet.val[j], matrix->p[i][1+j]);
    }
    for (k=j=Nv+1;j<nb_columns;j++) {
	if (bignum_is_new && j == Bg)
	  continue;
	entier_assign(p->row[current].objet.val[j], matrix->p[i][k++]);
    }
    for (j=0; j < Urs_parms; ++j) {
	int pos = nb_columns - Urs_parms + j;
	if (pos <= Nv)
	    --pos;
	if (pos <= Bg)
	    --pos;
	entier_oppose(p->row[current].objet.val[nb_columns+j],
		     p->row[current].objet.val[pos]);
    }
    entier_assign(p->row[current].objet.val[Nv], 
		 matrix->p[i][matrix->NbColumns-1]);
    if (Shift) {
      if (Shift < 0)
	entier_oppose(bignum, bignum);

      if (bignum_is_new)
	entier_assign(p->row[current].objet.val[Bg], bignum);
      else
	entier_addto(p->row[current].objet.val[Bg], 
		    p->row[current].objet.val[Bg], bignum);
    }
    
    /* Et ici lors de l'ajout de -p(x) >= 0 quand on traite une egalite. */
    if (!inequality) {
      decal ++ ;
      new = current + 1 ;
      Flag(p,new)= Unknown ;
      entier_set_si(Denom(p,new), 1);
      
      for (j=0;j<nb_columns+Urs_parms;j++)
	entier_oppose(p->row[new].objet.val[j], p->row[current].objet.val[j]);
    }
  }
  entier_clear(bignum);

  return(p);
}


char *Attr[] = {"Unit", "+", "-", "0", "*", "?"};

void tab_display(p, foo)
FILE *foo;
Tableau *p;
{

 int i, j, ff, fff, n;
 Entier x, d;
 #if defined(LINEAR_VALUE_IS_MP)
 mpz_init(d);
 #endif

 fprintf(foo, "%ld/[%d * %d]\n", cross_product, p->height, p->width);
 for(i = 0; i<p->height; i++){
   fff = ff = p->row[i].flags;
   /* if(fff ==0) continue; */
   #if defined(LINEAR_VALUE_IS_MP)
   mpz_set(d, Denom(p, i));
   #else
   d = Denom(p, i);
   #endif
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
       #if defined(LINEAR_VALUE_IS_MP)
       mpz_out_str(foo, 10, Index(p, i, j));
       putc(' ', foo);
       #else
       x = Index(p,i,j);
       fprintf(foo, FORMAT, x);
       fprintf(foo, " ");
       #endif
     }
   fprintf(foo, "]/");
   #if defined(LINEAR_VALUE_IS_MP)
   mpz_out_str(foo, 10, d);
   #else
   fprintf(foo, "%d", (int)d);
   #endif
   putc('\n', foo);
 }
 #if defined(LINEAR_VALUE_IS_MP)
 mpz_clear(d);
 #endif
}
