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

#include <piplib/piplib.h>

#define TAB_CHUNK 4096*sizeof(Entier)

static char *tab_free, *tab_top;
static struct A *tab_base;

extern int allocation;
extern long int cross_product, limit;
static int chunk_count;

int dgetc(FILE *);
#if defined(LINEAR_VALUE_IS_MP)
int dscanf(FILE *, int    *);
#else
int dscanf(FILE *, Entier *);
#endif

extern FILE * dump;

void tab_init(void)
{
 tab_free = malloc(sizeof (struct A));
 if(tab_free == NULL)
     {fprintf(stderr, "Your computer doesn't have enough memory\n");
      exit(1);
     }
 allocation = 1;
 tab_top = tab_free + sizeof (struct A);
 tab_base = (struct A *)tab_free;
 tab_free += sizeof(struct A);
 tab_base->precedent = NULL;
 tab_base->bout = tab_top;
 tab_base->free = tab_free;
 chunk_count = 1;
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
      p = (char *)tab_base + sizeof(struct A);
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
      d = taille + sizeof(struct A);
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
      tab_free += sizeof(struct A);
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
 #if defined(LINEAR_VALUE_IS_MP)
 int x ;
 #else
 Entier x;
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
	if(dscanf(foo, &x) < 0)
		return NULL;
        else
        #if defined(LINEAR_VALUE_IS_MP)
	mpz_set_si(p->row[i].objet.val[j], x);
        #else
	p->row[i].objet.val[j] = x;
        #endif
        }
      } 
      while((c = dgetc(foo)) != EOF)if(c == ']')break;
     
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
 * a dire en fait le nombre d'inconnues). Si Max vaut 0, on va rechercher
 * le minimum lexicographique, sinon on recherche le maximum. La fonction
 * met alors en place le bignum s'il n'y est pas deja et prepare les
 * contraintes au calcul du maximum lexicographique.
 * 27 juillet 2001 : Premiere version, Ced.
 * 30 juillet 2001 : Nombreuses modifications. Le calcul du nombre total
 *                   d'inequations (Nineq) se fait a present a l'exterieur.
 *  3 octobre 2001 : Pas mal d'ameliorations.
 * 18 octobre 2003 : Mise en place de la possibilite de calculer le
 *                   maximum lexicographique (parties 'if (Max)').
 */
Tableau * tab_Matrix2Tableau(matrix, Nineq, Nv, n, Max, Bg)
PipMatrix * matrix ;
int Nineq, Nv, n, Max, Bg ;
{ Tableau * p ;
  unsigned i, j, end, current, new, nb_columns, decal=0, bignum_is_new ;
  Entier * entier, inequality, bignum ;
  
  #if defined(LINEAR_VALUE_IS_MP)
  mpz_init(inequality) ;
  mpz_init(bignum) ;
  #endif
  nb_columns = matrix->NbColumns - 1 ;
  /* S'il faut un BigNum et qu'il n'existe pas, on lui reserve sa place. */
  if ((Max) && (bignum_is_new = (Bg > (matrix->NbColumns - 2))))
  nb_columns ++ ;

  p = tab_alloc(Nineq,nb_columns,n) ;
    
  /* La variable decal sert a prendre en compte les lignes supplementaires
   * issues des egalites.
   */
  end = matrix->NbRows + n ;
  for (i=n;i<end;i++)
  { current = i + decal ;
    Flag(p,current) = Unknown ;
    #if defined(LINEAR_VALUE_IS_MP)
    mpz_set_ui(Denom(p,current),1) ;
    if (Max)
    mpz_set_ui(bignum,0) ;
    #else
    Denom(p,current) = UN ;
    if (Max)
    bignum = 0 ;
    #endif
    entier = *(matrix->p + i - n) ;
    /* Pour passer l'indicateur d'egalite/inegalite. */
    #if defined(LINEAR_VALUE_IS_MP)
    mpz_set(inequality,*entier) ;
    #else
    inequality = *entier ;
    #endif
    entier ++ ;
         
    /* Dans le format de la polylib, l'element constant est place en
     * dernier. Dans le format de Pip, il se trouve apres la premiere
     * serie de variables (inconnues ou parametres). On remet donc les
     * choses dans l'ordre de Pip. Ici pour p(x) >= 0.
     */
    if (Max)
    { for (j=0;j<Nv;j++)
      {
        #if defined(LINEAR_VALUE_IS_MP)
        mpz_add(bignum,bignum,*entier) ;
	mpz_neg(*(p->row[current].objet.val + j),*entier++) ;
	#else
        bignum += *entier ;
        *(p->row[current].objet.val + j) = -(*entier++) ;
        #endif
      }

      for (j=Nv+1;j<Bg;j++)
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_set(*(p->row[current].objet.val + j),*entier++) ;
      #else
      *(p->row[current].objet.val + j) = *entier++ ;
      #endif

      if (bignum_is_new)
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_set(*(p->row[current].objet.val + Bg),bignum) ;
      else
      { mpz_add(*(p->row[current].objet.val + Bg),bignum,*entier++) ;
        for (j=Bg+1;j<nb_columns;j++)
        mpz_set(*(p->row[current].objet.val + j),*entier++) ;
      }
      #else
      *(p->row[current].objet.val + Bg) = bignum ;
      else
      { *(p->row[current].objet.val + Bg) = bignum + *entier++ ;
        for (j=Bg+1;j<nb_columns;j++)
        *(p->row[current].objet.val + j) = *entier++ ;
      }
      #endif
 
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_set(*(p->row[current].objet.val + Nv),*entier) ;
      #else
      *(p->row[current].objet.val + Nv) = *entier ;
      #endif
    }
    else
    { for (j=0;j<Nv;j++)
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_set(*(p->row[current].objet.val + j),*entier++) ;
      #else
      *(p->row[current].objet.val + j) = *entier++ ;
      #endif
      for (j=Nv+1;j<nb_columns;j++)
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_set(*(p->row[current].objet.val + j),*entier++) ;
      mpz_set(*(p->row[current].objet.val + Nv),*entier) ;
      #else
      *(p->row[current].objet.val + j) = *entier++ ;
      *(p->row[current].objet.val + Nv) = *entier ;
      #endif
    }
    
    /* Et ici lors de l'ajout de -p(x) >= 0 quand on traite une egalite. */
    #if defined(LINEAR_VALUE_IS_MP)
    if (mpz_sgn(inequality) == 0)
    #else
    if (!inequality)
    #endif
    { decal ++ ;
      new = current + 1 ;
      Flag(p,new)= Unknown ;
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_set(Denom(p,new),UN) ;
      #else
      Denom(p,new) = UN ;
      #endif
      
      for (j=0;j<nb_columns;j++)
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_neg(*(p->row[new].objet.val + j),*(p->row[current].objet.val + j)) ;
      #else
      *(p->row[new].objet.val + j) = -(*(p->row[current].objet.val + j)) ;
      #endif
    }
  }
  #if defined(LINEAR_VALUE_IS_MP)
  mpz_clear(inequality);
  mpz_clear(bignum);
  #endif

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
