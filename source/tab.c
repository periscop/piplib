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
 * Written by Paul Feautrier and Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/


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
int dscanf(FILE *, char *, Entier *);

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
 chunk_count = 1;
}

struct high_water_mark tab_hwm(void)
{struct high_water_mark p;
 p.chunk = chunk_count;
 p.top = tab_free;
 return p;
}

void tab_reset(struct high_water_mark p)

{struct A *g;
 while(chunk_count > p.chunk)
     {
      g = tab_base->precedent;
      free(tab_base);
      tab_base = g;
      tab_top = tab_base->bout;
      chunk_count--;
     }
 if(chunk_count > 0) tab_free = p.top;
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
 tp = (Tableau *)p;
 q = (Entier *)(p +  sizeof(struct T) + (h+n-1) * sizeof (struct L));
 tp->determinant[0] = (Entier) 1;
 tp->l_determinant = 1;
 for(i = 0; i<n ; i++){
   tp->row[i].flags = Unit;
   tp->row[i].objet.unit = i;
   Denom(tp, i) = UN ;
 }
 for(i = n; i < (h+n); i++){
   tp->row[i].flags = 0;
   tp->row[i].objet.val = q;
   for(j = 0; j < w; j++) *q++ = 0;
   Denom(tp, i) = ZERO ;
 }
 tp->height = h + n; tp->width = w;
 return(tp);
}

Tableau * tab_get(foo, h, w, n)
FILE * foo;
int h, w, n;
{
 Tableau *p;
 int i, j, c;
 Entier x;
 p = tab_alloc(h, w, n);
 while((c = dgetc(foo)) != EOF)
      if(c == '(')break;
 for(i = n; i<h+n; i++)
     {p->row[i].flags = Unknown;
      Denom(p, i) = UN;
      while((c = dgetc(foo)) != EOF)if(c == '[')break;
      for(j = 0; j<w; j++){
	if(dscanf(foo, FORMAT, &x) < 0)
		return NULL;
        else p->row[i].objet.val[j] = x;
        }
      while((c = dgetc(foo)) != EOF)if(c == ']')break;
     }
 while((c = dgetc(foo)) != EOF)if(c == ')')break;
 return((Tableau *) p);
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
 * a dire en fait le nombre d'inconnues).
 * 27 juillet 2001 : Premiere version, Ced.
 * 30 juillet 2001 : Nombreuses modifications. Le calcul du nombre total
 *                   d'inequations (Nineq) se fait a present a l'exterieur.
 *  3 octobre 2001 : Pas mal d'ameliorations. 
 */
Tableau * tab_Matrix2Tableau(PipMatrix * matrix, int Nineq, int Nv, int n)
{ Tableau * p ;
  unsigned i, j, end, current, new, nb_columns, decal=0 ;
  Entier * entier, inequality ;
  
  nb_columns = matrix->NbColumns - 1 ;
  p = tab_alloc(Nineq,nb_columns,n) ;
    
  /* La variable decal sert a prendre en compte les lignes supplementaires
   * issues des egalites.
   */
  end = matrix->NbRows + n ;
  for (i=n;i<end;i++)
  { current = i + decal ;
    Flag(p,current) = Unknown ;
    Denom(p,current) = UN ;
    entier = *(matrix->p + i - n) ;
    /* Pour passer l'indicateur d'egalite/inegalite. */
    inequality = *entier ;
    entier ++ ;
         
    /* Dans le format de la polylib, l'element constant est place en
     * dernier. Dans le format de Pip, il se trouve apres la premiere
     * serie de variables (inconnues ou parametres). On remet donc les
     * choses dans l'ordre de Pip. Ici pour p(x) >= 0.
     */
    for (j=0;j<Nv;j++)
    *(p->row[current].objet.val + j) = *entier++ ;
    for (j=Nv+1;j<nb_columns;j++)
    *(p->row[current].objet.val + j) = *entier++ ;
    *(p->row[current].objet.val + Nv) = *entier ;
    
    /* Et ici lors de l'ajout de -p(x) >= 0 quand on traite une egalite. */
    if (!inequality)
    { decal ++ ;
      new = current + 1 ;
      Flag(p,new)= Unknown ;
      Denom(p,new) = UN ;
      
      for (j=0;j<nb_columns;j++)
      *(p->row[new].objet.val + j) = -(*(p->row[current].objet.val + j)) ;
    }
  }
  return(p);
}


char *Attr[] = {"Unit", "+", "-", "0", "*", "?"};

void tab_display(p, foo)
FILE *foo;
Tableau *p;
{

 int i, j, ff, fff, n;
 Entier x, d;
 fprintf(foo, "%ld/[%d * %d]\n", cross_product, p->height, p->width);
 for(i = 0; i<p->height; i++){
   fff = ff = p->row[i].flags;
   /* if(fff ==0) continue; */
   d = Denom(p, i);
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
       x = Index(p,i,j);
       fprintf(foo, FORMAT, x);
       fprintf(foo, " ");
     }
   fprintf(foo, "]/");
   fprintf(foo, "%d", (int)d);
   putc('\n', foo);
 }
}
