/********************************************************/
/* Part of MultiPrecision PIP port (C) Zbigniew Chamski */
/* and Paul Feautrier, 2001.                            */
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

#define TAB_CHUNK 4096*sizeof(Entier)

static char *tab_free, *tab_top;
static struct A *tab_base;

extern int allocation;
extern long int cross_product, limit;
static int chunk_count;

extern FILE * dump;

/*   Structure of the heap after execution of tab_init.

     tab_base         tab_free         tab_top
        |                |                |
        |                |                |
        V                |                |
     struct A {          V                |
NULL <---- precedent;    .  <-------------.
           bout;
           } 
*/

void tab_init(void)
{
 tab_free = malloc(sizeof (struct A));
 if(tab_free == NULL)
     {fprintf(stdout, "Your computer doesn't have memory\n");
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

void tab_reset(struct high_water_mark by_the_mark)

{struct A *g;
 char *p;
 while(chunk_count > by_the_mark.chunk)
     {
      g = tab_base->precedent;
/* Before actually freeing the memory, one has to clear the included Tableaux.
   If this is not done, the GMP objects referenced in the Tableaux will be
   orphaned.
*/

/* Enumerate the included tableaux.
 */
      p = (char *)tab_base + sizeof(struct A);
      while(p < tab_free){
        Tableau *pt;
        pt = (Tableau *) p;
	tab_clear(pt);
        p += pt->taille;
      } 
      free(tab_base);
      tab_base = g;
      tab_free = tab_base->bout;
      chunk_count--;
     }
 if(chunk_count > 0) tab_free = by_the_mark.top;
 else {
     fprintf(stdout, "Syserr: tab_reset : error in memory allocation\n");
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
 unsigned long the_taille;
 int i, j;
 the_taille = sizeof(struct T) + (h+n-1) * sizeof (struct L)
	  + h * w * sizeof (Entier);
 if(tab_free + the_taille >= tab_top)
     {struct A * g;
      unsigned long d;
      d = the_taille + sizeof(struct A);
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
 tab_free += the_taille;
 tp = (Tableau *)p;
 q = (Entier *)(p +  sizeof(struct T) + (h+n-1) * sizeof (struct L));
 mpz_init_set_ui(tp->determinant, 1);
 for(i = 0; i<n ; i++){
   tp->row[i].flags = Unit;
   tp->row[i].objet.unit = i;
   mpz_init_set_ui(Denom(tp, i), 1);
 }
 for(i = n; i < (h+n); i++){
   tp->row[i].flags = 0;
   tp->row[i].objet.val = q;
   for(j = 0; j < w; j++)
     mpz_init_set_ui(*q++, 0);
   mpz_init_set_ui(Denom(tp, i), 0);
 }
 tp->height = h + n; tp->width = w; tp->taille = the_taille;
 return(tp);
}

Tableau * tab_get(foo, h, w, n)
FILE * foo;
int h, w, n;
{
 Tableau *p;
 int i, j, c;
 int x;
 p = tab_alloc(h, w, n);
 while((c = dgetc(foo)) != EOF)
      if(c == '(')break;
 for(i = n; i<h+n; i++){
   p->row[i].flags = Unknown;
   mpz_set_ui(Denom(p, i), 1);
   while((c = dgetc(foo)) != EOF)if(c == '[')break;
   for(j = 0; j<w; j++){
     if(dscanf(foo, FORMAT, &x) < 0)
       return NULL;
     else mpz_set_si(p->row[i].objet.val[j], x);
   }
 }
 while((c = dgetc(foo)) != EOF)if(c == ']')break;
 
 return(p);
}


char *Attr[] = {"Unit", "+", "-", "0", "*", "?"};

void tab_display(p, foo)
FILE *foo;
Tableau *p;
{

 int i, j, ff, fff, n;
 Entier d;
 mpz_init(d);
 fprintf(foo, "%ld/[%d * %d]\n", cross_product, p->height, p->width);
 for(i = 0; i<p->height; i++){
   fff = ff = p->row[i].flags;
   if(fff == 0) continue;
   mpz_set(d, Denom(p, i));
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
       mpz_out_str(foo, 10, Index(p, i, j));
       putc(' ', foo);
     }
   fprintf(foo, "]/");
   mpz_out_str(foo, 10, d);
   putc('\n', foo);
 }
 mpz_clear(d);
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
 * 24 octobre 2002 : Premiere version MP. 
 */
Tableau * tab_Matrix2Tableau(PipMatrix * matrix, int Nineq, int Nv, int n)
{ Tableau * p ;
  unsigned i, j, end, current, new, nb_columns, decal=0 ;
  Entier * entier, inequality ;
  
  mpz_init(inequality) ;
  nb_columns = matrix->NbColumns - 1 ;
  p = tab_alloc(Nineq,nb_columns,n) ;
  
  /* La variable decal sert a prendre en compte les lignes supplementaires
   * issues des egalites.
   */
  end = matrix->NbRows + n ;
  for (i=n;i<end;i++)
  { current = i + decal ;
    Flag(p,current) = Unknown ;
    mpz_set_ui(Denom(p,current),1) ;
    entier = *(matrix->p + i - n) ;
    /* Pour passer l'indicateur d'egalite/inegalite. */
    mpz_set(inequality,*entier) ;
    entier ++ ;
         
    /* Dans le format de la polylib, l'element constant est place en
     * dernier. Dans le format de Pip, il se trouve apres la premiere
     * serie de variables (inconnues ou parametres). On remet donc les
     * choses dans l'ordre de Pip. Ici pour p(x) >= 0.
     */
    for (j=0;j<Nv;j++)
    mpz_set(*(p->row[current].objet.val + j),*entier++) ;
    for (j=Nv+1;j<nb_columns;j++)
    mpz_set(*(p->row[current].objet.val + j),*entier++) ;
    mpz_set(*(p->row[current].objet.val + Nv),*entier) ;
    
    /* Et ici lors de l'ajout de -p(x) >= 0 quand on traite une egalite. */
    if (!inequality)
    { decal ++ ;
      new = current + 1 ;
      Flag(p,new)= Unknown ;
      mpz_set(Denom(p,new),UN) ;
      
      for (j=0;j<nb_columns;j++)
      mpz_neg(*(p->row[new].objet.val + j),*(p->row[current].objet.val + j)) ;
    }
  }
  mpz_clear(inequality);
  return(p);
}

