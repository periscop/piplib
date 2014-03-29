/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                   sol.h                                    *
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
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pip.h"

extern long int cross_product, limit;
extern int verbose;
extern FILE *dump;

struct S
    {int flags;
     piplib_int_t param1, param2;
    };

#define Free 0
#define Nil  1
#define If   2
#define List 3
#define Form 4
#define New  5
#define Div  6
#define Val  7
#define Error 8

struct S * sol_space;
static int sol_free;

void sol_init(void)
{
 sol_free = 0;
 sol_space = (struct S *)malloc(SOL_SIZE*sizeof(struct S)) ;
}

void sol_close(void)
{
 free(sol_space) ;
}

int sol_hwm()
{
 return(sol_free);
}

void sol_reset(p)
int p;
{int i;
 if(p<0 || p>=SOL_SIZE)
     {fprintf(stderr, "Syserr : sol_reset : Memory allocation error\n");
      exit(40);
     }
 for(i=p; i<sol_free; i++){
   piplib_int_clear(sol_space[i].param1);
   piplib_int_clear(sol_space[i].param2);
 }
 sol_free = p;
}

struct S *sol_alloc(void)
{struct S *r;
 r = sol_space + sol_free;
 r->flags = Free;
 piplib_int_init_set_si(r->param1,0);
 piplib_int_init_set_si(r->param2,0);
 sol_free++;
 if(sol_free >= SOL_SIZE)
     {fprintf(stderr, "The solution is too complex! : sol\n");
      exit(26);
     }
     return(r);
}

void sol_nil(void)
{
 struct S * r;
 r = sol_alloc();
 r -> flags = Nil;
 if(verbose > 0)
   {fprintf(dump, "\nNil");
    fflush(dump);
  }
}

void sol_error(int c)
{
 struct S *r;
 r = sol_alloc();
 r->flags = Nil;
 piplib_int_set_si(r->param1, c);
 if(verbose > 0) {
     fprintf(dump, "Erreur %d\n", c);
     fflush(dump);
     }
}

int is_not_Nil(p)
int p;
{
 return(sol_space[p].flags != Nil);
}

void sol_if(void)
{
 struct S *r;
 r = sol_alloc();
 r -> flags = If;
 if(verbose > 0) {
     fprintf(dump, "\nIf ");
     fflush(dump);
   }
}

void sol_list(n)
int n;
{struct S * r;
 r = sol_alloc();
 r->flags = List;
 piplib_int_set_si(r->param1, n);
 if(verbose > 0) {
     fprintf(dump, "\nList %d ", n);
     fflush(dump);
}
}

void sol_forme(l)
int l;
{
 struct S *r;
 r = sol_alloc();
 r -> flags = Form;
 piplib_int_set_si(r -> param1, l);
 if(verbose > 0) {
     fprintf(dump, "\nForme %d ", l);
     fflush(dump);
   }
}

void sol_new(k)
int k;
{
 struct S *r;
 r = sol_alloc();
 r -> flags = New;
 piplib_int_set_si(r -> param1, k);
 if(verbose > 0) {
     fprintf(dump, "New %d ", k);
     fflush(dump);
   }
}

void sol_div()
{
 struct S *r;
 r = sol_alloc();
 r -> flags = Div;
 if(verbose > 0) {
     fprintf(dump, "Div ");
     fflush(dump);
   }
}

void sol_val(n, d)
piplib_int_t n, d;
{
 struct S *r;
 r = sol_alloc();
 r -> flags = Val;
 piplib_int_assign(r->param1, n);
 piplib_int_assign(r->param2, d);
 if(verbose > 0) {
   fprintf(dump, "val(");
   piplib_int_print(dump, n);
   fprintf(dump, "/");
   piplib_int_print(dump, d);
   fprintf(dump, ") ");
   fflush(dump);
  }
}

int skip(int);

/* a` partir d'un point de la solution, sauter un objet
bien forme' ainsi qu'un e'ventuel New et pointer sur l'objet
suivant */

int skip_New (int i)
{
 if(sol_space[i].flags != New) return i;
 i = skip(i+1);      /* sauter le Div */
 return i;
}
/* au lancement, i indexe une cellule qui est la te^te d'un objet.
   la valeur retourne'e est la te^te de l'objet qui suit. Les
   objets de type New sont e'limine's                                */

int skip (int i)
{int n, f;
 while((f = sol_space[i].flags) == Free || f == Error) i++;
 switch (sol_space[i].flags) {
 case Nil : case Val : i++; break;
 case New : i = skip_New(i); break;
 case If : i = skip(i+1);        /* sauter le pre'dicat */
	   i = skip(i);          /* sauter le vrai */
	   i = skip(i); break;   /* sauter le faux */
 case List : case Form :
           n = piplib_int_get_si(sol_space[i].param1);
	   i++;
	   while(n--) i = skip(i);
	   break;
 case Div : i = skip(i+1);       /* sauter la forme */
	    i = skip(i);         /* sauter le diviseur */
	    break;
 default : fprintf(stderr,
	      "Syserr : skip : unknown %d\n", sol_space[i].flags);
 }
 return skip_New(i);
}
/* simplification de la solution : e'limination des constructions
   (if p () ()). N'est en service qu'en pre'sence de l'option -z */

void sol_simplify(int i)
{int j, k, l;
 if(sol_space[i].flags == If) {
     j = skip(i+1);        /* j : debut de la partie vraie */
     k = skip(j);          /* k : debut de la partie fausse */
     sol_simplify(k);
     sol_simplify(j);
     if(sol_space[j].flags == Nil && sol_space[k].flags == Nil) {
	 sol_space[i].flags = Nil;
	 if (k >= sol_free - 1) 
	    sol_reset(i+1);
	 else for(l = i+1; l<=k; l++) sol_space[l].flags = Free;
       }
   }

}
/* e'dition de la solution */

int sol_edit(FILE *foo, int i)
{int j, n;
 struct S *p;
 piplib_int_t N, D, d;
 piplib_int_init(N);
 piplib_int_init(D);
 piplib_int_init(d);
 
 p = sol_space + i;
 for(;;) {
   if(p->flags == Free) {
     p++;
     i++;
     continue;
   }
   if(p->flags == New) {
     n = piplib_int_get_si(p->param1);
     fprintf(foo, "(newparm %d ", n);
     if(verbose>0)fprintf(dump, "(newparm %d ", n);
     i = sol_edit(foo, ++i);
     p = sol_space +i;
     fprintf(foo, ")\n");
     if(verbose>0)fprintf(dump, ")\n");
     continue;
   }
   break;
 }
 switch(p->flags){
 case Nil : fprintf(foo, "()\n");
   if(verbose>0)fprintf(dump, "()\n");
   i++; break;
 case Error :
   fprintf(foo, "Error %d\n", piplib_int_get_si(p->param1));
   if(verbose>0)
   fprintf(dump, "Error %d\n", piplib_int_get_si(p->param1));
   i++; break;
 case If  : fprintf(foo, "(if ");
   if(verbose>0)fprintf(dump, "(if ");
   i = sol_edit(foo, ++i);
   i = sol_edit(foo, i);
   i = sol_edit(foo, i);
   fprintf(foo, ")\n");
   if(verbose>0)fprintf(dump, ")\n");
   break;
 case List: fprintf(foo, "(list ");
   if(verbose>0)fprintf(dump, "(list ");
   n = piplib_int_get_si(p->param1);
   i++;
   while(n--) i = sol_edit(foo, i);
   fprintf(foo, ")\n");
   if(verbose>0)fprintf(dump, ")\n");
   break;
 case Form: fprintf(foo, "#[");
   if(verbose>0)fprintf(dump, "#[");
   n = piplib_int_get_si(p->param1);
   for(j = 0; j<n; j++){
     i++; p++;
     piplib_int_assign(N, p->param1);
     piplib_int_assign(D, p->param2);
     piplib_int_gcd(d, N, D);
     if(piplib_int_eq(d, D)) {
       putc(' ', foo);
       piplib_int_div_exact(N, N, d);
       piplib_int_print(foo, N);
       if(verbose>0){
         putc(' ', dump);
         piplib_int_print(dump, N);
       }
     }
     else{
       piplib_int_div_exact(N, N, d);
       piplib_int_div_exact(D, D, d);
       putc(' ', foo);
       piplib_int_print(foo, N);
       putc('/', foo);
       piplib_int_print(foo, D);
       if(verbose>0){
         putc(' ', dump);
         piplib_int_print(dump, N);
         putc('/', dump);
         piplib_int_print(dump, D);
       }
     }
   }
   fprintf(foo, "]\n");
   if(verbose>0)fprintf(dump, "]\n");
   i++;
   break;
 case Div : fprintf(foo, "(div ");
   if(verbose>0)fprintf(dump, "(div ");
   i = sol_edit(foo, ++i);
   i = sol_edit(foo, i);
   fprintf(foo, ")\n");
   if(verbose>0)fprintf(dump, ")\n");
   break;
 case Val :
   piplib_int_assign(N, p->param1);
   piplib_int_assign(D, p->param2);
   piplib_int_gcd(d, N, D);
   if (piplib_int_eq(d, D)) {
     piplib_int_div_exact(N, N, d);
     putc(' ', foo);
     piplib_int_print(foo, N);
     if(verbose>0){
       putc(' ', dump);
       piplib_int_print(dump, N);
     }
   }
   else{
     piplib_int_div_exact(N, N, d);
     piplib_int_div_exact(D, D, d);
     putc(' ', foo);
     piplib_int_print(foo, N);
     fprintf(foo, "/");
     piplib_int_print(foo, D);
     if(verbose>0){
       putc(' ', dump);
       piplib_int_print(dump, N);
       fprintf(dump, "/");
       piplib_int_print(dump, D);
     }
   }
   i++;
   break;
 default  : fprintf(foo, "Inconnu : sol\n");
   if(verbose>0)fprintf(dump, "Inconnu : sol\n");
 }
 piplib_int_clear(d);
 piplib_int_clear(D);
 piplib_int_clear(N);
 return(i);
}


/* Fonction sol_vector_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a un Vector dans la grammaire dans une structure de type PipVector. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau sol_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipVector contenant les informations de ce Vector.
 * Premiere version : Ced. 20 juillet 2001. 
 */
PipVector * sol_vector_edit(int *i, int Bg, int Urs_p, int flags)
{ int j, k, n, unbounded  = 0, first_urs;
  struct S *p ;
  piplib_int_t N, D, d ;
  PipVector * vector ;

  piplib_int_init(N);
  piplib_int_init(D);
  piplib_int_init(d);
  
  vector = (PipVector *)malloc(sizeof(PipVector)) ;
  if (vector == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  p = sol_space + (*i) ;
  n = piplib_int_get_si(p->param1);
  if (flags & SOL_REMOVE)
    --n;
  n -= Urs_p;
  first_urs = Urs_p + (Bg >= 0);
  vector->nb_elements = n ;
  vector->the_vector = (piplib_int_t *)malloc(sizeof(piplib_int_t)*n) ;
  if (vector->the_vector == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  vector->the_deno = (piplib_int_t *)malloc(sizeof(piplib_int_t)*n) ;
  if (vector->the_deno == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  
  for (j=0, k=0; k < n; j++) {
    (*i)++ ;
    p++ ;

    piplib_int_assign(N, p->param1);
    piplib_int_assign(D, p->param2);
    piplib_int_gcd(d, N, D);

    if ((flags & SOL_SHIFT) && j == Bg) {
      piplib_int_sub(N, N, D);   /* subtract 1 */
      if (piplib_int_zero(N) == 0)
	unbounded = 1;
    }

    if ((flags & SOL_REMOVE) && j == Bg)
      continue;

    if (first_urs <= j && j < first_urs+Urs_p)
      continue;

    piplib_int_init(vector->the_vector[k]);
    piplib_int_div_exact(vector->the_vector[k], N, d);
    if (flags & SOL_NEGATE)
      piplib_int_oppose(vector->the_vector[k], vector->the_vector[k]);
    piplib_int_init(vector->the_deno[k]);
    if (piplib_int_eq(d, D))
      piplib_int_assign(vector->the_deno[k], UN);
    else
      piplib_int_div_exact(vector->the_deno[k], D, d);
    ++k;
  }
  if (unbounded)
    for (k=0; k < n; k++)
      piplib_int_assign(vector->the_deno[k], ZERO);
  (*i)++ ;

  piplib_int_clear(d);
  piplib_int_clear(D);
  piplib_int_clear(N);

  return(vector) ;
}


/* Fonction sol_newparm_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a un Newparm dans la grammaire dans une structure de type PipNewparm. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau sol_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipNewparm contenant les informations de ce Newparm.
 * Premiere version : Ced. 18 octobre 2001. 
 */
PipNewparm * sol_newparm_edit(int *i, int Bg, int Urs_p, int flags)
{ struct S * p ;
  PipNewparm * newparm, * newparm_first = NULL, * newparm_now = NULL;

  /* On place p au lieu de lecture. */
  p = sol_space + (*i) ;

  do {
    /* On passe le New et le Div pour aller a Form et lire le VECTOR. */
    (*i) += 2 ;

    newparm = (PipNewparm *)malloc(sizeof(PipNewparm)) ;
    if (newparm == NULL)
    { fprintf(stderr, "Memory Overflow.\n") ;
      exit(1) ;
    }
    newparm->vector = sol_vector_edit(i, Bg, Urs_p, flags);
    newparm->rank = piplib_int_get_si(p->param1);
    /* On met p a jour pour lire le denominateur (un Val de param2 UN). */
    p = sol_space + (*i) ;
    piplib_int_init(newparm->deno);
    piplib_int_assign(newparm->deno, p->param1);
    if (flags & SOL_REMOVE)
      newparm->rank--;
    newparm->rank -= Urs_p;
    newparm->next = NULL ;

    if (newparm_now)
      newparm_now->next = newparm;
    else
      newparm_first = newparm;
    newparm_now = newparm ;
    if (verbose > 0)
    { fprintf(dump,"\n(newparm ") ;
      fprintf(dump,piplib_int_format,newparm->rank) ;
      fprintf(dump," (div ") ;
      pip_vector_print(dump,newparm->vector) ;
      fprintf(dump," ") ;
      piplib_int_print(dump, newparm->deno);
      fprintf(dump,"))") ;
    }
  
    /* On passe aux elements suivants. */
    (*i) ++ ;
    p = sol_space + (*i) ;
  } while (p->flags == New);

  return newparm_first;
}


/* Fonction sol_list_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a une List dans la grammaire dans une structure de type PipList. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau sol_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipList contenant les informations de cette List.
 * Premiere version : Ced. 18 octobre 2001. 
 * 16 novembre 2005 : Ced. Prise en compte du cas 0 éléments, avant impossible.
 */
PipList * sol_list_edit(int *i, int nb_elements, int Bg, int Urs_p, int flags)
{ PipList * list, * list_new, * list_now ;
  
  /* Pour le premier element. */
  list = (PipList *)malloc(sizeof(PipList)) ;
  if (list == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  list->next = NULL ;
  
  if (nb_elements == 0)
  { list->vector = NULL ;
    return(list) ;
  }
  
  list->vector = sol_vector_edit(i, Bg, Urs_p, flags);

  list_now = list ;
  if (verbose > 0)
  { fprintf(dump,"\n(list ") ;
    pip_vector_print(dump,list->vector) ;
  }
  nb_elements-- ;

  /* Pour les elements suivants. */
  while (nb_elements--)
  { list_new = (PipList *)malloc(sizeof(PipList)) ;
    if (list_new == NULL)
    { fprintf(stderr, "Memory Overflow.\n") ;
      exit(1) ;
    }
    list_new->vector = sol_vector_edit(i, Bg, Urs_p, flags);
    list_new->next = NULL ;
		    
    if (verbose > 0)
    { fprintf(dump,"\n") ;
      pip_vector_print(dump,list_new->vector) ;
    }
    list_now->next = list_new ;
    list_now = list_now->next ;
  }
  if (verbose > 0)
  fprintf(dump,"\n)") ;
  
  return(list) ;
}


/* Fonction sol_quast_edit :
 * Cette fonction a pour but de placer les informations de la solution
 * (qui sont contenues dans le tableau sol_space) dans une structure de
 * type PipQuast en vue d'une utilisation directe de la solution par une
 * application exterieure. Elle prend en parametre un pointeur vers une
 * case memoire contenant le numero de cellule du tableau sol_space
 * a partir de laquelle on doit commencer la lecture des informations. Elle
 * recoit aussi l'adresse du PipQuast qui l'a appelle (pour le champ father).
 * Elle retourne un pointeur vers une structure de type PipQuast qui
 * contient toutes les informations sur la solution (sous forme d'arbre).
 * Remarques : cette fonction lit les informations comme elles doivent
 * se presenter a la fin du traitement. Elle respecte scrupuleusement
 * la grammaire attendue et n'accepte de passer des cellules a Free
 * qu'entre une des trois grandes formes (if, list ou suite de newparm).
 * 20  juillet 2001 : Premiere version, Ced. 
 * 31  juillet 2001 : Ajout du traitement de l'option verbose = code*2 :0( 
 * 18  octobre 2001 : Grands changements dus a l'eclatement de la structure
 *                    PipVector en PipVector, PipNewparm et PipList, et
 *                    eclatement de la fonction avec sol_newparm_edit et
 *                    sol_list_edit.
 * 16 novembre 2005 : (debug) Même si une liste est vide il faut la créer pour
 *                    afficher plus tard le (list), repéré par Sven Verdoolaege.
 */
PipQuast *sol_quast_edit(int *i, PipQuast *father, int Bg, int Urs_p, int flags)
{ int nb_elements ;
  struct S * p ;
  PipQuast * solution ;
  PipList * list_new, * list_now ;
  PipNewparm * newparm_new, * newparm_now ;
    
  /* On place p au lieu de lecture. */
  p = sol_space + (*i) ;
  /* En cas d'utilisation de l'option de simplification, une plage de
   * structures S peut avoir les flags a Free. On doit alors les passer.
   */
  while (p->flags == Free)
  { p ++ ;
    (*i) ++ ;
  }
  
  solution = (PipQuast *)malloc(sizeof(PipQuast)) ;
  if (solution == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  solution->newparm = NULL ;
  solution->list = NULL ;
  solution->condition = NULL ;
  solution->next_then = NULL ;
  solution->next_else = NULL ;
  solution->father = father ;
  
  /* On peut commencer par une chaine de nouveaux parametres... */
  if (p->flags == New)
  { solution->newparm = sol_newparm_edit(i, Bg, Urs_p, flags & SOL_REMOVE);
    p = sol_space + (*i) ;
  }
  
  /* ...ensuite soit par une liste (vide ou non) soit par un if. */
  (*i)++ ; /* Factorise de List, Nil et If. */
  switch (p->flags)
  { case List :
                nb_elements = piplib_int_get_si(p->param1) ;
                solution->list = sol_list_edit(i, nb_elements, Bg, Urs_p, flags);
		if (flags & SOL_DUAL)
		    solution->next_then = sol_quast_edit(i, solution, Bg, Urs_p, 0);
		break ;
    case Nil  : if (verbose > 0)
		fprintf(dump,"\n()") ;
                break ;
    case If   : solution->condition = 
			    sol_vector_edit(i, Bg, Urs_p, flags & SOL_REMOVE);
                if (verbose > 0)
		{ fprintf(dump,"\n(if ") ;
                  pip_vector_print(dump,solution->condition) ;
                }
		solution->next_then = sol_quast_edit(i, solution, Bg, Urs_p, flags);
                solution->next_else = sol_quast_edit(i, solution, Bg, Urs_p, flags);
                if (verbose > 0)
		fprintf(dump,"\n)") ;
                break ;
    default   : fprintf(stderr,"\nAie !!! Flag %d inattendu.\n",p->flags) ;
                if (verbose > 0)
		fprintf(dump,"\nAie !!! Flag %d inattendu.\n",p->flags) ;
                exit(1) ;
  }
  
  return(solution) ;
}
