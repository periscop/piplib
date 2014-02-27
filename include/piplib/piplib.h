/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 piplib.h                                   *
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
 * Written by Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/

#ifndef PIPLIB_H
#define PIPLIB_H

#include <stdio.h>

#include <osl/int.h>


#if defined(__cplusplus)
extern "C"
  {
#endif


/* Structure PipMatrix :
 * Structure de matrice au format PolyLib. Le premier element d'une ligne
 * indique quand il vaut 1 que la ligne decrit une inequation de la forme
 * p(x)>=0 et quand il vaut 0, que la ligne decrit une egalite de la forme
 * p(x)=0. Le dernier element de chaque ligne correspond au coefficient
 * constant. */
struct pipmatrix {
  unsigned NbRows, NbColumns ;
  osl_int_t** p;
  osl_int_t* p_Init;
  int p_Init_size;   /* Only for PolyLib compatibility under MP
                      * version: PolyLib makes sometimes
                      * overestimates on the size of the matrices,
                      * in order to go faster. Thus
                      * NbRows*NbColumns is not the number of
                      * allocated elements. With MP version, we
                      * have to think to osl_int_clear(PIPLIB_INT_PRECISION, &) all the
                      * initialized elements before freing, then
                      * we need to know the number of allocated
                      * elements: p_Init_size. */
};

typedef struct pipmatrix PipMatrix;


/* Structure PipVector :
 * Cette structure contient un Vector de 'nb_elements' la ieme composante de
 * ce vecteur vaut the_vector[i]/the_deno[i]. */
struct pipvector {
  int nb_elements;       /* Nombre d'elements du vecteur. */
  osl_int_t* the_vector; /* Numerateurs du vecteur. */
  osl_int_t* the_deno;   /* Denominateurs du vecteur. */
};

typedef struct pipvector PipVector ;


/* Structure PipNewparm :
 * Liste chainee de Newparm, les informations d'un newparm etant son rang, un
 * vecteur de coefficients et un denominateur. Le newparm est egal a la division
 * du vecteur par le denominateur. */
struct pipnewparm {
  int rank;                /* Rang du 'newparm'. */
  PipVector* vector;       /* Le vector decrivant le newparm. */
  osl_int_t deno;          /* Denominateur du 'newparm'. */
  struct pipnewparm* next; /* Pointeur vers le newparm suivant. */
};

typedef struct pipnewparm PipNewparm;


/* Structure PipList : Liste chainee de Vector. */
struct piplist {
  PipVector* vector;    /* Le vector contenant la partie de solution. */
  struct piplist* next; /* Pointeur vers l'element suivant. */
};

typedef struct piplist PipList;


/* Structure pipquast :
 * Arbre binaire. Conformement a la grammaire de sortie (voir mode d'emploi), un
 * noeud de l'arbre des solutions debute par une liste de 'newparm'. Il continue
 * ensuite soit par une 'list' (alors condition vaut null), soit par un 'if'
 * (alors le champ condition contient la condition). */
struct pipquast {
  PipNewparm* newparm;        /* Les 'newparm'. */
  PipList* list;              /* La 'list' si pas de 'if'. */
  PipVector* condition;       /* La condition si 'if'. */
  struct pipquast* next_then; /* Noeud si condition et si verifiee. */
  struct pipquast* next_else; /* Noeud si condition et si non verifiee. */
  struct pipquast* father;    /* Pointeur vers le quast pere. */
};

typedef struct pipquast PipQuast;


/* Structure pipoptions:
 * This structure contains each option that can be set to change the PIP
 * behaviour. */
struct pipoptions {
  int Nq;           /* 1 if an integer solution is needed, 0 otherwise. */
  int Verbose;      /* -1 -> absolute silence
                        0 -> relative silence
                        1 -> information on cuts
                             when an integer solution is needed
                        2 -> information sur les pivots et les d�terminants
                        3 -> information on arrays
                        Each option include the preceding. */
  int Simplify ;    /* Set to 1 to eliminate some trivial
                       solutions, 0 otherwise. */
  int Deepest_cut;  /* Set to 1 to include deepest cut algorithm. */
  int Maximize;     /* Set to 1 if maximum is needed. */
  int Urs_parms;    /* -1 -> all parameters may be negative
                        0 -> all parameters are non-negative */
  int Urs_unknowns; /* -1 -> all unknowns may be negative
                        0 -> all unknowns are non-negative */
  int Compute_dual;
};

typedef struct pipoptions PipOptions;

void* pip_options_print(FILE*, PipOptions*);


/* Prototypes des fonctions d'affichages des structures */
void pip_matrix_print(FILE*, PipMatrix*);
void pip_vector_print(FILE*, PipVector*);
void pip_newparm_print(FILE*, PipNewparm*, int);
void pip_list_print(FILE*, PipList*, int);
void pip_quast_print(FILE*, PipQuast*, int);


/* Prototypes des fonctions de liberation memoire des structures */
void pip_matrix_free(PipMatrix*);
void pip_vector_free(PipVector*);
void pip_newparm_free(PipNewparm*);
void pip_list_free(PipList*);
void pip_quast_free(PipQuast*);
void pip_options_free(PipOptions*);


/* Prototypes des fonctions d'acquisition de matrices de contraintes et
   options. */
PipMatrix* pip_matrix_alloc(unsigned, unsigned);
PipMatrix* pip_matrix_read(FILE*);
PipOptions* pip_options_init(void);
 

/* Initialization of pip library */
void pip_init();
void pip_init_sp();
void pip_init_dp();
void pip_init_mp();
void pip_close();


/* Prototype de la fonction de resolution :
 * pip_solve resoud le probleme qu'on lui passe en parametre, suivant les
 * options elles aussi en parametre. Elle renvoie la solution sous forme
 * d'un arbre de PipQuast. Parametres :
 * - probleme :
 * 1 PipMatrix  : systeme des inequations definissant le domaine des inconnues,
 * 2 PipMatrix  : systeme des inequations satisfaites par les parametres,
 * 3 int        : column rank of the bignum, or negative value if there
 *                is no big parameter.
 * 4 PipOptions : options for PIP. */ 
PipQuast* pip_solve(PipMatrix*, PipMatrix*, int, PipOptions*) ;

#define SOL_SHIFT  (1 << 0)                 /* Shift solution over -bigparam */
#define SOL_NEGATE (1 << 1)                 /* Negate solution */
#define SOL_REMOVE (1 << 2)                 /* Remove big parameter */
#define SOL_MAX    (SOL_SHIFT | SOL_NEGATE) /* Maximum was computed */
#define SOL_DUAL   (1 << 3)                 /* Create dual leaf */
PipQuast* sol_quast_edit(int*, PipQuast*, int, int, int);

#if defined(__cplusplus)
  }
#endif

#endif /* PIPLIB_H */
