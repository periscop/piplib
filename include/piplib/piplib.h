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
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this library; if not, write to the Free Software Foundation,    *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA         *
 *                                                                            *
 * Written by Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/

/* Premiere version du 18 septembre 2002. */


#ifndef PIPLIB_H
#define PIPLIB_H


#include <stdio.h>
#include <math.h>


// Copy from osl_int
long long int piplib_llgcd(long long int const, long long int const);
long long int piplib_llgcd_llabs(long long int const, long long int const);
size_t piplib_lllog2(long long int);
size_t piplib_lllog10(long long int);
long long int piplib_llmod(long long int const, long long int const);


#ifdef PIPLIB_INT_GMP

  #include <gmp.h>

  #define PIPLIB_ONE_DETERMINANT

  typedef mpz_t piplib_int_t;
  #define piplib_int_format "%d"

  #define piplib_int_init(i) (mpz_init(i))
  #define piplib_int_init_set(i, v) (mpz_init_set(i, v))
  #define piplib_int_init_set_si(i, v) (mpz_init_set_si(i, v))
  #define piplib_int_assign(r, i) (mpz_set(r, i))
  #define piplib_int_set_si(r, i) (mpz_set_si(r, i))
  #define piplib_int_clear(i) (mpz_clear(i))
  #define piplib_int_print(file, i) (mpz_out_str(file, 10, i))
  #define piplib_int_sscanf(string, i) (gmp_sscanf(string, "%lZd", i))

  #define piplib_int_get_si(i) ((int)mpz_get_si(i))
  #define piplib_int_get_d(i) (mpz_get_d(i))

  #define piplib_int_add(r, a, b) (mpz_add(r, a, b))
  #define piplib_int_sub(r, a, b) (mpz_sub(r, a, b))
  #define piplib_int_increment(r, i) (mpz_add_ui(r, i, 1))
  #define piplib_int_decrement(r, i) (mpz_sub_ui(r, i, 1))
  #define piplib_int_mul(r, a, b) (mpz_mul(r, a, b))
  #define piplib_int_div_exact(q, a, b) (mpz_divexact(q, a, b))
  #define piplib_int_floor_div_q(q, a, b) (mpz_fdiv_q(q, a, b))
  #define piplib_int_floor_div_r(r, a, b) (mpz_fdiv_r(r, a, b))
  #define piplib_int_floor_div_q_r(q, r, a, b) (mpz_fdiv_qr(q, r, a, b))
  #define piplib_int_mod(mod, a, b) (mpz_mod(mod, a, b))
  #define piplib_int_gcd(gcd, a, b) (mpz_gcd(gcd, a, b))
  #define piplib_int_oppose(r, i) (mpz_neg(r, i))
  #define piplib_int_size_in_base_2(i) (mpz_sizeinbase(i, 2))
  #define piplib_int_size_in_base_10(i) (mpz_sizeinbase(i, 10))

  #define piplib_int_eq(a, b) (mpz_cmp(a, b) == 0)
  #define piplib_int_ne(a, b) (mpz_cmp(a, b) != 0)
  #define piplib_int_zero(i) (mpz_sgn(i) == 0)
  #define piplib_int_one(i) (mpz_cmp_si(i, 1) == 0)
  #define piplib_int_pos(i) (mpz_sgn(i) > 0)
  #define piplib_int_neg(i) (mpz_sgn(i) < 0)

#elif PIPLIB_INT_OSL

  #include <osl/int.h>

  typedef osl_int_t piplib_int_t;
  #define piplib_int_format "%d"

  #define piplib_int_init(i) (osl_int_init(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_init_set(i, v) (osl_int_init_set(PIPLIB_INT_OSL_PRECISION, i, v))
  #define piplib_int_init_set_si(i, v) (osl_int_init_set_si(PIPLIB_INT_OSL_PRECISION, i, v))
  #define piplib_int_assign(r, i) (osl_int_assign(PIPLIB_INT_OSL_PRECISION, r, i))
  #define piplib_int_set_si(r, i) (osl_int_set_si(PIPLIB_INT_OSL_PRECISION, r, i))
  #define piplib_int_clear(i) (osl_int_clear(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_print(file, i) (osl_int_print(file, PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_sscanf(string, i) (osl_int_sscanf(string, PIPLIB_INT_OSL_PRECISION, i))

  #define piplib_int_get_si(i) (osl_int_get_si(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_get_d(i) (osl_int_get_d(PIPLIB_INT_OSL_PRECISION, i))

  #define piplib_int_add(r, a, b) (osl_int_add(PIPLIB_INT_OSL_PRECISION, r, a, b))
  #define piplib_int_sub(r, a, b) (osl_int_sub(PIPLIB_INT_OSL_PRECISION, r, a, b))
  #define piplib_int_increment(r, i) (osl_int_increment(PIPLIB_INT_OSL_PRECISION, r, i))
  #define piplib_int_decrement(r, i) (osl_int_decrement(PIPLIB_INT_OSL_PRECISION, r, i))
  #define piplib_int_mul(r, a, b) (osl_int_mul(PIPLIB_INT_OSL_PRECISION,r, a, b))
  #define piplib_int_div_exact(q, a, b) (osl_int_div_exact(PIPLIB_INT_OSL_PRECISION, q, a, b))
  #define piplib_int_floor_div_q(q, a, b) (osl_int_floor_div_q(PIPLIB_INT_OSL_PRECISION, q, a, b))
  #define piplib_int_floor_div_r(r, a, b) (osl_int_floor_div_r(PIPLIB_INT_OSL_PRECISION, r, a, b))
  #define piplib_int_floor_div_q_r(q, r, a, b) (osl_int_floor_div_q_r(q, r, a, b))
  #define piplib_int_mod(mod, a, b) (osl_int_mod(PIPLIB_INT_OSL_PRECISION, mod, a, b))
  #define piplib_int_gcd(gcd, a, b) (osl_int_gcd(PIPLIB_INT_OSL_PRECISION, gcd, a, b))
  #define piplib_int_oppose(r, i) (osl_int_oppose(PIPLIB_INT_OSL_PRECISION, r, i))
  #define piplib_int_size_in_base_2(i) (osl_int_size_in_base_2(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_size_in_base_10(i) (osl_int_size_in_base_10(PIPLIB_INT_OSL_PRECISION, i))

  #define piplib_int_eq(a, b) (osl_int_eq(PIPLIB_INT_OSL_PRECISION, a, b))
  #define piplib_int_ne(a, b) (osl_int_ne(PIPLIB_INT_OSL_PRECISION, a, b))
  #define piplib_int_zero(i) (osl_int_zero(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_one(i) (osl_int_one(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_pos(i) (osl_int_pos(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_neg(i) (osl_int_neg(PIPLIB_INT_OSL_PRECISION, i))

#elif PIPLIB_INT_SP

  typedef long int piplib_int_t;
  #define piplib_int_format "%ld"

#elif PIPLIB_INT_DP

  typedef long long int piplib_int_t;
  #define piplib_int_format "%lld"

#else

  #error "Define of PIPLIB_INT_SP or PIPLIB_INT_DP or PIPLIB_INT_GMP or \
          PIPLIB_INT_OSL not found"

#endif


#if defined(PIPLIB_INT_SP) || defined(PIPLIB_INT_DP)

  #define piplib_int_init(i) (i = 0)
  #define piplib_int_init_set(i, v) (i = v)
  #define piplib_int_init_set_si(i, v) (i = v)
  #define piplib_int_assign(r, i) (r = i)
  #define piplib_int_set_si(r, i) (r = (piplib_int_t)i)
  #define piplib_int_clear(i) do { } while (0)
  #define piplib_int_print(file, i) (fprintf(file, piplib_int_format, i))
  #define piplib_int_sscanf(string, i) (sscanf(string, piplib_int_format, &i))

  #define piplib_int_get_si(i) ((int)(i))
  #define piplib_int_get_d(i) ((double)(i))

  #define piplib_int_add(r, a, b) (r = a + b)
  #define piplib_int_sub(r, a, b) (r = a - b)
  #define piplib_int_increment(r, i) (r = i + 1)
  #define piplib_int_decrement(r, i) (r = i - 1)
  #define piplib_int_mul(r, a, b) (r = a * b)
  #define piplib_int_div_exact(q, a, b) (q = (a) / (b))
  //#define piplib_int_floor_div_q(q, a, b) (q = (piplib_int_t)(piplib_ll_floor_div_q(a, b)))
  #define piplib_int_floor_div_q(q, a, b) (q = (piplib_int_t)((a - piplib_llmod(a, b)) / (b)))
  #define piplib_int_floor_div_r(r, a, b) (r = (piplib_int_t)piplib_llmod(a, b))
  //#define piplib_int_floor_div_r(r, a, b) (r = (piplib_int_t)(piplib_ll_floor_div_r(a, b)))
  #define piplib_int_floor_div_q_r(q, r, a, b) do { piplib_int_floor_div_q(q, a, b); piplib_int_floor_div_r(r, a, b); } while (0)
  #define piplib_int_mod(mod, a, b) (mod = (piplib_int_t)(piplib_llmod(a, b)))
  #define piplib_int_gcd(gcd, a, b) (gcd = (piplib_int_t)(piplib_llgcd_llabs(a, b)))
  #define piplib_int_oppose(r, i) (r = - (i))
  #define piplib_int_size_in_base_2(i) (piplib_lllog2(i))
  #define piplib_int_size_in_base_10(i) (piplib_lllog10(i))

  #define piplib_int_eq(a, b) (a == b)
  #define piplib_int_ne(a, b) (a != b)
  #define piplib_int_zero(i) (i == 0)
  #define piplib_int_one(i) (i == 1)
  #define piplib_int_pos(i) (i > 0)
  #define piplib_int_neg(i) (i < 0)

#endif


#if defined(__cplusplus)
extern "C" 
  {
#endif 


/**
 * @brief Structure PipMatrix
 * 
 * Structure de matrice au format PolyLib. Le premier element d'une ligne
 * indique quand il vaut 1 que la ligne decrit une inequation de la forme
 * p(x)>=0 et quand il vaut 0, que la ligne decrit une egalite de la forme
 * p(x)=0. Le dernier element de chaque ligne correspond au coefficient
 * constant.
 */
struct pipmatrix {
  unsigned int NbRows;    /**< Number of rows */
  unsigned int NbColumns; /**< Number of columns */
  piplib_int_t** p;       /**< Data */
  piplib_int_t* p_Init;   /**< Init */
  int p_Init_size;        /**< Only for PolyLib compatibility under MP version:
                               PolyLib makes sometimes overestimates on the size
                               of the matrices, in order to go faster.
                               Thus NbRows*NbColumns is not the number of
                               allocated elements. With MP version, we have to
                               think to mpz_clear() all the initialized elements
                               before freing, then we need to know the number of
                               allocated elements: p_Init_size. */
};
typedef struct pipmatrix PipMatrix;


/**
 * @brief Structure PipVector
 * 
 * Cette structure contient un Vector de 'nb_elements' la ieme composante de
 * ce vecteur vaut the_vector[i]/the_deno[i].
 */
struct pipvector {
  int nb_elements;          /**< Nombre d'elements du vecteur. */
  piplib_int_t* the_vector; /**< Numerateurs du vecteur. */
  piplib_int_t* the_deno;   /**< Denominateurs du vecteur. */
};
typedef struct pipvector PipVector;


/**
 * @brief Structure PipNewparm
 * 
 * Liste chainee de Newparm, les informations d'un newparm etant son rang, un
 * vecteur de coefficients et un denominateur. Le newparm est egal a la division
 * du vecteur par le denominateur.
 */
struct pipnewparm {
  int rank;                /**< Rang du 'newparm'. */
  PipVector* vector;       /**< Le vector decrivant le newparm. */
  piplib_int_t deno;       /**< Denominateur du 'newparm'. */
  struct pipnewparm* next; /**< Pointeur vers le newparm suivant. */
};
typedef struct pipnewparm PipNewparm;


/**
 * @brief Structure PipList
 * 
 * Liste chainee de Vector.
 */
struct piplist {
  PipVector* vector;    /**< Le vector contenant la partie de solution. */
  struct piplist* next; /**< Pointeur vers l'element suivant. */
};
typedef struct piplist PipList;


/**
 * @brief Structure pipquast
 * 
 * Arbre binaire. Conformement a la grammaire de sortie (voir mode d'emploi), un
 * noeud de l'arbre des solutions debute par une liste de 'newparm'. Il continue
 * ensuite soit par une 'list' (alors condition vaut null), soit par un 'if'
 * (alors le champ condition contient la condition).
 */
struct pipquast {
  PipNewparm* newparm;        /**< Les 'newparm'. */
  PipList* list;              /**< La 'list' si pas de 'if'. */
  PipVector* condition;       /**< La condition si 'if'. */
  struct pipquast* next_then; /**< Noeud si condition et si verifiee. */
  struct pipquast* next_else; /**< Noeud si condition et si non verifiee. */
  struct pipquast* father;    /**< Pointeur vers le quast pere. */
};
typedef struct pipquast PipQuast;


/**
 * @brief Structure pipoptions
 * 
 * This structure contains each option that can be set to change the PIP
 * behaviour.
 */
struct pipoptions {
  /** 1 if an integer solution is needed, 0 otherwise. */
  int Nq;
  
  /**
   * @brief Verbose
   * 
   * -1 -> absolute silence @n
   *  0 -> relative silence @n
   *  1 -> information on cuts when an integer solution is needed @n
   *  2 -> information sur les pivots et les déterminants @n
   *  3 -> information on arrays
   * 
   * Each option include the preceding
   */
  int Verbose;
  
  /** Set to 1 to eliminate some trivial solutions, 0 otherwise */
  int Simplify;
  
  /** Set to 1 to include deepest cut algorithm */
  int Deepest_cut;
  
  /** Set to 1 if maximum is needed */
  int Maximize;
  
  /**
   * @brief Signs of parameters
   * 
   * -1 -> all parameters may be negative @n
   *  0 -> all parameters are non-negative
   */
  int Urs_parms;
  
  /**
   * @brief Signs of unknowns
   * 
   * -1 -> all unknowns may be negative @n
   *  0 -> all unknowns are non-negative
   */
  int Urs_unknowns;
  
  /** To compute the dual */
  int Compute_dual;
};
typedef struct pipoptions PipOptions;

void pip_options_print(FILE*, PipOptions*);


/* Fonctions d'affichages des structures de la PipLib. */
void pip_matrix_print(FILE*, PipMatrix*);
void pip_vector_print(FILE*, PipVector*);
void pip_newparm_print(FILE*, PipNewparm*, int);
void pip_list_print(FILE*, PipList*, int);
void pip_quast_print(FILE*, PipQuast*, int);


/* Fonctions de liberation memoire des structures de la PipLib.*/
void pip_matrix_free(PipMatrix*);
void pip_vector_free(PipVector*);
void pip_newparm_free(PipNewparm*);
void pip_list_free(PipList*);
void pip_quast_free(PipQuast*);
void pip_options_free(PipOptions*);


/* Fonctions d'acquisition de matrices de contraintes et options. */
PipMatrix* pip_matrix_alloc(unsigned int, unsigned int);
PipMatrix* pip_matrix_read(FILE*);
PipOptions* pip_options_init(void);


/* initialization of pip library */
void pip_init();
void pip_close();


/**
 * @brief Fonction de resolution
 * 
 * pip_solve resoud le probleme qu'on lui passe en parametre, suivant les
 * options elles aussi en parametre. Elle renvoie la solution sous forme
 * d'un arbre de PipQuast.
 * 
 * @param[in] domain     Inequations definissant le domaine des inconnues
 * @param[in] parameters Inequations satisfaites par les parametres
 * @param[in] bignum     Column rank of the bignum, or negative value if there is no big parameter
 * @param[in] options    PipLib options
 */ 
PipQuast* pip_solve(PipMatrix* domain, PipMatrix* parameters,
                    int bignum, PipOptions* options);


/**< Shift solution over -bigparam */
#define SOL_SHIFT  (1 << 0)

/**< Negate solution */
#define SOL_NEGATE (1 << 1)

/**< Remove big parameter */
#define SOL_REMOVE (1 << 2)

/**< Maximum was computed */
#define SOL_MAX    (SOL_SHIFT | SOL_NEGATE)

/**< Create dual leaf */
#define SOL_DUAL   (1 << 3)

/** sol_quast_edit */
PipQuast* sol_quast_edit(int* i, PipQuast* father,
                         int Bg, int Urs_p, int flags);


#if defined(__cplusplus)
  }
#endif

#endif /* define PIPLIB_H */
