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
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this library; if not, write to the Free Software Foundation,    *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA         *
 *                                                                            *
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 ******************************************************************************/

#ifndef TAB_H
#define TAB_H

#include <osl/int.h>


#if defined(__cplusplus)
extern "C" 
  {
#endif

struct A {
  struct A* precedent;
  char* bout;
  char* free;
};

struct L {
  int flags;
  osl_int_t d;
  float size;
  union {
    int unit;
    osl_int_t* val;
  } objet;
};

struct high_water_mark {
  int chunk;
  void* top;
};

#define Unit 1
#define Plus 2
#define Minus 4
#define Zero 8
#define Critic 16
#define Unknown 32

#define Sign 62

#define Index(p,i,j) (p)->row[i].objet.val[j]
#define Flag(p,i)    (p)->row[i].flags
#define Denom(p,i)   (p)->row[i].d

struct T {
  int height, width, taille;
  osl_int_t determinant;
  struct L row[1];
};

typedef struct T Tableau;

/* Ced : ajouts specifiques a la PipLib pour funcall. */
Tableau* tab_Matrix2Tableau(PipMatrix*, int, int, int, int, int, int);
Tableau* tab_Matrix2TableauMax(PipMatrix*, int, int, int, int) ;

int dscanf(FILE*, osl_int_p);

#if defined(__cplusplus)
  }
#endif

#endif /* TAB_H */
