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
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 ******************************************************************************/

struct A
    {struct A *precedent;
     char *bout;
    };

struct L
    {int flags;
     Entier d;
     float size;
     union { int unit;
	     Entier * val;
	   } objet;
    };

struct high_water_mark {
    int chunk;
    void * top;
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
#define MAX_DETERMINANT 4

struct T
    {int height, width;
     Entier determinant[MAX_DETERMINANT];
     int l_determinant;
     struct L row[1];
    };

typedef struct T Tableau;

