/*--------------------------------------------------------------------*/
/*                             T A B . H                              */
/*                                                                    */
/* Copyright Paul Feautrier, 1988, 1993, 1994                         */
/* This file is part of the PIP software                              */
/* PIP is NOT public domain. It can be                                */
/* used freely for research and teaching                              */
/* but NOT for commercial advantage.                                  */
/*--------------------------------------------------------------------*/
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

