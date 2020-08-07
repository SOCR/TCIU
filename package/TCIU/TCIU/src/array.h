#ifndef H_ARRAY_
#define H_ARRAY_

#include <stdio.h>
#include <stdlib.h>

/* ---------- 1D arrays ---------------------- */

#define MAKE_1ARRAY(a,n) do {                                                \
    (a) = malloc((n) * sizeof *(a));                                         \
    if ((a)==NULL) {}                                                        \
} while (0)                                                                  

#define FREE_1ARRAY(a)  do {                                                 \
    free(a);                                                                 \
    a = NULL;                                                                \
} while (0)

/* ---------- 2D arrays ---------------------- */

#define FREE_2ARRAY(a) do {                                                  \
    size_t ARRAY_H2RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H2RESERVED=0; (a)[ARRAY_H2RESERVED]!=NULL; ARRAY_H2RESERVED++)\
        FREE_1ARRAY((a)[ARRAY_H2RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)

/* note: parenthesize first arg because it may be given as `*a' */
#define MAKE_2ARRAY(a,m,n) do {                                              \
    size_t ARRAY_H2RESERVED;                                                 \
    MAKE_1ARRAY(a,(m)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[m] = NULL;                                                           \
    for (ARRAY_H2RESERVED=0; ARRAY_H2RESERVED<(m); ARRAY_H2RESERVED++) {     \
        MAKE_1ARRAY((a)[ARRAY_H2RESERVED],(n));                              \
        if ((a)[ARRAY_H2RESERVED]==NULL) {                                   \
            FREE_2ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

/* ---------- 3D arrays ---------------------- */

#define FREE_3ARRAY(a) do {                                                  \
    size_t ARRAY_H3RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H3RESERVED=0; (a)[ARRAY_H3RESERVED]!=NULL; ARRAY_H3RESERVED++)\
        FREE_2ARRAY((a)[ARRAY_H3RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)

#define MAKE_3ARRAY(a,p,q,r) do {                                            \
    size_t ARRAY_H3RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H3RESERVED=0; ARRAY_H3RESERVED<(p); ARRAY_H3RESERVED++) {     \
        MAKE_2ARRAY((a)[ARRAY_H3RESERVED],(q),(r));                          \
        if ((a)[ARRAY_H3RESERVED]==NULL) {                                   \
            FREE_3ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

/* ---------- 3D arrays ---------------------- */

#define FREE_4ARRAY(a) do {                                                  \
    size_t ARRAY_H4RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H4RESERVED=0; (a)[ARRAY_H4RESERVED]!=NULL; ARRAY_H4RESERVED++)\
        FREE_3ARRAY((a)[ARRAY_H4RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)

#define MAKE_4ARRAY(a,p,q,r,s) do {                                          \
    size_t ARRAY_H4RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H4RESERVED=0; ARRAY_H4RESERVED<(p); ARRAY_H4RESERVED++) {     \
        MAKE_3ARRAY((a)[ARRAY_H4RESERVED],(q),(r),(s));                      \
        if ((a)[ARRAY_H4RESERVED]==NULL) {                                   \
            FREE_4ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

/* ---------- synonyms ---------------------- */

#define MAKE_VECTOR(a,n)    MAKE_1ARRAY(a,n)
#define MAKE_MATRIX(a,m,n)  MAKE_2ARRAY(a,m,n)

#define FREE_VECTOR(a)      FREE_1ARRAY(a)
#define FREE_MATRIX(a)      FREE_2ARRAY(a)

#endif /* H_ARRAY_ */
