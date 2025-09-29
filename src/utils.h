/*
 * Copyright (c) 2025      High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 *
 * Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
 *
 */

#ifndef CLOCKTALK_UTILS_H__
#define CLOCKTALK_UTILS_H__

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define FREE_IF(x)                              \
if(NULL!= x) {                                  \
  free(x);                                      \
  x= NULL;                                      \
}

#define FREE_IF_2D(x)                           \
if(NULL!= x) {                                  \
  FREE_IF(x[0]);                                \
  free(x);                                      \
  x= NULL;                                      \
}

inline static double (**Alloc2d_double2(const long nr,
                                        const long *const nperrow,
                                        const long ntotal))[2]
{
  double (**x)[2]= (double (**)[2]) malloc(sizeof(double (*)[2])* nr);
  x[0]= (double (*)[2]) malloc(sizeof(double[2])* ntotal);
  for(long ir= 1; ir< nr; ++ir) {
    x[ir]= x[ir- 1]+ nperrow[ir- 1];
  }
  return x;
}

#define ALLOC_2D(Type,nr,nperrow,ntotal)                                \
inline static Type **Alloc2d_##Type(const long nr,                      \
                                    const long *const nperrow,          \
                                    const long ntotal)                  \
{                                                                       \
  Type **x= (Type **) malloc(sizeof(Type *)* nr);                       \
  x[0]= (Type *) malloc(sizeof(Type)* ntotal);                          \
  for(long ir= 1; ir< nr; ++ir) {                                       \
    x[ir]= x[ir- 1]+ nperrow[ir- 1];                                    \
  }                                                                     \
  return x;                                                             \
}
ALLOC_2D(double, nr, nperrow, ntotal);
ALLOC_2D(long, nr, nperrow, ntotal);
ALLOC_2D(int, nr, nperrow, ntotal);
/* ALLOC_2D(MsgType,nr,nperrow,ntotal); */

inline static char **Alloc2d(const long nr, const long *const nperrow,
                             const long ntotal, const size_t nb)
{
  char **x= (char **) malloc(sizeof(char *)* nr);
  x[0]= (char *) malloc(nb* ntotal);
  for(long ir= 1; ir< nr; ++ir) {
    x[ir]= x[ir- 1]+ nperrow[ir- 1]* nb;
  }
  return x;
}

#define ErrorIf(cond, ...) do { if(0!= GlOpts.show_opts.error&& (cond)) { printf("*** " __VA_ARGS__); } } while (false)
#define Error(...) do { if(0!= GlOpts.show_opts.error) { printf("*** " __VA_ARGS__); } } while (false)
#define ErrorNL(...) do { if(0!= GlOpts.show_opts.error) { printf("\n*** " __VA_ARGS__); } } while (false)

#define Debug(n, ...) do { if(n<= GlOpts.show_opts.diag) { printf("DEBUG: " __VA_ARGS__); } } while (false)
#define Debug1(...) Debug(1, __VA_ARGS__)
#define Debug2(...) Debug(2, __VA_ARGS__)
#define Debug3(...) Debug(3, __VA_ARGS__)

#define Log(n, ...) do { if(n<= GlOpts.show_opts.diag) { printf("LOG: " __VA_ARGS__); } } while (false)
#define Log1(...) Log(1, __VA_ARGS__)
#define Log2(...) Log(2, __VA_ARGS__)
#define Log3(...) Log(3, __VA_ARGS__)

inline static double Timer_s()
{
  struct timespec ts;
  if(0!= clock_gettime(CLOCK_MONOTONIC, &ts)) {
    Error("Error obtaining clock-value\n");
    return -1.0;
  }
  return ts.tv_sec+ (ts.tv_nsec* 1.0e-9);
}
inline static bool SameTime(const double t0, const double t1) { return fabs(t0- t1)< 0.1; }

#endif  /* CLOCKTALK_UTILS_H__ */
