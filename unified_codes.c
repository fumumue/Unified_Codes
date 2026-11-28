/*
 * unified_codes.c -- unified experimental implementation of algebraic-
 * geometry codes and classical error-correcting codes.
 *
 * The AG engine includes one-point and multi-point code experiments with
 * FR/BMS/VOTE/BASIC decoding paths. The classical engine includes BCH,
 * Golay, Goppa, alternant, Reed-Muller, Hadamard, GRS, Gabidulin, product,
 * and concatenated-code experiments.
 *
 * Build: cc -O2 -std=gnu11 -fopenmp -Wall -Wextra unified_codes.c -lm -o unified_codes
 * OpenMP is optional; omit -fopenmp when it is unavailable.
 * Run:   ./unified_codes --help
 */
#define _POSIX_C_SOURCE 200809L
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint16_t fe;

/* EXTENSION_POINT: algorithmic/memory capacities. Raising m alone is not
 * sufficient; audit these limits against genus, point count and M first. */
#define MAXB 400000      /* max number of basis functions (M+1)          */
#define MAXT 3000        /* max number of errors                         */
#define MAXGEN 100       /* max generators of the ideal E                */
#define MAXGENUS 1310006 /* max genus handled by the voting part         */
#define SGMAX (MAXB + MAXGENUS + 64)
#define MAXTERM 64 /* EXTENSION_POINT: sparse curve term capacity   */
#define MAXEDGE 32 /* EXTENSION_POINT: Newton polygon edge capacity */
#define MAXORB 128 /* EXTENSION_POINT: place/orbit capacity          */

static double now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ================= GF(2^m), supported m <= 12 =================
 * EXTENSION_POINT: increase this bound only after auditing field tables,
 * point enumeration, integer widths, and memory limits. */
static int FM, FS, FO;
static unsigned FPOLY;
static fe *gexp, *glog;

static int gf_init(int m, unsigned poly) {
  FM = m;
  FS = 1 << m;
  FO = FS - 1;
  FPOLY = poly;
  free(gexp);
  free(glog);
  gexp = malloc(sizeof(fe) * (2 * FO + 4));
  glog = malloc(sizeof(fe) * FS);
  unsigned x = 1;
  for (int i = 0; i < FO; ++i) {
    if (i > 0 && x == 1)
      return -1;
    gexp[i] = (fe)x;
    glog[x] = (fe)i;
    x <<= 1;
    if (x & (unsigned)FS)
      x ^= poly;
  }
  if (x != 1)
    return -1;
  for (int i = FO; i < 2 * FO + 4; ++i)
    gexp[i] = gexp[i - FO];
  glog[0] = 0;
  return 0;
}
static inline fe mul(fe a, fe b) {
  return (a && b) ? gexp[glog[a] + glog[b]] : 0;
}
static inline fe inv(fe a) { return gexp[FO - glog[a]]; }
static fe gpow(fe a, long long n) {
  if (n == 0)
    return 1;
  if (!a)
    return 0;
  long long e = n % FO;
  if (e < 0)
    e += FO;
  return gexp[(int)(((long long)glog[a] * e) % FO)];
}

/* ================= curve ================= */
static int T, TX[MAXTERM], TY[MAXTERM];
static fe TC[MAXTERM];
static int dx_deg;
/* Nonzero only for the built-in Hermitian family. */
static int hermitian_q = 0;
static fe coef_at(int px, int py) {
  for (int t = 0; t < T; ++t)
    if (TX[t] == px && TY[t] == py)
      return TC[t];
  return 0;
}
static int gcdi(int a, int b) {
  a = abs(a);
  b = abs(b);
  while (b) {
    int t = a % b;
    a = b;
    b = t;
  }
  return a;
}
static void egcd(int a, int b, int *x, int *y) {
  if (b == 0) {
    *x = (a >= 0) ? 1 : -1;
    *y = 0;
    return;
  }
  int x1, y1;
  egcd(b, a % b, &x1, &y1);
  *x = y1;
  *y = x1 - (a / b) * y1;
}
static long long cross(int ax, int ay, int bx, int by, int cx, int cy) {
  return (long long)(bx - ax) * (cy - ay) - (long long)(by - ay) * (cx - ax);
}
static int hull(int *hx, int *hy) {
  int idx[MAXTERM], n = T;
  for (int i = 0; i < n; i++)
    idx[i] = i;
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      if (TX[idx[j]] < TX[idx[i]] ||
          (TX[idx[j]] == TX[idx[i]] && TY[idx[j]] < TY[idx[i]])) {
        int t = idx[i];
        idx[i] = idx[j];
        idx[j] = t;
      }
  int h[2 * MAXTERM], k = 0;
  for (int i = 0; i < n; i++) {
    while (k >= 2 && cross(TX[h[k - 2]], TY[h[k - 2]], TX[h[k - 1]],
                           TY[h[k - 1]], TX[idx[i]], TY[idx[i]]) <= 0)
      k--;
    h[k++] = idx[i];
  }
  for (int i = n - 2, t = k + 1; i >= 0; i--) {
    while (k >= t && cross(TX[h[k - 2]], TY[h[k - 2]], TX[h[k - 1]],
                           TY[h[k - 1]], TX[idx[i]], TY[idx[i]]) <= 0)
      k--;
    h[k++] = idx[i];
  }
  k--;
  for (int i = 0; i < k; i++) {
    hx[i] = TX[h[i]];
    hy[i] = TY[h[i]];
  }
  return k;
}

/* trial mode: used by the curve normaliser to probe a model without exiting */
static int g_trial = 0, g_trial_fail = 0;

typedef struct {
  int vx, vy, ex, ey, len, wx, wy, a, b, w0, inf;
} Edge;
static Edge ed[MAXEDGE];
static int ne;

static int
build_edges(void) /* returns interior lattice points = genus (nondegenerate) */
{
  int hx[2 * MAXTERM], hy[2 * MAXTERM];
  int h = hull(hx, hy);
  if (h < 3) {
    if (g_trial) {
      g_trial_fail = 1;
      return 0;
    }
    fprintf(stderr, "Newton polygon is degenerate\n");
    exit(1);
  }
  ne = h;
  long long A2 = 0;
  int B = 0;
  for (int k = 0; k < h; k++) {
    int k2 = (k + 1) % h;
    int dx = hx[k2] - hx[k], dy = hy[k2] - hy[k];
    int g = gcdi(dx, dy);
    Edge *e = &ed[k];
    e->vx = hx[k];
    e->vy = hy[k];
    e->len = g;
    e->ex = dx / g;
    e->ey = dy / g;
    e->wx = -e->ey;
    e->wy = e->ex;
    e->inf = (e->wx < 0 || e->wy < 0);
    int x, y;
    egcd(e->wx, e->wy, &x, &y);
    e->a = x;
    e->b = y;
    if (e->wx * e->a + e->wy * e->b != 1) {
      e->a = -e->a;
      e->b = -e->b;
    }
    e->w0 = e->wx * e->vx + e->wy * e->vy;
    A2 += (long long)hx[k] * hy[k2] - (long long)hx[k2] * hy[k];
    B += g;
  }
  if (A2 < 0)
    A2 = -A2;
  return (int)((A2 - B + 2) / 2);
}

/* ================= polynomials over GF(2^m) (small degree) =================
 */
#define PMAXD 320 /* maximum supported Newton-edge polynomial degree */
typedef struct {
  int deg;
  fe c[PMAXD + 1];
} Poly;
static void ptrim(Poly *p) {
  while (p->deg > 0 && p->c[p->deg] == 0)
    p->deg--;
}
static fe peval(const Poly *p, fe x) {
  fe acc = 0;
  for (int i = p->deg; i >= 0; --i)
    acc = (fe)(mul(acc, x) ^ p->c[i]);
  return acc;
}
static void pmod(Poly *a, const Poly *b) /* a mod b, b monic or not */
{
  fe ib = inv(b->c[b->deg]);
  while (a->deg >= b->deg && !(a->deg == 0 && a->c[0] == 0)) {
    fe f = mul(a->c[a->deg], ib);
    int sh = a->deg - b->deg;
    for (int i = 0; i <= b->deg; ++i)
      a->c[i + sh] ^= mul(f, b->c[i]);
    a->c[a->deg] = 0;
    if (a->deg == 0)
      break;
    a->deg--;
    ptrim(a);
    if (a->deg < b->deg)
      break;
  }
}
static Poly pdiv_exact(const Poly *a, const Poly *b) /* assumes b | a */
{
  Poly r = *a, q;
  memset(&q, 0, sizeof q);
  q.deg = a->deg - b->deg;
  fe ib = inv(b->c[b->deg]);
  for (int i = q.deg; i >= 0; --i) {
    fe f = mul(r.c[i + b->deg], ib);
    q.c[i] = f;
    for (int j = 0; j <= b->deg; ++j)
      r.c[i + j] ^= mul(f, b->c[j]);
  }
  return q;
}
static int pdivides(const Poly *b, const Poly *a) {
  Poly r = *a;
  if (a->deg < b->deg)
    return 0;
  pmod(&r, b);
  for (int i = 0; i <= r.deg; ++i)
    if (r.c[i])
      return 0;
  return 1;
}
static Poly pmonic(Poly p) {
  fe i = inv(p.c[p.deg]);
  for (int k = 0; k <= p.deg; ++k)
    p.c[k] = mul(p.c[k], i);
  return p;
}

/* ================= places = Galois orbits of places at infinity
 * ================= */
typedef struct {
  int edge, wx, wy, a, b, w0;
  int d;   /* degree of the orbit (size of the Galois orbit) */
  fe p[8]; /* monic irreducible factor, p[0..d-1] (p[d] = 1) */
  int mult;
  int sigval; /* index of root among the orbit list (for printing) */
} Place;
static Place pl[MAXORB];
static int np;

/* split the edge polynomial of edge k into irreducible factors over GF(2^m) */
static int split_edge(int k, Poly *facs) {
  Edge *e = &ed[k];
  if (e->len > PMAXD) {
    if (g_trial) {
      g_trial_fail = 1;
      return 0;
    }
    fprintf(stderr, "edge %d: lattice length %d exceeds PMAXD = %d\n", k,
            e->len, PMAXD);
    exit(1);
  }
  Poly Q;
  memset(&Q, 0, sizeof Q);
  Q.deg = e->len;
  for (int t = 0; t <= e->len; ++t)
    Q.c[e->len - t] = coef_at(e->vx + t * e->ex, e->vy + t * e->ey);
  ptrim(&Q);
  if (Q.deg != e->len || Q.c[0] == 0) {
    if (g_trial) {
      g_trial_fail = 1;
      return 0;
    }
    fprintf(stderr, "edge polynomial has a zero end coefficient\n");
    exit(1);
  }
  Q = pmonic(Q);
  /* squarefree check: Q' = sum_{j odd} c_j x^{j-1};  gcd(Q,Q') = 1 */
  {
    Poly D;
    memset(&D, 0, sizeof D);
    for (int j = 1; j <= Q.deg; j += 2)
      D.c[j - 1] = Q.c[j];
    D.deg = Q.deg - 1;
    ptrim(&D);
    int allzero = 1;
    for (int i = 0; i <= D.deg; ++i)
      if (D.c[i])
        allzero = 0;
    if (allzero) {
      if (g_trial) {
        g_trial_fail = 1;
        return 0;
      }
      fprintf(stderr, "edge %d: polynomial is a p-th power (degenerate)\n", k);
      exit(1);
    }
    Poly A = Q, B2 = D;
    while (!(B2.deg == 0 && B2.c[0] == 0)) {
      Poly Rr = A;
      pmod(&Rr, &B2);
      A = B2;
      B2 = Rr;
      int z = 1;
      for (int i = 0; i <= B2.deg; ++i)
        if (B2.c[i])
          z = 0;
      if (z)
        break;
    }
    if (A.deg > 0) {
      if (g_trial) {
        g_trial_fail = 1;
        return 0;
      }
      fprintf(stderr, "edge %d: repeated factor (curve is degenerate)\n", k);
      exit(1);
    }
  }
  int nf = 0;
  /* linear factors */
  for (int x = 0; x < FS; ++x)
    if (peval(&Q, (fe)x) == 0) {
      Poly L;
      memset(&L, 0, sizeof L);
      L.deg = 1;
      L.c[1] = 1;
      L.c[0] = (fe)x;
      facs[nf++] = L;
      Poly q = pdiv_exact(&Q, &L);
      Q = q;
    }
  /* quadratic factors */
  if (Q.deg >= 4) {
    for (int a = 0; a < FS && Q.deg >= 4; ++a)
      for (int b = 1; b < FS && Q.deg >= 4; ++b) {
        Poly Qd;
        memset(&Qd, 0, sizeof Qd);
        Qd.deg = 2;
        Qd.c[2] = 1;
        Qd.c[1] = (fe)a;
        Qd.c[0] = (fe)b;
        if (pdivides(&Qd, &Q)) {
          facs[nf++] = Qd;
          Poly q = pdiv_exact(&Q, &Qd);
          Q = q;
        }
      }
  }
  /* EXTENSION_POINT: higher-degree edge factorization is not implemented. */
  if (Q.deg > 5) {
    if (g_trial) {
      g_trial_fail = 1;
      return 0;
    }
    fprintf(stderr, "edge %d: irreducible factor of degree > 5 not supported\n",
            k);
    exit(1);
  }
  if (Q.deg >= 1)
    facs[nf++] = Q;
  return nf;
}

static int ne_inf_orbits(void) {
  np = 0;
  for (int k = 0; k < ne; ++k) {
    if (!ed[k].inf)
      continue;
    Poly *facs = malloc(sizeof(*facs) * (PMAXD + 1));
    if (!facs) {
      fprintf(stderr, "out of memory while factoring edge polynomial\n");
      exit(1);
    }
    int nf = split_edge(k, facs);
    for (int f = 0; f < nf; ++f) {
      if (np >= MAXORB) {
        fprintf(stderr, "too many places\n");
        exit(1);
      }
      Place *P = &pl[np++];
      P->edge = k;
      P->wx = ed[k].wx;
      P->wy = ed[k].wy;
      P->a = ed[k].a;
      P->b = ed[k].b;
      P->w0 = ed[k].w0;
      P->d = facs[f].deg;
      P->mult = 0;
      for (int i = 0; i < 8; ++i)
        P->p[i] = (i <= facs[f].deg) ? facs[f].c[i] : 0;
    }
    free(facs);
  }
  return np;
}

/* ================= ring K_d = GF(2^m)[theta]/(p) and Laurent series over it
 * ================= */
static int KD;
static fe KP[8];
static uint64_t KQD; /* q^d */

static void k_setup(const Place *P) {
  KD = P->d;
  for (int i = 0; i < 8; ++i)
    KP[i] = P->p[i];
  KQD = 1;
  for (int i = 0; i < KD; ++i)
    KQD *= (uint64_t)FS;
}
static inline void k_mul_acc(fe *dst, const fe *a, const fe *b) {
  if (KD == 1) {
    dst[0] ^= mul(a[0], b[0]);
    return;
  }
  fe t[16] = {0};
  for (int i = 0; i < KD; ++i)
    if (a[i])
      for (int j = 0; j < KD; ++j)
        if (b[j])
          t[i + j] ^= mul(a[i], b[j]);
  for (int k = 2 * KD - 2; k >= KD; --k) {
    fe c = t[k];
    if (!c)
      continue;
    for (int l = 0; l < KD; ++l)
      t[k - KD + l] ^= mul(c, KP[l]);
    t[k] = 0;
  }
  for (int i = 0; i < KD; ++i)
    dst[i] ^= t[i];
}
static void k_mul(fe *dst, const fe *a, const fe *b) {
  for (int i = 0; i < KD; ++i)
    dst[i] = 0;
  k_mul_acc(dst, a, b);
}
static int k_zero(const fe *a) {
  for (int i = 0; i < KD; ++i)
    if (a[i])
      return 0;
  return 1;
}
static void k_pow_u64(fe *out, const fe *a, uint64_t n) {
  fe res[8] = {0}, base[8], tmp[8];
  res[0] = 1;
  for (int i = 0; i < KD; ++i)
    base[i] = a[i];
  while (n) {
    if (n & 1) {
      k_mul(tmp, res, base);
      for (int i = 0; i < KD; ++i)
        res[i] = tmp[i];
    }
    n >>= 1;
    if (n) {
      k_mul(tmp, base, base);
      for (int i = 0; i < KD; ++i)
        base[i] = tmp[i];
    }
  }
  for (int i = 0; i < KD; ++i)
    out[i] = res[i];
}
static void k_inv(fe *out, const fe *a) {
  if (KD == 1) {
    out[0] = inv(a[0]);
    return;
  }
  k_pow_u64(out, a, KQD - 2);
}

#define SC(arr, t) ((arr) + (size_t)(t) * KD)

static void s_mul(const fe *a, const fe *b, fe *r, int n) {
  memset(r, 0, (size_t)n * KD * sizeof(fe));
  if (KD == 1) {
    for (int i = 0; i < n; i++) {
      if (!a[i])
        continue;
      int la = glog[a[i]];
      for (int j = 0; i + j < n; j++)
        if (b[j])
          r[i + j] ^= gexp[la + glog[b[j]]];
    }
    return;
  }
  for (int i = 0; i < n; i++) {
    const fe *ai = SC(a, i);
    if (k_zero(ai))
      continue;
    for (int j = 0; i + j < n; j++) {
      const fe *bj = SC(b, j);
      if (k_zero(bj))
        continue;
      k_mul_acc(SC(r, i + j), ai, bj);
    }
  }
}
static void s_inv(const fe *a, fe *r, int n) {
  fe i0[8];
  k_inv(i0, a);
  for (int i = 0; i < KD; ++i)
    SC(r, 0)[i] = i0[i];
  for (int k = 1; k < n; k++) {
    fe s[8] = {0};
    for (int j = 1; j <= k; j++) {
      const fe *aj = SC(a, j), *rk = SC(r, k - j);
      if (k_zero(aj) || k_zero(rk))
        continue;
      k_mul_acc(s, aj, rk);
    }
    k_mul(SC(r, k), s, i0);
  }
}
static void s_pow(const fe *a, long long r, int n, fe *out) {
  size_t sz = (size_t)n * KD;
  fe *base = malloc(sz * sizeof(fe)), *res = calloc(sz, sizeof(fe)),
     *tmp = malloc(sz * sizeof(fe));
  if (r < 0) {
    s_inv(a, base, n);
    r = -r;
  } else
    memcpy(base, a, sz * sizeof(fe));
  res[0] = 1;
  while (r) {
    if (r & 1) {
      s_mul(res, base, tmp, n);
      memcpy(res, tmp, sz * sizeof(fe));
    }
    r >>= 1;
    if (r) {
      s_mul(base, base, tmp, n);
      memcpy(base, tmp, sz * sizeof(fe));
    }
  }
  memcpy(out, res, sz * sizeof(fe));
  free(base);
  free(res);
  free(tmp);
}

static fe *compute_sigma(const Place *P, int N) {
  size_t sz = (size_t)N * KD;
  fe *sg = calloc(sz, sizeof(fe)), *G = malloc(sz * sizeof(fe)),
     *Gs = malloc(sz * sizeof(fe));
  fe *pr = malloc(sz * sizeof(fe)), *iv = malloc(sz * sizeof(fe)),
     *dl = malloc(sz * sizeof(fe));
  if (KD == 1)
    sg[0] = P->p[0];
  else
    sg[1] = 1; /* root of the factor: theta (or -p0 = p0) */
  int iters = 1;
  while ((1 << iters) < N)
    iters++;
  iters += 2;
  for (int it = 0; it <= iters; it++) {
    memset(G, 0, sz * sizeof(fe));
    memset(Gs, 0, sz * sizeof(fe));
    for (int t = 0; t < T; t++) {
      int sh = P->wx * TX[t] + P->wy * TY[t] - P->w0;
      int r = -P->b * TX[t] + P->a * TY[t];
      if (sh >= N)
        continue;
      s_pow(sg, r, N, pr);
      for (int k = 0; k + sh < N; k++)
        for (int c = 0; c < KD; c++)
          SC(G, k + sh)[c] ^= mul(TC[t], SC(pr, k)[c]);
      if (r & 1) {
        s_pow(sg, r - 1, N, pr);
        for (int k = 0; k + sh < N; k++)
          for (int c = 0; c < KD; c++)
            SC(Gs, k + sh)[c] ^= mul(TC[t], SC(pr, k)[c]);
      }
    }
    if (it == iters) {
      for (size_t k = 0; k < sz; k++)
        if (G[k]) {
          fprintf(stderr, "Newton iteration did not converge\n");
          exit(1);
        }
      break;
    }
    if (k_zero(Gs)) {
      fprintf(stderr, "multiple root of the edge polynomial\n");
      exit(1);
    }
    s_inv(Gs, iv, N);
    s_mul(G, iv, dl, N);
    for (size_t k = 0; k < sz; k++)
      sg[k] ^= dl[k];
  }
  free(G);
  free(Gs);
  free(pr);
  free(iv);
  free(dl);
  return sg;
}

typedef struct {
  int vmin, N, rmin, rmax;
  fe *pw;
} PSer;

static void place_series(const Place *P, int ehi, int J, PSer *S) {
  k_setup(P);
  int c[4] = {0, P->wx * (dx_deg - 1), P->wy * J,
              P->wx * (dx_deg - 1) + P->wy * J};
  int vmin = c[0];
  for (int i = 1; i < 4; i++)
    if (c[i] < vmin)
      vmin = c[i];
  int rc[4] = {0, -P->b * (dx_deg - 1), P->a * J,
               -P->b * (dx_deg - 1) + P->a * J};
  int rmin = rc[0], rmax = rc[0];
  for (int i = 1; i < 4; i++) {
    if (rc[i] < rmin)
      rmin = rc[i];
    if (rc[i] > rmax)
      rmax = rc[i];
  }
  S->vmin = vmin;
  S->rmin = rmin;
  S->rmax = rmax;
  int N = ehi - vmin + 1;
  S->N = N > 0 ? N : 0;
  S->pw = NULL;
  if (N <= 0)
    return;
  size_t sz = (size_t)N * KD;
  fe *sg = compute_sigma(P, N), *si = malloc(sz * sizeof(fe));
  s_inv(sg, si, N);
  int R = rmax - rmin + 1;
  S->pw = calloc((size_t)R * sz, sizeof(fe));
  fe *z = S->pw + (size_t)(0 - rmin) * sz;
  z[0] = 1;
  for (int r = 1; r <= rmax; r++)
    s_mul(S->pw + (size_t)(r - 1 - rmin) * sz, sg,
          S->pw + (size_t)(r - rmin) * sz, N);
  for (int r = -1; r >= rmin; r--)
    s_mul(S->pw + (size_t)(r + 1 - rmin) * sz, si,
          S->pw + (size_t)(r - rmin) * sz, N);
  free(sg);
  free(si);
}

/* ================= L(D): D supported on places at infinity =================
 */
typedef struct {
  int dim, K, J;
  fe *B;
} LSp;

static int genus_g;
/* Experimental multi-place Newton ordering: 0=legacy column order, 1=lexicographic
 * valuation-vector order. This changes elimination pivot order, not the divisor
 * constraints or the resulting function space. */
static int multi_newton_order = 0;

static int monomial_order_less(int ca, int cb) {
  int ia = ca % dx_deg, ja = ca / dx_deg;
  int ib = cb % dx_deg, jb = cb / dx_deg;
  for (int q = 0; q < np; ++q) {
    long va = (long)pl[q].wx * ia + (long)pl[q].wy * ja;
    long vb = (long)pl[q].wx * ib + (long)pl[q].wy * jb;
    if (va != vb)
      return va < vb;
  }
  if (ja != jb) return ja < jb;
  return ia < ib;
}

/* 旧 solve_L() の名前を変更 */
static int solve_L_laurent(const int *mult, int Jfix, int Jstep, int Jmax,
                           int verbose, LSp *out) {
  int degD = 0, maxm = 0;
  for (int q = 0; q < np; q++) {
    degD += mult[q] * pl[q].d;
    if (mult[q] > maxm)
      maxm = mult[q];
  }
  int certifiable = degD >= 2 * genus_g - 1, expected = degD - genus_g + 1;
  int J = (Jfix >= 0) ? Jfix : (maxm / 4 + 16);
  if (J < 8)
    J = 8;
  fe *Bv = NULL;
  int dim = 0, K = 0, lastdim = -1, stable = 0;
  for (;;) {
    K = dx_deg * (J + 1);
    int ehi_q[MAXORB];
    PSer S[MAXORB];
    size_t E = 0;
    for (int q = 0; q < np; q++) {
      ehi_q[q] = -mult[q] - 1;
      place_series(&pl[q], ehi_q[q], J, &S[q]);
      E += (size_t)S[q].N * pl[q].d;
    }
    if ((double)E * K > 9e8) {
      fprintf(stderr, "matrix too large (E=%zu, K=%d)\n", E, K);
      exit(1);
    }
    fe *Mx = calloc((size_t)((E > 0 && K > 0) ? E * K : 1), sizeof(fe));
    int *col_order = malloc(sizeof(int) * (size_t)K);
    int *rank_of_col = malloc(sizeof(int) * (size_t)K);
    if (!col_order || !rank_of_col) {
      free(col_order); free(rank_of_col); free(Mx);
      for (int q = 0; q < np; ++q) free(S[q].pw);
      free(Bv);
      return -1;
    }
    for (int c = 0; c < K; ++c) col_order[c] = c;
    if (multi_newton_order) {
      /* Stable insertion sort by the vector of valuations at all places. */
      for (int c = 1; c < K; ++c) {
        int x = col_order[c], j = c;
        while (j > 0 && monomial_order_less(x, col_order[j - 1])) {
          col_order[j] = col_order[j - 1];
          --j;
        }
        col_order[j] = x;
      }
    }
    for (int c = 0; c < K; ++c) rank_of_col[col_order[c]] = c;
    size_t off = 0;
    for (int q = 0; q < np; q++) {
      int d = pl[q].d;
      if (S[q].N > 0) {
        size_t sz = (size_t)S[q].N * d;
        for (int j = 0; j <= J; j++)
          for (int i = 0; i < dx_deg; i++) {
            int col = j * dx_deg + i;
            int rank = rank_of_col[col];
            int v = pl[q].wx * i + pl[q].wy * j, r = -pl[q].b * i + pl[q].a * j;
            const fe *pw = S[q].pw + (size_t)(r - S[q].rmin) * sz;
            for (int t = 0; v + t <= ehi_q[q]; t++)
              for (int k = 0; k < d; k++) {
                size_t row = off + (size_t)(v + t - S[q].vmin) * d + k;
                Mx[row * K + rank] = pw[(size_t)t * d + k];
              }
          }
      }
      off += (size_t)S[q].N * d;
      free(S[q].pw);
    }
    int *pc = malloc(sizeof(int) * (K + 1));
    size_t rk = 0;
    int *lr = malloc(sizeof(int) * K);
    fe *tmprow = malloc(sizeof(fe) * K);
    for (int c = 0; c < K && rk < E; c++) {
      size_t p = rk;
      while (p < E && !Mx[p * K + c])
        p++;
      if (p == E)
        continue;
      if (p != rk) {
        memcpy(tmprow, Mx + p * K, K * sizeof(fe));
        memcpy(Mx + p * K, Mx + rk * K, K * sizeof(fe));
        memcpy(Mx + rk * K, tmprow, K * sizeof(fe));
      }
      fe *R = Mx + rk * K;
      fe iv = inv(R[c]);
      for (int col = c; col < K; col++)
        R[col] = mul(R[col], iv);
      for (int col = c; col < K; col++)
        lr[col] = R[col] ? glog[R[col]] : -1;
#pragma omp parallel for schedule(static)
      for (long long ii = (long long)rk + 1; ii < (long long)E; ii++) {
        fe *Dr = Mx + (size_t)ii * K;
        fe f = Dr[c];
        if (!f)
          continue;
        int lf = glog[f];
        for (int col = c; col < K; col++) {
          int x = lr[col];
          if (x >= 0)
            Dr[col] ^= gexp[lf + x];
        }
      }
      pc[rk++] = c;
    }
    char *isp = calloc(K, 1);
    for (size_t i = 0; i < rk; i++)
      isp[pc[i]] = 1;
    dim = K - (int)rk;
    free(Bv);
    Bv = calloc((size_t)(dim ? dim : 1) * K, sizeof(fe));
    int bi = 0;
    for (int f = 0; f < K; f++) {
      if (isp[f])
        continue;
      fe *x = Bv + (size_t)bi * K;
      x[f] = 1;
      for (long long i = (long long)rk - 1; i >= 0; i--) {
        int c = pc[i];
        const fe *R = Mx + (size_t)i * K;
        fe s = 0;
        for (int col = c + 1; col < K; col++)
          if (R[col] && x[col])
            s ^= mul(R[col], x[col]);
        x[c] = s;
      }
      /* Return coefficient vectors in the original j*dx_deg+i layout. */
      for (int c = 0; c < K; ++c) {
        int original_col = col_order[c];
        if (original_col != c) {
          /* x is temporarily kept in pivot-order coordinates; scatter below. */
        }
      }
      fe *scattered = calloc((size_t)K, sizeof(fe));
      if (!scattered) {
        free(col_order); free(rank_of_col); free(Mx); free(pc); free(lr);
        free(tmprow); free(isp); free(Bv);
        return -1;
      }
      for (int c = 0; c < K; ++c) scattered[col_order[c]] = x[c];
      memcpy(x, scattered, (size_t)K * sizeof(fe));
      free(scattered);
      bi++;
    }
    free(col_order);
    free(rank_of_col);
    free(Mx);
    free(pc);
    free(lr);
    free(tmprow);
    free(isp);
    if (verbose) {
      printf("#   J=%d: unknowns %d, equations %zu, dim L(D) = %d", J, K, E,
             dim);
      if (certifiable)
        printf("  (Riemann-Roch: %d)", expected);
      printf("\n");
      fflush(stdout);
    }
    if (Jfix >= 0)
      break;
    if (certifiable) {
      if (dim == expected)
        break;
      if (dim > expected) {
        fprintf(stderr,
                "dimension exceeds Riemann-Roch: curve assumptions violated\n");
        break;
      }
    } else {
      if (dim == lastdim) {
        stable++;
        if (stable >= 2)
          break;
      } else
        stable = 0;
      lastdim = dim;
    }
    if (J >= Jmax) {
      fprintf(stderr, "warning: Jmax reached, basis may be incomplete\n");
      break;
    }
    int deficit = certifiable ? expected - dim : 0;
    int st = Jstep;
    if (certifiable && deficit / 6 + 1 > st)
      st = deficit / 6 + 1;
    J += st;
    if (J > Jmax)
      J = Jmax;
  }
  out->dim = dim;
  out->K = K;
  out->J = J;
  out->B = Bv;
  if (certifiable && dim != expected)
    return -1;
  return 0;
}

/* ================= points ================= */
static int npts;
static fe *PXv, *PYv;

static void enumerate_points(void) {
  int cap = 1024;
  npts = 0;
  PXv = malloc(sizeof(fe) * cap);
  PYv = malloc(sizeof(fe) * cap);
  fe *ypw = malloc(sizeof(fe) * (size_t)T * FS);
  for (int t = 0; t < T; ++t)
    for (int y = 0; y < FS; ++y)
      ypw[(size_t)t * FS + y] = gpow((fe)y, TY[t]);
  for (int x = 0; x < FS; ++x) {
    fe xp[MAXTERM];
    for (int t = 0; t < T; ++t)
      xp[t] = mul(TC[t], gpow((fe)x, TX[t]));
    for (int y = 0; y < FS; ++y) {
      fe acc = 0;
      for (int t = 0; t < T; ++t)
        acc ^= mul(xp[t], ypw[(size_t)t * FS + y]);
      if (!acc) {
        if (npts == cap) {
          cap *= 2;
          PXv = realloc(PXv, sizeof(fe) * cap);
          PYv = realloc(PYv, sizeof(fe) * cap);
        }
        PXv[npts] = (fe)x;
        PYv[npts] = (fe)y;
        npts++;
      }
    }
  }
  free(ypw);
}

/* evaluation matrix: out[d*npts + p] = (d-th basis vector of L)(P_p) */
static fe *eval_space(const LSp *L) {
  int dim = L->dim, K = L->K, J = L->J;
  fe *out = malloc(sizeof(fe) * (size_t)(dim ? dim : 1) * npts);
  if (!out)
    return NULL;

  /* The Laurent/Newton basis vectors are very sparse.  Store only their
     nonzero monomial coefficients instead of scanning dim*K entries at
     every evaluation point. */
  size_t nnz = 0;
  for (size_t z = 0; z < (size_t)dim * K; ++z)
    nnz += (L->B[z] != 0);
  int *rowptr = malloc(sizeof(int) * (size_t)(dim + 1));
  int *sidx = malloc(sizeof(int) * (nnz ? nnz : 1));
  uint16_t *slog = malloc(sizeof(uint16_t) * (nnz ? nnz : 1));
  if (!rowptr || !sidx || !slog) {
    free(out);
    free(rowptr);
    free(sidx);
    free(slog);
    return NULL;
  }
  size_t nz = 0;
  rowptr[0] = 0;
  for (int d = 0; d < dim; ++d) {
    const fe *b = L->B + (size_t)d * K;
    for (int c = 0; c < K; ++c) {
      if (b[c]) {
        sidx[nz] = c;
        slog[nz] = glog[b[c]];
        ++nz;
      }
    }
    rowptr[d + 1] = (int)nz;
  }

#pragma omp parallel
  {
    fe *xp = malloc(sizeof(fe) * (size_t)dx_deg);
    fe *yp = malloc(sizeof(fe) * (size_t)(J + 1));
    int *ml = malloc(sizeof(int) * (size_t)K);
#pragma omp for schedule(static)
    for (int p = 0; p < npts; ++p) {
      xp[0] = 1;
      for (int i = 1; i < dx_deg; ++i)
        xp[i] = mul(xp[i - 1], PXv[p]);
      yp[0] = 1;
      for (int j = 1; j <= J; ++j)
        yp[j] = mul(yp[j - 1], PYv[p]);
      for (int c = 0; c < K; ++c) {
        fe v = mul(xp[c % dx_deg], yp[c / dx_deg]);
        ml[c] = v ? (int)glog[v] : -1;
      }
      for (int d = 0; d < dim; ++d) {
        fe acc = 0;
        for (int z = rowptr[d]; z < rowptr[d + 1]; ++z) {
          int lm = ml[sidx[z]];
          if (lm >= 0)
            acc ^= gexp[(int)slog[z] + lm];
        }
        out[(size_t)d * npts + p] = acc;
      }
    }
    free(xp);
    free(yp);
    free(ml);
  }
  free(rowptr);
  free(sidx);
  free(slog);
  return out;
}

static fe dotv(const fe *a, const fe *b, int n) {
  fe s = 0;
  for (int k = 0; k < n; ++k)
    if (a[k] && b[k])
      s ^= gexp[glog[a[k]] + glog[b[k]]];
  return s;
}

static int nullspace_gf(fe *A, int E, int K, fe **out) {
  int *pc = malloc(sizeof(int) * (K + 1));
  int rk = 0;
  int *lr = malloc(sizeof(int) * K);
  fe *tmp = malloc(sizeof(fe) * K);
  for (int c = 0; c < K && rk < E; ++c) {
    int p = rk;
    while (p < E && !A[(size_t)p * K + c])
      p++;
    if (p == E)
      continue;
    if (p != rk) {
      memcpy(tmp, A + (size_t)p * K, K * sizeof(fe));
      memcpy(A + (size_t)p * K, A + (size_t)rk * K, K * sizeof(fe));
      memcpy(A + (size_t)rk * K, tmp, K * sizeof(fe));
    }
    fe *R = A + (size_t)rk * K;
    fe iv = inv(R[c]);
    for (int col = c; col < K; ++col)
      R[col] = mul(R[col], iv);
    for (int col = c; col < K; ++col)
      lr[col] = R[col] ? glog[R[col]] : -1;
#pragma omp parallel for schedule(static)
    for (int ii = rk + 1; ii < E; ++ii) {
      fe *Dr = A + (size_t)ii * K;
      fe f = Dr[c];
      if (!f)
        continue;
      int lf = glog[f];
      for (int col = c; col < K; ++col) {
        int x = lr[col];
        if (x >= 0)
          Dr[col] ^= gexp[lf + x];
      }
    }
    pc[rk++] = c;
  }
  char *isp = calloc(K, 1);
  for (int i = 0; i < rk; ++i)
    isp[pc[i]] = 1;
  int dim = K - rk;
  fe *B = calloc((size_t)(dim ? dim : 1) * K, sizeof(fe));
  int bi = 0;
  for (int f = 0; f < K; ++f) {
    if (isp[f])
      continue;
    fe *x = B + (size_t)bi * K;
    x[f] = 1;
    for (int i = rk - 1; i >= 0; --i) {
      int c = pc[i];
      const fe *R = A + (size_t)i * K;
      fe s = 0;
      for (int col = c + 1; col < K; ++col)
        if (R[col] && x[col])
          s ^= mul(R[col], x[col]);
      x[c] = s;
    }
    bi++;
  }
  free(pc);
  free(lr);
  free(tmp);
  free(isp);
  *out = B;
  return dim;
}

/* right inverse: R (n x r) with U R = I_r, U is r x n of full row rank */
static int right_inverse(const fe *U, int r, int n, fe *R) {
  int W = n + r;
  fe *Aug = calloc((size_t)r * W, sizeof(fe));
  for (int i = 0; i < r; ++i) {
    memcpy(Aug + (size_t)i * W, U + (size_t)i * n, n * sizeof(fe));
    Aug[(size_t)i * W + n + i] = 1;
  }
  int *P = malloc(sizeof(int) * r);
  int rk = 0;
  int *lr = malloc(sizeof(int) * W);
  for (int c = 0; c < n && rk < r; ++c) {
    int p = rk;
    while (p < r && !Aug[(size_t)p * W + c])
      p++;
    if (p == r)
      continue;
    if (p != rk)
      for (int k = 0; k < W; ++k) {
        fe t = Aug[(size_t)p * W + k];
        Aug[(size_t)p * W + k] = Aug[(size_t)rk * W + k];
        Aug[(size_t)rk * W + k] = t;
      }
    fe *Rw = Aug + (size_t)rk * W;
    fe iv = inv(Rw[c]);
    for (int k = c; k < W; ++k)
      Rw[k] = mul(Rw[k], iv);
    for (int k = c; k < W; ++k)
      lr[k] = Rw[k] ? glog[Rw[k]] : -1;
#pragma omp parallel for schedule(static)
    for (int ii = 0; ii < r; ++ii) {
      if (ii == rk)
        continue;
      fe *Dr = Aug + (size_t)ii * W;
      fe f = Dr[c];
      if (!f)
        continue;
      int lf = glog[f];
      for (int k = c; k < W; ++k) {
        int x = lr[k];
        if (x >= 0)
          Dr[k] ^= gexp[lf + x];
      }
    }
    P[rk++] = c;
  }
  if (rk < r) {
    free(Aug);
    free(P);
    free(lr);
    return -1;
  }
  memset(R, 0, sizeof(fe) * (size_t)n * r);
  for (int a = 0; a < r; ++a)
    for (int b = 0; b < r; ++b)
      R[(size_t)P[a] * r + b] = Aug[(size_t)a * W + n + b];
  free(Aug);
  free(P);
  free(lr);
  return 0;
}

/* ================= one-point echelon basis at a rational place
 * ================= */
/* Turns the nullspace basis of L(mP) into functions f_s with pole order exactly
 * s and Laurent coefficient 1 at u^-s ("monic").  polev[d] = pole order of
 * basis vector d (0 for the constant).  Returns the number of constants (should
 * be 1).              */
static int semigroup_echelon(LSp *L, const Place *P, int m, int *polev) {
  int dim = L->dim, K = L->K, J = L->J;
  PSer S;
  place_series(P, -1, J, &S);
  fe *phi = calloc((size_t)dim * m, sizeof(fe));
#pragma omp parallel for schedule(static)
  for (int d = 0; d < dim; ++d) {
    const fe *b = L->B + (size_t)d * K;
    for (int col = 0; col < K; ++col) {
      fe x = b[col];
      if (!x)
        continue;
      int j = col / dx_deg, i = col % dx_deg;
      int v = P->wx * i + P->wy * j, r = -P->b * i + P->a * j;
      const fe *pw = S.pw + (size_t)(r - S.rmin) * S.N;
      for (int s = 1; s <= m; ++s) {
        int t = -s - v;
        if (t < 0)
          continue;
        if (t >= S.N)
          break;
        fe cf = pw[t];
        if (cf)
          phi[(size_t)d * m + s - 1] ^= mul(cf, x);
      }
    }
  }
  free(S.pw);
  char *used = calloc(dim, 1);
  for (int s = m; s >= 1; --s) {
    int p = -1;
    for (int d = 0; d < dim; ++d)
      if (!used[d] && phi[(size_t)d * m + s - 1]) {
        p = d;
        break;
      }
    if (p < 0)
      continue;
    fe iv = inv(phi[(size_t)p * m + s - 1]);
    for (int t = 0; t < m; ++t)
      phi[(size_t)p * m + t] = mul(phi[(size_t)p * m + t], iv);
    for (int c = 0; c < K; ++c)
      L->B[(size_t)p * K + c] = mul(L->B[(size_t)p * K + c], iv);
#pragma omp parallel for schedule(static)
    for (int d = 0; d < dim; ++d) {
      if (d == p || used[d])
        continue;
      fe f = phi[(size_t)d * m + s - 1];
      if (!f)
        continue;
      for (int t = 0; t < m; ++t)
        phi[(size_t)d * m + t] ^= mul(f, phi[(size_t)p * m + t]);
      for (int c = 0; c < K; ++c)
        L->B[(size_t)d * K + c] ^= mul(f, L->B[(size_t)p * K + c]);
    }
    used[p] = 1;
    polev[p] = s;
  }
  int nconst = 0;
  for (int d = 0; d < dim; ++d)
    if (!used[d]) {
      polev[d] = 0;
      nconst++;
    }
  free(phi);
  free(used);
  return nconst;
}

/* ================================================================== */
/*  One-point basis from the NEWTON POLYGON  (replaces solve_L +        */
/*  semigroup_echelon whenever it applies)                              */
/* ================================================================== */
/* Let P be the place of an edge e of Newton(F) of lattice length 1, with
 * primitive inward normal w_e.  A Laurent monomial x^i y^j has valuation
 * w_k . (i,j) at the places of edge k, so
 *      x^i y^j in L(M P)   <=>   w_k.(i,j) >= 0 (all k != e),  -w_e.(i,j) <= M,
 * its pole order at P being -w_e.(i,j).  Monomials of equal weight differ only
 * by lower pole orders (modulo F), so one representative per weight suffices.
 * If the weights found are exactly the M+1-g many pole orders that
 * Riemann-Roch demands (M >= 2g-1), they form a complete basis of L(M P):
 * the result is certified at run time, otherwise the caller falls back to
 * the Laurent-expansion method.                                              */
typedef struct {
  int n, laurent;
  int *rho, *ei, *ej; /* sorted by pole order rho */
} NBas;
static int use_newton = 1;
#define NCAND_MAX (MAXORB + MAXEDGE)

static long nb_fdiv(long a, long b) { /* floor(a/b), b > 0 */
  long q = a / b;
  if ((a % b) && (a < 0))
    --q;
  return q;
}

static int newton_basis(int e, int Mx, NBas *out) {
  if (ed[e].len != 1)
    return -1; /* a single place exists only on a primitive edge */
  long ca[MAXEDGE + 1], cb[MAXEDGE + 1], ck[MAXEDGE + 1];
  int nc = 0;
  for (int k = 0; k < ne; ++k)
    if (k != e) {
      ca[nc] = ed[k].wx;
      cb[nc] = ed[k].wy;
      ck[nc] = 0;
      nc++;
    }
  ca[nc] = ed[e].wx; /* pole order -w_e.(i,j) <= M  <=>  w_e.(i,j) + M >= 0 */
  cb[nc] = ed[e].wy;
  ck[nc] = Mx;
  nc++;
  int *ri = malloc(sizeof(int) * (Mx + 1)),
      *rj = malloc(sizeof(int) * (Mx + 1));
  char *have = calloc(Mx + 1, 1);
  int bad = 0;
  for (int dir = 0; dir < 2 && !bad; ++dir)
    for (int i = dir ? -1 : 0;; i += dir ? -1 : 1) {
      if (i > (1 << 20) || i < -(1 << 20)) {
        bad = 1;
        break;
      }
      long jl = -(1L << 40), jh = (1L << 40);
      double rl = -1e18, rh = 1e18;
      int empty = 0;
      for (int c = 0; c < nc && !empty; ++c) {
        long v = ca[c] * i + ck[c];
        if (cb[c] == 0) {
          if (v < 0)
            empty = 1;
        } else if (cb[c] > 0) {
          long l = -nb_fdiv(v, cb[c]);
          if (l > jl)
            jl = l;
          double r = -(double)v / (double)cb[c];
          if (r > rl)
            rl = r;
        } else {
          long h = nb_fdiv(v, -cb[c]);
          if (h < jh)
            jh = h;
          double r = (double)v / (double)(-cb[c]);
          if (r < rh)
            rh = r;
        }
      }
      if (empty || rl > rh + 1e-9)
        break; /* convex region: slices are contiguous in i */
      for (long j = jl; j <= jh; ++j) {
        int w = -(int)(ed[e].wx * (long)i + ed[e].wy * j);
        if (w < 0 || w > Mx) {
          bad = 1;
          break;
        }
        int neg = (i < 0 || j < 0);
        if (!have[w]) {
          have[w] = 1;
          ri[w] = i;
          rj[w] = (int)j;
        } else if (neg == 0 && (ri[w] < 0 || rj[w] < 0)) {
          ri[w] = i; /* prefer a representative without negative exponents */
          rj[w] = (int)j;
        }
      }
      if (bad)
        break;
    }
  int cnt = 0;
  for (int w = 0; w <= Mx; ++w)
    cnt += have[w];
  int ok =
      !bad && have[0] && ri[0] == 0 && rj[0] == 0 && cnt == Mx + 1 - genus_g;
  int laurent = 0;
  if (ok)
    for (int w = 0; w <= Mx; ++w)
      if (have[w] && (ri[w] < 0 || rj[w] < 0))
        laurent = 1;
  if (ok && laurent) {
    /* negative exponents need: the only point of the model on x=0 or y=0 is
     * the origin, i.e. F(0,y) and F(x,0) are single monomials */
    int tx0 = 0, ty0 = 0;
    for (int t = 0; t < T; ++t) {
      tx0 += (TX[t] == 0);
      ty0 += (TY[t] == 0);
    }
    if (tx0 > 1 || ty0 > 1)
      ok = 0;
  }
  if (ok) {
    out->n = cnt;
    out->laurent = laurent;
    out->rho = malloc(sizeof(int) * cnt);
    out->ei = malloc(sizeof(int) * cnt);
    out->ej = malloc(sizeof(int) * cnt);
    int k = 0;
    for (int w = 0; w <= Mx; ++w)
      if (have[w]) {
        out->rho[k] = w;
        out->ei[k] = ri[w];
        out->ej[k] = rj[w];
        k++;
      }
  }
  free(ri);
  free(rj);
  free(have);
  return ok ? 0 : -1;
}

/*
 * Newton-first L(D) solver.
 *
 * Supported:
 *   D = m P, where P is a rational place on a primitive Newton edge.
 *   The Newton basis must be non-Laurent and fit the existing
 *   rectangular monomial representation used by LSp/eval_space.
 *
 * Unsupported divisors are delegated to solve_L_laurent().
 */
static int solve_L(const int *mult, int Jfix, int Jstep, int Jmax, int verbose,
                   LSp *out) {
  if (!mult || !out)
    return -1;

  memset(out, 0, sizeof(*out));

  int qsel = -1;
  int degD = 0;

  for (int q = 0; q < np; ++q) {
    if (mult[q] < 0)
      return -1;

    degD += mult[q] * pl[q].d;

    if (mult[q] != 0) {
      if (qsel >= 0) {
        qsel = -2; /* multiple places */
        break;
      }
      qsel = q;
    }
  }

  /*
   * D = 0: L(0) = <1>.
   */
  if (qsel == -1) {
    out->dim = 1;
    out->J = 0;
    out->K = dx_deg;
    out->B = calloc((size_t)dx_deg, sizeof(fe));

    if (!out->B)
      return -1;

    out->B[0] = 1;
    return 0;
  }

  /*
   * Newton one-point basis does not directly handle
   * multiple places or non-rational places.
   */
  if (qsel < 0 || pl[qsel].d != 1 || pl[qsel].edge < 0 || pl[qsel].edge >= ne ||
      ed[pl[qsel].edge].len != 1) {
    return solve_L_laurent(mult, Jfix, Jstep, Jmax, verbose, out);
  }

  /*
   * Riemann-Roch certification is required here.
   */
  if (degD < 2 * genus_g - 1) {
    return solve_L_laurent(mult, Jfix, Jstep, Jmax, verbose, out);
  }

  NBas NB = {0};

  if (newton_basis(pl[qsel].edge, mult[qsel], &NB) != 0) {
    free(NB.rho);
    free(NB.ei);
    free(NB.ej);

    return solve_L_laurent(mult, Jfix, Jstep, Jmax, verbose, out);
  }

  int expected = degD - genus_g + 1;

  /*
   * Require the Riemann-Roch dimension and a basis that
   * can be represented by the existing LSp coefficient array.
   */
  int supported = (NB.n == expected && !NB.laurent);
  int maxj = 0;

  for (int k = 0; k < NB.n && supported; ++k) {
    if (NB.ei[k] < 0 || NB.ei[k] >= dx_deg || NB.ej[k] < 0) {
      supported = 0;
      break;
    }

    if (NB.ej[k] > maxj)
      maxj = NB.ej[k];
  }

  if (!supported) {
    free(NB.rho);
    free(NB.ei);
    free(NB.ej);

    return solve_L_laurent(mult, Jfix, Jstep, Jmax, verbose, out);
  }

  int J = maxj;
  int K = dx_deg * (J + 1);

  if (K <= 0 || (size_t)NB.n > SIZE_MAX / (size_t)K / sizeof(fe)) {
    free(NB.rho);
    free(NB.ei);
    free(NB.ej);
    return -1;
  }

  fe *B = calloc((size_t)NB.n * K, sizeof(fe));

  if (!B) {
    free(NB.rho);
    free(NB.ei);
    free(NB.ej);
    return -1;
  }

  /*
   * Convert each Newton monomial x^i y^j into
   * the existing LSp coefficient layout:
   *
   *   column = j * dx_deg + i
   */
  for (int k = 0; k < NB.n; ++k) {
    int col = NB.ej[k] * dx_deg + NB.ei[k];
    B[(size_t)k * K + col] = 1;
  }

  out->dim = NB.n;
  out->K = K;
  out->J = J;
  out->B = B;

  if (verbose) {
    printf("# Newton L(D): dim = %d "
           "(Riemann-Roch: %d), J = %d\n",
           NB.n, expected, J);
  }

  free(NB.rho);
  free(NB.ei);
  free(NB.ej);

  return 0;
}

/* value of x^i y^j at a point (i, j may be negative; then x, y != 0) */
static fe mono_eval(fe x, fe y, int i, int j) {
  if (i == 0 && j == 0)
    return 1;
  if ((i != 0 && !x) || (j != 0 && !y))
    return 0;
  long long ex =
      ((long long)i * (x ? glog[x] : 0) + (long long)j * (y ? glog[y] : 0)) %
      FO;
  if (ex < 0)
    ex += FO;
  return gexp[ex];
}

/* Fast row evaluation: out[k] = x_k^i y_k^j for k < n, using precomputed
 * log arrays (lx[k], ly[k] = -1 for a zero coordinate).  Same result as
 * mono_eval, but without per-element 64-bit division (Barrett reduction). */
static void mono_eval_row(fe *out, int n, const int *lx, const int *ly, int i,
                          int j) {
  if (i == 0 && j == 0) {
    for (int k = 0; k < n; ++k)
      out[k] = 1;
    return;
  }
  int ie = i % FO, je = j % FO;
  if (ie < 0)
    ie += FO;
  if (je < 0)
    je += FO;
  const uint64_t magic = ((uint64_t)1 << 40) / (uint64_t)FO + 1;
  for (int k = 0; k < n; ++k) {
    int a = lx[k], b = ly[k];
    if ((i != 0 && a < 0) || (j != 0 && b < 0)) {
      out[k] = 0;
      continue;
    }
    uint32_t ex = (uint32_t)(ie * (a < 0 ? 0 : a) + je * (b < 0 ? 0 : b));
    uint32_t q = (uint32_t)(((uint64_t)ex * magic) >> 40);
    out[k] = gexp[ex - q * (uint32_t)FO];
  }
}

/* d_FR(M) = min { nu(s) : s > M, s in H },  nu(s) = #{(a,b) in H^2, a+b = s} */
static int dfr_of(const int *inH_arr, int M, int span) {
  int best = 1 << 30;
  for (int s = M + 1; s <= M + span; ++s) {
    int c = 0;
    for (int a = 0; a <= s; ++a)
      if (inH_arr[a] && inH_arr[s - a])
        c++;
    if (c < best)
      best = c;
  }
  return best;
}

/* ================================================================== */
/*  One-point decoders (generalised from the Suzuki version)           */
/* ================================================================== */
static int GENUS, M_SYND, M_EXT, TAU, TAU_V, L_DEG, ROW_MAX, ASHIFT, NUM_ERRORS;
static int show_locator_equations = 0;

typedef struct {
  int rho;
  fe *val;
} Basis;
static Basis basis[MAXB];
static int nbasis;
static int in_H[SGMAX], idx_of_ord[SGMAX];
static int genH[MAXB], ngenH;

static int count_basis_upto(int rho) {
  int n = 0;
  for (int i = 0; i < nbasis; ++i)
    if (basis[i].rho <= rho)
      ++n;
  return n;
}

static void init_semigroup(void) {
  for (int x = 0; x < SGMAX; ++x) {
    in_H[x] = (x > M_EXT);
    idx_of_ord[x] = -1;
  }
  for (int i = 0; i < nbasis; ++i) {
    in_H[basis[i].rho] = 1;
    idx_of_ord[basis[i].rho] = i;
  }
  ngenH = 0;
  for (int s = 1; s <= M_EXT; ++s) {
    if (!in_H[s])
      continue;
    int gen = 1;
    for (int a = 1; a <= s / 2 && gen; ++a)
      if (in_H[a] && in_H[s - a])
        gen = 0;
    if (gen)
      genH[ngenH++] = s;
  }
}

/* ---- syndromes ---- */
static int use_poly;
static int compact_fr_mode;
static fe basis_value_at(int basis_index, int point_index);
static fe basis_syndrome(int basis_index);
static const fe *recv;
static int *nz_idx;
static int nz_cnt;

static fe synd(int ia, int ib, const fe *wt) {
  fe acc = 0;
  for (int n = 0; n < nz_cnt; ++n) {
    int k = nz_idx[n];
    fe term = mul(recv[k], mul(basis[ia].val[k], basis[ib].val[k]));
    if (wt)
      term = mul(term, wt[k]);
    acc ^= term;
  }
  return acc;
}
static fe synd1(int ia) {
  if (use_poly)
    return basis_syndrome(ia);
  fe acc = 0;
  for (int n = 0; n < nz_cnt; ++n)
    acc ^= mul(recv[nz_idx[n]], basis[ia].val[nz_idx[n]]);
  return acc;
}

typedef struct {
  int T, nfound, resolved, overflow, failed;
  int pos[MAXT + 1];
  fe val[MAXT + 1];
  uint8_t rowset[SGMAX], colset[SGMAX];
  int vote_steps, vote_wrong, vote_tie, vote_min_margin;
  double t_main, t_extract, t_resolve;
} Result;
static Result res_fr, res_bms, res_vote, res_frvote;

/* Public API bridge: caller-owned received word and output buffers. */
static const fe *api_received = NULL;
static size_t api_received_len = 0;
static size_t *api_error_positions = NULL;
static fe *api_error_values = NULL;
static size_t api_output_capacity = 0;
static size_t *api_error_count = NULL;
static int api_mode = 0;
static int api_completed = 0;
static fe *api_owned_received = NULL;
static size_t *api_position_map = NULL;

static fe *aug =
    NULL; /* resolve_by_values: allocated to the active row window */
static fe err_val[MAXT + 1];

/* Error values: solve  sum_k e_k h_a(P_k) = S_a  (a = 0 .. m-1).
 * Only ~n rows are needed to determine the n unknowns, so forward elimination
 * is done on the first R = n + slack rows only (O(n^3/3) instead of
 * O(m n^2/2)); the remaining rows are merely verified afterwards, O(m n).
 * The inner loop is branch-free: the pivot row is kept in log form and the
 * zero element is mapped to a sentinel log whose exp-table entry is 0. */
static int resolve_by_values(int *found, int *nfound) {
  int n = *nfound, m = count_basis_upto(M_SYND);
  if (n > MAXT)
    return 0;
  const int SENT = 2 * FO;
  fe *E = malloc(sizeof(fe) * (size_t)(3 * FO + 8));
  uint16_t *lg = malloc(sizeof(uint16_t) * (size_t)(n + 2));
  fe *x = malloc(sizeof(fe) * (size_t)(n + 1));
  int ok = 0;
  if (!E || !lg || !x)
    goto done;
  for (int i = 0; i < 3 * FO + 8; ++i)
    E[i] = i < 2 * FO ? gexp[i % FO] : 0;
  int R = n + 32 < m ? n + 32 : m;
  int aug_cap = R;
  aug = calloc((size_t)aug_cap * (size_t)(n + 1), sizeof(fe));
  if (!aug)
    goto done;
  for (;;) {
    for (int a = 0; a < R; ++a) {
      for (int k = 0; k < n; ++k)
        aug[(size_t)a * (size_t)(n + 1) + (size_t)k] =
            basis_value_at(a, found[k]);
      aug[(size_t)a * (size_t)(n + 1) + (size_t)n] = synd1(a);
    }
    int retry = 0;
    for (int c = 0; c < n; ++c) {
      int p = -1;
      for (int a = c; a < R; ++a)
        if (aug[(size_t)a * (size_t)(n + 1) + (size_t)c]) {
          p = a;
          break;
        }
      if (p < 0) {
        if (R < m) {
          retry = 1; /* rare: enlarge the row window and start over */
          break;
        }
        goto done;
      }
      if (p != c)
        for (int k = c; k <= n; ++k) {
          fe t = aug[(size_t)p * (size_t)(n + 1) + (size_t)k];
          aug[(size_t)p * (size_t)(n + 1) + (size_t)k] =
              aug[(size_t)c * (size_t)(n + 1) + (size_t)k];
          aug[(size_t)c * (size_t)(n + 1) + (size_t)k] = t;
        }
      fe iv = inv(aug[(size_t)c * (size_t)(n + 1) + (size_t)c]);
      aug[(size_t)c * (size_t)(n + 1) + (size_t)c] = 1;
      lg[c] = 0;
      for (int k = c + 1; k <= n; ++k) {
        aug[(size_t)c * (size_t)(n + 1) + (size_t)k] =
            mul(aug[(size_t)c * (size_t)(n + 1) + (size_t)k], iv);
        lg[k] = aug[(size_t)c * (size_t)(n + 1) + (size_t)k]
                    ? glog[aug[(size_t)c * (size_t)(n + 1) + (size_t)k]]
                    : (uint16_t)SENT;
      }
      for (int a = c + 1; a < R; ++a) {
        fe f = aug[(size_t)a * (size_t)(n + 1) + (size_t)c];
        if (!f)
          continue;
        const int lf = glog[f];
        fe *row = aug + (size_t)a * (size_t)(n + 1);
        for (int k = c + 1; k <= n; ++k)
          row[k] ^= E[lf + lg[k]];
        row[c] = 0;
      }
    }
    if (retry) {
      R = m;
      if (R > aug_cap) {
        fe *tmp = realloc(aug, (size_t)R * (size_t)(n + 1) * sizeof(fe));
        if (!tmp)
          goto done;
        aug = tmp;
        aug_cap = R;
      }
      continue;
    }
    break;
  }
  for (int a = n; a < R; ++a)
    if (aug[(size_t)a * (size_t)(n + 1) + (size_t)n])
      goto done; /* inconsistent */
  /* back substitution */
  for (int c = n - 1; c >= 0; --c) {
    fe acc = aug[(size_t)c * (size_t)(n + 1) + (size_t)n];
    for (int k = c + 1; k < n; ++k)
      acc ^= mul(aug[(size_t)c * (size_t)(n + 1) + (size_t)k], x[k]);
    x[c] = acc;
  }
  /* verify the rows that were not used for elimination */
  for (int a = R; a < m; ++a) {
    fe acc = synd1(a);
    for (int k = 0; k < n; ++k)
      acc ^= mul(x[k], basis_value_at(a, found[k]));
    if (acc)
      goto done;
  }
  {
    int keep = 0;
    for (int k = 0; k < n; ++k)
      if (x[k]) {
        found[keep] = found[k];
        err_val[keep] = x[k];
        ++keep;
      }
    *nfound = keep;
  }
  ok = 1;
done:
  free(E);
  free(lg);
  free(x);
  free(aug);
  aug = NULL;
  return ok;
}

static void resolve_and_store(Result *R) {
  double t0 = now();
  int resolved = 0;
  if (!R->overflow && R->nfound > 0) {
    if (R->nfound > R->T)
      printf("#   %d common zeros > t = %d: resolving by error values\n",
             R->nfound, R->T);
    resolved = resolve_by_values(R->pos, &R->nfound);
    if (!resolved)
      puts("#   could not compute the error values.");
  }
  R->resolved = resolved;
  for (int m = 0; m < R->nfound; ++m)
    R->val[m] = resolved ? err_val[m] : 0;
  R->t_resolve = now() - t0;
}

/* Experimental multivariate Forney-function evaluation.
 *
 * Construct functions f_i in the available check-function space such that
 *       f_i(P_j) = delta_ij.
 * Since S(f) = sum_j e_j f(P_j), this gives e_i = S(f_i).
 *
 * This is the interpolation/Forney-function identity, not the optimized
 * Kötter-update algorithm.  It is kept separate from resolve_by_values so
 * the two implementations can be compared on identical FR locations.
 */
static int resolve_by_forney_functions(int *found, int *nfound) {
  int n = *nfound, m = count_basis_upto(M_SYND);
  if (n <= 0 || n > MAXT || m < n)
    return 0;

  int R = n + 32 < m ? n + 32 : m;
  int cap = R;
  int *row_id = malloc((size_t)cap * sizeof(*row_id));
  fe *work = calloc((size_t)cap * (size_t)n, sizeof(*work));
  fe *aug2 = NULL;
  fe *sy = malloc((size_t)n * sizeof(*sy));
  fe *vals = malloc((size_t)n * sizeof(*vals));
  int ok = 0;
  if (!row_id || !work || !sy || !vals)
    goto done;

  for (;;) {
    for (int a = 0; a < R; ++a) {
      row_id[a] = a;
      for (int k = 0; k < n; ++k)
        work[(size_t)a * n + k] = basis_value_at(a, found[k]);
    }
    int rank = 0;
    for (int c = 0; c < n; ++c) {
      int pivot = -1;
      for (int a = rank; a < R; ++a)
        if (work[(size_t)a * n + c]) {
          pivot = a;
          break;
        }
      if (pivot < 0)
        continue;
      if (pivot != rank) {
        for (int k = 0; k < n; ++k) {
          fe t = work[(size_t)pivot * n + k];
          work[(size_t)pivot * n + k] = work[(size_t)rank * n + k];
          work[(size_t)rank * n + k] = t;
        }
        int ti = row_id[pivot];
        row_id[pivot] = row_id[rank];
        row_id[rank] = ti;
      }
      fe iv = inv(work[(size_t)rank * n + c]);
      for (int k = c; k < n; ++k)
        work[(size_t)rank * n + k] = mul(work[(size_t)rank * n + k], iv);
      for (int a = rank + 1; a < R; ++a) {
        fe f = work[(size_t)a * n + c];
        if (!f)
          continue;
        for (int k = c; k < n; ++k)
          work[(size_t)a * n + k] ^= mul(f, work[(size_t)rank * n + k]);
      }
      ++rank;
      if (rank == n)
        break;
    }
    if (rank == n)
      break;
    if (R == m)
      goto done;
    int next = m;
    fe *tmp = realloc(work, (size_t)next * (size_t)n * sizeof(*work));
    int *tmpid = realloc(row_id, (size_t)next * sizeof(*row_id));
    if (!tmp || !tmpid) {
      if (tmp)
        work = tmp;
      if (tmpid)
        row_id = tmpid;
      goto done;
    }
    work = tmp;
    row_id = tmpid;
    cap = next;
    R = next;
    (void)cap;
  }

  /* V[k,c] = h_{row_id[c]}(P_k).  Invert V^T to obtain the
   * coefficients of the interpolating functions f_i. */
  aug2 = calloc((size_t)n * (size_t)(2 * n), sizeof(*aug2));
  if (!aug2)
    goto done;
  for (int k = 0; k < n; ++k) {
    for (int c = 0; c < n; ++c)
      aug2[(size_t)k * (2 * n) + c] = basis_value_at(row_id[c], found[k]);
    aug2[(size_t)k * (2 * n) + n + k] = 1;
  }
  for (int c = 0; c < n; ++c) {
    int pivot = c;
    while (pivot < n && !aug2[(size_t)pivot * (2 * n) + c])
      ++pivot;
    if (pivot == n)
      goto done;
    if (pivot != c)
      for (int k = 0; k < 2 * n; ++k) {
        fe t = aug2[(size_t)pivot * (2 * n) + k];
        aug2[(size_t)pivot * (2 * n) + k] = aug2[(size_t)c * (2 * n) + k];
        aug2[(size_t)c * (2 * n) + k] = t;
      }
    fe iv = inv(aug2[(size_t)c * (2 * n) + c]);
    for (int k = 0; k < 2 * n; ++k)
      aug2[(size_t)c * (2 * n) + k] = mul(aug2[(size_t)c * (2 * n) + k], iv);
    for (int a = 0; a < n; ++a)
      if (a != c) {
        fe f = aug2[(size_t)a * (2 * n) + c];
        if (!f)
          continue;
        for (int k = 0; k < 2 * n; ++k)
          aug2[(size_t)a * (2 * n) + k] ^=
              mul(f, aug2[(size_t)c * (2 * n) + k]);
      }
  }
  for (int c = 0; c < n; ++c)
    sy[c] = synd1(row_id[c]);
  for (int i = 0; i < n; ++i) {
    fe e = 0;
    for (int c = 0; c < n; ++c)
      e ^= mul(aug2[(size_t)c * (2 * n) + n + i], sy[c]);
    vals[i] = e;
  }

  /* Verify the recovered values against the complete available syndrome
   * sequence before accepting them.  This prevents an underdetermined
   * interpolation from being mistaken for a successful Forney evaluation. */
  for (int a = 0; a < m; ++a) {
    fe check = 0;
    for (int i = 0; i < n; ++i)
      if (vals[i])
        check ^= mul(vals[i], basis_value_at(a, found[i]));
    if (check != synd1(a))
      goto done;
  }
  {
    int keep = 0;
    for (int i = 0; i < n; ++i)
      if (vals[i]) {
        found[keep] = found[i];
        err_val[keep] = vals[i];
        ++keep;
      }
    *nfound = keep;
  }
  ok = 1;
done:
  free(row_id);
  free(work);
  free(aug2);
  free(sy);
  free(vals);
  return ok;
}

static int use_forney_values = 0;
static int use_kotter_values = 0;
static int kotter_ready = 0;
static int resolve_by_kotter(int *found, int *nfound);
static void resolve_and_store_fr_forney(Result *R) {
  if (use_kotter_values) {
    double t0 = now();
    int resolved = 0;
    if (kotter_ready && !R->overflow && R->nfound > 0)
      resolved = resolve_by_kotter(R->pos, &R->nfound);
    if (resolved)
      puts("# FR: error values via Kötter generalized Forney update");
    else if (!R->overflow && R->nfound > 0)
      puts("# FR: Kötter update failed verification; no fallback was used.");
    R->resolved = resolved;
    for (int i = 0; i < R->nfound; ++i)
      R->val[i] = resolved ? err_val[i] : 0;
    R->t_resolve = now() - t0;
  } else if (use_forney_values) {
    double t0 = now();
    int resolved = 0;
    if (!R->overflow && R->nfound > 0)
      resolved = resolve_by_forney_functions(R->pos, &R->nfound);
    if (resolved)
      puts("# FR: error values via verified multivariate Forney functions");
    else if (!R->overflow && R->nfound > 0)
      puts("# FR: Forney-function evaluation failed verification; no fallback "
           "was used.");
    R->resolved = resolved;
    for (int i = 0; i < R->nfound; ++i)
      R->val[i] = resolved ? err_val[i] : 0;
    R->t_resolve = now() - t0;
  } else {
    /* Default: one augmented elimination.  The interpolation-based Forney
     * path needs a separate row-rank pass and a second inverse, and showed
     * no speed advantage in the current tests. */
    resolve_and_store(R);
  }
}

/* ================= FR: basic algorithm + characteristic polynomial
 * ================= */
static int slot_of_row[MAXB];
static fe *E0f;
static fe *pp_S;
static inline int pp_at(int i, int j);
static int *bas_ei, *bas_ej;
static int piv_of_row[MAXB], row_of_col[MAXB];
static fe pivinv[MAXB];
static fe **coefm; /* triangular locator coefficient rows, allocated to nl */
static int coefm_rows;
static fe *fr_coef_blk;
static size_t fr_coef_stride;
static fe *red_blk;
static size_t fr_stride;
static uint8_t **plr, **plc;
static void free_pivot_rows(void) {
  if (plr) {
    for (int i = 0; i < MAXT; ++i)
      free(plr[i]);
    free(plr);
    plr = NULL;
  }
  if (plc) {
    for (int i = 0; i < MAXT; ++i)
      free(plc[i]);
    free(plc);
    plc = NULL;
  }
}

static inline fe synd_fast(int a, int b, int shi, int shj) {
  if (!bas_ei || !bas_ej) {
    fprintf(stderr, "bas_ei/bas_ej not set\n");
    abort();
  }
  int t = pp_at(bas_ei[a] + bas_ei[b] + shi, bas_ej[a] + bas_ej[b] + shj);
  if (t < 0) {
    fprintf(stderr, "synd_fast: lattice point not in table (a=%d b=%d)\n", a,
            b);
    abort();
  }
  return pp_S[t];
}

static inline fe synd_fast2(int a, int b, int shi, int shj) {
  int t = pp_at(bas_ei[a] + bas_ei[b] + shi, bas_ej[a] + bas_ej[b] + shj);
  return pp_S[t];
}

#define FR_BLK 64
/* Log values are in [0, 2*FO].  FM+1 bits are sufficient, so store
 * pivot rows packed instead of spending 16 bits per log value. */
/* plain uint16 log rows (the bit-packed variant was ~1.7x slower in the hot
 * loop and saved only ~38% of a small table) */
static inline size_t packed_log_bytes(size_t n) { return n * sizeof(uint16_t); }
static inline void packed_log_set(uint8_t *p, size_t i, uint16_t v) {
  ((uint16_t *)p)[i] = v;
}
static inline uint16_t packed_log_get(const uint8_t *p, size_t i) {
  return ((const uint16_t *)p)[i];
}
static int fr_pos[MAXB], fr_ncol[MAXB];

/* ---- SIMD row update:  dst[0..n) ^= c * src[0..n)   (GF(2^m), m <= 12) ----
 * The hot loop of the Feng-Rao reduction is "reduce a row by a pivot row".
 * The scalar form (exp[lf + log[x]]) is a gather and cannot be vectorised.
 * Instead multiplication by the constant c is split into three 4-bit
 * nibbles; each nibble gives a 16-entry table (low byte / high byte of
 * c*nibble) that PSHUFB evaluates for 16 (AVX2) or 32 (AVX-512BW) field
 * elements at once.  Pivot rows are then stored as plain field values.
 * Runtime dispatch; falls back to the scalar log-table code.  Set the
 * environment variable AG_NOSIMD=1 to force the scalar path, AG_SIMD=2 to
 * force AVX2 even if AVX-512 is available. */
static int fr_simd = 0; /* 0 scalar (logs), 1 AVX2, 2 AVX-512BW, 3 SSSE3 */
/* fr_tab[c] = 6 nibble tables (3 low-byte, 3 high-byte) of the map x -> c*x */
static uint8_t (*fr_tab)[96];
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define FR_X86_SIMD 1
static inline void fr_nib_tables(fe c, uint8_t tl[3][16], uint8_t th[3][16]) {
  unsigned bs[12], x = c;
  for (int b = 0; b < 12; ++b) {
    bs[b] = x;
    x <<= 1;
    if (x & (unsigned)FS)
      x ^= FPOLY;
  }
  for (int k = 0; k < 3; ++k) {
    unsigned v[16];
    v[0] = 0;
    for (int n = 1; n < 16; ++n)
      v[n] = v[n & (n - 1)] ^ bs[4 * k + __builtin_ctz((unsigned)n)];
    for (int n = 0; n < 16; ++n) {
      tl[k][n] = (uint8_t)(v[n] & 0xFF);
      th[k][n] = (uint8_t)(v[n] >> 8);
    }
  }
}
__attribute__((target("avx2"))) static void fr_axpy_avx2(fe *dst, const fe *src,
                                                         fe c, size_t n) {
  const uint8_t (*tl)[16] = (const uint8_t (*)[16])fr_tab[c];
  const uint8_t (*th)[16] = tl + 3;
#define LD128(p)                                                               \
  _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(p)))
  const __m256i L0 = LD128(tl[0]), L1 = LD128(tl[1]), L2 = LD128(tl[2]);
  const __m256i H0 = LD128(th[0]), H1 = LD128(th[1]), H2 = LD128(th[2]);
#undef LD128
  const __m256i m4 = _mm256_set1_epi16(0x000F);
  size_t i = 0;
  for (; i + 16 <= n; i += 16) {
    __m256i x = _mm256_loadu_si256((const __m256i *)(src + i));
    __m256i n0 = _mm256_and_si256(x, m4);
    __m256i n1 = _mm256_and_si256(_mm256_srli_epi16(x, 4), m4);
    __m256i n2 = _mm256_srli_epi16(x, 8); /* values < 2^12 */
    __m256i lo = _mm256_xor_si256(_mm256_xor_si256(_mm256_shuffle_epi8(L0, n0),
                                                   _mm256_shuffle_epi8(L1, n1)),
                                  _mm256_shuffle_epi8(L2, n2));
    __m256i hi = _mm256_xor_si256(_mm256_xor_si256(_mm256_shuffle_epi8(H0, n0),
                                                   _mm256_shuffle_epi8(H1, n1)),
                                  _mm256_shuffle_epi8(H2, n2));
    __m256i p = _mm256_xor_si256(lo, _mm256_slli_epi16(hi, 8));
    __m256i d = _mm256_loadu_si256((const __m256i *)(dst + i));
    _mm256_storeu_si256((__m256i *)(dst + i), _mm256_xor_si256(d, p));
  }
  for (; i < n; ++i)
    dst[i] ^= mul(c, src[i]);
}
__attribute__((target("avx512f,avx512bw"))) static void
fr_axpy_avx512(fe *dst, const fe *src, fe c, size_t n) {
  const uint8_t (*tl)[16] = (const uint8_t (*)[16])fr_tab[c];
  const uint8_t (*th)[16] = tl + 3;
#define LD512(p) _mm512_broadcast_i32x4(_mm_loadu_si128((const __m128i *)(p)))
  const __m512i L0 = LD512(tl[0]), L1 = LD512(tl[1]), L2 = LD512(tl[2]);
  const __m512i H0 = LD512(th[0]), H1 = LD512(th[1]), H2 = LD512(th[2]);
#undef LD512
  const __m512i m4 = _mm512_set1_epi16(0x000F);
  size_t i = 0;
  for (; i + 32 <= n; i += 32) {
    __m512i x = _mm512_loadu_si512((const void *)(src + i));
    __m512i n0 = _mm512_and_si512(x, m4);
    __m512i n1 = _mm512_and_si512(_mm512_srli_epi16(x, 4), m4);
    __m512i n2 = _mm512_srli_epi16(x, 8);
    __m512i lo = _mm512_xor_si512(_mm512_xor_si512(_mm512_shuffle_epi8(L0, n0),
                                                   _mm512_shuffle_epi8(L1, n1)),
                                  _mm512_shuffle_epi8(L2, n2));
    __m512i hi = _mm512_xor_si512(_mm512_xor_si512(_mm512_shuffle_epi8(H0, n0),
                                                   _mm512_shuffle_epi8(H1, n1)),
                                  _mm512_shuffle_epi8(H2, n2));
    __m512i p = _mm512_xor_si512(lo, _mm512_slli_epi16(hi, 8));
    __m512i d = _mm512_loadu_si512((const void *)(dst + i));
    _mm512_storeu_si512((void *)(dst + i), _mm512_xor_si512(d, p));
  }
  for (; i < n; ++i)
    dst[i] ^= mul(c, src[i]);
}
/* Planar variant: a row is stored as a low-byte plane followed (at byte
 * offset `hoff`) by a high-byte plane.  One PSHUFB then handles 32 (AVX2) /
 * 64 (AVX-512BW) elements -- twice the elements per shuffle of the packed
 * kernels above, and the shuffle port is the bottleneck. */
__attribute__((target("avx2"))) static void
fr_axpyp_avx2(uint8_t *dl, uint8_t *dh, const uint8_t *sl, const uint8_t *sh,
              fe c, size_t n) {
  const uint8_t (*tl)[16] = (const uint8_t (*)[16])fr_tab[c];
  const uint8_t (*th)[16] = tl + 3;
#define LD128(p)                                                               \
  _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(p)))
  const __m256i L0 = LD128(tl[0]), L1 = LD128(tl[1]), L2 = LD128(tl[2]);
  const __m256i H0 = LD128(th[0]), H1 = LD128(th[1]), H2 = LD128(th[2]);
#undef LD128
  const __m256i m4 = _mm256_set1_epi8(0x0F);
  size_t i = 0;
  /* n is rounded up: rows carry >= 64 bytes of zero/garbage padding */
  for (; i < n; i += 32) {
    __m256i xl = _mm256_loadu_si256((const __m256i *)(sl + i));
    __m256i xh = _mm256_loadu_si256((const __m256i *)(sh + i)); /* <= 0x0F */
    __m256i n0 = _mm256_and_si256(xl, m4);
    __m256i n1 = _mm256_and_si256(_mm256_srli_epi16(xl, 4), m4);
    __m256i lo = _mm256_xor_si256(_mm256_xor_si256(_mm256_shuffle_epi8(L0, n0),
                                                   _mm256_shuffle_epi8(L1, n1)),
                                  _mm256_shuffle_epi8(L2, xh));
    __m256i hi = _mm256_xor_si256(_mm256_xor_si256(_mm256_shuffle_epi8(H0, n0),
                                                   _mm256_shuffle_epi8(H1, n1)),
                                  _mm256_shuffle_epi8(H2, xh));
    _mm256_storeu_si256(
        (__m256i *)(dl + i),
        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)(dl + i)), lo));
    _mm256_storeu_si256(
        (__m256i *)(dh + i),
        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)(dh + i)), hi));
  }
}
__attribute__((target("avx512f,avx512bw"))) static void
fr_axpyp_avx512(uint8_t *dl, uint8_t *dh, const uint8_t *sl, const uint8_t *sh,
                fe c, size_t n) {
  const uint8_t (*tl)[16] = (const uint8_t (*)[16])fr_tab[c];
  const uint8_t (*th)[16] = tl + 3;
#define LD512(p) _mm512_broadcast_i32x4(_mm_loadu_si128((const __m128i *)(p)))
  const __m512i L0 = LD512(tl[0]), L1 = LD512(tl[1]), L2 = LD512(tl[2]);
  const __m512i H0 = LD512(th[0]), H1 = LD512(th[1]), H2 = LD512(th[2]);
#undef LD512
  const __m512i m4 = _mm512_set1_epi8(0x0F);
  size_t i = 0;
  for (; i < n; i += 64) {
    __m512i xl = _mm512_loadu_si512((const void *)(sl + i));
    __m512i xh = _mm512_loadu_si512((const void *)(sh + i));
    __m512i n0 = _mm512_and_si512(xl, m4);
    __m512i n1 = _mm512_and_si512(_mm512_srli_epi16(xl, 4), m4);
    __m512i lo = _mm512_xor_si512(_mm512_xor_si512(_mm512_shuffle_epi8(L0, n0),
                                                   _mm512_shuffle_epi8(L1, n1)),
                                  _mm512_shuffle_epi8(L2, xh));
    __m512i hi = _mm512_xor_si512(_mm512_xor_si512(_mm512_shuffle_epi8(H0, n0),
                                                   _mm512_shuffle_epi8(H1, n1)),
                                  _mm512_shuffle_epi8(H2, xh));
    _mm512_storeu_si512(
        (void *)(dl + i),
        _mm512_xor_si512(_mm512_loadu_si512((const void *)(dl + i)), lo));
    _mm512_storeu_si512(
        (void *)(dh + i),
        _mm512_xor_si512(_mm512_loadu_si512((const void *)(dh + i)), hi));
  }
}
/* SSSE3 (128-bit PSHUFB): for CPUs without AVX2 (e.g. Athlon X4 / Steamroller)
 */
__attribute__((target("ssse3"))) static void
fr_axpyp_ssse3(uint8_t *dl, uint8_t *dh, const uint8_t *sl, const uint8_t *sh,
               fe c, size_t n) {
  const uint8_t (*tl)[16] = (const uint8_t (*)[16])fr_tab[c];
  const uint8_t (*th)[16] = tl + 3;
  const __m128i L0 = _mm_loadu_si128((const __m128i *)tl[0]),
                L1 = _mm_loadu_si128((const __m128i *)tl[1]),
                L2 = _mm_loadu_si128((const __m128i *)tl[2]);
  const __m128i H0 = _mm_loadu_si128((const __m128i *)th[0]),
                H1 = _mm_loadu_si128((const __m128i *)th[1]),
                H2 = _mm_loadu_si128((const __m128i *)th[2]);
  const __m128i m4 = _mm_set1_epi8(0x0F);
  for (size_t i = 0; i < n; i += 16) {
    __m128i xl = _mm_loadu_si128((const __m128i *)(sl + i));
    __m128i xh = _mm_loadu_si128((const __m128i *)(sh + i));
    __m128i n0 = _mm_and_si128(xl, m4);
    __m128i n1 = _mm_and_si128(_mm_srli_epi16(xl, 4), m4);
    __m128i lo = _mm_xor_si128(
        _mm_xor_si128(_mm_shuffle_epi8(L0, n0), _mm_shuffle_epi8(L1, n1)),
        _mm_shuffle_epi8(L2, xh));
    __m128i hi = _mm_xor_si128(
        _mm_xor_si128(_mm_shuffle_epi8(H0, n0), _mm_shuffle_epi8(H1, n1)),
        _mm_shuffle_epi8(H2, xh));
    _mm_storeu_si128(
        (__m128i *)(dl + i),
        _mm_xor_si128(_mm_loadu_si128((const __m128i *)(dl + i)), lo));
    _mm_storeu_si128(
        (__m128i *)(dh + i),
        _mm_xor_si128(_mm_loadu_si128((const __m128i *)(dh + i)), hi));
  }
}
__attribute__((target("ssse3"))) static void
fr_axpy_ssse3(fe *dst, const fe *src, fe c, size_t n) {
  const uint8_t (*tl)[16] = (const uint8_t (*)[16])fr_tab[c];
  const uint8_t (*th)[16] = tl + 3;
  const __m128i L0 = _mm_loadu_si128((const __m128i *)tl[0]),
                L1 = _mm_loadu_si128((const __m128i *)tl[1]),
                L2 = _mm_loadu_si128((const __m128i *)tl[2]);
  const __m128i H0 = _mm_loadu_si128((const __m128i *)th[0]),
                H1 = _mm_loadu_si128((const __m128i *)th[1]),
                H2 = _mm_loadu_si128((const __m128i *)th[2]);
  const __m128i m4 = _mm_set1_epi16(0x000F);
  size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    __m128i x = _mm_loadu_si128((const __m128i *)(src + i));
    __m128i n0 = _mm_and_si128(x, m4);
    __m128i n1 = _mm_and_si128(_mm_srli_epi16(x, 4), m4);
    __m128i n2 = _mm_srli_epi16(x, 8);
    __m128i lo = _mm_xor_si128(
        _mm_xor_si128(_mm_shuffle_epi8(L0, n0), _mm_shuffle_epi8(L1, n1)),
        _mm_shuffle_epi8(L2, n2));
    __m128i hi = _mm_xor_si128(
        _mm_xor_si128(_mm_shuffle_epi8(H0, n0), _mm_shuffle_epi8(H1, n1)),
        _mm_shuffle_epi8(H2, n2));
    __m128i p = _mm_xor_si128(lo, _mm_slli_epi16(hi, 8));
    _mm_storeu_si128(
        (__m128i *)(dst + i),
        _mm_xor_si128(_mm_loadu_si128((const __m128i *)(dst + i)), p));
  }
  for (; i < n; ++i)
    dst[i] ^= mul(c, src[i]);
}
static int fr_detect_simd(void) {
  if (getenv("AG_NOSIMD") || sizeof(fe) != 2 || FM > 12)
    return 0;
  __builtin_cpu_init();
  const char *f = getenv("AG_SIMD");
  int want512 = !(f && f[0] == '2');
  int forcessse = (f && f[0] == '3');
  if (!forcessse && want512 && __builtin_cpu_supports("avx512f") &&
      __builtin_cpu_supports("avx512bw"))
    return 2;
  if (!forcessse && __builtin_cpu_supports("avx2"))
    return 1;
  if (__builtin_cpu_supports("ssse3"))
    return 3;
  return 0;
}
static inline void fr_axpyp(uint8_t *dl, uint8_t *dh, const uint8_t *sl,
                            const uint8_t *sh, fe c, size_t n) {
  if (fr_simd == 2)
    fr_axpyp_avx512(dl, dh, sl, sh, c, n);
  else if (fr_simd == 1)
    fr_axpyp_avx2(dl, dh, sl, sh, c, n);
  else
    fr_axpyp_ssse3(dl, dh, sl, sh, c, n);
}
static inline void fr_axpy(fe *dst, const fe *src, fe c, size_t n) {
  if (fr_simd == 2)
    fr_axpy_avx512(dst, src, c, n);
  else if (fr_simd == 1)
    fr_axpy_avx2(dst, src, c, n);
  else
    fr_axpy_ssse3(dst, src, c, n);
}
#else
static int fr_detect_simd(void) { return 0; }
static inline void fr_axpyp(uint8_t *dl, uint8_t *dh, const uint8_t *sl,
                            const uint8_t *sh, fe c, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    fe p = mul(c, (fe)(sl[i] | (sh[i] << 8)));
    dl[i] ^= (uint8_t)(p & 0xFF);
    dh[i] ^= (uint8_t)(p >> 8);
  }
}
static inline void fr_axpy(fe *dst, const fe *src, fe c, size_t n) {
  for (size_t i = 0; i < n; ++i)
    dst[i] ^= mul(c, src[i]);
}
#endif

/* reduction rows are planar when fr_simd (low-byte plane, then the high-byte
 * plane at byte offset fr_stride); plain uint16 otherwise */
static inline fe fr_get(const fe *R, int j) {
  if (fr_simd) {
    const uint8_t *b = (const uint8_t *)R;
    return (fe)(b[j] | ((unsigned)b[(size_t)j + fr_stride] << 8));
  }
  return R[j];
}
static inline void fr_set(fe *R, int j, fe v) {
  if (fr_simd) {
    uint8_t *b = (uint8_t *)R;
    b[j] = (uint8_t)(v & 0xFF);
    b[(size_t)j + fr_stride] = (uint8_t)(v >> 8);
  } else
    R[j] = v;
}
/* 行 i を列 j0
 * から簡約する。ピボット列を持たない非零の列で止まり、その列を返す。
 * 最後まで消えたら ncol を返す。row_of_col / plr / plc は読むだけ。
 * fr_simd のときピボット行 plr は「値」(uint16)、そうでなければ「log」。 */
/* apply pivot p (column j) to one row; cc = Ri[j] / pivot */
static inline void fr_apply(fe *Ri, fe *Ci, int j, int ncol, int needc, int p,
                            fe cc) {
  int s = slot_of_row[p];
  int lf = glog[cc];
  if (fr_simd) {
    uint8_t *rb = (uint8_t *)Ri;
    const uint8_t *pb = (const uint8_t *)plr[s];
    size_t po = ((size_t)fr_ncol[p] + 63) & ~(size_t)63;
    fr_axpyp(rb + j, rb + fr_stride + j, pb + j, pb + po + j, cc,
             (size_t)(ncol - j));
    if (needc)
      fr_axpy(Ci, (const fe *)plc[s], cc, (size_t)p + 1);
    return;
  } else {
    const uint8_t *lr = plr[s];
    for (int jj = j; jj < ncol; ++jj)
      Ri[jj] ^= E0f[lf + packed_log_get(lr, (size_t)jj)];
  }
  if (needc) {
    const uint8_t *lc = plc[s];
    for (int k = 0; k <= p; ++k)
      Ci[k] ^= E0f[lf + packed_log_get(lc, (size_t)k)];
  }
}
static int fr_run(fe *Ri, fe *Ci, int ncol, int needc, int j0) {
  for (int j = j0; j < ncol; ++j) {
    fe v = fr_get(Ri, j);
    if (!v)
      continue;
    int p = row_of_col[j];
    if (p < 0)
      return j;
    fr_apply(Ri, Ci, j, ncol, needc, p, mul(v, pivinv[p]));
  }
  return ncol;
}

#define FR_GMAX 32
static int fr_group = 8; /* rows reduced together (cache blocking) */
/* Reduce rows ia..ib-1 against the finished pivots, TOGETHER and in pivot
 * column order: pivot row p is then read from memory once and reused for
 * every row of the group while it sits in L1, instead of being streamed
 * once per row.  Per row the sequence of operations is identical to fr_run
 * (rows are independent), so results are bit-identical.  fr_pos[i] receives
 * the column where row i stops (a non-zero entry without a pivot, or ncol). */
static void fr_run_group(int b0, int ia, int ib, int nl) {
  int n = ib - ia;
  fe *R[FR_GMAX], *C[FR_GMAX];
  int nc[FR_GMAX], needc[FR_GMAX], act[FR_GMAX], maxn = 0, nact = n;
  for (int r = 0; r < n; ++r) {
    int i = ia + r;
    R[r] = red_blk + (size_t)(i - b0) * fr_stride;
    needc[r] = (i < nl);
    C[r] = compact_fr_mode ? fr_coef_blk + (size_t)(i - b0) * fr_coef_stride
                           : (needc[r] ? coefm[i] : NULL);
    nc[r] = fr_ncol[i];
    act[r] = 1;
    if (nc[r] > maxn)
      maxn = nc[r];
  }
  for (int j = 0; j < maxn && nact; ++j) {
    int p = row_of_col[j];
    for (int r = 0; r < n; ++r) {
      if (!act[r])
        continue;
      if (j >= nc[r]) {
        act[r] = 0;
        nact--;
        fr_pos[ia + r] = nc[r];
        continue;
      }
      fe v = fr_get(R[r], j);
      if (!v)
        continue;
      if (p < 0) {
        act[r] = 0;
        nact--;
        fr_pos[ia + r] = j;
        continue;
      }
      fr_apply(R[r], C[r], j, nc[r], needc[r], p, mul(v, pivinv[p]));
    }
  }
  for (int r = 0; r < n; ++r)
    if (act[r])
      fr_pos[ia + r] = nc[r];
}

static int alloc_coefm(int nl) {
  if (coefm) {
    for (int i = 0; i < coefm_rows; ++i)
      free(coefm[i]);
    free(coefm);
    coefm = NULL;
  }
  coefm = calloc((size_t)nl, sizeof(*coefm));
  if (!coefm)
    return -1;
  coefm_rows = nl;
  for (int i = 0; i < nl; ++i) {
    coefm[i] = calloc((size_t)i + 1, sizeof(**coefm));
    if (!coefm[i]) {
      for (int j = 0; j < i; ++j)
        free(coefm[j]);
      free(coefm);
      coefm = NULL;
      coefm_rows = 0;
      return -1;
    }
  }
  return 0;
}

static void free_coefm(void) {
  if (coefm) {
    for (int i = 0; i < coefm_rows; ++i)
      free(coefm[i]);
    free(coefm);
    coefm = NULL;
  }
  coefm_rows = 0;
}

static int fr_locator_is_minimal(int basis_row) {
  int b = basis[basis_row].rho;
  for (int q = 0; q < ngenH; ++q) {
    int a = b - genH[q];
    if (a >= 0 && a < SGMAX && in_H[a]) {
      int ia = idx_of_ord[a];
      if (ia >= 0 && ia < basis_row && piv_of_row[ia] < 0)
        return 0;
    }
  }
  return 1;
}

static int feng_rao_reduce(int row_limit) {
  int nrows = count_basis_upto(ROW_MAX);
  if (row_limit > 0 && row_limit < nrows)
    nrows = row_limit;
  int nl = count_basis_upto(L_DEG);
  int nslot = 0;
  free_coefm();
  if (compact_fr_mode) {
    coefm = calloc((size_t)nl, sizeof(*coefm));
    if (!coefm)
      return -1;
    coefm_rows = nl;
    fr_coef_stride = (size_t)(nl > 0 ? nl : 1);
    free(fr_coef_blk);
    fr_coef_blk = calloc((size_t)FR_BLK * fr_coef_stride, sizeof(fe));
    if (!fr_coef_blk)
      return -1;
  } else if (alloc_coefm(nl) != 0)
    return -1;
  int maxcol = count_basis_upto(M_SYND - basis[0].rho);
  if (maxcol < 1)
    maxcol = 1;
  fr_stride = ((size_t)maxcol + 63) / 64 * 64 + 64; /* padded: SIMD tails */
  free(red_blk);
  red_blk = calloc((size_t)FR_BLK * fr_stride, sizeof(fe));
  if (!red_blk)
    return -1;
  if (!E0f) {
    E0f = calloc((size_t)4 * FO + 8, sizeof(fe));
    for (int k = 0; k < 2 * FO; ++k)
      E0f[k] = gexp[k];
    if (!E0f)
      return -1;
  }
  fr_simd = fr_detect_simd();
  if (fr_simd && !fr_tab) {
    fr_tab = malloc((size_t)FS * sizeof(*fr_tab));
    if (!fr_tab)
      fr_simd = 0;
    else
      for (int c = 0; c < FS; ++c) {
        uint8_t tl[3][16], th[3][16];
        fr_nib_tables((fe)c, tl, th);
        memcpy(fr_tab[c], tl, 48);
        memcpy(fr_tab[c] + 48, th, 48);
      }
  }
  if (getenv("AG_FRG"))
    fr_group = atoi(getenv("AG_FRG"));
  {
    static int said = 0;
    if (!said) {
      said = 1;
      printf("# SIMD for FR reduction: %s\n", fr_simd == 2   ? "AVX-512BW"
                                              : fr_simd == 1 ? "AVX2"
                                              : fr_simd == 3
                                                  ? "SSSE3"
                                                  : "none (scalar tables)");
    }
  }
  free_pivot_rows();
  plr = calloc((size_t)MAXT, sizeof(*plr));
  plc = calloc((size_t)MAXT, sizeof(*plc));
  if (!plr || !plc) {
    free_pivot_rows();
    return -1;
  }
  for (int j = 0; j < MAXB; ++j)
    row_of_col[j] = -1;

  for (int b0 = 0; b0 < nrows; b0 += FR_BLK) {
    int b1 = b0 + FR_BLK < nrows ? b0 + FR_BLK : nrows;

    /* フェーズA: 成分作り + 確定済みピボットでの簡約(並列) */
    {
      int G = fr_group < 1 ? 1 : (fr_group > FR_GMAX ? FR_GMAX : fr_group);
      int ngrp = (b1 - b0 + G - 1) / G;
#pragma omp parallel for schedule(dynamic, 1)
      for (int g = 0; g < ngrp; ++g) {
        int ia = b0 + g * G, ib = ia + G < b1 ? ia + G : b1;
        for (int i = ia; i < ib; ++i) {
          int ncol = count_basis_upto(M_SYND - basis[i].rho);
          int needc = (i < nl);
          fr_ncol[i] = ncol;
          piv_of_row[i] = -1;
          fe *Ri = red_blk + (size_t)(i - b0) * fr_stride;
          fe *Ci = compact_fr_mode
                       ? fr_coef_blk + (size_t)(i - b0) * fr_coef_stride
                       : (needc ? coefm[i] : NULL);
          for (int j = 0; j < ncol; ++j)
            fr_set(Ri, j, use_poly ? synd_fast(i, j, 0, 0) : synd(i, j, NULL));
          if (needc) {
            memset(Ci, 0, (size_t)(i + 1) * sizeof(fe));
            Ci[i] = 1;
          }
        }
        fr_run_group(b0, ia, ib, nl);
      }
    }

    /* フェーズB: ブロック内の新ピボットを使って続行(逐次) */
    for (int i = b0; i < b1; ++i) {
      int ncol = fr_ncol[i];
      int needc = (i < nl);
      fe *Ri = red_blk + (size_t)(i - b0) * fr_stride;
      fe *Ci = compact_fr_mode ? fr_coef_blk + (size_t)(i - b0) * fr_coef_stride
                               : (needc ? coefm[i] : NULL);
      int j = fr_run(Ri, Ci, ncol, needc, fr_pos[i]);
      if (j >= ncol) {
        /* Dependent row: FR-VOTE only retains coefficient vectors for the
         * minimal non-pivot generators that will actually become locators. */
        if (compact_fr_mode && needc && fr_locator_is_minimal(i)) {
          coefm[i] = malloc((size_t)(i + 1) * sizeof(fe));
          if (!coefm[i])
            return -1;
          memcpy(coefm[i], Ci, (size_t)(i + 1) * sizeof(fe));
        }
        continue;
      }
      if (nslot >= MAXT - 1)
        return -1; /* ピボット過多 */
      int s = nslot++;
      slot_of_row[i] = s;
      size_t po = ((size_t)ncol + 63) & ~(size_t)63;
      plr[s] = calloc(1, fr_simd ? 2 * po + 64
                                 : packed_log_bytes((size_t)(ncol ? ncol : 1)));
      if (!plr[s])
        return -1;
      if (fr_simd) {
        memcpy(plr[s], (const uint8_t *)Ri, (size_t)ncol);
        memcpy((uint8_t *)plr[s] + po, (const uint8_t *)Ri + fr_stride,
               (size_t)ncol);
      } else
        for (int jj = 0; jj < ncol; ++jj)
          packed_log_set(plr[s], (size_t)jj,
                         Ri[jj] ? glog[Ri[jj]] : (uint16_t)(2 * FO));
      if (needc) {
        plc[s] = calloc(1, packed_log_bytes((size_t)(i + 1)) + 64);
        if (!plc[s])
          return -1;
        if (fr_simd)
          memcpy(plc[s], Ci, (size_t)(i + 1) * sizeof(fe));
        else
          for (int k = 0; k <= i; ++k)
            packed_log_set(plc[s], (size_t)k,
                           Ci[k] ? glog[Ci[k]] : (uint16_t)(2 * FO));
      }
      piv_of_row[i] = j;
      row_of_col[j] = i;
      pivinv[i] = inv(fr_get(Ri, j));
    }
  }
  if (compact_fr_mode) {
    free(fr_coef_blk);
    fr_coef_blk = NULL;
    fr_coef_stride = 0;
  }
  return nrows;
}

static int feng_rao_reduce2(void) {
  int nrows = count_basis_upto(ROW_MAX);
  int maxcol = count_basis_upto(M_SYND - basis[0].rho);
  if (maxcol < 1)
    maxcol = 1;
  fe *red2 = calloc((size_t)nrows * (size_t)maxcol, sizeof(fe));
  if (!red2)
    return -1;
  int nl = count_basis_upto(L_DEG);
  int nslot = 0;
  if (alloc_coefm(nl) != 0) {
    free(red2);
    return -1;
  }
  if (!E0f) {
    E0f = calloc((size_t)4 * FO + 8, sizeof(fe));
    for (int k = 0; k < 2 * FO; ++k)
      E0f[k] = gexp[k];
    if (!E0f) {
      free(red2);
      return -1;
    }
  }
  free_pivot_rows();
  plr = calloc((size_t)MAXT, sizeof(*plr));
  plc = calloc((size_t)MAXT, sizeof(*plc));
  if (!plr || !plc) {
    free(red2);
    free_pivot_rows();
    return -1;
  }
  for (int j = 0; j < MAXB; ++j)
    row_of_col[j] = -1;

  for (int i = 0; i < nrows; ++i) {
    int ncol = count_basis_upto(M_SYND - basis[i].rho);
    int needc = (i < nl);
    for (int j = 0; j < ncol; ++j)
      red2[(size_t)i * (size_t)maxcol + (size_t)j] =
          use_poly ? synd_fast(i, j, 0, 0) : synd(i, j, NULL);
    if (needc) {
      memset(coefm[i], 0, (size_t)(i + 1) * sizeof(fe));
      coefm[i][i] = 1;
    }
    piv_of_row[i] = -1;

    for (int j = 0; j < ncol; ++j) {
      fe v = red2[(size_t)i * (size_t)maxcol + (size_t)j];
      if (!v)
        continue;
      int p = row_of_col[j];
      if (p >= 0) {
        int s = slot_of_row[p];
        int lf = glog[mul(v, pivinv[p])];
        const uint8_t *lr = plr[s];
        fe *Ri = red2 + (size_t)i * (size_t)maxcol;
        for (int jj = j; jj < ncol; ++jj)
          Ri[jj] ^= E0f[lf + packed_log_get(lr, (size_t)jj)];
        if (needc) {
          const uint8_t *lc = plc[s];
          fe *Ci = coefm[i];
          for (int k = 0; k <= p; ++k)
            Ci[k] ^= E0f[lf + packed_log_get(lc, (size_t)k)];
        }
      } else {
        if (nslot >= MAXT - 1)
          return -1; /* ピボット過多 */
        int s = nslot++;
        slot_of_row[i] = s;
        plr[s] = calloc(1, packed_log_bytes((size_t)(ncol ? ncol : 1)));
        if (!plr[s]) {
          free(red2);
          return -1;
        }
        for (int jj = 0; jj < ncol; ++jj)
          packed_log_set(
              plr[s], (size_t)jj,
              red2[(size_t)i * (size_t)maxcol + (size_t)jj]
                  ? glog[red2[(size_t)i * (size_t)maxcol + (size_t)jj]]
                  : (uint16_t)(2 * FO));
        if (needc) {
          plc[s] = calloc(1, packed_log_bytes((size_t)(i + 1)));
          if (!plc[s]) {
            free(red2);
            return -1;
          }
          for (int k = 0; k <= i; ++k)
            packed_log_set(plc[s], (size_t)k,
                           coefm[i][k] ? glog[coefm[i][k]]
                                       : (uint16_t)(2 * FO));
        }
        piv_of_row[i] = j;
        row_of_col[j] = i;
        pivinv[i] = inv(v);
        break;
      }
    }
  }
  free(red2);
  return nrows;
}

/* ================= FR: direct locator evaluation (no charpoly) ======= */
static int true_pos[MAXT];
static void report_fr_locator_zero_counts(int nl, const Result *R) {
  int locno = 0;
  puts("# FR locator zero counts on the evaluation-point set:");
  for (int i = 0; i < nl; ++i) {
    if (piv_of_row[i] >= 0)
      continue;
    int zeros = 0, true_hits = 0;
    for (int p = 0; p < npts; ++p) {
      fe acc = 0;
      for (int k = 0; k <= i; ++k)
        if (coefm[i][k])
          acc ^= mul(coefm[i][k], basis[k].val[p]);
      if (!acc) {
        ++zeros;
        for (int e = 0; e < NUM_ERRORS; ++e)
          if (true_pos[e] == p) {
            ++true_hits;
            break;
          }
      }
    }
    printf("#   L_%d: zeros=%d, injected-error points among them=%d\n", locno++,
           zeros, true_hits);
  }
  int common_true = 0;
  for (int e = 0; e < NUM_ERRORS; ++e)
    for (int k = 0; k < R->nfound; ++k)
      if (true_pos[e] == R->pos[k]) {
        ++common_true;
        break;
      }
  printf("# FR common zero set: %d points; injected-error points in it=%d/%d\n",
         R->nfound, common_true, NUM_ERRORS);
}

static void print_fr_locator_equations(int nl, int quiet) {
  if (!show_locator_equations || quiet)
    return;
  puts("# FR error-position equations (common-zero system):");
  int count = 0;
  for (int i = 0; i < nl; ++i) {
    if (piv_of_row[i] >= 0)
      continue;
    printf("#   L_%d(P) = ", count++);
    int first = 1;
    for (int k = 0; k <= i; ++k) {
      fe c = coefm[i][k];
      if (!c)
        continue;
      printf("%s0x%03X*h_%d(P)", first ? "" : " + ", (unsigned)c, k);
      first = 0;
    }
    if (first)
      printf("0");
    puts(" = 0");
  }
  puts("# Error positions are the common zeros of these locator equations.");
}

/* Zero set of the FR error locators, vectorised.
 * Instead of evaluating every locator at every point with scalar multiplies
 * and strided reads of basis[k].val[idx], (1) the first locator is evaluated
 * at all points by accumulating  acc += c_k * val_k  (contiguous SIMD axpy),
 * which leaves only the true error positions plus a few false zeros; (2) the
 * remaining locators are evaluated on those survivors only, again as axpy over
 * a compact [k][survivor] matrix.  Returns 0 on success, -1 if not applicable
 * (caller falls back to the scalar loop). */
static int fr_extract_simd(int nl, int *pos, int cap, int *nfound,
                           int *overflow) {
  if (!fr_simd || npts < 64)
    return -1;
  int i0 = -1, imax = -1;
  for (int i = 0; i < nl; ++i)
    if (piv_of_row[i] < 0) {
      if (i0 < 0)
        i0 = i;
      imax = i;
    }
  if (i0 < 0)
    return -1;
  double tx0 = now(), tx1, tx2;
  fe *acc = calloc((size_t)npts + 64, sizeof(fe));
  int *S = malloc(sizeof(int) * (size_t)npts);
  if (!acc || !S) {
    free(acc);
    free(S);
    return -1;
  }
  for (int k = 0; k <= i0; ++k)
    if (coefm[i0][k])
      fr_axpy(acc, basis[k].val, coefm[i0][k], (size_t)npts);
  int ns = 0;
  for (int idx = 0; idx < npts; ++idx)
    if (!acc[idx])
      S[ns++] = idx;
  free(acc);
  tx1 = now();
  if ((double)ns * (double)(imax + 1) > 6e7) { /* matrix too large */
    if (getenv("AG_XDEBUG"))
      fprintf(stderr, "[extract] FALLBACK scalar: ns=%d imax=%d\n", ns, imax);
    free(S);
    return -1;
  }
  fe *Vs = malloc((size_t)(imax + 1) * (size_t)(ns ? ns : 1) * sizeof(fe));
  fe *a2 = calloc((size_t)(ns ? ns : 1) + 64, sizeof(fe));
  char *alive = malloc((size_t)(ns ? ns : 1));
  if (!Vs || !a2 || !alive) {
    free(S);
    free(Vs);
    free(a2);
    free(alive);
    return -1;
  }
  for (int k = 0; k <= imax; ++k)
    for (int t = 0; t < ns; ++t)
      Vs[(size_t)k * ns + t] = basis[k].val[S[t]];
  tx2 = now();
  memset(alive, 1, (size_t)ns);
  int nalive = ns;
  for (int i = i0 + 1; i <= imax && nalive; ++i) {
    if (piv_of_row[i] >= 0)
      continue;
    memset(a2, 0, (size_t)ns * sizeof(fe));
    for (int k = 0; k <= i; ++k)
      if (coefm[i][k])
        fr_axpy(a2, Vs + (size_t)k * ns, coefm[i][k], (size_t)ns);
    for (int t = 0; t < ns; ++t)
      if (alive[t] && a2[t]) {
        alive[t] = 0;
        nalive--;
      }
  }
  int nf = 0, ov = 0;
  for (int t = 0; t < ns; ++t)
    if (alive[t]) {
      if (nf < cap)
        pos[nf++] = S[t];
      else
        ov = 1;
    }
  if (getenv("AG_XDEBUG"))
    fprintf(stderr,
            "[extract] simd mode %d: first row %.3f s (survivors ns=%d), "
            "matrix %.3f s (imax=%d), rows %.3f s\n",
            fr_simd, tx1 - tx0, ns, tx2 - tx1, imax, now() - tx2);
  *nfound = nf;
  *overflow = ov;
  free(S);
  free(Vs);
  free(a2);
  free(alive);
  return 0;
}

static int decode_fr(Result *R, int quiet) {
  int nl, nloc = 0;
  double t0 = now();
  memset(R, 0, sizeof *R);
  int nrows = feng_rao_reduce(0);
  if (nrows < 0) {
    R->failed = 1;
    return -1;
  }
  nl = count_basis_upto(L_DEG);
  for (int i = 0; i < nl; ++i)
    if (piv_of_row[i] < 0)
      ++nloc;
  for (int i = 0; i < nrows; ++i)
    if (piv_of_row[i] >= 0) {
      R->rowset[basis[i].rho] = 1;
      R->colset[basis[piv_of_row[i]].rho] = 1;
      ++R->T;
    }
  R->t_main = now() - t0;
  if (!quiet) {
    printf("# FR : %d error locators with rho <= %d\n", nloc, L_DEG);
    printf("# FR : staircase of %d elements -> estimated t = %d\n", R->T, R->T);
  }
  if (R->T == 0)
    return 0;
  if (nloc > 0)
    print_fr_locator_equations(nl, quiet);
  if (nloc == 0) {
    if (!quiet)
      puts("# FR : no error locator below the degree bound -> failed");
    R->failed = 1;
    return -1;
  }
  t0 = now();
  int nfound = 0, overflow = 0;
  if (fr_extract_simd(nl, R->pos, MAXT, &nfound, &overflow) != 0) {
    nfound = 0;
    overflow = 0;
    for (int idx = 0; idx < npts; ++idx) {
      int vanish = 1;
      for (int i = 0; i < nl && vanish; ++i) {
        if (piv_of_row[i] >= 0)
          continue;
        fe acc = 0;
        for (int k = 0; k <= i; ++k)
          if (coefm[i][k])
            acc ^= mul(coefm[i][k], basis[k].val[idx]);
        if (acc)
          vanish = 0;
      }
      if (vanish) {
        if (nfound < MAXT)
          R->pos[nfound++] = idx;
        else
          overflow = 1;
      }
    }
  }
  R->nfound = nfound;
  R->overflow = overflow;
  R->t_extract = now() - t0;
  if (quiet) {
    printf("# FR : %d positions found by direct locator evaluation "
           "(characteristic polynomial disabled)\n",
           nfound);
    report_fr_locator_zero_counts(nl, R);
  }
  resolve_and_store_fr_forney(R);
  return 0;
}

/* ================= BMS (Sakata) with value vectors and majority voting
 * ================= */
typedef struct {
  int b, d;
  fe *f;
} Gen;
static Gen gens[MAXGEN];
static int ngen;
/* FR-VOTE-only sparse locator supports.  The dense coefficient vectors remain
 * available for the common algebra, but the voting hot path visits only their
 * nonzero terms instead of scanning every lattice point. */
static int *frv_support_idx[MAXGEN];
static fe *frv_support_val[MAXGEN];
static int frv_support_n[MAXGEN];
static int frv_ancestor[SGMAX];
static int frv_support_ready;
static uint32_t *frv_vote_counts; /* ngen x (M_EXT-M_SYND), exact pair counts */
static int frv_vote_count_len;
static int frv_vote_conv_ready;
static fe *pivvec[SGMAX];
static fe pivd[SGMAX];
static uint8_t colused[SGMAX], is_piv[SGMAX];

static fe *vec_new(void) {
  fe *v = malloc(sizeof(fe) * npts);
  if (!v) {
    fprintf(stderr, "out of memory\n");
    exit(1);
  }
  return v;
}
static fe disc(const fe *f, int m) {
  const fe *ph = basis[idx_of_ord[m]].val;
  fe acc = 0;
  for (int n = 0; n < nz_cnt; ++n) {
    int k = nz_idx[n];
    acc ^= mul(recv[k], mul(f[k], ph[k]));
  }
  return acc;
}
static void vec_mul(fe *out, const fe *h, const fe *f) {
  for (int k = 0; k < npts; ++k)
    out[k] = mul(h[k], f[k]);
}
static void vec_axpy(fe *f, fe c, const fe *p) {
  for (int k = 0; k < npts; ++k)
    f[k] ^= mul(c, p[k]);
}

/* ---- majority voting machinery ---- */
static fe xhat[SGMAX];
static fe **Bw;
static int selp[MAXB];
static int nsel_dual;
static fe *Aug;
static int lgrow[MAXB + MAXGENUS];

static int build_dual(void) {
  int *perm = malloc(sizeof(int) * npts);
  const int NB = nbasis;
  if (Bw) {
    for (int v = 0; v < GENUS; ++v)
      free(Bw[v]);
    free(Bw);
    Bw = NULL;
  }
  Bw = calloc((size_t)GENUS, sizeof(*Bw));
  /* Keep the GENUS right-hand sides in Bw instead of appending them to Aug.
   * This changes Aug from NB x (NB+GENUS) to NB x NB and avoids a second
   * copy of the transformed right-hand sides. */
  Aug = calloc((size_t)NB * (size_t)NB, sizeof(fe));
  if (!Bw || !Aug) {
    free(perm);
    free(Aug);
    Aug = NULL;
    if (Bw) {
      free(Bw);
      Bw = NULL;
    }
    return -1;
  }
  for (int v = 0; v < GENUS; ++v) {
    Bw[v] = calloc((size_t)NB, sizeof(fe));
    if (!Bw[v]) {
      for (int z = 0; z < v; ++z)
        free(Bw[z]);
      free(Bw);
      Bw = NULL;
      free(perm);
      free(Aug);
      Aug = NULL;
      return -1;
    }
  }
  unsigned seed = 12345u;
  if (NB != M_SYND + 1 || basis[NB - GENUS].rho != M_SYND + 1) {
    fprintf(stderr,
            "dual basis: unexpected basis layout (need M >= conductor)\n");
    free(perm);
    free(Aug);
    Aug = NULL;
    for (int v = 0; v < GENUS; ++v)
      free(Bw[v]);
    free(Bw);
    Bw = NULL;
    return -1;
  }
  if (NB > npts) {
    fprintf(stderr, "dual basis: fewer points than basis functions\n");
    free(perm);
    free(Aug);
    Aug = NULL;
    for (int v = 0; v < GENUS; ++v)
      free(Bw[v]);
    free(Bw);
    Bw = NULL;
    return -1;
  }
  for (int attempt = 0; attempt < 8; ++attempt) {
    int singular = 0;
    for (int i = 0; i < npts; ++i)
      perm[i] = i;
    for (int i = 0; i < NB; ++i) {
      seed = seed * 1664525u + 1013904223u;
      int j = i + (int)((seed >> 8) % (unsigned)(npts - i));
      int t = perm[i];
      perm[i] = perm[j];
      perm[j] = t;
      selp[i] = perm[i];
    }
    nsel_dual = NB;
    for (int j = 0; j < NB; ++j) {
      for (int p = 0; p < NB; ++p)
        Aug[(size_t)j * (size_t)NB + (size_t)p] = basis[j].val[selp[p]];
      for (int v = 0; v < GENUS; ++v)
        Bw[v][j] = (j == NB - GENUS + v) ? 1 : 0;
    }
    for (int c = 0; c < NB && !singular; ++c) {
      int p = -1;
      for (int r = c; r < NB; ++r)
        if (Aug[(size_t)r * (size_t)NB + (size_t)c]) {
          p = r;
          break;
        }
      if (p < 0) {
        singular = 1;
        break;
      }
      if (p != c) {
        for (int k = 0; k < NB; ++k) {
          fe t = Aug[(size_t)p * (size_t)NB + (size_t)k];
          Aug[(size_t)p * (size_t)NB + (size_t)k] =
              Aug[(size_t)c * (size_t)NB + (size_t)k];
          Aug[(size_t)c * (size_t)NB + (size_t)k] = t;
        }
        for (int v = 0; v < GENUS; ++v) {
          fe t = Bw[v][p];
          Bw[v][p] = Bw[v][c];
          Bw[v][c] = t;
        }
      }
      fe iv = inv(Aug[(size_t)c * (size_t)NB + (size_t)c]);
      for (int k = c; k < NB; ++k)
        Aug[(size_t)c * (size_t)NB + (size_t)k] =
            mul(Aug[(size_t)c * (size_t)NB + (size_t)k], iv);
      for (int v = 0; v < GENUS; ++v)
        Bw[v][c] = mul(Bw[v][c], iv);
      for (int k = c + 1; k < NB; ++k)
        lgrow[k] = Aug[(size_t)c * (size_t)NB + (size_t)k]
                       ? (int)glog[Aug[(size_t)c * (size_t)NB + (size_t)k]]
                       : -1;
#pragma omp parallel for schedule(static)
      for (int r = c + 1; r < NB; ++r) {
        fe f = Aug[(size_t)r * (size_t)NB + (size_t)c];
        if (!f)
          continue;
        int lf = glog[f];
        fe *row = Aug + (size_t)r * (size_t)NB;
        for (int k = c + 1; k < NB; ++k)
          if (lgrow[k] >= 0)
            row[k] ^= gexp[lf + lgrow[k]];
        row[c] = 0;
        /* Each worker owns row r, so updating all RHS rows here is race-free
         * only if v is included in the same row operation. */
        for (int v = 0; v < GENUS; ++v)
          Bw[v][r] ^= mul(f, Bw[v][c]);
      }
    }
    if (singular)
      continue;
    for (int v = 0; v < GENUS; ++v) {
      /* Back substitution in place: Bw[v] already contains the transformed
       * right-hand side, and solved entries are overwritten from bottom up. */
      for (int i = NB - 1; i >= 0; --i) {
        fe acc = Bw[v][i];
        for (int k = i + 1; k < NB; ++k)
          if (Aug[(size_t)i * (size_t)NB + (size_t)k] && Bw[v][k])
            acc ^= mul(Aug[(size_t)i * (size_t)NB + (size_t)k], Bw[v][k]);
        Bw[v][i] = acc;
      }
    }
    free(perm);
    free(Aug);
    Aug = NULL;
    return 0;
  }
  fprintf(stderr, "dual basis: evaluation matrix singular\n");
  free(perm);
  free(Aug);
  Aug = NULL;
  return -1;
}

static int dual_selftest(void) {
  fe *coefv = malloc(sizeof(fe) * nbasis), *hs = malloc(sizeof(fe) * nbasis);
  for (int j = 0; j < nbasis; ++j)
    coefv[j] = (fe)(1 + rand() % FO);
  for (int p = 0; p < nsel_dual; ++p) {
    fe acc = 0;
    for (int j = 0; j < nbasis; ++j)
      acc ^= mul(coefv[j], basis[j].val[selp[p]]);
    hs[p] = acc;
  }
  int bad = 0;
  for (int v = 0; v < GENUS; ++v) {
    fe acc = 0;
    for (int p = 0; p < nsel_dual; ++p)
      acc ^= mul(Bw[v][p], hs[p]);
    if (acc != coefv[nbasis - GENUS + v])
      bad = 1;
  }
  free(coefv);
  free(hs);
  return bad ? -1 : 0;
}

/* branch-free GF multiply-accumulate in the log domain: zero is the log code
 * 2*FO, and the table E0 is 0 beyond 2*FO, so  E0[la + lb] = a*b  (0 if either
 * is zero).  LBw holds the logs of the dual rows Bw. */
static fe *E0;
static uint16_t **LBw;
static void init_dual_logs(void) {
  E0 = calloc((size_t)4 * FO + 8, sizeof(fe));
  for (int k = 0; k < 2 * FO; ++k)
    E0[k] = gexp[k];
  LBw = calloc((size_t)GENUS, sizeof *LBw);
  if (!LBw)
    return;
  for (int v = 0; v < GENUS; ++v) {
    LBw[v] = calloc((size_t)nsel_dual, sizeof(**LBw));
    if (!LBw[v])
      return;
    for (int p = 0; p < nsel_dual; ++p)
      LBw[v][p] = Bw[v][p] ? glog[Bw[v][p]] : (uint16_t)(2 * FO);
  }
}

/* ext_known results are cached per BMS step n (same function vector, same m):
 * the voting step and the discrepancy step ask for the same value. */
typedef struct {
  const fe *f;
  int m, n;
  fe known, an;
} ExtCache;
static ExtCache ecache[MAXGEN + 8];
static int necache, ecache_n = -1;

static fe ext_known(const fe *f, int m, int n, fe *an_out) {
  static fe hsel[MAXB];
  static uint16_t hl[MAXB];
  fe a[MAXGENUS + 1];
  const fe *ph = basis[idx_of_ord[m]].val;
  const fe *top[MAXGENUS + 1];
  const int nv = n - M_SYND;
  fe s = 0;
  if (ecache_n != n) {
    ecache_n = n;
    necache = 0;
  }
  for (int c = 0; c < necache; ++c)
    if (ecache[c].f == f && ecache[c].m == m) {
      *an_out = ecache[c].an;
      return ecache[c].known;
    }
  for (int p = 0; p < nsel_dual; ++p) {
    int k = selp[p];
    hsel[p] = mul(f[k], ph[k]);
    hl[p] = hsel[p] ? glog[hsel[p]] : (uint16_t)(2 * FO);
  }
  for (int v = 1; v <= nv; ++v) {
    const uint16_t *lb = LBw[v - 1];
    fe acc = 0;
    for (int p = 0; p < nsel_dual; ++p)
      acc ^= E0[hl[p] + lb[p]];
    a[v] = acc;
    top[v] = basis[idx_of_ord[M_SYND + v]].val;
  }
  for (int i = 0; i < nz_cnt; ++i) {
    int k = nz_idx[i];
    fe hv = mul(f[k], ph[k]);
    for (int v = 1; v <= nv; ++v)
      if (a[v])
        hv ^= mul(a[v], top[v][k]);
    s ^= mul(recv[k], hv);
  }
  for (int v = 1; v < nv; ++v)
    if (a[v])
      s ^= mul(a[v], xhat[M_SYND + v]);
  *an_out = a[nv];
  if (necache < MAXGEN + 8) {
    ecache[necache].f = f;
    ecache[necache].m = m;
    ecache[necache].n = n;
    ecache[necache].known = s;
    ecache[necache].an = a[nv];
    ++necache;
  }
  return s;
}
static fe xdisc(const fe *f, int m, int n) {
  fe an, known;
  if (n <= M_SYND)
    return disc(f, m);
  known = ext_known(f, m, n, &an);
  return (fe)(known ^ mul(an, xhat[n]));
}
static fe true_syndrome(int n) { return basis_syndrome(idx_of_ord[n]); }

/* ================================================================== */
/*  MONOMIAL SYNDROME TABLES (Newton polygon, NO dual basis)            */
/* ================================================================== */
/* All lattice points p=(i,j) of the cone {w_k.p>=0 (k!=e), pole order <= M+g}
 * are kept as monomials x^i y^j (there are only a few more than dim L((M+g)P)).
 * S(p) = sum_k e_k x_k^i y_k^j.  Pole orders <= M: S(p) is a parity check,
 * computed directly.  Above M: only the independent points ("reps", one per
 * pole order) are unknown (voted); a dependent point p = h + vA reduces via
 *      x^h * F = 0  :  S(h+vA) = S(h+vB) + sum_{other terms t} S(h+vt)
 * where vA, vB are the two terms of F on the pole edge (same pole order; the
 * other terms have strictly lower pole order).  So f is kept as a coefficient
 * vector over the lattice points: no evaluation, no dual basis, no elimination.
 */
static int use_poly = 0, use_poly_opt = 1;
static int NPT;
static int *pp_i, *pp_j, *pp_w, *pp_q, *pp_low, *pp_dep;
static int pp_nlow;
static int *pp_rep, *pp_ws, *pp_we, *pp_upto; /* by pole order 0..Mx */
static int pp_Mx;
static int pp_bi0, pp_bj0, pp_bw, pp_bh, *pp_box;
static fe *pp_S, *pp_lam, *pp_K;
static int *pivb;

static inline int pp_at(int i, int j) {
  i -= pp_bi0;
  j -= pp_bj0;
  if (i < 0 || i >= pp_bw || j < 0 || j >= pp_bh)
    return -1;
  return pp_box[(size_t)i * pp_bh + j];
}

/* The Newton path can evaluate a basis monomial directly at one point.  This
 * avoids retaining the full basis-by-position evaluation matrix when only
 * FR-VOTE is selected. */
static fe basis_value_at(int basis_index, int point_index) {
  if (use_poly && bas_ei && bas_ej)
    return mono_eval(PXv[point_index], PYv[point_index], bas_ei[basis_index],
                     bas_ej[basis_index]);
  return basis[basis_index].val[point_index];
}

static fe basis_syndrome(int basis_index) {
  if (use_poly && bas_ei && bas_ej) {
    int z = pp_at(bas_ei[basis_index], bas_ej[basis_index]);
    if (z >= 0)
      return pp_S[z];
  }
  fe acc = 0;
  for (int i = 0; i < nz_cnt; ++i) {
    int k = nz_idx[i];
    acc ^= mul(recv[k], basis_value_at(basis_index, k));
  }
  return acc;
}

/* Experimental decoder for a known common error value.
 * It solves the syndrome equations over GF(2), with one binary unknown per
 * evaluation position.  This is deliberately separate from FR/FR-VOTE.
 * The small-instance limits are intentional: minimum-weight binary syndrome
 * decoding is combinatorial, and this prototype enumerates the affine
 * solution space only when its dimension is at most 20. */
static int decode_common_value_positions(int common_value, int expected_weight,
                                         int verbose) {
  if (common_value <= 0 || common_value >= FS) {
    fprintf(stderr, "# COMMON-BINARY: requires a known nonzero common value\n");
    return -1;
  }
  const int n = npts;
  const int nb = count_basis_upto(M_SYND);
  const int nr = nb * FM;
  if (n > 4096 || nr > 8192) {
    printf("# COMMON-BINARY: prototype limit exceeded (n=%d, binary rows=%d; "
           "limits 4096/8192)\n",
           n, nr);
    return -1;
  }
  const int words = (n + 1 + 63) / 64;
  uint64_t *mat = calloc((size_t)nr * words, sizeof(uint64_t));
  int *pivot_col = malloc((size_t)(nr < n ? nr : n) * sizeof(int));
  if (!mat || !pivot_col) {
    free(mat);
    free(pivot_col);
    fprintf(stderr, "# COMMON-BINARY: allocation failed\n");
    return -1;
  }
  const fe ic = inv((fe)common_value);
  int row = 0;
  for (int bi = 0; bi < nb; ++bi) {
    if (basis[bi].rho > M_SYND)
      continue;
    fe rhs = mul(basis_syndrome(bi), ic);
    for (int bit = 0; bit < FM; ++bit, ++row) {
      uint64_t *r = mat + (size_t)row * words;
      for (int p = 0; p < n; ++p)
        if ((basis_value_at(bi, p) >> bit) & 1u)
          r[p >> 6] |= UINT64_C(1) << (p & 63);
      if ((rhs >> bit) & 1u)
        r[n >> 6] |= UINT64_C(1) << (n & 63);
    }
  }
  /* Reduced row-echelon form over GF(2). */
  int rank = 0;
  for (int col = 0; col < n && rank < row; ++col) {
    int sel = rank;
    while (sel < row &&
           !((mat[(size_t)sel * words + (col >> 6)] >> (col & 63)) & 1u))
      ++sel;
    if (sel == row)
      continue;
    if (sel != rank)
      for (int w = 0; w < words; ++w) {
        uint64_t t = mat[(size_t)rank * words + w];
        mat[(size_t)rank * words + w] = mat[(size_t)sel * words + w];
        mat[(size_t)sel * words + w] = t;
      }
    for (int rr = 0; rr < row; ++rr)
      if (rr != rank &&
          ((mat[(size_t)rr * words + (col >> 6)] >> (col & 63)) & 1u))
        for (int w = 0; w < words; ++w)
          mat[(size_t)rr * words + w] ^= mat[(size_t)rank * words + w];
    pivot_col[rank++] = col;
  }
  for (int rr = rank; rr < row; ++rr) {
    int any = 0;
    for (int w = 0; w < (n + 63) / 64; ++w)
      any |= mat[(size_t)rr * words + w] != 0;
    if (!any && ((mat[(size_t)rr * words + (n >> 6)] >> (n & 63)) & 1u)) {
      printf("# COMMON-BINARY: inconsistent equations for common value %d\n",
             common_value);
      free(mat);
      free(pivot_col);
      return 0;
    }
  }
  uint8_t *is_pivot = calloc((size_t)n, 1);
  int *free_col = malloc((size_t)(n - rank + 1) * sizeof(int));
  if (!is_pivot || !free_col) {
    free(mat);
    free(pivot_col);
    free(is_pivot);
    free(free_col);
    return -1;
  }
  for (int i = 0; i < rank; ++i)
    is_pivot[pivot_col[i]] = 1;
  int nf = 0;
  for (int c = 0; c < n; ++c)
    if (!is_pivot[c])
      free_col[nf++] = c;
  printf("# COMMON-BINARY: GF(2) system %d x %d, rank=%d, nullity=%d\n", row, n,
         rank, nf);
  if (nf > 20) {
    printf("# COMMON-BINARY: nullity %d exceeds exact-search cap 20; no "
           "position claim\n",
           nf);
    free(mat);
    free(pivot_col);
    free(is_pivot);
    free(free_col);
    return 0;
  }
  /* Diagnostic: exact minimum distance of the binary kernel. */
  int kernel_dmin = n + 1;
  uint8_t *kx = calloc((size_t)n, 1);
  if (!kx) {
    free(mat);
    free(pivot_col);
    free(is_pivot);
    free(free_col);
    return -1;
  }
  for (uint64_t km = 1; km < (UINT64_C(1) << nf); ++km) {
    memset(kx, 0, (size_t)n);
    for (int j = 0; j < nf; ++j)
      kx[free_col[j]] = (uint8_t)((km >> j) & 1u);
    for (int i = 0; i < rank; ++i) {
      const uint64_t *r = mat + (size_t)i * words;
      int v = 0;
      for (int j = 0; j < nf; ++j)
        if (kx[free_col[j]] &&
            ((r[free_col[j] >> 6] >> (free_col[j] & 63)) & 1u))
          v ^= 1;
      kx[pivot_col[i]] = (uint8_t)v;
    }
    int kw = 0;
    for (int p = 0; p < n; ++p)
      kw += kx[p];
    if (kw < kernel_dmin)
      kernel_dmin = kw;
  }
  printf("# COMMON-BINARY kernel minimum distance: %d; guaranteed unique "
         "radius: %d\n",
         kernel_dmin, (kernel_dmin - 1) / 2);
  free(kx);
  uint8_t *x = calloc((size_t)n, 1);
  uint8_t *bestx = calloc((size_t)n, 1);
  int bestwt = n + 1, best_count = 0;
  uint64_t total = UINT64_C(1) << nf;
  for (uint64_t mask = 0; mask < total; ++mask) {
    memset(x, 0, (size_t)n);
    for (int j = 0; j < nf; ++j)
      x[free_col[j]] = (uint8_t)((mask >> j) & 1u);
    for (int i = 0; i < rank; ++i) {
      const uint64_t *r = mat + (size_t)i * words;
      int v = (int)((r[n >> 6] >> (n & 63)) & 1u);
      for (int j = 0; j < nf; ++j)
        if (x[free_col[j]] &&
            ((r[free_col[j] >> 6] >> (free_col[j] & 63)) & 1u))
          v ^= 1;
      x[pivot_col[i]] = (uint8_t)v;
    }
    int wt = 0;
    for (int p = 0; p < n; ++p)
      wt += x[p];
    /* Certify the reconstructed support against every used field syndrome.
     * This remains valid when recv also contains an unknown codeword. */
    int ok = 1;
    for (int bi = 0; bi < nb && ok; ++bi)
      if (basis[bi].rho <= M_SYND) {
        fe acc = 0;
        for (int p = 0; p < n; ++p)
          if (x[p])
            acc ^= basis_value_at(bi, p);
        if (acc != mul(basis_syndrome(bi), ic))
          ok = 0;
      }
    if (ok) {
      if (wt < bestwt) {
        bestwt = wt;
        best_count = 1;
        memcpy(bestx, x, (size_t)n);
      } else if (wt == bestwt) {
        ++best_count;
      }
    }
  }
  if (best_count == 1) {
    if (verbose) {
      printf("# COMMON-BINARY positions (zero-based):");
      for (int p = 0; p < n; ++p)
        if (bestx[p])
          printf(" %d", p);
      putchar('\n');
    }
    int same = 0;
    for (int e = 0; e < expected_weight; ++e)
      if (bestx[true_pos[e]])
        ++same;
    printf("# COMMON-BINARY minimum-weight support: %d; injected-position "
           "comparison: %d/%d\n",
           bestwt, same, expected_weight);
  }
  printf("# COMMON-BINARY: %s (searched %llu affine solutions)\n",
         best_count == 1  ? "unique minimum-weight position set certified"
         : best_count > 1 ? "ambiguous minimum-weight position sets"
                          : "no certified solution found",
         (unsigned long long)total);
  free(mat);
  free(pivot_col);
  free(is_pivot);
  free(free_col);
  free(x);
  free(bestx);
  return best_count == 1 ? 1 : 0;
}
typedef struct {
  int i, j, w, depth;
} PPt;
static int pp_cmp(const void *a, const void *b) {
  const PPt *x = a, *y = b;
  if (x->w != y->w)
    return x->w - y->w;
  if (x->depth != y->depth)
    return x->depth - y->depth;
  return x->i - y->i;
}

static int poly_force_regular = 0;
static int poly_setup_impl(int e, int Mx, NBas *B);
/* weights of the single place at infinity (both x and y have poles there) */
static int single_inf_weights(int e, int *wx, int *wy) {
  int ninf = 0;
  for (int k = 0; k < ne; ++k)
    if (ed[k].inf)
      ninf++;
  if (ninf != 1 || !ed[e].inf || ed[e].len != 1 || ed[e].wx >= 0 ||
      ed[e].wy >= 0)
    return 0;
  *wx = -ed[e].wx;
  *wy = -ed[e].wy;
  return 1;
}
static int poly_setup(int e, int Mx, NBas *B) {
  int rc = poly_setup_impl(e, Mx, B);
  int wx, wy;
  if (rc != 0 && hermitian_q == 0 && single_inf_weights(e, &wx, &wy)) {
    poly_force_regular = 1;
    rc = poly_setup_impl(e, Mx, B);
    poly_force_regular = 0;
  }
  return rc;
}
static int poly_setup_impl(int e, int Mx, NBas *B) {
  /* 1. all lattice points of the cone */
  long ca[MAXEDGE + 1], cb[MAXEDGE + 1], ck[MAXEDGE + 1];
  int nc = 0;
  for (int k = 0; k < ne; ++k)
    if (k != e) {
      ca[nc] = ed[k].wx;
      cb[nc] = ed[k].wy;
      ck[nc] = 0;
      nc++;
    }
  ca[nc] = ed[e].wx;
  cb[nc] = ed[e].wy;
  ck[nc] = Mx;
  nc++;
  size_t cap = 1 << 14, n0 = 0;
  PPt *P = malloc(cap * sizeof *P);
  int bad = 0;
  if (hermitian_q > 0 || poly_force_regular) {
    /*
     * Hermitian affine coordinate ring:
     *   y^q + y + x^(q+1) = 0.
     * Its regular one-point monomials are x^i y^j (i >= 0, j >= 0),
     * and y^q is reduced by y^q = y + x^(q+1).  The generic toric cone
     * also contains Laurent monomials (for example x^i/y), which are not
     * regular on the affine curve.  Enumerate the actual polynomial
     * semigroup instead of that larger toric cone.
     */
    int wx = hermitian_q, wy = hermitian_q + 1;
    if (poly_force_regular)
      single_inf_weights(e, &wx, &wy);
    for (int i = 0; i <= Mx / wx; ++i) {
      int jmax = (Mx - wx * i) / wy;
      for (int j = 0; j <= jmax; ++j) {
        if (n0 == cap) {
          cap *= 2;
          PPt *tmp = realloc(P, cap * sizeof *P);
          if (!tmp) {
            free(P);
            return -1;
          }
          P = tmp;
        }
        P[n0].i = i;
        P[n0].j = j;
        P[n0].w = wx * i + wy * j;
        P[n0].depth = 0;
        ++n0;
      }
    }
  } else {
    for (int dir = 0; dir < 2 && !bad; ++dir)
      for (int i = dir ? -1 : 0;; i += dir ? -1 : 1) {
        if (i > (1 << 20) || i < -(1 << 20)) {
          bad = 1;
          break;
        }
        long jl = -(1L << 40), jh = (1L << 40);
        double rl = -1e18, rh = 1e18;
        int empty = 0;
        for (int c = 0; c < nc && !empty; ++c) {
          long v = ca[c] * i + ck[c];
          if (cb[c] == 0) {
            if (v < 0)
              empty = 1;
          } else if (cb[c] > 0) {
            long l = -nb_fdiv(v, cb[c]);
            if (l > jl)
              jl = l;
            double r = -(double)v / (double)cb[c];
            if (r > rl)
              rl = r;
          } else {
            long h = nb_fdiv(v, -cb[c]);
            if (h < jh)
              jh = h;
            double r = (double)v / (double)(-cb[c]);
            if (r < rh)
              rh = r;
          }
        }
        if (empty || rl > rh + 1e-9)
          break;
        for (long j = jl; j <= jh; ++j) {
          int w = -(int)(ed[e].wx * (long)i + ed[e].wy * j);
          if (w < 0 || w > Mx) {
            bad = 1;
            break;
          }
          if (n0 == cap) {
            cap *= 2;
            PPt *tmp = realloc(P, cap * sizeof *P);
            if (!tmp) {
              free(P);
              return -1;
            }
            P = tmp;
          }
          P[n0].i = i;
          P[n0].j = (int)j;
          P[n0].w = w;
          P[n0].depth = 0;
          n0++;
        }
        if (bad)
          break;
      }
  }
  if (bad || n0 == 0) {
    free(P);
    return -1;
  }
  /* 2. the two terms of F on the pole edge */
  long vmin = 1L << 40;
  for (int s = 0; s < T; ++s) {
    long v = ed[e].wx * (long)TX[s] + ed[e].wy * (long)TY[s];
    if (v < vmin)
      vmin = v;
  }
  int tA = -1, tB = -1;
  for (int s = 0; s < T; ++s)
    if (ed[e].wx * (long)TX[s] + ed[e].wy * (long)TY[s] == vmin) {
      if (tA < 0)
        tA = s;
      else if (tB < 0)
        tB = s;
      else
        tA = -2;
    }
  if (tA < 0 || tB < 0 || T > MAXTERM) {
    free(P);
    return -1;
  }
  const int dx = TX[tB] - TX[tA], dy = TY[tB] - TY[tA];
  /* 3. box index (unsorted) */
  int imin = 1 << 30, imax = -(1 << 30), jmin = 1 << 30, jmax = -(1 << 30);
  for (size_t p = 0; p < n0; ++p) {
    if (P[p].i < imin)
      imin = P[p].i;
    if (P[p].i > imax)
      imax = P[p].i;
    if (P[p].j < jmin)
      jmin = P[p].j;
    if (P[p].j > jmax)
      jmax = P[p].j;
  }
  pp_bi0 = imin;
  pp_bj0 = jmin;
  pp_bw = imax - imin + 1;
  pp_bh = jmax - jmin + 1;
  if ((double)pp_bw * pp_bh > 2e8) {
    free(P);
    return -1;
  }
  pp_box = malloc(sizeof(int) * (size_t)pp_bw * pp_bh);
  for (size_t x = 0; x < (size_t)pp_bw * pp_bh; ++x)
    pp_box[x] = -1;
  for (size_t p = 0; p < n0; ++p)
    pp_box[(size_t)(P[p].i - imin) * pp_bh + (P[p].j - jmin)] = (int)p;
  /* 4. dependent points and chain depth */
  char *dep = calloc(n0, 1);
  int *nb = malloc(sizeof(int) * n0);
  for (size_t p = 0; p < n0; ++p) {
    int hi = P[p].i - TX[tA], hj = P[p].j - TY[tA], ok = 1;
    for (int s = 0; s < T && ok; ++s)
      if (pp_at(hi + TX[s], hj + TY[s]) < 0)
        ok = 0;
    dep[p] = ok;
    nb[p] = ok ? pp_at(P[p].i + dx, P[p].j + dy) : -1;
    P[p].depth = ok ? -1 : 0;
  }
  for (int changed = 1; changed;) {
    changed = 0;
    for (size_t p = 0; p < n0; ++p)
      if (dep[p] && P[p].depth < 0 && nb[p] >= 0 && P[nb[p]].depth >= 0) {
        P[p].depth = P[nb[p]].depth + 1;
        changed = 1;
      }
  }
  int ok = 1;
  for (size_t p = 0; p < n0; ++p)
    if (P[p].depth < 0)
      ok = 0;
  /* 5. exactly one independent point per pole order, M+1-g in total */
  int *cnt = calloc(Mx + 1, sizeof(int)), nrep = 0;
  for (size_t p = 0; p < n0 && ok; ++p)
    if (!dep[p]) {
      if (++cnt[P[p].w] > 1)
        ok = 0;
      nrep++;
    }
  if (nrep != Mx + 1 - genus_g)
    ok = 0;
  free(cnt);
  if (ok) { /* axis points: negative exponents need x,y != 0 positions */
    int laurent = 0;
    for (size_t p = 0; p < n0; ++p)
      if (P[p].i < 0 || P[p].j < 0)
        laurent = 1;
    if (laurent) {
      int tx0 = 0, ty0 = 0;
      for (int s = 0; s < T; ++s) {
        tx0 += (TX[s] == 0);
        ty0 += (TY[s] == 0);
      }
      if (tx0 > 1 || ty0 > 1)
        ok = 0;
    }
    B->laurent = laurent;
  }
  free(dep);
  free(nb);
  free(pp_box);
  if (!ok) {
    free(P);
    return -1;
  }
  /* 6. sort by (pole order, chain depth) and rebuild everything */
  qsort(P, n0, sizeof *P, pp_cmp);
  NPT = (int)n0;
  pp_i = malloc(sizeof(int) * n0);
  pp_j = malloc(sizeof(int) * n0);
  pp_w = malloc(sizeof(int) * n0);
  pp_q = malloc(sizeof(int) * n0);
  pp_dep = malloc(sizeof(int) * n0);
  pp_nlow = T - 2;
  pp_low = malloc(sizeof(int) * n0 * (pp_nlow > 0 ? pp_nlow : 1));
  pp_box = malloc(sizeof(int) * (size_t)pp_bw * pp_bh);
  for (size_t x = 0; x < (size_t)pp_bw * pp_bh; ++x)
    pp_box[x] = -1;
  for (size_t p = 0; p < n0; ++p) {
    pp_i[p] = P[p].i;
    pp_j[p] = P[p].j;
    pp_w[p] = P[p].w;
    pp_box[(size_t)(P[p].i - imin) * pp_bh + (P[p].j - jmin)] = (int)p;
  }
  pp_Mx = Mx;
  pp_rep = malloc(sizeof(int) * (Mx + 1));
  pp_ws = malloc(sizeof(int) * (Mx + 2));
  pp_we = malloc(sizeof(int) * (Mx + 2));
  pp_upto = malloc(sizeof(int) * (Mx + 2));
  for (int w = 0; w <= Mx; ++w) {
    pp_rep[w] = -1;
    pp_ws[w] = pp_we[w] = 0;
  }
  for (size_t p = 0; p < n0; ++p) {
    int hi = pp_i[p] - TX[tA], hj = pp_j[p] - TY[tA], isdep = 1;
    for (int s = 0; s < T; ++s)
      if (pp_at(hi + TX[s], hj + TY[s]) < 0)
        isdep = 0;
    pp_dep[p] = isdep;
    if (!isdep)
      pp_rep[pp_w[p]] = (int)p;
    else {
      pp_q[p] = pp_at(pp_i[p] + dx, pp_j[p] + dy);
      int c = 0;
      for (int s = 0; s < T; ++s)
        if (s != tA && s != tB)
          pp_low[p * pp_nlow + c++] = pp_at(hi + TX[s], hj + TY[s]);
    }
    if (pp_ws[pp_w[p]] == 0 && pp_we[pp_w[p]] == 0)
      pp_ws[pp_w[p]] = (int)p;
    pp_we[pp_w[p]] = (int)p + 1;
  }
  {
    int cum = 0;
    for (int w = 0; w <= Mx; ++w) {
      if (pp_we[w] > 0)
        cum = pp_we[w];
      pp_upto[w] = cum;
    }
  }
  /* 7. the basis = the independent points; B in pole-order order */
  for (int k = 0; k < B->n; ++k) {
    int r = pp_rep[B->rho[k]];
    if (r < 0) {
      free(P);
      return -1; /* cannot happen: same pole-order set */
    }
    B->ei[k] = pp_i[r];
    B->ej[k] = pp_j[r];
  }
  pp_S = calloc(n0, sizeof(fe));
  pp_lam = calloc(n0, sizeof(fe));
  pp_K = calloc(n0, sizeof(fe));
  pivb = calloc(SGMAX, sizeof(int));
  free(P);
  return 0;
}

/* -------------------------------------------------------------------------
 * Kötter / O'Sullivan one-point error-value update (generic rational edge).
 * Functions are represented in the Apéry module basis x^a z_b.  The
 * multiplication table is derived from the same Newton-cone reductions used
 * by the syndrome engine; no interpolation matrix or per-error solve occurs.
 * ------------------------------------------------------------------------- */
typedef struct {
  int a, b;
  fe v;
} KTerm;
typedef struct {
  KTerm *t;
  int n;
} KPoly;
static int kotter_xmul_ready = 0;
static int kotter_kappa = 0, kotter_amax = 0, kotter_maxw = 0;
static int kotter_axis = 0; /* 0:Y^-1, 1:X, 2:X^-1, 3:Y */
static int kotter_lam[128], kotter_zi[128], kotter_zj[128];
static KPoly kotter_E[12000];
static KPoly kotter_C[128][128];
static fe *kotter_f = NULL, *kotter_g = NULL;
static int kotter_stride_a, kotter_stride_b;

static size_t kidx(int i, int a, int b) {
  return ((size_t)i * (size_t)kotter_stride_a + (size_t)a) *
             (size_t)kotter_stride_b +
         (size_t)b;
}
static int kpoly_from_dense(KPoly *p, const fe *d, int na, int k) {
  int n = 0;
  for (int a = 0; a < na; ++a)
    for (int b = 0; b < k; ++b)
      n += d[(size_t)a * k + b] != 0;
  p->t = n ? malloc((size_t)n * sizeof(*p->t)) : NULL;
  if (n && !p->t)
    return 0;
  p->n = n;
  int z = 0;
  for (int a = 0; a < na; ++a)
    for (int b = 0; b < k; ++b) {
      fe v = d[(size_t)a * k + b];
      if (v)
        p->t[z++] = (KTerm){a, b, v};
    }
  return 1;
}
static int kotter_nf(int p, uint8_t *mark, int *list, int *nout) {
  if (p < 0 || p >= NPT)
    return 0;
  memset(mark, 0, (size_t)NPT);
  mark[p] = 1;
  for (int r = p; r >= 0; --r)
    if (mark[r] && pp_dep[r]) {
      mark[r] = 0;
      int q = pp_q[r];
      if (q < 0 || q >= r)
        return 0;
      mark[q] ^= 1;
      for (int c = 0; c < pp_nlow; ++c) {
        int v = pp_low[(size_t)r * pp_nlow + c];
        if (v < 0 || v >= r)
          return 0;
        mark[v] ^= 1;
      }
    }
  int n = 0;
  for (int r = 0; r < NPT; ++r)
    if (mark[r]) {
      if (pp_dep[r])
        return 0;
      list[n++] = r;
    }
  *nout = n;
  return 1;
}
static void kotter_toggle_poly(fe *dense, const KPoly *p) {
  for (int i = 0; i < p->n; ++i) {
    const KTerm *t = &p->t[i];
    dense[(size_t)t->a * kotter_kappa + t->b] ^= t->v;
  }
}
static int kotter_prepare_module(int edge, NBas *B) {
  if (!use_poly || !B || edge < 0 || edge >= ne)
    return 0;
  /* Generic one-point Apéry module. Choose the inverse coordinate with
   * larger pole order as x; this keeps the module rank kappa manageable.
   * The curve polynomial is read from TX/TY/TC, not matched to fixed terms. */
  if (ed[edge].len != 1 || (ed[edge].wx == 0 && ed[edge].wy == 0))
    return 0;
  if (abs(ed[edge].wx) > abs(ed[edge].wy))
    kotter_axis = ed[edge].wx < 0 ? 1 : 2;
  else
    kotter_axis = ed[edge].wy > 0 ? 0 : 3;
  kotter_kappa = (kotter_axis == 1 || kotter_axis == 2) ? abs(ed[edge].wx)
                                                        : abs(ed[edge].wy);
  if (kotter_kappa < 1 || kotter_kappa > 128 ||
      B->n <= kotter_kappa) /* EXTENSION_POINT: module-rank cap */
    return 0;
  for (int b = 0; b < kotter_kappa; ++b) {
    int h = 0;
    while (h < B->n && B->rho[h] % kotter_kappa != b)
      ++h;
    if (h == B->n)
      return 0;
    kotter_lam[b] = B->rho[h];
    kotter_zi[b] = B->ei[h];
    kotter_zj[b] = B->ej[h];
  }
  int lmax = 0;
  for (int b = 0; b < kotter_kappa; ++b)
    if (kotter_lam[b] > lmax)
      lmax = kotter_lam[b];
  kotter_maxw = 2 * lmax;
  if (kotter_maxw >= (int)(sizeof(kotter_E) / sizeof(kotter_E[0])) ||
      kotter_maxw > pp_Mx)
    return 0;
  kotter_amax = M_SYND / kotter_kappa + 8;
  if (kotter_amax < kotter_maxw / kotter_kappa + 4)
    kotter_amax = kotter_maxw / kotter_kappa + 4;
  size_t dn = (size_t)(kotter_amax + 1) * kotter_kappa;
  fe *dense = calloc(dn, sizeof(fe));
  uint8_t *mark = malloc((size_t)NPT);
  int *list = malloc((size_t)NPT * sizeof(int));
  if (!dense || !mark || !list) {
    free(dense);
    free(mark);
    free(list);
    return 0;
  }

  /* E[w] converts the canonical Newton monomial of pole order w into the
   * Apéry module basis x^a z_b.  The recursion is triangular in pole order. */
  for (int w = 0; w <= kotter_maxw; ++w) {
    int pcan = pp_rep[w];
    if (pcan < 0)
      continue;
    int b = w % kotter_kappa;
    int a = (w - kotter_lam[b]) / kotter_kappa;
    if (a < 0) {
      free(dense);
      free(mark);
      free(list);
      return 0;
    }
    int I = kotter_zi[b], J = kotter_zj[b];
    if (kotter_axis == 0)
      J -= a;
    else if (kotter_axis == 1)
      I += a;
    else if (kotter_axis == 2)
      I -= a;
    else
      J += a;
    int p = pp_at(I, J), nn = 0;
    if (!kotter_nf(p, mark, list, &nn)) {
      free(dense);
      free(mark);
      free(list);
      return 0;
    }
    memset(dense, 0, dn * sizeof(fe));
    dense[(size_t)a * kotter_kappa + b] = 1;
    int saw_top = 0;
    for (int q = 0; q < nn; ++q) {
      int r = list[q], wr = pp_w[r];
      if (wr == w) {
        if (r != pcan) {
          free(dense);
          free(mark);
          free(list);
          return 0;
        }
        saw_top = 1;
      } else if (wr < w)
        kotter_toggle_poly(dense, &kotter_E[wr]);
      else {
        free(dense);
        free(mark);
        free(list);
        return 0;
      }
    }
    if (!saw_top ||
        !kpoly_from_dense(&kotter_E[w], dense, kotter_amax + 1, kotter_kappa)) {
      free(dense);
      free(mark);
      free(list);
      return 0;
    }
  }

  /* Multiplication table z_j*z_c in the Apéry basis. */
  for (int j = 0; j < kotter_kappa; ++j)
    for (int c = 0; c < kotter_kappa; ++c) {
      int p = pp_at(kotter_zi[j] + kotter_zi[c], kotter_zj[j] + kotter_zj[c]);
      int nn = 0;
      if (!kotter_nf(p, mark, list, &nn)) {
        free(dense);
        free(mark);
        free(list);
        return 0;
      }
      memset(dense, 0, dn * sizeof(fe));
      for (int q = 0; q < nn; ++q)
        kotter_toggle_poly(dense, &kotter_E[pp_w[list[q]]]);
      if (!kpoly_from_dense(&kotter_C[j][c], dense, kotter_amax + 1,
                            kotter_kappa)) {
        free(dense);
        free(mark);
        free(list);
        return 0;
      }
    }
  free(dense);
  free(mark);
  free(list);
  kotter_stride_a = kotter_amax + 1;
  kotter_stride_b = kotter_kappa;
  size_t fn = (size_t)kotter_kappa * kotter_stride_a * kotter_stride_b;
  kotter_f = calloc(fn, sizeof(fe));
  kotter_g = calloc(fn, sizeof(fe));
  if (!kotter_f || !kotter_g)
    return 0;
  for (int i = 0; i < kotter_kappa; ++i)
    kotter_f[kidx(i, 0, i)] = 1;
  kotter_xmul_ready = 1;
  return 1;
}
static int kotter_rho(const fe *f, int i) {
  int d = -1;
  for (int a = 0; a <= kotter_amax; ++a)
    for (int b = 0; b < kotter_kappa; ++b)
      if (f[kidx(i, a, b)]) {
        int w = a * kotter_kappa + kotter_lam[b];
        if (w > d)
          d = w;
      }
  return d;
}
static fe kotter_syndrome(int a, int b) {
  if (a < 0)
    return 0;
  int I = kotter_zi[b], J = kotter_zj[b];
  if (kotter_axis == 0)
    J -= a;
  else if (kotter_axis == 1)
    I += a;
  else if (kotter_axis == 2)
    I -= a;
  else
    J += a;
  int p = pp_at(I, J);
  if (p < 0 || pp_w[p] > M_SYND)
    return 0;
  return pp_S[p];
}
static fe kotter_mu(const fe *f, int i, int j, int r) {
  fe mu = 0;
  for (int a = 0; a <= kotter_amax; ++a)
    for (int c = 0; c < kotter_kappa; ++c) {
      fe fc = f[kidx(i, a, c)];
      if (!fc)
        continue;
      const KPoly *p = &kotter_C[j][c];
      for (int z = 0; z < p->n; ++z) {
        const KTerm *t = &p->t[z];
        fe s = kotter_syndrome(a + t->a + r, t->b);
        if (s)
          mu ^= mul(fc, mul(t->v, s));
      }
    }
  return mu;
}
static fe kotter_partial_curve(fe X, fe Y, int dx) {
  fe acc = 0;
  for (int t = 0; t < T; ++t) {
    int e = dx ? TX[t] : TY[t];
    if (e & 1) {
      int I = TX[t] - (dx ? 1 : 0);
      int J = TY[t] - (dx ? 0 : 1);
      acc ^= mul(TC[t], mono_eval(X, Y, I, J));
    }
  }
  return acc;
}

/* Differentiate a function on F(X,Y)=0 with respect to the selected
 * pole coordinate x (one of X, X^-1, Y, Y^-1). */
static int kotter_eval_derivative(const fe *f, int i, fe X, fe Y, fe *out) {
  if (!X || !Y)
    return 0;
  fe FX = kotter_partial_curve(X, Y, 1);
  fe FY = kotter_partial_curve(X, Y, 0);
  fe X2 = mul(X, X), Y2 = mul(Y, Y);
  fe dXdx = 0, dYdx = 0;
  switch (kotter_axis) {
  case 0: /* x=Y^-1 */
    if (!FX)
      return 0;
    dXdx = mul(mul(FY, Y2), inv(FX));
    dYdx = Y2;
    break;
  case 1: /* x=X */
    if (!FY)
      return 0;
    dXdx = 1;
    dYdx = mul(FX, inv(FY));
    break;
  case 2: /* x=X^-1 */
    if (!FY)
      return 0;
    dXdx = X2;
    dYdx = mul(mul(FX, X2), inv(FY));
    break;
  default: /* x=Y */
    if (!FX)
      return 0;
    dXdx = mul(FY, inv(FX));
    dYdx = 1;
    break;
  }
  fe acc = 0;
  for (int a = 0; a <= kotter_amax; ++a)
    for (int b = 0; b < kotter_kappa; ++b) {
      fe c = f[kidx(i, a, b)];
      if (!c)
        continue;
      int I = kotter_zi[b], J = kotter_zj[b];
      if (kotter_axis == 0)
        J -= a;
      else if (kotter_axis == 1)
        I += a;
      else if (kotter_axis == 2)
        I -= a;
      else
        J += a;
      if (I & 1)
        acc ^= mul(c, mul(mono_eval(X, Y, I - 1, J), dXdx));
      if (J & 1)
        acc ^= mul(c, mul(mono_eval(X, Y, I, J - 1), dYdx));
    }
  *out = acc;
  return 1;
}

static int resolve_by_kotter(int *found, int *nfound) {
  if (!kotter_xmul_ready || !found || !nfound || *nfound <= 0 || *nfound > MAXT)
    return 0;
  const int k = kotter_kappa, A = kotter_amax,
            stride = kotter_stride_a * kotter_stride_b;
  memset(kotter_f, 0, (size_t)k * stride * sizeof(fe));
  memset(kotter_g, 0, (size_t)k * stride * sizeof(fe));
  for (int i = 0; i < k; ++i)
    kotter_f[kidx(i, 0, i)] = 1;
  fe *oldf = malloc((size_t)stride * sizeof(fe));
  if (!oldf)
    return 0;
  fe *oldf2 = malloc((size_t)stride * sizeof(fe));
  if (!oldf2) {
    free(oldf);
    return 0;
  }
  for (int m = 0; m <= M_SYND; ++m) {
    for (int i = 0; i < k; ++i) {
      int j = (m - i) % k;
      if (j < 0)
        j += k;
      /* The paper's "pair i,j" is unordered.  Process {i,j} once,
       * applying the symmetric update to (f_i,g_j) and (f_j,g_i). */
      if (i > j)
        continue;
      int di = kotter_rho(kotter_f, i), dj = kotter_rho(kotter_f, j);
      if (di < 0 || dj < 0 || (m - di - kotter_lam[j]) % k ||
          (m - dj - kotter_lam[i]) % k) {
        free(oldf);
        free(oldf2);
        return 0;
      }
      int ri = (m - di - kotter_lam[j]) / k;
      int rj = (m - dj - kotter_lam[i]) / k;
      fe mui = kotter_mu(kotter_f, i, j, ri);
      fe muj = (i == j) ? mui : kotter_mu(kotter_f, j, i, rj);
      if (mui != muj) {
        free(oldf);
        free(oldf2);
        return 0;
      }
      int p = (di + dj - m) / k - 1;
      fe *fi = kotter_f + (size_t)i * stride;
      fe *gj = kotter_g + (size_t)j * stride;
      fe *fj = kotter_f + (size_t)j * stride;
      fe *gi = kotter_g + (size_t)i * stride;
      if (!mui || p >= 0) {
        if (mui) {
          for (int a = 0; a <= A; ++a) {
            int src = a - p;
            if (src < 0)
              continue;
            for (int b = 0; b < k; ++b) {
              fe v = gj[(size_t)src * k + b];
              if (v)
                fi[(size_t)a * k + b] ^= mul(mui, v);
              if (i != j) {
                v = gi[(size_t)src * k + b];
                if (v)
                  fj[(size_t)a * k + b] ^= mul(mui, v);
              }
            }
          }
        }
      } else {
        int sh = -p;
        memcpy(oldf, fi, (size_t)stride * sizeof(fe));
        if (i != j)
          memcpy(oldf2, fj, (size_t)stride * sizeof(fe));
        fe im = inv(mui);
        for (int a = 0; a <= A; ++a)
          for (int b = 0; b < k; ++b) {
            size_t z = (size_t)a * k + b;
            fe v = (a >= sh) ? oldf[(size_t)(a - sh) * k + b] : 0;
            if (gj[z])
              v ^= mul(mui, gj[z]);
            fi[z] = v;
            if (i != j) {
              v = (a >= sh) ? oldf2[(size_t)(a - sh) * k + b] : 0;
              if (gi[z])
                v ^= mul(mui, gi[z]);
              fj[z] = v;
            }
          }
        for (int a = 0; a <= A; ++a)
          for (int b = 0; b < k; ++b) {
            size_t z = (size_t)a * k + b;
            gj[z] = mul(im, oldf[z]);
            if (i != j)
              gi[z] = mul(im, oldf2[z]);
          }
      }
    }
  }
  free(oldf2);
  free(oldf);
  int n = *nfound;
  fe *vals = malloc((size_t)n * sizeof(fe));
  if (!vals)
    return 0;
  for (int u = 0; u < n; ++u) {
    fe X = PXv[found[u]], Y = PYv[found[u]], sum = 0;
    for (int i = 0; i < k; ++i) {
      fe fp = 0;
      if (!kotter_eval_derivative(kotter_f, i, X, Y, &fp)) {
        free(vals);
        return 0;
      }
      fe gp = 0;
      for (int a = 0; a <= kotter_amax; ++a)
        for (int b = 0; b < k; ++b) {
          fe c = kotter_g[kidx(i, a, b)];
          if (c) {
            int I = kotter_zi[b], J = kotter_zj[b];
            if (kotter_axis == 0)
              J -= a;
            else if (kotter_axis == 1)
              I += a;
            else if (kotter_axis == 2)
              I -= a;
            else
              J += a;
            gp ^= mul(c, mono_eval(X, Y, I, J));
          }
        }
      sum ^= mul(fp, gp);
    }
    if (!sum) {
      free(vals);
      return 0;
    }
    vals[u] = inv(sum);
  }
  /* Independent full syndrome verification. */
  int nb = count_basis_upto(M_SYND);
  for (int a = 0; a < nb; ++a) {
    fe check = 0;
    for (int u = 0; u < n; ++u)
      check ^= mul(vals[u], basis_value_at(a, found[u]));
    if (check != synd1(a)) {
      free(vals);
      return 0;
    }
  }
  for (int u = 0; u < n; ++u)
    err_val[u] = vals[u];
  free(vals);
  return 1;
}

/* syndromes of the received word: direct for pole order <= M, 0 above */
static void poly_syndromes(void) {
  for (int p = 0; p < NPT; ++p) {
    if (pp_w[p] > M_SYND) {
      pp_S[p] = 0;
      continue;
    }
    fe acc = 0;
    for (int z = 0; z < nz_cnt; ++z) {
      int k = nz_idx[z];
      acc ^= mul(recv[k], mono_eval(PXv[k], PYv[k], pp_i[p], pp_j[p]));
    }
    pp_S[p] = acc;
  }
}
/* Explicit BMS doubled-syndrome matrix diagnostic.
 * For basis functions f_i with rho(f_i) <= floor(M/2), form
 *     S2[i,j] = S(f_i f_j) = sum_k e_k f_i(P_k) f_j(P_k).
 * This is the finite square submatrix of the product-syndrome array
 * (the 2Sigma construction).  The normal BMS path accesses the same entries
 * implicitly through pp_shift()/poly_xdisc(); this routine materializes them
 * only for a small diagnostic run and checks them against pp_S.
 */
static int bms2_check_enabled = 0;
static int bms2_syndrome_check(void) {
  int ids[MAXB], r = 0, mismatch = 0, entries = 0;
  if (!use_poly || !pp_S || !bas_ei || !bas_ej) {
    puts("# BMS2: unavailable (Newton monomial syndrome table is not active)");
    return -1;
  }
  for (int bi = 0; bi < nbasis; ++bi)
    if (basis[bi].rho <= M_SYND / 2)
      ids[r++] = bi;
  if (r <= 0) {
    puts("# BMS2: empty product-syndrome matrix");
    return -1;
  }
  /* This explicit O(r^3) diagnostic is intentionally limited to small
   * instances; the production BMS path does not materialize this matrix. */
  if (r > 256) {
    printf("# BMS2: extraction skipped (matrix dimension %d exceeds diagnostic "
           "cap 256)\n",
           r);
    return -1;
  }
  fe *mat = calloc((size_t)r * (size_t)r, sizeof(fe));
  if (!mat) {
    puts("# BMS2: allocation failed");
    return -1;
  }
  for (int i = 0; i < r; ++i) {
    int bi = ids[i];
    for (int j = 0; j < r; ++j) {
      int bj = ids[j];
      fe v = 0;
      for (int z = 0; z < nz_cnt; ++z) {
        int k = nz_idx[z];
        fe vi = basis_value_at(bi, k);
        fe vj = basis_value_at(bj, k);
        v ^= mul(recv[k], mul(vi, vj));
      }
      mat[(size_t)i * r + j] = v;
      int t = pp_at(bas_ei[bi] + bas_ei[bj], bas_ej[bi] + bas_ej[bj]);
      if (t < 0 || pp_w[t] > M_SYND || pp_S[t] != v)
        ++mismatch;
      else
        ++entries;
    }
  }
  /* Reduced row-echelon form over GF(2^m).  We retain pivot columns so
   * that the full nullspace, not just the rank, is available below. */
  int *pivot_col = malloc((size_t)r * sizeof(*pivot_col));
  if (!pivot_col) {
    free(mat);
    puts("# BMS2: pivot allocation failed");
    return -1;
  }
  int rank = 0;
  for (int col = 0; col < r; ++col) {
    int piv = rank;
    while (piv < r && mat[(size_t)piv * r + col] == 0)
      ++piv;
    if (piv == r)
      continue;
    if (piv != rank)
      for (int j = 0; j < r; ++j) {
        fe tmp = mat[(size_t)rank * r + j];
        mat[(size_t)rank * r + j] = mat[(size_t)piv * r + j];
        mat[(size_t)piv * r + j] = tmp;
      }
    fe ip = inv(mat[(size_t)rank * r + col]);
    for (int j = 0; j < r; ++j)
      mat[(size_t)rank * r + j] = mul(mat[(size_t)rank * r + j], ip);
    for (int i = 0; i < r; ++i) {
      if (i == rank)
        continue;
      fe c = mat[(size_t)i * r + col];
      if (!c)
        continue;
      for (int j = 0; j < r; ++j)
        mat[(size_t)i * r + j] ^= mul(c, mat[(size_t)rank * r + j]);
    }
    pivot_col[rank++] = col;
  }
  int nullity = r - rank;
  printf("# BMS2: explicit product-syndrome matrix %d x %d, ", r, r);
  printf("rank=%d, nullity=%d, injected errors=%d\n", rank, nullity, nz_cnt);
  printf("# BMS2: checked entries=%d, mismatches=%d; entries are ", entries,
         mismatch);
  puts("S_ij = sum_e e_k f_i(P_k) f_j(P_k) = S(f_i f_j)");
  if (mismatch || nullity <= 0) {
    if (!nullity)
      puts("# BMS2: nullspace is trivial; cannot extract common zeros");
    free(pivot_col);
    free(mat);
    return mismatch ? -1 : rank;
  }
  if (nullity > 512 || r > 2048) {
    printf("# BMS2: extraction skipped (r=%d, nullity=%d; diagnostic cap "
           "r<=2048, nullity<=512)\n",
           r, nullity);
    free(pivot_col);
    free(mat);
    return rank;
  }
  uint8_t *is_pivot = calloc((size_t)r, 1);
  fe *ker = calloc((size_t)nullity * (size_t)r, sizeof(fe));
  if (!is_pivot || !ker) {
    free(is_pivot);
    free(ker);
    free(pivot_col);
    free(mat);
    puts("# BMS2: nullspace allocation failed");
    return -1;
  }
  for (int i = 0; i < rank; ++i)
    is_pivot[pivot_col[i]] = 1;
  int q = 0;
  for (int free_col = 0; free_col < r; ++free_col)
    if (!is_pivot[free_col]) {
      fe *v = ker + (size_t)q * r;
      v[free_col] = 1;
      for (int i = 0; i < rank; ++i)
        v[pivot_col[i]] = mat[(size_t)i * r + free_col];
      ++q;
    }
  /* A point is a common zero precisely when every nullspace function
   * evaluates to zero there.  For a full-row-rank error evaluation matrix,
   * these common zeros are the error support. */
  int found[MAXT + 1], nfound = 0, too_many = 0;
  for (int p = 0; p < npts; ++p) {
    int common_zero = 1;
    for (int h = 0; h < nullity && common_zero; ++h) {
      const fe *v = ker + (size_t)h * r;
      fe value = 0;
      for (int j = 0; j < r; ++j)
        if (v[j])
          value ^= mul(v[j], basis_value_at(ids[j], p));
      if (value)
        common_zero = 0;
    }
    if (common_zero) {
      if (nfound < MAXT)
        found[nfound++] = p;
      else
        too_many = 1;
    }
  }
  printf("# BMS2: nullspace common-zero candidates=%s%d\n",
         too_many ? ">=" : "", nfound);
  if (!too_many && nfound > 0 && resolve_by_values(found, &nfound)) {
    int same = 0;
    for (int i = 0; i < nfound; ++i)
      for (int j = 0; j < NUM_ERRORS; ++j)
        if (found[i] == true_pos[j] && err_val[i] == recv[true_pos[j]]) {
          ++same;
          break;
        }
    printf("# BMS2: value solve recovered=%d, exact position/value "
           "matches=%d/%d -> %s\n",
           nfound, same, NUM_ERRORS,
           (nfound == NUM_ERRORS && same == NUM_ERRORS) ? "SUCCESS"
                                                        : "MISMATCH");
  } else {
    puts("# BMS2: value solve failed (candidate set is not the exact error "
         "support)");
  }
  free(is_pivot);
  free(ker);
  free(pivot_col);
  free(mat);
  return mismatch ? -1 : rank;
}

/* relation check on the known range (dependent points <= M) */
static int poly_selfcheck(void) {
  int bad = 0, tested = 0;
  for (int p = 0; p < NPT; ++p)
    if (pp_dep[p] && pp_w[p] <= M_SYND) {
      fe v = pp_S[pp_q[p]];
      for (int c = 0; c < pp_nlow; ++c)
        v ^= pp_S[pp_low[p * pp_nlow + c]];
      tested++;
      if (v != pp_S[p])
        bad++;
    }
  (void)tested;
  return bad;
}
/* weight-n points: S(p) = lam*X_n + K, X_n the (unknown) syndrome of the rep */
static void poly_prepare(int n) {
  for (int p = pp_ws[n]; p < pp_we[n]; ++p) {
    if (!pp_dep[p]) {
      pp_lam[p] = 1;
      pp_K[p] = 0;
    } else {
      int q = pp_q[p];
      fe k = pp_K[q];
      for (int c = 0; c < pp_nlow; ++c)
        k ^= pp_S[pp_low[p * pp_nlow + c]];
      pp_lam[p] = pp_lam[q];
      pp_K[p] = k;
    }
  }
}
static void poly_finalize(int n, fe x) {
  for (int p = pp_ws[n]; p < pp_we[n]; ++p)
    pp_S[p] = (fe)(mul(pp_lam[p], x) ^ pp_K[p]);
}
static fe *vnew(void) {
  size_t n = (size_t)(npts > NPT ? npts : NPT);
  fe *v = calloc(n, sizeof(fe));
  if (!v) {
    fprintf(stderr, "out of memory\n");
    exit(1);
  }
  return v;
}
static inline int pp_shift(int p, int ri, int rj) {
  int t = pp_at(pp_i[p] + ri, pp_j[p] + rj);
  if (t < 0) {
    fprintf(stderr, "internal error: product leaves the lattice cone\n");
    exit(1);
  }
  return t;
}
/* out = x^(rep c) * f ,  f of pole order bsrc */
static void vmul_b(fe *out, int c, const fe *f, int bsrc) {
  if (!use_poly) {
    vec_mul(out, basis[idx_of_ord[c]].val, f);
    return;
  }
  int r = pp_rep[c], ri = pp_i[r], rj = pp_j[r];
  for (int p = 0; p < pp_upto[bsrc]; ++p)
    if (f[p])
      out[pp_shift(p, ri, rj)] = f[p];
}
static void vaxpy_b(fe *f, fe c, const fe *piv, int bpiv) {
  if (!use_poly) {
    vec_axpy(f, c, piv);
    return;
  }
  for (int p = 0; p < pp_upto[bpiv]; ++p)
    if (piv[p])
      f[p] ^= mul(c, piv[p]);
}
static fe poly_xdisc(const fe *f, int b, int m) {
  int r = pp_rep[m], ri = pp_i[r], rj = pp_j[r];
  fe acc = 0;
  for (int p = 0; p < pp_upto[b]; ++p)
    if (f[p])
      acc ^= mul(f[p], pp_S[pp_shift(p, ri, rj)]);
  return acc;
}
static fe poly_ext_known(const fe *f, int b, int m, int n, fe *an) {
  int r = pp_rep[m], ri = pp_i[r], rj = pp_j[r];
  fe known = 0, a = 0;
  for (int p = 0; p < pp_upto[b]; ++p)
    if (f[p]) {
      int t = pp_shift(p, ri, rj);
      if (pp_w[t] < n)
        known ^= mul(f[p], pp_S[t]);
      else {
        a ^= mul(f[p], pp_lam[t]);
        known ^= mul(f[p], pp_K[t]);
      }
    }
  *an = a;
  return known;
}
static fe ext_known_b(const fe *f, int b, int m, int n, fe *an) {
  return use_poly ? poly_ext_known(f, b, m, n, an) : ext_known(f, m, n, an);
}
static fe xdisc_b(const fe *f, int b, int m, int n) {
  if (use_poly)
    return poly_xdisc(f, b, m);
  return xdisc(f, m, n);
}
/* does the polynomial f (pole order b) not vanish at position k */
static int poly_nonzero(const fe *f, int b, int k) {
  fe acc = 0;
  for (int p = 0; p < pp_upto[b]; ++p)
    if (f[p])
      acc ^= mul(f[p], mono_eval(PXv[k], PYv[k], pp_i[p], pp_j[p]));
  return acc != 0;
}
/* Exact NTT convolution for the majority-vote candidate counts.
 * MOD is larger than M_EXT, so every coefficient is an exact integer count,
 * not merely a residue. */
#define FRV_NTT_MOD 998244353u
#define FRV_NTT_ROOT 3u
static uint32_t frv_mod_pow(uint32_t a, uint32_t e) {
  uint64_t r = 1, x = a;
  while (e) {
    if (e & 1u)
      r = (r * x) % FRV_NTT_MOD;
    x = (x * x) % FRV_NTT_MOD;
    e >>= 1;
  }
  return (uint32_t)r;
}
static void frv_ntt(uint32_t *a, int n, int invert) {
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j) {
      uint32_t t = a[i];
      a[i] = a[j];
      a[j] = t;
    }
  }
  for (int len = 2; len <= n; len <<= 1) {
    uint32_t wlen =
        frv_mod_pow(FRV_NTT_ROOT, (FRV_NTT_MOD - 1u) / (uint32_t)len);
    if (invert)
      wlen = frv_mod_pow(wlen, FRV_NTT_MOD - 2u);
    for (int i = 0; i < n; i += len) {
      uint64_t w = 1;
      int half = len >> 1;
      for (int j = 0; j < half; ++j) {
        uint32_t u = a[i + j];
        uint32_t v = (uint32_t)((uint64_t)a[i + j + half] * w % FRV_NTT_MOD);
        uint32_t x = u + v;
        if (x >= FRV_NTT_MOD)
          x -= FRV_NTT_MOD;
        uint32_t y = u >= v ? u - v : u + FRV_NTT_MOD - v;
        a[i + j] = x;
        a[i + j + half] = y;
        w = w * wlen % FRV_NTT_MOD;
      }
    }
  }
  if (invert) {
    uint32_t invn = frv_mod_pow((uint32_t)n, FRV_NTT_MOD - 2u);
    for (int i = 0; i < n; ++i)
      a[i] = (uint32_t)((uint64_t)a[i] * invn % FRV_NTT_MOD);
  }
}
static int frv_build_vote_convolutions(void) {
  const int L = M_EXT + 1;
  const int need = 2 * L - 1;
  int ntt_n = 1;
  while (ntt_n < need && ntt_n < (1 << 23))
    ntt_n <<= 1;
  if (ntt_n < need || (uint64_t)(M_EXT + 1) >= FRV_NTT_MOD)
    return -1;
  uint32_t *fb = calloc((size_t)ntt_n, sizeof(uint32_t));
  uint32_t *fa = calloc((size_t)ntt_n, sizeof(uint32_t));
  size_t count_len = (size_t)(M_EXT - M_SYND);
  uint32_t *counts = calloc((size_t)ngen * count_len, sizeof(uint32_t));
  if (!fb || !fa || !counts) {
    free(fb);
    free(fa);
    free(counts);
    return -1;
  }
  for (int j = 0; j <= M_EXT; ++j)
    if (in_H[j] && !colused[j])
      fb[j] = 1;
  frv_ntt(fb, ntt_n, 0);
  for (int g = 0; g < ngen; ++g) {
    memset(fa, 0, (size_t)ntt_n * sizeof(uint32_t));
    for (int i = 0; i <= M_EXT; ++i)
      if (frv_ancestor[i] == g)
        fa[i] = 1;
    frv_ntt(fa, ntt_n, 0);
    for (int k = 0; k < ntt_n; ++k)
      fa[k] = (uint32_t)((uint64_t)fa[k] * fb[k] % FRV_NTT_MOD);
    frv_ntt(fa, ntt_n, 1);
    uint32_t *dst = counts + (size_t)g * count_len;
    for (int n = M_SYND + 1; n <= M_EXT; ++n)
      dst[n - (M_SYND + 1)] = fa[n];
  }
  free(fb);
  free(fa);
  free(frv_vote_counts);
  frv_vote_counts = counts;
  frv_vote_count_len = (int)count_len;
  frv_vote_conv_ready = 1;
  return 0;
}

/* Build sparse supports once after the FR kernel has produced its locators.
 * The ancestor of a pole order i depends only on i and the locator ideal, not
 * on the syndrome order n; precomputing it removes the inner generator scan
 * from every majority-vote step. */
static int frv_build_sparse(void) {
  for (int g = 0; g < ngen; ++g) {
    int lim = pp_upto[gens[g].b], nnz = 0;
    for (int p = 0; p < lim; ++p)
      nnz += gens[g].f[p] != 0;
    frv_support_idx[g] = nnz ? malloc((size_t)nnz * sizeof(int)) : NULL;
    frv_support_val[g] = nnz ? malloc((size_t)nnz * sizeof(fe)) : NULL;
    if (nnz && (!frv_support_idx[g] || !frv_support_val[g]))
      return -1;
    frv_support_n[g] = nnz;
    int z = 0;
    for (int p = 0; p < lim; ++p)
      if (gens[g].f[p]) {
        frv_support_idx[g][z] = p;
        frv_support_val[g][z] = gens[g].f[p];
        ++z;
      }
  }
  for (int i = 0; i <= M_EXT; ++i) {
    frv_ancestor[i] = -1;
    if (!in_H[i] || is_piv[i])
      continue;
    for (int g = 0; g < ngen; ++g) {
      int a = i - gens[g].b;
      if (a >= 0 && in_H[a]) {
        frv_ancestor[i] = g;
        break;
      }
    }
  }
  frv_support_ready = 1;
  /* The convolution path replaces O(M_EXT^2) candidate scans by a small
   * number of O(N log N) exact convolutions.  If it cannot be allocated, the
   * original exact scan remains available. */
  (void)frv_build_vote_convolutions();
  return 0;
}
static fe poly_ext_known_sparse(int g, int m, int n, fe *an) {
  int r = pp_rep[m], ri = pp_i[r], rj = pp_j[r];
  fe known = 0, a = 0;
  for (int z = 0; z < frv_support_n[g]; ++z) {
    int p = frv_support_idx[g][z];
    fe c = frv_support_val[g][z];
    int t = pp_shift(p, ri, rj);
    if (pp_w[t] < n)
      known ^= mul(c, pp_S[t]);
    else {
      a ^= mul(c, pp_lam[t]);
      known ^= mul(c, pp_K[t]);
    }
  }
  *an = a;
  return known;
}
static int poly_nonzero_sparse(int g, int k) {
  fe acc = 0;
  for (int z = 0; z < frv_support_n[g]; ++z) {
    int p = frv_support_idx[g][z];
    acc ^=
        mul(frv_support_val[g][z], mono_eval(PXv[k], PYv[k], pp_i[p], pp_j[p]));
  }
  return acc != 0;
}

static int vote_position(int n, Result *R, int verbose) {
  fe pred[MAXGEN], dv[MAXGEN];
  int validg[MAXGEN], dc[MAXGEN], nd = 0, ncand = 0;
  if (use_poly)
    poly_prepare(n);
  for (int g = 0; g < ngen; ++g) {
    int m = n - gens[g].b;
    validg[g] = (m >= 0 && in_H[m]);
    if (validg[g]) {
      fe an, known = (use_poly && frv_support_ready)
                         ? poly_ext_known_sparse(g, m, n, &an)
                         : ext_known_b(gens[g].f, gens[g].b, m, n, &an);
      if (!an) {
        fprintf(stderr, "vote: zero leading coefficient (n=%d)\n", n);
        return -1;
      }
      pred[g] = mul(known, inv(an));
    }
  }
  if (frv_vote_conv_ready && frv_support_ready) {
    /* Count all (i,j) pairs with i+j=n by locator ancestor using exact
     * integer convolutions.  Group equal predicted field values afterwards. */
    for (int g = 0; g < ngen; ++g) {
      uint32_t cnt = frv_vote_counts[(size_t)g * frv_vote_count_len +
                                     (size_t)(n - M_SYND - 1)];
      if (!validg[g] || !cnt)
        continue;
      ncand += (int)cnt;
      int s = 0;
      while (s < nd && dv[s] != pred[g])
        ++s;
      if (s == nd) {
        dv[nd] = pred[g];
        dc[nd] = 0;
        ++nd;
      }
      dc[s] += (int)cnt;
    }
  } else {
    /* Exact reference path, retained as a fallback. */
    for (int i = 0; i <= n; ++i) {
      int j = n - i, anc = -1, s;
      if (!in_H[i] || is_piv[i] || !in_H[j] || colused[j])
        continue;
      anc = frv_support_ready ? frv_ancestor[i] : -1;
      if (!frv_support_ready) {
        for (int g = 0; g < ngen; ++g)
          if (gens[g].b <= i && in_H[i - gens[g].b]) {
            anc = g;
            break;
          }
      }
      if (anc < 0 || !validg[anc])
        continue;
      ++ncand;
      for (s = 0; s < nd; ++s)
        if (dv[s] == pred[anc])
          break;
      if (s == nd) {
        dv[nd] = pred[anc];
        dc[nd] = 0;
        ++nd;
      }
      ++dc[s];
    }
  }
  if (nd == 0) {
    fprintf(stderr, "vote: no candidates at n=%d\n", n);
    return -1;
  }
  int best = 0, second = 0;
  for (int s = 1; s < nd; ++s)
    if (dc[s] > dc[best])
      best = s;
  for (int s = 0; s < nd; ++s)
    if (s != best && dc[s] > second)
      second = dc[s];
  xhat[n] = dv[best];
  if (use_poly)
    poly_finalize(n, xhat[n]);
  int margin = dc[best] - second;
  if (frv_vote_conv_ready && frv_support_ready && margin <= 0) {
    /* The convolution gives exact vote counts, but not the order in which
     * equal-count predictions first appeared.  Re-run only tied steps with
     * the original i-ascending scan to preserve the previous tie-break. */
    nd = ncand = 0;
    for (int i = 0; i <= n; ++i) {
      int j = n - i, anc = frv_ancestor[i], ss;
      if (anc < 0 || !in_H[j] || colused[j] || !validg[anc])
        continue;
      ++ncand;
      for (ss = 0; ss < nd; ++ss)
        if (dv[ss] == pred[anc])
          break;
      if (ss == nd) {
        dv[nd] = pred[anc];
        dc[nd] = 0;
        ++nd;
      }
      ++dc[ss];
    }
    best = 0;
    second = 0;
    for (int ss = 1; ss < nd; ++ss)
      if (dc[ss] > dc[best])
        best = ss;
    for (int ss = 0; ss < nd; ++ss)
      if (ss != best && dc[ss] > second)
        second = dc[ss];
    margin = dc[best] - second;
  }
  ++R->vote_steps;
  if (margin < R->vote_min_margin)
    R->vote_min_margin = margin;
  if (margin <= 0)
    ++R->vote_tie;
  int right = true_syndrome(n) == xhat[n];
  if (!right)
    ++R->vote_wrong;
  if (verbose)
    printf("#   vote n=%d: %d candidates, %d predictions, winner %d votes, "
           "runner-up %d%s\n",
           n, ncand, nd, dc[best], second,
           right ? "  (= true syndrome)" : "  (WRONG)");
  return 0;
}

static void bms_cleanup(void);

/* ================= 純粋なFeng-Rao多数決復号 ==========================
 * 誤り位置消去イデアルは、Sakata/BMSの逐次更新ではなく、
 * Feng-Raoシンドローム行列の核（feng_rao_reduce）から構成する。
 * その最小半群生成元を使ってFeng-Rao多数決を行う。
 * Newton多項式シンドローム表現を必須とし、従来VOTEの補間双対基底は使わない。
 */
static int decode_fr_vote_pass(Result *R, int verbose, int quiet,
                               int row_limit) {
  double t0 = now();
  int nrows, nl, nsel = 0, nfound = 0, overflow = 0;
  int sel[MAXGEN];

  memset(R, 0, sizeof *R);
  R->vote_min_margin = 1 << 30;
  if (!use_poly) {
    fprintf(stderr, "FR-VOTE requires the Newton monomial syndrome path; "
                    "enable it with the default Newton setup (do not use "
                    "-nopoly/-nonewton).\n");
    R->failed = 1;
    return -1;
  }

  memset(xhat, 0, sizeof xhat);
  memset(is_piv, 0, sizeof is_piv);
  memset(colused, 0, sizeof colused);
  bms_cleanup();

  nrows = feng_rao_reduce(row_limit);
  if (nrows < 0) {
    fprintf(stderr, "FR-VOTE: Feng-Rao reduction failed.\n");
    R->failed = 1;
    return -1;
  }
  nl = count_basis_upto(L_DEG);

  /* Feng-Rao簡約のピボット階段を半群次数の集合へ移す。 */
  for (int i = 0; i < nrows; ++i) {
    if (piv_of_row[i] >= 0) {
      int b = basis[i].rho;
      int c = basis[piv_of_row[i]].rho;
      if (b < SGMAX) {
        is_piv[b] = 1;
        R->rowset[b] = 1;
      }
      if (c < SGMAX) {
        colused[c] = 1;
        R->colset[c] = 1;
      }
      ++R->T;
    }
  }
  if (R->T == 0) {
    R->resolved = 1;
    R->t_main = now() - t0;
    return 0;
  }

  /* 非ピボット行は誤り位置消去関数になる。
   * Weierstrass半群で割り切れない最小生成元だけを残す。 */
  for (int b = 0; b <= L_DEG; ++b) {
    int i, minimal = 1;
    if (!in_H[b] || is_piv[b])
      continue;
    i = idx_of_ord[b];
    if (i < 0 || i >= nl || i >= nrows || piv_of_row[i] >= 0)
      continue;
    for (int q = 0; q < ngenH; ++q) {
      int a = b - genH[q];
      if (a >= 0 && in_H[a] && !is_piv[a]) {
        minimal = 0;
        break;
      }
    }
    if (!minimal)
      continue;
    if (ngen >= MAXGEN) {
      fprintf(stderr, "FR-VOTE: too many locator generators.\n");
      R->failed = 1;
      bms_cleanup();
      return -1;
    }
    gens[ngen].b = b;
    gens[ngen].d = 0;
    gens[ngen].f = vnew();
    for (int k = 0; k <= i; ++k) {
      fe c = coefm[i][k];
      if (!c)
        continue;
      int z = pp_at(bas_ei[k], bas_ej[k]);
      if (z < 0 || z >= pp_upto[b]) {
        fprintf(stderr,
                "FR-VOTE: locator coefficient is outside Newton table "
                "(b=%d, basis=%d, index=%d).\n",
                b, k, z);
        R->failed = 1;
        bms_cleanup();
        return -1;
      }
      gens[ngen].f[z] ^= c;
    }
    ++ngen;
  }

  if (ngen == 0) {
    fprintf(
        stderr,
        "FR-VOTE: no locator generators were obtained from the FR kernel.\n");
    R->failed = 1;
    return -1;
  }
  if (frv_build_sparse() != 0) {
    fprintf(stderr, "FR-VOTE: cannot allocate sparse locator supports.\n");
    R->failed = 1;
    bms_cleanup();
    return -1;
  }

  /* 未知のg層のシンドロームをFeng-Rao多数決で順に復元する。 */
  for (int n = M_SYND + 1; n <= M_EXT; ++n) {
    if (vote_position(n, R, verbose) != 0) {
      R->failed = 1;
      bms_cleanup();
      return -1;
    }
  }
  R->t_main = now() - t0;

  /* Feng-Rao限界に従って必要な消去関数を選び、その共通零点を誤り位置とする。 */
  {
    int bmax = M_EXT - 2 * GENUS + 1 - R->T;
    for (int g = 0; g < ngen; ++g)
      if (gens[g].b <= bmax)
        sel[nsel++] = g;
  }
  if (!quiet) {
    printf("# FR-VOTE: %d FR-kernel locator generators, %d used (rho <= %d)\n",
           ngen, nsel, M_EXT - 2 * GENUS + 1 - R->T);
    printf("# FR-VOTE: staircase of %d elements; %d voted syndromes, "
           "smallest margin %d, %d wrong, %d ties\n",
           R->T, R->vote_steps, R->vote_min_margin, R->vote_wrong, R->vote_tie);
  }

  t0 = now();
  if (nsel == 0) {
    fprintf(stderr,
            "FR-VOTE: no locator generators satisfy the Feng-Rao bound.\n");
    R->failed = 1;
    bms_cleanup();
    return -1;
  }
  for (int k = 0; k < npts; ++k) {
    int vanish = 1;
    for (int q = 0; q < nsel && vanish; ++q)
      if (frv_support_ready ? poly_nonzero_sparse(sel[q], k)
                            : poly_nonzero(gens[sel[q]].f, gens[sel[q]].b, k))
        vanish = 0;
    if (vanish) {
      if (nfound < MAXT)
        R->pos[nfound++] = k;
      else
        overflow = 1;
    }
  }
  R->nfound = nfound;
  R->overflow = overflow;
  R->t_extract = now() - t0;
  if (!quiet) {
    puts("# FR-VOTE locator zero counts on the evaluation-point set:");
    for (int q = 0; q < nsel; ++q) {
      int zeros = 0, true_hits = 0;
      for (int p = 0; p < npts; ++p) {
        int iszero = frv_support_ready
                         ? !poly_nonzero_sparse(sel[q], p)
                         : !poly_nonzero(gens[sel[q]].f, gens[sel[q]].b, p);
        if (iszero) {
          ++zeros;
          for (int e = 0; e < NUM_ERRORS; ++e)
            if (true_pos[e] == p) {
              ++true_hits;
              break;
            }
        }
      }
      printf(
          "#   G_%d (rho=%d): zeros=%d, injected-error points among them=%d\n",
          q, gens[sel[q]].b, zeros, true_hits);
    }
    int common_true = 0;
    for (int e = 0; e < NUM_ERRORS; ++e)
      for (int k = 0; k < R->nfound; ++k)
        if (true_pos[e] == R->pos[k]) {
          ++common_true;
          break;
        }
    printf("# FR-VOTE common zero set: %d points; injected-error points in "
           "it=%d/%d\n",
           R->nfound, common_true, NUM_ERRORS);
  }
  bms_cleanup();
  resolve_and_store(R);
  return 0;
}

static void bms_cleanup(void) {
  for (int g = 0; g < MAXGEN; ++g) {
    free(frv_support_idx[g]);
    free(frv_support_val[g]);
    frv_support_idx[g] = NULL;
    frv_support_val[g] = NULL;
    frv_support_n[g] = 0;
  }
  frv_support_ready = 0;
  free(frv_vote_counts);
  frv_vote_counts = NULL;
  frv_vote_count_len = 0;
  frv_vote_conv_ready = 0;
  for (int g = 0; g < ngen; ++g)
    free(gens[g].f);
  ngen = 0;
  for (int m = 0; m < SGMAX; ++m) {
    free(pivvec[m]);
    pivvec[m] = NULL;
  }
}

/* Fast FR-VOTE first reduces only the rows that can produce locator
 * coefficient vectors.  This is a candidate solution, not an unchecked
 * truncation: the existing error-value solver verifies the recovered values
 * against the complete available syndrome sequence.  If it cannot certify
 * the result, retry with the original full row range. */
static int frv_known_syndromes_nonzero(void) {
  for (int p = 0; p < NPT; ++p)
    if (pp_w[p] <= M_SYND && pp_S[p])
      return 1;
  return 0;
}
static int decode_fr_vote(Result *R, int verbose, int quiet) {
  /* Small prefixes are tried first.  A result is accepted only when the
   * existing full-syndrome error-value check certifies it; otherwise enlarge
   * the prefix and ultimately fall back to the original full reduction. */
  int limits[] = {4096, 8192, 16384, count_basis_upto(L_DEG), 0};
  const int nlimits = (int)(sizeof(limits) / sizeof(limits[0]));
  int rc = -1, previous_limit = -1;
  for (int pass = 0; pass < nlimits; ++pass) {
    int lim = limits[pass];
    if (lim > 0 && lim <= previous_limit)
      continue;
    previous_limit = lim;
    rc = decode_fr_vote_pass(R, verbose, quiet, lim);
    if (rc == 0 && R->resolved && (R->T != 0 || !frv_known_syndromes_nonzero()))
      return 0;
    if (pass + 1 < nlimits) {
      fprintf(stderr,
              "# FR-VOTE: prefix %d was not certified; expanding reduction.\n",
              lim ? lim : count_basis_upto(ROW_MAX));
      bms_cleanup();
      poly_syndromes(); /* discard tentative unknown syndromes */
    }
  }
  if (rc == 0 && R->T == 0 && frv_known_syndromes_nonzero()) {
    R->resolved = 0;
    R->failed = 1;
    return -1;
  }
  return rc;
}

static int bms_run(Result *R, int nmax, int verbose) {
  int T_ = 0;
  memset(colused, 0, sizeof colused);
  memset(is_piv, 0, sizeof is_piv);
  memset(pivvec, 0, sizeof pivvec);
  ngen = 1;
  gens[0].b = 0;
  gens[0].d = 0;
  gens[0].f = vnew();
  if (use_poly)
    gens[0].f[pp_at(0, 0)] = 1;
  else
    for (int k = 0; k < npts; ++k)
      gens[0].f[k] = 1;

  for (int n = 0; n <= nmax; ++n) {
    int anyfail = 0;
    if (n > M_SYND && vote_position(n, R, verbose) != 0)
      return -1;
    for (int g = 0; g < ngen; ++g) {
      int m = n - gens[g].b;
      gens[g].d = (m >= 0 && in_H[m]) ? xdisc_b(gens[g].f, gens[g].b, m, n) : 0;
      if (gens[g].d)
        anyfail = 1;
    }
    if (!anyfail)
      continue;
    for (int a = 0; a <= n; ++a) {
      int m = n - a, src = -1;
      if (!in_H[a] || is_piv[a] || !in_H[m] || colused[m])
        continue;
      for (int g = 0; g < ngen; ++g)
        if (gens[g].d && gens[g].b <= a && in_H[a - gens[g].b]) {
          src = g;
          break;
        }
      if (src < 0)
        continue;
      fe *f = vnew();
      vmul_b(f, a - gens[src].b, gens[src].f, gens[src].b);
      fe d = xdisc_b(f, a, m, n);
      if (!d) {
        fprintf(stderr, "BMS: internal error (zero pivot at n=%d)\n", n);
        free(f);
        return -1;
      }
      pivvec[m] = f;
      if (use_poly)
        pivb[m] = a;
      pivd[m] = d;
      colused[m] = 1;
      is_piv[a] = 1;
      R->rowset[a] = 1;
      R->colset[m] = 1;
      ++T_;
    }
    Gen ng[MAXGEN];
    uint8_t reused[MAXGEN];
    int nn = 0;
    memset(reused, 0, sizeof reused);
    for (int b = 0; b <= nmax; ++b) {
      int minimal = 1, og = -1;
      if (!in_H[b] || is_piv[b])
        continue;
      for (int s = 0; s < ngenH; ++s)
        if (b - genH[s] >= 0 && in_H[b - genH[s]] && !is_piv[b - genH[s]]) {
          minimal = 0;
          break;
        }
      if (!minimal)
        continue;
      if (nn >= MAXGEN) {
        fprintf(stderr, "BMS: too many generators\n");
        return -1;
      }
      for (int g = 0; g < ngen; ++g)
        if (gens[g].b == b)
          og = g;
      if (og >= 0) {
        ng[nn] = gens[og];
        reused[og] = 1;
      } else {
        int src = -1, m = n - b;
        for (int g = 0; g < ngen; ++g)
          if (gens[g].b <= b && in_H[b - gens[g].b]) {
            src = g;
            break;
          }
        if (src < 0) {
          fprintf(stderr, "BMS: no source for generator %d\n", b);
          return -1;
        }
        ng[nn].b = b;
        ng[nn].f = vnew();
        vmul_b(ng[nn].f, b - gens[src].b, gens[src].f, gens[src].b);
        ng[nn].d = (m >= 0 && in_H[m]) ? xdisc_b(ng[nn].f, b, m, n) : 0;
      }
      {
        int m = n - ng[nn].b;
        if (m >= 0 && in_H[m] && ng[nn].d) {
          if (!colused[m] || !pivvec[m]) {
            fprintf(stderr, "BMS: inconsistent state (n=%d, b=%d)\n", n, b);
            return -1;
          }
          vaxpy_b(ng[nn].f, mul(ng[nn].d, inv(pivd[m])), pivvec[m],
                  use_poly ? pivb[m] : 0);
          ng[nn].d = 0;
        }
      }
      ++nn;
    }
    for (int g = 0; g < ngen; ++g)
      if (!reused[g])
        free(gens[g].f);
    for (int g = 0; g < nn; ++g)
      gens[g] = ng[g];
    ngen = nn;
  }
  R->T = T_;
  return 0;
}

static int decode_bms(Result *R, int use_vote, int verbose, int quiet) {
  const char *name = use_vote ? "VOTE" : "BMS";
  const int Mtop = use_vote ? M_EXT : M_SYND;
  double t0 = now();
  memset(R, 0, sizeof *R);
  R->vote_min_margin = 1 << 30;
  memset(xhat, 0, sizeof xhat);
  if (bms_run(R, Mtop, verbose) != 0) {
    R->failed = 1;
    bms_cleanup();
    return -1;
  }
  R->t_main = now() - t0;
  int bmax = use_vote ? Mtop - 2 * GENUS + 1 - R->T : L_DEG;
  int nsel = 0, sel[MAXGEN];
  for (int g = 0; g < ngen; ++g)
    if (gens[g].b <= bmax)
      sel[nsel++] = g;
  if (!quiet) {
    printf("# %-4s: %d generators of E, %d with rho <= %d\n", name, ngen, nsel,
           bmax);
    printf("# %-4s: staircase of %d elements -> estimated t = %d\n", name, R->T,
           R->T);
    if (use_vote)
      printf("# %-4s: %d voted syndromes, smallest margin %d, %d wrong (vs. "
             "simulated truth), %d ties\n",
             name, R->vote_steps, R->vote_min_margin, R->vote_wrong,
             R->vote_tie);
  }
  if (R->T == 0) {
    bms_cleanup();
    return 0;
  }
  t0 = now();
  int nfound = 0, overflow = 0;
  for (int k = 0; k < npts; ++k) {
    int vanish = nsel > 0;
    for (int s = 0; s < nsel && vanish; ++s)
      if (use_poly ? poly_nonzero(gens[sel[s]].f, gens[sel[s]].b, k)
                   : gens[sel[s]].f[k])
        vanish = 0;
    if (vanish) {
      if (nfound < MAXT)
        R->pos[nfound++] = k;
      else
        overflow = 1;
    }
  }
  R->nfound = nfound;
  R->overflow = overflow;
  R->t_extract = now() - t0;
  bms_cleanup();
  resolve_and_store(R);
  return 0;
}

/* ---- error simulation / verification ---- */
static fe true_val[MAXT];

/* Print the injected values next to the values recovered by one decoder. */
static void show_error_values(const char *name, const Result *R) {
  printf("# %s error details (%d injected, %d recovered)\n", name, NUM_ERRORS,
         R->nfound);
  puts("# pos       (X,Y)          error     before    after     corrected  "
       "match");
  for (int e = 0; e < NUM_ERRORS; ++e) {
    fe got = 0;
    int found = 0;
    int p = true_pos[e];
    for (int k = 0; k < R->nfound; ++k) {
      if (R->pos[k] == p) {
        got = R->val[k];
        found = 1;
        break;
      }
    }
    /* This simulator transmits the all-zero codeword.  Thus before=0,
       after=the injected symbol, and corrected=after XOR decoded error. */
    fe before = 0;
    fe after = recv ? recv[p] : true_val[e];
    fe corrected = found ? (fe)(after ^ got) : after;
    printf("# %-9d (%03X,%03X)   0x%03X    0x%03X    0x%03X    0x%03X    %s\n",
           p, (unsigned)PXv[p], (unsigned)PYv[p], (unsigned)true_val[e],
           (unsigned)before, (unsigned)after, (unsigned)corrected,
           found && got == true_val[e] && corrected == before
               ? "OK"
               : (found ? "MISMATCH" : "MISSING"));
  }
  puts("# Position condition: (X,Y) is a common zero of the displayed "
       "L_i(X,Y).");
}

static int verify(const char *name, const Result *R, int quiet) {
  int same = 0;
  if (R->failed) {
    if (!quiet)
      printf("# %-5s FAILED\n", name);
    return 0;
  }
  for (int m = 0; m < R->nfound; ++m)
    for (int e = 0; e < NUM_ERRORS; ++e)
      if (true_pos[e] == R->pos[m] && true_val[e] == R->val[m])
        ++same;
  int ok = !R->overflow && R->resolved && R->nfound == NUM_ERRORS &&
           same == NUM_ERRORS;
  if (!quiet)
    printf("# %-5s %d of %d errors recovered with the correct value, %d "
           "reported -> %s\n",
           name, same, NUM_ERRORS, R->nfound, ok ? "OK" : "MISMATCH");
  return ok;
}
static double total_time(const Result *R) {
  return R->t_main + R->t_extract + R->t_resolve;
}

static int fixed_error_value = 0; /* 0=random nonzero; otherwise force every
                                     injected error to this field element */
static uint64_t rng_s = 88172645463325252ULL;
static uint64_t rnd(void) {
  rng_s ^= rng_s << 13;
  rng_s ^= rng_s >> 7;
  rng_s ^= rng_s << 17;
  return rng_s;
}

/* draw NUM_ERRORS random errors; fills r (length npts), true_pos/true_val */
static void make_errors(fe *r) {
  memset(r, 0, sizeof(fe) * npts);
  for (int e = 0; e < NUM_ERRORS; ++e) {
    int p, dup;
    do {
      p = (int)(rnd() % (uint64_t)npts);
      dup = 0;
      for (int k = 0; k < e; ++k)
        if (true_pos[k] == p)
          dup = 1;
    } while (dup);
    true_pos[e] = p;
    if (fixed_error_value > 0) {
      r[p] = (fe)fixed_error_value;
    } else {
      do {
        r[p] = (fe)(rnd() % (uint64_t)FS);
      } while (!r[p]);
    }
    true_val[e] = r[p];
  }
}
static void set_received(fe *r) {
  recv = r;
  nz_cnt = 0;
  for (int k = 0; k < npts; ++k)
    if (r[k])
      nz_idx[nz_cnt++] = k;
  if (use_poly) {
    static int checked = 0;
    poly_syndromes();
    if (!checked) {
      checked = 1;
      int bad = poly_selfcheck();
      if (bad)
        fprintf(stderr,
                "# WARNING: monomial relation self-check failed at %d "
                "points\n",
                bad);
    }
  }
}

/* ================================================================== */
/*  multi-point codes: basic algorithm  S(h_a f_b) = sum_k e_k h_a f_b */
/* ================================================================== */
static int solve_values(const fe *VG, int rG, const fe *sigma, const int *Z,
                        int nz, fe *val) {
  int W = nz + 1;
  fe *A = malloc(sizeof(fe) * (size_t)rG * W);
  for (int s = 0; s < rG; ++s) {
    for (int c = 0; c < nz; ++c)
      A[(size_t)s * W + c] = VG[(size_t)s * npts + Z[c]];
    A[(size_t)s * W + nz] = sigma[s];
  }
  for (int c = 0; c < nz; ++c) {
    int p = -1;
    for (int a = c; a < rG; ++a)
      if (A[(size_t)a * W + c]) {
        p = a;
        break;
      }
    if (p < 0) {
      free(A);
      return 0;
    }
    if (p != c)
      for (int k = 0; k < W; ++k) {
        fe t = A[(size_t)p * W + k];
        A[(size_t)p * W + k] = A[(size_t)c * W + k];
        A[(size_t)c * W + k] = t;
      }
    fe iv = inv(A[(size_t)c * W + c]);
    for (int k = c; k < W; ++k)
      A[(size_t)c * W + k] = mul(A[(size_t)c * W + k], iv);
#pragma omp parallel for schedule(static)
    for (int a = 0; a < rG; ++a) {
      fe f = A[(size_t)a * W + c];
      if (a == c || !f)
        continue;
      for (int k = c; k < W; ++k)
        A[(size_t)a * W + k] ^= mul(f, A[(size_t)c * W + k]);
    }
  }
  for (int a = nz; a < rG; ++a)
    if (A[(size_t)a * W + nz]) {
      free(A);
      return 0;
    }
  for (int c = 0; c < nz; ++c)
    val[c] = A[(size_t)c * W + nz];
  free(A);
  return 1;
}

/* ---- what parameters does this build / this curve handle? ---- */
static int limit_Mmax(int g, int n) {
  int mm = n - g -
           2; /* M+g <= n-2: injective evaluation, one point may be dropped */
  if (mm > MAXB - 2)
    mm = MAXB - 2; /* basis array holds M+1 functions */
  return mm;
}
static int suggest_param_mode = 0;

/* Print a compact parameter frontier for the current one-point AG curve.
 * The bounds below are the same designed bounds used by print_limits(). */
static void print_suggest_params(int g, int n, int requested_t) {
  int Mmin = 2 * g;
  int Mmax = limit_Mmax(g, n);
  int cap = MAXT - 2;
  puts("# ---- suggested one-point AG parameters ----");
  printf("# exact curve data: genus g=%d, affine positions n=%d\n", g, n);
  printf("# build limits: M<=%d, errors<=%d\n", MAXB - 2, cap);
  if (Mmax < Mmin) {
    printf("# no feasible M: need M >= %d but injectivity/build limit gives M <= %d\n",
           Mmin, Mmax);
    puts("# condition: n > 3g+1 and M+g < n");
    return;
  }
  puts("# table is sorted by descending M (higher rate / larger syndrome bound first)");
  puts("#       M       k_RR      rate       t_FR/BMS    t_vote");
  int span = Mmax - Mmin;
  int step = (span + 98) / 99;
  if (step < 1) step = 1;
  int last = -1;
  for (int M = Mmax; M >= Mmin; --M) {
    if (M != Mmax && M != Mmin && ((Mmax - M) % step) != 0)
      continue;
    int k = M - g + 1;
    int tb = (M - 3 * g + 1) / 2;
    int tv = (M - 2 * g + 1) / 2;
    if (tb < 0) tb = 0;
    if (tv < 0) tv = 0;
    if (tb > cap) tb = cap;
    if (tv > cap) tv = cap;
    printf("# %7d %9d %10.6f %11d %10d\n", M, k,
           n > 0 ? (double)k / (double)n : 0.0, tb, tv);
    last = M;
  }
  if (last != Mmin) {
    int k = Mmin - g + 1;
    int tb = (Mmin - 3 * g + 1) / 2;
    int tv = (Mmin - 2 * g + 1) / 2;
    if (tb < 0) tb = 0;
    if (tv < 0) tv = 0;
    if (tb > cap) tb = cap;
    if (tv > cap) tv = cap;
    printf("# %7d %9d %10.6f %11d %10d\n", Mmin, k,
           n > 0 ? (double)k / (double)n : 0.0, tb, tv);
  }
  if (requested_t >= 0) {
    long long Mfr = 2LL * requested_t + 3LL * g - 1;
    long long Mvote = 2LL * requested_t + 2LL * g;
    printf("# for requested errors t=%d:\n", requested_t);
    printf("#   FR/BMS needs M >= %lld; ", Mfr);
    if (Mfr <= Mmax && Mfr <= MAXB - 2)
      printf("feasible, k=%lld, rate=%.6f\n", Mfr - g + 1,
             (double)(Mfr - g + 1) / (double)n);
    else
      puts("not feasible under current n/build limits");
    printf("#   voting needs M >= %lld; ", Mvote);
    if (Mvote <= Mmax && Mvote <= MAXB - 2)
      printf("feasible, k=%lld, rate=%.6f\n", Mvote - g + 1,
             (double)(Mvote - g + 1) / (double)n);
    else
      puts("not feasible under current n/build limits");
  }
  puts("# k_RR=M-g+1 is the Riemann-Roch dimension; actual decoding still depends on a certified basis.");
}

static void print_limits(int g, int n, int t) {
  int Mmax = limit_Mmax(g, n), cap = MAXT - 2;
  puts("# ---- parameter limits (this build, this curve) ----");
  printf("#   field GF(2^m): m <= 12;  here g = %d, n = %d positions, build "
         "caps: MAXB = %d, MAXT = %d\n",
         g, n, MAXB, MAXT);
  if (Mmax < 2 * g) {
    printf("#   NOT decodable here: needs M >= 2g = %d but M <= %d  (n must be "
           "> 3g+1 = %d and MAXB >= 3g+3)\n",
           2 * g, Mmax, 3 * g + 1);
    return;
  }
  int tb = (Mmax - 3 * g + 1) / 2, tv = (Mmax - 2 * g + 1) / 2;
  if (tb > cap)
    tb = cap;
  if (tv > cap)
    tv = cap;
  printf("#   syndrome bound M: %d .. %d\n", 2 * g, Mmax);
  printf("#   max errors at M = %d:  basic FR/BMS %d,  majority voting %d  "
         "(designed bounds)\n",
         Mmax, tb < 0 ? 0 : tb, tv);
  int Mb = 2 * t + 3 * g - 1, Mv = 2 * t + 2 * g;
  printf("#   for t = %d errors: needs M >= %d (FR/BMS) %s, M >= %d (voting) "
         "%s%s\n",
         t, Mb, Mb <= Mmax ? "[OK]" : "[too large]", Mv,
         Mv <= Mmax ? "[OK]" : "[too large]",
         t > cap ? "  [t exceeds MAXT-2]" : "");
}

static void usage(void) {
  puts(
      "noize [options] [seed]\n"
      "  --m m            base field GF(2^m), m <= 12             (default 6)\n"
      "  --poly 0x43      primitive polynomial of GF(2^m)       (default: "
      "standard)\n"
      "  --curve bbgs|hermitian|suzuki|klein  select a built-in curve family "
      "(default "
      "bbgs)\n"
      "  --q auto|power-of-2   BBGS parameter q (power of 2); Suzuki: "
      "q=2^(2n+1), default 8 (own engine; uses --M --errors -n)\n"
      "  --a A --b B      BBGS coprime positive parameters (default 1,2); "
      "c=a+b\n"
      "  BBGS field is GF(q^c), c=a+b; curve has c+1 terms\n"
      "  --terms \"16,0;0,1;15,4;12,16\"   custom plane curve instead (X^i "
      "Y^j, optionally \"i,j,coef\"); non-monic / unsupported models are "
      "normalised automatically (swap, Y->Y+cX^k, ...)\n"
      "  --errors t       number of injected errors              (default "
      "256)\n"
      "  --fixed-error-value v  use the same nonzero field value at every "
      "error\n"
      "  -common-binary    experimental GF(2) position-only decoder (small "
      "cases)\n"
      "  -bms2-check       build S(f_i f_j), compute its GF(2^m) nullspace,\n"
      "                    scan common zeros and solve error values (small "
      "cases)\n"
      "  --M M            syndrome bound / pole order (default 2t+2g)\n"
      "  --len L | -half  puncture: use only the first L points (or half of "
      "them) as positions\n"
      "  --place k        pole place (index from the list); default: automatic "
      "choice\n"
      "  --basis legacy|newton   multi-point elimination column order (default legacy)\n"
      "  --pole \"E2:140,E1:85\" or \"#3:100\"   multi-point divisor (E<k>: "
      "edge, #k: orbit)\n"
      "  -fr -bms -vote   one-point decoders to run (default: all)\n"
      "  -ecp             independent Error-Correcting-Pair decoder (one-point AG)\n"
      "  -frvote          pure Feng-Rao majority voting (Newton path; no BMS, "
      "no dual basis)\n"
      "  -forney          use multivariate AG Forney-function evaluation for "
      "FR\n"
      "  -kotter          use Kötter one-point update + generalized Forney "
      "formula\n"
      "  -show-values     print injected and decoded error values\n"
      "  -strict          skip FR/BMS when t exceeds the basic bound tau "
      "(default: run them anyway)\n"
      "  -limits          print the parameter limits of this build/curve and "
      "exit\n"
      "  --suggest-param  print a descending-M parameter/rate frontier and exit\n"
      "  -basic           (multi-point) Niederreiter syndrome decoding test\n"
      "                    (default: error vector -> syndrome -> decode)\n"
      "  -encode          additionally generate a random AG codeword, add "
      "errors,\n"
      "                    decode, and verify exact recovery (with -basic)\n"
      "  --concat-inner goppa:m:t[:n]  one-point Newton AG + binary Goppa concatenation (default);\n"
      "                                add --pole to opt in to multi-point AG concatenation\n"
      "  --concat-burst  inject a contiguous burst in the transmitted bit stream\n"
      "  --interleave D  interleave D adjacent inner blocks (default 1 = off)\n"
      "  -nonewton        disabled: this build requires Newton-polygon bases\n"
      "  --place k        also avoids probing every candidate place (saves "
      "setup time)\n"
      "  -n N             N random trials (summary only)       -v  votes     "
      "-list  only list places");
}

/* ================= Suzuki curve: dedicated one-point engine =================
 *   y^q + y = x^q0 (x^q + x),   q = 2^(2n+1),  q0 = 2^n,   genus g = q0 (q-1)
 * The Weierstrass semigroup at the single place at infinity is
 * <q, q+q0, q+2q0, q+2q0+1>; the plane monomials x^i y^j do not realise it, so
 * L(M P_inf) is spanned by the monomials x^a y^b z^c w^d with
 *   z = x^(2q0+1) + y^(2q0),   w = x y^(2q0) + z^(2q0)
 * of pole orders q, q+q0, q+2q0, q+2q0+1 (one monomial per semigroup element).
 * Decoding is the Feng-Rao basic algorithm on the semigroup syndrome matrix
 * S[i][j] = s(phi_i phi_j) (only function values at the points are needed),
 * followed by the univariate polynomial g(theta) = det(theta I + S0^-1 SX)
 * whose roots are the x-coordinates of the errors, and a linear solve for the
 * values. Guaranteed: t <= (M - 3g + 1) / 2.  Independent of the plane-curve
 * machinery of the rest of this file (own working arrays; uses only the field
 * above). */
static int sz_q, sz_q0, sz_g, sz_M, sz_tau, sz_ldeg, sz_rowmax;
static int sz_np;
static fe *sz_px, *sz_py;
static int *sz_first, *sz_count;
static int sz_nb;
static int *sz_rho, *sz_cnt;
static fe **sz_val;

static void sz_free_all(void) {
  free(sz_px);
  free(sz_py);
  free(sz_first);
  free(sz_count);
  free(sz_rho);
  free(sz_cnt);
  if (sz_val)
    for (int i = 0; i < sz_nb; ++i)
      free(sz_val[i]);
  free(sz_val);
  sz_px = sz_py = NULL;
  sz_first = sz_count = sz_rho = sz_cnt = NULL;
  sz_val = NULL;
  sz_nb = 0;
}
static int sz_count_upto(int rho) {
  if (rho < 0)
    return 0;
  if (rho > sz_M)
    rho = sz_M;
  return sz_cnt[rho];
}
static int sz_enumerate(void) {
  const int q = sz_q, q0 = sz_q0;
  fe *fib = malloc((size_t)FS * (size_t)q * sizeof(fe));
  int *fc = calloc((size_t)FS, sizeof(int));
  sz_first = malloc((size_t)FS * sizeof(int));
  sz_count = malloc((size_t)FS * sizeof(int));
  if (!fib || !fc || !sz_first || !sz_count) {
    free(fib);
    free(fc);
    return -1;
  }
  for (int y = 0; y < FS; ++y) {
    int v = gpow((fe)y, q) ^ y;
    if (fc[v] >= q) {
      free(fib);
      free(fc);
      return -1;
    }
    fib[(size_t)v * q + fc[v]++] = (fe)y;
  }
  size_t total = 0;
  for (int x = 0; x < FS; ++x) {
    int rhs = mul(gpow((fe)x, q0), (fe)(gpow((fe)x, q) ^ x));
    total += (size_t)fc[rhs];
  }
  sz_px = malloc((total + 1) * sizeof(fe));
  sz_py = malloc((total + 1) * sizeof(fe));
  if (!sz_px || !sz_py) {
    free(fib);
    free(fc);
    return -1;
  }
  size_t np = 0;
  for (int x = 0; x < FS; ++x) {
    int rhs = mul(gpow((fe)x, q0), (fe)(gpow((fe)x, q) ^ x));
    sz_first[x] = (int)np;
    sz_count[x] = fc[rhs];
    for (int j = 0; j < fc[rhs]; ++j) {
      sz_px[np] = (fe)x;
      sz_py[np] = fib[(size_t)rhs * q + j];
      ++np;
    }
  }
  sz_np = (int)np;
  free(fib);
  free(fc);
  return 0;
}
static int sz_build_basis(void) {
  const int q = sz_q, q0 = sz_q0, M = sz_M;
  const int wy = q + q0, wz = q + 2 * q0, ww = q + 2 * q0 + 1;
  sz_rho = malloc((size_t)(M + 1) * sizeof(int));
  sz_cnt = malloc((size_t)(M + 1) * sizeof(int));
  int *ea = malloc((size_t)(M + 1) * sizeof(int)),
      *eb = malloc((size_t)(M + 1) * sizeof(int)),
      *ec = malloc((size_t)(M + 1) * sizeof(int)),
      *ed = malloc((size_t)(M + 1) * sizeof(int));
  if (!sz_rho || !sz_cnt || !ea || !eb || !ec || !ed) {
    free(ea);
    free(eb);
    free(ec);
    free(ed);
    return -1;
  }
  int nb = 0;
  for (int s = 0; s <= M; ++s) {
    int found = 0, fa = 0, fb = 0, fc2 = 0, fd = 0;
    for (int d = 0; ww * d <= s && !found; ++d)
      for (int c = 0; ww * d + wz * c <= s && !found; ++c)
        for (int b = 0; ww * d + wz * c + wy * b <= s && !found; ++b) {
          int r = s - ww * d - wz * c - wy * b;
          if (r % q == 0) {
            found = 1;
            fa = r / q;
            fb = b;
            fc2 = c;
            fd = d;
          }
        }
    if (found) {
      sz_rho[nb] = s;
      ea[nb] = fa;
      eb[nb] = fb;
      ec[nb] = fc2;
      ed[nb] = fd;
      ++nb;
    }
    sz_cnt[s] = nb;
  }
  sz_nb = 0;
  sz_val = calloc((size_t)nb + 1, sizeof(fe *));
  if (!sz_val) {
    free(ea);
    free(eb);
    free(ec);
    free(ed);
    return -1;
  }
  sz_nb = nb;
  for (int i = 0; i < nb; ++i) {
    sz_val[i] = malloc((size_t)sz_np * sizeof(fe));
    if (!sz_val[i]) {
      free(ea);
      free(eb);
      free(ec);
      free(ed);
      return -1;
    }
  }
  for (int k = 0; k < sz_np; ++k) {
    fe x = sz_px[k], y = sz_py[k];
    fe y2 = gpow(y, 2 * q0);
    fe z = (fe)(gpow(x, 2 * q0 + 1) ^ y2);
    fe w = (fe)(mul(x, y2) ^ gpow(z, 2 * q0));
    for (int i = 0; i < nb; ++i)
      sz_val[i][k] = mul(mul(gpow(x, ea[i]), gpow(y, eb[i])),
                         mul(gpow(z, ec[i]), gpow(w, ed[i])));
  }
  free(ea);
  free(eb);
  free(ec);
  free(ed);
  return 0;
}
/* flat n x n matrices, stride n */
static int sz_mat_inverse(const fe *m, fe *out, int n) {
  fe *wk = malloc((size_t)n * n * sizeof(fe));
  if (!wk)
    return 0;
  memcpy(wk, m, (size_t)n * n * sizeof(fe));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      out[(size_t)i * n + j] = 0;
    out[(size_t)i * n + i] = 1;
  }
  for (int c = 0; c < n; ++c) {
    int p = -1;
    for (int r = c; r < n; ++r)
      if (wk[(size_t)r * n + c]) {
        p = r;
        break;
      }
    if (p < 0) {
      free(wk);
      return 0;
    }
    if (p != c)
      for (int i = 0; i < n; ++i) {
        fe t = wk[(size_t)p * n + i];
        wk[(size_t)p * n + i] = wk[(size_t)c * n + i];
        wk[(size_t)c * n + i] = t;
        t = out[(size_t)p * n + i];
        out[(size_t)p * n + i] = out[(size_t)c * n + i];
        out[(size_t)c * n + i] = t;
      }
    fe iv = inv(wk[(size_t)c * n + c]);
    for (int i = 0; i < n; ++i) {
      wk[(size_t)c * n + i] = mul(wk[(size_t)c * n + i], iv);
      out[(size_t)c * n + i] = mul(out[(size_t)c * n + i], iv);
    }
    for (int r = 0; r < n; ++r) {
      fe f = wk[(size_t)r * n + c];
      if (r == c || !f)
        continue;
      for (int i = 0; i < n; ++i) {
        wk[(size_t)r * n + i] ^= mul(f, wk[(size_t)c * n + i]);
        out[(size_t)r * n + i] ^= mul(f, out[(size_t)c * n + i]);
      }
    }
  }
  free(wk);
  return 1;
}
/* characteristic polynomial det(theta I + H); H (n x n) is destroyed */
static int sz_charpoly(fe *H, int n, fe *g) {
  for (int k = 0; k + 2 < n; ++k) {
    int p = -1;
    for (int i = k + 1; i < n; ++i)
      if (H[(size_t)i * n + k]) {
        p = i;
        break;
      }
    if (p < 0)
      continue;
    if (p != k + 1) {
      for (int j = 0; j < n; ++j) {
        fe t = H[(size_t)p * n + j];
        H[(size_t)p * n + j] = H[(size_t)(k + 1) * n + j];
        H[(size_t)(k + 1) * n + j] = t;
      }
      for (int i = 0; i < n; ++i) {
        fe t = H[(size_t)i * n + p];
        H[(size_t)i * n + p] = H[(size_t)i * n + k + 1];
        H[(size_t)i * n + k + 1] = t;
      }
    }
    fe iv = inv(H[(size_t)(k + 1) * n + k]);
    for (int i = k + 2; i < n; ++i) {
      if (!H[(size_t)i * n + k])
        continue;
      fe f = mul(H[(size_t)i * n + k], iv);
      for (int j = 0; j < n; ++j)
        H[(size_t)i * n + j] ^= mul(f, H[(size_t)(k + 1) * n + j]);
      for (int r = 0; r < n; ++r)
        H[(size_t)r * n + k + 1] ^= mul(f, H[(size_t)r * n + i]);
    }
  }
  fe **pol = calloc((size_t)n + 1, sizeof(fe *));
  if (!pol)
    return 0;
  for (int k = 0; k <= n; ++k) {
    pol[k] = calloc((size_t)k + 1, sizeof(fe));
    if (!pol[k]) {
      for (int j = 0; j < k; ++j)
        free(pol[j]);
      free(pol);
      return 0;
    }
  }
  pol[0][0] = 1;
  for (int k = 1; k <= n; ++k) {
    fe prod = 1;
    for (int d = 0; d < k; ++d) {
      pol[k][d + 1] ^= pol[k - 1][d];
      pol[k][d] ^= mul(H[(size_t)(k - 1) * n + k - 1], pol[k - 1][d]);
    }
    for (int i = 1; i <= k - 1; ++i) {
      prod = mul(prod, H[(size_t)(k - i) * n + k - i - 1]);
      if (!prod)
        break;
      fe c = mul(H[(size_t)(k - i - 1) * n + k - 1], prod);
      if (!c)
        continue;
      for (int d = 0; d <= k - i - 1; ++d)
        pol[k][d] ^= mul(c, pol[k - i - 1][d]);
    }
  }
  for (int d = 0; d <= n; ++d)
    g[d] = pol[n][d];
  for (int k = 0; k <= n; ++k)
    free(pol[k]);
  free(pol);
  return 1;
}
/* decode one received word r (length sz_np); returns 1 when every injected
 * error was recovered with its value (true_pos / r are the simulated truth) */
static int sz_decode(const fe *r, const int *true_pos, int nerr, int verbose,
                     int list_values, int *nfound_out) {
  const int np = sz_np, nb = sz_nb, M = sz_M, q = sz_q;
  int ok = 0;
  int *nz = malloc((size_t)np * sizeof(int));
  int nzc = 0;
  for (int k = 0; k < np; ++k)
    if (r[k])
      nz[nzc++] = k;
  int nrows = sz_count_upto(sz_rowmax);
  int nl = sz_count_upto(sz_ldeg);
  if (nl > nrows)
    nl = nrows;
  fe **red = calloc((size_t)nrows + 1, sizeof(fe *));
  fe **coef = calloc((size_t)nrows + 1, sizeof(fe *));
  int *piv = malloc(((size_t)nrows + 1) * sizeof(int));
  int *rowc = malloc(((size_t)nb + 1) * sizeof(int));
  fe *pivinv = malloc(((size_t)nrows + 1) * sizeof(fe));
  fe *u = malloc(((size_t)nzc + 1) * sizeof(fe));
  const int maxt = sz_tau + 2;
  int *R = malloc((size_t)(maxt + 1) * sizeof(int)),
      *C = malloc((size_t)(maxt + 1) * sizeof(int));
  fe *xv = calloc((size_t)np, sizeof(fe));
  int *found = malloc((size_t)(maxt + 1) * sizeof(int));
  fe *errv = malloc((size_t)(maxt + 1) * sizeof(fe));
  fe *S0 = NULL, *SX = NULL, *S0i = NULL, *Q = NULL, *gp = NULL, *aug = NULL;
  int *roots = NULL;
  int nfound = 0, overflow = 0, resolved = 0, tt = 0;
  if (!nz || !red || !coef || !piv || !rowc || !pivinv || !u || !R || !C ||
      !xv || !found || !errv)
    goto done;
  for (int j = 0; j < nb; ++j)
    rowc[j] = -1;
  for (int i = 0; i < nrows; ++i) {
    int ncol = sz_count_upto(M - sz_rho[i]);
    red[i] = malloc((size_t)(ncol + 1) * sizeof(fe));
    coef[i] = calloc((size_t)i + 1, sizeof(fe));
    if (!red[i] || !coef[i])
      goto done;
    for (int n = 0; n < nzc; ++n)
      u[n] = mul(r[nz[n]], sz_val[i][nz[n]]);
    for (int j = 0; j < ncol; ++j) {
      fe acc = 0;
      for (int n = 0; n < nzc; ++n)
        acc ^= mul(u[n], sz_val[j][nz[n]]);
      red[i][j] = acc;
    }
    coef[i][i] = 1;
    piv[i] = -1;
    for (int j = 0; j < ncol; ++j) {
      fe v = red[i][j];
      if (!v)
        continue;
      int p = rowc[j];
      if (p >= 0) {
        fe f = mul(v, pivinv[p]);
        for (int jj = j; jj < ncol; ++jj)
          red[i][jj] ^= mul(f, red[p][jj]);
        for (int k = 0; k <= p; ++k)
          coef[i][k] ^= mul(f, coef[p][k]);
      } else {
        piv[i] = j;
        rowc[j] = i;
        pivinv[i] = inv(v);
        break;
      }
    }
  }
  int nloc = 0;
  for (int i = 0; i < nl; ++i)
    if (piv[i] < 0)
      ++nloc;
  for (int i = 0; i < nrows; ++i)
    if (piv[i] >= 0) {
      if (tt >= maxt) {
        if (verbose)
          puts("# SUZUKI: too many pivots (more errors than the designed "
               "bound)");
        goto done;
      }
      R[tt] = i;
      C[tt] = piv[i];
      ++tt;
    }
  if (verbose) {
    printf("# SUZUKI: %d error locators with rho <= %d\n", nloc, sz_ldeg);
    printf("# SUZUKI: staircase of %d elements -> estimated number of errors "
           "t = %d\n",
           tt, tt);
  }
  if (tt == 0) {
    if (verbose)
      puts("# SUZUKI: no errors.");
    ok = (nerr == 0);
    goto done;
  }
  for (int a = 0; a < tt; ++a)
    for (int b = 0; b < tt; ++b)
      if (sz_rho[R[a]] + sz_rho[C[b]] + q > M) {
        if (verbose)
          puts("# SUZUKI: syndromes too short to form the shifted matrix.");
        goto done;
      }
  for (int n = 0; n < nzc; ++n)
    xv[nz[n]] = sz_px[nz[n]];
  S0 = malloc((size_t)tt * tt * sizeof(fe));
  SX = malloc((size_t)tt * tt * sizeof(fe));
  S0i = malloc((size_t)tt * tt * sizeof(fe));
  Q = malloc((size_t)tt * tt * sizeof(fe));
  gp = malloc(((size_t)tt + 1) * sizeof(fe));
  if (!S0 || !SX || !S0i || !Q || !gp)
    goto done;
  for (int a = 0; a < tt; ++a) {
    for (int n = 0; n < nzc; ++n)
      u[n] = mul(r[nz[n]], sz_val[R[a]][nz[n]]);
    for (int b = 0; b < tt; ++b) {
      fe a0 = 0, ax = 0;
      for (int n = 0; n < nzc; ++n) {
        fe t = mul(u[n], sz_val[C[b]][nz[n]]);
        a0 ^= t;
        ax ^= mul(t, xv[nz[n]]);
      }
      S0[(size_t)a * tt + b] = a0;
      SX[(size_t)a * tt + b] = ax;
    }
  }
  if (!sz_mat_inverse(S0, S0i, tt)) {
    if (verbose)
      puts("# SUZUKI: S0 is singular (rank deficient).");
    goto done;
  }
  for (int a = 0; a < tt; ++a)
    for (int b = 0; b < tt; ++b) {
      fe acc = 0;
      for (int k = 0; k < tt; ++k)
        acc ^= mul(S0i[(size_t)a * tt + k], SX[(size_t)k * tt + b]);
      Q[(size_t)a * tt + b] = acc;
    }
  if (!sz_charpoly(Q, tt, gp))
    goto done;
  roots = malloc((size_t)FS * sizeof(int));
  if (!roots)
    goto done;
  int nroots = 0;
  for (int th = 0; th < FS; ++th) {
    fe acc = gp[tt];
    for (int i = tt - 1; i >= 0; --i)
      acc = (fe)(mul(acc, (fe)th) ^ gp[i]);
    if (!acc)
      roots[nroots++] = th;
  }
  if (verbose)
    printf(
        "# SUZUKI: g(theta) of degree %d has %d distinct roots in GF(2^%d)\n",
        tt, nroots, FM);
  for (int m = 0; m < nroots; ++m) {
    int xr = roots[m];
    for (int j = 0; j < sz_count[xr]; ++j) {
      int idx = sz_first[xr] + j, vanish = 1;
      for (int i = 0; i < nl && vanish; ++i) {
        if (piv[i] >= 0)
          continue;
        fe acc = 0;
        for (int k = 0; k <= i; ++k)
          if (coef[i][k])
            acc ^= mul(coef[i][k], sz_val[k][idx]);
        if (acc)
          vanish = 0;
      }
      if (vanish) {
        if (nfound < maxt)
          found[nfound++] = idx;
        else
          overflow = 1;
      }
    }
  }
  if (!overflow && nfound > 0) {
    if (nfound > tt && verbose)
      printf("# SUZUKI: %d common zeros > t = %d: resolving by error values\n",
             nfound, tt);
    int n = nfound, mrows = nb;
    aug = malloc((size_t)mrows * ((size_t)n + 1) * sizeof(fe));
    if (aug) {
      const size_t st = (size_t)n + 1;
      for (int a = 0; a < mrows; ++a) {
        for (int k = 0; k < n; ++k)
          aug[(size_t)a * st + k] = sz_val[a][found[k]];
        fe s = 0;
        for (int i2 = 0; i2 < nzc; ++i2)
          s ^= mul(r[nz[i2]], sz_val[a][nz[i2]]);
        aug[(size_t)a * st + n] = s;
      }
      int good = 1;
      for (int c = 0; c < n && good; ++c) {
        int p = -1;
        for (int a = c; a < mrows; ++a)
          if (aug[(size_t)a * st + c]) {
            p = a;
            break;
          }
        if (p < 0) {
          good = 0;
          break;
        }
        if (p != c)
          for (int k = 0; k <= n; ++k) {
            fe t = aug[(size_t)p * st + k];
            aug[(size_t)p * st + k] = aug[(size_t)c * st + k];
            aug[(size_t)c * st + k] = t;
          }
        fe iv = inv(aug[(size_t)c * st + c]);
        for (int k = c; k <= n; ++k)
          aug[(size_t)c * st + k] = mul(aug[(size_t)c * st + k], iv);
        for (int a = 0; a < mrows; ++a) {
          fe f = aug[(size_t)a * st + c];
          if (a == c || !f)
            continue;
          for (int k = c; k <= n; ++k)
            aug[(size_t)a * st + k] ^= mul(f, aug[(size_t)c * st + k]);
        }
      }
      if (good)
        for (int a = n; a < mrows; ++a)
          if (aug[(size_t)a * st + n]) {
            good = 0;
            break;
          }
      if (good) {
        int keep = 0;
        for (int k = 0; k < n; ++k)
          if (aug[(size_t)k * st + n]) {
            found[keep] = found[k];
            errv[keep] = aug[(size_t)k * st + n];
            ++keep;
          }
        nfound = keep;
        resolved = 1;
      }
    }
    if (!resolved && verbose)
      puts("# SUZUKI: could not compute the error values.");
  }
  if (nfound_out)
    *nfound_out = nfound;
  {
    int same = 0;
    for (int m = 0; m < nfound && resolved; ++m)
      if (r[found[m]] && errv[m] == r[found[m]])
        ++same;
    ok = resolved && !overflow && nfound == nerr && same == nerr;
    if (verbose) {
      if (list_values)
        for (int m = 0; m < nfound; ++m)
          printf("#   decoded error %d: position %d (x=%u,y=%u)  value=%u\n",
                 m + 1, found[m], sz_px[found[m]], sz_py[found[m]],
                 resolved ? errv[m] : 0);
      printf("# SUZUKI %d of %d errors recovered with the correct value, %d "
             "reported -> %s\n",
             same, nerr, nfound, ok ? "OK" : "MISMATCH");
    }
  }
done:
  if (red)
    for (int i = 0; i < nrows; ++i)
      free(red[i]);
  if (coef)
    for (int i = 0; i < nrows; ++i)
      free(coef[i]);
  free(red);
  free(coef);
  free(piv);
  free(rowc);
  free(pivinv);
  free(u);
  free(R);
  free(C);
  free(xv);
  free(found);
  free(errv);
  free(S0);
  free(SX);
  free(S0i);
  free(Q);
  free(gp);
  free(aug);
  free(roots);
  free(nz);
  (void)true_pos;
  return ok;
}
static int suzuki_run(int q_req, int nerr, int Mreq, int trials, int verbose,
                      int fixed_val, int list_values, int only_list) {
  int q = q_req > 0 ? q_req : 8, ex = 0;
  while ((1 << ex) < q)
    ++ex;
  if ((1 << ex) != q || (ex % 2) == 0 || ex > 9) {
    fprintf(
        stderr,
        "--curve suzuki needs --q = 2^(2n+1) (2, 8, 32, 128, 512); got %d\n",
        q);
    return 1;
  }
  sz_q = q;
  sz_q0 = 1 << ((ex - 1) / 2);
  sz_g = sz_q0 * (q - 1);
  if (nerr < 0 || trials < 1) {
    fprintf(stderr, "invalid --errors or -n\n");
    return 1;
  }
  sz_M = Mreq > 0 ? Mreq : 2 * nerr + 3 * sz_g - 1;
  if (sz_M < 3 * sz_g + 1)
    sz_M = 3 * sz_g + 1;
  sz_tau = (sz_M - 3 * sz_g + 1) / 2;
  sz_ldeg = sz_tau + sz_g;
  sz_rowmax = sz_M - q;
  printf("# Suzuki curve: y^%d + y = x^%d (x^%d + x)  over GF(2^%d)  "
         "(q=%d, q0=%d, genus g=%d)\n",
         q, sz_q0, q, FM, q, sz_q0, sz_g);
  if (sz_enumerate() != 0) {
    fprintf(stderr, "Suzuki point enumeration failed\n");
    sz_free_all();
    return 1;
  }
  printf(
      "# %d affine rational points (+1 at infinity); semigroup <%d,%d,%d,%d>\n",
      sz_np, q, q + sz_q0, q + 2 * sz_q0, q + 2 * sz_q0 + 1);
  if (suggest_param_mode) {
    print_suggest_params(sz_g, sz_np, nerr);
    sz_free_all();
    return 0;
  }
  printf("# one-point code: M = %d, g = %d, designed t <= %d, locators in "
         "L(%d P_inf)\n",
         sz_M, sz_g, sz_tau, sz_ldeg);
  if (only_list) {
    sz_free_all();
    return 0;
  }
  if (nerr > sz_tau) {
    fprintf(stderr,
            "# %d errors exceed the designed bound t <= %d for M = %d; raise "
            "--M (need M >= %d) or lower --errors\n",
            nerr, sz_tau, sz_M, 2 * nerr + 3 * sz_g - 1);
    sz_free_all();
    return 1;
  }
  if (sz_np <= sz_M + sz_g) {
    fprintf(stderr,
            "# the curve has only n = %d rational points, but M + g = %d >= "
            "n: use a larger field (--m) or smaller M/errors\n",
            sz_np, sz_M + sz_g);
    sz_free_all();
    return 1;
  }
  if (nerr > sz_np || sz_build_basis() != 0) {
    fprintf(stderr, "Suzuki basis setup failed\n");
    sz_free_all();
    return 1;
  }
  printf("# basis: %d monomials x^a y^b z^c w^d, one per semigroup element "
         "<= M\n",
         sz_nb);
  fe *r = malloc((size_t)sz_np * sizeof(fe));
  int *tp = malloc(((size_t)nerr + 1) * sizeof(int));
  unsigned char *used = malloc((size_t)sz_np);
  int good = 0;
  for (int tr = 1; tr <= trials; ++tr) {
    memset(r, 0, (size_t)sz_np * sizeof(fe));
    memset(used, 0, (size_t)sz_np);
    for (int e = 0; e < nerr; ++e) {
      int p;
      do
        p = (int)(rnd() % (uint64_t)sz_np);
      while (used[p]);
      used[p] = 1;
      tp[e] = p;
      r[p] =
          fixed_val > 0 ? (fe)fixed_val : (fe)(1 + rnd() % (uint64_t)(FS - 1));
    }
    int vb = (trials == 1) || verbose;
    if (vb && trials > 1)
      printf("# trial %d\n", tr);
    if (vb)
      printf("# %d errors injected\n", nerr);
    int ok = sz_decode(r, tp, nerr, vb, list_values, NULL);
    if (ok)
      ++good;
    if (!vb)
      printf("# trial %d: SUZUKI %s\n", tr, ok ? "ok" : "FAIL");
  }
  if (trials > 1)
    printf("# summary of %d trials with %d errors: SUZUKI %d/%d\n", trials,
           nerr, good, trials);
  free(r);
  free(tp);
  free(used);
  sz_free_all();
  return good == trials ? 0 : 3;
}

/* =========================================================================
 * Automatic normalisation of a user-supplied plane model (--terms).
 *
 * The decoder needs a model F(X,Y)=0 that (a) is monic in X and (b) has a
 * rational place whose Weierstrass semigroup the Newton-polygon engine can
 * certify (a primitive edge, lattice length 1).  Both are properties of the
 * COORDINATES, not of the curve.  Affine-plane automorphisms keep the
 * affine rational points (= code positions, up to relabelling) and the ring
 * of regular functions, so we search for one that makes the model supported:
 *    swap X<->Y;   Y -> Y + c X^k;   X -> X + c Y^k   (k >= 0, c in GF(2^m))
 * The shear Y -> Y + aX makes any curve monic in X (top form non-zero at a);
 * the substitutions with c a root of an edge polynomial split a place off a
 * long edge (this is a Tschirnhaus step: the place lands on a shorter edge).
 * Supported models are never touched.  Not covered: Newton-degenerate models
 * (tangent line at infinity, repeated edge roots), edge roots of degree > 1
 * with |ex|,|ey| >= 2, and curves with no rational point at all.
 * ========================================================================= */
typedef struct {
  int T;
  int X[MAXTERM], Y[MAXTERM];
  fe C[MAXTERM];
} Curve;
static void cv_load(Curve *c) {
  c->T = T;
  for (int t = 0; t < T; ++t) {
    c->X[t] = TX[t];
    c->Y[t] = TY[t];
    c->C[t] = TC[t];
  }
}
static void cv_store(const Curve *c) {
  T = c->T;
  for (int t = 0; t < c->T; ++t) {
    TX[t] = c->X[t];
    TY[t] = c->Y[t];
    TC[t] = c->C[t];
  }
}
static void cv_swap(const Curve *in, Curve *out) {
  out->T = in->T;
  for (int t = 0; t < in->T; ++t) {
    out->X[t] = in->Y[t];
    out->Y[t] = in->X[t];
    out->C[t] = in->C[t];
  }
}
static int cv_xmonic(const Curve *c) {
  int dx = 0, top = 0, pure = 0;
  for (int t = 0; t < c->T; ++t)
    if (c->X[t] > dx)
      dx = c->X[t];
  for (int t = 0; t < c->T; ++t)
    if (c->X[t] == dx) {
      top++;
      if (c->Y[t] == 0)
        pure = 1;
    }
  return top == 1 && pure;
}
/* onY=1:  Y -> Y + cc X^k      onY=0:  X -> X + cc Y^k     (binomials mod 2) */
static int cv_subst(const Curve *in, Curve *out, int onY, fe cc, int k) {
  enum { CAP = 4096 };
  static int X[CAP], Y[CAP];
  static fe C[CAP];
  int n = 0;
  for (int t = 0; t < in->T; ++t) {
    int i = in->X[t], j = in->Y[t], e = onY ? j : i;
    for (int l = 0; l <= e; ++l) {
      if ((e & l) != l)
        continue; /* C(e,l) is odd iff l is a bit-subset of e */
      fe cf = mul(in->C[t], gpow(cc, e - l));
      int ni = onY ? i + k * (e - l) : l;
      int nj = onY ? l : j + k * (e - l);
      int f = 0;
      for (int z = 0; z < n; ++z)
        if (X[z] == ni && Y[z] == nj) {
          C[z] ^= cf;
          f = 1;
          break;
        }
      if (!f) {
        if (n >= CAP)
          return 0;
        X[n] = ni;
        Y[n] = nj;
        C[n] = cf;
        n++;
      }
    }
  }
  out->T = 0;
  for (int z = 0; z < n; ++z)
    if (C[z]) {
      if (out->T >= MAXTERM)
        return 0;
      out->X[out->T] = X[z];
      out->Y[out->T] = Y[z];
      out->C[out->T] = C[z];
      out->T++;
    }
  return out->T >= 2;
}
/* can the existing Newton engine certify a one-point code on this model? */
static int cv_trial(const Curve *c) {
  Curve save;
  cv_load(&save);
  cv_store(c);
  g_trial = 1;
  g_trial_fail = 0;
  int ok = 0;
  if (cv_xmonic(c)) {
    int g = build_edges();
    if (!g_trial_fail) {
      genus_g = g;
      ne_inf_orbits();
      if (!g_trial_fail && use_newton)
        for (int k = 0; k < ne && !ok; ++k)
          if (ed[k].len == 1) {
            NBas nb;
            memset(&nb, 0, sizeof nb);
            int rcnb = newton_basis(k, 2 * g + 6, &nb);
            if (getenv("AG_CVDEBUG"))
              fprintf(stderr,
                      "[cv] edge %d len1 inf=%d w=(%d,%d) newton_basis=%d\n", k,
                      ed[k].inf, ed[k].wx, ed[k].wy, rcnb);
            if (rcnb == 0) {
              ok = 1;
              free(nb.rho);
              free(nb.ei);
              free(nb.ej);
            }
          }
    }
  }
  if (getenv("AG_CVDEBUG")) {
    fprintf(stderr, "[cv] trial ok=%d monic=%d g=%d terms:", ok, cv_xmonic(c),
            genus_g);
    for (int t = 0; t < c->T; ++t)
      fprintf(stderr, " %d,%d,%u", c->X[t], c->Y[t], (unsigned)c->C[t]);
    fprintf(stderr, "\n");
    for (int k = 0; k < ne; ++k)
      fprintf(stderr, "[cv]   E%d v=(%d,%d) e=(%d,%d) len=%d inf=%d\n", k,
              ed[k].vx, ed[k].vy, ed[k].ex, ed[k].ey, ed[k].len, ed[k].inf);
  }
  g_trial = 0;
  cv_store(&save);
  return ok;
}
typedef struct {
  int kind; /* 0 swap, 1 Y->Y+cX^k, 2 X->X+cY^k */
  fe c;
  int k;
} CVStep;
static CVStep cvs[10];
static int cvn, cv_budget;
static int cv_search(const Curve *cur, int depth, int need_place, Curve *res) {
  int monic = cv_xmonic(cur);
  if (monic && (!need_place || cv_trial(cur))) {
    *res = *cur;
    return 1;
  }
  if (depth == 0 || --cv_budget < 0)
    return 0;
  Curve nx;
  if (!monic) {
    cv_swap(cur, &nx);
    if (cv_xmonic(&nx)) {
      cvs[cvn++] = (CVStep){0, 0, 0};
      if (cv_search(&nx, depth - 1, need_place, res))
        return 1;
      cvn--;
    }
    int tried = 0;
    for (int a = 1; a < FS && tried < 2; ++a)
      if (cv_subst(cur, &nx, 1, (fe)a, 1) && cv_xmonic(&nx)) {
        ++tried;
        cvs[cvn++] = (CVStep){1, (fe)a, 1};
        if (cv_search(&nx, depth - 1, need_place, res))
          return 1;
        cvn--;
      }
    return 0;
  }
  /* monic but no certified place: split places off long edges */
  Curve save;
  cv_load(&save);
  cv_store(cur);
  g_trial = 1;
  g_trial_fail = 0;
  int g = build_edges();
  typedef struct {
    int onY, k;
    fe c;
  } CVCandidate;
  CVCandidate cand[64];
  int nc = 0;
  if (!g_trial_fail) {
    genus_g = g;
    for (int k = 0; k < ne && nc < 60; ++k) {
      if (!ed[k].inf || ed[k].len < 2)
        continue;
      Poly *facs = malloc(sizeof(*facs) * (PMAXD + 1));
      if (!facs)
        continue;
      int nf = split_edge(k, facs);
      for (int f = 0; f < nf && !g_trial_fail && nc < 60; ++f) {
        if (facs[f].deg != 1 || !facs[f].c[0])
          continue;
        fe u0 = inv(facs[f].c[0]); /* x^ex y^ey = u0 at this place */
        int ex = ed[k].ex, ey = ed[k].ey;
        if (ey == 1 && ex <= 0)
          cand[nc++] = (CVCandidate){1, -ex, u0};
        else if (ey == -1 && ex >= 0)
          cand[nc++] = (CVCandidate){1, ex, inv(u0)};
        else if (ex == 1 && ey <= 0)
          cand[nc++] = (CVCandidate){0, -ey, u0};
        else if (ex == -1 && ey >= 0)
          cand[nc++] = (CVCandidate){0, ey, inv(u0)};
      }
      free(facs);
    }
  }
  g_trial = 0;
  cv_store(&save);
  for (int z = 0; z < nc; ++z)
    if (cv_subst(cur, &nx, cand[z].onY, cand[z].c, cand[z].k)) {
      cvs[cvn++] = (CVStep){cand[z].onY ? 1 : 2, cand[z].c, cand[z].k};
      if (cv_search(&nx, depth - 1, need_place, res))
        return 1;
      cvn--;
    }
  return 0;
}
/* returns 1 if the model was replaced by an equivalent supported one */
static int cv_normalize(int need_place) {
  Curve orig, res;
  cv_load(&orig);
  if (orig.T < 2)
    return 0;
  cvn = 0;
  cv_budget = 400;
  if (cv_xmonic(&orig) && (!need_place || cv_trial(&orig)))
    return 0; /* already supported: leave untouched */
  if (!cv_search(&orig, 5, need_place, &res)) {
    /* at least try to reach a monic model, so the usual diagnostics apply */
    cvn = 0;
    cv_budget = 400;
    if (cv_xmonic(&orig) || !cv_search(&orig, 3, 0, &res))
      return 0;
  }
  printf("# curve normalisation: the given model is not directly supported by "
         "the Newton engine;\n#   using an equivalent model (affine-plane "
         "automorphism, same rational points up to relabelling):\n");
  for (int z = 0; z < cvn; ++z) {
    if (cvs[z].kind == 0)
      printf("#   swap X <-> Y\n");
    else
      printf("#   %s -> %s + %u * %s^%d   (c as field element, integer "
             "representation)\n",
             cvs[z].kind == 1 ? "Y" : "X", cvs[z].kind == 1 ? "Y" : "X",
             (unsigned)cvs[z].c, cvs[z].kind == 1 ? "X" : "Y", cvs[z].k);
  }
  printf("#   new model (i,j,coef):");
  for (int t = 0; t < res.T; ++t)
    printf(" %d,%d,%u;", res.X[t], res.Y[t], (unsigned)res.C[t]);
  printf("\n");
  cv_store(&res);
  return 1;
}

static int any_newton_place(int g) {
  for (int k = 0; k < ne; ++k)
    if (ed[k].len == 1) {
      NBas nb;
      memset(&nb, 0, sizeof nb);
      if (newton_basis(k, 2 * g + 6, &nb) == 0) {
        free(nb.rho);
        free(nb.ei);
        free(nb.ej);
        return 1;
      }
    }
  return 0;
}

/* -------------------------------------------------------------------------
 * Error-correcting-pair decoder for one-point evaluation codes.
 * C = ev(L(MP)); A = ev(L(aP)), a=t+g; B=(A*C)^perp.
 * Since A*C is contained in ev(L((M+a)P)), the usual Goppa bounds certify
 * d(A)>=n-a, d(C)>=n-M and d(B^perp)>=n-M-a.  The ECP inequalities are
 * checked before decoding.  This path is deliberately independent of FR/BMS.
 * The one-point simulator transmits the zero codeword, so recovered error
 * values are the received symbols at the ECP common-zero positions.
 */
static int ecp_fit_support(int n, int k, const fe *G, const fe *y,
                           const int *support, int ns, fe *err) {
  int rows=n-ns, W=k+1;
  if (rows<k) return 0;
  fe *A=calloc((size_t)rows*W,sizeof(fe));
  int *pivot=malloc((size_t)k*sizeof(int));
  if (!A || !pivot) { free(A); free(pivot); return 0; }
  int rr=0;
  for (int p=0;p<n;p++) {
    int skip=0; for (int q=0;q<ns;q++) if (support[q]==p) {skip=1;break;}
    if (skip) continue;
    for (int j=0;j<k;j++) A[(size_t)rr*W+j]=G[(size_t)j*n+p];
    A[(size_t)rr*W+k]=y[p]; rr++;
  }
  int rank=0;
  for (int c=0;c<k;c++) {
    int p=rank; while (p<rows && !A[(size_t)p*W+c]) p++;
    if (p==rows) { free(A); free(pivot); return 0; }
    if (p!=rank) for (int j=c;j<W;j++) {fe z=A[(size_t)p*W+j];A[(size_t)p*W+j]=A[(size_t)rank*W+j];A[(size_t)rank*W+j]=z;}
    fe iv=inv(A[(size_t)rank*W+c]);
    for (int j=c;j<W;j++) A[(size_t)rank*W+j]=mul(A[(size_t)rank*W+j],iv);
    for (int q=0;q<rows;q++) if(q!=rank) {
      fe f=A[(size_t)q*W+c]; if(!f) continue;
      for(int j=c;j<W;j++) A[(size_t)q*W+j]^=mul(f,A[(size_t)rank*W+j]);
    }
    pivot[rank++]=c;
  }
  for(int q=k;q<rows;q++) if(A[(size_t)q*W+k]) {free(A);free(pivot);return 0;}
  fe *msg=calloc((size_t)k,sizeof(fe));
  if(!msg){free(A);free(pivot);return 0;}
  for(int q=0;q<k;q++) msg[pivot[q]]=A[(size_t)q*W+k];
  for(int q=0;q<ns;q++) {
    int pos=support[q]; fe cw=0;
    for(int j=0;j<k;j++) if(msg[j]) cw^=mul(msg[j],G[(size_t)j*n+pos]);
    err[q]=y[pos]^cw;
    if(!err[q]) {free(msg);free(A);free(pivot);return 0;}
  }
  free(msg);free(A);free(pivot);return 1;
}
static int ecp_support_search(int n,int k,const fe *G,const fe *y,
                              const int *cand,int nc,int t,int start,int depth,
                              int *work,int *best,fe *workv,fe *bestv,
                              unsigned long *tries) {
  if(depth==t) {
    if(++*tries>100000UL) return 0;
    if(ecp_fit_support(n,k,G,y,work,t,workv)) {
      memcpy(best,work,(size_t)t*sizeof(int));
      memcpy(bestv,workv,(size_t)t*sizeof(fe)); return 1;
    }
    return 0;
  }
  for(int i=start;i<=nc-(t-depth);i++) {
    work[depth]=cand[i];
    if(ecp_support_search(n,k,G,y,cand,nc,t,i+1,depth+1,work,best,workv,bestv,tries)) return 1;
    if(*tries>100000UL) return 0;
  }
  return 0;
}

static int ecp_ag_decode(int n, int M, int g, int t, const fe *received,
                         int *out_pos, fe *out_val, int *out_n) {
  if (!received || !out_pos || !out_val || !out_n || t < 0)
    return 0;
  int a = t + g;
  if (M < 2*g-1 || a + M >= n || n - M - a <= t) {
    printf("# ECP: certified Goppa bounds do not fit (need M+a<n and "
           "n-M-a>t; M=%d, a=%d, n=%d, t=%d)\n", M, a, n, t);
    return 0;
  }
  int dc = count_basis_upto(M), da = count_basis_upto(a);
  if (dc != M + 1 - g || da != a + 1 - g || da <= t || dc <= 0) {
    printf("# ECP: Newton basis dimensions do not match Riemann-Roch "
           "(dim C=%d expected %d; dim A=%d expected %d)\n",
           dc, M+1-g, da, a+1-g);
    return 0;
  }
  if ((size_t)da * (size_t)dc * (size_t)n > 20000000u) {
    puts("# ECP: product-space row cap exceeded in this build");
    return 0;
  }
  fe *P = calloc((size_t)da * dc * n, sizeof(fe));
  if (!P) return 0;
  for (int i=0;i<da;i++) for (int j=0;j<dc;j++) {
    fe *row=P+((size_t)i*dc+j)*n;
    if (!basis[i].val || !basis[j].val) { free(P); return 0; }
    for (int x=0;x<n;x++) row[x]=mul(basis[i].val[x],basis[j].val[x]);
  }
  fe *B=NULL;
  int db=nullspace_gf(P, da*dc, n, &B);
  free(P);
  if (db <= 0 || !B) { free(B); puts("# ECP: B construction produced no nonzero space"); return 0; }
  fe *S=calloc((size_t)db*da,sizeof(fe));
  if (!S) { free(B); return 0; }
  for (int b=0;b<db;b++) for (int i=0;i<da;i++) {
    fe acc=0;
    for (int x=0;x<n;x++) {
      fe z=mul(mul(B[(size_t)b*n+x],received[x]),basis[i].val[x]);
      acc ^= z;
    }
    S[(size_t)b*da+i]=acc;
  }
  fe *K=NULL;
  int dk=nullspace_gf(S,db,da,&K);
  free(S); free(B);
  if (dk <= 0 || !K) { free(K); puts("# ECP: locator kernel is trivial"); return 0; }
  int nz=0, *cand=malloc(sizeof(int)*(size_t)n);
  for (int x=0;x<n;x++) {
    int common=1;
    for (int q=0;q<dk && common;q++) {
      fe z=0;
      for (int i=0;i<da;i++) if (K[(size_t)q*da+i])
        z ^= mul(K[(size_t)q*da+i],basis[i].val[x]);
      if (z) common=0;
    }
    if (common) cand[nz++]=x;
  }
  free(K);
  if (nz<t || nz>MAXT) { free(cand); printf("# ECP: common-zero candidate set has %d positions for t=%d\n",nz,t); return 0; }
  fe *C=calloc((size_t)dc*n,sizeof(fe));
  int *work=malloc(sizeof(int)*(size_t)t), *best=malloc(sizeof(int)*(size_t)t);
  fe *wv=malloc(sizeof(fe)*(size_t)t), *bv=malloc(sizeof(fe)*(size_t)t);
  if(!C||!work||!best||!wv||!bv){free(cand);free(C);free(work);free(best);free(wv);free(bv);return 0;}
  for(int j=0;j<dc;j++) memcpy(C+(size_t)j*n,basis[j].val,sizeof(fe)*(size_t)n);
  unsigned long tries=0;
  int found=ecp_support_search(n,dc,C,received,cand,nz,t,0,0,work,best,wv,bv,&tries);
  free(cand);free(C);free(work);free(wv);
  if(!found){free(best);free(bv);printf("# ECP: no support passed the codeword consistency check (%lu subsets tested)\n",tries);return 0;}
  memcpy(out_pos,best,sizeof(int)*(size_t)t); memcpy(out_val,bv,sizeof(fe)*(size_t)t);
  free(best);free(bv); *out_n=t;
  printf("# ECP: A=L(%dP) dim=%d, C=L(%dP) dim=%d, B dim=%d; "
         "bounds d(A)>=%d, d(C)>=%d, d(B^perp)>%d\n",
         a,da,M,dc,db,n-a,n-M,t);
  return 1;
}

static int ag_concat_run(const fe *Gen, const fe *VG, int npts, int kcode,
                         int rG, int outer_m, const char *inner_spec,
                         unsigned poly, int trials, int nerr, int verbose,
                         int showv, int burst, int interleave);
static int ag_concat_frvote_run(int npts, int kcode, int outer_m,
                                const char *inner_spec, unsigned poly,
                                int trials, int nerr, int verbose,
                                int showv, int burst, int interleave);

int noize_run(int argc, char **argv) {
  const char *terms = "16,0;0,1;15,4;12,16", *polespec = NULL,
             *concat_inner = NULL;
  const char *curve_name = "bbgs";
  int len_req = -1, half = 0, q_req = -1, terms_set = 0, m_bits = 6, Mreq = -1,
      bbgs_a = 1, bbgs_b = 2, place_req = -1, trials = 1, verbose = 0,
      show_values = 0, only_list = 0, Jstep = 16, Jmax = 700;
  unsigned poly = 0, seed = 1;
  int run_fr = 1, run_bms = 1, run_vote = 1, run_frvote = 0, only = 0,
      nerr = 256, force = 1, fixed_value_req = 0, run_common_binary = 0,
      encode_requested = 0, run_ecp = 0, concat_burst = 0,
      concat_interleave = 1;
  for (int i = 1; i < argc; ++i) {
    { /* accept --force / --fr / --bms / --vote / --basic / --list as well as
         -force ... */
      static const char *flags[] = {
          "half",   "force",  "fr",     "bms",      "vote",   "frvote",
          "basic",  "list",   "v",      "nonewton", "strict", "limits",
          "nopoly", "forney", "kotter", "encode", "ecp", NULL};
      if (argv[i][0] == '-' && argv[i][1] == '-')
        for (int f = 0; flags[f]; ++f)
          if (!strcmp(argv[i] + 2, flags[f])) {
            argv[i]++;
            break;
          }
    }
    if (!strcmp(argv[i], "--curve") && i + 1 < argc)
      curve_name = argv[++i];
    else if (!strcmp(argv[i], "--m") && i + 1 < argc)
      m_bits = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--poly") && i + 1 < argc)
      poly = (unsigned)strtoul(argv[++i], NULL, 0);
    else if (!strcmp(argv[i], "--terms") && i + 1 < argc) {
      terms = argv[++i];
      terms_set = 1;
    } else if (!strcmp(argv[i], "--len") && i + 1 < argc)
      len_req = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--q") && i + 1 < argc) {
      ++i;
      q_req = strcmp(argv[i], "auto") ? atoi(argv[i]) : -1;
    } else if (!strcmp(argv[i], "--a") && i + 1 < argc)
      bbgs_a = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--b") && i + 1 < argc)
      bbgs_b = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--errors") && i + 1 < argc)
      nerr = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--suggest-param"))
      suggest_param_mode = 1;
    else if (!strcmp(argv[i], "--fixed-error-value") && i + 1 < argc)
      fixed_value_req = (int)strtol(argv[++i], NULL, 0);
    else if (!strcmp(argv[i], "-common-binary"))
      run_common_binary = 1;
    else if (!strcmp(argv[i], "-bms2-check"))
      bms2_check_enabled = 1;
    else if (!strcmp(argv[i], "--M") && i + 1 < argc)
      Mreq = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--place") && i + 1 < argc)
      place_req = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--pole") && i + 1 < argc)
      polespec = argv[++i];
    else if (!strcmp(argv[i], "--basis") && i + 1 < argc) {
      const char *b = argv[++i];
      if (!strcmp(b, "legacy")) multi_newton_order = 0;
      else if (!strcmp(b, "newton")) multi_newton_order = 1;
      else { fprintf(stderr, "--basis must be legacy or newton\n"); return 2; }
    }
    else if (!strcmp(argv[i], "--concat-inner") && i + 1 < argc) {
      concat_inner = argv[++i];
      encode_requested = 1;
    } else if (!strcmp(argv[i], "--concat-burst")) {
      concat_burst = 1;
    } else if (!strcmp(argv[i], "--interleave") && i + 1 < argc) {
      concat_interleave = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "--Jstep") && i + 1 < argc)
      Jstep = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--Jmax") && i + 1 < argc)
      Jmax = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-n") && i + 1 < argc)
      trials = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-v"))
      verbose = 1;
    else if (!strcmp(argv[i], "-show-values")) {
      show_values = 1;
      show_locator_equations = 1;
    } else if (!strcmp(argv[i], "-list"))
      only_list = 1;
    else if (!strcmp(argv[i], "-force"))
      force = 1;
    else if (!strcmp(argv[i], "-strict"))
      force = 0;
    else if (!strcmp(argv[i], "-nopoly"))
      use_poly_opt = 0;
    else if (!strcmp(argv[i], "-forney"))
      use_forney_values = 1;
    else if (!strcmp(argv[i], "-kotter"))
      use_kotter_values = 1;
    else if (!strcmp(argv[i], "-encode"))
      encode_requested = 1;
    else if (!strcmp(argv[i], "-limits"))
      only_list = 1;
    else if (!strcmp(argv[i], "-nonewton")) {
      fprintf(stderr, "-nonewton is disabled: this build uses Newton-polygon "
                      "bases only.\n");
      return 2;
    } else if (!strcmp(argv[i], "-half"))
      half = 1;
    else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      usage();
      return 0;
    } else if (!strcmp(argv[i], "-fr") || !strcmp(argv[i], "-bms") ||
               !strcmp(argv[i], "-vote") || !strcmp(argv[i], "-frvote") ||
               !strcmp(argv[i], "-basic") || !strcmp(argv[i], "-ecp")) {
      if (!only)
        run_fr = run_bms = run_vote = run_frvote = 0;
      only = 1;
      if (!strcmp(argv[i], "-fr"))
        run_fr = 1;
      else if (!strcmp(argv[i], "-bms"))
        run_bms = 1;
      else if (!strcmp(argv[i], "-vote"))
        run_vote = 1;
      else if (!strcmp(argv[i], "-frvote"))
        run_frvote = 1;
      else if (!strcmp(argv[i], "-ecp"))
        run_ecp = 1;
    } else if (argv[i][0] == '-') {
      fprintf(stderr, "unknown option: %s  (see --help)\n", argv[i]);
      return 1;
    } else
      seed = (unsigned)strtoul(argv[i], NULL, 10);
  }
  /* The one-point AG/Goppa concatenation is decoded as one outer AG
   * syndrome problem after the binary inner decoders.  Select FR-VOTE
   * automatically; the legacy multi-point concatenation path remains separate. */
  if (concat_inner && !polespec) {
    run_fr = run_bms = run_vote = 0;
    run_frvote = 1;
  }
  if (run_frvote && run_fr) {
    fprintf(stderr, "-frvote cannot be combined with -fr: FR-VOTE keeps no "
                    "basis-value table and only the minimal locator rows.\n");
    return 1;
  }
  compact_fr_mode = run_frvote;
  if (api_mode &&
      (polespec != NULL || run_ecp || !run_frvote || run_fr || run_bms || run_vote)) {
    fprintf(stderr, "API mode supports only one-point Newton FR-VOTE; use "
                    "-frvote alone and do not set --pole.\n");
    return 1;
  }
  if (m_bits < 1 || m_bits > 12) {
    fprintf(stderr, "--m must be in 1..12 in this build\n");
    return 1;
  }
  if (!poly) {
    static const unsigned dp[13] = {0,     3,     7,     0xB,   0x13,
                                    0x25,  0x43,  0x89,  0x11D, 0x211,
                                    0x409, 0x805, 0x1053};
    poly = dp[m_bits];
  }
  if (gf_init(m_bits, poly) != 0) {
    fprintf(stderr, "polynomial 0x%X is not primitive of degree %d\n", poly,
            m_bits);
    return 1;
  }
  rng_s ^= (uint64_t)seed * 0x9E3779B97F4A7C15ULL;
  for (int i = 0; i < 8; ++i)
    rnd();
  srand(seed);
  if (!terms_set && !strcmp(curve_name, "suzuki")) {
    if (api_mode) {
      fprintf(stderr, "--curve suzuki is not available in API mode\n");
      return 1;
    }
    return suzuki_run(q_req, nerr, Mreq, trials, verbose, fixed_value_req,
                      show_values, only_list);
  }

  /* ---- curve family ----
   * F(X,Y)=0 is stored as a sparse polynomial.  Both built-in families are
   * handled by the same point enumerator, Newton-polygon code and FR-VOTE.
   */
  static char tbuf[4096];
  int qfam = 0;
  if (!terms_set) {
    int q = q_req;
    if (!strcmp(curve_name, "bbgs")) {
      int a = bbgs_a, b = bbgs_b, c = a + b;
      if (a < 1 || b < 1 || gcdi(a, b) != 1 || c > 10) {
        /* EXTENSION_POINT: c<=10 is a structural guard; revisit together
         * with qpow[] sizing and curve/point-count complexity. */
        fprintf(stderr,
                "BBGS requires positive coprime --a and --b, with a+b <= 10\n");
        return 1;
      }
      if (q < 0) {
        if (m_bits % c != 0) {
          fprintf(stderr, "For --a %d --b %d, --m must be divisible by c=%d\n",
                  a, b, c);
          return 1;
        }
        int r = m_bits / c;
        if (r < 1 || r > 10)
          return 1;
        q = 1 << r;
      }
      if (q < 2 || (q & (q - 1))) {
        fprintf(stderr, "--q must be a power of 2 (>= 2)\n");
        return 1;
      }
      /* EXTENSION_POINT: fixed exponent table; current m<=12 makes every
       * accepted q^c <= 512, but retain a checked bound for future growth. */
      int qpow[16];
      qpow[0] = 1;
      for (int k = 1; k <= c; ++k) {
        if (qpow[k - 1] > 8192 / q) { /* EXTENSION_POINT: exponent ceiling */
          fprintf(stderr, "BBGS exponents exceed this build's limits\n");
          return 1;
        }
        qpow[k] = qpow[k - 1] * q;
      }
      if (qpow[c] != FS) {
        int r = 0, qq = q;
        while (qq > 1) {
          ++r;
          qq >>= 1;
        }
        fprintf(stderr,
                "BBGS field mismatch: need --m %d for GF(q^c)=GF(2^%d), got "
                "--m %d\n",
                r * c, r * c, m_bits);
        return 1;
      }
      int used = snprintf(tbuf, sizeof tbuf, "%d,0", qpow[c - 1]);
      for (int k = 0; k < b; ++k)
        used += snprintf(tbuf + used, sizeof tbuf - (size_t)used, ";%d,%d",
                         qpow[c - 1] - qpow[k], qpow[a + k]);
      for (int k = 0; k < a; ++k)
        used += snprintf(tbuf + used, sizeof tbuf - (size_t)used, ";%d,%d",
                         qpow[c - 1] - qpow[b + k], qpow[k]);
      if (used < 0 || (size_t)used >= sizeof tbuf) {
        fprintf(stderr, "BBGS term list exceeds buffer\n");
        return 1;
      }
      terms = tbuf;
      if (a == 1 && b == 2)
        qfam = q;
      printf("# BBGS curve (a=%d,b=%d,c=%d,q=%d), %d terms:\n#   X^%d", a, b, c,
             q, c + 1, qpow[c - 1]);
      for (int k = 0; k < b; ++k)
        printf(" + X^%d Y^%d", qpow[c - 1] - qpow[k], qpow[a + k]);
      for (int k = 0; k < a; ++k)
        printf(" + X^%d Y^%d", qpow[c - 1] - qpow[b + k], qpow[k]);
      puts(" = 0");
    } else if (!strcmp(curve_name, "hermitian")) {
      if (q < 0)
        q = (m_bits % 3 == 0) ? (1 << (m_bits / 3)) : 4;
      if (q < 2 || (q & (q - 1))) {
        fprintf(stderr, "--q must be a power of 2 (>= 2)\n");
        return 1;
      }
      snprintf(tbuf, sizeof tbuf, "0,%d;%d,0;0,1", q, q + 1);
      terms = tbuf;
      hermitian_q = q;
      printf("# Hermitian curve: Y^%d + Y + X^%d = 0 (q=%d)\n", q, q + 1, q);
    } else if (!strcmp(curve_name, "klein")) {
      /* Klein quartic X^3 Y + Y^3 Z + Z^3 X = 0.
       * In the affine chart Z=1 this is x^3 y + y^3 + x = 0.
       * The common curve engine performs its usual model normalisation. */
      terms = "3,1;0,3;1,0";
      terms_set = 1; /* run the custom-model normaliser for this affine form */
      printf("# Klein quartic: X^3 Y + Y^3 Z + Z^3 X = 0\n"
             "# affine chart Z=1: x^3 y + y^3 + x = 0\n");
    } else {
      fprintf(stderr,
              "--curve must be bbgs, hermitian, suzuki or klein (or use --terms)\n");
      return 1;
    }
  }
  {
    char *s = strdup(terms);
    if (!s) {
      fprintf(stderr, "out of memory while parsing --terms\n");
      return 1;
    }
    char *tok = strtok(s, ";");
    while (tok) {
      int a, b, c = 1;
      int k = sscanf(tok, "%d,%d,%d", &a, &b, &c);
      if (k < 2) {
        fprintf(stderr, "--terms: \"i,j[,coef];...\"\n");
        free(s);
        return 1;
      }
      if (c <= 0 || c >= FS) {
        fprintf(stderr, "coefficient out of range\n");
        free(s);
        return 1;
      }
      if (T >= MAXTERM) {
        fprintf(stderr, "too many curve terms (maximum %d)\n", MAXTERM);
        free(s);
        return 1;
      }
      TX[T] = a;
      TY[T] = b;
      TC[T] = (fe)c;
      T++;
      tok = strtok(NULL, ";");
    }
    free(s);
  }
  if (terms_set && !hermitian_q)
    cv_normalize(polespec == NULL && place_req < 0 && use_newton);
  dx_deg = 0;
  for (int t = 0; t < T; ++t)
    if (TX[t] > dx_deg)
      dx_deg = TX[t];
  /* --terms that spell the Hermitian curve X^(q+1) + Y^q + Y (all coefficients
   * 1) get the same dedicated polynomial-semigroup path as --curve hermitian */
  if (hermitian_q == 0 && T == 3) {
    int q = 0, hx = 0, hy = 0, h1 = 0, ok = 1;
    for (int t = 0; t < T; ++t)
      if (TC[t] != 1)
        ok = 0;
    for (int t = 0; t < T && ok; ++t) {
      if (TY[t] == 0 && TX[t] > 1) {
        hx++;
        q = TX[t] - 1;
      } else if (TX[t] == 0 && TY[t] == 1)
        h1++;
      else if (TX[t] == 0 && TY[t] > 1)
        hy++;
    }
    if (ok && hx == 1 && h1 == 1 && hy == 1 && q >= 2 && !(q & (q - 1))) {
      int ycnt = 0;
      for (int t = 0; t < T; ++t)
        if (TX[t] == 0 && TY[t] == q)
          ycnt++;
      if (ycnt == 1) {
        hermitian_q = q;
        printf("# --terms recognised as the Hermitian curve (q=%d)\n", q);
      }
    }
  }
  {
    int top = 0, pure = 0;
    for (int t = 0; t < T; ++t)
      if (TX[t] == dx_deg) {
        top++;
        if (TY[t] == 0)
          pure = 1;
      }
    if (top != 1 || !pure) {
      fprintf(stderr, "only curves that are monic in X are supported (swap X "
                      "and Y in --terms)\n");
      return 1;
    }
  }
  genus_g = build_edges();
  int ng_orb = ne_inf_orbits();
  GENUS = genus_g;
  double ts = now();
  enumerate_points();
  int nfull = npts;
  if (half)
    len_req = npts / 2;
  else if (len_req < 0 && npts > 8192) {
    len_req = 8192; /* default code length cap (override with --len) */
    printf("# code length capped at 8192 (use --len L to change)\n");
  }
  if (len_req > 0 && len_req < npts)
    npts = len_req; /* puncturing: keep the first points only */
  printf("# curve over GF(2^%d) [poly 0x%X]: X-degree %d, genus g = %d, n = %d "
         "affine rational points\n",
         FM, FM ? poly : 0, dx_deg, genus_g, nfull);
  if (npts != nfull)
    printf("# punctured: only the first %d of %d points are used as code "
           "positions\n",
           npts, nfull);
  for (int k = 0; k < ne; ++k)
    printf("#   edge E%d: from (%d,%d) dir (%d,%d) length %d  w=(%d,%d) %s\n",
           k, ed[k].vx, ed[k].vy, ed[k].ex, ed[k].ey, ed[k].len, ed[k].wx,
           ed[k].wy, ed[k].inf ? "[infinity]" : "[affine side]");
  int nrat = 0, ratidx[MAXORB];
  for (int q = 0; q < ng_orb; ++q) {
    printf("#   place #%d: edge E%d, degree %d (%s)\n", q, pl[q].edge, pl[q].d,
           pl[q].d == 1 ? "rational" : "Galois orbit");
    if (pl[q].d == 1)
      ratidx[nrat++] = q;
  }
  if (qfam) {
    long long nopt = 1 + ((long long)qfam * qfam * qfam - 1) * qfam * qfam;
    if ((long long)qfam * qfam * qfam == (long long)FS)
      printf("# GF(2^%d) = GF(q^3): maximal-structure point count 1 + (q^3-1) "
             "q^2 = %lld, found n = %d%s\n",
             FM, nopt, nfull, nopt == nfull ? "  (matches)" : "  (differs)");
    else
      printf("# note: 2^%d != q^3, so the q = %d curve has no special point "
             "structure here\n",
             FM, qfam);
  }
  print_limits(GENUS, npts, concat_inner ? 0 : nerr);
  if (suggest_param_mode) {
    print_suggest_params(GENUS, npts, nerr);
    return 0;
  }
  if (only_list)
    return 0;
  /*
  if (GENUS > MAXGENUS - 1 || 2 * GENUS + 2 > MAXB - GENUS) {
      fprintf(stderr, "# genus g = %d: decoding needs M >= 2g = %d syndromes (>
  this build's limit: g <= %d).\n"
                      "# The dense linear algebra of this program is not
  feasible at this size; -list works, decoding does not.\n", GENUS, 2 * GENUS,
  MAXGENUS - 1); return 1;
  }
  */

  NUM_ERRORS = nerr;
  if (fixed_value_req < 0 || fixed_value_req >= FS) {
    fprintf(stderr, "--fixed-error-value must be in 1..%d (or omit it)\n",
            FS - 1);
    return 1;
  }
  fixed_error_value = fixed_value_req;
  if (run_common_binary) {
    run_fr = run_bms = run_vote = run_frvote = 0;
  }
  if ((concat_inner ? nerr < 0 : nerr < 1) ||
      (!concat_inner && nerr > MAXT - 2)) {
    if (concat_inner)
      fprintf(stderr, "--errors must be nonnegative for concatenation\n");
    else
      fprintf(stderr, "--errors out of range (1..%d)\n", MAXT - 2);
    return 1;
  }
  int mult[MAXORB];
  memset(mult, 0, sizeof mult);

  /* ===================== choose the code ===================== */
  /* Multi-point construction is opt-in only.  The default path is the
   * one-point Newton-polygon basis; never silently switch algorithms merely
   * because this curve has no Newton-certifiable rational pole. */
  int multi = (polespec != NULL);
  int nvirt = 0;
  if (use_newton)
    for (int k = 0; k < ne; ++k)
      if (!ed[k].inf && ed[k].len == 1)
        nvirt++;
  if (!multi && nrat == 0 && nvirt == 0) {
    fprintf(stderr,
            "ERROR: no rational one-point place is available for the Newton "
            "basis on this curve.\n"
            "       Multi-point construction is opt-in; specify --pole "
            "explicitly to use it.\n");
    return 1;
  }
  if (encode_requested && !multi && !concat_inner) {
    fprintf(stderr, "ERROR: -encode is available only for the multi-point "
                    "-basic path (use --pole).\n");
    return 1;
  }

  if (multi) {
    /* ---------- multi-point: basic Niederreiter syndrome decoding ----------
     */
    /* For concatenation, nerr counts binary channel errors, not AG-symbol
       errors.  The AG divisor only needs to support the outer parity-check
       construction; do not size its Riemann--Roch split using bit errors. */
    int ag_t = concat_inner ? 0 : nerr;
    int degG = 0;
    if (polespec) {
      char *s = strdup(polespec), *tok = strtok(s, ",");
      while (tok) {
        int id, mm;
        if (sscanf(tok, "E%d:%d", &id, &mm) == 2) {
          if (mm <= 0) {
            fprintf(stderr, "ERROR: pole multiplicities must be positive\n");
            free(s);
            return 1;
          }
          if (id < 0 || id >= ne || !ed[id].inf) {
            fprintf(stderr, "E%d is not an edge at infinity\n", id);
            free(s);
            return 1;
          }
          for (int q = 0; q < ng_orb; ++q)
            if (pl[q].edge == id)
              mult[q] = mm;
        } else if (sscanf(tok, "#%d:%d", &id, &mm) == 2) {
          if (mm <= 0) {
            fprintf(stderr, "ERROR: pole multiplicities must be positive\n");
            free(s);
            return 1;
          }
          if (id < 0 || id >= ng_orb) {
            fprintf(stderr, "no place #%d\n", id);
            free(s);
            return 1;
          }
          mult[id] = mm;
        } else {
          fprintf(stderr, "--pole: \"E1:10,#3:5\"\n");
          return 1;
        }
        tok = strtok(NULL, ",");
      }
      free(s);
    } else {
      int best = 0;
      for (int q = 1; q < ng_orb; ++q)
        if (pl[q].d < pl[best].d)
          best = q;
      int need = 2 * ag_t + 3 * GENUS + 4;
      mult[best] = (need + pl[best].d - 1) / pl[best].d;
    }
    for (int q = 0; q < ng_orb; ++q)
      degG += mult[q] * pl[q].d;
    int need = 2 * ag_t + 3 * GENUS - 1;
    printf("# multi-point divisor G: deg G = %d, designed distance >= %d; "
           "%s needs deg G >= %d for t = %d\n",
           degG, degG - 2 * GENUS + 2,
           concat_inner ? "AG/Goppa outer construction" : "basic algorithm",
           need, ag_t);
    if (degG < need) {
      fprintf(
          stderr,
          "deg G too small for %d AG errors with the basic algorithm (tau = %d)\n",
          ag_t, (degG - 3 * GENUS + 1) / 2);
      return 1;
    }
    /* split G = F + H with deg F >= t+g and deg H >= t+2g-1 */
    int lo = ag_t + GENUS, hi = degG - ag_t - 2 * GENUS + 1,
        target = (lo + hi) / 2;
    int fm[MAXORB], hm[MAXORB], degF = 0;
    for (int q = 0; q < ng_orb; ++q) {
      fm[q] = (int)(((long long)mult[q] * target + degG / 2) / degG);
      if (fm[q] > mult[q])
        fm[q] = mult[q];
      degF += fm[q] * pl[q].d;
    }
    for (int tries = 0; tries < 400 && (degF < lo || degF > hi); ++tries) {
      int q = tries % ng_orb;
      if (degF < lo && fm[q] < mult[q]) {
        fm[q]++;
        degF += pl[q].d;
      } else if (degF > hi && fm[q] > 0) {
        fm[q]--;
        degF -= pl[q].d;
      }
    }
    for (int q = 0; q < ng_orb; ++q)
      hm[q] = mult[q] - fm[q];
    if (degF < lo || degF > hi) {
      fprintf(stderr,
              "could not split G into F + H (deg F = %d, need %d..%d)\n", degF,
              lo, hi);
      return 1;
    }
    printf("# splitting G = F + H: deg F = %d (>= t+g = %d), deg H = %d (>= "
           "t+2g-1 = %d)\n",
           degF, lo, degG - degF, ag_t + 2 * GENUS - 1);
    LSp LG = {0}, LF = {0}, LH = {0};
    double t0 = now();
    if (solve_L(mult, -1, Jstep, Jmax, 0, &LG) != 0) {
      fprintf(stderr, "ERROR: could not certify a complete basis of L(G); "
                      "multi-point construction aborted\n");
      free(LG.B);
      return 1;
    }
    if (solve_L(fm, -1, Jstep, Jmax, 0, &LF) != 0) {
      fprintf(stderr, "ERROR: could not certify a complete basis of L(F); "
                      "multi-point construction aborted\n");
      free(LG.B);
      free(LF.B);
      return 1;
    }
    if (solve_L(hm, -1, Jstep, Jmax, 0, &LH) != 0) {
      fprintf(stderr, "ERROR: could not certify a complete basis of L(H); "
                      "multi-point construction aborted\n");
      free(LG.B);
      free(LF.B);
      free(LH.B);
      return 1;
    }
    printf("# dim L(G) = %d, dim L(F) = %d, dim L(H) = %d  (%.2f s)\n", LG.dim,
           LF.dim, LH.dim, now() - t0);
    fe *VG = eval_space(&LG), *VF = eval_space(&LF), *VH = eval_space(&LH);
    int rG = LG.dim, dF = LF.dim, dH = LH.dim;
    if (!VG || !VF || !VH) {
      fprintf(stderr,
              "ERROR: could not allocate one or more evaluation matrices\n");
      free(LG.B);
      free(LF.B);
      free(LH.B);
      free(VG);
      free(VF);
      free(VH);
      return 1;
    }
    if (rG <= 0 || rG >= npts || dF <= 0 || dH <= 0) {
      fprintf(
          stderr,
          "ERROR: invalid evaluation dimensions (rG=%d, dF=%d, dH=%d, n=%d)\n",
          rG, dF, dH, npts);
      free(LG.B);
      free(LF.B);
      free(LH.B);
      free(VG);
      free(VF);
      free(VH);
      return 1;
    }
    fe *R0 = malloc(sizeof(fe) * (size_t)npts * rG);
    if (!R0 || right_inverse(VG, rG, npts, R0) != 0) {
      fprintf(stderr, "ERROR: evaluation map L(G) -> GF(q)^n is not injective "
                      "or right-inverse construction failed\n");
      free(LG.B);
      free(LF.B);
      free(LH.B);
      free(VG);
      free(VF);
      free(VH);
      free(R0);
      return 1;
    }
    /* Niederreiter mode needs only the parity-check matrix VG and its
       right inverse.  Construct a generator matrix only when -encode was
       explicitly requested. */
    fe *Gen = NULL;
    int kcode = npts - rG;
    if (encode_requested) {
      fe *VGwork = malloc(sizeof(fe) * (size_t)rG * npts);
      if (!VGwork) {
        fprintf(stderr,
                "ERROR: out of memory while constructing generator matrix\n");
        free(R0);
        free(LG.B);
        free(LF.B);
        free(LH.B);
        free(VG);
        free(VF);
        free(VH);
        return 1;
      }
      memcpy(VGwork, VG, sizeof(fe) * (size_t)rG * npts);
      kcode = nullspace_gf(VGwork, rG, npts, &Gen);
      free(VGwork);
      if (kcode != npts - rG || !Gen) {
        fprintf(stderr,
                "ERROR: generator dimension mismatch (got %d, expected %d)\n",
                kcode, npts - rG);
        free(Gen);
        free(R0);
        free(LG.B);
        free(LF.B);
        free(LH.B);
        free(VG);
        free(VF);
        free(VH);
        return 1;
      }
      for (int s = 0; s < rG; ++s)
        for (int j = 0; j < kcode; ++j)
          if (dotv(VG + (size_t)s * npts, Gen + (size_t)j * npts, npts)) {
            fprintf(
                stderr,
                "ERROR: generated basis fails parity-check row %d, vector %d\n",
                s, j);
            free(Gen);
            free(R0);
            free(LG.B);
            free(LF.B);
            free(LH.B);
            free(VG);
            free(VF);
            free(VH);
            return 1;
          }
    }
    if (concat_inner) {
      int rc = ag_concat_run(Gen, VG, npts, kcode, rG, m_bits,
                             concat_inner, poly, trials, nerr, verbose,
                             show_values, concat_burst, concat_interleave);
      free(Gen);
      free(R0);
      free(LG.B);
      free(LF.B);
      free(LH.B);
      free(VG);
      free(VF);
      free(VH);
      return rc;
    }
    printf("# code: n = %d, check rows = %d, k = %d, designed distance >= %d, "
           "setup %.2f s\n",
           npts, rG, npts - rG, degG - 2 * GENUS + 2, now() - ts);
    fe *err = calloc((size_t)npts, sizeof(fe));
    fe *codeword = encode_requested ? calloc((size_t)npts, sizeof(fe)) : NULL;
    fe *r = calloc((size_t)npts, sizeof(fe));
    fe *sigma = malloc(sizeof(fe) * rG), *estar = malloc(sizeof(fe) * npts);
    fe *Mt = malloc(sizeof(fe) * (size_t)dH * dF),
       *AE = malloc(sizeof(fe) * (size_t)dF * npts);
    if (!err || (encode_requested && !codeword) || !r || !sigma || !estar ||
        !Mt || !AE) {
      fprintf(stderr,
              "ERROR: out of memory in multi-point Niederreiter decoder\n");
      free(err);
      free(codeword);
      free(r);
      free(sigma);
      free(estar);
      free(Mt);
      free(AE);
      free(Gen);
      free(R0);
      free(LG.B);
      free(LF.B);
      free(LH.B);
      free(VG);
      free(VF);
      free(VH);
      return 1;
    }
    int okc = 0;
    for (int tr = 0; tr < trials; ++tr) {
      make_errors(err);
      if (encode_requested) {
        memset(codeword, 0, sizeof(fe) * (size_t)npts);
        /* Optional conventional AG-codeword round-trip test. */
        for (int j = 0; j < kcode; ++j) {
          fe msg = (fe)(rnd() % (uint64_t)FS);
          if (!msg)
            continue;
          for (int k = 0; k < npts; ++k)
            codeword[k] ^= mul(msg, Gen[(size_t)j * npts + k]);
        }
      }
      for (int k = 0; k < npts; ++k)
        r[k] = (encode_requested ? codeword[k] : 0) ^ err[k];
      double t1 = now();
      /* Syndrome of the actual received word c+e, not a pre-filled error
         syndrome.  Hc=0, hence this must equal He. */
#pragma omp parallel
      {
#pragma omp for schedule(static)
        for (int s = 0; s < rG; ++s)
          sigma[s] = dotv(VG + (size_t)s * npts, r, npts);
#pragma omp for schedule(static)
        for (int k = 0; k < npts; ++k) {
          fe acc = 0;
          for (int s = 0; s < rG; ++s)
            if (sigma[s] && R0[(size_t)k * rG + s])
              acc ^= mul(R0[(size_t)k * rG + s], sigma[s]);
          estar[k] = acc;
        }
#pragma omp for schedule(static)
        for (int a = 0; a < dF; ++a)
          for (int k = 0; k < npts; ++k)
            AE[(size_t)a * npts + k] = mul(estar[k], VF[(size_t)a * npts + k]);
#pragma omp for schedule(static)
        for (int z = 0; z < dH * dF; ++z) {
          int b = z / dF, a = z - b * dF;
          Mt[z] = dotv(AE + (size_t)a * npts, VH + (size_t)b * npts, npts);
        }
      }
      fe *ker;
      int kd = nullspace_gf(Mt, dH, dF, &ker);
      int *Z = malloc(sizeof(int) * npts), nz = 0;
      for (int k = 0; k < npts; ++k)
        Z[nz++] = k;
      for (int v = 0; v < kd && nz > 0; ++v) {
        int nn = 0;
        for (int zi = 0; zi < nz; ++zi) {
          int k = Z[zi];
          fe acc = 0;
          for (int a = 0; a < dF; ++a)
            if (ker[(size_t)v * dF + a])
              acc ^= mul(ker[(size_t)v * dF + a], VF[(size_t)a * npts + k]);
          if (!acc)
            Z[nn++] = k;
        }
        nz = nn;
      }
      free(ker);
      fe *val = malloc(sizeof(fe) * (nz + 1));
      int ok = 0, nfound = 0;
      if (nz > 0 && nz <= MAXT && solve_values(VG, rG, sigma, Z, nz, val)) {
        int same = 0;
        for (int c = 0; c < nz; ++c)
          if (val[c]) {
            nfound++;
            for (int e = 0; e < NUM_ERRORS; ++e)
              if (true_pos[e] == Z[c] && true_val[e] == val[c])
                same++;
          }
        ok = (same == NUM_ERRORS && nfound == NUM_ERRORS);
        if (ok) {
          /* Correct the received word.  In default Niederreiter mode the
             transmitted vector is zero; with -encode it is a random codeword.
           */
          for (int c = 0; c < nz; ++c)
            r[Z[c]] ^= val[c];
          if (encode_requested) {
            if (memcmp(r, codeword, sizeof(fe) * (size_t)npts) != 0)
              ok = 0;
          } else {
            for (int k = 0; k < npts; ++k)
              if (r[k] != 0) {
                ok = 0;
                break;
              }
          }
        }
      }
      okc += ok;
      if (trials == 1 || verbose)
        printf("# BASIC trial %d: %d kernel vectors of dim L(F) = %d, %d "
               "common zeros, %d errors reported -> %s (%.2f s)\n",
               tr + 1, kd, dF, nz, nfound, ok ? "OK" : "FAILED", now() - t1);
      free(Z);
      free(val);
    }
    if (encode_requested)
      printf("# BASIC: %d of %d optional encode -> corrupt -> decode -> "
             "exact-recovery trials succeeded (%d errors)\n",
             okc, trials, NUM_ERRORS);
    else
      printf("# BASIC: %d of %d Niederreiter syndrome -> decode -> "
             "zero-error-vector recovery trials succeeded (%d errors)\n",
             okc, trials, NUM_ERRORS);
    free(err);
    free(codeword);
    free(r);
    free(sigma);
    free(estar);
    free(Mt);
    free(AE);
    free(Gen);
    free(R0);
    free(LG.B);
    free(LF.B);
    free(LH.B);
    free(VG);
    free(VF);
    free(VH);
    return okc == trials ? 0 : 1;
  }

  /* ---------- one-point codes ---------- */
  int M = Mreq;
  if (concat_inner && M <= 0) {
    /* For a one-point concatenated outer code, leave 2g+1 parity checks
     * (Riemann--Roch gives k=M+1-g).  This is a conservative default with
     * useful erasure capacity; users may tune it explicitly with --M. */
    M = npts - GENUS - 2;
    if (M < 2 * GENUS - 1) {
      fprintf(stderr,
              "ERROR: one-point AG/Goppa concatenation needs more rational "
              "evaluation points: n=%d, g=%d; default M=%d is below "
              "2g-1=%d. Increase the field/point count or use --pole for "
              "the optional multi-point construction.\n",
              npts, GENUS, M, 2 * GENUS - 1);
      return 1;
    }
    printf("# one-point concatenation default: M=%d (outer parity checks "
           "about %d)\n", M, npts - (M + 1 - GENUS));
  }
  if (run_ecp && M <= 0)
    M = 2 * GENUS; /* choose the smallest one-point divisor for ECP bounds */
  if (M <= 0) {
    /* default: the smallest M for which FR/BMS (basic bound) AND voting all
     * correct nerr errors; if that does not fit, the voting minimum */
    int Mb = 2 * nerr + 3 * GENUS - 1;
    if (Mb <= limit_Mmax(GENUS, npts))
      M = Mb;
    else {
      M = 2 * nerr + 2 * GENUS;
      printf("# M = %d (voting minimum): the basic bound would need M = %d, "
             "which does not fit; FR/BMS will be run beyond their proven "
             "bound\n",
             M, Mb);
    }
  }
  if (M < 2 * GENUS)
    M = 2 * GENUS;
  M_SYND = M;
  M_EXT = M + GENUS;
  TAU = (M - 3 * GENUS + 1) / 2;
  TAU_V = (M - 2 * GENUS + 1) / 2;
  /* FR-VOTE needs locator functions through the voting radius, not the
   * smaller basic-decoding radius used by FR/BMS. */
  L_DEG = (run_frvote ? TAU_V : TAU) + GENUS;
  if (M + 1 > MAXB - 1) {
    fprintf(stderr, "M too large for this build\n");
    print_limits(GENUS, npts, nerr);
    return 1;
  }
  if (!concat_inner && nerr > TAU_V && !run_common_binary && !run_ecp) {
    fprintf(stderr,
            "%d errors exceed tau_vote = %d for M = %d (need M >= %d)\n", nerr,
            TAU_V, M, 2 * nerr + 2 * GENUS);
    print_limits(GENUS, npts, nerr);
    return 1;
  }
  if (concat_inner)
    printf("# one-point AG/Goppa outer setup: M = %d, g = %d\n", M, GENUS);
  else if (run_ecp)
    printf("# one-point code: M = %d, g = %d, ECP requested for t = %d\n",
           M, GENUS, NUM_ERRORS);
  else
    printf("# one-point code: M = %d, g = %d, basic algorithm t <= %d, with "
           "majority voting t <= %d\n",
           M, GENUS, TAU, TAU_V);
  if (npts <= M_EXT) {
    fprintf(stderr,
            "# the curve has only n = %d rational points, but M + g = %d >= n: "
            "the evaluation map is not injective.\n"
            "# Use a larger field (--m) or a smaller M/number of errors (need "
            "n > M + g).\n",
            npts, M_EXT);
    print_limits(GENUS, npts, nerr);
    return 1;
  }

  /* Candidate poles are accepted only when the Newton polygon certifies the
   * complete one-point basis.  No Laurent-series fallback is performed. */
  LSp Lc[NCAND_MAX];
  int *pv[NCAND_MAX];
  NBas NB[NCAND_MAX];
  int via_newton[NCAND_MAX];
  memset(via_newton, 0, sizeof via_newton);
  memset(NB, 0, sizeof NB);
  int best = -1, bestdfr = -1, besta = 1 << 30;
  int cand[NCAND_MAX], ncand = 0;
  if (place_req >= 0) {
    if (place_req < ng_orb && pl[place_req].d == 1)
      cand[ncand++] = place_req;
    else
      printf("# requested place #%d is not a rational place at infinity -> "
             "searching for a substitute\n",
             place_req);
  }
  if (ncand == 0) {
    for (int i = 0; i < nrat; ++i)
      cand[ncand++] = ratidx[i];
    if (use_newton)
      for (int k = 0; k < ne; ++k)
        if (!ed[k].inf && ed[k].len == 1)
          cand[ncand++] = MAXORB + k;
  }
  printf("# candidate poles: %d\n", ncand);
  for (int pass = 0; pass < 1; ++pass) {
    for (int ci = 0; ci < ncand; ++ci) {
      int id = cand[ci], virt = id >= MAXORB;
      double t0 = now();
      int dimq = 0, nconst = 1;
      const int *rhos = NULL;
      if (pass == 0) {
        int e = virt ? id - MAXORB : pl[id].edge;
        if (!use_newton || newton_basis(e, M_EXT, &NB[id]) != 0)
          continue;
        via_newton[id] = 1;
        dimq = NB[id].n;
        rhos = NB[id].rho;
      } else {
        /* Laurent-series basis construction is intentionally disabled. */
        continue;
      }
      static int inHt[SGMAX];
      memset(inHt, 0, sizeof inHt);
      for (int x = M_EXT + 1; x < SGMAX; ++x)
        inHt[x] = 1;
      for (int d = 0; d < dimq; ++d)
        inHt[rhos[d]] = 1;
      int ngaps = 0;
      for (int x = 1; x <= M_EXT; ++x)
        if (!inHt[x])
          ngaps++;
      int a1 = 1 << 30;
      for (int x = 1; x <= M_EXT; ++x)
        if (inHt[x]) {
          a1 = x;
          break;
        }
      int span = SGMAX - M - 2;
      if (span > 400)
        span = 400;
      int dfr = dfr_of(inHt, M, span);
      char lab[64];
      if (virt)
        snprintf(lab, sizeof lab, "place of edge E%d (%s)", id - MAXORB,
                 "affine side");
      else
        snprintf(lab, sizeof lab, "place #%d", id);
      printf("#   %s: first non-zero pole order %d, %d gaps (g = %d)%s, "
             "d_FR(M) = %d -> t = %d  [basis via %s]  (%.3f s)\n",
             lab, a1, ngaps, GENUS,
             ngaps == GENUS && nconst == 1 ? "" : "  [semigroup check FAILED]",
             dfr, (dfr - 1) / 2, "Newton polygon", now() - t0);
      if (ngaps == GENUS && nconst == 1 &&
          (dfr > bestdfr || (dfr == bestdfr && a1 < besta))) {
        best = id;
        bestdfr = dfr;
        besta = a1;
      }
    }
  }
  if (best < 0) {
    fprintf(stderr, "no usable rational place for a one-point code\n");
    return 1;
  }
  /* free the Newton bases of the candidates that were not chosen */
  for (int ci = 0; ci < ncand; ++ci) {
    int id = cand[ci];
    if (id != best && via_newton[id]) {
      free(NB[id].rho);
      free(NB[id].ei);
      free(NB[id].ej);
      memset(&NB[id], 0, sizeof NB[id]);
      via_newton[id] = 0;
    }
  }
  if (best >= MAXORB)
    printf("# chosen pole: place of edge E%d  (d_FR = %d, designed t = %d)\n",
           best - MAXORB, bestdfr, (bestdfr - 1) / 2);
  else
    printf("# chosen pole: place #%d  (d_FR = %d, designed t = %d)\n", best,
           bestdfr, (bestdfr - 1) / 2);
  if (via_newton[best] && use_poly_opt) {
    int e = best >= MAXORB ? best - MAXORB : pl[best].edge;
    if (poly_setup(e, M_EXT + (use_kotter_values ? 128 : 0), &NB[best]) == 0) {
      use_poly = 1;
      if (use_kotter_values) {
        kotter_ready = kotter_prepare_module(e, &NB[best]);
        if (!kotter_ready) {
          fprintf(stderr, "Kötter path unsupported for this curve/place or "
                          "module setup failed.\n");
          return 2;
        }
      }
      bas_ei = NB[best].ei; /* poly_setup の後 */
      bas_ej = NB[best].ej;
      if (hermitian_q > 0)
        printf("# Hermitian monomial syndrome tables on %d polynomial "
               "semigroup points (no dual basis)\n",
               NPT);
      else
        printf("# monomial syndrome tables on %d lattice points of the Newton "
               "cone (no dual basis, f kept as monomial coefficients)\n",
               NPT);
    }
  }

  if (via_newton[best]) {
    /* ---- basis straight from the lattice points: monomials x^i y^j ---- */
    NBas *B = &NB[best];
    if (B->laurent) {
      /* Laurent monomials: x = 0 or y = 0 are not allowed as positions.
       * The only axis point of the model is the origin (checked), which is a
       * point where the chosen functions are regular but whose x^-1 / y^-1
       * values are undefined, so it is dropped from the positions. */
      int original_npts = npts;
      int k = 0;
      if (api_mode && api_received_len == (size_t)original_npts) {
        api_owned_received = malloc(sizeof(fe) * (size_t)original_npts);
        api_position_map = malloc(sizeof(size_t) * (size_t)original_npts);
        if (!api_owned_received || !api_position_map) {
          fprintf(stderr,
                  "API allocation failed while puncturing Laurent positions\n");
          return 1;
        }
      }
      for (int p = 0; p < original_npts; ++p)
        if (PXv[p] && PYv[p]) {
          PXv[k] = PXv[p];
          PYv[k] = PYv[p];
          if (api_mode && api_received_len == (size_t)original_npts) {
            api_owned_received[k] = api_received[p];
            api_position_map[k] = (size_t)p;
          }
          k++;
        }
      if (api_mode && api_received_len == (size_t)original_npts) {
        api_received = api_owned_received;
        api_received_len = (size_t)k;
      }
      printf("# Laurent monomials used: %d axis points dropped, n = %d\n",
             npts - k, k);
      npts = k;
      if (npts <= M_EXT) {
        fprintf(stderr, "# only n = %d positions left, but M + g = %d >= n.\n",
                npts, M_EXT);
        return 1;
      }
    }
    nbasis = B->n;
    if (run_frvote) {
      /* FR-VOTE's Newton syndrome path only needs monomial exponents and
       * point coordinates.  Do not allocate nbasis*npts field elements;
       * selected error positions are evaluated on demand. */
      for (int i = 0; i < B->n; ++i) {
        basis[i].rho = B->rho[i];
        basis[i].val = NULL;
      }
    } else {
      int *lxs = malloc(sizeof(int) * (size_t)npts);
      int *lys = malloc(sizeof(int) * (size_t)npts);
      for (int k = 0; k < npts; ++k) {
        lxs[k] = PXv[k] ? (int)glog[PXv[k]] : -1;
        lys[k] = PYv[k] ? (int)glog[PYv[k]] : -1;
      }
      for (int i = 0; i < B->n; ++i) {
        basis[i].rho = B->rho[i];
        basis[i].val = malloc(sizeof(fe) * npts);
        mono_eval_row(basis[i].val, npts, lxs, lys, B->ei[i], B->ej[i]);
      }
      free(lxs);
      free(lys);
    }
  } else {
    LSp *L = &Lc[best];
    fe *V = eval_space(L);
    { /* basis sorted by pole order */
      int dim = L->dim;
      int *ord = malloc(sizeof(int) * dim);
      for (int i = 0; i < dim; ++i)
        ord[i] = i;
      for (int i = 1; i < dim; ++i) {
        int x = ord[i], j = i - 1;
        while (j >= 0 && pv[best][ord[j]] > pv[best][x]) {
          ord[j + 1] = ord[j];
          j--;
        }
        ord[j + 1] = x;
      }
      nbasis = dim;
      for (int i = 0; i < dim; ++i) {
        basis[i].rho = pv[best][ord[i]];
        /* V is already a contiguous dim x npts evaluation matrix.  Point
         * directly into its rows instead of allocating/copying every row.
         * This also fixes the former leak of V and halves this phase's peak
         * evaluation-storage requirement. */
        basis[i].val = V + (size_t)ord[i] * (size_t)npts;
      }
      free(ord);
    }
  }
  if (concat_inner) {
    /* One-point AG/Goppa concatenation: use the Newton basis functions as
     * outer parity checks and let FR-VOTE decode their direct syndromes.
     * No outer generator matrix or nullspace is constructed in this path. */
    int kcode = 0;
    for (int i = 0; i < nbasis; ++i)
      if (basis[i].rho <= M)
        kcode++;
    int expected_k = M + 1 - GENUS;
    if (kcode != expected_k || kcode <= 0 || kcode >= npts) {
      fprintf(stderr,
              "ERROR: one-point AG basis dimension mismatch: got k=%d, "
              "Riemann-Roch expects %d (M=%d, g=%d).\n"
              "       This curve/place does not provide a certified "
              "one-point Newton basis for the requested divisor.\n",
              kcode, expected_k, M, GENUS);
      return 2;
    }
    printf("# one-point AG outer code: [n=%d,k=%d] over GF(2^%d), "
           "Newton check functions; designed distance >= %d\n",
           npts, npts - kcode, m_bits, M - 2 * GENUS + 2);
    int rc = ag_concat_frvote_run(npts, kcode, m_bits, concat_inner, poly,
                                  trials, nerr, verbose, show_values,
                                  concat_burst, concat_interleave);
    for (int i = 0; i < nbasis; ++i)
      free(basis[i].val);
    if (via_newton[best]) {
      free(NB[best].rho);
      free(NB[best].ei);
      free(NB[best].ej);
    }
    return rc;
  }

  init_semigroup();
  ASHIFT = basis[1].rho;
  ROW_MAX = M_SYND - ASHIFT;
  nz_idx = malloc(sizeof(int) * npts);
  printf("# setup (curve, places, basis up to M+g, evaluation): %.2f s;  "
         "semigroup has %d generators, smallest pole %d\n",
         now() - ts, ngenH, ASHIFT);
  if (run_vote && !use_poly) {
    double t0 = now();
    if (build_dual() != 0 || dual_selftest() != 0) {
      fprintf(stderr, "dual basis construction failed\n");
      return 1;
    }
    init_dual_logs();
    printf("# setup (dual basis for the top g coefficients, once per code): "
           "%.2f s\n",
           now() - t0);
  }
  fe *r = malloc(sizeof(fe) * npts);
  int cnt_fr = 0, cnt_bms = 0, cnt_vote = 0, cnt_frvote = 0, cnt_ecp = 0;
  int quiet = trials > 1;
  for (int tr = 0; tr < trials; ++tr) {
    if (api_mode) {
      if (!run_frvote || run_fr || run_bms || run_vote) {
        fprintf(stderr, "API mode requires -frvote alone (do not enable "
                        "-fr/-bms/-vote).\n");
        return 1;
      }
      if (api_received_len != (size_t)npts) {
        fprintf(stderr, "API received length mismatch: got %zu, expected %d\n",
                api_received_len, npts);
        return 1;
      }
      for (int k = 0; k < npts; ++k) {
        if (api_received[k] >= FS) {
          fprintf(stderr, "API received[%d]=%u is outside GF(2^%d)\n", k,
                  (unsigned)api_received[k], FM);
          return 1;
        }
      }
      memcpy(r, api_received, sizeof(fe) * (size_t)npts);
      set_received(r);
      if (decode_fr_vote(&res_frvote, verbose, 0) != 0 || res_frvote.failed ||
          !res_frvote.resolved) {
        fprintf(stderr,
                "API FR-VOTE decoding failed or error values unresolved.\n");
        return 2;
      }
      if ((size_t)res_frvote.nfound > api_output_capacity) {
        fprintf(stderr, "API output buffer too small: need %d, capacity %zu\n",
                res_frvote.nfound, api_output_capacity);
        return 3;
      }
      for (int k = 0; k < res_frvote.nfound; ++k) {
        api_error_positions[k] = api_position_map
                                     ? api_position_map[res_frvote.pos[k]]
                                     : (size_t)res_frvote.pos[k];
        api_error_values[k] = res_frvote.val[k];
      }
      *api_error_count = (size_t)res_frvote.nfound;
      api_completed = 1;
      break;
    }
    make_errors(r);
    set_received(r);
    if (bms2_check_enabled)
      (void)bms2_syndrome_check();
    if (run_common_binary) {
      if (!fixed_error_value) {
        fprintf(stderr, "# COMMON-BINARY requires --fixed-error-value v (known "
                        "common value)\n");
        return 1;
      }
      decode_common_value_positions(fixed_error_value, NUM_ERRORS, verbose);
    }
    if (!quiet) {
      printf("# %d errors injected (seed %u)%s\n", NUM_ERRORS, seed,
             NUM_ERRORS > TAU ? "  -- beyond the basic bound" : "");
      if (fixed_error_value)
        printf("# injected error values: all fixed to %u\n",
               (unsigned)fixed_error_value);
      else
        puts("# injected error values: random nonzero field elements");
    }
    if (force && NUM_ERRORS > TAU && !quiet && (run_fr || run_bms))
      puts("# FR/BMS: t exceeds the proven basic bound; running anyway "
           "(expected to "
           "fail)");
    int okf = 1, okb = 1, okv = 1, okfv = 1;
    int okecp = 1;
    if (run_ecp) {
      int ep[MAXT], en=0; fe ev[MAXT];
      double te=now();
      okecp=ecp_ag_decode(npts,M,GENUS,NUM_ERRORS,recv,ep,ev,&en);
      int same=0;
      if (okecp) for (int q=0;q<en;q++) for (int e=0;e<NUM_ERRORS;e++)
        if (ep[q]==true_pos[e] && ev[q]==true_val[e]) same++;
      okecp = okecp && en==NUM_ERRORS && same==NUM_ERRORS;
      cnt_ecp += okecp;
      if (!quiet) printf("# ECP  %d of %d errors recovered with correct values -> %s (%.3f s)\n",
                         same,NUM_ERRORS,okecp?"OK":"FAILED",now()-te);
    }
    if (run_frvote) {
      decode_fr_vote(&res_frvote, verbose, quiet);
      okfv = verify("FR-VOTE", &res_frvote, quiet);
      cnt_frvote += okfv;
    }
    if (run_fr) {
      if (NUM_ERRORS > TAU && !force) {
        okf = 0;
        if (!quiet)
          puts("# FR   skipped: beyond the basic bound");
      } else {
        decode_fr(&res_fr, quiet);
        okf = verify("FR", &res_fr, quiet);
        cnt_fr += okf;
      }
    }
    if (run_bms) {
      if (NUM_ERRORS > TAU && !force) {
        okb = 0;
        if (!quiet)
          puts("# BMS  skipped: beyond the basic bound");
      } else {
        decode_bms(&res_bms, 0, verbose, quiet);
        okb = verify("BMS", &res_bms, quiet);
        cnt_bms += okb;
      }
    }
    if (run_vote) {
      decode_bms(&res_vote, 1, verbose, quiet);
      okv = verify("VOTE", &res_vote, quiet);
      cnt_vote += okv;
    }
    if (show_values && !quiet) {
      if (run_fr && (NUM_ERRORS <= TAU || force))
        show_error_values(
            use_kotter_values
                ? "FR/Kotter"
                : (use_forney_values ? "FR/Forney" : "FR/solve_L"),
            &res_fr);
      if (run_bms && (NUM_ERRORS <= TAU || force))
        show_error_values("BMS", &res_bms);
      if (run_vote)
        show_error_values("VOTE", &res_vote);
      if (run_frvote)
        show_error_values("FR-VOTE", &res_frvote);
    }
    if (!quiet) {
      if (run_fr && (NUM_ERRORS <= TAU || force))
        printf("# time FR   : reduction %.3f s, position extraction %.3f s, "
               "values "
               "%.3f s\n",
               res_fr.t_main, res_fr.t_extract, res_fr.t_resolve);
      if (run_bms && (NUM_ERRORS <= TAU || force))
        printf("# time BMS  : loop %.3f s, scan %.3f s, values %.3f s\n",
               res_bms.t_main, res_bms.t_extract, res_bms.t_resolve);
      if (run_vote)
        printf("# time VOTE : BMS+voting %.3f s, scan %.3f s, values %.3f s\n",
               res_vote.t_main, res_vote.t_extract, res_vote.t_resolve);
      if (run_frvote)
        printf("# time FR-VOTE: FR-kernel+voting %.3f s, scan %.3f s, values "
               "%.3f s\n",
               res_frvote.t_main, res_frvote.t_extract, res_frvote.t_resolve);
      if (run_fr && run_vote && (NUM_ERRORS <= TAU || force) && !res_fr.failed)
        printf("# FR / VOTE time ratio: %.1f x\n",
               total_time(&res_fr) / total_time(&res_vote));
    } else {
      int bo = NUM_ERRORS <= TAU || force;
      printf("# trial %d: FR %s  BMS %s  VOTE %s\n", tr + 1,
             (run_fr && bo) ? (okf ? "ok" : "FAIL") : "n/a",
             (run_bms && bo) ? (okb ? "ok" : "FAIL") : "n/a",
             run_vote ? (okv ? "ok" : "FAIL") : "n/a");
      if (run_frvote)
        printf("# trial %d: FR-VOTE %s\n", tr + 1, okfv ? "ok" : "FAIL");
      if (run_ecp)
        printf("# trial %d: ECP %s\n", tr + 1, okecp ? "ok" : "FAIL");
    }
  }
  if (trials > 1) {
    char sf[32], sb[32], sv[32];
    int basic_ok = NUM_ERRORS <= TAU || force;
    if (run_fr && basic_ok)
      snprintf(sf, sizeof sf, "%d/%d", cnt_fr, trials);
    else
      snprintf(sf, sizeof sf, "n/a");
    if (run_bms && basic_ok)
      snprintf(sb, sizeof sb, "%d/%d", cnt_bms, trials);
    else
      snprintf(sb, sizeof sb, "n/a");
    if (run_vote)
      snprintf(sv, sizeof sv, "%d/%d", cnt_vote, trials);
    else
      snprintf(sv, sizeof sv, "n/a");
    printf("# summary of %d trials with %d errors: FR %s, BMS %s, VOTE %s\n",
           trials, NUM_ERRORS, sf, sb, sv);
    if (run_frvote)
      printf("# FR-VOTE: %d/%d\n", cnt_frvote, trials);
    if (run_ecp)
      printf("# ECP: %d/%d\n", cnt_ecp, trials);
  }
  if (via_newton[best])
    free(NB[best].rho); /* ei/ej stay referenced by bas_ei/bas_ej */
  return 0;
}
int noize_decode(int argc, char **argv, const uint16_t *received,
                 size_t received_len, size_t *error_positions,
                 uint16_t *error_values, size_t output_capacity,
                 size_t *error_count) {
  if (!argv || argc < 1 || !received || !error_count ||
      (output_capacity && (!error_positions || !error_values)))
    return -1;
  *error_count = 0;
  api_received = (const fe *)received;
  api_received_len = received_len;
  api_error_positions = error_positions;
  api_error_values = (fe *)error_values;
  api_output_capacity = output_capacity;
  api_error_count = error_count;
  api_mode = 1;
  api_completed = 0;
  api_owned_received = NULL;
  api_position_map = NULL;
  int rc = noize_run(argc, argv);
  free(api_owned_received);
  free(api_position_map);
  api_owned_received = NULL;
  api_position_map = NULL;
  api_mode = 0;
  api_received = NULL;
  api_error_positions = NULL;
  api_error_values = NULL;
  api_error_count = NULL;
  if (rc != 0)
    return rc;
  return api_completed ? 0 : -10;
}

#ifndef noize_NO_MAIN
/* Unified dispatcher is defined at the end of this file. */
#endif
/* ================= classical-code module (namespaced from et.c) ============
 */
#define CODES_NO_MAIN
#define fe et_fe
#define Poly et_Poly
#define rng_s et_rng_s
#define now et_now
#define peval et_peval
#define ptrim et_ptrim
#define rnd et_rnd
#define usage et_usage
#define gpow et_gpow
/*
 * codes.c -- encoders and decoders for classical algebraic code families
 *            (same conventions as noizy2.c: GF(2^m) with m <= 12, xorshift
 *            RNG with a seed, "# ..." report lines, --errors / -n / -v)
 *
 *   bch      narrow-sense / any-b binary BCH code (primitive, cyclic n | 2^m-1,
 *            or shortened), decoded by Berlekamp-Massey + Chien search
 *   golay    binary Golay [23,12,7] and extended [24,12,8], syndrome table
 *            (perfect code: every error pattern of weight <= 3 is corrected)
 *   goppa    binary classical Goppa code Gamma(L,g), g irreducible, corrects
 *            t = deg g errors: extended Euclid on the key equation modulo g^2
 *            (Patterson-equivalent, no square roots needed)
 *   alt      random binary alternant code (alpha_i, y_i random), t = floor(r/2)
 *   srivastava  original Srivastava alternant (multiplicity 1, z_i=alpha_i^mu)
 *   gs       generalized Srivastava alternant, with s poles and multiplicity t
 *            (designed distance >= s*t+1); BM by default, generic EEA via -eea
 *   ldpc     regular LDPC code (n, dv, dc, seed), systematic encoding, soft
 *            sum-product / normalized min-sum belief-propagation decoding
 *            (also usable as the inner code of product/concat)
 *   rm       binary Reed-Muller code RM(r,m) = [2^m, sum_{i<=r} C(m,i),
 *            2^(m-r)], Reed majority-logic decoding (corrects 2^(m-r-1)-1)
 *   hadamard Hadamard code [2^m, m, 2^(m-1)]; with -fht also RM(1,m): maximum-
 *            likelihood decoding by the fast Walsh-Hadamard transform,
 *            O(n log n), beyond the bounded-distance radius
 *   grs      generalized Reed-Solomon code over GF(2^m); errors AND erasures
 *            (2e + s <= n-k), Berlekamp-Massey + Forney
 *   rank     Gabidulin code (rank metric, linearized polynomials), n <= m,
 *            corrects t = (n-k)/2 rank errors, linearized key equation
 *   product  product code C1 x C2 of two binary codes (bch/golay/goppa/alt
 *            specs), iterative row/column hard-decision decoding
 *   concat   concatenated code: outer GRS over GF(2^m), inner binary code
 *            (any spec with k_in >= m, shortened if k_in > m); inner decoding
 *            failures become erasures for the outer decoder
 *
 * Binary alternant-type codes (bch, alt, goppa) are encoded through the exact
 * binary kernel of the GF(2) parity-check matrix (systematic on the free
 * columns), so k is always the true dimension, not just the lower bound.
 *
 * Build: cc -O2 -std=gnu11 -o codes codes.c
 * Run:   ./codes --help
 */
#define _POSIX_C_SOURCE 200809L
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint16_t fe;

#define MAXM 12 /* GF(2^m), m <= 12 (fe is 16 bit)                      */
#define ET_PMAXD                                                               \
  1100 /* maximum polynomial degree (EXTENSION_POINT)           */

static double now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}
static void die(const char *f, ...) {
  va_list ap;
  va_start(ap, f);
  fprintf(stderr, "error: ");
  vfprintf(stderr, f, ap);
  fputc('\n', stderr);
  va_end(ap);
  exit(1);
}

/* ================= random numbers (same generator as noizy2.c) ============ */
static uint64_t rng_s = 88172645463325252ULL;
static uint64_t rnd(void) {
  rng_s ^= rng_s << 13;
  rng_s ^= rng_s >> 7;
  rng_s ^= rng_s << 17;
  return rng_s;
}
/* choose cnt distinct positions in [0,n) that are not flagged in excl */
static void pick(int n, int cnt, int *dst, const uint8_t *excl) {
  int avail = n;
  if (excl)
    for (int i = 0; i < n; ++i)
      avail -= excl[i] != 0;
  if (cnt > avail)
    die("cannot place %d distinct positions among %d available", cnt, avail);
  for (int e = 0; e < cnt; ++e) {
    int p, dup;
    do {
      p = (int)(rnd() % (uint64_t)n);
      dup = excl && excl[p];
      for (int k = 0; k < e; ++k)
        if (dst[k] == p)
          dup = 1;
    } while (dup);
    dst[e] = p;
  }
}

/* ================= GF(2^m), supported m <= 12 ================= */
typedef struct {
  int m, size, ord;
  unsigned poly;
  fe *exp, *log;
} GF;

static const unsigned DEFPOLY[13] = {
    0, 0, 0x7, 0xB, 0x13, 0x25, 0x43, 0x89, 0x11D, 0x211, 0x409, 0x805, 0x1053};

static GF *gf_new(int m, unsigned poly) {
  if (m < 2 || m > MAXM)
    die("--m must be in 2..%d", MAXM);
  if (!poly)
    poly = DEFPOLY[m];
  GF *F = calloc(1, sizeof *F);
  F->m = m;
  F->size = 1 << m;
  F->ord = F->size - 1;
  F->poly = poly;
  F->exp = malloc(sizeof(fe) * (2 * F->ord + 4));
  F->log = malloc(sizeof(fe) * F->size);
  unsigned x = 1;
  for (int i = 0; i < F->ord; ++i) {
    if (i > 0 && x == 1)
      die("polynomial 0x%X is not primitive of degree %d", poly, m);
    F->exp[i] = (fe)x;
    F->log[x] = (fe)i;
    x <<= 1;
    if (x & (unsigned)F->size)
      x ^= poly;
  }
  if (x != 1)
    die("polynomial 0x%X is not primitive of degree %d", poly, m);
  for (int i = F->ord; i < 2 * F->ord + 4; ++i)
    F->exp[i] = F->exp[i - F->ord];
  F->log[0] = 0;
  return F;
}
static inline fe gmul(const GF *F, fe a, fe b) {
  return (a && b) ? F->exp[F->log[a] + F->log[b]] : 0;
}
static inline fe ginv(const GF *F, fe a) { return F->exp[F->ord - F->log[a]]; }
static inline fe gdiv(const GF *F, fe a, fe b) {
  return gmul(F, a, ginv(F, b));
}
static fe gpow(const GF *F, fe a, long long n) {
  if (n == 0)
    return 1;
  if (!a)
    return 0;
  long long e = n % F->ord;
  if (e < 0)
    e += F->ord;
  return F->exp[(int)(((long long)F->log[a] * e) % F->ord)];
}
static fe gfrob(const GF *F, fe a, int i) { /* a^(2^i) = a^[i] */
  for (; i > 0; --i)
    a = gmul(F, a, a);
  return a;
}
static fe grand_nz(const GF *F) { return (fe)(1 + rnd() % (uint64_t)F->ord); }

/* ================= polynomials over GF(2^m); deg = -1 is the zero poly =====
 */
typedef struct {
  int deg;
  fe c[ET_PMAXD + 1];
} Poly;
static void ptrim(Poly *p) {
  while (p->deg >= 0 && p->c[p->deg] == 0)
    p->deg--;
}
static fe peval(const GF *F, const Poly *p, fe x) {
  fe acc = 0;
  for (int i = p->deg; i >= 0; --i)
    acc = (fe)(gmul(F, acc, x) ^ p->c[i]);
  return acc;
}
static void pcopy(Poly *d, const Poly *s) {
  d->deg = s->deg;
  if (s->deg >= 0)
    memcpy(d->c, s->c, sizeof(fe) * (size_t)(s->deg + 1));
}
static void pmul(const GF *F, const Poly *a, const Poly *b, Poly *r) {
  if (a->deg < 0 || b->deg < 0) {
    r->deg = -1;
    return;
  }
  if (a->deg + b->deg > ET_PMAXD)
    die("polynomial degree exceeds ET_PMAXD=%d", ET_PMAXD);
  r->deg = a->deg + b->deg;
  memset(r->c, 0, sizeof(fe) * (size_t)(r->deg + 1));
  for (int i = 0; i <= a->deg; ++i)
    if (a->c[i])
      for (int j = 0; j <= b->deg; ++j)
        r->c[i + j] ^= gmul(F, a->c[i], b->c[j]);
  ptrim(r);
}
static void padd(const Poly *a, const Poly *b, Poly *r) {
  int d = a->deg > b->deg ? a->deg : b->deg;
  fe tmp[ET_PMAXD + 1];
  for (int i = 0; i <= d; ++i)
    tmp[i] = (fe)((i <= a->deg ? a->c[i] : 0) ^ (i <= b->deg ? b->c[i] : 0));
  r->deg = d;
  if (d >= 0)
    memcpy(r->c, tmp, sizeof(fe) * (size_t)(d + 1));
  ptrim(r);
}
/* r = a mod b, q = a div b (q may be NULL); b != 0; r may alias a */
static void pdivmod(const GF *F, const Poly *a, const Poly *b, Poly *q,
                    Poly *r) {
  if (r != a)
    pcopy(r, a);
  if (q) {
    q->deg = a->deg >= b->deg ? a->deg - b->deg : -1;
    if (q->deg >= 0)
      memset(q->c, 0, sizeof(fe) * (size_t)(q->deg + 1));
  }
  fe ib = ginv(F, b->c[b->deg]);
  while (r->deg >= b->deg) {
    fe f = gmul(F, r->c[r->deg], ib);
    int sh = r->deg - b->deg;
    if (q)
      q->c[sh] = f;
    for (int i = 0; i <= b->deg; ++i)
      r->c[i + sh] ^= gmul(F, f, b->c[i]);
    ptrim(r);
  }
}
static void pderiv(const Poly *a, Poly *r) { /* characteristic 2 */
  r->deg = a->deg - 1;
  for (int i = 0; i <= r->deg; ++i)
    r->c[i] = (i % 2 == 0) ? a->c[i + 1] : 0;
  ptrim(r);
}
static void pgcd(const GF *F, const Poly *a, const Poly *b, Poly *g) {
  Poly x, y, r;
  pcopy(&x, a);
  pcopy(&y, b);
  while (y.deg >= 0) {
    pdivmod(F, &x, &y, NULL, &r);
    pcopy(&x, &y);
    pcopy(&y, &r);
  }
  pcopy(g, &x);
}
static void psqr_mod(const GF *F, const Poly *p, const Poly *g, Poly *out) {
  Poly t;
  if (p->deg < 0) {
    out->deg = -1;
    return;
  }
  t.deg = 2 * p->deg;
  memset(t.c, 0, sizeof(fe) * (size_t)(t.deg + 1));
  for (int i = 0; i <= p->deg; ++i)
    t.c[2 * i] = gmul(F, p->c[i], p->c[i]);
  pdivmod(F, &t, g, NULL, out);
}
/* g (monic, deg >= 1) is irreducible over GF(2^m) iff it has no irreducible
 * factor of degree <= deg/2, i.e. gcd(x^(Q^i) - x, g) = 1 for i <= deg/2 */
static int pirreducible(const GF *F, const Poly *g) {
  Poly P, x, d, t;
  x.deg = 1;
  x.c[0] = 0;
  x.c[1] = 1;
  if (g->deg == 1)
    return 1;
  pdivmod(F, &x, g, NULL, &P);
  for (int i = 1; i <= g->deg / 2; ++i) {
    for (int s = 0; s < F->m; ++s) { /* P <- P^(2^m) = P^Q */
      psqr_mod(F, &P, g, &t);
      pcopy(&P, &t);
    }
    padd(&P, &x, &d);
    if (d.deg < 0)
      return 0;
    pgcd(F, g, &d, &t);
    if (t.deg > 0)
      return 0;
  }
  return 1;
}

/* ================= linear algebra helpers ================= */
/* solve A x = b over GF(2^m), A is N x N row-major, destroyed; x -> b */
static int gf_solve(const GF *F, fe *A, fe *b, int N) {
  for (int c = 0; c < N; ++c) {
    int p = -1;
    for (int r = c; r < N; ++r)
      if (A[r * N + c]) {
        p = r;
        break;
      }
    if (p < 0)
      return -1;
    if (p != c) {
      for (int j = 0; j < N; ++j) {
        fe t = A[p * N + j];
        A[p * N + j] = A[c * N + j];
        A[c * N + j] = t;
      }
      fe t = b[p];
      b[p] = b[c];
      b[c] = t;
    }
    fe iv = ginv(F, A[c * N + c]);
    for (int j = 0; j < N; ++j)
      A[c * N + j] = gmul(F, A[c * N + j], iv);
    b[c] = gmul(F, b[c], iv);
    for (int r = 0; r < N; ++r)
      if (r != c && A[r * N + c]) {
        fe f = A[r * N + c];
        for (int j = 0; j < N; ++j)
          A[r * N + j] ^= gmul(F, f, A[c * N + j]);
        b[r] ^= gmul(F, f, b[c]);
      }
  }
  return 0;
}
/* GF(2) basis with combination tracking (vectors < 2^16, combos < 2^16) */
typedef struct {
  int has[16];
  unsigned vec[16], comb[16];
} B2;
static int hibit(unsigned v) { return 31 - __builtin_clz(v); }
/* insert v (tagged comb). returns 1 if independent; else 0 and *dep = the
 * combination of tags that adds up to zero */
static int b2_insert(B2 *b, unsigned v, unsigned comb, unsigned *dep) {
  while (v) {
    int hb = hibit(v);
    if (!b->has[hb]) {
      b->has[hb] = 1;
      b->vec[hb] = v;
      b->comb[hb] = comb;
      return 1;
    }
    v ^= b->vec[hb];
    comb ^= b->comb[hb];
  }
  if (dep)
    *dep = comb;
  return 0;
}
/* reduce v; returns the remainder (0 iff v in span), *mask = tags used */
static unsigned b2_reduce(const B2 *b, unsigned v, unsigned *mask) {
  *mask = 0;
  while (v) {
    int hb = hibit(v);
    if (!b->has[hb])
      return v;
    v ^= b->vec[hb];
    *mask ^= b->comb[hb];
  }
  return 0;
}

/* ================= Berlekamp-Massey ================= */
/* shortest LFSR for u[0..N-1]; C[0..N] = connection poly (C[0]=1); returns L */
static int bm(const GF *F, const fe *u, int N, fe *C) {
  fe *B = calloc((size_t)N + 2, sizeof(fe)),
     *Tt = calloc((size_t)N + 2, sizeof(fe));
  memset(C, 0, sizeof(fe) * (size_t)(N + 1));
  C[0] = 1;
  B[0] = 1;
  int L = 0, m = 1;
  fe b = 1;
  for (int k = 0; k < N; ++k) {
    fe d = u[k];
    for (int i = 1; i <= L; ++i)
      d ^= gmul(F, C[i], u[k - i]);
    if (!d) {
      ++m;
      continue;
    }
    fe coef = gdiv(F, d, b);
    if (2 * L <= k) {
      memcpy(Tt, C, sizeof(fe) * (size_t)(N + 1));
      for (int i = 0; i + m <= N; ++i)
        C[i + m] ^= gmul(F, coef, B[i]);
      L = k + 1 - L;
      memcpy(B, Tt, sizeof(fe) * (size_t)(N + 1));
      b = d;
      m = 1;
    } else {
      for (int i = 0; i + m <= N; ++i)
        C[i + m] ^= gmul(F, coef, B[i]);
      ++m;
    }
  }
  free(B);
  free(Tt);
  return L;
}

/* ================= generic alternant / GRS decoder =================
 * Parity checks   s_j = sum_i y_i c_i alpha_i^j   (j = 0..r-1, alpha_i != 0
 * distinct).  Corrects e errors + s erasures whenever 2e + s <= r:
 *   erasure locator Gamma, modified syndromes T = Gamma*S mod x^r, BM on
 *   T[s..r-1] -> Lambda, Psi = Gamma*Lambda, Omega = Psi*S mod x^r, Chien
 *   search, Forney  f_i = alpha_i Omega(1/alpha_i) / Psi'(1/alpha_i),
 *   error value e_i = f_i / y_i.
 * rcv[i] must be 0 at erased positions.  Returns the number of (nonzero)
 * error values, positions in epos[], values in eval[] (c_i = rcv_i + e_i);
 * -1 on decoding failure.  The result is re-verified against the syndromes. */
static int alt_decode(const GF *F, int n, const fe *al, const fe *y, int r,
                      const fe *rcv, const uint8_t *era, int *epos, fe *eval) {
  int res = -1;
  fe *S = calloc((size_t)r + 1, sizeof(fe)),
     *T = calloc((size_t)r + 1, sizeof(fe)),
     *C = calloc((size_t)r + 2, sizeof(fe));
  Poly *Gam = malloc(sizeof(Poly)), *Lam = malloc(sizeof(Poly)),
       *Psi = malloc(sizeof(Poly)), *Om = malloc(sizeof(Poly)),
       *dPsi = malloc(sizeof(Poly));
  if (r > ET_PMAXD / 2)
    die("too many parity checks (r=%d) for ET_PMAXD", r);
  for (int i = 0; i < n; ++i)
    if (rcv[i]) {
      fe t = gmul(F, y[i], rcv[i]);
      for (int j = 0; j < r; ++j) {
        S[j] ^= t;
        t = gmul(F, t, al[i]);
      }
    }
  int nz = 0;
  for (int j = 0; j < r; ++j)
    nz |= S[j];
  if (!nz) {
    res = 0;
    goto done;
  }
  int nu = 0;
  Gam->deg = 0;
  Gam->c[0] = 1;
  if (era)
    for (int i = 0; i < n; ++i)
      if (era[i]) {
        if (++nu > r)
          goto done;
        Gam->c[Gam->deg + 1] = 0;
        for (int j = Gam->deg + 1; j >= 1; --j)
          Gam->c[j] ^= gmul(F, al[i], Gam->c[j - 1]);
        Gam->deg++;
      }
  for (int j = 0; j < r; ++j) {
    fe acc = 0;
    for (int l = 0; l <= Gam->deg && l <= j; ++l)
      acc ^= gmul(F, Gam->c[l], S[j - l]);
    T[j] = acc;
  }
  int N = r - nu;
  int L = bm(F, T + nu, N, C);
  if (2 * L > N)
    goto done;
  Lam->deg = L;
  memcpy(Lam->c, C, sizeof(fe) * (size_t)(L + 1));
  ptrim(Lam);
  pmul(F, Gam, Lam, Psi);
  Om->deg = r - 1;
  for (int j = 0; j < r; ++j) {
    fe acc = 0;
    for (int l = 0; l <= Psi->deg && l <= j; ++l)
      acc ^= gmul(F, Psi->c[l], S[j - l]);
    Om->c[j] = acc;
  }
  ptrim(Om);
  pderiv(Psi, dPsi);
  int cnt = 0, nroots = 0;
  for (int i = 0; i < n; ++i) {
    fe x = ginv(F, al[i]);
    if (peval(F, Psi, x))
      continue;
    ++nroots;
    fe den = peval(F, dPsi, x);
    if (!den)
      goto done;
    fe f = gmul(F, al[i], gdiv(F, peval(F, Om, x), den));
    fe e = gdiv(F, f, y[i]);
    if (e) {
      epos[cnt] = i;
      eval[cnt++] = e;
    }
  }
  if (nroots != Psi->deg)
    goto done;
  for (int q = 0; q < cnt; ++q) { /* verify: corrected word has zero syndrome */
    fe t = gmul(F, y[epos[q]], eval[q]);
    for (int j = 0; j < r; ++j) {
      S[j] ^= t;
      t = gmul(F, t, al[epos[q]]);
    }
  }
  for (int j = 0; j < r; ++j)
    if (S[j])
      goto done;
  res = cnt;
done:
  free(S);
  free(T);
  free(C);
  free(Gam);
  free(Lam);
  free(Psi);
  free(Om);
  free(dPsi);
  return res;
}

/* ================= binary alternant-type code (bch / alt / goppa) =========
 * parity check  H = (y_i alpha_i^j), j < r, expanded over GF(2); the code is
 * the exact binary kernel.  goppa=1: y_i = 1/g(alpha_i), r = deg g, decoded
 * by the Euclidean algorithm modulo g^2. */
typedef struct {
  GF *F;
  int n, r, goppa;
  fe *al, *y, *invG2, *rb, *val;
  int *pos;
  Poly g, G2;
  int k, W, *info;
  uint64_t *gen; /* k rows of W words */
} AltCode;

static AltCode *altc_alloc(GF *F, int n, int r) {
  AltCode *A = calloc(1, sizeof *A);
  A->F = F;
  A->n = n;
  A->r = r;
  A->al = calloc((size_t)n, sizeof(fe));
  A->y = calloc((size_t)n, sizeof(fe));
  A->rb = calloc((size_t)n, sizeof(fe));
  A->val = calloc((size_t)n, sizeof(fe));
  A->pos = calloc((size_t)n, sizeof(int));
  return A;
}
static void altc_build(AltCode *A) {
  GF *F = A->F;
  int n = A->n, r = A->r, m = F->m, W = (n + 63) / 64, R = m * r;
  A->W = W;
  uint64_t *H = calloc((size_t)R * (size_t)W, sizeof(uint64_t));
  for (int i = 0; i < n; ++i) {
    fe t = A->y[i];
    for (int j = 0; j < r; ++j) {
      for (int b = 0; b < m; ++b)
        if ((t >> b) & 1)
          H[(size_t)(j * m + b) * W + i / 64] |= 1ULL << (i % 64);
      t = gmul(F, t, A->al[i]);
    }
  }
  int *piv = malloc(sizeof(int) * (size_t)(n + 1)), rank = 0;
  for (int c = 0; c < n && rank < R; ++c) {
    int p = -1;
    for (int q = rank; q < R; ++q)
      if ((H[(size_t)q * W + c / 64] >> (c % 64)) & 1) {
        p = q;
        break;
      }
    if (p < 0)
      continue;
    if (p != rank)
      for (int w = 0; w < W; ++w) {
        uint64_t t = H[(size_t)p * W + w];
        H[(size_t)p * W + w] = H[(size_t)rank * W + w];
        H[(size_t)rank * W + w] = t;
      }
    for (int q = 0; q < R; ++q)
      if (q != rank && ((H[(size_t)q * W + c / 64] >> (c % 64)) & 1))
        for (int w = c / 64; w < W; ++w)
          H[(size_t)q * W + w] ^= H[(size_t)rank * W + w];
    piv[rank++] = c;
  }
  uint8_t *isp = calloc((size_t)n, 1);
  for (int q = 0; q < rank; ++q)
    isp[piv[q]] = 1;
  A->k = n - rank;
  if (A->k <= 0)
    die("the code has dimension 0 (parameters too strong for n=%d)", n);
  A->info = malloc(sizeof(int) * (size_t)A->k);
  A->gen = calloc((size_t)A->k * (size_t)W, sizeof(uint64_t));
  int u = 0;
  for (int f = 0; f < n; ++f)
    if (!isp[f]) {
      A->info[u] = f;
      uint64_t *row = A->gen + (size_t)u * W;
      row[f / 64] |= 1ULL << (f % 64);
      for (int q = 0; q < rank; ++q)
        if ((H[(size_t)q * W + f / 64] >> (f % 64)) & 1)
          row[piv[q] / 64] |= 1ULL << (piv[q] % 64);
      ++u;
    }
  free(H);
  free(piv);
  free(isp);
}
static void altc_encode(const AltCode *A, const uint8_t *msg, uint8_t *cw) {
  int W = A->W;
  uint64_t *acc = calloc((size_t)W, sizeof(uint64_t));
  for (int u = 0; u < A->k; ++u)
    if (msg[u])
      for (int w = 0; w < W; ++w)
        acc[w] ^= A->gen[(size_t)u * W + w];
  for (int i = 0; i < A->n; ++i)
    cw[i] = (uint8_t)((acc[i / 64] >> (i % 64)) & 1);
  free(acc);
}
static int altc_decode_bm(const AltCode *A, const uint8_t *rcv,
                          const uint8_t *era, uint8_t *out) {
  int n = A->n;
  for (int i = 0; i < n; ++i) {
    int er = era && era[i];
    A->rb[i] = er ? 0 : rcv[i];
    out[i] = A->rb[i];
  }
  int c = alt_decode(A->F, n, A->al, A->y, A->r, A->rb, era, A->pos, A->val);
  if (c < 0)
    return -1;
  for (int q = 0; q < c; ++q) {
    if (A->val[q] != 1)
      return -1; /* not a binary error pattern */
    out[A->pos[q]] ^= 1;
  }
  return c;
}
/* Generic binary alternant decoder using the extended Euclidean algorithm.
 * For S(x)=sum_j s_j x^j, solve Lambda*S = Omega (mod x^r), with
 * deg(Lambda)<=floor(r/2), then Chien-search Lambda and apply Forney.
 * This is the EEA counterpart of alt_decode()'s Berlekamp-Massey path.
 * Erasures are deliberately rejected here; the BM path supports them. */
static int altc_decode_eea(const AltCode *A, const uint8_t *rcv,
                           uint8_t *out) {
  const GF *F = A->F;
  const int n = A->n, r = A->r, t = r / 2;
  if (r < 1 || r > ET_PMAXD / 2)
    return -1;
  memcpy(out, rcv, (size_t)n);
  Poly S, mod, r0, r1, t0, t1, q, rem, qt, t2, Lam, Om, dLam;
  S.deg = r - 1;
  memset(S.c, 0, sizeof(fe) * (size_t)r);
  for (int i = 0; i < n; ++i) {
    if (!rcv[i]) continue;
    fe v = gmul(F, A->y[i], rcv[i]);
    for (int j = 0; j < r; ++j) {
      S.c[j] ^= v;
      v = gmul(F, v, A->al[i]);
    }
  }
  ptrim(&S);
  if (S.deg < 0) return 0;

  /* EEA: r0=x^r, r1=S; t0=0, t1=1. */
  mod.deg = r;
  memset(mod.c, 0, sizeof(fe) * (size_t)(r + 1));
  mod.c[r] = 1;
  pcopy(&r0, &mod);
  pcopy(&r1, &S);
  t0.deg = -1;
  t1.deg = 0;
  t1.c[0] = 1;
  while (r1.deg >= (r + 1) / 2) {
    if (r1.deg < 0) return -1;
    pdivmod(F, &r0, &r1, &q, &rem);
    pmul(F, &q, &t1, &qt);
    padd(&t0, &qt, &t2);
    pcopy(&r0, &r1);
    pcopy(&r1, &rem);
    pcopy(&t0, &t1);
    pcopy(&t1, &t2);
    if (t1.deg > t) return -1;
  }
  if (t1.deg < 1 || t1.deg > t || !t1.c[0]) return -1;
  fe scale = ginv(F, t1.c[0]);
  for (int j = 0; j <= t1.deg; ++j)
    t1.c[j] = gmul(F, t1.c[j], scale);
  pcopy(&Lam, &t1);

  /* Omega = Lambda*S mod x^r. */
  pmul(F, &Lam, &S, &Om);
  if (Om.deg >= r) Om.deg = r - 1;
  ptrim(&Om);
  pderiv(&Lam, &dLam);

  int cnt = 0, roots = 0;
  for (int i = 0; i < n; ++i) {
    if (!A->al[i]) return -1;
    fe x = ginv(F, A->al[i]);
    if (peval(F, &Lam, x)) continue;
    ++roots;
    fe den = peval(F, &dLam, x);
    if (!den) return -1;
    fe num = gmul(F, A->al[i], peval(F, &Om, x));
    fe e = gdiv(F, gdiv(F, num, den), A->y[i]);
    if (e != 1) return -1; /* this module is a binary code */
    A->pos[cnt] = i;
    A->val[cnt] = e;
    out[i] ^= 1;
    ++cnt;
  }
  if (roots != Lam.deg || cnt != Lam.deg) return -1;

  /* Independent post-correction syndrome check. */
  for (int j = 0; j < r; ++j) {
    fe sum = 0;
    for (int i = 0; i < n; ++i)
      if (out[i]) sum ^= gmul(F, A->y[i], gpow(F, A->al[i], j));
    if (sum) return -1;
  }
  return cnt;
}
/* S(z) = sum_{w_i=1} (z - alpha_i)^(-1) mod g^2 */
static void goppa_syndrome(const AltCode *A, const uint8_t *w, Poly *S) {
  const GF *F = A->F;
  int d = A->G2.deg;
  fe q[ET_PMAXD + 1];
  S->deg = d - 1;
  memset(S->c, 0, sizeof(fe) * (size_t)d);
  for (int i = 0; i < A->n; ++i)
    if (w[i]) {
      q[d - 1] = A->G2.c[d];
      for (int j = d - 1; j >= 1; --j)
        q[j - 1] = (fe)(A->G2.c[j] ^ gmul(F, A->al[i], q[j]));
      for (int j = 0; j < d; ++j)
        S->c[j] ^= gmul(F, A->invG2[i], q[j]);
    }
  ptrim(S);
}
static int altc_decode_goppa(const AltCode *A, const uint8_t *rcv,
                             uint8_t *out) {
  const GF *F = A->F;
  int t = A->g.deg;
  memcpy(out, rcv, (size_t)A->n);
  Poly *S = malloc(sizeof(Poly)), *r0 = malloc(sizeof(Poly)),
       *r1 = malloc(sizeof(Poly)), *t0 = malloc(sizeof(Poly)),
       *t1 = malloc(sizeof(Poly)), *q = malloc(sizeof(Poly)),
       *rem = malloc(sizeof(Poly)), *qt = malloc(sizeof(Poly)),
       *t2 = malloc(sizeof(Poly));
  int res = -1;
  goppa_syndrome(A, rcv, S);
  if (S->deg < 0) {
    res = 0;
    goto done;
  }
  pcopy(r0, &A->G2);
  pcopy(r1, S);
  t0->deg = -1;
  t1->deg = 0;
  t1->c[0] = 1;
  while (r1->deg >= t) { /* Euclid until deg r < t: sigma = t1, omega = r1 */
    pdivmod(F, r0, r1, q, rem);
    pmul(F, q, t1, qt);
    padd(t0, qt, t2);
    pcopy(r0, r1);
    pcopy(r1, rem);
    pcopy(t0, t1);
    pcopy(t1, t2);
  }
  if (t1->deg < 1 || t1->deg > t)
    goto done;
  int cnt = 0;
  for (int i = 0; i < A->n; ++i)
    if (!peval(F, t1, A->al[i])) {
      out[i] ^= 1;
      ++cnt;
    }
  if (cnt != t1->deg)
    goto done;
  goppa_syndrome(A, out, S);
  if (S->deg >= 0)
    goto done;
  res = cnt;
done:
  free(S);
  free(r0);
  free(r1);
  free(t0);
  free(t1);
  free(q);
  free(rem);
  free(qt);
  free(t2);
  return res;
}

/* ---- constructors ---- */
static AltCode *make_bch(int m, unsigned poly, int t, int n, int b) {
  GF *F = gf_new(m, poly);
  if (n <= 0)
    n = F->ord;
  if (n > F->ord)
    die("BCH length n=%d exceeds 2^m-1=%d (use a larger --m)", n, F->ord);
  if (t < 1 || 2 * t >= n)
    die("BCH needs 1 <= t and 2t < n");
  /* n | 2^m-1, n < 2^m-1: cyclic code of length n (beta of order n);
   * otherwise primitive or shortened primitive BCH code */
  fe beta = (n < F->ord && F->ord % n == 0) ? F->exp[F->ord / n] : F->exp[1];
  AltCode *A = altc_alloc(F, n, 2 * t);
  for (int i = 0; i < n; ++i) {
    A->al[i] = gpow(F, beta, i);
    A->y[i] = gpow(F, A->al[i], b);
  }
  altc_build(A);
  return A;
}
static AltCode *make_alt(int m, unsigned poly, int t, int n) {
  GF *F = gf_new(m, poly);
  if (n <= 0)
    n = F->ord;
  if (n > F->ord)
    die("alternant length n=%d exceeds 2^m-1=%d (alpha_i must be nonzero)", n,
        F->ord);
  if (t < 1 || 2 * t >= n)
    die("alternant needs 1 <= t and 2t < n");
  AltCode *A = altc_alloc(F, n, 2 * t);
  fe *perm = malloc(sizeof(fe) * (size_t)F->ord);
  for (int i = 0; i < F->ord; ++i)
    perm[i] = F->exp[i];
  for (int i = F->ord - 1; i > 0; --i) { /* Fisher-Yates */
    int j = (int)(rnd() % (uint64_t)(i + 1));
    fe tmp = perm[i];
    perm[i] = perm[j];
    perm[j] = tmp;
  }
  for (int i = 0; i < n; ++i) {
    A->al[i] = perm[i];
    A->y[i] = grand_nz(F);
  }
  free(perm);
  altc_build(A);
  return A;
}
static AltCode *make_goppa(int m, unsigned poly, int t, int n) {
  GF *F = gf_new(m, poly);
  if (t < 1 || t > ET_PMAXD / 4)
    die("Goppa degree t must be in 1..%d", ET_PMAXD / 4);
  Poly g;
  int tries = 0;
  do { /* random monic irreducible g of degree t */
    g.deg = t;
    for (int i = 0; i < t; ++i)
      g.c[i] = (fe)(rnd() % (uint64_t)F->size);
    g.c[t] = 1;
    if (++tries > 200000)
      die("could not find an irreducible polynomial of degree %d", t);
  } while (!pirreducible(F, &g));
  fe *L = malloc(sizeof(fe) * (size_t)F->size);
  int nl = 0;
  for (int x = 0; x < F->size; ++x)
    if (peval(F, &g, (fe)x))
      L[nl++] = (fe)x;
  if (n <= 0)
    n = nl;
  if (n > nl)
    die("Goppa length n=%d exceeds |L|=%d", n, nl);
  if (n <= 2 * t)
    die("Goppa needs n > 2t");
  AltCode *A = altc_alloc(F, n, t);
  A->goppa = 1;
  pcopy(&A->g, &g);
  pmul(F, &g, &g, &A->G2);
  A->invG2 = malloc(sizeof(fe) * (size_t)n);
  for (int i = 0; i < n; ++i) {
    A->al[i] = L[i];
    fe gv = peval(F, &g, L[i]);
    A->y[i] = ginv(F, gv);
    A->invG2[i] = ginv(F, gmul(F, gv, gv));
  }
  free(L);
  altc_build(A);
  return A;
}

/* Generalized Srivastava / Srivastava binary alternant construction.
 * GS parity checks are z_i/(alpha_i-w_j)^ell, j=1..s, ell=1..mult.
 * In standard alternant form they have r=s*mult rows with
 * y_i=z_i/prod_j(alpha_i-w_j)^mult and a suitable polynomial basis.
 * The original Srivastava family is mult=1, z_i=alpha_i^mu. */
static AltCode *make_srivastava(int m, unsigned poly, int s, int mult,
                                int n, int mu, int original) {
  GF *F = gf_new(m, poly);
  if (s < 1 || mult < 1 || mult > ET_PMAXD / 2 ||
      s > (ET_PMAXD / 2) / mult)
    die("Srivastava requires s>=1, multiplicity>=1, and s*multiplicity<=%d",
        ET_PMAXD / 2);
  if (original && mult != 1)
    die("the original Srivastava construction uses multiplicity 1");
  if (s >= F->ord)
    die("need fewer poles than nonzero field elements (s < 2^m-1)");

  /* Choose distinct nonzero poles, then support points disjoint from them.
   * Nonzero alpha values are required by the current Forney/Chien decoder. */
  fe *pool = malloc(sizeof(fe) * (size_t)F->ord);
  if (!pool) die("out of memory creating Srivastava support");
  for (int i = 0; i < F->ord; ++i) pool[i] = F->exp[i];
  for (int i = F->ord - 1; i > 0; --i) {
    int j = (int)(rnd() % (uint64_t)(i + 1));
    fe z = pool[i]; pool[i] = pool[j]; pool[j] = z;
  }
  fe *w = malloc(sizeof(fe) * (size_t)s);
  if (!w) die("out of memory creating Srivastava poles");
  memcpy(w, pool, sizeof(fe) * (size_t)s);
  int maxn = F->ord - s;
  if (n <= 0) n = maxn;
  if (n < 2 || n > maxn)
    die("Srivastava length n must be in 2..%d (nonzero support excluding s poles)", maxn);
  int r = s * mult;
  if (r >= n)
    die("Srivastava requires s*multiplicity < n (got %d >= %d)", r, n);

  AltCode *A = altc_alloc(F, n, r);
  for (int i = 0; i < n; ++i) {
    fe alpha = pool[s + i];
    A->al[i] = alpha;
    fe denom = 1;
    for (int j = 0; j < s; ++j) {
      fe d = (fe)(alpha ^ w[j]);
      for (int ell = 0; ell < mult; ++ell)
        denom = gmul(F, denom, d);
    }
    fe z = original ? gpow(F, alpha, mu) : grand_nz(F);
    A->y[i] = gdiv(F, z, denom);
  }
  free(pool);
  free(w);
  altc_build(A);
  return A;
}

/* ================= Golay codes ================= */
#define GOLAY_G 0xC75u /* x^11+x^10+x^6+x^5+x^4+x^2+1 divides x^23 + 1 */
static uint32_t gol_rem(uint32_t v) {
  for (int i = 22; i >= 11; --i)
    if ((v >> i) & 1)
      v ^= GOLAY_G << (i - 11);
  return v;
}
static uint32_t *golay_table(void) {
  uint32_t *tab = malloc(2048 * sizeof(uint32_t));
  for (int i = 0; i < 2048; ++i)
    tab[i] = 0xFFFFFFFFu;
  tab[0] = 0;
  for (int a = 0; a < 23; ++a) {
    tab[gol_rem(1u << a)] = 1u << a;
    for (int b = a + 1; b < 23; ++b) {
      tab[gol_rem((1u << a) | (1u << b))] = (1u << a) | (1u << b);
      for (int c = b + 1; c < 23; ++c) {
        uint32_t e = (1u << a) | (1u << b) | (1u << c);
        tab[gol_rem(e)] = e;
      }
    }
  }
  for (int i = 0; i < 2048; ++i)
    if (tab[i] == 0xFFFFFFFFu)
      die("Golay generator polynomial check failed (not a perfect code)");
  return tab;
}

/* ================= regular LDPC code (Gallager type, soft BP decoder) ======
 * The parity-check matrix is built from a seed with a private splitmix64
 * stream, so the same (n, dv, dc, seed) always gives the same code.  Encoding
 * is systematic on the free columns of the reduced H; decoding is
 * sum-product or normalized min-sum belief propagation on the Tanner graph
 * of the original sparse H.  A bit-flip channel is fed in as LLR = +-4,
 * erasures as LLR = 0.  There is no guaranteed t: "t" below is nominal. */
#define LDPC_CLAMP 50.0
typedef struct {
  int n, m, dv, dc, k, rank;
  uint32_t seed;
  uint8_t *H;     /* m x n; reduced row-echelon form (used for encoding) */
  int *edge_var;  /* m*dc: variable of edge (check c, slot j) = c*dc+j   */
  int *var_edges; /* n*dv: edge ids of each variable                     */
  int *pivot_col, *free_col;
  uint8_t *is_pivot;
  double *llr, *q, *rmsg, *tv, *pre, *suf; /* decoder scratch */
} LdpcCode;
static int g_ldpc_iters = 50, g_ldpc_minsum = 0;
static double g_ldpc_norm = 0.8;

static uint64_t ldpc_next(uint64_t *st) {
  uint64_t z = (*st += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static void ldpc_free(LdpcCode *g) {
  if (!g) return;
  free(g->H); free(g->edge_var); free(g->var_edges); free(g->pivot_col);
  free(g->free_col); free(g->is_pivot); free(g->llr); free(g->q);
  free(g->rmsg); free(g->tv); free(g->pre); free(g->suf);
  free(g);
}
/* true if two columns of the sparse H coincide (a weight-2 codeword) */
static int ldpc_has_equal_columns(const LdpcCode *g) {
  int n = g->n, dv = g->dv;
  size_t sz = 1;
  while (sz < (size_t)2 * (size_t)n) sz <<= 1;
  int *slot = malloc(sz * sizeof(int)), *key = malloc((size_t)n * (size_t)dv * sizeof(int));
  if (!slot || !key) die("out of memory (ldpc)");
  for (size_t i = 0; i < sz; ++i) slot[i] = -1;
  int dup = 0;
  for (int v = 0; v < n && !dup; ++v) {
    int *kv = key + (size_t)v * (size_t)dv;
    for (int j = 0; j < dv; ++j) kv[j] = g->var_edges[v * dv + j] / g->dc;
    for (int a = 1; a < dv; ++a) {
      int x = kv[a], b = a - 1;
      while (b >= 0 && kv[b] > x) { kv[b + 1] = kv[b]; --b; }
      kv[b + 1] = x;
    }
    uint64_t h = 1469598103934665603ULL;
    for (int j = 0; j < dv; ++j) h = (h ^ (uint64_t)kv[j]) * 1099511628211ULL;
    size_t pos = (size_t)(h ^ (h >> 29)) & (sz - 1);
    while (slot[pos] >= 0) {
      const int *ku = key + (size_t)slot[pos] * (size_t)dv;
      int same = 1;
      for (int j = 0; j < dv; ++j) if (ku[j] != kv[j]) { same = 0; break; }
      if (same) { dup = 1; break; }
      pos = (pos + 1) & (sz - 1);
    }
    if (!dup) slot[pos] = v;
  }
  free(slot); free(key);
  return dup;
}
static LdpcCode *ldpc_new(int n, int dv, int dc, uint32_t seed) {
  if (n <= 0 || dv < 1 || dc < 2 || dv >= n || dc > n || ((long long)n * dv) % dc)
    die("invalid regular LDPC parameters (need dc>=2, dv<n, dc<=n, dc | n*dv)");
  LdpcCode *g = calloc(1, sizeof *g);
  if (!g) die("out of memory (ldpc)");
  g->n = n; g->dv = dv; g->dc = dc; g->m = (int)((long long)n * dv / dc); g->seed = seed;
  if ((double)g->m * (double)n > 4e9) die("LDPC matrix too large (m*n > 4e9)");
  size_t E = (size_t)g->m * (size_t)dc;
  g->H = calloc((size_t)g->m * (size_t)n, 1);
  g->edge_var = malloc(E * sizeof(int));
  g->var_edges = malloc((size_t)n * (size_t)dv * sizeof(int));
  g->pivot_col = malloc((size_t)g->m * sizeof(int));
  g->free_col = malloc((size_t)n * sizeof(int));
  g->is_pivot = calloc((size_t)n, 1);
  g->llr = malloc((size_t)n * sizeof(double));
  g->q = malloc(E * sizeof(double));
  g->rmsg = malloc(E * sizeof(double));
  g->tv = malloc(((size_t)dc + 1) * sizeof(double));
  g->pre = malloc(((size_t)dc + 1) * sizeof(double));
  g->suf = malloc(((size_t)dc + 1) * sizeof(double));
  int *perm = malloc((size_t)n * sizeof(int)), *vcount = calloc((size_t)n, sizeof(int));
  if (!g->H || !g->edge_var || !g->var_edges || !g->pivot_col || !g->free_col ||
      !g->is_pivot || !g->llr || !g->q || !g->rmsg || !g->tv || !g->pre ||
      !g->suf || !perm || !vcount)
    die("out of memory (ldpc)");
  uint64_t st = (uint64_t)seed * 0x2545F4914F6CDD1DULL + 12345;
  /* Regenerate (continuing the same seeded stream, so it stays reproducible)
     while two columns of H coincide: that is a weight-2 codeword, and not
     even a single bit error could then be corrected reliably. */
  for (int attempt = 0;; ++attempt) {
    memset(g->H, 0, (size_t)g->m * (size_t)n);
    memset(vcount, 0, (size_t)n * sizeof(int));
    for (int b = 0; b < dv; ++b) {
      for (int i = 0; i < n; ++i) perm[i] = i;
      for (int i = n - 1; i > 0; --i) {
        int j = (int)(ldpc_next(&st) % (uint64_t)(i + 1));
        int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
      }
      for (int r = 0; r < n / dc; ++r) {
        int c = b * (n / dc) + r;
        for (int j = 0; j < dc; ++j) {
          int v = perm[r * dc + j], e = c * dc + j;
          g->H[(size_t)c * n + v] = 1;
          g->edge_var[e] = v;
          g->var_edges[v * dv + vcount[v]++] = e;
        }
      }
    }
    if (attempt >= 200 || !ldpc_has_equal_columns(g)) break;
  }
  free(perm); free(vcount);
  /* Gauss-Jordan on H (in place): pivot columns become an identity. */
  int row = 0;
  for (int col = 0; col < n && row < g->m; ++col) {
    int p = row;
    while (p < g->m && !g->H[(size_t)p * n + col]) ++p;
    if (p == g->m) continue;
    if (p != row)
      for (int j = 0; j < n; ++j) {
        uint8_t t = g->H[(size_t)row * n + j];
        g->H[(size_t)row * n + j] = g->H[(size_t)p * n + j];
        g->H[(size_t)p * n + j] = t;
      }
    for (int r = 0; r < g->m; ++r)
      if (r != row && g->H[(size_t)r * n + col])
        for (int j = col; j < n; ++j) g->H[(size_t)r * n + j] ^= g->H[(size_t)row * n + j];
    g->pivot_col[row] = col; g->is_pivot[col] = 1; ++row;
  }
  g->rank = row; g->k = n - row;
  int z = 0;
  for (int c = 0; c < n; ++c) if (!g->is_pivot[c]) g->free_col[z++] = c;
  return g;
}
static void ldpc_encode(const LdpcCode *g, const uint8_t *msg, uint8_t *word) {
  memset(word, 0, (size_t)g->n);
  for (int i = 0; i < g->k; ++i) word[g->free_col[i]] = msg[i] & 1;
  for (int r = 0; r < g->rank; ++r) {
    const uint8_t *h = g->H + (size_t)r * g->n;
    uint8_t x = 0;
    for (int i = 0; i < g->k; ++i) x ^= (uint8_t)(h[g->free_col[i]] & word[g->free_col[i]]);
    word[g->pivot_col[r]] = x;
  }
}
static int ldpc_syndrome_ok(const LdpcCode *g, const uint8_t *w) {
  for (int c = 0; c < g->m; ++c) {
    uint8_t x = 0;
    for (int j = 0; j < g->dc; ++j) x ^= w[g->edge_var[c * g->dc + j]];
    if (x) return 0;
  }
  return 1;
}
/* g->llr holds the channel LLRs; hard decisions go to out[0..n).  Returns 1
 * if a zero syndrome was reached within the iteration limit. */
static int ldpc_decode_llr(LdpcCode *g, uint8_t *out, int *used) {
  int E = g->m * g->dc;
  double *q = g->q, *rmsg = g->rmsg, *llr = g->llr;
  for (int e = 0; e < E; ++e) { q[e] = llr[g->edge_var[e]]; rmsg[e] = 0.0; }
  for (int it = 1; it <= g_ldpc_iters; ++it) {
    for (int c = 0; c < g->m; ++c) {
      int base = c * g->dc, dc = g->dc;
      if (!g_ldpc_minsum) {
        g->pre[0] = 1.0; g->suf[dc] = 1.0;
        for (int l = 0; l < dc; ++l) {
          g->tv[l] = tanh(q[base + l] * 0.5);
          g->pre[l + 1] = g->pre[l] * g->tv[l];
        }
        for (int l = dc - 1; l >= 0; --l) g->suf[l] = g->suf[l + 1] * g->tv[l];
        for (int j = 0; j < dc; ++j) {
          double prod = g->pre[j] * g->suf[j + 1];
          if (prod > 0.999999999999) prod = 0.999999999999;
          if (prod < -0.999999999999) prod = -0.999999999999;
          rmsg[base + j] = 2.0 * atanh(prod);
        }
      } else {
        int s = 1, a1 = -1;
        double m1 = 1e100, m2 = 1e100;
        for (int l = 0; l < dc; ++l) {
          double a = fabs(q[base + l]);
          if (q[base + l] < 0) s = -s;
          if (a < m1) { m2 = m1; m1 = a; a1 = l; } else if (a < m2) m2 = a;
        }
        for (int j = 0; j < dc; ++j) {
          int sj = q[base + j] < 0 ? -s : s;
          rmsg[base + j] = sj * g_ldpc_norm * (j == a1 ? m2 : m1);
        }
      }
    }
    for (int v = 0; v < g->n; ++v) {
      double sum = llr[v];
      for (int j = 0; j < g->dv; ++j) sum += rmsg[g->var_edges[v * g->dv + j]];
      if (sum > LDPC_CLAMP) sum = LDPC_CLAMP;
      if (sum < -LDPC_CLAMP) sum = -LDPC_CLAMP;
      out[v] = (uint8_t)(sum < 0);
      for (int j = 0; j < g->dv; ++j) {
        int e = g->var_edges[v * g->dv + j];
        double x = sum - rmsg[e];
        if (x > LDPC_CLAMP) x = LDPC_CLAMP;
        if (x < -LDPC_CLAMP) x = -LDPC_CLAMP;
        q[e] = x;
      }
    }
    if (ldpc_syndrome_ok(g, out)) { if (used) *used = it; return 1; }
  }
  if (used) *used = g_ldpc_iters;
  return 0;
}

/* ================= unified binary code object ================= */
static int g_fht =
    0; /* -fht: RM(1,m) is decoded by the Walsh-Hadamard transform */
static int g_alt_eea = 0; /* -eea: generic alternant EEA instead of BM */
typedef enum {
  BC_ALT = 0,
  BC_GOPPA = 1,
  BC_GOLAY = 2,
  BC_RM = 3,
  BC_HAD = 4,
  BC_LDPC = 5,
  BC_DUP = 6
} BcKind;
typedef struct {
  char name[128];
  BcKind kind;
  int n, k, t, dmin, ext;
  int *info; /* message bit j sits at codeword position info[j] */
  AltCode *alt;
  uint32_t *gtab;
  int rm_r, rm_m;
  int dup_count; /* number of identical RM(1,7) copies */
  unsigned *mon; /* RM: monomial x_S (bit mask S), degree <= r */
  uint8_t *gsys; /* RM: systematic generator, k rows of n bytes */
  LdpcCode *ldpc; /* LDPC: Tanner graph, reduced H, decoder scratch */
} BinCode;

/* Release a binary-code object and all storage it owns.  For alternant/Goppa
   objects B->info aliases A->info, so that array is freed exactly once. */
static void bc_free(BinCode *B) {
  if (!B) return;
  if (B->alt) {
    AltCode *A = B->alt;
    if (A->F) {
      free(A->F->exp);
      free(A->F->log);
      free(A->F);
    }
    free(A->al);
    free(A->y);
    free(A->invG2);
    free(A->rb);
    free(A->val);
    free(A->pos);
    free(A->info);
    free(A->gen);
    free(A);
  } else {
    free(B->info);
  }
  free(B->gtab);
  free(B->mon);
  free(B->gsys);
  ldpc_free(B->ldpc);
  free(B);
}

static BinCode *bc_from_alt(AltCode *A, const char *name, int t) {
  BinCode *B = calloc(1, sizeof *B);
  snprintf(B->name, sizeof B->name, "%s", name);
  B->kind = A->goppa ? BC_GOPPA : BC_ALT;
  B->n = A->n;
  B->k = A->k;
  B->t = t;
  B->dmin = 2 * t + 1;
  B->info = A->info;
  B->alt = A;
  return B;
}
static BinCode *bc_ldpc(int n, int dv, int dc, int t, uint32_t seed) {
  LdpcCode *g = ldpc_new(n, dv, dc, seed);
  if (g->k <= 0) die("LDPC code has no information bits");
  BinCode *B = calloc(1, sizeof *B);
  if (!B) die("out of memory");
  snprintf(B->name, sizeof B->name,
           "LDPC[%d,%d] regular dv=%d dc=%d seed=%u (nominal t, no guaranteed d)",
           g->n, g->k, dv, dc, (unsigned)seed);
  B->kind = BC_LDPC;
  B->n = g->n;
  B->k = g->k;
  B->t = t;
  B->dmin = 1;
  B->info = malloc((size_t)g->k * sizeof(int));
  if (!B->info) die("out of memory");
  memcpy(B->info, g->free_col, (size_t)g->k * sizeof(int));
  B->ldpc = g;
  return B;
}
static BinCode *bc_golay(int ext) {
  BinCode *B = calloc(1, sizeof *B);
  B->kind = BC_GOLAY;
  B->ext = ext;
  B->n = ext ? 24 : 23;
  B->k = 12;
  B->t = 3;
  B->dmin = ext ? 8 : 7;
  snprintf(B->name, sizeof B->name, "Golay[%d,12,%d]", B->n, B->dmin);
  B->info = malloc(sizeof(int) * 12);
  for (int j = 0; j < 12; ++j)
    B->info[j] = 11 + j;
  B->gtab = golay_table();
  return B;
}
/* Reed-Muller RM(r,m): point p in [0,2^m) has coordinates = bits of p,
 * monomial x_S takes the value 1 at p iff (p & S) == S. */
static BinCode *bc_rm(int r, int m) {
  if (m < 1 || m > MAXM || r < 0 || r > m)
    die("Reed-Muller needs 0 <= r <= m <= %d", MAXM);
  int n = 1 << m, k = 0;
  BinCode *B = calloc(1, sizeof *B);
  B->kind = BC_RM;
  B->rm_r = r;
  B->rm_m = m;
  B->n = n;
  for (int S = 0; S < n; ++S)
    if (__builtin_popcount(S) <= r)
      ++k;
  B->k = k;
  B->dmin = 1 << (m - r);
  B->t = (B->dmin - 1) / 2;
  snprintf(B->name, sizeof B->name, "RM(%d,%d)[%d,%d,%d]", r, m, n, k, B->dmin);
  B->mon = malloc(sizeof(unsigned) * (size_t)k);
  int u = 0;
  for (int d = r; d >= 0; --d) /* highest degree first */
    for (int S = 0; S < n; ++S)
      if (__builtin_popcount(S) == d)
        B->mon[u++] = (unsigned)S;
  uint8_t *G = malloc((size_t)k * (size_t)n);
  for (int i = 0; i < k; ++i)
    for (int pt = 0; pt < n; ++pt)
      G[(size_t)i * n + pt] = ((unsigned)pt & B->mon[i]) == B->mon[i];
  B->info = malloc(sizeof(int) * (size_t)k);
  int rank = 0;
  for (int c = 0; c < n && rank < k; ++c) { /* Gauss-Jordan -> systematic */
    int pv = -1;
    for (int q = rank; q < k; ++q)
      if (G[(size_t)q * n + c]) {
        pv = q;
        break;
      }
    if (pv < 0)
      continue;
    if (pv != rank)
      for (int j = 0; j < n; ++j) {
        uint8_t t = G[(size_t)pv * n + j];
        G[(size_t)pv * n + j] = G[(size_t)rank * n + j];
        G[(size_t)rank * n + j] = t;
      }
    for (int q = 0; q < k; ++q)
      if (q != rank && G[(size_t)q * n + c])
        for (int j = 0; j < n; ++j)
          G[(size_t)q * n + j] ^= G[(size_t)rank * n + j];
    B->info[rank++] = c;
  }
  if (rank != k)
    die("Reed-Muller generator has rank %d != %d", rank, k);
  B->gsys = G;
  return B;
}
/* Duplicate one RM(1,7) codeword several times (HQC-style dRM). */
static BinCode *bc_drm(int copies) {
  if (copies < 2 || copies > 16)
    die("dRM repetition count must be between 2 and 16");
  BinCode *B = bc_rm(1, 7);
  const int base_n = B->n;
  uint8_t *G = malloc((size_t)B->k * (size_t)base_n * (size_t)copies);
  if (!G) die("out of memory");
  for (int u = 0; u < B->k; ++u)
    for (int c = 0; c < copies; ++c)
      memcpy(G + (size_t)u * (base_n * copies) + (size_t)c * base_n,
             B->gsys + (size_t)u * base_n, (size_t)base_n);
  free(B->gsys);
  B->gsys = G;
  B->kind = BC_DUP;
  B->dup_count = copies;
  B->n = base_n * copies;
  B->dmin *= copies;
  B->t = (B->dmin - 1) / 2;
  snprintf(B->name, sizeof B->name, "dRM(1,7)x%d[%d,%d,%d]",
           copies, B->n, B->k, B->dmin);
  return B;
}
/* Reed majority-logic decoding: for d = r..0 every coefficient a_S, |S| = d,
 * is the majority of the 2^(m-d) disjoint check sums  sum_{x_S free} w(x)
 * (the other coordinates fixed); decoded degree-d terms are then removed. */
static int rm_decode(const BinCode *B, const uint8_t *rcv, uint8_t *out) {
  int n = B->n, m = B->rm_m, k = B->k;
  uint8_t *w = malloc((size_t)n), *a = calloc((size_t)k, 1);
  memcpy(w, rcv, (size_t)n);
  int res = -1, i = 0;
  for (int d = B->rm_r; d >= 0; --d) {
    int i0 = i;
    for (; i < k && __builtin_popcount(B->mon[i]) == d; ++i) {
      unsigned S = B->mon[i], C = ~S & (unsigned)(n - 1);
      long ones = 0, checks = 1L << (m - d);
      for (unsigned t = C;; t = (t - 1) & C) {
        int par = 0;
        for (unsigned s = S;; s = (s - 1) & S) {
          par ^= w[t | s];
          if (!s)
            break;
        }
        ones += par;
        if (!t)
          break;
      }
      if (2 * ones == checks)
        goto done; /* tie: failure */
      a[i] = 2 * ones > checks;
    }
    for (int q = i0; q < i; ++q)
      if (a[q])
        for (int pt = 0; pt < n; ++pt)
          if (((unsigned)pt & B->mon[q]) == B->mon[q])
            w[pt] ^= 1;
  }
  memset(out, 0, (size_t)n);
  for (int q = 0; q < k; ++q)
    if (a[q])
      for (int pt = 0; pt < n; ++pt)
        if (((unsigned)pt & B->mon[q]) == B->mon[q])
          out[pt] ^= 1;
  int cnt = 0;
  for (int pt = 0; pt < n; ++pt)
    cnt += out[pt] != rcv[pt];
  if (cnt <= B->t)
    res = cnt;
done:
  free(w);
  free(a);
  return res;
}
/* in-place fast Walsh-Hadamard transform (n = 2^m), n log n additions */
static void fwht(int *a, int n) {
  for (int h = 1; h < n; h <<= 1)
    for (int i = 0; i < n; i += h << 1)
      for (int j = i; j < i + h; ++j) {
        int x = a[j], y = a[j + h];
        a[j] = x + y;
        a[j + h] = x - y;
      }
}
/* ML decoding of first-order words  c(x) = a0 + a.x  via
 *   F[a] = sum_x (-1)^(r(x)+a.x) = n - 2 d(r, a.x):
 * the largest |F| (RM(1,m)) or largest F (Hadamard code, a0 = 0 forced) gives
 * the nearest codeword; ties = detected failure.  Returns distance or -1. */
static int fht_decode(const BinCode *B, const uint8_t *rcv, uint8_t *out,
                      int allow_complement) {
  int n = B->n;
  int *F = malloc(sizeof(int) * (size_t)n);
  for (int x = 0; x < n; ++x)
    F[x] = rcv[x] ? -1 : 1;
  fwht(F, n);
  int best = -1, bv = -1, ties = 0;
  for (int a = 0; a < n; ++a) {
    int v = allow_complement ? abs(F[a]) : F[a];
    if (v > bv) {
      bv = v;
      best = a;
      ties = 1;
    } else if (v == bv)
      ++ties;
  }
  int a0 = allow_complement && F[best] < 0;
  int cnt = 0;
  for (int x = 0; x < n; ++x) {
    out[x] = (uint8_t)((__builtin_popcount((unsigned)(best & x)) & 1) ^ a0);
    cnt += out[x] != rcv[x];
  }
  free(F);
  return ties > 1 ? -1 : cnt;
}
/* Soft-combine repeated copies, then ML-decode RM(1,7) with a 128-point
 * Walsh-Hadamard transform. */
static int drm_fht_decode(const BinCode *B, const uint8_t *rcv,
                          uint8_t *out) {
  const int base_n = 128;
  int *F = malloc(sizeof(int) * base_n);
  if (!F) return -1;
  for (int x = 0; x < base_n; ++x) {
    int sum = 0;
    for (int c = 0; c < B->dup_count; ++c)
      sum += (rcv[(size_t)c * base_n + x] & 1) ? -1 : 1;
    F[x] = sum;
  }
  fwht(F, base_n);
  int best = -1, bv = -1, ties = 0;
  for (int a = 0; a < base_n; ++a) {
    int v = abs(F[a]);
    if (v > bv) { bv = v; best = a; ties = 1; }
    else if (v == bv) ++ties;
  }
  int a0 = F[best] < 0;
  int dist = 0;
  for (int c = 0; c < B->dup_count; ++c)
    for (int x = 0; x < base_n; ++x) {
      uint8_t bit = (uint8_t)((__builtin_popcount((unsigned)(best & x)) & 1) ^ a0);
      out[(size_t)c * base_n + x] = bit;
      dist += bit != (rcv[(size_t)c * base_n + x] & 1);
    }
  free(F);
  return ties > 1 ? -1 : dist;
}
static BinCode *bc_had(int m) {
  if (m < 2 || m > MAXM)
    die("Hadamard code needs 2 <= m <= %d", MAXM);
  BinCode *B = calloc(1, sizeof *B);
  B->kind = BC_HAD;
  B->n = 1 << m;
  B->k = m;
  B->rm_m = m;
  B->dmin = 1 << (m - 1);
  B->t = (B->dmin - 1) / 2;
  snprintf(B->name, sizeof B->name, "Hadamard[%d,%d,%d]", B->n, m, B->dmin);
  B->info = malloc(sizeof(int) * (size_t)m);
  for (int j = 0; j < m; ++j)
    B->info[j] = 1 << j; /* c(e_j) = a_j: systematic */
  return B;
}
static void bc_encode(const BinCode *B, const uint8_t *msg, uint8_t *cw) {
  if (B->kind == BC_LDPC) {
    ldpc_encode(B->ldpc, msg, cw);
  } else if (B->kind == BC_GOLAY) {
    uint32_t m = 0;
    for (int j = 0; j < 12; ++j)
      m |= (uint32_t)(msg[j] & 1) << j;
    uint32_t v = m << 11;
    v |= gol_rem(v);
    for (int i = 0; i < 23; ++i)
      cw[i] = (uint8_t)((v >> i) & 1);
    if (B->ext)
      cw[23] = (uint8_t)(__builtin_popcount(v) & 1);
  } else if (B->kind == BC_HAD) {
    unsigned a = 0;
    for (int j = 0; j < B->k; ++j)
      a |= (unsigned)(msg[j] & 1) << j;
    for (int x = 0; x < B->n; ++x)
      cw[x] = (uint8_t)(__builtin_popcount(a & (unsigned)x) & 1);
  } else if (B->kind == BC_RM || B->kind == BC_DUP) {
    memset(cw, 0, (size_t)B->n);
    for (int u = 0; u < B->k; ++u)
      if (msg[u])
        for (int j = 0; j < B->n; ++j)
          cw[j] ^= B->gsys[(size_t)u * B->n + j];
  } else
    altc_encode(B->alt, msg, cw);
}
/* corrected-error count, or -1 (detected failure).  era may be NULL;
 * erasures are supported by BCH/alternant codes only. */
static int bc_decode(const BinCode *B, const uint8_t *rcv, const uint8_t *era,
                     uint8_t *out) {
  if (B->kind == BC_LDPC) {
    LdpcCode *g = B->ldpc;
    for (int i = 0; i < g->n; ++i)
      g->llr[i] = (era && era[i]) ? 0.0 : ((rcv[i] & 1) ? -4.0 : 4.0);
    if (!ldpc_decode_llr(g, out, NULL)) return -1;
    int cnt = 0;
    for (int i = 0; i < g->n; ++i) cnt += (out[i] != (rcv[i] & 1));
    return cnt;
  }
  if (B->kind == BC_GOLAY) {
    uint32_t v = 0;
    for (int i = 0; i < 23; ++i)
      v |= (uint32_t)(rcv[i] & 1) << i;
    uint32_t e = B->gtab[gol_rem(v)], c = v ^ e;
    int cnt = __builtin_popcount(e);
    for (int i = 0; i < 23; ++i)
      out[i] = (uint8_t)((c >> i) & 1);
    if (B->ext) {
      int par = __builtin_popcount(c) & 1;
      out[23] = (uint8_t)par;
      cnt += (rcv[23] & 1) != par;
      if (cnt > 3)
        return -1; /* 4 errors: detected */
    }
    return cnt;
  }
  if (B->kind == BC_HAD)
    return fht_decode(B, rcv, out, 0);
  if (B->kind == BC_DUP)
    return drm_fht_decode(B, rcv, out);
  if (B->kind == BC_RM)
    return (g_fht && B->rm_r == 1) ? fht_decode(B, rcv, out, 1)
                                   : rm_decode(B, rcv, out);
  if (B->kind == BC_GOPPA)
    return altc_decode_goppa(B->alt, rcv, out);
  if (g_alt_eea) {
    if (era) return -1;
    return altc_decode_eea(B->alt, rcv, out);
  }
  return altc_decode_bm(B->alt, rcv, era, out);
}

/* parse BCH/Goppa/Alternant and Srivastava-family specifications */
static BinCode *bc_parse(const char *spec, unsigned poly) {
  char buf[128], *tok[8];
  int nt = 0;
  snprintf(buf, sizeof buf, "%s", spec);
  for (char *p = strtok(buf, ":"); p && nt < 8; p = strtok(NULL, ":"))
    tok[nt++] = p;
  if (!nt)
    die("empty code spec");
  int v[8] = {0};
  for (int i = 1; i < nt; ++i)
    v[i] = atoi(tok[i]);
  char name[128];
  if (!strcmp(tok[0], "golay"))
    return bc_golay(0);
  if (!strcmp(tok[0], "golay24"))
    return bc_golay(1);
  if (!strcmp(tok[0], "had") && nt >= 2)
    return bc_had(v[1]);
  if (!strcmp(tok[0], "ldpc") && nt >= 4)
    return bc_ldpc(v[1], v[2], v[3], nt > 4 && v[4] > 0 ? v[4] : 2,
                   nt > 5 ? (uint32_t)strtoul(tok[5], NULL, 10)
                          : ((uint32_t)(rnd() >> 16) | 1u));
  if (!strcmp(tok[0], "rm") && nt >= 3)
    return bc_rm(v[1], v[2]);
  if (!strcmp(tok[0], "drm") && nt >= 2)
    return bc_drm(v[1]);
  if (!strcmp(tok[0], "bch") && nt >= 3) {
    int n = nt > 3 ? v[3] : 0, b = nt > 4 ? v[4] : 1;
    AltCode *A = make_bch(v[1], poly, v[2], n, b);
    snprintf(name, sizeof name, "BCH[%d,%d,>=%d] GF(2^%d) b=%d", A->n, A->k,
             2 * v[2] + 1, v[1], b);
    return bc_from_alt(A, name, v[2]);
  }
  if (!strcmp(tok[0], "goppa") && nt >= 3) {
    AltCode *A = make_goppa(v[1], poly, v[2], nt > 3 ? v[3] : 0);
    snprintf(name, sizeof name, "Goppa[%d,%d,>=%d] GF(2^%d) deg g=%d", A->n,
             A->k, 2 * v[2] + 1, v[1], v[2]);
    return bc_from_alt(A, name, v[2]);
  }
  if (!strcmp(tok[0], "alt") && nt >= 4) {
    AltCode *A = make_alt(v[1], poly, v[3], v[2]);
    snprintf(name, sizeof name, "alternant[%d,>=%d,>=%d] GF(2^%d)", A->n, A->k,
             2 * v[3] + 1, v[1]);
    return bc_from_alt(A, name, v[3]);
  }
  if (!strcmp(tok[0], "srivastava") && nt >= 3) {
    int mm=v[1], ss=v[2], nn=nt>3?v[3]:0, mu=nt>4?v[4]:0;
    AltCode *A=make_srivastava(mm,poly,ss,1,nn,mu,1);
    int r=ss, terr=r/2;
    snprintf(name,sizeof name,"Srivastava[%d,%d,>=%d] GF(2^%d) s=%d mu=%d",
             A->n,A->k,r+1,mm,ss,mu);
    BinCode *B=bc_from_alt(A,name,terr);
    B->dmin=r+1;
    return B;
  }
  if (!strcmp(tok[0], "gs") && nt >= 4) {
    int mm=v[1], ss=v[2], mult=v[3], nn=nt>4?v[4]:0;
    AltCode *A=make_srivastava(mm,poly,ss,mult,nn,0,0);
    int r=ss*mult, terr=r/2;
    snprintf(name,sizeof name,"GS[%d,%d,>=%d] GF(2^%d) s=%d t=%d",
             A->n,A->k,r+1,mm,ss,mult);
    BinCode *B=bc_from_alt(A,name,terr);
    B->dmin=r+1;
    return B;
  }
  die("bad code spec '%s' (bch:m:t[:n[:b]] | golay | golay24 | goppa:m:t[:n] | "
      "alt:m:n:t | srivastava:m:s[:n[:mu]] | gs:m:s:t[:n] | rm:r:m | drm:copies | had:m | "
      "ldpc:n:dv:dc[:t[:seed]])",
      spec);
  return NULL;
}


/* ================= same-rate benchmark =================
 * Compare binary code families supported by bc_parse().  Each candidate is
 * tested with the same PRNG initial state; timings are per block. */
static int benchmark_rate(const char *specs, double target, double tol,
                          int trials, int fixed_errors, uint64_t seed,
                          unsigned poly) {
  if (!specs || !*specs || !(target > 0.0 && target <= 1.0) ||
      tol < 0.0 || trials < 1) {
    fprintf(stderr, "invalid --compare / --rate / --rate-tol / -n value\n");
    return 2;
  }
  char *list = strdup(specs);
  if (!list) die("out of memory");
  uint64_t saved_rng = rng_s;
  uint64_t base_rng = saved_rng ^ (seed * 0x9E3779B97F4A7C15ULL);
  printf("# same-rate benchmark: target=%.6f tolerance=%.6f trials=%d\n",
         target, tol, trials);
  printf("# spec | code | n | k | rate | enc_ms/block | dec_ms/block | failures | fail_rate\n");
  int tested = 0;
  char *p = list;
  while (*p) {
    char *comma = strchr(p, ',');
    if (comma) *comma = '\0';
    while (*p == ' ' || *p == '\t') ++p;
    char *end = p + strlen(p);
    while (end > p && (end[-1] == ' ' || end[-1] == '\t')) *--end = '\0';
    if (*p) {
      rng_s = base_rng;
      BinCode *B = bc_parse(p, poly);
      double rate = (double)B->k / (double)B->n;
      if (rate < target - tol || rate > target + tol) {
        printf("# skip %s: rate %.6f outside [%.6f, %.6f]\n",
               p, rate, target - tol, target + tol);
      } else {
        int errors = fixed_errors >= 0 ? fixed_errors : B->t;
        if (errors < 0 || errors > B->n) {
          fprintf(stderr, "invalid error count %d for %s (n=%d)\n",
                  errors, p, B->n);
          bc_free(B); free(list); rng_s = saved_rng; return 2;
        }
        uint8_t *msg = malloc((size_t)B->k);
        uint8_t *cw = malloc((size_t)B->n);
        uint8_t *rv = malloc((size_t)B->n);
        uint8_t *dec = malloc((size_t)B->n);
        int *perm = malloc(sizeof(int) * (size_t)B->n);
        if (!msg || !cw || !rv || !dec || !perm)
          die("out of memory in rate benchmark");
        uint64_t trial_rng = rng_s;
        double enc_total = 0.0, dec_total = 0.0;
        int failures = 0;
        for (int tr = 0; tr < trials; ++tr) {
          for (int j = 0; j < B->k; ++j)
            msg[j] = (uint8_t)(rnd() & 1u);
          double t0 = now();
          bc_encode(B, msg, cw);
          enc_total += now() - t0;
          memcpy(rv, cw, (size_t)B->n);
          for (int j = 0; j < B->n; ++j) perm[j] = j;
          for (int j = 0; j < errors; ++j) {
            int q = j + (int)(rnd() % (uint64_t)(B->n - j));
            int z = perm[j]; perm[j] = perm[q]; perm[q] = z;
            rv[perm[j]] ^= 1;
          }
          t0 = now();
          int rc = bc_decode(B, rv, NULL, dec);
          dec_total += now() - t0;
          if (rc < 0 || memcmp(dec, cw, (size_t)B->n) != 0)
            ++failures;
        }
        rng_s = trial_rng;
        printf("%s | %s | %d | %d | %.6f | %.6f | %.6f | %d/%d | %.6f\n",
               p, B->name, B->n, B->k, rate,
               1000.0 * enc_total / trials, 1000.0 * dec_total / trials,
               failures, trials, (double)failures / trials);
        ++tested;
        free(msg); free(cw); free(rv); free(dec); free(perm);
      }
      bc_free(B);
    }
    if (!comma) break;
    p = comma + 1;
  }
  free(list);
  rng_s = saved_rng;
  if (!tested)
    fprintf(stderr, "no candidate matched the requested rate window\n");
  return tested ? 0 : 1;
}

/* ================= automatic exact same-rate search =================
 * Enumerate bounded RM/dRM/BCH parameter sets, group by exact k/n using
 * cross multiplication, then benchmark every candidate at the selected rate.
 * If --errors is omitted, use the smallest guaranteed correction capability
 * among the selected candidates so every code receives the same error load. */
typedef struct {
  char spec[64];
  int n, k, t, family; /* 0=dRM, 1=RM, 2=BCH */
} AutoRateCand;

static int auto_rate_equal(const AutoRateCand *a, const AutoRateCand *b) {
  return (int64_t)a->k * b->n == (int64_t)b->k * a->n;
}

static int auto_add_candidate(AutoRateCand *v, int *nv, int cap,
                              const char *spec, int family, unsigned poly) {
  if (*nv >= cap) return -1;
  BinCode *B = bc_parse(spec, poly);
  AutoRateCand *c = &v[(*nv)++];
  snprintf(c->spec, sizeof c->spec, "%s", spec);
  c->n = B->n; c->k = B->k; c->t = B->t; c->family = family;
  bc_free(B);
  return 0;
}

static int auto_same_rate(const char *families, double target, double tol,
                          int trials, int fixed_errors, uint64_t seed,
                          unsigned poly) {
  enum { AUTO_CAP = 512 };
  AutoRateCand cand[AUTO_CAP];
  int nc = 0, want[3] = {0, 0, 0};
  char fbuf[128];
  if (!families || !*families || strlen(families) >= sizeof fbuf) {
    fprintf(stderr, "--auto-same-rate needs a family list: drm,rm,bch\n");
    return 2;
  }
  snprintf(fbuf, sizeof fbuf, "%s", families);
  for (char *q = strtok(fbuf, ","); q; q = strtok(NULL, ",")) {
    while (*q == ' ' || *q == '\t') ++q;
    char *e = q + strlen(q);
    while (e > q && (e[-1] == ' ' || e[-1] == '\t')) *--e = '\0';
    if (!strcmp(q, "drm")) want[0] = 1;
    else if (!strcmp(q, "rm")) want[1] = 1;
    else if (!strcmp(q, "bch")) want[2] = 1;
    else {
      fprintf(stderr, "unsupported auto-search family '%s' (use drm,rm,bch)\n", q);
      return 2;
    }
  }
  if (want[0] + want[1] + want[2] < 2) {
    fprintf(stderr, "--auto-same-rate needs at least two distinct families\n");
    return 2;
  }
  int requested_families = want[0] + want[1] + want[2];
  uint64_t saved_rng = rng_s;
  uint64_t base_rng = saved_rng ^ (seed * 0x9E3779B97F4A7C15ULL);
  char spec[64];
  if (want[0]) for (int c = 2; c <= 16; ++c) {
    snprintf(spec, sizeof spec, "drm:%d", c);
    rng_s = base_rng;
    if (auto_add_candidate(cand, &nc, AUTO_CAP, spec, 0, poly)) goto oom;
  }
  if (want[1]) for (int m = 2; m <= MAXM; ++m)
    for (int r = 0; r < m; ++r) {
      snprintf(spec, sizeof spec, "rm:%d:%d", r, m);
      rng_s = base_rng;
      if (auto_add_candidate(cand, &nc, AUTO_CAP, spec, 1, poly)) goto oom;
    }
  if (want[2]) for (int m = 3; m <= MAXM; ++m) {
    int n = (1 << m) - 1;
    int tmax = (n - 1) / 2;
    if (tmax > 32) tmax = 32; /* bounded search; keeps enumeration practical */
    for (int t = 1; t <= tmax; ++t) {
      snprintf(spec, sizeof spec, "bch:%d:%d", m, t);
      rng_s = base_rng;
      if (auto_add_candidate(cand, &nc, AUTO_CAP, spec, 2, poly)) goto oom;
    }
  }
  rng_s = saved_rng;
  if (nc == 0) {
    fprintf(stderr, "auto-search generated no candidates\n");
    return 1;
  }

  int best = -1, best_coverage = 0;
  double best_delta = 1e300, best_rate = 0.0;
  for (int i = 0; i < nc; ++i) {
    int seen[3] = {0, 0, 0};
    for (int j = 0; j < nc; ++j)
      if (auto_rate_equal(&cand[i], &cand[j])) seen[cand[j].family] = 1;
    int coverage = seen[0] + seen[1] + seen[2];
    if (coverage < 2) continue;
    double rate = (double)cand[i].k / cand[i].n;
    double delta = target > 0.0 ? fabs(rate - target) : 0.0;
    if (target > 0.0 && delta > tol) continue;
    if (coverage > best_coverage ||
        (coverage == best_coverage && target > 0.0 && delta < best_delta) ||
        (coverage == best_coverage && (target <= 0.0 || delta == best_delta) && rate > best_rate)) {
      best = i; best_coverage = coverage; best_delta = delta; best_rate = rate;
    }
  }
  if (best < 0) {
    fprintf(stderr, "no exact common rate found for requested families\n");
    return 1;
  }

  int min_t = INT32_MAX, selected = 0;
  size_t listcap = 1, used = 0;
  for (int i = 0; i < nc; ++i) if (auto_rate_equal(&cand[best], &cand[i])) {
    listcap += strlen(cand[i].spec) + 1;
    if (cand[i].t < min_t) min_t = cand[i].t;
    ++selected;
  }
  char *list = malloc(listcap);
  if (!list) goto oom;
  list[0] = '\0';
  for (int i = 0; i < nc; ++i) if (auto_rate_equal(&cand[best], &cand[i])) {
    if (used) list[used++] = ',';
    size_t z = strlen(cand[i].spec);
    memcpy(list + used, cand[i].spec, z); used += z; list[used] = '\0';
  }
  int errors = fixed_errors >= 0 ? fixed_errors : min_t;
  printf("# auto same-rate search: families=%s candidates=%d selected=%d\n",
         families, nc, selected);
  printf("# exact common rate = %d/%d = %.9f; family coverage=%d/%d\n",
         cand[best].k, cand[best].n, best_rate, best_coverage, requested_families);
  printf("# common error load = %d bit(s) per block%s\n", errors,
         fixed_errors >= 0 ? " (user-specified)" : " (minimum guaranteed capability)");
  int rc = benchmark_rate(list, best_rate, 0.0, trials, errors, seed, poly);
  free(list);
  return rc;

oom:
  rng_s = saved_rng;
  fprintf(stderr, "out of memory or too many auto-search candidates\n");
  return 2;
}

/* ================= GRS ================= */
typedef struct {
  GF *F;
  int n, k;
  fe *al, *v, *y, *rb, *val;
  int *pos;
  int found; /* number of error values of the last decode */
} GRS;

static GRS *grs_new(GF *F, int n, int k, int shuffle, int randv) {
  if (n > F->ord)
    die("GRS length n=%d exceeds 2^m-1=%d", n, F->ord);
  if (k < 1 || k >= n)
    die("GRS needs 1 <= k < n");
  GRS *G = calloc(1, sizeof *G);
  G->F = F;
  G->n = n;
  G->k = k;
  G->al = malloc(sizeof(fe) * (size_t)n);
  G->v = malloc(sizeof(fe) * (size_t)n);
  G->y = malloc(sizeof(fe) * (size_t)n);
  G->rb = malloc(sizeof(fe) * (size_t)n);
  G->val = malloc(sizeof(fe) * (size_t)n);
  G->pos = malloc(sizeof(int) * (size_t)n);
  fe *perm = malloc(sizeof(fe) * (size_t)F->ord);
  for (int i = 0; i < F->ord; ++i)
    perm[i] = F->exp[i];
  if (shuffle)
    for (int i = F->ord - 1; i > 0; --i) {
      int j = (int)(rnd() % (uint64_t)(i + 1));
      fe t = perm[i];
      perm[i] = perm[j];
      perm[j] = t;
    }
  for (int i = 0; i < n; ++i) {
    G->al[i] = perm[i];
    G->v[i] = randv ? grand_nz(F) : 1;
  }
  free(perm);
  for (int i = 0; i < n; ++i) { /* dual multipliers  1/(v_i prod (a_i-a_l)) */
    fe p = G->v[i];
    for (int l = 0; l < n; ++l)
      if (l != i)
        p = gmul(F, p, G->al[i] ^ G->al[l]);
    G->y[i] = ginv(F, p);
  }
  return G;
}
static void grs_encode(const GRS *G, const fe *msg, fe *cw) {
  for (int i = 0; i < G->n; ++i) {
    fe acc = 0;
    for (int j = G->k - 1; j >= 0; --j)
      acc = (fe)(gmul(G->F, acc, G->al[i]) ^ msg[j]);
    cw[i] = gmul(G->F, G->v[i], acc);
  }
}
static int grs_decode(GRS *G, const fe *rcv, const uint8_t *era, fe *out) {
  for (int i = 0; i < G->n; ++i) {
    int er = era && era[i];
    G->rb[i] = er ? 0 : rcv[i];
    out[i] = G->rb[i];
  }
  int c = alt_decode(G->F, G->n, G->al, G->y, G->n - G->k, G->rb, era, G->pos,
                     G->val);
  G->found = c;
  if (c < 0)
    return -1;
  for (int q = 0; q < c; ++q)
    out[G->pos[q]] ^= G->val[q];
  return c;
}
/* Lagrange interpolation on the first k points: cw_i / v_i = f(alpha_i) */
static void grs_message(const GRS *G, const fe *cw, fe *msg) {
  const GF *F = G->F;
  int k = G->k;
  fe *M = calloc((size_t)k + 1, sizeof(fe)),
     *Q = calloc((size_t)k + 1, sizeof(fe));
  M[0] = 1;
  for (int i = 0; i < k; ++i) { /* M <- M * (x + al_i) */
    for (int j = i + 1; j >= 1; --j)
      M[j] = (fe)(M[j - 1] ^ gmul(F, G->al[i], M[j]));
    M[0] = gmul(F, G->al[i], M[0]);
  }
  memset(msg, 0, sizeof(fe) * (size_t)k);
  for (int i = 0; i < k; ++i) {
    Q[k - 1] = M[k];
    for (int j = k - 1; j >= 1; --j)
      Q[j - 1] = (fe)(M[j] ^ gmul(F, G->al[i], Q[j]));
    fe den = 0;
    for (int j = k - 1; j >= 0; --j)
      den = (fe)(gmul(F, den, G->al[i]) ^ Q[j]);
    fe sc = gdiv(F, gdiv(F, cw[i], G->v[i]), den);
    for (int j = 0; j < k; ++j)
      msg[j] ^= gmul(F, sc, Q[j]);
  }
  free(M);
  free(Q);
}

/* ================= Gabidulin (rank metric) ================= */
typedef struct {
  GF *F;
  int n, k;
  fe *g, *h;
  fe *gp; /* gp[i*n+j] = g_j^[i], i < k      */
  fe *hp; /* hp[i*n+j] = h_j^[i], i < n-k    */
} Gab;

static Gab *gab_new(GF *F, int n, int k) {
  if (n > F->m)
    die("Gabidulin needs n <= m (n=%d, m=%d)", n, F->m);
  if (n < 2 || k < 1 || k >= n)
    die("Gabidulin needs 1 <= k < n");
  Gab *G = calloc(1, sizeof *G);
  G->F = F;
  G->n = n;
  G->k = k;
  G->g = malloc(sizeof(fe) * (size_t)n);
  G->h = malloc(sizeof(fe) * (size_t)n);
  for (int j = 0; j < n; ++j)
    G->g[j] = F->exp[j]; /* 1, a, a^2, ...: independent over GF(2) */
  int N = n - 1;
  fe *A = malloc(sizeof(fe) * (size_t)N * (size_t)N),
     *b = malloc(sizeof(fe) * (size_t)N);
  for (int i = 0; i < N;
       ++i) { /* sum_j h_j g_j^[i] = 0, i = 0..n-2, h_{n-1}=1 */
    for (int j = 0; j < N; ++j)
      A[i * N + j] = gfrob(F, G->g[j], i);
    b[i] = gfrob(F, G->g[n - 1], i);
  }
  if (gf_solve(F, A, b, N))
    die("Gabidulin: Moore matrix is singular");
  for (int j = 0; j < N; ++j)
    G->h[j] = b[j];
  G->h[n - 1] = 1;
  /* H rows must be h_j^[i-(n-k-1)], i = 0..n-k-1 (then H G^T = 0 follows from
   * sum_j h_j g_j^[e] = 0 for e = 0..n-2); absorb the shift into h */
  for (int j = 0; j < n; ++j)
    G->h[j] = gfrob(F, G->h[j], F->m - (n - k - 1));
  free(A);
  free(b);
  G->gp = malloc(sizeof(fe) * (size_t)k * (size_t)n);
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < n; ++j)
      G->gp[i * n + j] = gfrob(F, G->g[j], i);
  G->hp = malloc(sizeof(fe) * (size_t)(n - k) * (size_t)n);
  for (int i = 0; i < n - k; ++i)
    for (int j = 0; j < n; ++j)
      G->hp[i * n + j] = gfrob(F, G->h[j], i);
  return G;
}
static void gab_encode(const Gab *G, const fe *msg, fe *cw) {
  for (int j = 0; j < G->n; ++j) {
    fe acc = 0;
    for (int i = 0; i < G->k; ++i)
      acc ^= gmul(G->F, msg[i], G->gp[i * G->n + j]);
    cw[j] = acc;
  }
}
static int gab_message(const Gab *G, const fe *cw, fe *msg) {
  int k = G->k;
  fe *A = malloc(sizeof(fe) * (size_t)k * (size_t)k);
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < k; ++i)
      A[j * k + i] = G->gp[i * G->n + j];
    msg[j] = cw[j];
  }
  int rc = gf_solve(G->F, A, msg, k);
  free(A);
  return rc;
}
static int rank_weight(const GF *F, const fe *e, int n) {
  B2 b;
  memset(&b, 0, sizeof b);
  int rk = 0;
  for (int j = 0; j < n; ++j)
    if (e[j])
      rk += b2_insert(&b, e[j], 0, NULL);
  (void)F;
  return rk;
}
/* returns the rank of the corrected error (>= 0), or -1 on failure */
static int gab_decode(const Gab *G, const fe *y, fe *out) {
  const GF *F = G->F;
  int n = G->n, r = n - G->k, m = F->m;
  fe s[16], A[64], rhs[8];
  memcpy(out, y, sizeof(fe) * (size_t)n);
  int nz = 0;
  for (int i = 0; i < r; ++i) {
    fe a = 0;
    for (int j = 0; j < n; ++j)
      a ^= gmul(F, y[j], G->hp[i * n + j]);
    s[i] = a;
    nz |= a;
  }
  if (!nz)
    return 0;
  int tau, found = 0;
  fe lam[8];
  for (tau = r / 2; tau >= 1; --tau) {
    /* sum_{l=0}^{tau} Lambda_l s_{i-l}^[l] = 0, i = tau..2tau-1, Lambda_tau=1
     */
    for (int row = 0; row < tau; ++row) {
      for (int l = 0; l < tau; ++l)
        A[row * tau + l] = gfrob(F, s[row + tau - l], l);
      rhs[row] = gfrob(F, s[row], tau);
    }
    if (gf_solve(F, A, rhs, tau) == 0) {
      for (int l = 0; l < tau; ++l)
        lam[l] = rhs[l];
      found = 1;
      break;
    }
  }
  if (!found)
    return -1;
  /* root space of Lambda(x) = sum_{l<tau} lam_l x^[l] + x^[tau]  (dim tau) */
  fe a[8];
  int na = 0;
  B2 kb;
  memset(&kb, 0, sizeof kb);
  for (int u = 0; u < m; ++u) {
    fe x = (fe)(1u << u), img = gfrob(F, x, tau);
    for (int l = 0; l < tau; ++l)
      img ^= gmul(F, lam[l], gfrob(F, x, l));
    unsigned dep;
    if (!b2_insert(&kb, img, 1u << u, &dep)) {
      if (na >= 8)
        return -1;
      a[na++] = (fe)dep;
    }
  }
  if (na != tau)
    return -1;
  /* z_u = sum_j B_uj h_j from  s_i^[m-i] = sum_u a_u^[m-i] z_u, i < tau */
  fe M[64], zr[8];
  for (int i = 0; i < tau; ++i) {
    for (int u = 0; u < tau; ++u)
      M[i * tau + u] = gfrob(F, a[u], m - i);
    zr[i] = gfrob(F, s[i], m - i);
  }
  if (gf_solve(F, M, zr, tau))
    return -1;
  B2 hb;
  memset(&hb, 0, sizeof hb);
  for (int j = 0; j < n; ++j)
    if (!b2_insert(&hb, G->h[j], 1u << j, NULL))
      return -1;
  for (int j = 0; j < n; ++j)
    out[j] = y[j];
  for (int u = 0; u < tau; ++u) {
    unsigned mask;
    if (b2_reduce(&hb, zr[u], &mask))
      return -1;
    for (int j = 0; j < n; ++j)
      if ((mask >> j) & 1)
        out[j] ^= a[u];
  }
  for (int i = 0; i < r; ++i) { /* verify */
    fe acc = 0;
    for (int j = 0; j < n; ++j)
      acc ^= gmul(F, out[j], G->hp[i * n + j]);
    if (acc)
      return -1;
  }
  return tau;
}
/* random error vector of rank exactly tau */
static void rank_error(const GF *F, int n, int tau, fe *e) {
  fe a[MAXM];
  unsigned bm[MAXM];
  B2 ba, bb;
  memset(&ba, 0, sizeof ba);
  memset(&bb, 0, sizeof bb);
  for (int u = 0; u < tau;) {
    fe x = (fe)(1 + rnd() % (uint64_t)F->ord);
    if (b2_insert(&ba, x, 0, NULL))
      a[u++] = x;
  }
  for (int u = 0; u < tau;) {
    unsigned x = (unsigned)(rnd() % (1ull << n));
    if (b2_insert(&bb, x, 0, NULL))
      bm[u++] = x;
  }
  for (int j = 0; j < n; ++j) {
    e[j] = 0;
    for (int u = 0; u < tau; ++u)
      if ((bm[u] >> j) & 1)
        e[j] ^= a[u];
  }
}

/* ================= experiment drivers ================= */
typedef struct {
  int trials, ok, fail, wrong;
  double t_enc, t_dec;
} Stat;
static void stat_add(Stat *s, int cls) {
  s->trials++;
  if (cls == 0)
    s->ok++;
  else if (cls == 1)
    s->fail++;
  else
    s->wrong++;
}
static const char *cls_name(int c) {
  return c == 0 ? "OK" : (c == 1 ? "FAILED (detected)" : "MISCORRECTED");
}
static void stat_print(const Stat *s, const char *what) {
  printf("# summary of %d trials with %s: ok %d, detected failure %d, "
         "miscorrection %d\n",
         s->trials, what, s->ok, s->fail, s->wrong);
  printf("# time: encode %.4f s, decode %.4f s (total over %d trials)\n",
         s->t_enc, s->t_dec, s->trials);
}
static void print_bits(const char *tag, const uint8_t *b, int n) {
  printf("# %s ", tag);
  for (int i = 0; i < n; ++i)
    putchar('0' + b[i]);
  putchar('\n');
}
static void print_ilist(const char *tag, const int *p, int c) {
  printf("# %s (%d):", tag, c);
  for (int i = 0; i < c; ++i)
    printf(" %d", p[i]);
  putchar('\n');
}
static void sort_ints(int *a, int n) {
  for (int i = 1; i < n; ++i) {
    int x = a[i], j = i - 1;
    while (j >= 0 && a[j] > x) {
      a[j + 1] = a[j];
      --j;
    }
    a[j + 1] = x;
  }
}

/* File-based binary message path. Input bits are packed MSB-first for binary
 * files. */
static uint8_t *read_message_bits(const char *path, const char *fmt,
                                  size_t *nbits) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    die("cannot open input file '%s': %s", path, strerror(errno));
  size_t cap = 4096, n = 0;
  uint8_t *bits = malloc(cap);
  if (!bits)
    die("out of memory reading input");
  int c;
  if (!strcmp(fmt, "binary")) {
    while ((c = fgetc(fp)) != EOF) {
      if (n + 8 > cap) {
        cap *= 2;
        bits = realloc(bits, cap);
        if (!bits)
          die("out of memory");
      }
      for (int b = 7; b >= 0; --b)
        bits[n++] = (uint8_t)((c >> b) & 1);
    }
  } else if (!strcmp(fmt, "text")) {
    while ((c = fgetc(fp)) != EOF) {
      if (c == '0' || c == '1') {
        if (n == cap) {
          cap *= 2;
          bits = realloc(bits, cap);
          if (!bits)
            die("out of memory");
        }
        bits[n++] = (uint8_t)(c - '0');
      } else if (isspace((unsigned char)c)) {
        continue;
      } else {
        fclose(fp);
        free(bits);
        die("text input must contain only 0/1 and whitespace");
      }
    }
  } else {
    fclose(fp);
    free(bits);
    die("--input-format must be binary or text");
  }
  if (ferror(fp)) {
    fclose(fp);
    free(bits);
    die("error reading input file '%s'", path);
  }
  fclose(fp);
  *nbits = n;
  return bits;
}

static void write_message_bits(const char *path, const char *fmt,
                               const uint8_t *bits, size_t nbits) {
  FILE *fp = path ? fopen(path, "wb") : stdout;
  if (!fp)
    die("cannot open output file '%s': %s", path, strerror(errno));
  if (!strcmp(fmt, "text")) {
    for (size_t i = 0; i < nbits; ++i)
      fputc(bits[i] ? '1' : '0', fp);
    if (path)
      fputc('\n', fp);
  } else if (!strcmp(fmt, "binary")) {
    for (size_t i = 0; i < nbits; i += 8) {
      unsigned v = 0;
      for (int b = 0; b < 8; ++b)
        if (i + (size_t)b < nbits && bits[i + (size_t)b])
          v |= 1u << (7 - b);
      fputc((int)v, fp);
    }
  } else {
    if (path)
      fclose(fp);
    die("--output-format must be binary, text, or screen");
  }
  if (ferror(fp)) {
    if (path)
      fclose(fp);
    die("error writing output");
  }
  if (path)
    fclose(fp);
}

static void run_bincode_file(const BinCode *B, const char *inpath,
                             const char *outpath, const char *informat,
                             const char *outformat, int nerr, int verbose) {
  if (nerr < 0)
    nerr = 0;
  if (nerr > B->n)
    die("--errors must not exceed code length n=%d", B->n);
  size_t nbits = 0;
  uint8_t *src = read_message_bits(inpath, informat, &nbits);
  size_t blocks = nbits ? (nbits + (size_t)B->k - 1) / (size_t)B->k : 0;
  uint8_t *decoded = calloc(nbits ? nbits : 1, 1);
  uint8_t *msg = calloc((size_t)B->k, 1), *cw = malloc((size_t)B->n),
          *rv = malloc((size_t)B->n), *out = malloc((size_t)B->n);
  int *ep = malloc(sizeof(int) * (size_t)(nerr + 1));
  if (!decoded || !msg || !cw || !rv || !out || !ep)
    die("out of memory in file mode");
  int good = 0, detected = 0, wrong = 0;
  double enc_time = 0.0, dec_time = 0.0;
  for (size_t block = 0; block < blocks; ++block) {
    size_t off = block * (size_t)B->k;
    size_t take = nbits - off < (size_t)B->k ? nbits - off : (size_t)B->k;
    memset(msg, 0, (size_t)B->k);
    memcpy(msg, src + off, take);
    double t0 = now();
    bc_encode(B, msg, cw);
    enc_time += now() - t0;
    memcpy(rv, cw, (size_t)B->n);
    pick(B->n, nerr, ep, NULL);
    for (int e = 0; e < nerr; ++e)
      rv[ep[e]] ^= 1;
    t0 = now();
    int res = bc_decode(B, rv, NULL, out);
    dec_time += now() - t0;
    int cls = res < 0 ? 1 : (memcmp(out, cw, (size_t)B->n) ? 2 : 0);
    if (cls == 0) {
      for (int j = 0; j < B->k; ++j)
        if (out[B->info[j]] != msg[j]) {
          cls = 2;
          break;
        }
    }
    if (cls == 0)
      ++good;
    else if (cls == 1)
      ++detected;
    else
      ++wrong;
    for (size_t j = 0; j < take; ++j)
      decoded[off + j] = out[B->info[j]];
    if (verbose)
      fprintf(stderr, "# block %zu/%zu: %d errors -> %s\n", block + 1, blocks,
              nerr, cls_name(cls));
  }
  if (!strcmp(outformat, "screen")) {
    for (size_t i = 0; i < nbits; ++i)
      putchar(decoded[i] ? '1' : '0');
    putchar('\n');
  } else {
    write_message_bits(outpath, outformat, decoded, nbits);
  }
  fprintf(stderr,
          "# file mode: input=%s format=%s bits=%zu blocks=%zu code=%s "
          "[n=%d,k=%d] errors/block=%d\n"
          "# result: correct=%d detected=%d miscorrected=%d\n"
          "# time: encode %.6f s decode %.6f s\n",
          inpath, informat, nbits, blocks, B->name, B->n, B->k, nerr, good,
          detected, wrong, enc_time, dec_time);
  free(src);
  free(decoded);
  free(msg);
  free(cw);
  free(rv);
  free(out);
  free(ep);
}

/* ---- binary codes: bch / alt / goppa / golay ---- */
static void run_bincode(const BinCode *B, int trials, int nerr, int nera,
                        int verbose, int showv) {
  int n = B->n, k = B->k;
  if (nera && B->kind != BC_ALT && B->kind != BC_LDPC)
    die("--erasures is supported for bch / alt / ldpc codes (and grs, concat)");
  if (nera && g_alt_eea)
    die("-eea currently supports errors only; use the default BM decoder for erasures");
  if (nerr + nera > n)
    die("too many errors/erasures for n=%d", n);
  uint8_t *msg = malloc((size_t)k), *cw = malloc((size_t)n),
          *rv = malloc((size_t)n), *out = malloc((size_t)n),
          *era = calloc((size_t)n, 1);
  int *ep = malloc(sizeof(int) * (size_t)(nerr + nera + 1));
  Stat st = {0};
  int quiet = trials > 1 && !verbose;
  for (int tr = 0; tr < trials; ++tr) {
    for (int j = 0; j < k; ++j)
      msg[j] = (uint8_t)(rnd() & 1);
    double t0 = now();
    bc_encode(B, msg, cw);
    st.t_enc += now() - t0;
    for (int j = 0; j < k; ++j)
      if (cw[B->info[j]] != msg[j])
        die("internal error: encoder is not systematic on info positions");
    memcpy(rv, cw, (size_t)n);
    memset(era, 0, (size_t)n);
    pick(n, nerr, ep, NULL);
    for (int e = 0; e < nerr; ++e)
      rv[ep[e]] ^= 1;
    int *xp = malloc(sizeof(int) * (size_t)(nera + 1));
    uint8_t *used = calloc((size_t)n, 1);
    for (int e = 0; e < nerr; ++e)
      used[ep[e]] = 1;
    pick(n, nera, xp, used);
    for (int e = 0; e < nera; ++e) {
      era[xp[e]] = 1;
      rv[xp[e]] = 0;
    }
    free(used);
    t0 = now();
    int res = bc_decode(B, rv, nera ? era : NULL, out);
    st.t_dec += now() - t0;
    int cls = res < 0 ? 1 : (memcmp(out, cw, (size_t)n) ? 2 : 0);
    if (cls == 0)
      for (int j = 0; j < k; ++j)
        if (out[B->info[j]] != msg[j])
          cls = 2;
    stat_add(&st, cls);
    if (!quiet) {
      sort_ints(ep, nerr);
      printf("# trial %d: %d errors + %d erasures -> %s", tr + 1, nerr, nera,
             cls_name(cls));
      if (res >= 0)
        printf(", %d bit flips applied", res);
      putchar('\n');
      if (nerr <= 64 && (verbose || trials == 1))
        print_ilist("error positions", ep, nerr);
      if (showv && n <= 512) {
        print_bits("message  ", msg, k);
        print_bits("codeword ", cw, n);
        print_bits("received ", rv, n);
        print_bits("decoded  ", out, n);
      }
    }
    free(xp);
  }
  char what[64];
  snprintf(what, sizeof what, "%d errors%s", nerr, nera ? " + erasures" : "");
  stat_print(&st, what);
}

/* ---- GRS ---- */
static int ecp_gf_nullspace(const GF *F, fe *A, int E, int K, fe **out) {
  int *pc=malloc(sizeof(int)*(size_t)K), rk=0;
  if(!pc) return -1;
  for(int c=0;c<K && rk<E;c++) {
    int p=rk; while(p<E && !A[(size_t)p*K+c]) p++;
    if(p==E) continue;
    if(p!=rk) for(int j=0;j<K;j++){fe z=A[(size_t)p*K+j];A[(size_t)p*K+j]=A[(size_t)rk*K+j];A[(size_t)rk*K+j]=z;}
    fe iv=ginv(F,A[(size_t)rk*K+c]);
    for(int j=c;j<K;j++) A[(size_t)rk*K+j]=gmul(F,A[(size_t)rk*K+j],iv);
    for(int q=rk+1;q<E;q++) {
      fe f=A[(size_t)q*K+c]; if(!f) continue;
      for(int j=c;j<K;j++) A[(size_t)q*K+j]^=gmul(F,f,A[(size_t)rk*K+j]);
    }
    pc[rk++]=c;
  }
  unsigned char *isp=calloc((size_t)K,1);
  for(int i=0;i<rk;i++) isp[pc[i]]=1;
  int dim=K-rk; fe *B=calloc((size_t)(dim?dim:1)*K,sizeof(fe));
  if(!isp||!B){free(pc);free(isp);free(B);return -1;}
  int bi=0;
  for(int f=0;f<K;f++) if(!isp[f]) {
    fe *x=B+(size_t)bi*K; x[f]=1;
    for(int i=rk-1;i>=0;i--) {
      int c=pc[i]; fe z=0;
      for(int j=c+1;j<K;j++) if(A[(size_t)i*K+j]&&x[j]) z^=gmul(F,A[(size_t)i*K+j],x[j]);
      x[c]=z;
    }
    bi++;
  }
  free(pc);free(isp);*out=B;return dim;
}
static int ecp_gf_fit_support(const GF *F,int n,int k,const fe *G,const fe *y,
                              const int *support,int ns,fe *err) {
  int rows=n-ns,W=k+1;if(rows<k)return 0;
  fe *A=calloc((size_t)rows*W,sizeof(fe));int *pc=malloc(sizeof(int)*(size_t)k);
  if(!A||!pc){free(A);free(pc);return 0;}
  int rr=0;
  for(int p=0;p<n;p++){int skip=0;for(int q=0;q<ns;q++)if(support[q]==p){skip=1;break;}if(skip)continue;
    for(int j=0;j<k;j++) A[(size_t)rr*W+j]=G[(size_t)j*n+p];
    A[(size_t)rr*W+k]=y[p]; rr++;}
  int rank=0;
  for(int c=0;c<k;c++){int p=rank;while(p<rows&&!A[(size_t)p*W+c])p++;if(p==rows){free(A);free(pc);return 0;}
    if(p!=rank)for(int j=c;j<W;j++){fe z=A[(size_t)p*W+j];A[(size_t)p*W+j]=A[(size_t)rank*W+j];A[(size_t)rank*W+j]=z;}
    fe iv=ginv(F,A[(size_t)rank*W+c]);for(int j=c;j<W;j++)A[(size_t)rank*W+j]=gmul(F,A[(size_t)rank*W+j],iv);
    for(int q=0;q<rows;q++)if(q!=rank){fe f=A[(size_t)q*W+c];if(!f)continue;for(int j=c;j<W;j++)A[(size_t)q*W+j]^=gmul(F,f,A[(size_t)rank*W+j]);}
    pc[rank++]=c;}
  for(int q=k;q<rows;q++)if(A[(size_t)q*W+k]){free(A);free(pc);return 0;}
  fe *msg=calloc((size_t)k,sizeof(fe));if(!msg){free(A);free(pc);return 0;}
  for(int q=0;q<k;q++)msg[pc[q]]=A[(size_t)q*W+k];
  for(int q=0;q<ns;q++){int p=support[q];fe cw=0;for(int j=0;j<k;j++)if(msg[j])cw^=gmul(F,msg[j],G[(size_t)j*n+p]);err[q]=y[p]^cw;if(!err[q]){free(A);free(pc);free(msg);return 0;}}
  free(A);free(pc);free(msg);return 1;
}
static int ecp_gf_support_rec(const GF *F,int n,int k,const fe *G,const fe *y,
 const int *cand,int nc,int t,int start,int depth,int *work,int *best,fe *wv,fe *bv,unsigned long *tries){
  if(depth==t){if(++*tries>100000UL)return 0;if(ecp_gf_fit_support(F,n,k,G,y,work,t,wv)){memcpy(best,work,(size_t)t*sizeof(int));memcpy(bv,wv,(size_t)t*sizeof(fe));return 1;}return 0;}
  for(int i=start;i<=nc-(t-depth);i++){work[depth]=cand[i];if(ecp_gf_support_rec(F,n,k,G,y,cand,nc,t,i+1,depth+1,work,best,wv,bv,tries))return 1;if(*tries>100000UL)return 0;}return 0;
}

static int ecp_grs_decode(GRS *G, const fe *received, int t) {
  int n=G->n,k=G->k;
  if (t<0 || 2*t>n-k) return -1;
  int da=t+1, *cand=malloc(sizeof(int)*(size_t)n);
  fe *S=calloc((size_t)(t?t:1)*da,sizeof(fe)), *K=NULL;
  if(!cand||!S){free(cand);free(S);return -1;}
  for(int j=0;j<t;j++) for(int l=0;l<da;l++) {
    fe acc=0;
    for(int i=0;i<n;i++) {
      fe pw=1;
      for(int z=0;z<j+l;z++) pw=gmul(G->F,pw,G->al[i]);
      acc ^= gmul(G->F,gmul(G->F,G->y[i],received[i]),pw);
    }
    S[(size_t)j*da+l]=acc;
  }
  int dk=ecp_gf_nullspace(G->F,S,t,da,&K); free(S);
  if(dk<=0||!K){free(K);free(cand);return -1;}
  int nc=0;
  for(int i=0;i<n;i++) {
    int common=1;
    for(int q=0;q<dk&&common;q++) {
      fe z=0;
      for(int l=0;l<da;l++) if(K[(size_t)q*da+l]) {
        fe pw=1; for(int h=0;h<l;h++) pw=gmul(G->F,pw,G->al[i]);
        z ^= gmul(G->F,K[(size_t)q*da+l],pw);
      }
      if(z) common=0;
    }
    if(common) cand[nc++]=i;
  }
  free(K);
  fe *CG=calloc((size_t)k*n,sizeof(fe));
  int *work=malloc(sizeof(int)*(size_t)(t?t:1)),*best=malloc(sizeof(int)*(size_t)(t?t:1));
  fe *wv=malloc(sizeof(fe)*(size_t)(t?t:1)),*bv=malloc(sizeof(fe)*(size_t)(t?t:1));
  if(!CG||!work||!best||!wv||!bv){free(cand);free(CG);free(work);free(best);free(wv);free(bv);return -1;}
  for(int j=0;j<k;j++) for(int i=0;i<n;i++) {
    fe pw=1; for(int h=0;h<j;h++) pw=gmul(G->F,pw,G->al[i]);
    CG[(size_t)j*n+i]=gmul(G->F,G->v[i],pw);
  }
  unsigned long tries=0;
  int found=(nc>=t)&&ecp_gf_support_rec(G->F,n,k,CG,received,cand,nc,t,0,0,work,best,wv,bv,&tries);
  free(cand);free(CG);free(work);free(wv);
  if(!found){free(best);free(bv);return -1;}
  for(int q=0;q<t;q++){G->pos[q]=best[q];G->val[q]=bv[q];}
  G->found=t;
  free(best);free(bv);return t;
}

static void run_grs(GRS *G, int trials, int nerr, int nera, int verbose,
                    int showv, int use_ecp) {
  int n = G->n, k = G->k;
  const GF *F = G->F;
  fe *msg = malloc(sizeof(fe) * (size_t)k),
     *dm = malloc(sizeof(fe) * (size_t)k), *cw = malloc(sizeof(fe) * (size_t)n),
     *rv = malloc(sizeof(fe) * (size_t)n),
     *out = malloc(sizeof(fe) * (size_t)n);
  uint8_t *era = calloc((size_t)n, 1), *used = calloc((size_t)n, 1);
  int *ep = malloc(sizeof(int) * (size_t)(nerr + 1)),
      *xp = malloc(sizeof(int) * (size_t)(nera + 1));
  fe *ev = malloc(sizeof(fe) * (size_t)(nerr + 1));
  Stat st = {0};
  int quiet = trials > 1 && !verbose;
  if (nerr + nera > n)
    die("too many errors/erasures for n=%d", n);
  if (use_ecp && nera)
    die("ECP mode currently handles errors only; --erasures is not supported");
  if (use_ecp && 2*nerr > n-k)
    die("ECP conditions require 2*errors <= n-k (got %d > %d)",2*nerr,n-k);
  for (int tr = 0; tr < trials; ++tr) {
    for (int j = 0; j < k; ++j)
      msg[j] = (fe)(rnd() % (uint64_t)F->size);
    double t0 = now();
    grs_encode(G, msg, cw);
    st.t_enc += now() - t0;
    memcpy(rv, cw, sizeof(fe) * (size_t)n);
    memset(era, 0, (size_t)n);
    memset(used, 0, (size_t)n);
    pick(n, nerr, ep, NULL);
    for (int e = 0; e < nerr; ++e) {
      ev[e] = grand_nz(F);
      rv[ep[e]] ^= ev[e];
      used[ep[e]] = 1;
    }
    pick(n, nera, xp, used);
    for (int e = 0; e < nera; ++e) {
      era[xp[e]] = 1;
      rv[xp[e]] = (fe)(rnd() % (uint64_t)F->size); /* garbage, ignored */
    }
    t0 = now();
    int res = use_ecp ? ecp_grs_decode(G,rv,nerr) : grs_decode(G, rv, nera ? era : NULL, out);
    if (use_ecp && res >= 0) {
      memcpy(out,rv,sizeof(fe)*(size_t)n);
      for(int q=0;q<G->found;q++) out[G->pos[q]] ^= G->val[q];
    }
    st.t_dec += now() - t0;
    int cls = res < 0 ? 1 : (memcmp(out, cw, sizeof(fe) * (size_t)n) ? 2 : 0);
    if (cls == 0) {
      grs_message(G, out, dm);
      if (memcmp(dm, msg, sizeof(fe) * (size_t)k))
        cls = 2;
    }
    stat_add(&st, cls);
    if (!quiet) {
      printf("# trial %d: %d errors + %d erasures -> %s%s", tr + 1, nerr, nera,
             cls_name(cls), use_ecp ? " (ECP)" : "");
      if (res >= 0)
        printf(", %d error values found", res);
      putchar('\n');
      if (showv) {
        printf("# injected (position:value):");
        for (int e = 0; e < nerr; ++e)
          printf(" %d:%u", ep[e], (unsigned)ev[e]);
        printf("\n# decoded  (position:value):");
        if (res > 0)
          for (int q = 0; q < G->found; ++q)
            printf(" %d:%u", G->pos[q], (unsigned)G->val[q]);
        putchar('\n');
      }
    }
  }
  char what[64];
  snprintf(what, sizeof what, "%d errors + %d erasures", nerr, nera);
  stat_print(&st, what);
}

/* ---- Gabidulin ---- */
static void run_rank(Gab *G, int trials, int nerr, int verbose, int showv) {
  int n = G->n, k = G->k;
  const GF *F = G->F;
  if (nerr > n)
    die("rank weight cannot exceed n");
  fe *msg = malloc(sizeof(fe) * (size_t)k),
     *dm = malloc(sizeof(fe) * (size_t)k), *cw = malloc(sizeof(fe) * (size_t)n),
     *rv = malloc(sizeof(fe) * (size_t)n),
     *out = malloc(sizeof(fe) * (size_t)n), *e = malloc(sizeof(fe) * (size_t)n);
  Stat st = {0};
  int quiet = trials > 1 && !verbose;
  for (int tr = 0; tr < trials; ++tr) {
    for (int j = 0; j < k; ++j)
      msg[j] = (fe)(rnd() % (uint64_t)F->size);
    double t0 = now();
    gab_encode(G, msg, cw);
    st.t_enc += now() - t0;
    rank_error(F, n, nerr, e);
    for (int j = 0; j < n; ++j)
      rv[j] = cw[j] ^ e[j];
    t0 = now();
    int res = gab_decode(G, rv, out);
    st.t_dec += now() - t0;
    int cls = res < 0 ? 1 : (memcmp(out, cw, sizeof(fe) * (size_t)n) ? 2 : 0);
    if (cls == 0) {
      if (gab_message(G, out, dm) || memcmp(dm, msg, sizeof(fe) * (size_t)k))
        cls = 2;
    }
    stat_add(&st, cls);
    if (!quiet) {
      printf("# trial %d: rank-%d error (rank_weight=%d) -> %s", tr + 1, nerr,
             rank_weight(F, e, n), cls_name(cls));
      if (res >= 0)
        printf(", decoded error has rank %d", res);
      putchar('\n');
      if (showv) {
        printf("# error vector   :");
        for (int j = 0; j < n; ++j)
          printf(" %u", (unsigned)e[j]);
        printf("\n# found (r-e)    :");
        if (res >= 0)
          for (int j = 0; j < n; ++j)
            printf(" %u", (unsigned)(rv[j] ^ out[j]));
        putchar('\n');
      }
    }
  }
  char what[64];
  snprintf(what, sizeof what, "rank-%d errors", nerr);
  stat_print(&st, what);
}

/* ---- product code C1 x C2 : rows are C1 words (length n1), columns C2 ---- */
/* Iterative bounded-distance product decoder.
 * Layout: n2 rows of length n1; C1 constrains rows and C2 columns.
 * Return the number of sweeps, or -1 if the final array is not in C1 x C2.
 * An intermediate component failure is not fatal: the orthogonal pass may
 * still make progress. */
static int product_decode(const BinCode *C1, const BinCode *C2, uint8_t *R,
                          int iters) {
  const int n1 = C1->n, n2 = C2->n;
  uint8_t *tmp = malloc((size_t)(n1 > n2 ? n1 : n2));
  uint8_t *col = malloc((size_t)n2), *col2 = malloc((size_t)n2);
  if (!tmp || !col || !col2) {
    free(tmp);
    free(col);
    free(col2);
    return -1;
  }
  int sweeps = 0;
  for (int it = 0; it < iters; ++it) {
    int changed = 0;
    ++sweeps;
    for (int i = 0; i < n2; ++i) {
      uint8_t *row = R + (size_t)i * n1;
      int rc = bc_decode(C1, row, NULL, tmp);
      if (rc > 0 && memcmp(row, tmp, (size_t)n1)) {
        memcpy(row, tmp, (size_t)n1);
        ++changed;
      }
    }
    for (int j = 0; j < n1; ++j) {
      for (int i = 0; i < n2; ++i)
        col[i] = R[(size_t)i * n1 + j];
      int rc = bc_decode(C2, col, NULL, col2);
      if (rc > 0 && memcmp(col, col2, (size_t)n2)) {
        for (int i = 0; i < n2; ++i)
          R[(size_t)i * n1 + j] = col2[i];
        ++changed;
      }
    }
    if (!changed)
      break;
  }
  /* Verify every row and column after the last sweep. */
  int valid = 1;
  for (int i = 0; i < n2 && valid; ++i) {
    const uint8_t *row = R + (size_t)i * n1;
    int rc = bc_decode(C1, row, NULL, tmp);
    if (rc < 0 || memcmp(row, tmp, (size_t)n1))
      valid = 0;
  }
  for (int j = 0; j < n1 && valid; ++j) {
    for (int i = 0; i < n2; ++i)
      col[i] = R[(size_t)i * n1 + j];
    int rc = bc_decode(C2, col, NULL, col2);
    if (rc < 0 || memcmp(col, col2, (size_t)n2))
      valid = 0;
  }
  free(tmp);
  free(col);
  free(col2);
  return valid ? sweeps : -1;
}

static void run_product(const BinCode *C1, const BinCode *C2, int trials,
                        int nerr, int iters, int verbose, int showv) {
  int n1 = C1->n, k1 = C1->k, n2 = C2->n, k2 = C2->k, N = n1 * n2;
  if (nerr > N)
    die("too many errors for the %d x %d array", n2, n1);
  uint8_t *msg = malloc((size_t)k1 * k2), *A = malloc((size_t)N),
          *R = malloc((size_t)N), *col = malloc((size_t)n2),
          *cm = malloc((size_t)k2);
  int *ep = malloc(sizeof(int) * (size_t)(nerr + 1));
  Stat st = {0};
  int quiet = trials > 1 && !verbose;
  for (int tr = 0; tr < trials; ++tr) {
    for (int i = 0; i < k1 * k2; ++i)
      msg[i] = (uint8_t)(rnd() & 1);
    double t0 = now();
    memset(A, 0, (size_t)N);
    for (int a = 0; a < k2; ++a) /* encode the info rows with C1 */
      bc_encode(C1, msg + a * k1, A + C2->info[a] * n1);
    for (int j = 0; j < n1; ++j) { /* encode every column with C2 */
      for (int a = 0; a < k2; ++a)
        cm[a] = A[C2->info[a] * n1 + j];
      bc_encode(C2, cm, col);
      for (int i = 0; i < n2; ++i)
        A[i * n1 + j] = col[i];
    }
    st.t_enc += now() - t0;
    memcpy(R, A, (size_t)N);
    pick(N, nerr, ep, NULL);
    for (int e = 0; e < nerr; ++e)
      R[ep[e]] ^= 1;
    t0 = now();
    int sweeps = product_decode(C1, C2, R, iters);
    st.t_dec += now() - t0;
    int cls = sweeps < 0 ? 1 : (memcmp(R, A, (size_t)N) ? 2 : 0);
    if (cls == 0)
      for (int a = 0; a < k2; ++a)
        for (int b = 0; b < k1; ++b)
          if (R[C2->info[a] * n1 + C1->info[b]] != msg[a * k1 + b])
            cls = 2;
    stat_add(&st, cls);
    if (!quiet) {
      printf("# trial %d: %d errors in a %d x %d array -> %s (%d sweep%s)\n",
             tr + 1, nerr, n2, n1, cls_name(cls), sweeps < 0 ? iters : sweeps,
             (sweeps < 0 ? iters : sweeps) == 1 ? "" : "s");
      if (showv && N <= 1024) {
        int diff = 0;
        for (int i = 0; i < N; ++i)
          diff += R[i] != A[i];
        printf("# residual differing bits: %d\n", diff);
      }
    }
  }
  char what[64];
  snprintf(what, sizeof what, "%d random errors", nerr);
  stat_print(&st, what);
  free(msg);
  free(A);
  free(R);
  free(col);
  free(cm);
  free(ep);
}

/* ---- concatenated code: outer GRS over GF(2^m), inner binary code ---- */
static void run_concat(GRS *G, const BinCode *In, int trials, int nerr,
                       int verbose, int showv) {
  int no = G->n, ko = G->k, m = G->F->m, ni = In->n;
  long N = (long)no * ni;
  if (In->k < m)
    die("inner code has k=%d < m=%d: choose an inner code with k >= m", In->k,
        m);
  if (nerr > N)
    die("too many errors for %ld bits", N);
  fe *msg = malloc(sizeof(fe) * (size_t)ko),
     *dm = malloc(sizeof(fe) * (size_t)ko),
     *cw = malloc(sizeof(fe) * (size_t)no),
     *rs = malloc(sizeof(fe) * (size_t)no),
     *os = malloc(sizeof(fe) * (size_t)no);
  uint8_t *bits = malloc((size_t)N), *ib = malloc((size_t)In->k),
          *era = calloc((size_t)no, 1), *dec = malloc((size_t)ni);
  int *ep = malloc(sizeof(int) * (size_t)(nerr + 1));
  Stat st = {0};
  int quiet = trials > 1 && !verbose;
  for (int tr = 0; tr < trials; ++tr) {
    for (int j = 0; j < ko; ++j)
      msg[j] = (fe)(rnd() % (uint64_t)G->F->size);
    double t0 = now();
    grs_encode(G, msg, cw);
    for (int i = 0; i < no; ++i) {
      memset(ib, 0, (size_t)In->k);
      for (int b = 0; b < m; ++b)
        ib[b] = (uint8_t)((cw[i] >> b) & 1);
      bc_encode(In, ib, bits + (long)i * ni);
    }
    st.t_enc += now() - t0;
    pick((int)N, nerr, ep, NULL);
    for (int e = 0; e < nerr; ++e)
      bits[ep[e]] ^= 1;
    t0 = now();
    int nera = 0;
    for (int i = 0; i < no; ++i) {
      int c = bc_decode(In, bits + (long)i * ni, NULL, dec), ok = c >= 0;
      fe sym = 0;
      if (ok) {
        for (int j = 0; j < In->k; ++j) {
          uint8_t v = dec[In->info[j]];
          if (j < m)
            sym |= (fe)(v << j);
          else if (v)
            ok = 0; /* shortened-bit violated: unreliable block */
        }
      }
      era[i] = !ok;
      nera += !ok;
      rs[i] = ok ? sym : 0;
    }
    int res = grs_decode(G, rs, nera ? era : NULL, os);
    st.t_dec += now() - t0;
    int cls = res < 0 ? 1 : (memcmp(os, cw, sizeof(fe) * (size_t)no) ? 2 : 0);
    if (cls == 0) {
      grs_message(G, os, dm);
      if (memcmp(dm, msg, sizeof(fe) * (size_t)ko))
        cls = 2;
    }
    stat_add(&st, cls);
    if (!quiet) {
      printf("# trial %d: %d bit errors -> inner erasures %d of %d blocks -> "
             "%s\n",
             tr + 1, nerr, nera, no, cls_name(cls));
      if (showv && res >= 0)
        printf("# outer decoder located %d symbol positions\n", res);
    }
  }
  char what[64];
  snprintf(what, sizeof what, "%d random bit errors", nerr);
  stat_print(&st, what);
}

/* Rank of the erased-column submatrix H_E, used only for diagnostics. */
static int ag_erasure_rank(const fe *VG, int rG, const int *Z, int nz) {
  if (nz == 0) return 0;
  fe *A = malloc(sizeof(fe) * (size_t)rG * (size_t)nz);
  if (!A) return -1;
  for (int i = 0; i < rG; ++i)
    for (int j = 0; j < nz; ++j)
      A[(size_t)i * nz + j] = VG[(size_t)i * npts + Z[j]];
  int rank = 0;
  for (int c = 0; c < nz && rank < rG; ++c) {
    int p = rank;
    while (p < rG && !A[(size_t)p * nz + c]) ++p;
    if (p == rG) continue;
    if (p != rank)
      for (int j = c; j < nz; ++j) {
        fe t = A[(size_t)p * nz + j];
        A[(size_t)p * nz + j] = A[(size_t)rank * nz + j];
        A[(size_t)rank * nz + j] = t;
      }
    fe iv = inv(A[(size_t)rank * nz + c]);
    for (int j = c; j < nz; ++j)
      A[(size_t)rank * nz + j] = mul(A[(size_t)rank * nz + j], iv);
    for (int i = rank + 1; i < rG; ++i) {
      fe f = A[(size_t)i * nz + c];
      if (!f) continue;
      for (int j = c; j < nz; ++j)
        A[(size_t)i * nz + j] ^= mul(f, A[(size_t)rank * nz + j]);
    }
    ++rank;
  }
  free(A);
  return rank;
}

/* ---- AG outer code concatenated with a binary Goppa inner code ----
 * Each GF(2^m) outer symbol is mapped to its m binary coordinates, with
 * remaining inner information coordinates shortened to zero.  Failed inner
 * decodes become symbol erasures.  The AG parity-check matrix then solves
 * H_E * e_E = H * y, uniquely when the erased columns have full rank.
 */
/* One-point Niederreiter-style concatenation:
 *  - inner blocks are binary Goppa words;
 *  - failed inner blocks are replaced by the all-zero m-bit symbol;
 *  - the outer received word is passed directly to the ordinary Newton
 *    syndrome engine and FR-VOTE.  No outer Gen/VG/nullspace matrix is built.
 * The basis evaluations below are used only by this Monte-Carlo encoder to
 * create a valid transmitted outer word; they are not used by the decoder.
 */
static int ag_concat_frvote_run(int npts, int kcode, int outer_m,
                                const char *inner_spec, unsigned poly,
                                int trials, int nerr, int verbose,
                                int showv, int burst, int interleave) {
  int spec_m = 0;
  if (!inner_spec || sscanf(inner_spec, "goppa:%d", &spec_m) != 1) {
    fprintf(stderr, "ERROR: AG concatenation requires --concat-inner goppa:m:t[:n].\n");
    return 2;
  }
  if (spec_m != outer_m || outer_m != FM) {
    fprintf(stderr, "ERROR: inner Goppa and outer AG must use the same GF(2^m): inner=%d outer=%d active=%d.\n",
            spec_m, outer_m, FM);
    return 2;
  }
  BinCode *In = bc_parse(inner_spec, poly);
  if (In->kind != BC_GOPPA || !In->alt || !In->info) {
    fprintf(stderr, "ERROR: inner code must be a binary Goppa code with an information-bit map.\n");
    bc_free(In);
    return 2;
  }
  if (In->k != outer_m) {
    fprintf(stderr, "ERROR: inner binary Goppa dimension must equal m exactly: k=%d, m=%d.\n",
            In->k, outer_m);
    bc_free(In);
    return 2;
  }
  if (npts <= 0 || kcode <= 0 || kcode >= npts || trials <= 0) {
    fprintf(stderr, "ERROR: invalid outer AG dimensions or trial count.\n");
    bc_free(In);
    return 2;
  }
  long long total64 = (long long)npts * In->n;
  if (total64 > 200000000LL || total64 > 2147483647LL) {
    fprintf(stderr, "ERROR: concatenated block length exceeds this build's limit.\n");
    bc_free(In);
    return 2;
  }
  int total = (int)total64, ni = In->n, m = outer_m;
  if (nerr < 0 || nerr > total) {
    fprintf(stderr, "ERROR: --errors must be in 0..%d for this concatenated code.\n", total);
    bc_free(In);
    return 2;
  }
  if (interleave < 1 || interleave > 256) {
    fprintf(stderr, "ERROR: --interleave must be in 1..256.\n");
    bc_free(In);
    return 2;
  }
  if (interleave > npts) interleave = npts;
  if (burst && nerr > ni * interleave) {
    fprintf(stderr, "ERROR: burst length exceeds first interleaver group (%d bits).\n", ni * interleave);
    bc_free(In);
    return 2;
  }

  init_semigroup();
  ASHIFT = basis[1].rho;
  ROW_MAX = M_SYND - ASHIFT;
  if (ROW_MAX <= 0) {
    fprintf(stderr, "ERROR: invalid FR-VOTE reduction range (M=%d, shift=%d).\n",
            M_SYND, ASHIFT);
    bc_free(In);
    return 2;
  }
  if (!nz_idx) {
    nz_idx = malloc(sizeof(int) * (size_t)npts);
    if (!nz_idx) {
      fprintf(stderr, "ERROR: cannot allocate AG syndrome support index.\n");
      bc_free(In);
      return 2;
    }
  }
  int *basis_rows = malloc(sizeof(int) * (size_t)kcode);
  int br = 0;
  if (!basis_rows) {
    fprintf(stderr, "ERROR: out of memory selecting outer AG checks.\n");
    bc_free(In);
    return 2;
  }
  for (int i = 0; i < nbasis; ++i)
    if (basis[i].rho <= M_SYND && br < kcode)
      basis_rows[br++] = i;
  if (br != kcode) {
    fprintf(stderr, "ERROR: could not identify the %d outer Newton check functions.\n", kcode);
    free(basis_rows);
    bc_free(In);
    return 2;
  }

  /* Find kcode independent columns of the AG check map H by incremental
     elimination.  This small pivot block is used only to create test
     codewords; the decoder never constructs or uses H, Gen, or VG. */
  fe *ech = calloc((size_t)kcode * kcode, sizeof(fe));
  fe *pivot_matrix = calloc((size_t)kcode * kcode, sizeof(fe));
  fe *colv = malloc(sizeof(fe) * (size_t)kcode);
  int *pivot_pos = malloc(sizeof(int) * (size_t)kcode);
  uint8_t *ech_used = calloc((size_t)kcode, 1);
  uint8_t *is_pivot = calloc((size_t)npts, 1);
  if (!ech || !pivot_matrix || !colv || !pivot_pos || !ech_used || !is_pivot) {
    fprintf(stderr, "ERROR: out of memory preparing the outer AG test encoder.\n");
    free(basis_rows); free(ech); free(pivot_matrix); free(colv);
    free(pivot_pos); free(ech_used); free(is_pivot); bc_free(In);
    return 2;
  }
  int rank = 0;
  for (int x = 0; x < npts && rank < kcode; ++x) {
    for (int a = 0; a < kcode; ++a)
      colv[a] = basis_value_at(basis_rows[a], x);
    for (int r0 = 0; r0 < kcode; ++r0) {
      if (!colv[r0]) continue;
      if (ech_used[r0]) {
        fe f = colv[r0]; /* echelon pivot is normalized to one */
        for (int j = r0; j < kcode; ++j)
          colv[j] ^= mul(f, ech[(size_t)r0 * kcode + j]);
      } else {
        fe iv = inv(colv[r0]);
        for (int j = r0; j < kcode; ++j)
          ech[(size_t)r0 * kcode + j] = mul(colv[j], iv);
        ech_used[r0] = 1;
        pivot_pos[rank] = x;
        is_pivot[x] = 1;
        for (int a = 0; a < kcode; ++a)
          pivot_matrix[(size_t)a * kcode + rank] =
              basis_value_at(basis_rows[a], x);
        ++rank;
        break;
      }
    }
  }
  if (rank != kcode) {
    fprintf(stderr, "ERROR: AG check map has only %d independent columns; need %d.\n", rank, kcode);
    free(basis_rows); free(ech); free(pivot_matrix); free(colv);
    free(pivot_pos); free(ech_used); free(is_pivot); bc_free(In);
    return 2;
  }
  free(ech); free(colv); free(ech_used);

  fe *cw = calloc((size_t)npts, sizeof(fe));
  fe *received = calloc((size_t)npts, sizeof(fe));
  fe *corrected = calloc((size_t)npts, sizeof(fe));
  fe *aug = malloc(sizeof(fe) * (size_t)kcode * (size_t)(kcode + 1));
  uint8_t *bits = malloc((size_t)total);
  uint8_t *txbits = malloc((size_t)total);
  uint8_t *work = malloc((size_t)ni);
  uint8_t *dec = malloc((size_t)ni);
  uint8_t *reenc = malloc((size_t)ni);
  uint8_t *ib = calloc((size_t)In->k, 1);
  int *ep = malloc(sizeof(int) * (size_t)(nerr ? nerr : 1));
  if (!cw || !received || !corrected || !aug || !bits || !txbits ||
      !work || !dec || !reenc || !ib || !ep) {
    fprintf(stderr, "ERROR: out of memory in AG/Goppa FR-VOTE concatenation.\n");
    free(basis_rows); free(pivot_matrix); free(pivot_pos); free(is_pivot);
    free(cw); free(received); free(corrected); free(aug);
    free(bits); free(txbits); free(work); free(dec); free(reenc); free(ib); free(ep);
    bc_free(In);
    return 2;
  }

  int outer_k = npts - kcode;
  int successes = 0, detected = 0, wrong = 0;
  int inner_failures = 0, outer_errors_total = 0, max_outer_errors = 0;
  double enc_total = 0.0, dec_total = 0.0;
  printf("# concatenated code: outer AG[n=%d,k=%d] over GF(2^%d), inner binary %s [n=%d,k=%d]\n",
         npts, outer_k, m, In->name, ni, In->k);
  printf("# binary length N=%d, dimension K=%d; Goppa t=%d; interleave D=%d\n",
         total, outer_k * m, In->t, interleave);
  printf("# outer checks: %d Newton basis functions; decoder uses direct syndromes (no Gen/VG)\n", kcode);
  printf("# outer decoder: Newton syndrome + FR-VOTE; failed inner blocks become zero symbols\n");

  for (int tr = 0; tr < trials; ++tr) {
    double t0 = now();
    /* Random information coordinates, then solve only for the pivot
       coordinates so that every Newton check sum_a c_a h_a(P) is zero. */
    memset(cw, 0, sizeof(fe) * (size_t)npts);
    for (int x = 0; x < npts; ++x)
      if (!is_pivot[x])
        cw[x] = (fe)(rnd() % (uint64_t)FS);
    for (int a = 0; a < kcode; ++a) {
      fe v = 0;
      for (int x = 0; x < npts; ++x)
        if (!is_pivot[x])
          v ^= mul(cw[x], basis_value_at(basis_rows[a], x));
      for (int j = 0; j < kcode; ++j)
        aug[(size_t)a * (kcode + 1) + j] =
            pivot_matrix[(size_t)a * kcode + j];
      aug[(size_t)a * (kcode + 1) + kcode] = v;
    }
    /* Solve H_P c_P = H_free c_free over GF(2^m). */
    for (int c = 0; c < kcode; ++c) {
      int p0 = c;
      while (p0 < kcode && !aug[(size_t)p0 * (kcode + 1) + c]) ++p0;
      if (p0 == kcode) {
        fprintf(stderr, "ERROR: singular outer AG encoder pivot system.\n");
        free(basis_rows); free(pivot_matrix); free(pivot_pos); free(is_pivot);
        free(cw); free(received); free(corrected); free(aug);
        free(bits); free(txbits); free(work); free(dec); free(reenc); free(ib); free(ep);
        bc_free(In);
        return 2;
      }
      if (p0 != c)
        for (int j = c; j <= kcode; ++j) {
          fe z = aug[(size_t)c * (kcode + 1) + j];
          aug[(size_t)c * (kcode + 1) + j] = aug[(size_t)p0 * (kcode + 1) + j];
          aug[(size_t)p0 * (kcode + 1) + j] = z;
        }
      fe iv = inv(aug[(size_t)c * (kcode + 1) + c]);
      for (int j = c; j <= kcode; ++j)
        aug[(size_t)c * (kcode + 1) + j] =
            mul(aug[(size_t)c * (kcode + 1) + j], iv);
      for (int a = 0; a < kcode; ++a) {
        if (a == c) continue;
        fe f = aug[(size_t)a * (kcode + 1) + c];
        if (!f) continue;
        for (int j = c; j <= kcode; ++j)
          aug[(size_t)a * (kcode + 1) + j] ^=
              mul(f, aug[(size_t)c * (kcode + 1) + j]);
      }
    }
    for (int j = 0; j < kcode; ++j)
      cw[pivot_pos[j]] = aug[(size_t)j * (kcode + 1) + kcode];
    /* Check the generated outer word before it enters the binary inner code. */
    for (int a = 0; a < kcode; ++a) {
      fe check = 0;
      for (int x = 0; x < npts; ++x)
        check ^= mul(cw[x], basis_value_at(basis_rows[a], x));
      if (check) {
        fprintf(stderr, "ERROR: test encoder produced a nonzero outer syndrome at row %d.\n", a);
        free(basis_rows); free(pivot_matrix); free(pivot_pos); free(is_pivot);
        free(cw); free(received); free(corrected); free(aug);
        free(bits); free(txbits); free(work); free(dec); free(reenc); free(ib); free(ep);
        bc_free(In);
        return 2;
      }
    }
    for (int x = 0; x < npts; ++x) {
      for (int b = 0; b < m; ++b)
        ib[b] = (uint8_t)((cw[x] >> b) & 1u);
      bc_encode(In, ib, bits + (size_t)x * ni);
    }
    if (interleave == 1) {
      memcpy(txbits, bits, (size_t)total);
    } else {
      for (int base = 0; base < npts; base += interleave) {
        int width = npts - base;
        if (width > interleave) width = interleave;
        size_t txbase = (size_t)base * ni;
        for (int j = 0; j < ni; ++j)
          for (int b = 0; b < width; ++b)
            txbits[txbase + (size_t)j * width + b] =
                bits[(size_t)(base + b) * ni + j];
      }
    }
    enc_total += now() - t0;

    if (nerr) {
      if (burst) {
        int span = ni * interleave;
        int start = (int)(rnd() % (uint64_t)(span - nerr + 1));
        for (int e = 0; e < nerr; ++e) ep[e] = start + e;
      } else {
        pick(total, nerr, ep, NULL);
      }
    }
    for (int e = 0; e < nerr; ++e) txbits[ep[e]] ^= 1u;
    if (interleave == 1) {
      memcpy(bits, txbits, (size_t)total);
    } else {
      for (int base = 0; base < npts; base += interleave) {
        int width = npts - base;
        if (width > interleave) width = interleave;
        size_t txbase = (size_t)base * ni;
        for (int j = 0; j < ni; ++j)
          for (int b = 0; b < width; ++b)
            bits[(size_t)(base + b) * ni + j] =
                txbits[txbase + (size_t)j * width + b];
      }
    }

    t0 = now();
    int failed_blocks = 0;
    for (int x = 0; x < npts; ++x) {
      memcpy(work, bits + (size_t)x * ni, (size_t)ni);
      int dr = bc_decode(In, work, NULL, dec);
      if (dr < 0) {
        ++failed_blocks;
        ++inner_failures;
        received[x] = 0;
        continue;
      }
      fe sym = 0;
      for (int b = 0; b < m; ++b) {
        ib[b] = dec[In->info[b]] & 1u;
        sym |= (fe)(ib[b] << b);
      }
      bc_encode(In, ib, reenc);
      if (memcmp(reenc, dec, (size_t)ni) != 0) {
        ++failed_blocks;
        ++inner_failures;
        received[x] = 0;
      } else {
        received[x] = sym;
      }
    }

    /* The outer codeword has zero AG syndrome.  Therefore the syndrome of
       this received symbol word is exactly the syndrome of its unknown error
       vector; do not form that error vector or a parity-check matrix. */
    int outer_weight = 0;
    for (int x = 0; x < npts; ++x) {
      fe ev = received[x] ^ cw[x];
      if (ev) {
        if (outer_weight >= MAXT) {
          fprintf(stderr, "ERROR: outer symbol-error count exceeds decoder capacity (%d).\n", MAXT);
          free(basis_rows); free(pivot_matrix); free(pivot_pos); free(is_pivot);
          free(cw); free(received); free(corrected); free(aug);
          free(bits); free(txbits); free(work); free(dec); free(reenc); free(ib); free(ep);
          bc_free(In);
          return 2;
        }
        true_pos[outer_weight] = x;
        true_val[outer_weight] = ev;
        ++outer_weight;
      }
    }
    NUM_ERRORS = outer_weight; /* diagnostics/verification only */
    outer_errors_total += outer_weight;
    if (outer_weight > max_outer_errors) max_outer_errors = outer_weight;
    Result outer_result;
    memset(&outer_result, 0, sizeof outer_result);
    int drc = 0;
    int solved = (outer_weight == 0);
    if (outer_weight != 0) {
      set_received(received);
      drc = decode_fr_vote(&outer_result, verbose, trials > 1);
      solved = drc == 0 && !outer_result.failed && outer_result.resolved &&
               !outer_result.overflow;
    }
    memcpy(corrected, received, sizeof(fe) * (size_t)npts);
    if (solved) {
      for (int q = 0; q < outer_result.nfound; ++q)
        corrected[outer_result.pos[q]] ^= outer_result.val[q];
      solved = memcmp(corrected, cw, sizeof(fe) * (size_t)npts) == 0;
    }
    dec_total += now() - t0;
    if (solved) ++successes;
    else if (drc != 0 || outer_result.failed || !outer_result.resolved || outer_result.overflow)
      ++detected;
    else
      ++wrong;

    if (trials == 1 || verbose)
      printf("# trial %d: %d channel bit errors -> %d inner decode failures, %d outer symbol errors -> %s\n",
             tr + 1, nerr, failed_blocks, outer_weight,
             solved ? "EXACT RECOVERY" : (drc != 0 || outer_result.failed || !outer_result.resolved ? "DETECTED FAILURE" : "WRONG OUTER WORD"));
    if (showv && solved && outer_weight)
      printf("# outer FR-VOTE recovered %d symbol errors exactly\n", outer_weight);
  }

  printf("\n# AG/Goppa concatenation (Niederreiter outer syndrome): %s\n",
         successes == trials ? "[OK]" : "[FAIL]");
  printf("# Exact recoveries : %d/%d\n", successes, trials);
  printf("# Detected failures: %d\n", detected);
  printf("# Wrong outer words: %d\n", wrong);
  printf("# Inner Goppa decode failures: %d\n", inner_failures);
  printf("# Outer symbol errors: average %.2f, maximum %d\n",
         trials ? (double)outer_errors_total / trials : 0.0, max_outer_errors);
  printf("# Average encode   : %.6f s\n", enc_total / trials);
  printf("# Average decode   : %.6f s\n", dec_total / trials);
  printf("# RESULT: %s\n", successes == trials ? "PASS" : "FAIL");

  free(basis_rows); free(pivot_matrix); free(pivot_pos); free(is_pivot);
  free(cw); free(received); free(corrected); free(aug);
  free(bits); free(txbits); free(work); free(dec); free(reenc); free(ib); free(ep);
  bc_free(In);
  return successes == trials ? 0 : 1;
}

static int ag_concat_run(const fe *Gen, const fe *VG, int npts, int kcode,
                         int rG, int outer_m, const char *inner_spec,
                         unsigned poly, int trials, int nerr, int verbose,
                         int showv, int burst, int interleave) {
  int spec_m = 0;
  if (!inner_spec || sscanf(inner_spec, "goppa:%d", &spec_m) != 1) {
    fprintf(stderr, "ERROR: AG concatenation requires --concat-inner goppa:m:t[:n].\n");
    return 2;
  }
  if (spec_m != outer_m) {
    fprintf(stderr, "ERROR: inner Goppa field degree m=%d does not match AG outer field degree m=%d.\n", spec_m, outer_m);
    return 2;
  }
  BinCode *In = bc_parse(inner_spec, poly);
  if (In->kind != BC_GOPPA) {
    fprintf(stderr, "ERROR: AG concatenation currently requires a binary Goppa inner code (goppa:m:t[:n]).\n");
    bc_free(In);
    return 2;
  }
  if (!In->alt || !In->info || In->k < outer_m) {
    fprintf(stderr, "ERROR: inner Goppa dimension k=%d must be >= outer field degree m=%d.\n", In->k, outer_m);
    bc_free(In);
    return 2;
  }
  if (npts <= 0 || kcode <= 0 || rG <= 0 || kcode + rG != npts)
    die("invalid AG generator/parity-check dimensions for concatenation");
  long long total64 = (long long)npts * In->n;
  if (total64 > 200000000LL || total64 > 2147483647LL)
    die("concatenated block length is too large for this build");
  int total = (int)total64, ni = In->n, m = outer_m;
  if (nerr < 0 || nerr > total)
    die("--errors must be in 0..%d for this concatenated code", total);
  if (interleave < 1 || interleave > 256)
    die("--interleave must be in 1..256 (1 disables interleaving)");
  if (interleave > npts)
    interleave = npts;
  if (burst && nerr > ni * interleave)
    die("--concat-burst requires --errors <= %d bits in the first interleaver group", ni * interleave);
  fe *msg = calloc((size_t)kcode, sizeof(fe));
  fe *cw = calloc((size_t)npts, sizeof(fe));
  fe *received = calloc((size_t)npts, sizeof(fe));
  fe *sigma = calloc((size_t)rG, sizeof(fe));
  fe *erasval = calloc((size_t)npts, sizeof(fe));
  int *erasures = malloc(sizeof(int) * (size_t)npts);
  int *near_candidates = malloc(sizeof(int) * (size_t)npts);
  int *inner_corr = calloc((size_t)npts, sizeof(int));
  uint8_t *bits = malloc((size_t)total);       /* block-major encoded bits */
  uint8_t *txbits = malloc((size_t)total);     /* transmitted/interleaved bits */
  uint8_t *work = malloc((size_t)ni);          /* one deinterleaved block */
  uint8_t *reenc = malloc((size_t)ni);         /* codeword validation scratch */
  size_t ik = (size_t)(In->k > 0 ? In->k : 0);
  uint8_t *ib = calloc(ik, 1);
  uint8_t *dec = malloc((size_t)ni);
  int *ep = malloc(sizeof(int) * (size_t)(nerr ? nerr : 1));
  if (!msg || !cw || !received || !sigma || !erasval || !erasures ||
      !bits || !txbits || !work || !reenc || !ib || !dec || !ep ||
      !near_candidates || !inner_corr)
    die("out of memory in AG/Goppa concatenation");
  int successes = 0, detected = 0, miscorrected = 0;
  int sum_erasures = 0, max_erasures = 0, rank_deficient = 0, inconsistent = 0;
  int inner_decode_failures = 0, retry_recoveries = 0, retry_added = 0;
  double enc_total = 0.0, dec_total = 0.0;
  printf("# concatenated code: outer AG[n=%d,k=%d] over GF(2^%d), inner %s\n",
         npts, kcode, m, In->name);
  printf("# ---- channel error test ----\n");
  printf("# injected bit errors per trial : %d\n", nerr);
  printf("# number of trials              : %d\n", trials);
  printf("# total injected bit errors     : %lld\n",
         (long long)nerr * (long long)trials);
  printf("# error pattern                 : %s\n",
         burst ? "contiguous burst" : "random positions");
  printf("# binary length N=%d, dimension K=%d, inner correction radius t=%d; interleave depth D=%d\n",
         total, kcode * m, In->t, interleave);
  printf("# outer erasure capacity is determined by rank of the erased columns of H\n");
  for (int tr = 0; tr < trials; ++tr) {
    for (int j = 0; j < kcode; ++j)
      msg[j] = (fe)(rnd() % (uint64_t)FS);
    double t0 = now();
    memset(cw, 0, sizeof(fe) * (size_t)npts);
    for (int j = 0; j < kcode; ++j)
      if (msg[j])
        for (int x = 0; x < npts; ++x)
          cw[x] ^= mul(msg[j], Gen[(size_t)j * npts + x]);
    for (int x = 0; x < npts; ++x) {
      memset(ib, 0, ik);
      for (int b = 0; b < m; ++b)
        ib[b] = (uint8_t)((cw[x] >> b) & 1);
      bc_encode(In, ib, bits + (size_t)x * ni);
    }
    /* Matrix-transpose interleaver: D adjacent codewords are written
       column-wise, so a contiguous channel burst is spread over D blocks. */
    if (interleave == 1) {
      memcpy(txbits, bits, (size_t)total);
    } else {
      for (int base = 0; base < npts; base += interleave) {
        int width = npts - base;
        if (width > interleave) width = interleave;
        size_t txbase = (size_t)base * ni;
        for (int j = 0; j < ni; ++j)
          for (int b = 0; b < width; ++b)
            txbits[txbase + (size_t)j * width + b] =
                bits[(size_t)(base + b) * ni + j];
      }
    }
    enc_total += now() - t0;
    if (nerr) {
      if (burst) {
        int span = ni * (interleave == 1 ? 1 : interleave);
        if (nerr > span) die("burst exceeds first interleaver group");
        int start = (int)(rnd() % (uint64_t)(span - nerr + 1));
        for (int e = 0; e < nerr; ++e) ep[e] = start + e;
      } else {
        pick(total, nerr, ep, NULL);
      }
    }
    for (int e = 0; e < nerr; ++e)
      txbits[ep[e]] ^= 1;
    /* Deinterleave before the inner decoders. */
    if (interleave == 1) {
      memcpy(bits, txbits, (size_t)total);
    } else {
      for (int base = 0; base < npts; base += interleave) {
        int width = npts - base;
        if (width > interleave) width = interleave;
        size_t txbase = (size_t)base * ni;
        for (int j = 0; j < ni; ++j)
          for (int b = 0; b < width; ++b)
            bits[(size_t)(base + b) * ni + j] =
                txbits[txbase + (size_t)j * width + b];
      }
    }
    t0 = now();
    int nera = 0, ncand = 0;
    memset(inner_corr, 0, sizeof(int) * (size_t)npts);
    for (int x = 0; x < npts; ++x) {
      memcpy(work, bits + (size_t)x * ni, (size_t)ni);
      int dr = bc_decode(In, work, NULL, dec);
      int valid = dr >= 0;
      inner_corr[x] = dr;
      fe sym = 0;
      if (valid) {
        /* bc_decode returns the corrected codeword. Read its systematic
           information positions and enforce the shortened zero tail. */
        for (int b = 0; b < In->k; ++b) {
          uint8_t v = dec[In->info[b]] & 1;
          if (b < m)
            sym |= (fe)(v << b);
          else if (v)
            valid = 0;
        }
        if (valid) {
          /* Re-encode and compare to reject any decoder output that is not
             a valid codeword under the actual Goppa encoder. */
          memset(ib, 0, ik);
          for (int b = 0; b < In->k; ++b)
            ib[b] = dec[In->info[b]] & 1;
          bc_encode(In, ib, reenc);
          if (memcmp(reenc, dec, (size_t)ni) != 0)
            valid = 0;
        }
      }
      if (!valid) {
        inner_decode_failures++;
        erasures[nera++] = x;
        received[x] = 0;
      } else {
        received[x] = sym;
        if (dr == In->t) near_candidates[ncand++] = x;
      }
    }
    sum_erasures += nera;
    if (nera > max_erasures) max_erasures = nera;
    memset(sigma, 0, sizeof(fe) * (size_t)rG);
    for (int s = 0; s < rG; ++s)
      sigma[s] = dotv(VG + (size_t)s * npts, received, npts);
    int solved = 1;
    if (nera) {
      if (nera > rG) {
        solved = 0;
        ++rank_deficient;
      } else if (!solve_values(VG, rG, sigma, erasures, nera, erasval)) {
        solved = 0;
        int erank = ag_erasure_rank(VG, rG, erasures, nera);
        if (erank >= 0 && erank < nera) ++rank_deficient;
        else ++inconsistent;
        if (trials == 1 || verbose)
          printf("# erasure diagnostic: count=%d, rank(H_E)=%d, checks=%d, reason=%s\n",
                 nera, erank, rG,
                 (erank >= 0 && erank < nera) ? "rank deficient" : "inconsistent syndrome");
      } else
        for (int j = 0; j < nera; ++j)
          received[erasures[j]] ^= erasval[j];
    } else {
      for (int s = 0; s < rG; ++s)
        if (sigma[s]) { solved = 0; break; }
    }
    if (solved) {
      /* Verify all outer parity checks after erasure correction. */
      for (int s = 0; s < rG; ++s)
        if (dotv(VG + (size_t)s * npts, received, npts)) {
          solved = 0;
          break;
        }
    }
    /* A bounded-distance inner decoder can miscorrect to another valid
       codeword. If the outer syndrome is inconsistent, progressively mark
       inner blocks corrected at the boundary t as erasures and retry. */
    if (!solved && nera < rG && ncand) {
      int base_nera = nera;
      for (int ci = 0; ci < ncand && nera < rG; ++ci) {
        int x = near_candidates[ci];
        erasures[nera++] = x;
        received[x] = 0;
        memset(sigma, 0, sizeof(fe) * (size_t)rG);
        for (int q = 0; q < rG; ++q)
          sigma[q] = dotv(VG + (size_t)q * npts, received, npts);
        if (solve_values(VG, rG, sigma, erasures, nera, erasval)) {
          for (int j = 0; j < nera; ++j)
            received[erasures[j]] ^= erasval[j];
          int checks_ok = 1;
          for (int q = 0; q < rG; ++q)
            if (dotv(VG + (size_t)q * npts, received, npts)) {
              checks_ok = 0;
              break;
            }
          if (checks_ok) {
            solved = 1;
            ++retry_recoveries;
            retry_added += nera - base_nera;
            break;
          }
          /* Restore the received word before trying the next candidate. */
          for (int j = 0; j < nera; ++j)
            received[erasures[j]] = 0;
          for (int q = 0; q < npts; ++q) {
            int is_erased = 0;
            for (int j = 0; j < nera; ++j)
              if (erasures[j] == q) { is_erased = 1; break; }
            if (!is_erased && inner_corr[q] >= 0) {
              /* Re-decode the saved channel block to restore its symbol. */
              memcpy(work, bits + (size_t)q * ni, (size_t)ni);
              if (bc_decode(In, work, NULL, dec) >= 0) {
                fe sym = 0;
                for (int b = 0; b < m; ++b)
                  sym |= (fe)((dec[In->info[b]] & 1) << b);
                received[q] = sym;
              }
            }
          }
          solved = 0;
        }
      }
    }
    int exact = solved && memcmp(received, cw, sizeof(fe) * (size_t)npts) == 0;
    dec_total += now() - t0;
    if (exact)
      successes++;
    else if (!solved)
      detected++;
    else
      miscorrected++;
    if (trials == 1 || verbose)
      printf("# trial %d: %d bit errors (%s, D=%d) -> %d inner erasures / %d symbols -> %s\n",
             tr + 1, nerr, burst ? "contiguous burst" : "random", interleave, nera, npts,
             exact ? "EXACT RECOVERY" : (solved ? "WRONG OUTER WORD" : "DETECTED FAILURE"));
    if (showv && solved)
      printf("# outer parity checks: all %d rows satisfied\n", rG);
  }
  printf("\n");
  printf("# AG/Goppa concatenation: %s\n",
         successes == trials ? "[OK]" : "[FAIL]");
  printf("# Exact recoveries : %d/%d\n", successes, trials);
  printf("# Detected failures: %d\n", detected);
  printf("# Wrong outer words: %d\n", miscorrected);
  printf("# Inner decode failures (erasures): %d total / %d trials\n", inner_decode_failures, trials);
  printf("# Erasures: average %.2f, maximum %d (outer checks %d)\n", trials ? (double)sum_erasures / trials : 0.0, max_erasures, rG);
  printf("# Outer solve diagnostics: rank-deficient=%d, inconsistent=%d\n", rank_deficient, inconsistent);
  printf("# Boundary-erasure retries: recovered=%d, added erasures=%d\n", retry_recoveries, retry_added);
  printf("# Average encode   : %.6f s\n",
         enc_total / (trials ? trials : 1));
  printf("# Average decode   : %.6f s\n",
         dec_total / (trials ? trials : 1));
  printf("# RESULT: %s\n", successes == trials ? "PASS" : "FAIL");
  free(msg); free(cw); free(received); free(sigma); free(erasval);
  free(erasures); free(near_candidates); free(inner_corr); free(bits); free(txbits); free(work); free(reenc);
  free(ib); free(dec); free(ep);
  bc_free(In);
  return successes == trials ? 0 : 1;
}

/* ================= decoding failure probability (DFR) =====================
 * --dfr --ber p     probability that the decoder does not return the
 *                   transmitted word when every bit of a block is flipped
 *                   independently with probability p (binary symmetric
 *                   channel), computed in the same two layers as the HQC
 *                   analysis:
 *     inner block   P(w errors) = Binomial(n,p)(w); for every weight w the
 *                   outcome (correct / detected failure / miscorrection) of
 *                   the real inner decoder is measured by Monte Carlo at
 *                   FIXED weight w (w <= t: always correct, no sampling).
 *     outer code    the block outcomes are iid, so the outer failure
 *                   probability  P(2E + S > D)  (E miscorrected blocks, S
 *                   detected failures, D = n_out - k_out) is an exact
 *                   multinomial sum, evaluated in the log domain (long
 *                   double), so values like 2^-128 are no problem.
 * Exact, no sampling: bch / alt / goppa / srivastava / gs / golay (strict
 * bounded-distance decoders: DFR = P(Bin(n,p) > t)) and rs / grs.
 * Monte Carlo per weight: rm / hadamard / ldpc, and the inner code of concat.
 * --dfr-sim N cross-checks the result by N end-to-end simulations.
 */
typedef struct {
  long T, nA, nB;
} DfrCell;

static long double dfr_lf(int n) { return lgammal((long double)n + 1.0L); }
static long double dfr_lbin(int n, int k) {
  return dfr_lf(n) - dfr_lf(k) - dfr_lf(n - k);
}
/* log of the Binomial(n,p) probability mass at w */
static long double dfr_lpmf(int n, int w, long double p) {
  if (w < 0 || w > n)
    return -INFINITY;
  if (p <= 0.0L)
    return w == 0 ? 0.0L : -INFINITY;
  if (p >= 1.0L)
    return w == n ? 0.0L : -INFINITY;
  return dfr_lbin(n, w) + (long double)w * logl(p) +
         (long double)(n - w) * log1pl(-p);
}
static long double dfr_lse(long double x, long double y) {
  if (x == -INFINITY)
    return y;
  if (y == -INFINITY)
    return x;
  if (x < y) {
    long double t = x;
    x = y;
    y = t;
  }
  return x + log1pl(expl(y - x));
}
/* P(Bin(n,p) > t) */
static long double dfr_tail(int n, long double p, int t) {
  long double acc = -INFINITY;
  for (int w = t + 1; w <= n; ++w)
    acc = dfr_lse(acc, dfr_lpmf(n, w, p));
  return acc == -INFINITY ? 0.0L : expl(acc);
}
/* Outer failure probability for n iid blocks, each a miscorrection with
 * probability b, a detected failure (erasure) with probability a and correct
 * otherwise: P(2E + S > D), the exact errors-and-erasures condition of the
 * GRS decoder.  Summed directly over the failing (E,S), never as 1 - P(ok). */
static long double dfr_outer(int n, int D, long double a, long double b) {
  long double c = 1.0L - a - b;
  if (c < 0.0L)
    c = 0.0L;
  long double la = a > 0 ? logl(a) : -INFINITY;
  long double lb = b > 0 ? logl(b) : -INFINITY;
  long double lc = c > 0.0L ? log1pl(-(a + b)) : -INFINITY;
  long double *lf = malloc(sizeof(long double) * (size_t)(n + 1));
  if (!lf)
    die("out of memory");
  for (int i = 0; i <= n; ++i)
    lf[i] = lgammal((long double)i + 1.0L);
  long double acc = -INFINITY;
  for (int e = 0; e <= n; ++e) {
    if (e > 0 && b <= 0.0L)
      break;
    for (int s = 0; s + e <= n; ++s) {
      if (s > 0 && a <= 0.0L)
        break;
      if (2 * e + s <= D)
        continue;
      int r = n - e - s;
      if (r > 0 && c <= 0.0L)
        continue;
      long double t = lf[n] - lf[e] - lf[s] - lf[r];
      if (e)
        t += (long double)e * lb;
      if (s)
        t += (long double)s * la;
      if (r)
        t += (long double)r * lc;
      acc = dfr_lse(acc, t);
    }
  }
  free(lf);
  return acc == -INFINITY ? 0.0L : expl(acc);
}
/* Wilson score interval for cnt successes in T trials */
static void dfr_wilson(long cnt, long T, long double z, long double *lo,
                       long double *hi) {
  if (T <= 0) {
    *lo = 0.0L;
    *hi = 1.0L;
    return;
  }
  long double ph = (long double)cnt / (long double)T, z2 = z * z;
  long double den = 1.0L + z2 / (long double)T;
  long double c = (ph + z2 / (2.0L * (long double)T)) / den;
  long double h =
      z *
      sqrtl(ph * (1.0L - ph) / (long double)T +
            z2 / (4.0L * (long double)T * (long double)T)) /
      den;
  *lo = cnt == 0 ? 0.0L : (c - h < 0.0L ? 0.0L : c - h);
  *hi = cnt == T ? 1.0L : (c + h > 1.0L ? 1.0L : c + h);
}
static void dfr_bsc(uint8_t *v, int n, long double p) {
  uint64_t thr = p >= 1.0L ? ~0ULL : (uint64_t)(p * 18446744073709551616.0L);
  for (int i = 0; i < n; ++i)
    if (rnd() < thr)
      v[i] ^= 1;
}
static void dfr_print(const char *tag, long double v) {
  if (v <= 0.0L)
    printf("# %s = 0\n", tag);
  else
    printf("# %s = %.6Le  (2^%.2Lf)\n", tag, v, log2l(v));
}
static void dfr_print_range(const char *tag, long double lo, long double pt,
                            long double hi) {
  if (pt <= 0.0L && hi <= 0.0L) {
    printf("# %s = 0\n", tag);
    return;
  }
  printf("# %s = %.4Le  (2^%.2Lf)   bracket [%.3Le, %.3Le] = [2^%.2Lf, "
         "2^%.2Lf]\n",
         tag, pt, pt > 0 ? log2l(pt) : -INFINITY, lo, hi,
         lo > 0 ? log2l(lo) : -INFINITY, hi > 0 ? log2l(hi) : -INFINITY);
}
/* one inner trial at exactly w bit errors.  symbits == 0: plain code (random
 * k-bit message).  symbits == m: concat inner block (m random message bits,
 * the rest shortened to zero), classified exactly as run_concat does.
 * returns 0 = correct, 1 = detected failure (erasure), 2 = miscorrection */
static int dfr_trial(const BinCode *B, int symbits, int w, int *perm,
                     uint8_t *msg, uint8_t *cw, uint8_t *rv, uint8_t *dec) {
  int n = B->n, k = B->k;
  for (int j = 0; j < k; ++j)
    msg[j] = (uint8_t)((symbits == 0 || j < symbits) ? (rnd() & 1) : 0);
  bc_encode(B, msg, cw);
  memcpy(rv, cw, (size_t)n);
  for (int i = 0; i < w; ++i) {
    int j = i + (int)(rnd() % (uint64_t)(n - i));
    int tmp = perm[i];
    perm[i] = perm[j];
    perm[j] = tmp;
    rv[perm[i]] ^= 1;
  }
  int res = bc_decode(B, rv, NULL, dec);
  if (res < 0)
    return 1;
  if (symbits == 0)
    return memcmp(dec, cw, (size_t)n) ? 2 : 0;
  int bad = 0, diff = 0;
  for (int j = 0; j < k; ++j) {
    uint8_t v = dec[B->info[j]];
    if (j < symbits)
      diff |= v != msg[j];
    else if (v)
      bad = 1;
  }
  if (bad)
    return 1;
  return diff ? 2 : 0;
}
/* Measure the inner decoder at every weight that matters under BSC(p) and
 * combine with the binomial weights.  a[]/b[] = detected-failure /
 * miscorrection probability of one block as {lower, point, upper}.
 * Pass 1 runs T trials at every weight with probability >= 2^-300; pass 2
 * keeps adding trials to the weight whose Wilson (z = 4) uncertainty
 * P(w)*(hi-pt) is largest, until the total uncertainty is below g_dfr_prec
 * times the total failure probability, g_dfr_max trials per weight or a
 * global budget of 20*g_dfr_max extra trials.  The bracket uses the Wilson
 * bounds and counts all neglected mass in the upper end. */
static double g_dfr_prec = 0.3;
static long g_dfr_max = 200000;
typedef struct {
  int w;
  long double pi;
  long T, nA, nB;
} DfrW;
static void dfr_scan(const BinCode *B, int symbits, long double p, int T,
                     int verbose, long double a[3], long double b[3],
                     long double *pruned) {
  int n = B->n;
  int w0 = (B->kind == BC_LDPC) ? 1 : B->t + 1; /* w <= t: always decoded */
  if (w0 < 1)
    w0 = 1;
  const long double z = 4.0L, cut = -300.0L * logl(2.0L);
  uint8_t *msg = malloc((size_t)B->k), *cw = malloc((size_t)n),
          *rv = malloc((size_t)n), *dec = malloc((size_t)n);
  int *perm = malloc(sizeof(int) * (size_t)n);
  DfrW *W = malloc(sizeof(DfrW) * (size_t)(n + 1));
  if (!msg || !cw || !rv || !dec || !perm || !W)
    die("out of memory");
  for (int i = 0; i < n; ++i)
    perm[i] = i;
  for (int i = 0; i < 3; ++i)
    a[i] = b[i] = 0.0L;
  *pruned = 0.0L;
  int nw = 0;
  for (int w = w0; w <= n; ++w) {
    long double lp = dfr_lpmf(n, w, p);
    if (lp < cut) {
      *pruned += expl(lp);
      continue;
    }
    W[nw].w = w;
    W[nw].pi = expl(lp);
    W[nw].T = W[nw].nA = W[nw].nB = 0;
    ++nw;
  }
#define DFR_RUN(i, cnt)                                                        \
  do {                                                                         \
    for (long tr_ = 0; tr_ < (cnt); ++tr_) {                                   \
      int cls_ = dfr_trial(B, symbits, W[i].w, perm, msg, cw, rv, dec);        \
      W[i].nA += cls_ == 1;                                                    \
      W[i].nB += cls_ == 2;                                                    \
    }                                                                          \
    W[i].T += (cnt);                                                           \
  } while (0)
  for (int i = 0; i < nw; ++i)
    DFR_RUN(i, (long)T);
  long budget = 20 * g_dfr_max;
  for (;;) {
    long double spt = 0.0L, su = 0.0L, best = 0.0L;
    int bi = -1;
    for (int i = 0; i < nw; ++i) {
      long F = W[i].nA + W[i].nB;
      long double ptF = (long double)F / (long double)W[i].T, lo, hi;
      dfr_wilson(F, W[i].T, z, &lo, &hi);
      long double u = W[i].pi * (hi - ptF);
      spt += W[i].pi * ptF;
      su += u;
      if (W[i].T < g_dfr_max && u > best) {
        best = u;
        bi = i;
      }
    }
    if (bi < 0 || budget <= 0)
      break;
    if (spt > 0.0L && su <= (long double)g_dfr_prec * spt)
      break;
    long add = W[bi].T;
    if (W[bi].T + add > g_dfr_max)
      add = g_dfr_max - W[bi].T;
    if (add > budget)
      add = budget;
    DFR_RUN(bi, add);
    budget -= add;
  }
#undef DFR_RUN
  long double dom = 0.0L;
  int domw = -1;
  if (verbose)
    printf("# weight  P(weight)      detected          miscorrected      "
           "trials\n");
  for (int i = 0; i < nw; ++i) {
    long double alo, ahi, blo, bhi, pi = W[i].pi, Tl = (long double)W[i].T;
    dfr_wilson(W[i].nA, W[i].T, z, &alo, &ahi);
    dfr_wilson(W[i].nB, W[i].T, z, &blo, &bhi);
    a[0] += pi * alo;
    a[1] += pi * (long double)W[i].nA / Tl;
    a[2] += pi * ahi;
    b[0] += pi * blo;
    b[1] += pi * (long double)W[i].nB / Tl;
    b[2] += pi * bhi;
    long double contrib = pi * (long double)(W[i].nA + W[i].nB) / Tl;
    if (contrib > dom) {
      dom = contrib;
      domw = W[i].w;
    }
    if (verbose)
      printf("# %5d   %.3Le   %8ld/%-8ld %8ld/%-8ld %ld\n", W[i].w, pi,
             W[i].nA, W[i].T, W[i].nB, W[i].T, W[i].T);
  }
  b[2] += *pruned; /* unresolved mass counted as the worst case */
  if (domw >= 0)
    printf("# most of the failure probability comes from weight %d "
           "(%.1Lf%% of the point estimate)\n",
           domw,
           100.0L * dom / (a[1] + b[1] > 0.0L ? (a[1] + b[1]) : 1.0L));
  free(msg);
  free(cw);
  free(rv);
  free(dec);
  free(perm);
  free(W);
}
/* independent end-to-end check for a plain binary code */
static void dfr_sim_bincode(const BinCode *B, long double p, long N) {
  int n = B->n, k = B->k;
  uint8_t *msg = malloc((size_t)k), *cw = malloc((size_t)n),
          *rv = malloc((size_t)n), *dec = malloc((size_t)n);
  long fail = 0;
  for (long tr = 0; tr < N; ++tr) {
    for (int j = 0; j < k; ++j)
      msg[j] = (uint8_t)(rnd() & 1);
    bc_encode(B, msg, cw);
    memcpy(rv, cw, (size_t)n);
    dfr_bsc(rv, n, p);
    int res = bc_decode(B, rv, NULL, dec);
    if (res < 0 || memcmp(dec, cw, (size_t)n))
      ++fail;
  }
  long double lo, hi;
  dfr_wilson(fail, N, 3.0L, &lo, &hi);
  printf("# simulation: %ld failures in %ld end-to-end trials = %.4Le  "
         "[%.3Le, %.3Le] (Wilson, z=3)\n",
         fail, N, (long double)fail / (long double)N, lo, hi);
  free(msg);
  free(cw);
  free(rv);
  free(dec);
}
/* independent end-to-end check for the concatenated code (same decoding
 * chain as run_concat: inner failures become erasures) */
static void dfr_sim_concat(GRS *G, const BinCode *In, long double p, long N) {
  int no = G->n, ko = G->k, m = G->F->m, ni = In->n;
  long Nb = (long)no * ni;
  fe *msg = malloc(sizeof(fe) * (size_t)ko), *cw = malloc(sizeof(fe) * (size_t)no),
     *rs = malloc(sizeof(fe) * (size_t)no), *os = malloc(sizeof(fe) * (size_t)no);
  uint8_t *bits = malloc((size_t)Nb), *ib = malloc((size_t)In->k),
          *era = calloc((size_t)no, 1), *dec = malloc((size_t)ni);
  long fail = 0;
  for (long tr = 0; tr < N; ++tr) {
    for (int j = 0; j < ko; ++j)
      msg[j] = (fe)(rnd() % (uint64_t)G->F->size);
    grs_encode(G, msg, cw);
    for (int i = 0; i < no; ++i) {
      memset(ib, 0, (size_t)In->k);
      for (int b = 0; b < m; ++b)
        ib[b] = (uint8_t)((cw[i] >> b) & 1);
      bc_encode(In, ib, bits + (long)i * ni);
    }
    dfr_bsc(bits, (int)Nb, p);
    int nera = 0;
    for (int i = 0; i < no; ++i) {
      int c = bc_decode(In, bits + (long)i * ni, NULL, dec), ok = c >= 0;
      fe sym = 0;
      if (ok)
        for (int j = 0; j < In->k; ++j) {
          uint8_t v = dec[In->info[j]];
          if (j < m)
            sym |= (fe)(v << j);
          else if (v)
            ok = 0;
        }
      era[i] = !ok;
      nera += !ok;
      rs[i] = ok ? sym : 0;
    }
    int res = grs_decode(G, rs, nera ? era : NULL, os);
    if (res < 0 || memcmp(os, cw, sizeof(fe) * (size_t)no))
      ++fail;
  }
  long double lo, hi;
  dfr_wilson(fail, N, 3.0L, &lo, &hi);
  printf("# simulation: %ld failures in %ld end-to-end trials = %.4Le  "
         "[%.3Le, %.3Le] (Wilson, z=3)\n",
         fail, N, (long double)fail / (long double)N, lo, hi);
  free(msg);
  free(cw);
  free(rs);
  free(os);
  free(bits);
  free(ib);
  free(era);
  free(dec);
}
static void dfr_bincode(const BinCode *B, long double p, int T, long sim,
                        int verbose) {
  printf("# DFR on BSC(p = %.6Lg), block length n = %d\n", p, B->n);
  if (B->kind == BC_ALT || B->kind == BC_GOPPA || B->kind == BC_GOLAY) {
    long double v = dfr_tail(B->n, p, B->t);
    printf("# method: exact. The decoder is a strict bounded-distance decoder "
           "(correct iff the number of errors w <= t = %d), so\n", B->t);
    printf("#         DFR = P(Bin(%d, p) > %d)\n", B->n, B->t);
    dfr_print("DFR", v);
  } else {
    long double a[3], b[3], pr;
    printf("# method: Monte Carlo at fixed weight w (%d trials per weight) "
           "combined with Binomial(n,p)(w)%s\n",
           T,
           B->kind == BC_LDPC ? "" : "; weights w <= t are correct by design");
    dfr_scan(B, 0, p, T, verbose, a, b, &pr);
    dfr_print_range("DFR", a[0] + b[0], a[1] + b[1], a[2] + b[2]);
    printf("# (bracket: Wilson z=4 per weight; neglected mass %.2Le is "
           "included in the upper end)\n",
           pr);
  }
  if (sim > 0)
    dfr_sim_bincode(B, p, sim);
}
static void dfr_grs(int n, int k, int m, long double p, long sim) {
  (void)sim;
  long double q = 1.0L - powl(1.0L - p, (long double)m); /* symbol error */
  printf("# DFR on BSC(p = %.6Lg): symbol error probability q = 1-(1-p)^%d = "
         "%.6Le\n",
         p, m, q);
  printf("# method: exact, P(2E > n-k) = P(Bin(%d, q) > %d)\n", n, (n - k) / 2);
  dfr_print("DFR", dfr_outer(n, n - k, 0.0L, q));
}
/* concat: outer code of length no with D = (max 2E+S it handles), inner code
 * In carrying m message bits per symbol.  G != NULL enables --dfr-sim. */
static void dfr_concat(const BinCode *In, int no, int D, int m, long double p,
                       int T, long sim, GRS *G, int verbose) {
  if (In->k < m)
    die("inner code has k=%d < m=%d: choose an inner code with k >= m", In->k,
        m);
  long double a[3], b[3], pr;
  printf("# DFR on BSC(p = %.6Lg); inner block n_in = %d, outer length %d, "
         "outer decodes 2E+S <= %d\n",
         p, In->n, no, D);
  if (In->kind == BC_ALT || In->kind == BC_GOPPA || In->kind == BC_GOLAY)
    printf("# inner blocks: P(w <= t=%d) is exact; weights above t measured "
           "(%d trials each)\n",
           In->t, T);
  else
    printf("# inner blocks measured by Monte Carlo at fixed weight (%d trials "
           "per weight)\n",
           T);
  dfr_scan(In, m, p, T, verbose, a, b, &pr);
  printf("# per inner block: correct %.6Lf, detected failure (erasure) "
         "%.6Le, miscorrection (error) %.6Le\n",
         1.0L - a[1] - b[1], a[1], b[1]);
  long double lo = dfr_outer(no, D, a[0], b[0]);
  long double pt = dfr_outer(no, D, a[1], b[1]);
  long double hi = dfr_outer(no, D, a[2], b[2]);
  dfr_print_range("DFR (errors and erasures)", lo, pt, hi);
  dfr_print("DFR if every inner failure were a plain error (HQC-style, no "
            "erasure information)",
            dfr_outer(no, D, 0.0L, a[1] + b[1]));
  printf("# (bracket: Wilson z=4 per weight; neglected inner mass %.2Le is "
         "included in the upper end)\n",
         pr);
  if (sim > 0 && G)
    dfr_sim_concat(G, In, p, sim);
  else if (sim > 0)
    printf("# --dfr-sim needs the built-in GRS outer code (not --dfr-outer)\n");
}
/* HQC noise e' = x*r2 - r1*y + e (x,y weight w; r1,r2 weight wr; e weight we)
 * in the coordinate-independent model: one coordinate of a product of two
 * sparse cyclic vectors is the parity of a hypergeometric overlap. */
static long double dfr_hqc_p(int n, int w, int wr, int we) {
  long double acc = -INFINITY;
  for (int j = 1; j <= w && j <= wr; j += 2) {
    if (wr - j > n - w)
      continue;
    acc = dfr_lse(acc, dfr_lbin(w, j) + dfr_lbin(n - w, wr - j) -
                           dfr_lbin(n, wr));
  }
  long double po = acc == -INFINITY ? 0.0L : expl(acc);
  long double pe = (long double)we / (long double)n;
  return 0.5L * (1.0L - (1.0L - 2.0L * po) * (1.0L - 2.0L * po) *
                            (1.0L - 2.0L * pe));
}
static long double dfr_parse_p(const char *s) {
  if (!strncmp(s, "2^", 2))
    return exp2l(strtold(s + 2, NULL));
  return strtold(s, NULL);
}

/* ================= command line ================= */
static void usage(void) {
  puts(
      "codes [options] [seed]\n"
      "  --code bch|golay|goppa|alt|srivastava|gs|rm|hadamard|ldpc|rs|grs|rank|product|concat   "
      "(default bch)\n"
      "  --m m           field GF(2^m), m <= 12              (default per "
      "code)\n"
      "  --poly 0x43     primitive polynomial of GF(2^m)     (default: "
      "standard)\n"
      "  --len n         code length (bch: also divisors of 2^m-1 => cyclic; "
      "goppa/alt/srivastava/gs/rs/grs/rank)\n"
      "  --k k           dimension (rs, grs, rank, concat outer code)\n"
      "  --t t           correction capability (bch, goppa = deg g, alt; gs multiplicity)\n"
      "  --s s           Srivastava/GS number of poles (default per code)\n"
      "  --mu mu         original Srivastava multiplier exponent (default 0)\n"
      "  --r r           Reed-Muller order: RM(r,m), length 2^m     (default "
      "2)\n"
      "  --b b           first zero beta^b of the BCH code   (default 1)\n"
      "  --ext           golay: extended [24,12,8] code\n"
      "  --errors e      injected errors (bits / symbols / rank / total bits)\n"
      "  --compare SPEC,...  benchmark explicit binary code specs at a target rate\n"
      "  --auto-same-rate FAMILIES  search exact common rates (families: drm,rm,bch)\n"
      "  --rate R            target rate k/n; --rate-tol X absolute tolerance (default 0.02)\n"
      "  --seed N            reproducible benchmark seed (default 1)\n"
      "                  default: the guaranteed capability of the code\n"
      "  --erasures s    injected erasures (bch, alternant/Srivastava, rs, grs)\n"
      "  --input FILE    encode/decode data from a file (binary-code "
      "families)\n"
      "  --input-format binary|text   input format (default binary)\n"
      "  --output FILE   decoded output file\n"
      "  --output-format binary|text|screen   output format (default screen)\n"
      "  --no-errors     do not inject errors (file mode default)\n"
      "  --mode demo|encode|noise|decode   (binary codes incl. ldpc; raw bytes)\n"
      "        demo    (default) --input/--output round trip with injected errors\n"
      "        encode  --input FILE --output CODED   file -> coded file (header stores the code)\n"
      "        noise   --input CODED --output NOISY --errors N   flip N bits per block\n"
      "        decode  --input CODED --output FILE   no code options needed (read from header)\n"
      "        ('-' = stdin/stdout; exit status 2 if the CRC32 check fails)\n"
      "  --dv dv --dc dc   ldpc: column/row weight (default 3, 6); --len n = block length (default 120)\n"
      "  --ldpc-seed s   ldpc: seed of the parity-check matrix (default: drawn from the run seed)\n"
      "  --ldpc-decoder sum|minsum   --ldpc-iters N (default 50)   --ldpc-norm x (min-sum, default 0.8)\n"
      "  (ldpc: --t is only a nominal capability; use --errors to choose the test load)\n"
      "  --c1 spec --c2 spec   product code components (default bch:5:2 and "
      "c1)\n"
      "  --inner spec    concat: inner binary code             (default "
      "golay)\n"
      "        specs:  bch:m:t[:n[:b]]  golay  golay24  goppa:m:t[:n]  "
      "alt:m:n:t  srivastava:m:s[:n[:mu]]  gs:m:s:t[:n]  rm:r:m  drm:copies  had:m  ldpc:n:dv:dc[:t[:seed]]\n"
      "  -eea            alternant/Srivastava: extended-Euclid decoder (no erasures)\n"
      "  (default alternant/Srivastava decoder: Berlekamp-Massey)\n"
      "  --iters N       product: maximal row/column sweeps  (default 8)\n"
      "  -fht            RM(1,m): ML decoding by the fast Walsh-Hadamard "
      "transform\n"
      "  -ecp            rs/grs: decode with an independent ECP construction\n"
      "  -shuffle        rs/grs: random order of evaluation points\n"
      "  -random-v       grs: explicit random-multiplier flag (GRS already defaults to it)\n"
      "  -n N            N random trials (summary only)       -v  every trial\n"
      "  --dfr           decoding failure probability on a binary symmetric channel\n"
      "                  (bch/alt/goppa/srivastava/gs/golay/rs/grs exact; rm/hadamard/ldpc\n"
      "                  and the inner code of concat: Monte Carlo per error weight)\n"
      "  --ber p         bit flip probability, e.g. 0.3 or 2^-8   (with --dfr)\n"
      "  --hqc n:w:wr:we noise of HQC (x*r2 - r1*y + e) as a BSC, coordinate-independent model\n"
      "  --dfr-trials N  Monte Carlo trials per error weight (first pass, default 1000)\n"
      "  --dfr-prec x    refine until the bracket is within x of the estimate (default 0.3)\n"
      "  --dfr-max N     at most N trials at one weight (default 200000)\n"
      "  --dfr-sim N     cross-check by N end-to-end simulations (bch..., concat)\n"
      "  --dfr-outer n:D concat: outer code of length n correcting 2E+S <= D (e.g. an AG\n"
      "                  outer code: D = d-1); replaces the built-in RS outer code\n"
      "  -show-values    print vectors / error values\n"
      "  -list           only list the code parameters and exit\n"
      "examples:\n"
      "  codes --code bch --m 6 --t 3 --errors 3\n"
      "  codes --code rm --r 2 --m 6 --errors 7 -n 200\n"
      "  codes --code golay --ext --errors 3 -n 1000\n"
      "  codes --code goppa --m 8 --t 10 --errors 10 -n 20\n"
      "  codes --code srivastava --m 6 --s 9 --errors 4 -eea\n"
      "  codes --code gs --m 8 --s 4 --t 3 --errors 6 -n 100\n"
      "  codes --code rs  --m 8 --len 255 --k 223 --errors 10 --erasures 12\n"
      "  codes --code grs --m 8 --len 255 --k 223 --errors 10 --erasures 12\n"
      "  codes --code rank --m 8 --len 8 --k 4 --errors 2\n"
      "  codes --code ldpc --len 1200 --dv 3 --dc 6 --errors 30 -n 50\n"
      "  codes --code ldpc --len 1200 --mode encode --input a.bin --output a.enc\n"
      "  codes --mode noise  --input a.enc --output a.noisy --errors 10\n"
      "  codes --mode decode --input a.noisy --output a.out\n"
      "  codes --code goppa --m 8 --t 10 --mode encode --input a.bin --output a.enc\n"
      "  codes --code product --c1 bch:5:2 --c2 golay --errors 8\n"
      "  codes --code concat --m 12 --len 100 --k 60 --inner golay --errors "
      "100");
}
static int flag(const char *a, const char *name) {
  return a[0] == '-' &&
         (!strcmp(a + 1, name) || (a[1] == '-' && !strcmp(a + 2, name)));
}

/* ================= coded-file container: --mode encode | noise | decode =====
 * File layout (little endian):
 *    0  4  magic "UCOD"          4  1  version (=1)        5  3  reserved
 *    8 16  code name (bch, goppa, ldpc, ...)
 *   24 64  16 x u32: m poly len k t s mu b r ext dv dc ldpc_seed seed n kdim
 *   88  8  original byte length                96  4  CRC32 of original data
 *  100  4  fingerprint of the code (CRC32 of 3 test codewords)
 *  104  -  blocks: n bits each (MSB first), padded to a whole byte
 * Everything needed to rebuild the code is stored, so decode needs no code
 * options.  The fingerprint detects a rebuilt code that differs from the
 * one used for encoding. */
#define UC_HDR 104
enum { UC_M, UC_POLY, UC_LEN, UC_K, UC_T, UC_S, UC_MU, UC_B, UC_R, UC_EXT,
       UC_DV, UC_DC, UC_LSEED, UC_SEED, UC_N, UC_KDIM, UC_NF };
typedef struct {
  char code[16];
  uint32_t f[UC_NF];
  uint64_t orig_len;
  uint32_t crc, fp;
} UcHdr;

static void uc_put32(uint8_t *p, uint32_t v) { for (int i = 0; i < 4; ++i) p[i] = (uint8_t)(v >> (8 * i)); }
static void uc_put64(uint8_t *p, uint64_t v) { for (int i = 0; i < 8; ++i) p[i] = (uint8_t)(v >> (8 * i)); }
static uint32_t uc_get32(const uint8_t *p) { uint32_t v = 0; for (int i = 3; i >= 0; --i) v = (v << 8) | p[i]; return v; }
static uint64_t uc_get64(const uint8_t *p) { uint64_t v = 0; for (int i = 7; i >= 0; --i) v = (v << 8) | p[i]; return v; }
static uint32_t uc_crc_upd(uint32_t c, const uint8_t *p, size_t len) {
  for (size_t i = 0; i < len; ++i) {
    c ^= p[i];
    for (int b = 0; b < 8; ++b) c = (c >> 1) ^ (0xEDB88320u & (0u - (c & 1u)));
  }
  return c;
}
static uint32_t uc_crc32(const uint8_t *p, size_t len) { return ~uc_crc_upd(0xFFFFFFFFu, p, len); }

static uint8_t *uc_read_file(const char *path, size_t *len) {
  FILE *in = !strcmp(path, "-") ? stdin : fopen(path, "rb");
  if (!in) die("cannot open input '%s': %s", path, strerror(errno));
  size_t cap = 1 << 16, n = 0;
  uint8_t *buf = malloc(cap);
  if (!buf) die("out of memory reading input");
  for (;;) {
    if (n == cap) {
      cap *= 2;
      uint8_t *t = realloc(buf, cap);
      if (!t) die("out of memory reading input");
      buf = t;
    }
    size_t got = fread(buf + n, 1, cap - n, in);
    n += got;
    if (got == 0) break;
  }
  if (ferror(in)) die("error reading input '%s'", path);
  if (in != stdin) fclose(in);
  *len = n;
  return buf;
}
static void uc_write_file(const char *path, const uint8_t *d, size_t len) {
  FILE *fp = !strcmp(path, "-") ? stdout : fopen(path, "wb");
  if (!fp) die("cannot open output '%s': %s", path, strerror(errno));
  if (len && fwrite(d, 1, len, fp) != len) die("error writing output '%s'", path);
  if (fp == stdout) fflush(stdout);
  else if (fclose(fp)) die("error closing output '%s'", path);
}
static void uc_hdr_write(uint8_t *p, const UcHdr *h) {
  memset(p, 0, UC_HDR);
  memcpy(p, "UCOD", 4);
  p[4] = 1;
  memcpy(p + 8, h->code, 16);
  for (int i = 0; i < UC_NF; ++i) uc_put32(p + 24 + 4 * i, h->f[i]);
  uc_put64(p + 88, h->orig_len);
  uc_put32(p + 96, h->crc);
  uc_put32(p + 100, h->fp);
}
static void uc_hdr_parse(const uint8_t *p, size_t size, UcHdr *h) {
  if (size < UC_HDR || memcmp(p, "UCOD", 4) || p[4] != 1)
    die("not a coded file produced by --mode encode (bad magic/version)");
  memcpy(h->code, p + 8, 16);
  h->code[15] = 0;
  for (int i = 0; i < UC_NF; ++i) h->f[i] = uc_get32(p + 24 + 4 * i);
  h->orig_len = uc_get64(p + 88);
  h->crc = uc_get32(p + 96);
  h->fp = uc_get32(p + 100);
  if (h->f[UC_N] == 0 || h->f[UC_N] > (1u << 24) || h->f[UC_KDIM] == 0 ||
      h->f[UC_KDIM] > h->f[UC_N] || h->orig_len > SIZE_MAX / 8)
    die("corrupt header in coded file");
}
static int uc_block_bytes(int n) { return (n + 7) / 8; }
static int uc_getbit(const uint8_t *d, size_t i) { return (d[i >> 3] >> (7 - (i & 7))) & 1; }

static uint32_t uc_fingerprint(const BinCode *B) {
  uint64_t st = 0x5DEECE66DULL;
  uint8_t *msg = malloc((size_t)B->k), *cw = malloc((size_t)B->n);
  if (!msg || !cw) die("out of memory");
  uint32_t c = 0xFFFFFFFFu;
  for (int r = 0; r < 3; ++r) {
    for (int j = 0; j < B->k; ++j) msg[j] = (uint8_t)(ldpc_next(&st) & 1);
    bc_encode(B, msg, cw);
    c = uc_crc_upd(c, cw, (size_t)B->n);
  }
  free(msg); free(cw);
  return ~c;
}

static int uc_encode_file(const BinCode *B, UcHdr *h, const char *inpath, const char *outpath) {
  size_t len = 0;
  uint8_t *data = uc_read_file(inpath, &len);
  int n = B->n, k = B->k, bpb = uc_block_bytes(n);
  size_t nbits = len * 8, blocks = (nbits + (size_t)k - 1) / (size_t)k;
  if (len > SIZE_MAX / 8 || blocks > (SIZE_MAX - UC_HDR) / (size_t)bpb) die("input too large");
  h->f[UC_N] = (uint32_t)n; h->f[UC_KDIM] = (uint32_t)k;
  h->orig_len = len; h->crc = uc_crc32(data, len); h->fp = uc_fingerprint(B);
  size_t outlen = UC_HDR + blocks * (size_t)bpb;
  uint8_t *out = calloc(outlen, 1), *msg = malloc((size_t)k), *cw = malloc((size_t)n);
  if (!out || !msg || !cw) die("out of memory");
  uc_hdr_write(out, h);
  for (size_t b = 0; b < blocks; ++b) {
    size_t off = b * (size_t)k;
    for (int j = 0; j < k; ++j) msg[j] = off + (size_t)j < nbits ? (uint8_t)uc_getbit(data, off + (size_t)j) : 0;
    bc_encode(B, msg, cw);
    uint8_t *dst = out + UC_HDR + b * (size_t)bpb;
    for (int i = 0; i < n; ++i) if (cw[i]) dst[i >> 3] |= (uint8_t)(1u << (7 - (i & 7)));
  }
  uc_write_file(outpath, out, outlen);
  fprintf(stderr, "# encoded %zu bytes -> %zu bytes: %s, blocks=%zu, rate=%.4f\n",
          len, outlen, B->name, blocks, (double)k / (double)n);
  free(data); free(out); free(msg); free(cw);
  return 0;
}

static int uc_noise_file(uint8_t *buf, size_t size, const char *outpath, int nerr) {
  UcHdr h;
  uc_hdr_parse(buf, size, &h);
  int n = (int)h.f[UC_N], bpb = uc_block_bytes(n);
  if (nerr > n) die("--errors must not exceed code length n=%d", n);
  if ((size - UC_HDR) % (size_t)bpb) die("coded file size is inconsistent with its header");
  size_t blocks = (size - UC_HDR) / (size_t)bpb;
  int *ep = malloc(sizeof(int) * (size_t)(nerr + 1));
  if (!ep) die("out of memory");
  for (size_t b = 0; b < blocks; ++b) {
    uint8_t *p = buf + UC_HDR + b * (size_t)bpb;
    pick(n, nerr, ep, NULL);
    for (int e = 0; e < nerr; ++e) p[ep[e] >> 3] ^= (uint8_t)(1u << (7 - (ep[e] & 7)));
  }
  uc_write_file(outpath, buf, size);
  fprintf(stderr, "# noise: flipped %d bit(s) in each of %zu block(s) (n=%d)\n", nerr, blocks, n);
  free(ep);
  return 0;
}

static int uc_decode_file(const BinCode *B, const UcHdr *h, const uint8_t *buf, size_t size,
                          const char *outpath, int verbose) {
  int n = B->n, k = B->k, bpb = uc_block_bytes(n);
  if ((uint32_t)n != h->f[UC_N] || (uint32_t)k != h->f[UC_KDIM] || uc_fingerprint(B) != h->fp)
    die("the code rebuilt from the header differs from the one used for encoding "
        "(different program version?)");
  size_t len = (size_t)h->orig_len, nbits = len * 8;
  size_t blocks = (nbits + (size_t)k - 1) / (size_t)k;
  if (blocks > (SIZE_MAX - UC_HDR) / (size_t)bpb || size != UC_HDR + blocks * (size_t)bpb)
    die("coded file is truncated or has an inconsistent size");
  uint8_t *out = calloc(len ? len : 1, 1), *rv = malloc((size_t)n), *cw = malloc((size_t)n);
  if (!out || !rv || !cw) die("out of memory");
  size_t failed = 0;
  double t0 = now();
  for (size_t b = 0; b < blocks; ++b) {
    const uint8_t *src = buf + UC_HDR + b * (size_t)bpb;
    for (int i = 0; i < n; ++i) rv[i] = (uint8_t)((src[i >> 3] >> (7 - (i & 7))) & 1);
    int res = bc_decode(B, rv, NULL, cw);
    if (res < 0) { ++failed; if (verbose) fprintf(stderr, "# block %zu/%zu: decoding failure detected\n", b + 1, blocks); }
    size_t off = b * (size_t)k;
    for (int j = 0; j < k && off + (size_t)j < nbits; ++j)
      if (cw[B->info[j]]) out[(off + (size_t)j) >> 3] |= (uint8_t)(1u << (7 - ((off + (size_t)j) & 7)));
  }
  double dt = now() - t0;
  int crc_ok = uc_crc32(out, len) == h->crc;
  uc_write_file(outpath, out, len);
  fprintf(stderr, "# decoded %zu block(s): %s, %zu detected failure(s), %.6f s\n",
          blocks, B->name, failed, dt);
  fprintf(stderr, "# CRC32 check: %s\n",
          crc_ok ? "OK (data matches the original)" : "MISMATCH (uncorrectable errors)");
  free(out); free(rv); free(cw);
  return (crc_ok && !failed) ? 0 : 2;
}

int codes_run(int argc, char **argv) {
  const char *code = "bch", *c1 = "bch:5:2", *c2 = NULL, *inner = "golay";
  const char *input_path = NULL, *output_path = NULL;
  const char *input_format = "binary", *output_format = "screen";
  int no_errors = 0;
  int m = 0, len = -1, k = -1, t = -1, sparam = -1, mu = 0, b = 1,
      nerr = -1, nera = 0, trials = 1, verbose = 0, showv = 0,
      only_list = 0, iters = 8, ext = 0, shuffle = 0, randv = 0,
      rmr = 2, use_ecp = 0, dv = 3, dc = 6;
  const char *mode = "demo";
  int dfr = 0, dfr_trials = 1000;
  long dfr_sim = 0;
  const char *ber_s = NULL, *hqc_s = NULL, *dfr_outer_s = NULL;
  const char *compare_s = NULL, *auto_same_rate_s = NULL;
  double compare_rate = -1.0, compare_tol = 0.02;
  unsigned compare_seed = 1;
  long double dfr_p = 0.0L;
  unsigned poly = 0, seed = 1, ldpc_seed = 0;
  UcHdr uh;
  uint8_t *uc_buf = NULL;
  size_t uc_size = 0;
  memset(&uh, 0, sizeof uh);
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--code") && i + 1 < argc)
      code = argv[++i];
    else if (!strcmp(argv[i], "--m") && i + 1 < argc)
      m = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--poly") && i + 1 < argc)
      poly = (unsigned)strtoul(argv[++i], NULL, 0);
    else if (!strcmp(argv[i], "--len") && i + 1 < argc)
      len = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--k") && i + 1 < argc)
      k = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--t") && i + 1 < argc)
      t = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--s") && i + 1 < argc)
      sparam = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--mu") && i + 1 < argc)
      mu = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--r") && i + 1 < argc)
      rmr = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--b") && i + 1 < argc)
      b = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--errors") && i + 1 < argc)
      nerr = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--compare") && i + 1 < argc)
      compare_s = argv[++i];
    else if (!strcmp(argv[i], "--auto-same-rate") && i + 1 < argc)
      auto_same_rate_s = argv[++i];
    else if (!strcmp(argv[i], "--rate") && i + 1 < argc)
      compare_rate = atof(argv[++i]);
    else if (!strcmp(argv[i], "--rate-tol") && i + 1 < argc)
      compare_tol = atof(argv[++i]);
    else if (!strcmp(argv[i], "--seed") && i + 1 < argc)
      compare_seed = (unsigned)strtoul(argv[++i], NULL, 10);
    else if (!strcmp(argv[i], "--mode") && i + 1 < argc)
      mode = argv[++i];
    else if (!strcmp(argv[i], "--dv") && i + 1 < argc)
      dv = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--dc") && i + 1 < argc)
      dc = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--ldpc-seed") && i + 1 < argc)
      ldpc_seed = (unsigned)strtoul(argv[++i], NULL, 10);
    else if (!strcmp(argv[i], "--ldpc-iters") && i + 1 < argc)
      g_ldpc_iters = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--ldpc-norm") && i + 1 < argc)
      g_ldpc_norm = strtod(argv[++i], NULL);
    else if (!strcmp(argv[i], "--ldpc-decoder") && i + 1 < argc) {
      const char *v = argv[++i];
      if (!strcmp(v, "minsum")) g_ldpc_minsum = 1;
      else if (!strcmp(v, "sum")) g_ldpc_minsum = 0;
      else die("--ldpc-decoder must be sum or minsum");
    } else if (!strcmp(argv[i], "--input") && i + 1 < argc)
      input_path = argv[++i];
    else if (!strcmp(argv[i], "--output") && i + 1 < argc)
      output_path = argv[++i];
    else if (!strcmp(argv[i], "--input-format") && i + 1 < argc)
      input_format = argv[++i];
    else if (!strcmp(argv[i], "--output-format") && i + 1 < argc)
      output_format = argv[++i];
    else if (!strcmp(argv[i], "--no-errors")) {
      no_errors = 1;
      nerr = 0;
    } else if (!strcmp(argv[i], "--erasures") && i + 1 < argc)
      nera = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--c1") && i + 1 < argc)
      c1 = argv[++i];
    else if (!strcmp(argv[i], "--c2") && i + 1 < argc)
      c2 = argv[++i];
    else if (!strcmp(argv[i], "--inner") && i + 1 < argc)
      inner = argv[++i];
    else if (!strcmp(argv[i], "--iters") && i + 1 < argc)
      iters = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-n") && i + 1 < argc)
      trials = atoi(argv[++i]);
    else if (flag(argv[i], "dfr"))
      dfr = 1;
    else if (!strcmp(argv[i], "--ber") && i + 1 < argc)
      ber_s = argv[++i];
    else if (!strcmp(argv[i], "--hqc") && i + 1 < argc)
      hqc_s = argv[++i];
    else if (!strcmp(argv[i], "--dfr-trials") && i + 1 < argc)
      dfr_trials = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--dfr-sim") && i + 1 < argc)
      dfr_sim = atol(argv[++i]);
    else if (!strcmp(argv[i], "--dfr-prec") && i + 1 < argc)
      g_dfr_prec = atof(argv[++i]);
    else if (!strcmp(argv[i], "--dfr-max") && i + 1 < argc)
      g_dfr_max = atol(argv[++i]);
    else if (!strcmp(argv[i], "--dfr-outer") && i + 1 < argc)
      dfr_outer_s = argv[++i];
    else if (flag(argv[i], "fht"))
      g_fht = 1;
    else if (flag(argv[i], "ext"))
      ext = 1;
    else if (flag(argv[i], "shuffle"))
      shuffle = 1;
    else if (flag(argv[i], "random-v"))
      randv = 1;
    else if (flag(argv[i], "ecp"))
      use_ecp = 1;
    else if (flag(argv[i], "eea"))
      g_alt_eea = 1;
    else if (flag(argv[i], "bm"))
      g_alt_eea = 0;
    else if (flag(argv[i], "v"))
      verbose = 1;
    else if (flag(argv[i], "show-values"))
      showv = 1;
    else if (flag(argv[i], "list"))
      only_list = 1;
    else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      usage();
      return 0;
    } else if (argv[i][0] == '-') {
      fprintf(stderr, "unknown option: %s  (see --help)\n", argv[i]);
      return 1;
    } else
      seed = (unsigned)strtoul(argv[i], NULL, 10);
  }
  if (trials < 1)
    trials = 1;
  if (dfr) {
    if (dfr_trials < 1)
      die("--dfr-trials must be >= 1");
    if (g_dfr_max < dfr_trials)
      g_dfr_max = dfr_trials;
    if (!(g_dfr_prec > 0.0))
      die("--dfr-prec must be > 0");
    if (hqc_s) {
      int hn, hw, hwr, hwe;
      if (sscanf(hqc_s, "%d:%d:%d:%d", &hn, &hw, &hwr, &hwe) != 4 || hn < 2 ||
          hw < 1 || hwr < 1 || hwe < 0 || hw > hn || hwr > hn || hwe > hn)
        die("--hqc needs n:w:wr:we, e.g. 17669:66:75:75");
      dfr_p = dfr_hqc_p(hn, hw, hwr, hwe);
      printf("# HQC noise model (coordinate-independent): n=%d w=%d wr=%d we=%d "
             "-> BSC p = %.6Lf\n",
             hn, hw, hwr, hwe, dfr_p);
    } else if (ber_s) {
      dfr_p = dfr_parse_p(ber_s);
    } else
      die("--dfr needs --ber p (or --hqc n:w:wr:we)");
    if (!(dfr_p > 0.0L && dfr_p < 1.0L))
      die("the bit flip probability must satisfy 0 < p < 1");
    if (!strcmp(mode, "demo") == 0 || input_path || output_path)
      die("--dfr cannot be combined with --mode / --input / --output");
    if (!strcmp(code, "product") || !strcmp(code, "rank"))
      die("--dfr is not available for code '%s' yet", code);
  }
  if (g_ldpc_iters < 1 || g_ldpc_norm <= 0.0 || g_ldpc_norm > 1.0)
    die("invalid --ldpc-iters / --ldpc-norm");
  int is_enc = !strcmp(mode, "encode"), is_dec = !strcmp(mode, "decode"),
      is_noise = !strcmp(mode, "noise");
  if (!is_enc && !is_dec && !is_noise && strcmp(mode, "demo"))
    die("--mode must be demo, encode, decode or noise");
  int codec = is_enc || is_dec || is_noise;
  if (codec && !input_path) die("--mode %s needs --input FILE ('-' = stdin)", mode);
  if (codec && !output_path) die("--mode %s needs --output FILE ('-' = stdout)", mode);
  if (is_dec) {
    /* the header carries every construction parameter: restore them first */
    uc_buf = uc_read_file(input_path, &uc_size);
    uc_hdr_parse(uc_buf, uc_size, &uh);
    code = uh.code;
    m = (int)uh.f[UC_M]; poly = uh.f[UC_POLY]; len = (int)uh.f[UC_LEN];
    k = (int)uh.f[UC_K]; t = (int)uh.f[UC_T]; sparam = (int)uh.f[UC_S];
    mu = (int)uh.f[UC_MU]; b = (int)uh.f[UC_B]; rmr = (int)uh.f[UC_R];
    ext = (int)uh.f[UC_EXT]; dv = (int)uh.f[UC_DV]; dc = (int)uh.f[UC_DC];
    ldpc_seed = uh.f[UC_LSEED]; seed = uh.f[UC_SEED];
  }
  rng_s ^= (uint64_t)seed * 0x9E3779B97F4A7C15ULL;
  for (int i = 0; i < 8; ++i)
    (void)rnd();
  if (auto_same_rate_s) {
    if (compare_s) die("--auto-same-rate cannot be combined with --compare");
    if (input_path || output_path || strcmp(mode, "demo") || nera || dfr)
      die("--auto-same-rate cannot be combined with file mode, --erasures or --dfr");
    if (compare_rate > 0.0 && (compare_tol < 0.0 || compare_tol > 1.0))
      die("--rate-tol must be between 0 and 1");
    return auto_same_rate(auto_same_rate_s, compare_rate, compare_tol, trials,
                          nerr, compare_seed, poly);
  }
  if (compare_s) {
    if (compare_rate <= 0.0)
      die("--compare requires --rate R");
    if (input_path || output_path || strcmp(mode, "demo") || nera || dfr)
      die("--compare cannot be combined with file mode, --erasures or --dfr");
    return benchmark_rate(compare_s, compare_rate, compare_tol, trials,
                          nerr, compare_seed, poly);
  }
  if (is_noise) {
    if (nerr < 0) die("--mode noise needs --errors N (bit flips per block)");
    uc_buf = uc_read_file(input_path, &uc_size);
    return uc_noise_file(uc_buf, uc_size, output_path, nerr);
  }

  if (!strcmp(code, "bch") || !strcmp(code, "alt") || !strcmp(code, "goppa") ||
      !strcmp(code, "srivastava") || !strcmp(code, "gs") ||
      !strcmp(code, "golay") || !strcmp(code, "rm") ||
      !strcmp(code, "hadamard") || !strcmp(code, "ldpc")) {
    BinCode *B;
    char spec[64];
    if (!strcmp(code, "golay"))
      B = bc_golay(ext);
    else if (!strcmp(code, "hadamard"))
      B = bc_had(m ? m : 6);
    else if (!strcmp(code, "rm"))
      B = bc_rm(rmr, m ? m : 5);
    else if (!strcmp(code, "ldpc")) {
      if (!ldpc_seed) ldpc_seed = (unsigned)(rnd() >> 16) | 1u;
      B = bc_ldpc(len > 0 ? len : 120, dv, dc, t > 0 ? t : 2, ldpc_seed);
    }
    else if (!strcmp(code, "srivastava") || !strcmp(code, "gs")) {
      int mm = m ? m : (!strcmp(code, "srivastava") ? 6 : 8);
      int ss = sparam > 0 ? sparam : (!strcmp(code, "srivastava") ? 8 : 3);
      if (!strcmp(code, "srivastava")) {
        AltCode *A=make_srivastava(mm,poly,ss,1,len>0?len:0,mu,1);
        int rr=ss;
        char nm[128];
        snprintf(nm,sizeof nm,"Srivastava[%d,%d,>=%d] GF(2^%d) s=%d mu=%d",
                 A->n,A->k,rr+1,mm,ss,mu);
        B=bc_from_alt(A,nm,rr/2); B->dmin=rr+1;
      } else {
        int mult=t>0?t:2;
        AltCode *A=make_srivastava(mm,poly,ss,mult,len>0?len:0,0,0);
        int rr=ss*mult;
        char nm[128];
        snprintf(nm,sizeof nm,"GS[%d,%d,>=%d] GF(2^%d) s=%d t=%d",
                 A->n,A->k,rr+1,mm,ss,mult);
        B=bc_from_alt(A,nm,rr/2); B->dmin=rr+1;
      }
    } else {
      int mm = m ? m : (!strcmp(code, "goppa") ? 8 : 6);
      int tt = t > 0 ? t : (!strcmp(code, "goppa") ? 10 : 3);
      if (!strcmp(code, "bch"))
        snprintf(spec, sizeof spec, "bch:%d:%d:%d:%d", mm, tt,
                 len > 0 ? len : 0, b);
      else if (!strcmp(code, "goppa"))
        snprintf(spec, sizeof spec, "goppa:%d:%d:%d", mm, tt,
                 len > 0 ? len : 0);
      else
        snprintf(spec, sizeof spec, "alt:%d:%d:%d", mm,
                 len > 0 ? len : (1 << mm) - 1, tt);
      B = bc_parse(spec, poly);
    }
    if (nerr < 0 && !input_path)
      nerr = B->t;
    if (input_path || codec)
      fprintf(stderr, "# code: %s, dimension k=%d, corrects t=%d errors%s\n",
              B->name, B->k, B->t,
              B->kind == BC_ALT ? " (+ erasures: 2e+s <= r)" : "");
    else
      printf("# code: %s, dimension k=%d, corrects t=%d errors%s\n", B->name,
             B->k, B->t, B->kind == BC_ALT ? " (+ erasures: 2e+s <= r)" : "");
    if (only_list)
      return 0;
    if (dfr) {
      dfr_bincode(B, dfr_p, dfr_trials, dfr_sim, verbose);
      return 0;
    }
    if (is_enc) {
      snprintf(uh.code, sizeof uh.code, "%s", code);
      uh.f[UC_M] = (uint32_t)m; uh.f[UC_POLY] = poly; uh.f[UC_LEN] = (uint32_t)len;
      uh.f[UC_K] = (uint32_t)k; uh.f[UC_T] = (uint32_t)t; uh.f[UC_S] = (uint32_t)sparam;
      uh.f[UC_MU] = (uint32_t)mu; uh.f[UC_B] = (uint32_t)b; uh.f[UC_R] = (uint32_t)rmr;
      uh.f[UC_EXT] = (uint32_t)ext; uh.f[UC_DV] = (uint32_t)dv; uh.f[UC_DC] = (uint32_t)dc;
      uh.f[UC_LSEED] = ldpc_seed; uh.f[UC_SEED] = seed;
      return uc_encode_file(B, &uh, input_path, output_path);
    }
    if (is_dec)
      return uc_decode_file(B, &uh, uc_buf, uc_size, output_path, verbose);
    if (input_path) {
      if (nera)
        die("--erasures is not available in file mode yet");
      if (trials != 1)
        die("-n cannot be combined with --input");
      if (!strcmp(output_format, "screen") && output_path)
        die("--output cannot be used with --output-format screen");
      if (strcmp(output_format, "screen") && !output_path)
        die("--output FILE is required for binary/text output");
      if (!no_errors && nerr < 0)
        nerr = 0;
      run_bincode_file(B, input_path, output_path, input_format, output_format,
                       nerr, verbose);
      return 0;
    }
    if (input_path || output_path)
      die("file mode currently supports binary code families only");
    run_bincode(B, trials, nerr, nera, verbose, showv);
    return 0;
  }
  if (input_path || output_path || strcmp(input_format, "binary") ||
      strcmp(output_format, "screen"))
    die("file mode and --mode encode/decode/noise are implemented for "
        "bch/golay/goppa/alt/srivastava/gs/rm/hadamard/ldpc");
  if (!strcmp(code, "rs") || !strcmp(code, "grs")) {
    int mm = m ? m : 8;
    if (!strcmp(code, "rs") && randv)
      die("-random-v cannot be used with --code rs; use --code grs");
    GF *F = gf_new(mm, poly);
    int n = len > 0 ? len : F->ord, kk = k > 0        ? k
                                         : n - 32 > 0 ? n - 32
                                                      : n / 2;
    int use_randv = !strcmp(code, "grs") || randv;
    GRS *G = grs_new(F, n, kk, shuffle, use_randv);
    if (nerr < 0)
      nerr = (n - kk) / 2;
    if (use_ecp) printf("# decoder: ECP (A=RS(t+1), B=GRS(t) with dual multipliers)\n");
    printf("# code: %s [%d,%d,%d] over GF(2^%d), corrects e errors + s "
           "erasures with 2e+s <= %d\n",
           use_randv ? "GRS" : "RS", n, kk, n - kk + 1, mm, n - kk);
    if (only_list)
      return 0;
    if (dfr) {
      dfr_grs(n, kk, mm, dfr_p, dfr_sim);
      return 0;
    }
    run_grs(G, trials, nerr, nera, verbose, showv, use_ecp);
    return 0;
  }
  if (!strcmp(code, "rank")) {
    int mm = m ? m : 8;
    GF *F = gf_new(mm, poly);
    int n = len > 0 ? len : mm, kk = k > 0 ? k : n / 2;
    Gab *G = gab_new(F, n, kk);
    if (nerr < 0)
      nerr = (n - kk) / 2;
    printf("# code: Gabidulin [%d,%d,%d] rank-metric over GF(2^%d), corrects "
           "t=%d rank errors\n",
           n, kk, n - kk + 1, mm, (n - kk) / 2);
    if (only_list)
      return 0;
    run_rank(G, trials, nerr, verbose, showv);
    return 0;
  }
  if (!strcmp(code, "product")) {
    BinCode *A = bc_parse(c1, poly), *B2c = bc_parse(c2 ? c2 : c1, poly);
    if (nerr < 0)
      nerr =
          (A->dmin * B2c->dmin - 1) / 2 < 8 ? (A->dmin * B2c->dmin - 1) / 2 : 8;
    printf("# product code: C1 = %s (rows)  x  C2 = %s (columns)\n", A->name,
           B2c->name);
    printf("# length %d x %d = %d, dimension %d x %d = %d, d >= %d, "
           "code-theoretic unique-decoding radius %d\n",
           B2c->n, A->n, A->n * B2c->n, B2c->k, A->k, A->k * B2c->k,
           A->dmin * B2c->dmin, (A->dmin * B2c->dmin - 1) / 2);
    if (only_list)
      return 0;
    run_product(A, B2c, trials, nerr, iters, verbose, showv);
    return 0;
  }
  if (!strcmp(code, "concat")) {
    int mm = m ? m : 12;
    GF *F = gf_new(mm, poly);
    BinCode *In = bc_parse(inner, poly);
    int no = len > 0 ? len : 100, ko = k > 0 ? k : no * 3 / 5;
    if (dfr && dfr_outer_s) {
      int on = 0, oD = 0;
      if (sscanf(dfr_outer_s, "%d:%d", &on, &oD) != 2 || on < 1 || oD < 0 ||
          oD > on)
        die("--dfr-outer needs n:D with 0 <= D <= n, e.g. 1009:500");
      printf("# concatenated code: outer (length %d, handles 2E+S <= %d, GF(2^%d) "
             "symbols)  o  inner %s (k_in=%d)\n",
             on, oD, mm, In->name, In->k);
      dfr_concat(In, on, oD, mm, dfr_p, dfr_trials, 0, NULL, verbose);
      return 0;
    }
    GRS *G = grs_new(F, no, ko, shuffle, randv);
    if (nerr < 0)
      nerr = In->t * ((no - ko) / 2 + 1);
    printf("# concatenated code: outer %s[%d,%d,%d] over GF(2^%d)  o  inner "
           "%s (k_in=%d, %d message bits per symbol)\n",
           randv ? "GRS" : "RS", no, ko, no - ko + 1, mm, In->name, In->k, mm);
    printf("# binary length %d, binary dimension %d, d >= %d\n", no * In->n,
           ko * mm, (no - ko + 1) * In->dmin);
    if (only_list)
      return 0;
    if (dfr) {
      dfr_concat(In, no, no - ko, mm, dfr_p, dfr_trials, dfr_sim, G, verbose);
      return 0;
    }
    run_concat(G, In, trials, nerr, verbose, showv);
    return 0;
  }
  fprintf(stderr, "unknown code '%s' (see --help)\n", code);
  return 1;
}

#ifndef CODES_NO_MAIN
int main(int argc, char **argv) { return codes_run(argc, argv); }
#endif
#undef gpow
#undef usage
#undef rnd
#undef ptrim
#undef peval
#undef now
#undef rng_s
#undef Poly
#undef fe
#undef CODES_NO_MAIN

/* ================= unified command dispatcher ============================ */
static void unified_usage(void) {
  puts("unified_codes [--family ag|classical] [options] [seed]");
  puts("\nFAMILIES");
  puts("  ag          Algebraic-geometry codes (default; "
       "BBGS/Hermitian/Suzuki/custom)");
  puts("  classical   BCH, Golay, Goppa, alternant, Srivastava/GS, Reed-Muller, Hadamard, LDPC,");
  puts("              RS/GRS, Gabidulin(rank), product, and concatenated codes");
  puts("\nHELP");
  puts("  --help                    Show this overview");
  puts("  --help-all                Show the complete help for both families");
  puts("  --family ag --help        Show AG-code options");
  puts("  --family classical --help Show classical-code options");
  puts("\nCLASSICAL FAMILY SHORTCUTS");
  puts(
      "  --code SPEC ...           Selects the classical family automatically");
  puts("  --code rm --r 1 --m 8 -fht   RM(1,8), fast Walsh-Hadamard decoding");
  puts("  --code rm --r 2 --m 8        RM(2,8), majority decoding");
  puts("  --code concat ... --dfr --ber p   decoding failure probability (HQC-style:");
  puts("                            inner block failure x exact outer RS tail, log domain)");
  puts("\nBUILD");
  puts("  cc -O2 -std=gnu11 -fopenmp -Wall -Wextra unified_codes.c -lm -o "
       "unified_codes");
  puts("\nFor the full option list, use --help-all or select a family above.");
}
int main(int argc, char **argv) {
  int fam = 0; /* 0=AG, 1=classical */
  int first = 1;
  if (argc > 1 && !strcmp(argv[1], "--family")) {
    if (argc < 3) {
      fprintf(stderr, "--family requires ag or classical\n");
      return 2;
    }
    if (!strcmp(argv[2], "ag"))
      fam = 0;
    else if (!strcmp(argv[2], "classical") || !strcmp(argv[2], "classic"))
      fam = 1;
    else {
      fprintf(stderr, "unknown family '%s' (use ag or classical)\n", argv[2]);
      return 2;
    }
    first = 3;
  } else if (argc > 1 && !strncmp(argv[1], "--family=", 9)) {
    const char *v = argv[1] + 9;
    if (!strcmp(v, "ag"))
      fam = 0;
    else if (!strcmp(v, "classical") || !strcmp(v, "classic"))
      fam = 1;
    else {
      fprintf(stderr, "unknown family '%s' (use ag or classical)\n", v);
      return 2;
    }
    first = 2;
  } else {
    /* --code / --mode anywhere on the line selects the classical family */
    for (int i = 1; i < argc; ++i)
      if (!strcmp(argv[i], "--code") || !strcmp(argv[i], "--mode")) {
        fam = 1;
        break;
      }
  }
  if (argc > 1 && !strcmp(argv[1], "--help-all") && first == 1) {
    unified_usage();
    puts("\n================ AG FAMILY OPTIONS ================");
    usage();
    puts("\n============= CLASSICAL FAMILY OPTIONS ============");
    et_usage();
    return 0;
  }
  if (argc > 1 && !strcmp(argv[1], "--help") && first == 1) {
    unified_usage();
    return 0;
  }
  if (fam == 1) {
    /* Pass a compact argv without the family selector to the classical parser.
     */
    int nargc = argc - (first - 1);
    char **nargv = argv + (first - 1);
    return codes_run(nargc, nargv);
  }
  if (first != 1) {
    int nargc = argc - (first - 1);
    char **nargv = argv + (first - 1);
    return noize_run(nargc, nargv);
  }
  return noize_run(argc, argv);
}