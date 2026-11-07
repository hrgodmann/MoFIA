// MoFIA: per-person marginal likelihood of a moderated factor model (categorical
// and continuous items) by (adaptive) product Gauss-Hermite quadrature for any
// number of factors D, plus the posterior-weighted derivative sums from which R
// assembles the per-person scores.
//
// Inputs are person-level quantities precomputed in R (cheap, O(N P Q)):
//   NU  (N x P)      item shifts nu_ij (incl. intercepts of continuous items)
//   LAM (N x P*D)    loadings lambda_ijd, column j + P*d (0 if not loading)
//   MU  (N x D)      factor means
//   CL  (N x T)      Cholesky factor of the factor covariance matrix, row-wise
//                    lower triangle, T = D(D+1)/2, entry (d,e) at d(d+1)/2 + e
//   AM, ALC          adaptive quadrature: posterior means / Cholesky (same layout)
//   YC, RV, cont     continuous items: values, residual variances, indicator
// Thresholds are passed flattened (tau_flat, tau_off).
//
// Outputs: ll (per person), posterior moments (PM, PC), and for scores
//   PG  = E[g_ij], PGE = E[g_ij eta_d], DLO/DHI (thresholds), PH (residual
//   variances), DMU = E[A_d], MZ = E[A_d z_e] (column d*D + e), where
//   g = d log p / d a and A_d = sum_j lambda_jd g_j, expectations over the
//   posterior node weights.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

static inline int lt(int d, int e) { return d * (d + 1) / 2 + e; }   // 0-based, d >= e

// ---- category probabilities ---------------------------------------------------
// For a category with bounds lo < hi and linear predictor a, p = F(u) - F(v) with
// u = hi - a, v = lo - a (u = +Inf / v = -Inf for the open ends). To avoid
// cancellation, p is computed in the tail closest to the interval: if u + v > 0
// (interval mostly in the upper half), p = F(-v) - F(-u). Differences u, v are
// used directly, so no exponential of the predictor can overflow. If p underflows,
// log p and the derivative ratios f/p are computed on the log scale.
static inline void logis_Ff(double x, double &F, double &f) {   // F(x), f(x); x finite
  double e = std::exp(-std::fabs(x)), d = 1.0 + e;
  F = x >= 0 ? 1.0 / d : e / d; f = e / (d * d);
}
static inline double logF_(double x, int link) {
  if (x == R_PosInf) return 0.0;
  if (x == R_NegInf) return R_NegInf;
  if (link == 0) return x >= 0 ? -std::log1p(std::exp(-x)) : x - std::log1p(std::exp(x));
  return R::pnorm(x, 0.0, 1.0, 1, 1);
}
static inline double logf_(double x, int link) {
  if (!R_FINITE(x)) return R_NegInf;
  if (link == 0) { double ax = std::fabs(x); return -ax - 2.0 * std::log1p(std::exp(-ax)); }
  return R::dnorm(x, 0.0, 1.0, 1);
}
static inline double log1mexp_(double x) {      // log(1 - exp(x)), x <= 0
  if (x == R_NegInf) return 0.0;
  return x > -M_LN2 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
}
// returns true if p is returned for the running product (p > 1e-30, so the product
// cannot underflow before it is logged), false if logp is returned;
// rlo = f(v)/p, rhi = f(u)/p in both cases
static inline bool cat_prob(double u, double v, int link, double &p, double &logp, double &rlo, double &rhi) {
  const bool up = (u + v) > 0;
  const double b = up ? -v : u, c = up ? -u : v;        // p = F(b) - F(c), b > c, b finite
  double Fb, fb, Fc = 0.0, fc = 0.0;
  if (link == 0) { logis_Ff(b, Fb, fb); if (c != R_NegInf) logis_Ff(c, Fc, fc); }
  else { Fb = R::pnorm(b, 0.0, 1.0, 1, 0); fb = R::dnorm(b, 0.0, 1.0, 0);
         if (c != R_NegInf) { Fc = R::pnorm(c, 0.0, 1.0, 1, 0); fc = R::dnorm(c, 0.0, 1.0, 0); } }
  const double fu = up ? fc : fb, fv = up ? fb : fc;    // densities are symmetric
  p = Fb - Fc;
  if (p > 1e-250) {
    rlo = fv / p; rhi = fu / p;
    if (p > 1e-30) return true;
    logp = std::log(p); return false;
  }
  const double lFb = logF_(b, link), lFc = (c == R_NegInf) ? R_NegInf : logF_(c, link);
  double d = lFc - lFb; if (d > 0) d = 0;
  logp = lFb + log1mexp_(d);
  rlo = std::exp(logf_(v, link) - logp); rhi = std::exp(logf_(u, link) - logp);
  return false;
}

// [[Rcpp::export]]
List mofia_core_cpp(IntegerMatrix Y, NumericMatrix NU, NumericMatrix LAM, IntegerVector loadflag,
                    NumericVector tau_flat, IntegerVector tau_off,
                    NumericMatrix MU, NumericMatrix CL,
                    NumericMatrix Tn, NumericVector logw,
                    bool adaptive, NumericMatrix AM, NumericMatrix ALC,
                    int link, bool want_scores, bool want_post, int nthreads,
                    NumericMatrix YC, NumericMatrix RV, IntegerVector cont) {
  const int N = Y.nrow(), P = Y.ncol(), D = MU.ncol(), M = Tn.nrow(), T = D * (D + 1) / 2;
  NumericVector ll(N);
  NumericMatrix PM(want_post ? N : 0, D), PC(want_post ? N : 0, T);
  NumericMatrix PG(want_scores ? N : 0, P), PGE(want_scores ? N : 0, P * D),
                DLO(want_scores ? N : 0, P), DHI(want_scores ? N : 0, P), PH(want_scores ? N : 0, P),
                DMU(want_scores ? N : 0, D), MZ(want_scores ? N : 0, D * D);
  int bad = 0;
  const int *y = Y.begin(), *toff = tau_off.begin(), *lf = loadflag.begin(), *ct = cont.begin();
  const double *nu = NU.begin(), *lam = LAM.begin(), *mu = MU.begin(), *cl = CL.begin(), *tn = Tn.begin(),
               *lw0 = logw.begin(), *tf = tau_flat.begin(), *am = AM.begin(), *alc = ALC.begin(),
               *yc = YC.begin(), *rv = RV.begin();
  double *pll = ll.begin();

#ifdef _OPENMP
  if (nthreads < 1) nthreads = 1;
  #pragma omp parallel num_threads(nthreads) reduction(+:bad)
#endif
  {
    std::vector<double> z(M * D), lwv(M), E(M * D), logL(M), post(M), prodv(M), la(D), t(D);
    std::vector<double> gb, lob, hib, hcb, A;
    if (want_scores) { gb.assign(P * M, 0.0); lob.assign(P * M, 0.0); hib.assign(P * M, 0.0);
                       hcb.assign(P * M, 0.0); A.assign(D * M, 0.0); }
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int i = 0; i < N; ++i) {
      // ---- nodes, log weights, factor values ----
      double logdet = 0.0;
      if (adaptive) for (int d = 0; d < D; ++d) logdet += std::log(alc[i + (size_t)N * lt(d, d)]);
      for (int q = 0; q < M; ++q) {
        double w = lw0[q];
        for (int d = 0; d < D; ++d) t[d] = tn[q + (size_t)M * d];
        for (int d = 0; d < D; ++d) {
          double zd;
          if (adaptive) {
            zd = am[i + (size_t)N * d];
            for (int e = 0; e <= d; ++e) zd += alc[i + (size_t)N * lt(d, e)] * t[e];
            w += -0.5 * zd * zd + 0.5 * t[d] * t[d];
          } else zd = t[d];
          z[q + (size_t)M * d] = zd;
        }
        if (adaptive) w += logdet;
        lwv[q] = w;
        for (int d = 0; d < D; ++d) {
          double ed = mu[i + (size_t)N * d];
          for (int e = 0; e <= d; ++e) ed += cl[i + (size_t)N * lt(d, e)] * z[q + (size_t)M * e];
          E[q + (size_t)M * d] = ed;
        }
        logL[q] = 0.0; prodv[q] = 1.0;
      }
      // ---- items ----
      for (int j = 0; j < P; ++j) {
        double nij = nu[i + (size_t)N * j];
        for (int d = 0; d < D; ++d) la[d] = lf[j + P * d] ? lam[i + (size_t)N * (j + P * d)] : 0.0;
        if (ct[j]) {                                   // continuous item: normal density
          double yv = yc[i + (size_t)N * j];
          if (ISNAN(yv)) continue;
          double th = rv[i + (size_t)N * j], lc = -0.5 * std::log(2.0 * M_PI * th);
          for (int q = 0; q < M; ++q) {
            double a = nij; for (int d = 0; d < D; ++d) a += la[d] * E[q + (size_t)M * d];
            double r = yv - a;
            logL[q] += lc - r * r / (2.0 * th);
            if (want_scores) { gb[j * M + q] = r / th; hcb[j * M + q] = -0.5 + r * r / (2.0 * th);
                               lob[j * M + q] = 0.0; hib[j * M + q] = 0.0; }
          }
          continue;
        }
        int yy = y[i + (size_t)N * j];
        if (yy == NA_INTEGER) continue;
        const int K = toff[j + 1] - toff[j] + 1;
        const double lo = (yy == 0) ? R_NegInf : tf[toff[j] + yy - 1];
        const double hi = (yy == K - 1) ? R_PosInf : tf[toff[j] + yy];
        for (int q = 0; q < M; ++q) {
          double a = nij; for (int d = 0; d < D; ++d) a += la[d] * E[q + (size_t)M * d];
          double p, logp, rlo, rhi;
          if (cat_prob(hi - a, lo - a, link, p, logp, rlo, rhi)) {
            prodv[q] *= p;
            if (prodv[q] < 1e-250) { logL[q] += std::log(prodv[q]); prodv[q] = 1.0; }
          } else logL[q] += logp;
          if (want_scores) { gb[j * M + q] = rlo - rhi; lob[j * M + q] = -rlo; hib[j * M + q] = rhi; }
        }
      }
      for (int q = 0; q < M; ++q) logL[q] += std::log(prodv[q]);
      // ---- posterior over nodes ----
      double mx = R_NegInf;
      for (int q = 0; q < M; ++q) { post[q] = logL[q] + lwv[q]; if (post[q] > mx) mx = post[q]; }
      double sum = 0.0;
      for (int q = 0; q < M; ++q) sum += std::exp(post[q] - mx);
      double lli = mx + std::log(sum);
      if (!std::isfinite(lli)) bad += 1;
      pll[i] = lli;
      for (int q = 0; q < M; ++q) post[q] = std::exp(post[q] - lli);

      if (want_post) {
        for (int d = 0; d < D; ++d) {
          double md = 0; for (int q = 0; q < M; ++q) md += post[q] * z[q + (size_t)M * d];
          PM(i, d) = md;
        }
        for (int d = 0; d < D; ++d) for (int e = 0; e <= d; ++e) {
          double s2 = 0; for (int q = 0; q < M; ++q) s2 += post[q] * z[q + (size_t)M * d] * z[q + (size_t)M * e];
          PC(i, lt(d, e)) = s2 - PM(i, d) * PM(i, e);
        }
      }
      if (want_scores) {
        std::fill(A.begin(), A.end(), 0.0);
        for (int j = 0; j < P; ++j) {
          if (ct[j]) { if (ISNAN(yc[i + (size_t)N * j])) continue; }
          else { int yy = y[i + (size_t)N * j]; if (yy == NA_INTEGER) continue; }
          for (int d = 0; d < D; ++d) la[d] = lf[j + P * d] ? lam[i + (size_t)N * (j + P * d)] : 0.0;
          double spg = 0, sdlo = 0, sdhi = 0, sph = 0;
          for (int q = 0; q < M; ++q) {
            double g = gb[j * M + q], wq = post[q];
            spg += wq * g; sdlo += wq * lob[j * M + q]; sdhi += wq * hib[j * M + q];
            if (ct[j]) sph += wq * hcb[j * M + q];
            for (int d = 0; d < D; ++d) A[q + (size_t)M * d] += la[d] * g;
          }
          PG(i, j) = spg; DLO(i, j) = sdlo; DHI(i, j) = sdhi; if (ct[j]) PH(i, j) = sph;
          for (int d = 0; d < D; ++d) if (lf[j + P * d]) {
            double s = 0; for (int q = 0; q < M; ++q) s += post[q] * gb[j * M + q] * E[q + (size_t)M * d];
            PGE(i, j + P * d) = s;
          }
        }
        for (int d = 0; d < D; ++d) {
          double dm = 0; for (int q = 0; q < M; ++q) dm += post[q] * A[q + (size_t)M * d];
          DMU(i, d) = dm;
          for (int e = 0; e < D; ++e) {
            double s = 0; for (int q = 0; q < M; ++q) s += post[q] * A[q + (size_t)M * d] * z[q + (size_t)M * e];
            MZ(i, d * D + e) = s;
          }
        }
      }
    }
  }
  return List::create(_["ll"] = ll, _["bad"] = bad, _["PM"] = PM, _["PC"] = PC,
                      _["PG"] = PG, _["PGE"] = PGE, _["DLO"] = DLO, _["DHI"] = DHI, _["PH"] = PH,
                      _["DMU"] = DMU, _["MZ"] = MZ);
}

// [[Rcpp::export]]
int mofia_openmp_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 0;
#endif
}
