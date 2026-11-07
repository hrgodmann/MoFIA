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

static inline double plogis_(double x) {
  if (x >= 0) { double e = std::exp(-x); return 1.0 / (1.0 + e); }
  double e = std::exp(x); return e / (1.0 + e);
}
static inline double Fcdf(double x, int link) {
  if (x == R_PosInf) return 1.0;
  if (x == R_NegInf) return 0.0;
  return link == 0 ? plogis_(x) : R::pnorm(x, 0.0, 1.0, 1, 0);
}
static inline double fpdf(double x, int link) {
  if (!R_FINITE(x)) return 0.0;
  if (link == 0) { double e = std::exp(-std::fabs(x)); return e / ((1.0 + e) * (1.0 + e)); }
  return R::dnorm(x, 0.0, 1.0, 0);
}
static inline int lt(int d, int e) { return d * (d + 1) / 2 + e; }   // 0-based, d >= e

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
        int K = toff[j + 1] - toff[j] + 1;
        bool lo_inf = (yy == 0), hi_inf = (yy == K - 1);
        double lo = lo_inf ? R_NegInf : tf[toff[j] + yy - 1];
        double hi = hi_inf ? R_PosInf : tf[toff[j] + yy];
        if (link == 0) {
          // logit: with u = exp(a), F(hi - a) = ehi / (ehi + u), F(lo - a) = elo / (elo + u),
          // p = u (ehi - elo) / ((ehi + u)(elo + u)); f(x) = F(x)(1 - F(x)). One exp per node.
          double elo = lo_inf ? 0.0 : std::exp(lo), ehi = hi_inf ? 0.0 : std::exp(hi);
          for (int q = 0; q < M; ++q) {
            double a = nij; for (int d = 0; d < D; ++d) a += la[d] * E[q + (size_t)M * d];
            if (a > 700.0) a = 700.0; else if (a < -700.0) a = -700.0;
            double u = std::exp(a), p, flo = 0.0, fhi = 0.0;
            if (lo_inf && hi_inf) p = 1.0;
            else if (lo_inf)      { double dd = ehi + u; p = ehi / dd; fhi = ehi * u / (dd * dd); }
            else if (hi_inf)      { double dd = elo + u; p = u / dd;   flo = elo * u / (dd * dd); }
            else { double dh = ehi + u, dl = elo + u; p = u * (ehi - elo) / (dh * dl);
                   fhi = ehi * u / (dh * dh); flo = elo * u / (dl * dl); }
            if (!(p > 1e-300)) p = 1e-300;
            prodv[q] *= p;
            if (prodv[q] < 1e-250) { logL[q] += std::log(prodv[q]); prodv[q] = 1.0; }
            if (want_scores) { gb[j * M + q] = (flo - fhi) / p; lob[j * M + q] = -flo / p; hib[j * M + q] = fhi / p; }
          }
        } else {
          for (int q = 0; q < M; ++q) {
            double a = nij; for (int d = 0; d < D; ++d) a += la[d] * E[q + (size_t)M * d];
            double p = Fcdf(hi - a, link) - Fcdf(lo - a, link);
            if (!(p > 1e-300)) p = 1e-300;
            prodv[q] *= p;
            if (prodv[q] < 1e-250) { logL[q] += std::log(prodv[q]); prodv[q] = 1.0; }
            if (want_scores) {
              double flo = fpdf(lo - a, link), fhi = fpdf(hi - a, link);
              gb[j * M + q] = (flo - fhi) / p; lob[j * M + q] = -flo / p; hib[j * M + q] = fhi / p;
            }
          }
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
