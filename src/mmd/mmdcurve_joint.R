Rcpp::cppFunction('
double mmdcurve_vstat_c(Rcpp::NumericMatrix y, Rcpp::NumericMatrix z,
                        Rcpp::NumericVector t_y, Rcpp::NumericVector t_z,
                        double ell_t, double ell_x) {

  int ny  = y.nrow();
  int nz  = z.nrow();
  int dim = y.ncol();

  if (z.ncol() != dim)                           Rcpp::stop("y and z dimensions differ");
  if (t_y.size() != ny)                          Rcpp::stop("t_y length mismatch");
  if (t_z.size() != nz)                          Rcpp::stop("t_z length mismatch");
  if (ell_t <= 0.0 || ell_x <= 0.0)              Rcpp::stop("ell_t and ell_x must be positive");

  // 1/(2*ell^2)
  const double inv2lt2 = 1.0 / (2.0 * ell_t * ell_t);
  const double inv2lx2 = 1.0 / (2.0 * ell_x * ell_x);

  double sum_yy = 0.0;
  double sum_zz = 0.0;
  double sum_yz = 0.0;

  // YY
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < ny; ++j) {
      double vert = 0.0;
      for (int d = 0; d < dim; ++d) {
        double diff = y(i,d) - y(j,d);
        vert += diff * diff;
      }
      double dt = t_y[i] - t_y[j];
      double k  = std::exp( - dt*dt * inv2lt2 - vert * inv2lx2 );
      sum_yy += k;
    }
  }

  // ZZ
  for (int i = 0; i < nz; ++i) {
    for (int j = 0; j < nz; ++j) {
      double vert = 0.0;
      for (int d = 0; d < dim; ++d) {
        double diff = z(i,d) - z(j,d);
        vert += diff * diff;
      }
      double dt = t_z[i] - t_z[j];
      double k  = std::exp( - dt*dt * inv2lt2 - vert * inv2lx2 );
      sum_zz += k;
    }
  }

  // YZ
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nz; ++j) {
      double vert = 0.0;
      for (int d = 0; d < dim; ++d) {
        double diff = y(i,d) - z(j,d);
        vert += diff * diff;
      }
      double dt = t_y[i] - t_z[j];
      double k  = std::exp( - dt*dt * inv2lt2 - vert * inv2lx2 );
      sum_yz += k;
    }
  }

  // V-statistics
  double mmd2 = (sum_yy / ( (double)ny * ny ))
              - 2.0 * (sum_yz / ( (double)ny * nz ))
              + (sum_zz / ( (double)nz * nz ));

  return mmd2;
}
')

