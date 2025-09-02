Rcpp::cppFunction('
  double mmd_V_delay(Rcpp::NumericMatrix y, Rcpp::NumericMatrix z, double epsilon){
    int n1 = y.nrow();
    int n2 = z.nrow();
    int d  = y.ncol();

    if (z.ncol() != d) Rcpp::stop("y and z must have the same number of columns");
    if (n1 <= 0 || n2 <= 0) Rcpp::stop("y and z must have at least one row");
    if (!(epsilon > 0.0))   Rcpp::stop("epsilon must be positive");

    const double inv2eps2 = 1.0 / (2.0 * epsilon * epsilon);

    double sum_yy = 0.0;
    double sum_zz = 0.0;
    double sum_yz = 0.0;

    // YY
    for (int i = 0; i < n1; ++i){
      for (int j = 0; j < n1; ++j){
        double dist2 = 0.0;
        for (int p = 0; p < d; ++p){
          double diff = y(i,p) - y(j,p);
          dist2 += diff * diff;
        }
        sum_yy += std::exp(- dist2 * inv2eps2);
      }
    }

    // ZZ
    for (int i = 0; i < n2; ++i){
      for (int j = 0; j < n2; ++j){
        double dist2 = 0.0;
        for (int p = 0; p < d; ++p){
          double diff = z(i,p) - z(j,p);
          dist2 += diff * diff;
        }
        sum_zz += std::exp(- dist2 * inv2eps2);
      }
    }

    // YZ
    for (int i = 0; i < n1; ++i){
      for (int j = 0; j < n2; ++j){
        double dist2 = 0.0;
        for (int p = 0; p < d; ++p){
          double diff = y(i,p) - z(j,p);
          dist2 += diff * diff;
        }
        sum_yz += std::exp(- dist2 * inv2eps2);
      }
    }

    // V-statistics 
    double mmd2 = (sum_yy / ( (double) n1 * n1 ))
                - 2.0 * (sum_yz / ( (double) n1 * n2 ))
                + (sum_zz / ( (double) n2 * n2 ));

    return mmd2;
  }
')
