# Function to compute the MMD^2 by the V-statistic estimator
Rcpp::cppFunction("
  double mmd_V(NumericMatrix y, NumericMatrix z, double epsilon){
    int nobs      = y.rows();
    int dimension = y.cols();
    double yy = 0.0, zz = 0.0, yz = 0.0;

    for (int i1 = 0; i1 < nobs; i1++){
      for (int i2 = 0; i2 < nobs; i2++){
        double cost_yy = 0.0, cost_zz = 0.0, cost_yz = 0.0;
        for (int j = 0; j < dimension; j++){
          cost_yy += std::pow(y(i1, j) - y(i2, j), 2.0);
          cost_zz += std::pow(z(i1, j) - z(i2, j), 2.0);
          cost_yz += std::pow(y(i1, j) - z(i2, j), 2.0);
        }
        yy += std::exp(- cost_yy / (2.0 * epsilon * epsilon));
        zz += std::exp(- cost_zz / (2.0 * epsilon * epsilon));
        yz += std::exp(- cost_yz / (2.0 * epsilon * epsilon));
      }
    }

    double norm = static_cast<double>(nobs) * static_cast<double>(nobs);
    double result = yy / norm + zz / norm - 2.0 * yz / norm;
    return result;
  }
")
