# Function to compute the MMD^2 U-estimator of curve matching
Rcpp::cppFunction(
  "double mmdcurve_c(NumericMatrix y, NumericMatrix z,
                      NumericVector t_y, NumericVector t_z,
                      double sigma, double lambda){
   int nobs = y.rows();
   int dim  = y.cols();
   if (z.rows() != nobs || z.cols() != dim)  stop(\"y and z dimensions differ\");
   if (t_y.size() != nobs || t_z.size() != nobs) stop(\"Time index length mismatch\");

   double yy = 0.0, zz = 0.0, yz = 0.0;

   for (int i1 = 0; i1 < nobs; ++i1){
     for (int i2 = 0; i2 < nobs; ++i2){
       // Compute squared spatial distance
       double vert_yy = 0.0, vert_zz = 0.0, vert_yz = 0.0;
       for (int j = 0; j < dim; ++j){
         if (i1 != i2){
           vert_yy += std::pow(y(i1, j) - y(i2, j), 2.0);
           vert_zz += std::pow(z(i1, j) - z(i2, j), 2.0);
         }
         vert_yz += std::pow(y(i1, j) - z(i2, j), 2.0);
       }

       // Compute squared temporal distance scaled by lambda
       double dt_yy = t_y[i1] - t_y[i2];
       double horiz2_yy = lambda * dt_yy * dt_yy;  // λ·|Δt|^2

       double dt_zz = t_z[i1] - t_z[i2];
       double horiz2_zz = lambda * dt_zz * dt_zz;

       double dt_yz = t_y[i1] - t_z[i2];
       double horiz2_yz = lambda * dt_yz * dt_yz;

       // Gaussian kernel
       double k_yz = std::exp(-(vert_yz + horiz2_yz) / (2.0 * sigma * sigma));
       yz += k_yz;

       if (i1 != i2){
         double k_yy = std::exp(-(vert_yy + horiz2_yy) / (2.0 * sigma * sigma));
         double k_zz = std::exp(-(vert_zz + horiz2_zz) / (2.0 * sigma * sigma));
         yy += k_yy;
         zz += k_zz;
       }
     }
   }

   double denom = static_cast<double>(nobs) * (nobs - 1);
   double result = yy / denom + zz / denom - 2.0 * yz / (nobs * nobs);
   return result;
}"
)
