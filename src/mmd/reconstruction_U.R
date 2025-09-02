Rcpp::cppFunction('
  double mmd_U_delay(NumericMatrix y, NumericMatrix z, double epsilon){
    int n1 = y.nrow();
    int n2 = z.nrow();
    int d  = y.ncol();
    if (z.ncol() != d) stop("y and z must have the same number of columns");

    double yy = 0.0, zz = 0.0, yz = 0.0;

    // YY 
    for (int i=0;i<n1;i++){
      for (int j=0;j<n1;j++){
        if (i==j) continue;
        double dist2 = 0.0;
        for (int p=0;p<d;p++){
          double diff = y(i,p) - y(j,p);
          dist2 += diff*diff;
        }
        yy += exp(- dist2 / (2.0*epsilon*epsilon));
      }
    }

    // ZZ 
    for (int i=0;i<n2;i++){
      for (int j=0;j<n2;j++){
        if (i==j) continue;
        double dist2 = 0.0;
        for (int p=0;p<d;p++){
          double diff = z(i,p) - z(j,p);
          dist2 += diff*diff;
        }
        zz += exp(- dist2 / (2.0*epsilon*epsilon));
      }
    }

    // YZ 
    for (int i=0;i<n1;i++){
      for (int j=0;j<n2;j++){
        double dist2 = 0.0;
        for (int p=0;p<d;p++){
          double diff = y(i,p) - z(j,p);
          dist2 += diff*diff;
        }
        yz += exp(- dist2 / (2.0*epsilon*epsilon));
      }
    }
    
    // U-statistics
    double mmd2 = yy / ( (double)n1*(n1-1) )
                + zz / ( (double)n2*(n2-1) )
                - 2.0 * yz / ( (double)n1*n2 );
    return mmd2;
  }
')

