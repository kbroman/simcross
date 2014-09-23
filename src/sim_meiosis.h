NumericVector fromR_sim_crossovers(const double L, const int m, const double p,
                                   const bool obligate_chiasma, const double Lstar);
NumericVector cpp_sim_crossovers(const double L, const int m, const double p,
                                 const bool obligate_chiasma, const double Lstar);
List fromR_sim_meiosis(const List parent, const int m, const double p,
                       const bool obligate_chiasma, const double Lstar);

