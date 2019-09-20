

auto equal = [](auto x, auto y, double tol = 1e-6){ return x == Approx(y).epsilon(tol); };
