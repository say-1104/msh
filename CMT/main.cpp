#include "main.hpp"

#include "LI.hpp"
// #include "Structure.hpp"

int main(int argc, char* argv[]) {
  auto output_err = [](int step, double cur_z, std::complex<double> ab0,
                       std::complex<double> ab1) -> void { // ラムダ式
    std::cerr << std::fixed << std::setprecision(10);
    std::cerr << step << '\t';
    std::cerr << cur_z << '\t';
    std::cerr << abs(ab0) * abs(ab0) << '\t';
    std::cerr << abs(ab1) * abs(ab1) << '\n';
  };

  
  Eigen::Vector2cd ab;
  ab << std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0);

  LI<double, double> beta_even;
  LI<double, double> beta_odd;
  LI<double, double> beta_1;
  LI<double, double> beta_2;

  int N;
  std::cin >> N;
  for (int _ = 0; _ < N; _++) {
    double w;
    double b1, b2, b3, b4;
    std::cin >> w >> b1 >> b2 >> b3 >> b4;
    beta_even.append(w, b1);
    beta_odd.append(w, b2);
    beta_1.append(w, b3);
    beta_2.append(w, b4);
  }

  double Lc = 36.90;
  double dz = 0.010;
  double n_div = Lc / dz;
  double wst = 0.794;
  double wfi = 0.836 + 0.042;

  std::cerr << "step\tposition_z\tabs(a)\tabs(b)\n";


  std::vector<std::pair<double, double>> ZtoW(n_div);
  for (int step = 0; step < n_div; step++) {
    double cur_z = step * dz;
    double cur_w = wst + (cur_z / Lc) * (wfi - wst);
    ZtoW[step] = std::make_pair(cur_z, cur_w);
  }

  for (int step = 0; step < n_div; step++) {
    auto z = ZtoW[step].first;
    auto w = ZtoW[step].second;

    double be = beta_even[w];
    double bo = beta_odd[w];
    double b1 = beta_1[w];
    double b2 = beta_2[w];

    auto calc_CMT = [&]() -> void { // ラムダ式
      double beta_ave = (b1 + b2) / 2.0;
      double delta = (b1 - b2) / 2.0;
      double q = (be - bo) / 2.0;
      double kappa =
          ((q * q - delta * delta >= 0.0) ? sqrt(q * q - delta * delta) : 0.0);

      Eigen::Matrix2d P;
      P << delta + q, kappa, kappa, -(delta + q);
      P *= 1.0 / (sqrt(kappa * kappa + (delta + q) * (delta + q)));

      Eigen::Matrix2cd mid;
      mid << exp(-cj * (beta_ave + q) * dz), std::complex<double>(0.0, 0.0),
          std::complex<double>(0.0, 0.0), exp(-cj * (beta_ave - q) * dz);

      ab = P * ab;
      ab = mid * ab;
      ab = P * ab;
    };

    output_err(step, z, ab(0), ab(1));
    calc_CMT();
  }
  output_err(n_div, Lc, ab(0), ab(1));
}