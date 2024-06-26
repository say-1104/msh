#include "main.hpp"

#include "LI.hpp"
// #include "Structure.hpp"

int main(int argc, char* argv[]) {
  int i, j, k;

  auto output_err = [](int step, double cur_z, std::complex<double> ab0,
                       std::complex<double> ab1) -> void { // ラムダ式
    std::cerr << std::fixed << std::setprecision(10);
    std::cerr << step << '\t';
    std::cerr << cur_z << '\t';
    std::cerr << abs(ab0) * abs(ab0) << '\t';
    std::cerr << abs(ab1) * abs(ab1) << '\n';
  };

  //input_[波長].pre の読み込み
  int N_Mat;
  std::cin >> N_Mat;
  
  int N;
  LI<double, double> beta_even;
  LI<double, double> beta_odd;
  LI<double, double> beta_1;
  LI<double, double> beta_2;
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
  
  LI<double, double> beta_even_;
  LI<double, double> beta_odd_;
  LI<double, double> beta_1_;
  LI<double, double> beta_2_;

  for(int n = 0; n < N_Mat; n++){
    
  }
  

  if(str_name == "pcmDC"){
    std::cin >> N_;
    for (int _ = 0; _ < N_; _++) {
      double w;
      double b1, b2, b3, b4;
      std::cin >> w >> b1 >> b2 >> b3 >> b4;
      beta_even_.append(w, b1);
      beta_odd_.append(w, b2);
      beta_1_.append(w, b3);
      beta_2_.append(w, b4);
    }
  }

  //structure.cfgの読み込み
  std::ifstream ifs("structure.cfg");
  if (!ifs) {
    std::cerr << "Error: cannot open strcture.cfg" << std::endl;
    exit(1);
  }
  std::string buff;        //必要ない文字列の格納場所
  std::string str_name;   //構造名
  ifs >> str_name;
  double Lc, dz, wst, Leff;
  int str_N = 0;
std::vector<int> Mat;
  if(str_name == "DC") {
    ifs >> buff >> Lc >> buff >> dz >> buff >> wst;
    std::cerr << "Lc: " << Lc << "\ndz: " << dz << "\nw_start: " << wst << std::endl;
    str_N++;
    Mat.push_back();
  } else if(str_name == "pcmDC") {
    ifs >> buff >> Lc >> buff >> dz >> buff >> wst >> buff >> Leff;
    str_N += 2;
  }
  double n_div = Lc / dz;




  
  for(int _=0; _++; _<str_N){
    ifs >> str_name;
    if (str_name == "DC"){
      

      

      std::vector<std::pair<double, double>> ZtoW(n_div);
      for (int step = 0; step < n_div; step++) {
        double cur_z = step * dz;
        double cur_w = wst;
        ZtoW[step] = std::make_pair(cur_z, cur_w);
      }
    }

  }


  
  
  Eigen::Vector2cd ab;
  //ab << std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0);




  

  double Lc = 36.90;
  double dz = 0.010;
  double n_div = Lc / dz;
  double wst = 0.794;
  double wfi = 0.836 + 0.042;

  std::cerr << "step\tposition_z\tabs(a)\tabs(b)\n";


  std::vector<std::pair<double, double>> ZtoW(n_div);

  for (int str = 0; str < str_N; str++){
    for (int step = 0; step < n_div[str]; step++) {
      auto z = ZtoW[str][step].first;
      auto w = ZtoW[str][step].second;

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
  }
  output_err(n_div, Lc, ab(0), ab(1));
}