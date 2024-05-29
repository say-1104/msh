#include "main.hpp"

#include "LI.hpp"
// #include "Structure.hpp"

int main(int argc, char* argv[]) {
  

  Eigen::Vector2cd ab;
  auto output_err = [](int step, double cur_z, std::complex<double> ab0,
                      std::complex<double> ab1) -> void { // ラムダ式
    std::cerr << std::fixed << std::setprecision(10);
    std::cerr << step << '\t';
    std::cerr << cur_z << '\t';
    std::cerr << abs(ab0) * abs(ab0) << '\t';
    std::cerr << abs(ab1) * abs(ab1) << '\n';
  };
  LI<double, double> beta_even;
  LI<double, double> beta_odd;
  LI<double, double> beta_1;
  LI<double, double> beta_2;
  int i,j;
  int stage;
  double dz,w_bus,w_acc;
  int Sbend_number,coupler_number;
  double **Str;
  double l;
  int *stagetype;
  int Sbend_Str_n = 3;
  int Coupler_Str_n = 2;
  int Str_n;
  int N;
  int all_step = 0;
  double all_z;
  std::ofstream writing_file;
  char filename[256];


  char string[256];
  FILE *fp;

  if ((fp = fopen("Structure", "r")) == NULL){
    fprintf(stderr,"cannnot open file\n");
    exit(0);
  }
  // std::cerr << beta_1_Sbend[0.400] << '\n';

  fscanf(fp, "%*s = %d", &stage);
  fscanf(fp, "%*s = %lf", &dz);
  fscanf(fp, "%*s = %lf", &w_bus);
  fscanf(fp, "%*s = %lf", &w_acc);
  fscanf(fp, "%*s = %d", &Sbend_number);
  fscanf(fp, "%*s = %d", &coupler_number);

  // std::cerr << stage << "\t" << dz << "\t" << w_bus << "\t" << w_acc << "\t" << Sbend_number << "\t" << coupler_number << "\n";

  // std::cerr << "stagetype" << "\n";
  ALLOCATION(stagetype, int, stage);
  for (i = 0; i < stage; i++){
    fscanf(fp, "%*s = %d", &stagetype[i]);
    // std::cerr << stagetype[i] << "\n";
  }
  // std::cerr << "Str" << "\n";
  ALLOCATION(Str, double *,stage);
  for (i = 0; i < stage; i++){
    fscanf(fp, "\n");
    if((stagetype[i] == 1) || (stagetype[i] == 3)){
      Str_n = Sbend_Str_n;
    }else if((stagetype[i] == 2)){
      Str_n = Coupler_Str_n;
    }
    // std::cerr << "Str" << i << "\n";
    ALLOCATION(Str[i], double, Str_n);
    for(j = 0; j < Str_n; j++){
      fscanf(fp, "%*s = %lf", &Str[i][j]);
      // std::cerr << Str[i][j] << "\n";
    }

  }
  fclose(fp);
  std::cerr << "step\tposition_z\tabs(a)\tabs(b)\n";
  for(i = 0; i < stage; i++){
    if ((stagetype[i] == 1) || (stagetype[i] == 3)){
      sprintf(string, "./data_Sbend/input_%s.pre",argv[1]);
    }else if(stagetype[i] == 2){
      sprintf(string, "./data/input_%s.pre", argv[1]);
    }
    if ((fp = fopen(string, "r")) == NULL){
    std::cerr << "can't open file" << string << "\n";
    exit(1);
    }

    fscanf(fp, "%d", &N);
    for (int _ = 0; _ < N; _++) {
      double w;
      double b1, b2, b3, b4;
      fscanf(fp, "%lf\t%lf\t%lf\t%lf\t%lf", &w, &b1, &b2, &b3, &b4);
      // std::cerr << w << '\t' << b1 << '\t' << b2 << '\t' << b3 << '\t' << b4 << '\n';
      if((stagetype[i] == 1) || (stagetype[i] == 3)){
        if(w == 400){
          w = 0.400;
        }
        w = w * 0.001;
      }
      beta_even.append(w, b1);
      beta_odd.append(w, b2);
      beta_1.append(w, b3);
      beta_2.append(w, b4);
      // std::cerr << w << '\t' << beta_even_Sbend[w] << '\n';
    }
    fclose(fp);
    
    int n_div = Str[i][0]/dz;
    std::vector<std::pair<double, double>> ZtoW(n_div);

    // 構造を記述する
    for (int step = 0; step < n_div; step++){
      double cur_z = step*dz;
      double cur_w;
      if (stagetype[i] == 1){
        double r,bufc11, theta;
        double ld = Str[i][0];
        double g = Str[i][1];
        double sd = Str[i][2];
        l = ld;

        r = (ld*ld+sd*sd)/(4.0*sd);
        bufc11 = g / 2.0 + w_bus/2.0+sd;
        if (cur_z < ld / 2.0){
          theta = asin(cur_z / r);
          cur_w = g + 2*sd-2*r+2*r*cos(theta);
        }else if(cur_z < ld){
          theta=asin((ld-cur_z)/r);
          cur_w=g+2*r-2*r*cos(theta);
        }
      }else if ( stagetype[i] == 2){
        double Lc = Str[i][0];
        double dw = Str[i][1];
        l = Lc;
        if ( cur_z < Lc / 4.0){
          cur_w = w_bus - 4 * dw * cur_z / Lc;
        }else if ( cur_z < 3*Lc/4.0){
          cur_w = w_bus - 2 * dw + 4*dw*(cur_z)/ Lc;
        }else if ( cur_z < Lc){
          cur_w = w_bus - 4*dw*cur_z/Lc + 4 * dw;
        }
      }else if( stagetype[i] == 3){
        double r,bufc11, theta;
        double ld = Str[i][0];
        double g = Str[i][1];
        double sd = Str[i][2];
        l = ld;

        r = (ld*ld+sd*sd)/(4.0*sd);
        bufc11 = g / 2.0 + w_bus/2.0+sd;
        if (cur_z < ld / 2.0){
          theta=asin((cur_z)/r);
          cur_w=g+2*r-2*r*cos(theta);
        }else if(cur_z < ld){
          theta=asin((ld-cur_z)/r);
          cur_w = g + 2*sd-2*r+2*r*cos(theta);
        }
      }
      ZtoW[step] = std::make_pair(cur_z, cur_w);
    }
    if (i == 0){
      ab << std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0);
    }else{
      double a_re,a_im,b_re,b_im;
      sprintf(string, "./result/com_amp_end.dat");
      if ((fp = fopen(string,"r")) == NULL){
        std::cerr << "can't open file" << string << "\n";
        exit(1);
      }
      while (fscanf(fp, "%*s\t%*s\t(%lf,%lf)\t(%lf,%lf)",&a_re,&a_im,&b_re,&b_im) != EOF){
        // std::cerr << a_re << "\t" << a_im << "\t" << b_re << "\t" << b_im << "\n";
      }
      ab << std::complex<double>(a_re,a_im), std::complex<double>(b_re,b_im);
    }

    for ( int step = 0; step < n_div; step++){
      double z,w;
      std::tie(z,w) = ZtoW[step];
      // std::cerr << all_z << "\t" << w << "\n";
      double be,bo,b1,b2;

      be = beta_even[w];
      bo = beta_odd[w];
      b1 = beta_1[w];
      b2 = beta_2[w];

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

      // sprintf(filename, "./result/com_amp_%d.dat",i);
      // // filename << "./result/com_amp_" << i << ".dat";
      // writing_file.open(filename, std::ios::app);
      // if(writing_file == NULL){
      //   std::cerr << "can't open file" << filename << "\n";
      //   exit(1);
      // }

      // writing_file << step << "\t" << z << "\t" << ab(0) << "\t" << ab(1) << "\n";
      // writing_file.close();

      // output_err(all_step,all_z,ab(0),ab(1));
      all_step = 1 + all_step;
      all_z = all_step*dz;
      calc_CMT();
    }

    sprintf(filename, "./result/com_amp_%d.dat",i);
    writing_file.open(filename, std::ios::app);
    if(writing_file == NULL){
      std::cerr << "can't open file" << filename << "\n";
      exit(1);
    }
    writing_file << n_div << "\t" << l << "\t" << ab(0) << "\t" << ab(1) << "\n" ;
    writing_file.close();

    // output_err(all_step,all_z,ab(0),ab(1));

    sprintf(filename, "./result/com_amp_end.dat");
    // filename = "./result/com_amp_end.dat";
    writing_file.open(filename,std::ios::app);
    if(writing_file == NULL){
      std::cerr << "can't open file" << filename << "\n";
      exit(1);
    }
    writing_file << n_div << "\t" << l << "\t" << ab(0) << "\t" << ab(1) << "\n" ;
    writing_file.close();
  }
  output_err(all_step,all_z,ab(0),ab(1));
}