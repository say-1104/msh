typedef struct _Param {
    double wl;  //波長
    int N_dset; //データセットの個数
    double dz;  //伝搬方向の分割数
    double Leff;    //cPCMの長さ
    double wst;     //始点のwidth
    int N_taper;    //テーパの接続個数
    std::vector<std::pair<double, double>> taper;   //テーパの構造パラメータ
    std::vector<std::tuple<double, double, int>> ZtoW;
} Param;

typedef struct _Dataset {
    int N;  //widthの個数
    LI<double, double> beta_even;
    LI<double, double> beta_odd;
    LI<double, double> beta_1;
    LI<double, double> beta_2;
} Dataset;

typedef struct _Flag {
    int pcm;    //0なら堆積無し、1なら片側堆積、2なら両側堆積
} Flag;

typedef struct _DataTable {
    Flag flag;
    Param par;
    Dataset *dset;
} DataTable;