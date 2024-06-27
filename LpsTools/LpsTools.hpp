#pragma once
#include <fstream>
#include <string>
#include <vector>

enum class Shape {
    point, line, surface, volume
};

enum class axis {
    x, y, z
};

enum class Function {
    point, line, surface, volume, copy, rotation, material, unstr, str, none
};

std::vector<int> Plusv(std::vector<int> a, int n);

class LpsTools {
    private:
        std::string lpsname;        //lpsファイル名
        Function currentfunc;       //現在展開中の機能

        std::string getFuncname(Function func);

    public:
        std::fstream ofs;

        double y_offset;
        int step_offset;
        std::vector<std::vector<int> > step = std::vector<std::vector<int> >(4, std::vector<int>(1, 0));
        double unstr, trans, gene;
        

        LpsTools(const char* name, double unstr, double trans, double gene);

        void SwitchingFunc(Function func);
        std::vector<int> Atovec(const char* text);      //文字列をベクタに変換
        void Appendstep(int p, int l, int s, int v);    //step値を更新
        void Printstep(void);

        void Point(double x, double y, double z);
        void Point(double x, double z);
        void Line(double x1, double y1, double z1, double x2, double y2, double z2);
        void Line(double x1, double z1, double x2, double z2);
        void Surface(const char* text);
        void Surface(std::vector<int> v);
        void Surface(int st, int end, int inc=1);
        void Volume(const char* text);
        void Volume(std::vector<int> v);
        void Volume(int st, int end, int inc=1);

        void Copy(Shape shape, double x, double y, double z, const char* text);
        void Copy(Shape shape, double x, double y, double z, std::vector<int> v);
        void Copy(Shape shape, double x, double y, double z, int st, int end, int inc=1);

        void Mat2D(int mat, const char* text);
        void Mat2D(int mat, std::vector<int> v);
        void Mat3D(int mat, const char* text);
        void Mat3D(int mat, std::vector<int> v);

        void Unstr(Shape shape, double unstr, const char* text);
        void Unstr(Shape shape, double unstr, std::vector<int> v);

        void TransGene();
        void Fileclose();
};