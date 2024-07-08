#include "main.hpp"

int main(int argc, char* argv[]){
    DataTable data;
    //引数の確認
    checkArgument(&data, argc, argv);
    checkConfig(&data);

    inputData(&data);

    makeZtoW(&data);

    calcCMT(&data);
    return 0;
}