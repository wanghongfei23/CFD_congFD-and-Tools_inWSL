#include <iostream>
#include <fstream>
using namespace std;
int main(int, char**){
    fstream in,out;
    in.open("data.txt",ios::in);

    int n,nLine;
    in>>n;
    in>>nLine;

    double** data=new double* [nLine];
    for(int j=0;j<nLine;j++)
    {
        data[j]=new double[n];
        for(int i=0;i<n;i++)
        {
            in>>data[j][i];
        }
    }

    
    for(int i=1;i<nLine;i++)
    {
        string str1="./tec/tec";
        string str2=".dat";
        out.open(str1+to_string(i)+str2,ios::out);
        out<<"Title=\"tecplot dataformat"<<i<<"\"\n";
        out<<"Variables=\"x\",\"u\"\n";
        out<<"Zone I="<<n<<" F=POINT\n";
        for(int j=0;j<n;j++)
        {
            out<<data[0][j]<<' ';
            out<<data[i][j]<<'\n';
        }
        out.close();
    }

    for(int j=0;j<nLine;j++)
    {
        delete[] data[j];
    }
    delete[] data;

    return 0;

}
