#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    ifstream in("input.txt");
    ofstream out("output.txt");
    int nx, ny;

    nx = 101;
    ny = 51;

    double T[nx*ny][3];
    double a[ny][nx];


    for (int i = 0; i < nx*ny; i++)
    {
        double b;
        in >> b;
        T[i][0] = b;
        in >> b;
        T[i][1] = b;
        in >> b;
        T[i][2] = b;
    }


    int k;

    do{

       k = 0;

       for (int i = 0; i < nx*ny - 1; i++)
       {
           if (T[i][0] > T[i + 1][0])
           {
              double t;
              t = T[i][0];
              T[i][0] = T[i + 1][0];
              T[i + 1][0] = t;

              t = T[i][1];
              T[i][1] = T[i + 1][1];
              T[i + 1][1] = t;

              t = T[i][2];
              T[i][2] = T[i + 1][2];
              T[i + 1][2] = t;
           }
           else
           {
              k++;
           }
       }


    }while(k != (nx*ny - 1));




    for (int j = 0; j < nx; j++)
    {

      do{

       k = 0;

       for (int i = j*ny; i < (j+1)*ny - 1; i++)
       {
           if (T[i][1] > T[i + 1][1])
           {
              double t;
              t = T[i][0];
              T[i][0] = T[i + 1][0];
              T[i + 1][0] = t;

              t = T[i][1];
              T[i][1] = T[i + 1][1];
              T[i + 1][1] = t;

              t = T[i][2];
              T[i][2] = T[i + 1][2];
              T[i + 1][2] = t;
           }
           else
           {
              k++;
           }
       }


       }while(k != (ny - 1));

    }





    for (int j = 0; j < nx; j++)
        for (int i = 0; i < ny; i++)
        {
            a[i][j] = T[i + j*ny][2];
        }



    for (int i = ny - 1; i > 0; i--)
    {

        for (int j = 0; j < nx; j++)
        {
            out << a[i][j] << " ";
        }

        out << endl;
    }

    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            out << a[i][j] << " ";
        }

        out << endl;

    }

    in.close();
    out.close();
    delete [] a,T;

    return 0;
}
