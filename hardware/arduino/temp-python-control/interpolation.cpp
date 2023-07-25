#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main()
{
double v, xb, yb,p;
int n=0, i=0, j;
fstream fil, doc;
fil.open("measured-energy.dat", ios::in);

do
{
    fil >> v;
    n++;
}while(!fil.eof());
n=n-1;
fil.close();
double x[n/2];
double y[n/2];
double f[n/2][n/2];

fil.open("measured-energy.dat", ios::in);

//read in values
do
{
    fil >> x[i];
    fil >> y[i];
    i++;
}while(i<=n/2);


fil.close();

//set up starting values for matrix
i=0;
do
{   f[i][0]=y[i];
    i++;
}while(i<=n/2);

//calculate coefficients
j=1;
do{
    i=j;
    do
    {   f[i][j]=(f[i][j-1]-f[i-1][j-1])/(x[i]-x[i-j]);
        i++;
    }while(i<=n/2);
    j++;
}while(j<=n/2);


doc.open("interpolation.dat", ios::out);

//return file which contains real voltages respective to requested voltages
//lowestVoltage = ?
//highestVoltage = ?
//stepSize = ?
xb=lowestVoltage;
do{i=0;

    do{j=0;
        p=1;
        while(j<i){
            p=p*(xb-x[j]);
            j++;
        }
        yb=yb+p*f[i][i];
        i++;
    }while(i<=n/2);
    doc << xb << " " << yb << endl;
    xb=xb+stepSize;
    yb=0;
}while(xb<=highestVoltage);
doc.close();
return 0;
}
