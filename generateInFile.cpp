#include<iostream>
#include <fstream>
#include<cstdlib>

using namespace std;

main()
{
    ofstream myfile;
    myfile.open ("query_15_1.txt");

    int n; //for number of data points
    int dim;

    cout << "enter no of transaction" << endl ;
    cin >> n;

    cout << "enter dimension" << endl;
    cin >> dim;
    myfile << dim << endl;
    myfile << n << endl;


    for ( int i = 1 ; i <= n ; i++ )
    {
        for(int j =1; j<=dim; j++)
        {
            double num = ((double) rand() / (RAND_MAX));
            myfile << num << ",";
        }
        myfile <<endl;
    }

    myfile.close();
    return 0;
}
