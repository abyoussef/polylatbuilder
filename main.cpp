#include <iostream>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2X.h>
#include <string>


#include "tools.h"

using namespace std;
using namespace NTL;

int main(){
    long w,m,s,iter;

    string answer,yon; 
    cout << "Enter the precision :" << "\n"; 
    cin >>w ;  
    cout << "Enter the power m (such as the number of points is 2^m):" << "\n"; 
    cin >>m ;  
    cout << "Enter the dimension :" << "\n"; 
    cin >>s ;  
    cout << "Enter the iteration :" << "\n"; 
    cin >>iter ;
    cout << "Enter the search type : "<< "\n"; 
    cout << "(CBC for component by component search, Korobov for Korobov PLR, Random for random search)"<< "\n";
    cin >>  answer;  
    cout << "Do you want to output points ? (y/n)"<<"\n";      
    cin >> yon;
    Vec<vec_RR> plr;
    plr = tools::gen(w,m,s,iter,answer);
    tools::outputPLR(plr,w,yon) ;
}
