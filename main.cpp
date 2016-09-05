#include <iostream>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2X.h>

#include "tools.h"

using namespace std;
using namespace NTL;

int main(){
    long w,m,s,iter;
    cout << "Enter the precision :" << "\n"; 
    cin >>w ;  
    cout << "Enter the power m (such as the number of points is 2^m):" << "\n"; 
    cin >>m ;  
    cout << "Enter the dimension :" << "\n"; 
    cin >>s ;  
    cout << "Enter the iteration :" << "\n"; 
    cin >>iter ;  
    Vec<vec_RR> plr;
    plr = tools::genCBC(w,m,s,iter);
    tools::outputPLR(plr,w) ;
}