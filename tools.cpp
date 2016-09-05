#include "tools.h"
#include <iostream>

// NTL libraries used :
#include <NTL/GF2X.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/vec_RR.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/vector.h>
using namespace std;
using namespace NTL;

vec_GF2     tools::coeffVector( const GF2X& q, long n){
    vec_GF2 v;
    v.SetLength(n);
    for(int i=1; i <= n; ++i ){
        v(i) = coeff(q,n-i);
    }

    return v;
};

mat_GF2     tools::coeffMatrix( const GF2X& p){
    long n;
    n = deg(p);
    mat_GF2 C;
    C = ident_mat_GF2(n);
    for(int i = 2; i <= n ; ++i ){
        for(int j = 1; j < i ; ++j){
            C(i,j) = coeff(p,n - i + j);
        }
    }
    return C;
};

vec_RR      tools::coord(const mat_GF2& C){
    long r,c;
    vec_RR res ;
    RR b;

    r = C.NumRows();
    c = C.NumCols();
    res.SetLength(r);

    for(long i = 0; i < r ; ++i ){
        b=0.5;
        for(long j = 0; j < c ; ++j ){
            res[i] += b * rep(C[i][j]);
            b /= 2;
        }
    }

    return res;
};

vec_GF2X    tools::genk( long m){
    long N;
    N = pow(2,m);

    vec_GF2X k;
    k.SetLength(N);
    k[0] = GF2X(INIT_MONO, 0, 0);
    k[1] = GF2X(INIT_MONO, 0, 1);
    long j = 2;
    GF2X p;
    for(int i = 1; i < m  ; ++i  ){
        p = GF2X(INIT_MONO,i);
        for( long l = j; l< 2*j ; ++l ){
           k[l]= k[l-j]+p;
        }
        j*=2;
    }
    p.kill();
    return k;
};



void        tools::productMod(vec_GF2X& x, const GF2X& q, const vec_GF2X& v, const GF2X& f){
    GF2XModulus F(f);
    x.SetLength(v.length());

    for (long i = 0; i < v.length(); i++){
            MulMod(x[i], q, v[i], F);
    }
};

vec_GF2     tools::laurentVec(const GF2X& q, const GF2X& p, long w ){
    long m;
    mat_GF2 C;
    vec_GF2 b;
    GF2 d;
    vec_GF2 x;

    m  =  deg(p);
    C  =  coeffMatrix(p);
    b  =  coeffVector(q,deg(p));

    solve( d, C, x , b);

    if( w <= m ){
        return x;
    }else{
        vec_GF2 y;
        y.SetLength(w);
        for(int k=0; k < x.length(); ++k){
            y[k]=x[k];
        }
        for( int i = m + 1 ; i < w ; ++i){
            for(int j = 1; j <= m ; ++j){
                y[i] += y[i-j] * coeff(p,m-j);
            }
        }
        return y;
    }
};

mat_GF2     tools::laurentMat(const vec_GF2X& k, const GF2X& p, long w ){
    long m,N;
    mat_GF2 x;

    m  =  deg(p);
    N = k.length() ;
    if( w <= m ){
        x.SetDims( N, m);
    }else{
        x.SetDims( N, w);
    }
    for(long i = 0; i < N ; ++i){
        x[i] = laurentVec( k[i] , p,  w );
    }
    return x;
};

RR          tools::wal( const RR& x){
    RR deux (2);
    return  12 * (1/6 - pow( deux , floor(log(x)/log(deux)) - 1 ));
};

RR          tools::walVec(const vec_RR& v){
    RR res(1);
    for(long i = 0; i < v.length(); ++i ){
        if(v[i] > 0)
               	res *= wal(v[i]);
    }
    return res;
};



vec_GF2X       	tools::extract(const vec_GF2X& v , long n) {
       	vec_GF2X res;
       	long l(n);
       	if( l > v.length() ) {
       		l = v.length();
       	}
       	res.SetLength(l);
       	for(long i=0; i< l; ++ i){
       		res[i] = v[i];
       	}
       	return res;
}

//void 		tools::kill(const vec_GF2X& v){
//     	for(long i=0; i<v.length() ; ++i ){
//     		v[i].kill();
//     	}
//}

// Figures of Merit
RR          tools::ODWs( const Vec<vec_RR>& plr){
    RR res(0);
    long N;

    N = plr.length();
    for( long i = 1; i < N ; ++i ){
        res += walVec(plr[i]);
    }
    return abs(res/N);
};

RR          tools::ODW1( const Vec<vec_RR>& plr, const vec_RR& w){
    RR res(0);
    RR temp(0);
    long N,s;

    N = plr.length();
    s = w.length();

    for(long l = 0; l < s ; ++l ){
        for(long i = 1 ; i < N ; ++i ){
            temp += wal(plr[i][l]);
        }
        res += temp * w[l];
        temp = 0;
    }
    return res/N;
};


// Generation algorithms of PLR

Vec<vec_RR> tools::genplr( const vec_GF2X& z, const vec_GF2X& k, const GF2X& p, long w){
    long N,s;
    mat_GF2 U;
    Vec<vec_RR> y;
    s = z.length();
    N = k.length();

    y.SetLength(N);


    for(long i=0; i < N; ++i ){
        y[i].SetLength(s);
    }

    vec_GF2X x;
    mat_GF2 C;
    vec_RR coords;
    for(long l = 0; l < s ; ++l){
        productMod( x, z[l] , k, p);
        C = laurentMat( x, p,  w ) ;
        coords = coord(C);
        for(long i=0; i < N;++i){
            y[i][l] = coords[i];
        }
    }

    for (long i = 0; i < x.length(); i++){
        x[i].kill();
    }
    coords.kill();
    C.kill();

    return y;
};


Vec<vec_RR> tools::genRandom( long w, long m, long s , long iter){
    vec_GF2X z,x;
    Vec<vec_RR> plr;
    vec_GF2X k = tools::genk(m);
    GF2X p = BuildSparseIrred_GF2X(m);
    z.SetLength(s);
    RR best(0);
    RR odws(0);

    for(long i = 0; i < iter ; ++i){
        for(long l = 0; l<s ; ++l){
            z[l] = random_GF2X( m );
            //temp = random_GF2X( m );
            //div(z[l], temp, GCD( temp, p) );
        }
        plr = genplr(  z,  k,  p, w);
        odws = ODWs(plr);
        // cout << "Step " << i << "  "<< odws << "\n";
        if ( (odws <=  best) || (i==0) ) {
            best =odws;
            x=z;
        }
        cout << best << "\n";
    }
    cout << "Polynomial P is " << p << "\n";
    cout << "Generating vector z in "<<s<<" dimensions :" << "\n";
    for(long l =0; l< s ; ++l ) {
        cout <<"        Poly["<<l<<"] = "<<z[l]<<"\n";
    }
    return genplr(x,k,p,w);
};

Vec<vec_RR> tools::genKorobov( long w, long m, long s , long iter){
    vec_GF2X z,x;
    Vec<vec_RR> plr;
    vec_GF2X k = tools::genk(m);
    GF2X p = BuildSparseIrred_GF2X(m);
    z.SetLength(s);
    RR best(0);
    RR odws(0);
    for(long i = 0; i < iter ; ++i){
        set(z[0]);
        z[1] = random_GF2X( m );
        for(long l = 2; l<s ; ++l){
            z[l] = z[l-1] * z[1];
            //temp = random_GF2X( m );
            //div(z[l], temp, GCD( temp, p) );
        }
        plr = genplr(  z,  k,  p, w);
        odws = ODWs(plr);
        //cout << "Step " << i << "  "<< odws << "\n";
        if ( (odws <=  best) || (i==0) ) {
            best =odws;
            x=z;
        }
    }
    //cout << best << "\n";
    cout << "Polynomial P is " << p << "\n";
    cout << "Generating vector z in "<<s<<" dimensions :" << "\n";
    for(long l =0; l< s ; ++l ) {
        cout <<"        Poly["<<l<<"] = "<<x[l]<<"\n";
    }
    return genplr(x,k,p,w);
};

Vec<vec_RR> tools::genCBC( long w, long m, long s , long iter){
    vec_GF2X z;
    GF2X x;
    Vec<vec_RR> plr;
    vec_GF2X k = tools::genk(m);
    GF2X p = BuildSparseIrred_GF2X(m);
    RR best(0);
    RR odws(0);
    GF2X fin;
    z.SetLength(s);
    for(long l=0; l<s;++l){
       	z[l] = 0 ;
    }
    for(long l = 0; l < s ; ++l){
       	best = 0;
       	odws = 0;
        for(long i = 0 ; i < iter ; ++i ){
       		x = random_GF2X(m);
       		z[l] = x ;
       		plr = genplr(extract(z,l+1),k,p,w);
       		odws = ODWs(plr);
       		if ( (odws <=  best) || (i==0) ) {
               		best = odws;
               		fin = x ;
               	}
       	}
        z[l] = fin ;
        // cout << "Step " << l << "  "<< odws << "\n";
    }
    cout << "Polynomial P is " << p << "\n";
    cout << "Generating vector z in "<<s<<" dimensions :" << "\n";
    for(long l =0; l< s ; ++l ) {
       	cout <<"       	Poly["<<l<<"] = "<<z[l]<<"\n";
    }
    return genplr(z,k,p,w);
};

void        tools::outputPLR( const Vec<vec_RR>& plr, long w){
    long N,s;

    RR::SetOutputPrecision(w);
    N = plr.length();
    s = plr[0].length();

    for(long l = 0; l< s ; ++l){
        for(long i = 0; i<plr.length() ; ++i){
            cout << plr[i][l] << "  ";
        }
        cout << "\n";
    }
};
