#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>

// NTL libraries used :
#include <NTL/GF2X.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/vec_RR.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vector.h>

using namespace NTL;
using namespace std;

class tools {
public:
    static mat_GF2     coeffMatrix( const GF2X& p);
    static vec_GF2     coeffVector( const GF2X& q, long n);
    static vec_RR      coord(const mat_GF2& C);
    static vec_GF2X    genk( long m);
    static void        productMod(vec_GF2X& x, const GF2X& q, const vec_GF2X& v, const GF2X& f);
    static vec_GF2     laurentVec(const GF2X& q, const GF2X& p, long w );
    static mat_GF2     laurentMat(const vec_GF2X& k, const GF2X& p, long w );
    static vec_GF2X    extract(const vec_GF2X& v , long n );
    //static void        kill(const vec_GF2X& v);

    static RR          wal( const RR& x);
    static RR          walVec(const vec_RR& v);

    static RR          ODWs( const Vec<vec_RR>& plr);
    static RR          ODW1( const Vec<vec_RR>& plr, const vec_RR& w);

    static Vec<vec_RR> genplr( const vec_GF2X& z, const vec_GF2X& k, const GF2X& p, long w);
    static Vec<vec_RR> genRandom( long w, long m, long s, long iter );
    static Vec<vec_RR> genKorobov( long w, long m, long s, long iter );
    static Vec<vec_RR> genCBC( long w, long m, long s, long iter);

    static void outputPLR( const Vec<vec_RR>& plr, long w);
};

#endif