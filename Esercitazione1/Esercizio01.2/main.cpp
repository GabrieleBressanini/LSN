/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <fstream>

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  
    int n = 10000;      //numero di "variabili somma" che voglio ottenere per ogni casistica
    int N;              //numero di variabili casuali che sommo per ottenere la variabile somma
    double sum_uni, sum_exp, sum_lor;       //variabili somma relative alle diverse distribuzioni
    
    sum_uni = 0;
    sum_exp = 0;
    sum_lor = 0;
    
    N=1;
    
    ofstream out_uni;  //output dado uniforme, N=1
    out_uni.open("uniforme1.txt");
    ofstream out_exp;  //output dado esponenziale, N=1
    out_exp.open("esponenziale1.txt");
    ofstream out_lor;  //output dado lorentziano, N=1
    out_lor.open("lorentziana1.txt");
    
    for(int i=0; i<n; i++){
        
        sum_uni = 0;
        sum_exp = 0;
        sum_lor = 0;
        
        for(int j=0; j<N; j++){
            sum_uni += rnd.Rannyu();
            sum_exp += rnd.Exp(1.);
            sum_lor += rnd.Lorentz(1.,0.);
        }
        
        out_uni << sum_uni/N << endl;
        out_exp << sum_exp/N << endl;
        out_lor << sum_lor/N << endl;
    }
    
    out_uni.close();
    out_uni.clear();
    out_exp.close();
    out_exp.clear();
    out_lor.close();
    out_lor.clear();
    
    
    N=2;
    
    out_uni.open("uniforme2.txt");
    out_exp.open("esponenziale2.txt");
    out_lor.open("lorentziana2.txt");
    
    for(int i=0; i<n; i++){
        
        sum_uni = 0;
        sum_exp = 0;
        sum_lor = 0;
        
        for(int j=0; j<N; j++){
            sum_uni += rnd.Rannyu();
            sum_exp += rnd.Exp(1.);
            sum_lor += rnd.Lorentz(1.,0.);
        }
        
        out_uni << sum_uni/N << endl;
        out_exp << sum_exp/N << endl;
        out_lor << sum_lor/N << endl;
    }
    
    out_uni.close();
    out_uni.clear();
    out_exp.close();
    out_exp.clear();
    out_lor.close();
    out_lor.clear();
    
    
    N=10;
    
    out_uni.open("uniforme10.txt");
    out_exp.open("esponenziale10.txt");
    out_lor.open("lorentziana10.txt");
    
    for(int i=0; i<n; i++){
        
        sum_uni = 0;
        sum_exp = 0;
        sum_lor = 0;
        
        for(int j=0; j<N; j++){
            sum_uni += rnd.Rannyu();
            sum_exp += rnd.Exp(1.);
            sum_lor += rnd.Lorentz(1.,0.);
        }
        
        out_uni << sum_uni/N << endl;
        out_exp << sum_exp/N << endl;
        out_lor << sum_lor/N << endl;
    }
    
    out_uni.close();
    out_uni.clear();
    out_exp.close();
    out_exp.clear();
    out_lor.close();
    out_lor.clear();
    
    
    N=100;
    
    out_uni.open("uniforme100.txt");
    out_exp.open("esponenziale100.txt");
    out_lor.open("lorentziana100.txt");
    
    for(int i=0; i<n; i++){
        
        sum_uni = 0;
        sum_exp = 0;
        sum_lor = 0;
        
        for(int j=0; j<N; j++){
            sum_uni += rnd.Rannyu();
            sum_exp += rnd.Exp(1.);
            sum_lor += rnd.Lorentz(1.,0.);
        }
        
        out_uni << sum_uni/N << endl;
        out_exp << sum_exp/N << endl;
        out_lor << sum_lor/N << endl;
    }
    
    out_uni.close();
    out_uni.clear();
    out_exp.close();
    out_exp.clear();
    out_lor.close();
    out_lor.clear();
    
   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
