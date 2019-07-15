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
#include <math.h>
#include <fstream>

using namespace std;

double error(double,double,int);
 
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
    
    // PARTE 1
    
    int M = 100000;     //definisco il numero totale di lanci
    int N = 100;        //definisco il numero di blocchi in cui soddivido i lanci
    int L = M/N;        //lunghezza del blocco
    double sum;
    double ave[N],av2[N];
    double sum_prog[N],su2_prog[N],err_prog[N]; //prog sta per "progressivo"
    
    for(int i=0; i<N; i++){         //NB: sempre inizializzare i vettori in questo modo!
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    ofstream out;
    out.open("results.txt");
    
    for(int i=0; i<N; i++){
        sum = 0;
        for(int j=0; j<L; j++){
            sum = sum + rnd.Rannyu();
        }
        ave[i] = sum / L;
        av2[i] = ave[i] * ave[i];
    }
    
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            sum_prog[i] = sum_prog[i] + ave[j];
            su2_prog[i] = su2_prog[i] + av2[j];
        }
        sum_prog[i] = sum_prog[i]/(i+1);
        su2_prog[i] = su2_prog[i]/(i+1);
        err_prog[i] = error(sum_prog[i],su2_prog[i],i);
        out << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
    }
    
    out.close();
    out.clear();
    
    // PARTE 2
    
    out.open("results2.txt");
    
    for(int i=0; i<N; i++){
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    double x = 0;
    
    for(int i=0; i<N; i++){
        sum = 0;
        for(int j=0; j<L; j++){
            x = rnd.Rannyu();
            sum = sum + (x - 0.5)*(x - 0.5);
        }
        ave[i] = sum / L;
        av2[i] = ave[i] * ave[i];
    }
    
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            sum_prog[i] = sum_prog[i] + ave[j];
            su2_prog[i] = su2_prog[i] + av2[j];
        }
        sum_prog[i] = sum_prog[i]/(i+1);
        su2_prog[i] = su2_prog[i]/(i+1);
        err_prog[i] = error(sum_prog[i],su2_prog[i],i);
        out << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
    }
    
    
    out.close();
    out.clear();
    
    //PARTE 3
    
    out.open("results3.txt");
    
    M = 100;                //numero di intervalli in cui divido l'intervallo [0,1]
    int n = 10000;          //numeri casuali distribuiti uniformemente generati ad ogni "blocco"
    double posizione[M];    //i-esimo elemento indica il riempimento dell'i-esimo intervallo
    int pos = 0;
    
    for(int i=0; i<M; i++){
        posizione[i]=0;
    }
    
    for(int k=0; k<100; k++){
        sum = 0;
        
        for(int i=0; i<M; i++){
            posizione[i]=0;
        }
        
        for(int i=0; i<n; i++){             //riempimento del vettore posizione
            pos = int(rnd.Rannyu()*100);    //si noti il modo in cui questo avviene
            posizione[pos] += 1;
        }

        for(int j=0; j<M; j++){             //calcolo del X^2
            sum = sum + (posizione[j]-n/M)*(posizione[j]-n/M);
        }
        sum = sum /(n/M);
        
        out << k << "  " << sum << endl;
    }
    
    
    out.close();
    out.clear();
    
   rnd.SaveSeed();
   return 0;
}


double error(double sum_prog, double su2_prog, int n){
    double err;
    if(n==0){
        return 0;
    }
    err = sqrt(1./(n) * (su2_prog - sum_prog * sum_prog));
    return err;
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
