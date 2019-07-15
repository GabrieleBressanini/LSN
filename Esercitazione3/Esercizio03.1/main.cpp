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
#include <fstream>
#include <math.h>

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

    double S0 = 100.;        //asset price at t=0;
    double T = 1.;           //delivery time
    double k = 100.;         //strike price
    double r = 0.1;         //risk-free interest rate
    double sigma = 0.25;    //volatility
    
    int N = 100;            //number of blocks
    int M = 10000;          //number of "throws" in total
    
    double w = 0.;               //stochastic process
    
    //Task: compute at time t=0 via Monte Carlo the European call-option price C[S(0),0] and put-option price P[S(0),0]
    
    //First Method: direct sampling of the final asset price
    
    double ST;              //asset price at t=T, initial value S0
    double call[N];         //vettore che conterrà le 100 stime della call (una da ogni blocco)
    double put[N];          //vettore che conterrà le 100 stime della put (una da ogni blocco)
    double call2[N];        //salvo i valori al quadrato delle medie per il futuro calcolo dell'errore
    double put2[N];         //salvo i valori al quadrato delle media per il futuro calcolo dell'errore
    
    double c=0,p=0;             //variabili di comodo per il calcolo della media su ogni blocco delle call e put rispettivamente
   
    double sum_prog[N],su2_prog[N],err_prog[N];
    
    for(int i=0; i<N; i++){
        call[i] = put[i] = put2[i] = call2[i] = 0;
    }
    
    for(int j=0; j<N; j++){
        c = 0;
        p = 0;
        for(int i=0; i<(M/N); i++){
            w = rnd.Gauss(0.,1.);
            ST = S0 * exp((r - 0.5 * sigma * sigma) * T + sigma * w);
            c += exp(-r*T) * max(0.,ST-k);
            p += exp(-r*T) * max(0.,k-ST);
        }
        
        call[j] = c / double(M/N);
        put[j] = p / double(M/N);
        put2[j] = put[j] * put[j];
        call2[j] = call[j] * call[j];
    }
    
    //calcolo dell'errore "progressivo"

    for(int i=0; i<N; i++){
        sum_prog[i] = su2_prog[i] = err_prog[i] = 0;
    }
    
    ofstream out;
    out.open("CallDirect.txt");
    
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            sum_prog[i] = sum_prog[i] + call[j];
            su2_prog[i] = su2_prog[i] + call2[j];
        }
        sum_prog[i] = sum_prog[i]/double(i+1);
        su2_prog[i] = su2_prog[i]/double(i+1);
        err_prog[i] = error(sum_prog[i],su2_prog[i],i);
        out << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
    }
    
    out.close();
    out.clear();
    
    out.open("PutDirect.txt");
    
    for(int i=0; i<N; i++){
        sum_prog[i] = su2_prog[i] = err_prog[i] = 0;
    }
    
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            sum_prog[i] = sum_prog[i] + put[j];
            su2_prog[i] = su2_prog[i] + put2[j];
        }
        sum_prog[i] = sum_prog[i]/double(i+1);
        su2_prog[i] = su2_prog[i]/double(i+1);
        err_prog[i] = error(sum_prog[i],su2_prog[i],i);
        out << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
    }
    
    out.close();
    out.clear();
  
    //Second Method: discretized sampling of the final asset price
    
    for(int i=0; i<N; i++){
        call[i] = put[i] = put2[i] = call2[i] = 0;
    }
    
    for(int j=0; j<N; j++){
        c = 0;
        p = 0;
        for(int i=0; i<(M/N); i++){
            ST = S0;
            for(int m=0; m<100; m++){   //divido l'intervallo 0-T in 100 parti
                w = rnd.Gauss(0.,1.);
                ST = ST * exp((r - 0.5 * sigma * sigma) * 0.01 + sigma * w * sqrt(0.01));
            }
            
            c += exp(-r*T) * max(0.,ST-k);
            p += exp(-r*T) * max(0.,k-ST);
        }
        
        call[j] = c / double(M/N);
        put[j] = p / double(M/N);
        put2[j] = put[j] * put[j];
        call2[j] = call[j] * call[j];
    }
    
    //calcolo dell'errore "progressivo"
    
    for(int i=0; i<N; i++){
        sum_prog[i] = su2_prog[i] = err_prog[i] = 0;
    }
    
    out.open("CallDiscretized.txt");
    
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            sum_prog[i] = sum_prog[i] + call[j];
            su2_prog[i] = su2_prog[i] + call2[j];
        }
        sum_prog[i] = sum_prog[i]/double(i+1);
        su2_prog[i] = su2_prog[i]/double(i+1);
        err_prog[i] = error(sum_prog[i],su2_prog[i],i);
        out << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
    }
    
    out.close();
    out.clear();
    
    out.open("PutDiscretized.txt");
    
    for(int i=0; i<N; i++){
        sum_prog[i] = su2_prog[i] = err_prog[i] = 0;
    }
    
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            sum_prog[i] = sum_prog[i] + put[j];
            su2_prog[i] = su2_prog[i] + put2[j];
        }
        sum_prog[i] = sum_prog[i]/double(i+1);
        su2_prog[i] = su2_prog[i]/double(i+1);
        err_prog[i] = error(sum_prog[i],su2_prog[i],i);
        out << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
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
