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

    ofstream out;
    out.open("results.txt");
    
    int M = 100000;        //numero di lanci
    int N = 100;           //numero di blocchi
    double L = 0.8;       //lunghezza dell'ago
    double d = 1;          //distanza tra le linee parallele
    double x;              //distanza tra il centro dell'ago e la linea più vicina (compreso tra 0 e d/2)
    double px,py;          //coordinate per generazione angolo tra 0 e pi/2 (distribuzione uniforme)
    int Nhit=0,Ntot=0;         //contatori per il numero di hit ed il numero totale di aghi lanciati
    double p[N];                  //vettore che contiene le N stime di pi
    
    for(int i=0; i<N; i++){             //ciclo sui blocchi
        Nhit=0;
        Ntot=0;
        for(int k=0; k<M/N; k++){       //ciclo sul numero di lanci per blocco
            x = d/2 * rnd.Rannyu();     //estraggo uniformemente la distanza del centro dell'ago dalla linea più vicina
            do{
                px = rnd.Rannyu();      //due punti uniformemente tra 0 ed 1. Se cadono all'interno della circonferenza
                py = rnd.Rannyu();      //unitaria li uso per calcolare il seno dell'angolo (in questo modo la distrubuzione è uniforme)
            }while(px*px + py*py > 1);
           
            if(x <= L/2 * py/(sqrt(px*px + py*py))){       //si noti che py/sqrt(px*px + py*py) è il seno di un angolo distribuito uniformemente nel primo quadrante
                Nhit++;
            }
            Ntot++;
        }
        p[i] = 2 * L * Ntot / (d * Nhit);
    }
    
    double sum_prog[N], su2_prog[N];        //somme progressive
    double err_prog[N];
    
    for(int i=0; i<N; i++){
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    for(int i=0; i<N; i++){
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
        for(int k=0; k<i+1; k++){
            sum_prog[i] += p[k];
            su2_prog[i] += p[k] * p[k];
        }
        sum_prog[i] = sum_prog[i] / (i+1);
        su2_prog[i] = su2_prog[i] / (i+1);
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
