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
    
    //PRIMA PARTE: Calcolo dell'integrale con sampling uniforme (metodo della media)
    
    ofstream out;
    out.open("SamplingUniforme.txt");
    
    int M = 100000;         //numero di lanci
    int N = 100;            //numero di blocchi
    double integrale[N];    //vettore che contiene le N stime dell'integrale
    double x;               //conterrà la variabile indipendente della funzione
    
    for(int i=0; i<N; i++){
        integrale[i]=0;
    }
    
    
    for(int i=0; i<N; i++){
        for(int k=0; k<(M/N); k++){
            x = rnd.Rannyu();                               //sampling uniforme
            integrale[i] += M_PI/2 * cos(M_PI/2 * x);       //valuto la funzione nel punto estratto
        }
        integrale[i] /= (M/N);                              //valuto l'integrale con il teorema della media
    }
    
    double sum_prog[N], su2_prog[N];        //somme progressive di medie ed errori
    double err_prog[N];
    
    for(int i=0; i<N; i++){
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    
    //Calcolo della media e dell'errore "a blocchi"
    for(int i=0; i<N; i++){
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
        for(int k=0; k<i+1; k++){
            sum_prog[i] += integrale[k];
            su2_prog[i] += integrale[k] * integrale[k];
        }
        sum_prog[i] = sum_prog[i] / (i+1);
        su2_prog[i] = su2_prog[i] / (i+1);
        err_prog[i] = error(sum_prog[i],su2_prog[i],i);
        out << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
    }
    
    out.close();
    out.clear();
    
    //SECONDA PARTE: IMPORTANCE SAMPLING
    
    //Utilizzo come distribuzione di probabilità (occhio che la funzione deve essere non negativa e normalizzata!) p(x)=2(1-x)
    
    out.open("ImportanceSampling.txt");
    
    double y;       //sarà la variabile che genero uniformemente tra 0 ed 1
    
    for(int i=0; i<N; i++){
        integrale[i]=0;
    }
    
    for(int i=0; i<N; i++){
        for(int k=0; k<(M/N); k++){
            y = rnd.Rannyu();           //estraggo uniformemente tra 0 ed 1
            x = 1 - sqrt(1-y);          //estraggo secondo la distribuzione pa (Metodo della Trasformata : tutto analitico!)
            integrale[i] += M_PI/2 * cos(M_PI/2 * x) / (2*(1-x));       //somma della funzione valutata in x, pesata sulla distribuzione!
        }
        integrale[i] /= (M/N);
        cout << integrale[i] << endl;
    }
    
    for(int i=0; i<N; i++){
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    
    //Calcolo della media e dell'errore "a blocchi"
    for(int i=0; i<N; i++){
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
        for(int k=0; k<i+1; k++){
            sum_prog[i] += integrale[k];
            su2_prog[i] += integrale[k] * integrale[k];
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
