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

    //PARTE UNO: Random Walk in spazio discreto.
    //NB: il punto di partenza di tutti i cammini sarà sempre l'origine
    
    ofstream out;
    out.open("RandomWalkDiscreto.txt");
    
    double a = 1.;          //lato dei cubetti in cui divido lo spazio
    int N = 10000;          //numero di iterazioni del random walk
    int M = 100;            //numero di step compiuti ad ogni iterazione
    double pos[6];          //vettore che contiene le "coordinate" (su,giu,destra,sinistra,avanti,indietro)
    
    //matrice che contiene la distanza del walker dall'origine all'iterazione i - primo indice - , step j - secondo indice -
    double distanza[N][M];
    double passo[6];        //vettore che contiene le 6 possibili direzioni di moto prima indicate
    
    for(int i=0; i<6; i++){
        passo[i] = a;
        pos[i] = 0;
    }
    
    for(int i=0; i<N; i++){
        for(int j=0; j<M; j++){
            distanza[i][j] = 0;
        }
    }
    
    int index;      //indice casuale tra 0 e 5 - probabilità uniforme - che mi dirà la direzione in cui si muove il walker
    
    for(int i=0; i<N; i++){                                         //ciclo sulle iterazioni del RW
        pos[0] = pos[1] = pos[2] = pos[3] = pos[4] = pos[5] = 0;    //posiziono il walker nell'origine all'inizio di un nuovo RW
        for(int j=0; j<M; j++){                                     //ciclo sugli step del singolo RW
            
            //in questo modo genero un numero casuale uniformemente distribuito tra 0 ... 5
            index = int(rnd.Rannyu(0.,6.));
            
            pos[index] += passo[index];
            
            //calcolo la distanza: NB pos[0]-pos[1] = "destra" - "sinistra" = coordinata x. Gli altri sono analoghi
            distanza[i][j] = sqrt(pow(pos[0]-pos[1],2) + pow(pos[2]-pos[3],2) + pow(pos[4]-pos[5],2));
        }
    }
    
    //ora ho, ad ogni step, N=10000 valori della distanza. Faccio la media a blocchi di M=100 per ottenere N/M medie.
    
    double media[N/M][M];       //matrice in cui salvo le medie a blocchi e lo step a cui fanno riferimento - secondo indice -
    
    for(int i=0; i<(N/M); i++){
        for(int j=0; j<M; j++){
            media[i][j] = 0;
        }
    }
    
    for(int k=0; k<M; k++){                             //ciclo sugli step
        for(int i=0; i<(N/M); i++){                     //ciclo sul numero di medie che otterrò per ogni step
            for(int j=0; j<M; j++){                     //ciclo sul numero di distanze che uso per calcolare la singola media
                media[i][k] += distanza[i*M + j][k];    //i*M + j è fatto per far si che la somma continui all'elemento successivo al termine del ciclo
            }
            media[i][k] /= M;
        }
    }
    
    //Calcolo quindi la media delle medie - per ogni step - e la relativa deviazione standard della media
    
    double sum[M], su2[M];        //variabili somma e somma dei quadrati che uso per calcolare medie ed errori
    double err[M];
    
    for(int i=0; i<M; i++){
        sum[i] = 0;
        su2[i] = 0;
        err[i] = 0;
    }
    
    for(int i=0; i<M; i++){             //ciclo sugli step
        for(int k=0; k<(N/M); k++){     //ciclo sul numero di medie che ho ad ogni step
            sum[i] += media[k][i];
            su2[i] += media[k][i] * media[k][i];
        }
        sum[i] /= (N/M);
        su2[i] /= (N/M);
        err[i] = error(sum[i],su2[i],N/M - 1);
        out << i <<"   "<< sum[i] <<"   "<< err[i] << endl;
    }
    
    //OSS: gli errori mi escono dell'ordine dei centesimi, e sono pertanto troppo piccoli per essere visti sul grafico
    
    out.close();
    out.clear();
    
    //SECONDA PARTE: come prima, ma ora la direzione di movimento è casuale (in 3D ovviamente)
    
    out.open("RandomWalkContinuo.txt");
    
    double r[3];            //questo è il mio nuovo vettore posizione
    double phi,theta;       //angolo polare e azimutale
    
    for(int i=0; i<3; i++){
        r[i] = 0;
    }
    
    for(int i=0; i<N; i++){
        for(int j=0; j<M; j++){
            distanza[i][j] = 0;
        }
    }
    
    for(int i=0; i<N; i++){
        r[0] = r[1] = r[2] = 0;    //posiziono il walker nell'origine
        for(int j=0; j<M; j++){
            phi = 2 * M_PI * rnd.Rannyu();          //angolo polare generato uniformemente tra 0 e 2pi
            theta = acos(1 - 2 * rnd.Rannyu());     //angolo azimutale generato in modo da garantireuna direzione casuale uniforme - Lecture 1
            r[0] += sin(theta)*cos(phi);            //coordinate sferiche
            r[1] += sin(theta)*sin(phi);
            r[0] += cos(theta);
            distanza[i][j] = sqrt(pow(r[0],2) + pow(r[1],2) + pow(r[2],2));
        }
    }
    
    for(int i=0; i<(N/M); i++){
        for(int j=0; j<M; j++){
            media[i][j] = 0;
        }
    }

    for(int k=0; k<M; k++){                             //ciclo sugli step
        for(int i=0; i<(N/M); i++){                     //ciclo sul numero di medie che otterrò per ogni step
            for(int j=0; j<M; j++){                     //ciclo sul numero di distanze che uso per calcolare la media
                media[i][k] += distanza[i*M + j][k];    //i*M + j è fatto per far si che la somma continui all'elemento successivo al termine del ciclo
            }
            media[i][k] /= M;
        }
    }
    
    for(int i=0; i<M; i++){
        sum[i] = 0;
        su2[i] = 0;
        err[i] = 0;
    }
    
    for(int i=0; i<M; i++){
        for(int k=0; k<(N/M); k++){
            sum[i] += media[k][i];
            su2[i] += media[k][i] * media[k][i];
        }
        sum[i] /= (N/M);
        su2[i] /= (N/M);
        err[i] = error(sum[i],su2[i],N/M - 1);
        out << i <<"   "<< sum[i] <<"   "<< err[i] << endl;
    }
    
    //OSS: gli errori, ancora una volta, mi escono dell'ordine dei centesimi, e sono pertanto troppo piccoli per essere visti
    
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
