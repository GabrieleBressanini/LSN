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

    // *** Le lunghezze sono espresse in unità del Raggio di Bohr ***
    
    //Distribuzione Uniforme
    
    ofstream out,graph;
    out.open("1s.txt");
    graph.open("graph1s.txt");
    
    double sum = 0.;
    double r = 0.;      //raggio
    double x,y,z;
    int M = 1000000;    //numero di lanci
    int N = 100;        //numero di blocchi
    double step = 1.;   //step del campionamento (in unità di a0). Provo lo step che mi assicura una p=50%
    
    //per lo stato 1s pongo il punto iniziale nell'origine degli assi
    x = y = z = 0;
    double xnew,ynew,znew;
    xnew = ynew = znew = 0;
    double A = 0,rnew = 0;
    
    double ave[N],av2[N];
    double sum_prog[N],su2_prog[N],err_prog[N]; //prog sta per "progressivo"
    
    for(int i=0; i<N; i++){         //NB: sempre inizializzare i vettori in questo modo!
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }

    for(int j=0; j<N; j++){
        sum = 0;
        for(int i=0; i<(M/N); i++){
            r = sqrt(x*x + y*y + z*z);
            xnew = x + step * rnd.Rannyu(-1.,1.);
            ynew = y + step * rnd.Rannyu(-1.,1.);
            znew = z + step * rnd.Rannyu(-1.,1.);
            rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
            A = min(1.,exp(-2 * (rnew - r)));
            if(rnd.Rannyu() <= A){
                sum += rnew;
                x = xnew;
                y = ynew;
                z = znew;
                graph << x <<" "<< y <<" "<< z << endl;
            }else{
                sum += r;
            }
            ave[j] = sum/(M/N);
            av2[j] = ave[j]*ave[j];
        }
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
    graph.close();
    graph.clear();
    
    out.open("2p.txt");
    graph.open("graph2p.txt");
    
    for(int i=0; i<N; i++){         //NB: sempre inizializzare i vettori in questo modo!
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    x = y = 0;
    z = 1;
    step = 1.2;
    
    for(int j=0; j<N; j++){
        sum = 0;
        for(int i=0; i<(M/N); i++){
            r = sqrt(x*x + y*y + z*z);
            xnew = x + step * rnd.Rannyu(-1.,1.);
            ynew = y + step * rnd.Rannyu(-1.,1.);
            znew = z + step * rnd.Rannyu(-1.,1.);
            rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
            A = min(1.,pow((rnew/r),2) * exp(-(rnew - r)) * pow(znew/z * r/rnew,2));
            if(rnd.Rannyu() <= A){
                sum += rnew;
                x = xnew;
                y = ynew;
                z = znew;
                graph << x <<" "<< y <<" "<< z << endl;
            }else{
                sum += r;
            }
            ave[j] = sum/(M/N);
            av2[j] = ave[j]*ave[j];
        }
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
    graph.close();
    graph.clear();
    
    //Distribuzione Gaussiana
    
    out.open("1sG.txt");
    
    //per lo stato 1s pongo il punto iniziale nell'origine degli assi
    x = y = z = 0;
    step = 0.75;
    
    for(int i=0; i<N; i++){         //NB: sempre inizializzare i vettori in questo modo!
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    for(int j=0; j<N; j++){
        sum = 0;
        for(int i=0; i<(M/N); i++){
            r = sqrt(x*x + y*y + z*z);
            xnew = x + step * rnd.Gauss(0.,1.);
            ynew = y + step * rnd.Gauss(0.,1.);
            znew = z + step * rnd.Gauss(0.,1.);
            rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
            A = min(1.,exp(-2 * (rnew - r)));
            if(rnd.Rannyu() <= A){
                sum += rnew;
                x = xnew;
                y = ynew;
                z = znew;
            }else{
                sum += r;
            }
            ave[j] = sum/(M/N);
            av2[j] = ave[j]*ave[j];
        }
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

    out.open("2pG.txt");

    for(int i=0; i<N; i++){         //NB: sempre inizializzare i vettori in questo modo!
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    x = y = 0;
    z = 1;
    step = 1.88;
    
    for(int j=0; j<N; j++){
        sum = 0;
        for(int i=0; i<(M/N); i++){
            r = sqrt(x*x + y*y + z*z);
            xnew = x + step * rnd.Gauss(0.,1.);
            ynew = y + step * rnd.Gauss(0.,1.);
            znew = z + step * rnd.Gauss(0.,1.);
            rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
            A = min(1.,pow((rnew/r),2) * exp(-(rnew - r)) * pow(znew/z * r/rnew,2));
            if(rnd.Rannyu() <= A){
                sum += rnew;
                x = xnew;
                y = ynew;
                z = znew;
            }else{
                sum += r;
            }
            ave[j] = sum/(M/N);
            av2[j] = ave[j]*ave[j];
        }
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
