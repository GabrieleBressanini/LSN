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
double psi(double,double,double);  //funzione d'onda, prende x e i due parametri mu e sigma
double pot(double);                //potenziale, prende x
double psi2(double,double,double); //derivata seconda, prende x e i due parametri mu e sigma
 
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
    out.open("energia.txt");
    
    double sum = 0.;
    double sigma;
    double mu;
    double x;
    int M = 10000000;    //numero di lanci
    int N = 100;        //numero di blocchi
    double step = 2.2;   //step del campionamento. Provo lo step che mi assicura una p=50%
    
    //pongo il punto iniziale nell'origine degli assi
    x = 0;
    double xnew;
    xnew = 0;
    double A;   //accettazione del Metropolis
    
    double ave[N],av2[N];
    double sum_prog[N],su2_prog[N],err_prog[N]; //prog sta per "progressivo"
    
    for(int i=0; i<N; i++){         //NB: sempre inizializzare i vettori in questo modo!
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    
    int tentativi = 0;
    int accettati = 0;
   
    double minimo_energia = 999999.;
    double mu_opt, sigma_opt;   //valori ottimali dei due parametri, da stimare
    mu_opt = sigma_opt = 0;
    
    for(double mu = 0.75 ; mu <= 0.81 ; mu = mu + 0.01 ){
        for(double sigma = 0.55; sigma <= 0.65; sigma = sigma + 0.01){
            sum = 0;
            accettati = 0;
            tentativi = 0;
            for(int i=0; i<M; i++){
                xnew = x + step * rnd.Rannyu(-1.,1.);
                double rapporto;
                rapporto = pow(psi(xnew,mu,sigma),2.) / pow(psi(x,mu,sigma),2.) ;
                A = min(1.,rapporto);
                if(rnd.Rannyu() <= A){
                    x = xnew;
                    accettati ++;
                }
                sum += (-0.5 * psi2(x,mu,sigma) + pot(x)*psi(x,mu,sigma))/(psi(x,mu,sigma));
                tentativi++;
            }
            cout <<"Sigma:   "<< sigma <<"  mu:   "<< mu << endl;
            cout << "L'accettazione è:  " << (double)(accettati)/(double)(tentativi) << endl;
            sum /= M;
            if(sum < minimo_energia){
                minimo_energia = sum;
                mu_opt = mu;
                sigma_opt = sigma;
            }
        }
        
    }
    
    cout << "Miglior valore di mu: "<< mu_opt << endl;
    cout << "Miglior valore di sigma: "<< sigma_opt << endl;
    mu = mu_opt;
    sigma = sigma_opt;
    
    cout << "STIMA DELL'ENERGIA CON I PARAMETRI OTTIMALI E DATA BLOCKING" << endl << endl;
    
    tentativi = 0;
    accettati = 0;
    x = 0;
    int nbins = 100;
    int istogrammi[nbins][N];   //questa matrice contiene in ogni riga l'istogramma riferito ad un singolo blocco
    
    for(int i=0; i<nbins; i++){
        for(int j=0; j<N; j++){
            istogrammi[i][j] = 0;
        }
    }
    
    //per riempire gli istogrammi vado da -2.5 a 2.5 con 100 intervalli, ognuno largo 0.05
    
    for(int j=0; j<N; j++){
        sum = 0;
        for(int i=0; i<(M/N); i++){
            xnew = x + step * rnd.Rannyu(-1.,1.);
            double rapporto;
            rapporto = pow(psi(xnew,mu,sigma),2.) / pow(psi(x,mu,sigma),2.) ;
            A = min(1.,rapporto);
            if(rnd.Rannyu() <= A){
                x = xnew;
                accettati++;
            }
            sum += (-0.5 * psi2(x,mu,sigma) + pot(x)*psi(x,mu,sigma))/(psi(x,mu,sigma));
            tentativi++;
            
            int index; //indice che indicherà il bin corretto dell'istogramma
            index = int((x+2.5)/0.05);  //0.05 è la larghezza del bin
            istogrammi[index][j]++;
        }
        ave[j] = sum/(M/N);
        av2[j] = ave[j]*ave[j];
    }
    
    cout << "L'accettazione è:  " << (double)(accettati)/(double)(tentativi) << endl;
    
    //Data blocking per l'energia
    
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
    
    //Data blocking per l'istogramma
    
    out.open("istogramma.txt");
    
    for(int i=0; i<N; i++){
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }

    cout << "ESEGUO DATA BLOCKING ISTOGRAMMA" << endl << endl;
    
    for(int ibin=0; ibin<nbins; ibin++){
        for(int i=0; i<N; i++){
            sum_prog[i] = 0;
            su2_prog[i] = 0;
            err_prog[i] = 0;
        }
        for(int i=0; i<N; i++){
            for(int j=0; j<i+1; j++){
                sum_prog[i] = sum_prog[i] + istogrammi[ibin][j];
                su2_prog[i] = su2_prog[i] + istogrammi[ibin][j]*istogrammi[ibin][j];
            }
            
            sum_prog[i] = sum_prog[i]/(double)(i+1);
            su2_prog[i] = su2_prog[i]/(double)(i+1);
            err_prog[i] = error(sum_prog[i],su2_prog[i],i);
            //stampo solo media ed errore relativi all'ultimo blocco! Come ascissa uso valore centrale al bin.
            if(i == N-1){
                out << (-2.5+0.025)+(ibin*0.05) <<"   "<< sum_prog[i]/(double(M/N)*0.05) <<"   "<< err_prog[i]/(double(M/N)*0.05) << endl;
                //occhio alla normalizzazione: divido per il numero di lanci in un blocco e per la larghezza del bin
            }
        }
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

double psi(double x,double mu, double sigma){
    double result;
    result = exp(-pow(x-mu,2.)/(2*sigma*sigma)) + exp(-pow(x+mu,2.)/(2*sigma*sigma));
    return result;
}

double pot(double x){
    double result;
    result = pow(x,4.)-2.5*pow(x,2.);
    return result;
}

double psi2(double x, double mu, double sigma){
    double result;
    result = exp(-pow(x-mu,2.)/(2*sigma*sigma))*(pow((x-mu)/(sigma*sigma),2.)-1./(sigma*sigma)) + exp(-pow(x+mu,2.)/(2*sigma*sigma))*(pow((x+mu)/(sigma*sigma),2.)-1./(sigma*sigma));
    return result;
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
