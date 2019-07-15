/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Gabriele Bressanini
_/    _/  _/_/_/  _/_/_/_/ email: gabriele.bressanini@studenti.unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <fstream>
#include "city.h"

using namespace std;

double error(double,double,int);
 
int main (int argc, char *argv[]){

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
    
    ofstream best;
    best.open("bestcost.txt");  //stampo il costo migliore ad ogni generazione
    ofstream path;
    path.open("bestpath.txt");  //stampo il percorso migliore
    ofstream avg;
    avg.open("bestavg.txt");    //stampo il costo medio, ad ogni generazione, della metà migliore della popolazione di salesmen
    
    //Genero le 30 città, a cui assegno posizione (x,y) e il loro indice, che va da 0 a 29;
    //le posizioni sono punti sulla circonferenza goniometrica oppure nel quadrato di lato 1 con vertice nell'origine.
    City *cities = new City[30];
    for(int i=0; i<30; i++){
        double theta;
        theta = 2*M_PI*rnd.Rannyu();
        cities[i].setx(cos(theta));
        cities[i].sety(sin(theta));
        //cities[i].setx(rnd.Rannyu());
        //cities[i].sety(rnd.Rannyu());
        cities[i].setindex(i);
    }
    
    //Genero i 900 Salesmen
    Salesman *salesmen = new Salesman[900];
    
    //Verifico, come dovrò fare ogni volta che genero dei percorsi, che questi rispettino i vincoli
    for(int i=0; i<900; i++){
        salesmen[i].check();
    }
    
    
    for(int gen = 0; gen<200; gen++){
        //MUTAZIONI
    
        //Permutazione di una coppia
        double prob1 = 0.05; //Probabilità di effettuare una permutazione tra una coppia di città
        for(int i=0; i<900; i++){
            if(rnd.Rannyu()<prob1){
                int i1,i2;
                i1 = (int)(30 * rnd.Rannyu());  //genero gli indici delle città da scambiare
                i2 = (int)(30 * rnd.Rannyu());  //se per caso genero due indici uguali pazienza
               
                int scambio; //variabile che uso per effettuare lo scambio
                scambio = salesmen[i].getroute(i1);
                salesmen[i].setroute(salesmen[i].getroute(i2),i1);
                salesmen[i].setroute(scambio,i2);
                
                salesmen[i].check();
            }
        }
    
    
        //Shifting di TUTTE le città di un posto in avanti
        double prob2 = 0.05;
        for(int i=0; i<900; i++){
            if(rnd.Rannyu()<prob2){
                int primo;
                primo = salesmen[i].getroute(0);    //salvo il punto iniziale, che viene sovrascritto immediatamente nel ciclo
                
                for(int k=1; k<30; k++){        //comincio dal secondo elemento e vado a sostituire nell'elemento precedente
                    salesmen[i].setroute(salesmen[i].getroute(k),k-1);
                }
                
                salesmen[i].setroute(primo,29);
                
                salesmen[i].check();
            }
        }
    
    
        //Shifting delle prime 10 città di 2 posti in avanti
        double prob3 = 0.05;
        for(int i=0; i<900; i++){
            if(rnd.Rannyu()<prob3){
                int undicesima, dodicesima; //salvo la città 11 e 12
                undicesima = salesmen[i].getroute(10);
                dodicesima = salesmen[i].getroute(11);
               
                for(int k=0; k<10; k++){
                    salesmen[i].setroute(salesmen[i].getroute(9-k),9-k+2);
                }
                
                salesmen[i].setroute(undicesima,0);
                salesmen[i].setroute(dodicesima,1);
                
                salesmen[i].check();
            }
        }
    
    
        //Scambio le città dalla 10 alla 14 con quelle dalla 20 alla 24
        double prob4 = 0.05;
        for(int i=0; i<900; i++){
            if(rnd.Rannyu()<prob4){
                for(int k=0; k<5; k++){
                    int scambio;
                    scambio = salesmen[i].getroute(20+k);
                    salesmen[i].setroute(salesmen[i].getroute(10+k),20+k);
                    salesmen[i].setroute(scambio,10+k);
                }
                
                salesmen[i].check();
            }
        }
    
    
        //Inversione delle ultime 10 città
        double prob5 = 0.05;
        for(int i=0; i<900; i++){
            if(rnd.Rannyu()<prob5){
                for(int k=0; k<5; k++){
                    int scambio;
                    scambio = salesmen[i].getroute(20+k);
                    salesmen[i].setroute(salesmen[i].getroute(29-k),20+k);
                    salesmen[i].setroute(scambio,29-k);
                }
                salesmen[i].check();
            }
        }
    
    
        //Calcolo il costo per ogni salesman
        for(int i=0; i<900; i++){
            salesmen[i].cost(cities);
        }
    
        //Ordino il vettore dei salesmen in orine CRESCENTE di costo
        Salesman comodo;  //variabile di comodo che serve per gli scambi
        for(int i=0; i<899; i++){
            for(int j=i+1; j<900; j++){
                if(salesmen[i].getcost() > salesmen[j].getcost()){
                    comodo = salesmen[i];
                    salesmen[i] = salesmen[j];
                    salesmen[j] = comodo;
                }
            }
        }
    
    
        //Selezione (voglio che vengano scelti con più frequenza quelli con costi più piccoli)
        //Dopo averli selezionati essi fanno crossover, con probabilità p6, per generare due figli.
        //Se non effettuano crossover mando avanti i due genitori alla generazione successiva.
    
        double prob6 = 0.80;
        double P = 5; //potenza del numero casuale tra 0 ed 1. Vedere JN del professore per spiegazione
    
        Salesman *new_salesmen = new Salesman[900];
    
        //faccio 450 crossovers per creare 900 figli che faranno parte della nuova generazione
        for(int i=0; i<450; i++){
            int i1, i2; //indici dei canditati genitori
            i1 = (int)(900 * pow(rnd.Rannyu(),P));
            i2 = (int)(900 * pow(rnd.Rannyu(),P));
            if(rnd.Rannyu() < prob6){
                int cut;
                cut = (int)(rnd.Rannyu(1.,29.9));     //scelgo casualmente il sito dove avviene il taglio
                
                //copio nei due figli il percorso di mamma(papà) fino al taglio. NB: se taglio è 0 allora solo la prima città viene copiata etc.
                for(int k=0; k<cut; k++){
                    new_salesmen[i].setroute(salesmen[i1].getroute(k),k);
                    new_salesmen[899-i].setroute(salesmen[i2].getroute(k),k);
                }
                
                //ora devo aggiungere le città mancanti, prendendole IN ORDINE dal consorte
                //PRIMO FIGLIO
                for(int k=cut; k<30; k++){      //ciclo sugli elementi mancanti del figlio
                    for(int j=0; j<30; j++){    //ciclo sugli elementi del genitore "COMPLEMENTARE"
                        int trovato = 0;
                        int m=0;
                        do{
                            if(salesmen[i2].getroute(j) == new_salesmen[i].getroute(m)){
                                trovato = 1;
                            }
                            m++;
                        }while(trovato == 0 && m<k);
                        if(trovato == 0){
                            new_salesmen[i].setroute(salesmen[i2].getroute(j),k);
                            break;
                        }
                    }
                }
                
                //SECONDO FIGLIO
                for(int k=cut; k<30; k++){      //ciclo sugli elementi mancanti del figlio
                    for(int j=0; j<30; j++){    //ciclo sugli elementi del genitore "COMPLEMENTARE"
                        int trovato = 0;
                        int m=0;
                        do{
                            if(salesmen[i1].getroute(j) == new_salesmen[899-i].getroute(m)){
                                trovato = 1;
                            }
                            m++;
                        }while(trovato == 0 && m<k);
                        if(trovato == 0){
                            new_salesmen[899-i].setroute(salesmen[i1].getroute(j),k);
                            break;
                        }
                    }
                }
                new_salesmen[i].check();
                new_salesmen[899-i].check();
                
                
            }else{
                //faccio le copie dei genitori nei figli
                new_salesmen[i] = salesmen[i1];
                new_salesmen[i].check();
                new_salesmen[899-i] = salesmen[i2];
                new_salesmen[899-i].check();
                
            }
        }
    
        //ora copio new salesmen in salesmen, e il ciclo riparte...
        for(int k=0; k<900; k++){
            salesmen[k] = new_salesmen[k];
            salesmen[k].check();
        }
        delete[]new_salesmen;
        
        //Calcolo il costo per ogni salesman
        for(int i=0; i<900; i++){
            salesmen[i].cost(cities);
        }
        
        //Ordino il vettore dei salesmen in orine CRESCENTE di costo
        for(int i=0; i<899; i++){
            for(int j=i+1; j<900; j++){
                if(salesmen[i].getcost() > salesmen[j].getcost()){
                    comodo = salesmen[i];
                    salesmen[i] = salesmen[j];
                    salesmen[j] = comodo;
                }
            }
        }
         
        
        
        double sum = 0;     //media del costo della metà migliore della popolazione
        for(int l=0; l<450; l++){
            sum += salesmen[l].getcost();
        }
        sum /= 450.;
        
        avg << gen << "   " << sum << endl;
        
        //Stampo in un file il costo minore ad ogni generazione.
        best << gen <<"   "<< salesmen[0].getcost() << endl;
         
    }
    
    for(int k=0; k<30; k++){
        double x,y;
        x = cities[salesmen[0].getroute(k)].getx();
        y = cities[salesmen[0].getroute(k)].gety();
        path << x << "   " << y << endl;
    }
    //ora ristampo nuovamente le coordinate della prima città, in modo che nel grafico sia collegata all'ultima
    path << cities[salesmen[0].getroute(0)].getx()<<"   "<< cities[salesmen[0].getroute(0)].gety() << endl;
    
    best.close();
    best.clear();
    path.close();
    path.clear();
    avg.close();
    avg.clear();
    
    delete[]cities;
    delete[]salesmen;
   
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
 _/    _/       _/ _/       Gabriele Bressanini
_/    _/  _/_/_/  _/_/_/_/ email: gabriele.bressanini@gmail.com
*****************************************************************
*****************************************************************/
