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
#include <string>
#include "city.h"
#include "mpi.h"

//DISCLAIMER : Se si vuole stampare l'accettazione delle 3 mosse (mutazioni) bisogna scommentare una parte di codice in fondo.

using namespace std;
 
int main (int argc, char *argv[]){
    
    MPI::Init(argc,argv);

    int size = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();
    string r = to_string(rank);
    
    ofstream best;
    best.open(("bestcost"+r+".txt").c_str());  //stampo il costo migliore ad ogni temperatura
    ofstream path;
    path.open(("bestpath"+r+".txt").c_str());  //stampo il percorso migliore
    
    
    //Inizializzazione generatore di numeri casuali. Uso lo stesso per tutti i nodi, in modo tale da avere sempre
    //la stessa disposizione delle città. Successivamente reinizializzo in modo da dare ad ognuno due Primes diversi.
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
    
    //Genero le 30 città, a cui assegno posizione (x,y) e il loro indice, che va da 0 a 29;
    //le posizioni sono punti sulla circonferenza goniometrica oppure nel quadrato di lato 1 con vertice nell'origine.
    City *cities = new City[30];
    for(int i=0; i<30; i++){
        //double theta;
        //theta = 2*M_PI*rnd.Rannyu();
        //cities[i].setx(cos(theta));
        //cities[i].sety(sin(theta));
        cities[i].setx(rnd.Rannyu());
        cities[i].sety(rnd.Rannyu());
        cities[i].setindex(i);
    }
    
    //Reinizializzo il generatore di numeri casuali, in modo da dare ad ogni nodo una coppia di primes diversi.
    Primes.open("Primes");
    if (Primes.is_open()){
        for(int i=0; i<rank+1; i++){    //leggo primes differenti per ogni nodo
            Primes >> p1 >> p2 ;
        }
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    input.open("seed.in");
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
    
    //Genero un solo individuo
    Salesman salesmen;
    
    //Nel best_salesmen salvo progressivamente l'individuo avente il percorso ottimale fino a quel momento
    Salesman best_salesmen;
    
    //Inizialmente il migliore è il primo percorso che genero
    best_salesmen = salesmen;
    
    //Verifico, come dovrò fare ogni volta che genero dei percorsi (o "salesmen"), che questi rispettino i vincoli
    salesmen.check();
    
    double accepted1, accepted2, accepted3; //accettazione delle tre mutazioni proposte
    double beta; //nelle unità utilizzate è l'inverso della temperatura
    double nstep = 1000; //numero di step che compio ad ogni temperatura. Ad ogni step provo ad effettuare tutte e tre le mutazioni
    Salesman new_salesmen; //new_salesmen viene usato per "testare" la mutazione e calcolare la A del Metropolis
    double prob1, prob2, prob3; //probabilità di chiamata della i-esima mutazione
    prob1 = prob2 = prob3 = 0.5;
    double A = 0;   //accettazione dei Metropolis
    
    //Ciclo sulle temperature, che vengono portate gradualmente a zero
    for(double T = 1; T>0; T=T-0.001){
        
        accepted1 = accepted2 = accepted3 = 0;
        beta = 1./T;
        
        for(int step = 0; step<nstep; step++){
            
            //MUTAZIONE1 : Permutazione di una coppia
            
            new_salesmen = salesmen;
            if(rnd.Rannyu()<prob1){
                int i1,i2;
                i1 = (int)(30 * rnd.Rannyu());  //genero gli indici delle città da scambiare
                i2 = (int)(30 * rnd.Rannyu());  //se per caso genero due indici uguali pazienza
                
                int scambio; //variabile che uso per effettuare lo scambio
                scambio = new_salesmen.getroute(i1);
                new_salesmen.setroute(new_salesmen.getroute(i2),i1);
                new_salesmen.setroute(scambio,i2);
                
                new_salesmen.check();
            }
       
            new_salesmen.cost(cities);
            salesmen.cost(cities);
            
            A = min(1.,exp(-beta*(new_salesmen.getcost() - salesmen.getcost())));
            
            if(rnd.Rannyu() < A){
                salesmen = new_salesmen;
                accepted1++;
            }
            
            salesmen.cost(cities);
            best_salesmen.cost(cities);
            
            if(salesmen.getcost() < best_salesmen.getcost()){
                best_salesmen = salesmen;
            }
            
            
            //MUTAZIONE2 : Scambio le 5 città dalla 10 alla 15 e le 5 città dalla 20 alla 25;
            
            new_salesmen = salesmen;
            if(rnd.Rannyu()<prob2){
                for(int k=0; k<5; k++){
                    int scambio;
                    scambio = new_salesmen.getroute(20+k);
                    new_salesmen.setroute(new_salesmen.getroute(10+k),20+k);
                    new_salesmen.setroute(scambio,10+k);
                }
                
                new_salesmen.check();
            }
            
            new_salesmen.cost(cities);
            salesmen.cost(cities);
            
            A = min(1.,exp(-beta*(new_salesmen.getcost() - salesmen.getcost())));
            
            if(rnd.Rannyu() < A){
                salesmen = new_salesmen;
                accepted2++;
            }
            
            salesmen.cost(cities);
            best_salesmen.cost(cities);
            
            if(salesmen.getcost() < best_salesmen.getcost()){
                best_salesmen = salesmen;
            }
            
            
            //MUTAZIONE3 : Inversione delle ultime 10 città
            
            new_salesmen = salesmen;
            if(rnd.Rannyu()<prob3){
                for(int k=0; k<5; k++){
                    int scambio;
                    scambio = new_salesmen.getroute(20+k);
                    new_salesmen.setroute(new_salesmen.getroute(29-k),20+k);
                    new_salesmen.setroute(scambio,29-k);
                }
                new_salesmen.check();
            }
           
            new_salesmen.cost(cities);
            salesmen.cost(cities);
            
            A = min(1.,exp(-beta*(new_salesmen.getcost() - salesmen.getcost())));
            
            if(rnd.Rannyu() < A){
                salesmen = new_salesmen;
                accepted3++;
            }
            
            salesmen.cost(cities);
            best_salesmen.cost(cities);
            
            if(salesmen.getcost() < best_salesmen.getcost()){
                best_salesmen = salesmen;
            }
            
            best_salesmen.cost(cities);
            
        }
        /*
        cout << "T = " << T << endl;
        cout << "Accettazione Mutazione 1 : " << (double)accepted1/nstep << endl;
        cout << "Accettazione Mutazione 2 : " << (double)accepted2/nstep << endl;
        cout << "Accettazione Mutazione 3 : " << (double)accepted3/nstep << endl << endl;
        */
        best << T << "   " << best_salesmen.getcost() << endl;
        
    }
    
    for(int k=0; k<30; k++){
        double x,y;
        x = cities[best_salesmen.getroute(k)].getx();
        y = cities[best_salesmen.getroute(k)].gety();
        path << x << "   " << y << endl;
    }

    
    //ora ristampo nuovamente le coordinate della prima città, in modo che nel grafico sia collegata all'ultima
    path << cities[best_salesmen.getroute(0)].getx()<<"   "<< cities[best_salesmen.getroute(0)].gety() << endl;
    

    double migliori[size];  //in questo vettore salvo la lunghezza del percorso migliore di ogni nodo
    for(int i=0; i<size; i++){
            migliori[i] = 0;
    }
    
    //Tutti i nodi comunicano al nodo zero il loro costo migliore
    double costo;
    costo = best_salesmen.getcost();    //costo migliore di ogni nodo
    MPI_Gather(&costo,1,MPI_DOUBLE_PRECISION,migliori,1,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);
    
    //al nodo zero stampo i valori dei percorsi migliori e trovo il migliore in assoluto
    if(rank == 0){
        double min = 99999.;
        int rank_minimo = 0;
        cout << "Lunghezza dei migliori percorsi di ogni nodo:"<< endl;
        for(int i=0; i<size; i++){
            cout << "Nodo " << i <<", lunghezza miglior percorso: " << migliori[i] << endl;
            if(migliori[i] < min){
                min = migliori[i];
                rank_minimo = i;
            }
        }
        cout << "Il nodo che ha trovato il percorso migliore è il numero " << rank_minimo << endl;
        cout << "Il percorso migliore in assoluto trovato ha costo pari a " << min << endl;
    }
    
    best.close();
    best.clear();
    path.close();
    path.clear();
    
    MPI::Finalize();
   
    return 0;
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
