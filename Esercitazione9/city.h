#ifndef City_h
#define City_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include "random.h"

using namespace std;

int seed[4];
Random rnd;


class City{
public:
    City(){
        m_x = 0;
        m_y = 0;
        m_index = 0;
    }
    
    City(double x, double y, int index){
        m_x = x;
        m_y = y;
        m_index = index;
    }
    
    void setx(double x){
        m_x = x;
        return;
    }
    
    void sety(double y){
        m_y = y;
        return;
    }
    
    void setindex(int index){
        m_index = index;
        return;
    }
    
    double getx()const{
        return m_x;
    }
    
    double gety()const{
        return m_y;
    }
    
    double getindex()const{
        return m_index;
    }
    
private:
    double m_x,m_y;
    int m_index;
};


class Salesman{
public:
    Salesman(){
        m_cost = 0;
        
        //inizialmente genero il percorso ordinato
        for(int i=0; i<30; i++){
            m_route[i] = i;
        }
        
        //permuto il percorso: ciclo su tutte le posizioni, scegliendo il secondo indice casualmente
        for(int i=0; i<30; i++){
            int index = 0;
            index = (int)(30 * rnd.Rannyu());
            swap(m_route[i],m_route[index]);
        }
    }
    
    //setta l'iesima città del percorso
    void setroute(int citta,int i){
        m_route[i] = citta;
        return;
    }
    
    //questa funzione restituisce l'i-esima città del percorso
    int getroute(int i){
        return m_route[i];
    }
    
    //A questa funzione passo il vettore delle città. Verrà quindi calcolato il costo del percorso del Salesman
    void cost(City *c){
        double dist = 0;
        for(int i=0; i<29; i++){
            dist += sqrt( pow(c[m_route[i]].getx()-c[m_route[i+1]].getx(),2.) + pow(c[m_route[i]].gety()-c[m_route[i+1]].gety(),2.) );
        }
        //esterno al ciclo sommo anche la distanza tra il primo e l'ultimo punto
        dist += sqrt( pow(c[m_route[0]].getx()-c[m_route[29]].getx(),2.) + pow(c[m_route[0]].gety()-c[m_route[29]].gety(),2.) );
        m_cost = dist;
        return;
    }
    
    double getcost(){
        return m_cost;
    }
    
    //questa funzione controlla se il percorso del salesman rispetta i vincoli
    //richiedo che la somma degli indici delle città (da 0 a 29) faccia esattamente 435
    void check(){
        int sum=0;
        for(int i=0; i<30; i++){
            sum += m_route[i];
        }
        if(sum != 435){
            cout << "ERRORE, UN PERCORSO NON RISPETTA I VINCOLI" << endl;
        }
        return;
    }
    
private:
    double m_cost;
    int m_route[30]; //salvo il percorso del salesman in questo vettore
    
};

#endif
