/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Gabriele Bressanini
_/    _/  _/_/_/  _/_/_/_/ email: gabriele.bressanini@studenti.unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number (i will use my own library: random.h)
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

#include "random.h"

using namespace std;

double error(double,double,int);
int start_from_previous = 0;
int cont = 0;   //contatore nececessario per salvare nei vettori le stime ad ogni step di energie e temperature
int nblocks = 30;  //definisco il numero dei blocchi per il calcolo dei valori medi

//vettori che contengono le stime delle grandezze ad ogni step. RICORDO CHE MISURO OGNI 10.
static double potenziale[3000], cinetica[3000], totale[3000], temperatura[3000];


int main(){
    Input();             //Inizialization
    int nconf = 1;
    for(int istep=1; istep <= nstep; ++istep){
        Move();           //Move particles with Verlet algorithm
        if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        if(istep%10 == 0){
            Measure();     //Properties measurement
            //ConfXYZ(nconf);//Write actual configuration in XYZ format (only for visualization via Ovito)
            nconf += 1;
        }
    }
    ConfFinal();         //Write final configuration to restart
    ConfPrevious();    //Write second to last position reached by the simulation
    return 0;
}


void Input(void){ //Prepare all stuff for the simulation
    
    //Initialization of random object
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
    
    //Actual input preparation
    
    ifstream ReadInput,ReadConf;
    double ep, ek, pr, et, vir;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;
  
    ReadInput.open("input.dat"); //Read input
    ReadInput >> start_from_previous;     //if start_from_previous = 1 the programm will start the simulation from config.final
    ReadInput >> temp;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;

    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    ReadInput.close();

    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    n_props = 4; //Number of observables

    //Read second to last configuration of the previous run in order to scale velocities to match the target temperature
    
    if(start_from_previous == 1){
        cout << "Read initial configuration from file config.final " << endl << endl;
        ReadConf.open("config.final");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();
    
        cout << "Read configuration from file config.prev in order to scale velocities" << endl << endl;
        ReadConf.open("config.prev");
        for (int i=0; i<npart; ++i){
            ReadConf >> xold[i] >> yold[i] >> zold[i];
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        ReadConf.close();
        
        //find x(t+dt) with a single step of verlet algorithm and calculate velocity
     
     //NB: ho fatto diventare le xnew,ynew e znew dei VETTORI
     
        double xnew[m_part], ynew[m_part], znew[m_part], fx[m_part], fy[m_part], fz[m_part];
        
        for(int i=0; i<npart; ++i){ //Force acting on particle i
            fx[i] = Force(i,0);
            fy[i] = Force(i,1);
            fz[i] = Force(i,2);
        }
        
        double sumv[3] = {0.0, 0.0, 0.0};
        
        for(int i=0; i<npart; ++i){ //Verlet integration scheme
            
            xnew[i] = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
            ynew[i] = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
            znew[i] = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
            
            vx[i] = Pbc(xnew[i] - xold[i])/(2.0 * delta);
            vy[i] = Pbc(ynew[i] - yold[i])/(2.0 * delta);
            vz[i] = Pbc(znew[i] - zold[i])/(2.0 * delta);
            
            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        
        //scale velocities to match the target temperature
        
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
        double sumv2 = 0.0, fs;
        for (int i=0; i<npart; ++i){
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];
            
            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;
     
     //Ricalcolo delle posizioni "vecchie" (sto usando la formula al PRIMO ordine)
     
            xold[i] = x[i] - vx[i]  * delta;
            yold[i] = y[i] - vy[i]  * delta;
            zold[i] = z[i] - vz[i]  * delta;
            
        }
        return;
    }else{
        
        //Read initial configuration
        cout << "Read initial configuration from file config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();

        //Prepare initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i=0; i<npart; ++i){
            vx[i] = rnd.Rannyu() - 0.5;
            vy[i] = rnd.Rannyu() - 0.5;
            vz[i] = rnd.Rannyu() - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
        double sumv2 = 0.0, fs;
        for (int i=0; i<npart; ++i){
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];
        
            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;

        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = x[i] - vx[i] * delta;
            yold[i] = y[i] - vy[i] * delta;
            zold[i] = z[i] - vz[i] * delta;
        }
    }
        return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
    
  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
    
    //output for the mean values (these are overwritten every run)
    ofstream Epotave, Ekinave, Etotave, Tempave;
    
    Epotave.open("ave_epot.dat");
    Ekinave.open("ave_ekin.dat");
    Tempave.open("ave_temp.dat");
    Etotave.open("ave_etot.dat");

  v = 0.0; //reset observables
  t = 0.0;
    
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery
    
    potenziale[cont] = stima_pot;
    cinetica[cont] = stima_kin;
    totale[cont] = stima_etot;
    temperatura[cont] = stima_temp;

    //se è l'ultimo step allora calcolo le media a blocchi con le relative incertezze.
    //il diviso 10 arriva dal fatto che la misura viene fatta ogni 10 step
    
    double sum = 0.;
    
    if(cont == nstep/10 - 1){
        cout << "Evaluate averages and errors as functions of blocks." << endl << endl;
        
        double ave[nblocks],av2[nblocks];
        double sum_prog[nblocks],su2_prog[nblocks],err_prog[nblocks];
        
        for(int i=0; i<nblocks; i++){
            ave[i] = av2[i] = sum_prog[i] = su2_prog[i] = err_prog[i] = 0;
        }
        
        //numero di misure per ogni blocco
        int M;
        M = (nstep/10)/nblocks;
        
        //Medie ed errori a blocchi per l'energia potenziale
        for(int i=0; i<nblocks; i++){
            sum = 0;
            for(int j=0; j<M; j++){
                sum = sum + potenziale[i * M + j];
            }
            ave[i] = sum / double(M);
            av2[i] = ave[i] * ave[i];
        }
        
        for(int i=0; i<nblocks; i++){
            for(int j=0; j<i+1; j++){
                sum_prog[i] = sum_prog[i] + ave[j];
                su2_prog[i] = su2_prog[i] + av2[j];
            }
            sum_prog[i] = sum_prog[i]/(i+1);
            su2_prog[i] = su2_prog[i]/(i+1);
            err_prog[i] = error(sum_prog[i],su2_prog[i],i);
            Epotave << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
        }
        
        for(int i=0; i<nblocks; i++){
            ave[i] = av2[i] = sum_prog[i] = su2_prog[i] = err_prog[i] = 0;
        }
        
        //Medie ed errori a blocchi per l'energia cinetica
        for(int i=0; i<nblocks; i++){
            sum = 0;
            for(int j=0; j<M; j++){
                sum = sum + cinetica[i * M + j];
            }
            ave[i] = sum / double(M);
            av2[i] = ave[i] * ave[i];
        }
        
        for(int i=0; i<nblocks; i++){
            for(int j=0; j<i+1; j++){
                sum_prog[i] = sum_prog[i] + ave[j];
                su2_prog[i] = su2_prog[i] + av2[j];
            }
            sum_prog[i] = sum_prog[i]/(i+1);
            su2_prog[i] = su2_prog[i]/(i+1);
            err_prog[i] = error(sum_prog[i],su2_prog[i],i);
            Ekinave << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
        }
        
        for(int i=0; i<nblocks; i++){
            ave[i] = av2[i] = sum_prog[i] = su2_prog[i] = err_prog[i] = 0;
        }
        
        //Medie ed errori a blocchi per l'energia totale
        for(int i=0; i<nblocks; i++){
            sum = 0;
            for(int j=0; j<M; j++){
                sum = sum + totale[i * M + j];
            }
            ave[i] = sum / double(M);
            av2[i] = ave[i] * ave[i];
        }
        
        for(int i=0; i<nblocks; i++){
            for(int j=0; j<i+1; j++){
                sum_prog[i] = sum_prog[i] + ave[j];
                su2_prog[i] = su2_prog[i] + av2[j];
            }
            sum_prog[i] = sum_prog[i]/(i+1);
            su2_prog[i] = su2_prog[i]/(i+1);
            err_prog[i] = error(sum_prog[i],su2_prog[i],i);
            Etotave << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
        }
        
        for(int i=0; i<nblocks; i++){
            ave[i] = av2[i] = sum_prog[i] = su2_prog[i] = err_prog[i] = 0;
        }
        
        
        //Medie ed errori a blocchi per la temperatura
        for(int i=0; i<nblocks; i++){
            sum = 0;
            for(int j=0; j<M; j++){
                sum = sum + temperatura[i * M + j];
            }
            ave[i] = sum / double(M);
            av2[i] = ave[i] * ave[i];
        }
        
        for(int i=0; i<nblocks; i++){
            for(int j=0; j<i+1; j++){
                sum_prog[i] = sum_prog[i] + ave[j];
                su2_prog[i] = su2_prog[i] + av2[j];
            }
            sum_prog[i] = sum_prog[i]/(i+1);
            su2_prog[i] = su2_prog[i]/(i+1);
            err_prog[i] = error(sum_prog[i],su2_prog[i],i);
            Tempave << i <<"   "<< sum_prog[i] <<"   "<< err_prog[i] << endl;
        }
        
    }
    
    //aggiorno il contatore
    cont++;
    
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    
    Epotave.close();
    Ekinave.close();
    Tempave.close();
    Etotave.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}


void ConfPrevious(void){
    ofstream WriteConf;
    
    cout << "Print second to last configuration to file config.prev " << endl << endl;
    WriteConf.open("config.prev");
    
    for (int i=0; i<npart; ++i){
        WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    }
    WriteConf.close();
    return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
_/    _/  _/_/_/  _/_/_/_/ email: gabriele.bressanini@studenti.unimi.it
*****************************************************************
*****************************************************************/
