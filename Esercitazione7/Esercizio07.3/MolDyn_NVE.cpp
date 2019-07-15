/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Gabriele Bressanini
_/    _/  _/_/_/  _/_/_/_/ email: gabriele.bressanini@studenti.unimi.it
****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number (i will use my own library: random.h)
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

double Error(double,double,int);
int start_from_previous = 0;
int nblk = 30;  //definisco il numero dei blocchi per il calcolo dei valori medi

int main()
{
    Input(); //Inizialization
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
        cout << "Blocco:" << iblk << endl << endl;
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
            Move();
            Measure();
            Accumulate(); //Update block averages
        }
        Averages(iblk);   //Print results for current block
    }
    GetGave();
    ConfFinal(); //Write final configuration
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

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;
  
    ReadInput.open("input.gas"); //Read input
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
    
    //measurement of g(r)
    igofr = 4;
    nbins = 100;
    n_props = n_props + nbins;
    bin_size = (box/2.0)/((double)(nbins));
    cout << "Bin Size = " << bin_size << endl;

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

  double t;
  double dx, dy, dz, dr;
  ofstream Temp;
    
    //reset the hystogram of g(r)
    for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  Temp.open("output_temp.dat",ios::app);

  t = 0.0; //reset observables
    
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

        dx = Pbc( x[i] - x[j] );
        dy = Pbc( y[i] - y[j] );
        dz = Pbc( z[i] - z[j] );

        dr = dx*dx + dy*dy + dz*dz;
        dr = sqrt(dr);
        
        //update of the histogram of g(r)
        for (int i_bin = 0; i_bin < nbins; i_bin++){
            if(dr>(bin_size*i_bin) && dr<=(bin_size*(i_bin+1))) walker[igofr+i_bin]+=2;
        }
    }
  }
  

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature

  
    //nei file di output salvo le variabili in unità LJ. (utile per la prima fase si equilibrazione)
    Temp << stima_temp << endl;

    Temp.close();
    return;
}

void Reset(int iblk) //Reset block averages
{
    
    if(iblk == 1)
    {
        for(int i=0; i<n_props; ++i)
        {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }
    
    for(int i=0; i<n_props; ++i)
    {
        blk_av[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}

void Accumulate(void) //Update block averages
{
    
    for(int i=0; i<n_props; ++i)
    {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
    
    ofstream Gofr, Gave;
    
    //g(r)
    
    for (int i_bin = 0; i_bin < nbins; i_bin++){
        double delta_V = (4.*M_PI/3.)*(pow((bin_size*(i_bin+1)),3)-pow((bin_size*i_bin),3));
        stima_gofr = blk_av[igofr+i_bin]/blk_norm/(npart*delta_V*rho);
        glob_av[igofr+i_bin] += stima_gofr;
        glob_av2[igofr+i_bin] += stima_gofr*stima_gofr;
    }
    
    cout << "----------------------------" << endl << endl;
}

void GetGave(void){
    const int wd=12;
    ofstream Gave;
    Gave.open("output.gave.0");
    for (int i_bin = 0; i_bin < nbins; i_bin++){
        err_gdir=Error(glob_av[igofr+i_bin],glob_av2[igofr+i_bin],nblk);
        Gave << setw(wd)<< i_bin << setw(wd)  << glob_av[igofr+i_bin]/(double)nblk <<  setw(wd) << err_gdir <<  endl;
    }
    Gave.close();
    Gave.clear();
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

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
