#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include "../../common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <numeric>

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{
private:
    // Existing private members of Engine...
  const double pi=3.1415926535897932384626433832795028841971e0;

  // definition des variables

  double g, m, L, Omega, r, kappa;         // accélération gravitationnelle, masse, longueur, fréquence angulaire, rayon, coefficient de frottement


  double theta;
  double thetadot;
  double F_theta;
  

  double t;  // Temps courant pas de temps
  double tf;          // Temps final
  double dt;      // Intervalle de temps
  int N_excit;  // Nombre de périodes d'excitation
  int nsteps_per; // Nombre de pas de temps par période d'excitation

  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  void printOut(bool write)
  {

    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      double emec = Emec(theta, thetadot, t); // TODO: Evaluer l'energie mecanique
      double pnc = Pnonc(theta, thetadot, t); // TODO: Evaluer la puissance des forces non conservatives
      *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // TODO definir l'énergie mechanique
  double Emec(double theta, double thetadot, double t_)
  {

      return 0.5*m*pow(L*thetadot,2) + m*g*L*(1-cos(theta));
  }

  // TODO definir la puissance des forces non conservatives
  double Pnonc(double theta, double thetadot, double t_)
  {

      return - kappa*(pow(L*thetadot,2) + L*r*thetadot*Omega*cos(Omega*t_ - theta)) + m*L*r*thetadot*pow(Omega,2)*sin(Omega*t_-theta);
  }

  // TODO écrire la fonction pour l'acceleration (theta_doubledot)
  double compute_acc(double theta, double thetadot, double t_)
  {
      double acc = -g/L*sin(theta) - kappa/m*(thetadot + r*Omega*cos(Omega*t_ - theta)/L) + r*Omega*Omega*sin(Omega*t_ - theta)/L;


      return acc;
  }
  // TODO implementer le schéma Velocity Verlet pour une accélération dependante du theta, thetadot et t.

  void step()
  {



    double F_thetadot = -kappa/m*thetadot;

    double thetadot_half = thetadot + 0.5*(F_theta + F_thetadot)*dt;

    double F_thetadot_half = -kappa/m*thetadot_half;
    double theta_1 = theta + thetadot_half*dt;

    double F_theta_1 = compute_acc(theta_1, thetadot, t+dt) + kappa/m*thetadot;

    thetadot = thetadot + dt*F_theta/2 + dt*F_theta_1/2 + dt*F_thetadot_half/2; 
    F_theta = F_theta_1;
    theta = theta_1;

    t += dt;
  }

public:
    // Modified constructor
    Engine(ConfigFile configFile)
    {
      // Stockage des parametres de simulation dans les attributs de la classe
      tf     = configFile.get<double>("tf",tf);	        // t final (overwritten if N_excit >0)
      g     = configFile.get<double>("g", g);         // lire l'acceleration de gravite
      m     = configFile.get<double>("m", m);         // lire la masse
      L     = configFile.get<double>("L", L);         // lire la longueur
      Omega = configFile.get<double>("Omega", Omega); // lire la frequence angulaire
      r     = configFile.get<double>("r", r);         // lire le rayon
      kappa = configFile.get<double>("kappa", kappa); // lire le coefficient de frottement
      theta    = configFile.get<double>("theta0", theta);    // lire la condition initiale en theta
      thetadot = configFile.get<double>("thetadot0", thetadot); // lire la condition initiale en thetadot

      N_excit  = configFile.get<int>("N");            // number of periods of excitation
      nsteps_per= configFile.get<int>("nsteps");        // number of time step per period
      sampling = configFile.get<unsigned int>("sampling",sampling); // lire le nombre de pas de temps entre chaque ecriture des diagnostics

      // Ouverture du fichier de sortie
      outputFile = new ofstream(configFile.get<string>("output").c_str());
      outputFile->precision(15);
      if(N_excit>0){
        tf = N_excit*(2*pi/Omega);
        dt   = (2*pi/Omega)/nsteps_per;
      }
      else{
        dt = tf/nsteps_per;
      }
      
    };


    // Destructeur virtuel
    virtual ~Engine()
    {
      outputFile->close();
      delete outputFile;
    };
      // Simulation complete
    void run()
    {
      t = 0.;
      last = 0;

      F_theta = compute_acc(theta, thetadot, t) + kappa/m*thetadot;

      printOut(true);

      while( t < tf-0.5*dt )
      {
        step();
        printOut(false);
      }
      printOut(true);
    };
};

// programme
int main(int argc, char* argv[])
{
  // Existing main function implementation
  // ...
  string inputPath("configuration.in.example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

  Engine* engine;

  // Create an instance of Engine instead of EngineEuler
  engine = new Engine(configFile);

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}


