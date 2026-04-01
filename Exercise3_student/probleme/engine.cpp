#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include "../../common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <numeric>
#include <valarray>

using namespace std;

class Engine {
    private:

    static constexpr double RT = 6378.1*1e3; // rayon de la Terre
    static constexpr double G = 6.674*1e-11; // constante gravitationnele
    static constexpr double lambda = 7232.2; //épaisseur caractéristique
    static constexpr double Cx = 0.3; // coefficient de traînée
    static constexpr double mT =5,972*1e24; //masse de la Terre
    static constexpr double pi 3.1415926535897932384626433832795028841971e0; 


    //parametres de simulation
    double rho0; //densité de l'air au niveau de la mer
    double dA; //diamètre de la sonde
    double mA; //masse de la sonde
    double mL; //masse de la Lune

    double epsilon; //précision
    double t; //temps courant
    double dt; //pas de temps
    double tf; //temps de la simulation

    //variable d'itération contenant les positions et vitesses des 3 corps
    valarray<double> y = std::valarray<double>(0.e0, 12);

    unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
    unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
    ofstream *outputFile;    // Pointeur vers le fichier de sortie


    void printOut(bool write){}

    double Emec(){} //Energie mécanique

    double p() {} //quantité de mouvement

    void step(){}

    public:

    Engine(ConfigFile configFile) {};

    virtual ~Engine(){
        outputFile->close();
        delete outputFile;
    };

    void run() {};
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






