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
    static constexpr double dA = 5.02; //diamètre de la sonde
    static constexpr double pi 3.1415926535897932384626433832795028841971e0; 
    static constexpr size_t Terre = 0; //indice de la Terre dans y
    static constexpr size_t Lune = 1; //indice de la Lune dans y
    static constexpr size_t Art = 2; //indice de Artemis dans y
    static constexpr numbodies = 3; //nombre de corps
    static constexpr dimension = 2; //dimension du problème



    //parametres de simulation
    double rho0; //densité de l'air au niveau de la mer
    double mA; //masse de Artemis
    double mL; //masse de la Lune

    double epsilon; //précision
    double t; //temps courant
    double dt; //pas de temps
    double tf; //temps de la simulation
    bool adaptatif;

    //variable d'itération contenant les positions et vitesses des 3 corps
    valarray<double> y = std::valarray<double>(0.e0, numbodies*dimension*2);

    unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
    unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
    ofstream *outputFile;    // Pointeur vers le fichier de sortie


    void printOut(bool write){
        if((!write && last>=sampling) || (write && last!=1))
        {
            double emec = Emec(); 
            double p = p(); 
            *outputFile << t << " " <<  y[ix(Art)] << " " <<  y[iy(Art)] << " " << emec << " " << p << endl;
            last = 1;
        }
        else
        {
            last++;
        }
    }

    double Emec() const {} //Energie mécanique

    double p() const {return mA*sqrt(pow(y[ivx(Art)],2) + pow(y[ivy(Art)]))} //quantité de mouvement

    size_t ix(size_t i) const { return 2 * i; }
    size_t iy(size_t i) const { return 2 * i + 1; }
    size_t ivx(size_t i) const { return 2 * numBodies + 2 * i; }
    size_t ivy(size_t i) const { return 2 * numBodies + 2 * i + 1; }

    void compute_f(valarray<double>& f)const {
        //evolution de la lune
        valarray<double> r L= valarray<double> (y[ix(1)], y[iy(1)]);
        valarray<double> vL = valarray<double> (y[ivx(1)], y[ivy(1)]);

        vL_after = -G*mT*rL/(np.linalg.norm(rL)**3);

        f[[ix(1)]] = vL[0];
        f[[iy(1)]] = vL[1];
        f[ivx(1)] = vl_after[0];
        f[ivx(1)] = vl_after[1];

        //evolution d'Artemis
        valarray<double> rA = valarray<double> (y[ix(1)], y[iy(1)]);
        valarray<double> vA = valarray<double> (y[ivx(1)], y[ivy(1)]);

        vA_after = -G*mT*rA/(np.linalg.norm(rA)**3);

        f[[ix(1)]] = vA[0];
        f[[iy(1)]] = vA[1];
        f[ivx(1)] = vA_after[0];
        f[ivx(1)] = vA_after[1];

    }

    void step(){
        valarray<double> f =valarray<double>(0.e0,12); 

        compute_f(f);
        valarray<double> k1 = f;

        y += 0.5*k1;
        compute_f(f);
        valarray<double> k2 = f;

        y += -0.5*k1 + 0.5*k2;
        compute_f(f);
        valarray<double> k3 = f;

        y += -0.5*k2 + k3
        compute_f(f);
        valarray<double> k4 = f;

        y += -k3 + dt*(1.0/6)*(k1+2*k2+2*k3+k4);
    }

    public:

    Engine(ConfigFile configFile) {
        rho0 = configFile.get<double>("rho0", rho0);
        mA = configFile.get<double>("mA", mA);
        mL = configFile.get<double>("mL", mL);
        dt = configFile.get<double>("dt",dt);
        tf = configFile.get<double>("tf",tf);
        y[ix(Terre)] = configFile.get<double>("xT",xT);
        y[iy(Terre)] = configFile.get<double>("yT",yT);
        y[ivx(Terre)] = configFile.get<double>("vxT",vxT);
        y[ivy(Terre)] = configFile.get<double>("vyT",vyT);
        y[ix(Lune)] = configFile.get<double>("xL",xL);
        y[iy(Lune)] = configFile.get<double>("yL",yL);
        y[ivx(Lune)] = configFile.get<double>("vxL",vxL);
        y[ivy(Lune)] = configFile.get<double>("vyL",vyL);
        y[ix(Art)] = configFile.get<double>("xA",xA);
        y[iy(Art)] = configFile.get<double>("yA",yA);
        y[ivx(Art)] = configFile.get<double>("vxA",vxA);
        y[ivy(Art)] = configFile.get<double>("vyA",vyA);
        sampling = configFile.get<double>("sampling", sampling)

        //Ouverture du fichier de sortie
        outputFile = new ofstream(configFile.get<string>("output").c_str());
        outputFile->precision(15);


    };

    virtual ~Engine(){
        outputFile->close();
        delete outputFile;
    };

    void run() {
        t = 0;
        last = 0;

        printOut(true);

        while(t < tf-0.5*dt){
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






