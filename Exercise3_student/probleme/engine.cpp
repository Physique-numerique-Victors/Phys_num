#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include "../common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <array>

using namespace std;


class Engine {
    private:

    static const array<double, 3> R; // rayon de la Terre, de la Lune et de la sonde
    static constexpr double G = 6.674*1e-11; // constante gravitationnele
    static constexpr double lambda = 7238.2; //épaisseur caractéristique
    static constexpr double Cx = 0.3; // coefficient de traînée
    static constexpr double pi = 3.1415926535897932384626433832795028841971e0; 
    static constexpr size_t Terre = 0; //indice de la Terre dans y
    static constexpr size_t Lune = 1; //indice de la Lune dans y
    static constexpr size_t Art = 2; //indice de Artemis dans y
    static constexpr unsigned int numbodies = 3; //nombre de corps
    static constexpr unsigned int dimension = 2; //dimension du problème

    //parametres de simulation
    double rho0; //densité de l'air au niveau de la mer
    array<double, 3> m; //masse de la Terre, de la Lune et de la sonde

    double epsilon; //précision du pas de temps adaptatif
    double t; //temps courant
    double dt; //pas de temps
    double tf; //temps de la simulation
    bool adaptatif;

    //variable d'itération contenant les positions et vitesses des 3 corps
    valarray<double> y = std::valarray<double>(0.e0, numbodies * dimension * 2);

    unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
    unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
    ofstream *outputFile;    // Pointeur vers le fichier de sortie

    size_t ix(size_t i) const { return 2 * i; }
    size_t iy(size_t i) const { return 2 * i + 1; }
    size_t ivx(size_t i) const { return 2 * numbodies + 2 * i; }
    size_t ivy(size_t i) const { return 2 * numbodies + 2 * i + 1; }

    valarray<double> r(size_t i, const valarray<double>& state) const {return state[slice(ix(i), dimension, 1)];}

    valarray<double> v(size_t i, const valarray<double>& state) const { return state[slice(ivx(i), dimension, 1)];}   

    void printOut(bool write){
        if((!write && last>=sampling) || (write && last!=1))
        {
            *outputFile << t << " " <<  y[ix(Art)] << " " <<  y[iy(Art)] << " " << y[ivx(Art)] << " " << y[ivy(Art)] << " " << dt << " " << norm(acc_Artemis()) << " " << Pt_Art() << endl;
            last = 1;
        }
        else
        {
            last++;
        }
    }
    double norm(const valarray<double>& V) const {return sqrt(pow(V[0] , 2) + pow(V[1] , 2));} //norme de vecteurs de deux dimensions

    double EmecArt() const { //Energie mécanique d'Artemis
        return 0.5 * m[Art] * pow(norm(v(Art, y) - v(Terre, y)), 2) - G * m[Terre] * m[Art] / norm(r(Art, y) - r(Terre, y)) - G * m[Lune] * m[Art] / norm(r(Art, y) - r(Lune, y));
    } 

    double p(const valarray<double>& v, double m) const { //quantité de mouvement normée
        return m*norm(v);
    } 

    valarray<double> acc_Artemis() const {
        valarray<double> rT = r(Terre, y);
        valarray<double> vT = v(Terre, y);
        valarray<double> rL = r(Lune, y);
        valarray<double> rA = r(Art, y);
        valarray<double> vA = v(Art, y);

        return acc_grav(rA, rT, m[Terre]) + Ft(rA-rT, vA-vT)/m[Art] + acc_grav(rA, rL, m[Lune]);
    }

    double Pt_Art() const {
        valarray<double> rT = r(Terre, y);
        valarray<double> rA = r(Art, y);
        valarray<double> vrel = v(Art, y) - v(Terre, y);
        valarray<double> Fdrag = Ft(rA - rT, vrel);

        return (Fdrag * vrel).sum();
    }

    //méthode pour initialiser la position et la vitesse des corps dans y
    void set_body(const ConfigFile& configFile, size_t body, const string& x, const string& y_, const string& vx, const string& vy){
    y[slice(ix(body), 2, 1)] = valarray<double>{configFile.get<double>(x, 0.0), configFile.get<double>(y_, 0.0)};

    y[slice(ivx(body), 2, 1)] = valarray<double>{configFile.get<double>(vx, 0.0), configFile.get<double>(vy, 0.0)};
    }


    valarray<double> acc_grav(const valarray<double>& r, const valarray<double>& r_a, double m) const { return -G * m *(r-r_a)/(pow(norm(r-r_a) , 3)); } //Force gravitationnelle/m

    double rho(const valarray<double>& r) const {return rho0 * exp(-(norm(r) - R[Terre])/lambda);}

    valarray<double> Ft(const valarray<double>& r, const valarray<double>& v) const { return -0.5 * rho(r) * pi * pow(R[Art], 2) * Cx * norm(v) * v;} //Force de trainée aérodynamique

    void compute_f(const valarray<double>& y_, valarray<double>& df) const {
        valarray<double> rT = r(Terre, y_);
        valarray<double> vT = v(Terre, y_);
        valarray<double> rL = r(Lune, y_);
        valarray<double> vL = v(Lune, y_);
        valarray<double> rA = r(Art, y_);
        valarray<double> vA = v(Art, y_);

        //evolution de la Terre

        valarray<double> vT_after = acc_grav(rT, rL, m[Lune]) + acc_grav(rT, rA, m[Art]);

        df[ix(Terre)] = vT[0];
        df[iy(Terre)] = vT[1];
        df[ivx(Terre)] = vT_after[0];
        df[ivy(Terre)] = vT_after[1];


        //evolution de la Lune

        valarray<double> vL_after = acc_grav(rL, rT, m[Terre]) + acc_grav(rL, rA, m[Art]);

        df[ix(Lune)] = vL[0];
        df[iy(Lune)] = vL[1];
        df[ivx(Lune)] = vL_after[0];
        df[ivy(Lune)] = vL_after[1];

        //evolution d'Artemis
        
        valarray<double> vA_after= acc_grav(rA, rT, m[Terre]) + Ft(rA-rT, vA-vT)/m[Art] + acc_grav(rA, rL, m[Lune]);

        df[ix(Art)] = vA[0];
        df[iy(Art)] = vA[1];
        df[ivx(Art)] = vA_after[0];
        df[ivy(Art)] = vA_after[1];

    }

valarray<double> rk4Step(double step, const valarray<double>& y0){
        valarray<double> y_ = y0;
        valarray<double> f = valarray<double>(0.e0, numbodies * dimension * 2); 

        compute_f(y_, f);
        valarray<double> k1 = f;

        y_ = y0 + 0.5 * k1 * step;
        compute_f(y_, f);
        valarray<double> k2 = f;

        y_ = y0 + 0.5 * k2 * step;
        compute_f(y_, f);
        valarray<double> k3 = f;

        y_ = y0 + k3 * step;
        compute_f(y_, f);
        valarray<double> k4 = f;

        y_ = y0 + step * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
        
        return y_;
    }

    bool checkCollisions() const { //retourne true si il y a une collision entre les corps et false sinon

        if ((norm(r(Terre,y) - r(Art,y)) >= R[Terre] + R[Art]) and (norm(r(Terre,y) - r(Lune,y)) >= R[Terre] + R[Lune]) and (norm(r(Lune,y) - r(Art,y)) >= R[Lune] + R[Art])) { return false;}
        else {return true;}
    }

    public:

    Engine(const ConfigFile& configFile) {
        rho0 = configFile.get<double>("rho0", rho0);
        m[Art] = configFile.get<double>("mA", 0);
        m[Lune] = configFile.get<double>("mL", 0);
        m[Terre] = configFile.get<double>("mT", 0);
        epsilon = configFile.get<double>("epsilon", epsilon);
        dt = configFile.get<double>("dt",dt);
        tf = configFile.get<double>("tf",tf);
        adaptatif = configFile.get<bool>("adaptatif",adaptatif);

        set_body(configFile, Terre, "xT", "yT", "vxT", "vyT");
        set_body(configFile, Lune,  "xL", "yL", "vxL", "vyL");
        set_body(configFile, Art,   "xA", "yA", "vxA", "vyA");

        sampling = configFile.get<unsigned int>("sampling", sampling);

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

        while((t < tf-0.5*dt) and !(checkCollisions())){

            if(adaptatif){

                bool step_accepted = false;

                while (!step_accepted) {
                     double error = 0.0;

                    valarray<double> y1 = rk4Step(dt, y);
                    valarray<double> y2 = rk4Step(dt*0.5, rk4Step(dt*0.5, y));

                    for(size_t i = 0; i < y.size(); ++i){
                        error += pow(y1[i] - y2[i], 2);
                    }

                    error =sqrt(error);
                    if (error <= epsilon) {
                        y = y2;
                        step_accepted = true;
                        t += dt;

                        // Si l'erreur est très petite, on peut augmenter le pas de temps pour accélérer la simulation
                        if (error < epsilon * 0.1) {
                            dt *= 2.0;
                        }

                    } else {
                        dt *= 0.5;
                    }
                }
                
            } else {
                y = rk4Step(dt, y);
                t += dt;
            }

            printOut(false);       
            
        }
        
        printOut(true);
    };
};

const std::array<double, 3> Engine::R = {6378.1e3, 1737e3, 5.02 * 0.5};

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






