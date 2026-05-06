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
            valarray<double> y_inG = in_G();
            *outputFile << t << " ";

            for (int i=0; i <3; i++) {

                  *outputFile <<  y_inG[ix(i)] << " " <<  y_inG[iy(i)] << " " << y_inG[ivx(i)] << " " << y_inG[ivy(i)] << " ";

            }
            *outputFile << dt << " " << norm(acc_Artemis_inG()) << " " << Pt_Art_inG()  << endl;
           
            last = 1;
        }
        else
        {
            last++;
        }
    }
    double norm(const valarray<double>& V) const {return sqrt(pow(V[0] , 2) + pow(V[1] , 2));} //norme de vecteurs de deux dimensions
    
    valarray<double> in_G() const{ //convertit le vecteur y dans le référentiel du centre de masse G M
        valarray<double> rG = (m[Terre]*r(Terre, y)+m[Lune]*r(Lune, y))/(m[Terre]+m[Lune]);
        valarray<double> vG = (m[Terre]*v(Terre, y)+m[Lune]*v(Lune, y))/(m[Terre]+m[Lune]);
        valarray<double> rT_inG = r(Terre, y) - rG; valarray<double> rL_inG = r(Lune, y) - rG; valarray<double> rA_inG = r(Art, y) - rG;
        valarray<double> vT_inG = v(Terre, y) - vG; valarray<double> vL_inG = v(Lune, y) - vG; valarray<double> vA_inG = v(Art, y) - vG;

        valarray<double> y_inG(numbodies * dimension * 2);

        y_inG[slice(0, dimension, 1)] = rT_inG; y_inG[slice(2, dimension, 1)] = rL_inG; y_inG[slice(4, dimension, 1)] = rA_inG;
        y_inG[slice(6, dimension, 1)] = vT_inG; y_inG[slice(8, dimension, 1)] = vL_inG; y_inG[slice(10, dimension, 1)] = vA_inG;

        return y_inG;
    }

    double EmecArt_inG() const { //Energie mécanique d'Artemis
        valarray<double> y_inG = in_G();
        return 0.5 * m[Art] * pow(norm(v(Art, y_inG)), 2) - G * m[Terre] * m[Art] / norm(r(Art, y_inG) - r(Terre, y_inG)) - G * m[Lune] * m[Art] / norm(r(Art, y_inG) - r(Lune, y_inG));
    } 

    double p(const valarray<double>& v, double m) const { //quantité de mouvement normée
        return m*norm(v);
    }

    valarray<double> acc_Artemis_inG() const {
        valarray<double> y_inG = in_G();

        valarray<double> rT = r(Terre, y_inG);
        valarray<double> vT = v(Terre, y_inG);
        valarray<double> rL = r(Lune, y_inG);
        valarray<double> rA = r(Art, y_inG);
        valarray<double> vA = v(Art, y_inG);

        return acc_grav(rA, rT, m[Terre]) + Ft(rA-rT, vA-vT)/m[Art] + acc_grav(rA, rL, m[Lune]);
    }

    double Pt_Art_inG() const {
        valarray<double> y_inG = in_G();

        valarray<double> rT = r(Terre, y_inG);
        valarray<double> rA = r(Art, y_inG);
        valarray<double> vrel = v(Art, y_inG) - v(Terre, y_inG);
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
        valarray<valarray<double>> tab = {r(Terre, y_), v(Terre, y_), r(Lune, y_), v(Lune, y_), r(Art, y_), v(Art, y_)};

        for (int i = 0; i<=2; i++) {
            valarray<double> v_after = acc_grav(tab[i*2], tab[((i+1)%3*2)], m[(i+1)%3]) + acc_grav(tab[i*2], tab[((i+2)%3)*2], m[(i+2)%3]); // calcule v_after correspondant à l'indice i

            if (i == Art) {
                v_after += Ft(tab[i*2]-tab[0], tab[i*2+1]-tab[1])/m[Art];
            }

            df[ix(i)] = tab[i*2+1][0];
            df[iy(i)] = tab[i*2+1][1];
            df[ivx(i)] = v_after[0];
            df[ivy(i)] = v_after[1];
        }
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
                    valarray<double> y1 = rk4Step(dt, y);
                    valarray<double> y2 = rk4Step(dt*0.5, rk4Step(dt*0.5, y));

                    double d = 0;

                    for(size_t i = 0; i < y.size(); ++i){
                        d += pow(y1[i]-y2[i], 2);
                        //error += pow(y1[i] - y2[i], 2);
                    }
                    d = sqrt(d);
                    double f = 0.95;
                    double dt_new;

                    if (d < 1e-6*epsilon) {
                        dt_new = dt*2;
                    }
                    else {
                        dt_new = dt * pow(epsilon/d, 1.0/5.0);
                    }
                    if (d <= epsilon) {
                        y = y2;
                        t += dt;
                        dt = dt_new;
                        step_accepted = true;
                    } else {
                        dt = f*dt_new;
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
