#ifndef Particle_H
#define Particle_H
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>

using namespace std;

class Particle
{
 public:
  double mass, energy, theta, phi, beta;
  TLorentzVector LorentzVector;
  
  // Constructeurs
  Particle();
  Particle(TVector3 momentum, int pid);
  Particle(TVector3 momentum, TVector3 position, int pid, int chi2pid, int status, int charge);
  Particle(TVector3 momentum, TVector3 position, int pid, int chi2pid, int status, int charge,double beta);  
  //Accesseurs et mutateurs
  void setmomentum(TVector3 momentum);
  void setposition(TVector3 position);
  void setpid(int pid);
  void setchi2pid(int chi2pid);
  void setstatus(int status);
  void setcharge(int charge);
  TVector3 getmomentum() const;
  TVector3 getposition() const;
  int getpid() const;
  int getchi2pid() const;
  int getstatus() const;
  int getcharge() const;

  //Other functions
  void compute_energy();
  void compute_theta();
  void compute_phi();
  void get_mass();
  void compute_LorentzVector();
  void construct_all_properties();
  
 private: 
  TVector3 momentum;
  TVector3 position;
  int pid, chi2pid, status, charge;
};

#endif
