#include "Particle.h"
using namespace std;

//----------------------------------------------------------------------------------------Constructeurs--------------------------------------------------------------------//
Particle::Particle() : momentum(), pid()
{}
Particle::Particle(TVector3 momentum, TVector3 position, int pid, int chi2pid, int status, int charge) : momentum(momentum), position(position), pid(pid), chi2pid(chi2pid), status(status), charge(charge)
{}

Particle::Particle(TVector3 momentum, TVector3 position, int pid, int chi2pid, int status, int charge, double beta) : momentum(momentum), position(position), pid(pid), chi2pid(chi2pid), status(status), charge(charge), beta(beta)
{}


Particle::Particle(TVector3 momentum, int pid) : momentum(momentum), pid(pid)
{}

//----------------------------------------------------------------------------------------Setters----------------------------------------------------------------------------------------//
void Particle::setmomentum(TVector3 momentum)
{
this->momentum = momentum;
}

void Particle::setposition(TVector3 position)
{
  this->position = position;
}

void Particle::setpid(int pid)
{
  this->pid = pid;
}

void Particle::setchi2pid(int chi2pid)
{
  this->chi2pid = chi2pid;
}

void Particle::setstatus(int status)
{
  this->status = status;
}

void Particle::setcharge(int charge)
{
  this->charge = charge;
}

//----------------------------------------------------------------------------------------Getters----------------------------------------------------------------------------------------//
TVector3 Particle::getmomentum() const
{
  return this->momentum;
}

TVector3 Particle::getposition() const
{
  return this->position;
}

int Particle::getpid() const
{
  return this->pid;
}

int Particle::getchi2pid() const
{
  return this->chi2pid;
}

int Particle::getstatus() const
{
  return this->status;
}

int Particle::getcharge() const
{
  return this->charge;
}

//----------------------------------------------------------------------------------------Other methods----------------------------------------------------------------------------------------//

void Particle::get_mass()
{
  if(this->pid==11) this->mass=0.51099895000/1000.;
  else if(this->pid==2212) this->mass=938.27208816/1000.;
  else cout << "PID unknown : " << this->pid << endl;
}

void Particle::compute_energy()
{
  double E =  this->momentum.Mag2() + this->mass*this->mass;
  this->energy = TMath::Sqrt(E);
}

void Particle::compute_theta()
{
  this->theta = (180. / TMath::Pi()) * this->momentum.Theta();
}

void Particle::compute_phi()
{
  this->phi = (180. / TMath::Pi()) * this->momentum.Phi();
}

void Particle::compute_LorentzVector()
{
  TLorentzVector V(this->momentum, this->energy);
  this->LorentzVector = V;
}

void Particle::construct_all_properties()
{
  this->get_mass();
  this->compute_energy();
  this->compute_LorentzVector();
  this->compute_phi();
  this->compute_theta();
}
