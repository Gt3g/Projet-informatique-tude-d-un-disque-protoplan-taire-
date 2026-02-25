#include <iostream>
#include <math.h>
#include<stdlib.h>
#include<random>

#include "nuage.h"
#include "particule.h"

using namespace std;



void Nuage::init_random(int n, double R_e, double m_e, double m_p, double q_e, double q_p , double L){

  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<double> dist(-L,L);
  
  Particule p;
  p.init(0,0,0,0,0,0,m_e,m_p);

  vecteur_de_part.push_back(p);

  for (int i=0 ; i<n;i++){

    double x,y,z;
    do{
      x = dist(gen);
      y = dist(gen);
      z = dist(gen);
    }
    while (x*x + y*y + z*z <R_e*R_e);

    Particule p;
    p.init(x,y,z,0,0,0,m_p,q_p);

    vecteur_de_part.push_back(p);
  }

}

void Nuage::affiche_nuage(){
  int i;
  for(i=0;i<vecteur_de_part.size();i++){
    cout<< "part" << i << " ";
    vecteur_de_part[i].afficher();
    cout<<endl;
  }
}
    
  
