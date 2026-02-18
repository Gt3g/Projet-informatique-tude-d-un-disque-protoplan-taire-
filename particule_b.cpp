#include <iostream>
#include <cmath>
#include "particule_b.h"


using namespace std;

void Particule::init(double x , double y , double z , double v_x, double v_y , double v_z , double m , double q){
  this->coord_x = x;
  this->coord_y = y;
  this->coord_z = z;
  this->vit_x = v_x;
  this->vit_y = v_y;
  this->vit_z = v_z;
  this->masse = m;
  this->charge = q;

}

void Particule::afficher(){
  cout<< " x=" << this->coord_x << " y=" << this->coord_y << "z="<< this->coord_x <<endl;
  cout<< " v_x=" << this->vit_x << " v_y=" << this->vit_y << "v_z="<< this->vit_x <<endl;
  cout<< " m=" << this->masse  <<endl;
  cout<< " q=" << this->charge  <<endl;


}

double Particule::dist_origine(){
  return sqrt((this->coord_x)*(this->coord_x) + (this->coord_y)*(this->coord_y) + (this->coord_z)*(this->coord_z)
);
}

double Particule::dist_axe(){
  return sqrt((this->coord_x)*(this->coord_x) + (this->coord_y)*(this->coord_y) 
);
}
