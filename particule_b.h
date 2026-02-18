#ifndef PARTICULE_H
#define PARTICULE_H


class Particule{
 public:

  double coord_x;
  double coord_y;
  double coord_z;

  double vit_x;
  double vit_y;
  double vit_z;

  double masse;

  double charge;
  

  void init(double coord_x , double coord_y , double coord_z ,double vit_x , double vit_y , double vit_z , double masse , double charge );
  void afficher();
  double dist_origine();
  double dist_axe();
};

#endif
