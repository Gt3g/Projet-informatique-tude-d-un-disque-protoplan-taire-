# Projet-informatique-tude-d-un-disque-protoplan-taire-
Projet info Ludivine Baptiste Nils


DESCRIPTION PROGRAMMES:

-particules.h/.cpp : classe particule avec position (x,y,z), vitesse(vx,vy,vz), masse et charge(pas utilisé,utilisé comme moment_cin mais finallement non)

-nuage.h/.cpp : classe nuage (vecteur de particule), initialisation avec(nbr de particule,rayon_particule,masse_part,rayon nuage,masse_étoile,rayon étoile,vitesse de rotation de nuage)
                les part initialisées dans un volume de rayon ray_nuage exclue du cylindre centré sur l'axe z et de rayon_étoile, vitesse initial = distance à l'axe*vitesse de rotation.

-nuage_solver.h/.cpp : évolution du système en fonction du temps
                       -position et vitesse de chaque part dans une unique liste pour passer dans rk4
                       -système avec les différentes forces (gravit avec l'étoile, frottement selon z (sinon le système se stabilise pas)
                       -résolution rk4 et écriture dans le fichier evolution_nuage

-main_nuage_etoile.cpp : déclaration des constantes et des différents parramètres
                         calcul du temps de simulation en fonction du temps de chute libre
                         définition des paramètre temporelle ( t0, tEnd prop à tff, pas en fraction de tff)
                         initialisation du nuage
                         simulation rk4
                         
