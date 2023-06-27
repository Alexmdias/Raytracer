#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "scene.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>
using namespace std;

// Profondeur de récursion maximum des rayons réfléchis/réfractés/etc. pour chaque pixel.
#define MAX_RAY_RECURSION 10


class Raytracer 
{
public:
    // Rend la scène donnée en utilisant le lancer de rayon.
    // Sauvegarde l'image résulta et la carte de profondeur dans les fichiers indiqués.
    static void render(const char *filename, const char *depth_filename, 
                       Scene const &scene);

private:
    // Lancer un rayon dans la scène et calcule ses valeurs de couleur/profondeur.
    // Paramètres :
    //   rayon -- le rayon lancé à travers la scène
    //   rayDepth -- la profondeur de récursion du rayon actuellement lancé
    //   scene -- la scène dans laquelle le rayon est lancé
    //   outColor -- rempli avec la couleur calculée pour le rayon à la sortie
    //   depth -- rempli avec la z-profondeur de la plus proche intersection, si elle existe
    //      NOTE : La profondeur passée avec ce paramètre est traitée comme une limite
    //             supérieure. Les objets plus lointain que cette profondeur doivent être
    //             ignorés.
    // Renvoie true ssi le rayon intersecte un objet de la scène.
    static bool trace(Ray const &ray, 
                      int &rayDepth,
                      Scene const &scene,
                      Vector &outColor, double &depth);

    // Calcule l'ombrage pour un point et une normale donnés.
    // Utilise le matériau donné ainsi que les sources de lumières
    // et les autres objets de la scène pour les ombres et réflections.
    // Params:
    // Paramètres :
    //   rayon -- le rayon lancé à travers la scène
    //   rayDepth -- la profondeur de récursivité du rayon actuellement lancé
    //   intersection -- informations à propos de l'intersection rayon-surface
    //   material -- le matériau de l'objet au point d'intersection
    //   scene -- la scène dans laquelle le rayon est lancé
    // Renvoie la couleur calculée.
	static Vector shade(Ray const &ray,
                        int &rayDepth,
                        Intersection const &intersection,
                        Material const &material,
                        Scene const &scene);
};


#endif
