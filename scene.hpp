#ifndef SCENE_H
#define SCENE_H

#include "object.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>


struct Camera
{
	Camera(void) : fovy(45), aspect(1.0), zNear(1.0), zFar(10000.0),
		position(0.0, 0.0, 0.0), center(0.0, 0.0, 1.0), up(0.0, 1.0, 0.0) {}

	double fovy;
	double aspect;
	double zNear;
	double zFar;
	Vector position;
	Vector center;
	Vector up;
};


// Une classe pour encapsuler tous les paramètres d'une lumière.
class PointLight 
{
public:
	// Constructeurs
	PointLight();
	PointLight(Vector const &position, ParamList &params) : position(position) { init(params); }

	// Initialise les attributs de la lumière avec la liste de paramètres donnée.
	void init(ParamList &params)
	{
		#define SET_VECTOR(_name) _name = params[#_name];
		SET_VECTOR(ambient)
		SET_VECTOR(diffuse)
		SET_VECTOR(specular)
		SET_VECTOR(attenuation)
	}

	Vector position; // Position de la lumière.

	// Couleurs ambiante/diffuse/speculaire de la lumière.
	Vector ambient;
	Vector diffuse;
	Vector specular;

	// Coefficients d'attenuation.
	// attenuation[0] = CONSTANT
	// attenuation[1] = LINEAR
	// attenuation[2] = QUADRATIC
	Vector attenuation;
};


// Une classe qui stocke tous les paramètres, matériaux et objets
// dans une scène que l'on cherche à rendre.
class Scene {
public:
	Scene(void)
	{
		resolution[0] = resolution[1] = 640;
	}
public:
    // Résolution (largeur/hauteur) de l'image de sortie, en pixels.
    int resolution[2]; 

    // La caméra utilisée durant le rendu de la scène.
    Camera camera;
    
    // Map les noms de matériaux aux matériaux eux-mêmes.
    std::map<std::string, Material> materials;

    // Liste des lumières ponctuelles.
    std::vector<PointLight> lights; 

    // Liste des pointeurs vers les objets de la scène.
    // Notez que la classe Object est abstraite, donc les items pointeront réellement
	// vers des objets Spheres, Planes, Mehses, etc.
    std::vector<Object*> objects;
};


#endif
