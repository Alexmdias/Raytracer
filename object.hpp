#ifndef OBJECT_H
#define OBJECT_H

#include "basic.hpp"

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>

// Le type d'une "liste de paramètres", e.g. une map de strings vers des listes de nombres.
typedef std::map<std::string, std::vector<double> > ParamList;


// Une classe pour encapsuler tous les paramètres pour un matériau.
class Material
{
public:
	// Constructeurs
	Material() {};
	Material(ParamList &params) { init(params); }

	void init(ParamList &params)
	{
		#define SET_VECTOR(_name) _name = params[#_name];
		SET_VECTOR(ambient)
		SET_VECTOR(diffuse)
		SET_VECTOR(specular)
		SET_VECTOR(emission)

		#define SET_FLOAT(_name) _name = params[#_name].size() ? params[#_name][0] : 0;
		SET_FLOAT(shininess)
		SET_FLOAT(reflect)
	}

	// Couleurs ambiante/diffuse/speculaire/emissive.
	Vector ambient;
	Vector diffuse;
	Vector specular;
	Vector emission;

	// Facteur de "brillance" (specular exponent).
	double shininess;

	// Coefficient de reflection coefficient [0 -> aucune reflection, 1 -> reflection totale]
	// Toutes les valeurs entre sont mélangées par ce facteur avec la couleur de la surface.
	double reflect;
};


// Classe abstraite de base pour les objets.
class Object 
{
public:
    Matrix transform;   // Transformation de l'espace global à l'espace de l'objet (global --> local).
    Matrix i_transform; // Transformation de l'espace de l'objet à l'espace global (local --> global).
    Matrix n_transform; // Trasnformation vers l'espace global pour les normales (local --> global).

    // Mets en place les 3 transformations à partir de la transformation (global-vers-objet) donnée.
    void setup_transform(Matrix const &m)
	{
		transform = m;
		m.invert(i_transform);
		n_transform = i_transform.transpose();
	}

    // Intersecte l'objet avec le rayon donné, dans l'espace global.
	// Retourne true s'il y a eu une intersection, hit est alors mis à jour avec les paramètres.
    bool intersect(Ray ray, Intersection &hit) const;

    // Intersecte l'objet avec le rayon donné, dans l'espace de l'objet.
	// Cette fonction est spécifique à chaque sous-type d'objet.
	// Retourne true s'il y a eu une intersection, hit est alors mis à jour avec les paramètres.
    virtual bool localIntersect(Ray const &ray, Intersection &hit) const = 0;

    Material material; // Matériau de l'objet.
};


// Une sphère centrée autour de son origine locale, avec un certain rayon.
class Sphere : public Object 
{
public:
    double radius;
    
    bool localIntersect(Ray const &ray, Intersection &hit) const;
};


// Un plan à l'origine, utilisant Z+ comme normale dans l'espace objet.
class Plane : public Object 
{
public:
    virtual bool localIntersect(Ray const &ray, Intersection &hit) const;
};


// Un cone le long de l'axe Z+, limité le long de Z par zMin et zMax,
// avec les rayons radius1 et radius2.
class Conic : public Object 
{
public:
    double radius1, radius2;
    double zMin, zMax;

    bool localIntersect(Ray const &ray, Intersection &hit) const;
};


// Une classe pour représenter un seul sommet d'un polygone. Les entiers stockés
// sont des indices pour les vecteurs positions/texCoords/normals/colors de
// l'objet auquel le sommet appartient.
class Vertex {
public:
	// Indices dans les vecteurs positions, texCoods, normals, et colors.
	int pi, ti, ni, ci;

	Vertex() : pi(-1), ti(-1), ni(-1), ci(-1) {}

	Vertex(int pi, int ti, int ni, int ci) :
		pi(pi),
		ti(ti),
		ni(ni),
		ci(ci)
	{}
};


class Triangle
{
public:
	Vertex v[3];

	Triangle(Vertex const &v0, Vertex const &v1, Vertex const &v2)
	{
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
	}

	Vertex& operator[](int i) { return v[i]; }
	const Vertex& operator[](int i) const { return v[i]; }
};


class Mesh : public Object {
public:
	// Contenant pour les positions, coordonnées de texture, normales et couleurs. Recherche par indice.
	std::vector<Vector> positions;
	std::vector<Vector> texCoords;
	std::vector<Vector> normals;
	std::vector<Vector> colors;

	// Les triangles sont des triplets de sommets.
	std::vector<Triangle> triangles;

	// Boite englobante de l'objet.
	Vector bboxMin, bboxMax;

	// Lis les données OBJ d'un fichier donné.
	bool readOBJ(std::string const &filename)
	{
		// Essaie d'ouvrir le fichier.
		std::ifstream file(filename.c_str());
		if (!file.good()) {
			std::cerr << "Unable to open OBJ file \"" << filename << "\"" << std::endl;
			return false;
		}

		// Continue de récupérer les codes opérationnel et de les analyser. Nous supposons
		// qu'il n'y a qu'une opération par ligne.
		while (file.good()) {

			std::string opString;
			std::getline(file, opString);

			std::stringstream opStream(opString);
			std::string opCode;
			opStream >> opCode;

			// Saute les lignes blanches et les commentaires.
			if (!opCode.size() || opCode[0] == '#') {
				continue;
			}

			// Ignore les groupes.
			if (opCode[0] == 'g') {
				std::cerr << "ignored OBJ opCode '" << opCode << "'" << std::endl;
			} // Données de sommet.
			else if (opCode[0] == 'v') {

				// Lis jusqu'à 4 doubles.
				Vector vec;
				for (int i = 0; opStream.good() && i < 4; i++) {
					opStream >> vec[i];
				}

				// Stocke cette donnée dans le bon vecteur.
				switch (opCode.size() > 1 ? opCode[1] : 'v') {
				case 'v':
					positions.push_back(vec);
					break;
				case 't':
					texCoords.push_back(vec);
					break;
				case 'n':
					normals.push_back(vec);
					break;
				case 'c':
					colors.push_back(vec);
					break;
				default:
					std::cerr << "unknown vertex type '" << opCode << "'" << std::endl;
					break;
				}
			} // Un polygone (ou face).
			else if (opCode == "f") {
				std::vector<Vertex> polygon;
				// Limite à 4 sommets, puisque nous ne gérons que les triangles ou les quads.
				for (int i = 0; opStream.good() && i < 4; i++) {

					// Récupère la spécification complète d'un sommet.
					std::string vertexString;
					opStream >> vertexString;

					if (!vertexString.size()) {
						break;
					}

					// Analyse le sommet en un set d'indices pour les positions, coordonnées de tetxure,
					// normales et couleurs, respectivement.
					std::stringstream vertexStream(vertexString);
					std::vector<int> indices;
					for (int j = 0; vertexStream.good() && j < 4; j++) {
						// Saute les slashes.
						if (vertexStream.peek() == '/') {
							vertexStream.ignore(1);
						}
						int index;
						if (vertexStream >> index)
							indices.push_back(index);
					}

					// Transforme les données récupérées en un véritable sommet, et l'ajoute au polygone.
					if (indices.size()) {
						indices.resize(4, 0);
						polygon.push_back(Vertex(
							indices[0] - 1,
							indices[1] - 1,
							indices[2] - 1,
							indices[3] - 1
							));
					}

				}

				// On n'accepte que les triangles...
				if (polygon.size() == 3) {
					triangles.push_back(Triangle(polygon[0],
						polygon[1],
						polygon[2]));
				} // ... et les quads...
				else if (polygon.size() == 4) {
					// ... mais on les décompose en triangle.
					triangles.push_back(Triangle(polygon[0],
						polygon[1],
						polygon[2]));
					triangles.push_back(Triangle(polygon[0],
						polygon[2],
						polygon[3]));
				}

				// Tous les autres codes op sont ignorés.
			}
			else {
				std::cerr << "unknown opCode '" << opCode << "'" << std::endl;
			}
		}

		updateBBox();

		return true;
	}

	// Construit la boite englobante des sommets.
	// Ceci doit être appelé après que les données du maillages aient été
	// initialisées et avant que le lancer de rayons ne commence.
	void updateBBox()
	{
		bboxMin = Vector(DBL_MAX, DBL_MAX, DBL_MAX);
		bboxMax = Vector(-DBL_MAX, -DBL_MAX, -DBL_MAX);
		for (std::vector<Vector>::iterator pItr = positions.begin();
			pItr != positions.end(); ++pItr) {
			Vector const& p = *pItr;
			if (p[0] < bboxMin[0]) bboxMin[0] = p[0];
			if (p[0] > bboxMax[0]) bboxMax[0] = p[0];
			if (p[1] < bboxMin[1]) bboxMin[1] = p[1];
			if (p[1] > bboxMax[1]) bboxMax[1] = p[1];
			if (p[2] < bboxMin[2]) bboxMin[2] = p[2];
			if (p[2] > bboxMax[2]) bboxMax[2] = p[2];
		}
	}

	// Intersections !
	bool localIntersect(Ray const &ray, Intersection &hit) const;

private:
	// Calcule le résultat de l'équation de ligne implicite en 2D
	// pour un point donné et une ligne avec les points de terminaison donnés.
	double implicitLineEquation(double p_x, double p_y,
		double e1_x, double e1_y,
		double e2_x, double e2_y) const;

	// Trouve le point d'intersection entre le rayon donné et le maillage triangulaire.
	// Renvoie true ssi une intersection existe, et remplit les données de
	// la structure hit avec les bonnes informations.
	bool intersectTriangle(Ray const &ray,
		Triangle const &tri,
		Intersection &hit) const;
};


#endif
