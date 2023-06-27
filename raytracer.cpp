#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>

#include "raytracer.hpp"
#include "image.hpp"


void Raytracer::render(const char *filename, const char *depth_filename,
                       Scene const &scene)
{
    // Alloue les deux images qui seront sauvegardées à la fin du programme.
    Image colorImage(scene.resolution[0], scene.resolution[1]);
    Image depthImage(scene.resolution[0], scene.resolution[1]);
    
    // Crée le zBuffer.
    double *zBuffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        zBuffer[i] = DBL_MAX;
    }
	//direction
	Vector direction = (scene.camera.center - scene.camera.position).normalized();

	Vector u = (scene.camera.up).normalized(); //u = upvector
	Vector v = direction.cross(u).normalized(); //v = right vector

	double distance = scene.camera.zNear; //position du znear de la camera

	//plan image
	double planeTop = distance * tan(deg2rad(scene.camera.fovy) / 2); //t
	double planeRight = planeTop * scene.camera.aspect; //r
	double planeLeft = -planeRight; //l
	double planeBottom = -planeTop; //b


	double deltaU = (planeRight - planeLeft) / scene.resolution[0]; // (r-l)/Nx
	double deltaV = (planeTop - planeBottom) / scene.resolution[1]; // (t-b)/Ny

	Vector planO = (scene.camera.position + (scene.camera.zNear * direction)) + (planeLeft * v) + (planeBottom * u); // O

    // Itère sur tous les pixels de l'image.
    for(int y = 0; y < scene.resolution[1]; y++) {
        for(int x = 0; x < scene.resolution[0]; x++) {

            // Génère le rayon approprié pour ce pixel.
			Ray ray;
			if (scene.objects.empty())
			{
				// Pas d'objet dans la scène --> on rend la scène par défaut.
				// Pour celle-ci, le plan de vue est à z = 640 avec une largeur et une hauteur toute deux à 640 pixels.
				ray = Ray(scene.camera.position, (Vector(-320, -320, 640) + Vector(x + 0.5, y + 0.5, 0) - scene.camera.position).normalized());
			}
			else
			{

				//Pij = O + (i+0.5) * deltaU*v + (j+0.5) * deltaV*u
				Vector pxy = planO + (x + 0.5) * deltaU * v + (y + 0.5) * deltaV * u;

				Vector rayDirection = (pxy - scene.camera.position).normalized();

				ray = Ray(scene.camera.position, rayDirection);
				
			}

            // Initialise la profondeur de récursivité du rayon.
            int rayDepth = 0;
           
            // Notre lancer de rayons récursif calculera la couleur et la z-profondeur.
            Vector color;

            // Ceci devrait être la profondeur maximum, correspondant à l'arrière plan.
            // NOTE : Ceci suppose que la direction du rayon est de longueur unitaire (normalisée)
			//        et que l'origine du rayon est à la position de la caméra.
            double depth = scene.camera.zFar;

            // Calcule la valeur du pixel en lançant le rayon dans la scène.
            trace(ray, rayDepth, scene, color, depth);

            // Test de profondeur
            if(depth >= scene.camera.zNear && depth <= scene.camera.zFar && 
                depth < zBuffer[x + y*scene.resolution[0]]) {
                zBuffer[x + y*scene.resolution[0]] = depth;

                // Met à jour la couleur de l'image (et sa profondeur)
                colorImage.setPixel(x, y, color);
                depthImage.setPixel(x, y, (depth-scene.camera.zNear) / 
                                        (scene.camera.zFar-scene.camera.zNear));
            }
        }

		// Affiche les informations de l'étape
		if (y % 100 == 0)
		{
			printf("Row %d pixels finished.\n", y);
		}
    }

	// Sauvegarde l'image
    colorImage.writeBMP(filename);
    depthImage.writeBMP(depth_filename);

	printf("Ray tracing finished with images saved.\n");
    delete[] zBuffer;

}


bool Raytracer::trace(Ray const &ray, 
                 int &rayDepth,
                 Scene const &scene,
                 Vector &outColor, double &depth)
{
    // Incrémente la profondeur du rayon.
    rayDepth++;

    // - itérer sur tous les objets en appelant   Object::intersect.
    // - ne pas accepter les intersections plus lointaines que la profondeur donnée.
    // - appeler Raytracer::shade avec l'intersection la plus proche.
    // - renvoyer true ssi le rayon intersecte un objet.
	if (scene.objects.empty())
	{
		// Pas d'objet dans la scène --> on rend la scène par défaut :
		// Par défaut, un cube est centré en (0, 0, 1280 + 160) avec une longueur de côté de 320, juste en face de la caméra.
		// Test d'intersection :
		double x = 1280 / ray.direction[2] * ray.direction[0] + ray.origin[0];
		double y = 1280 / ray.direction[2] * ray.direction[1] + ray.origin[1];
		if ((x <= 160) && (x >= -160) && (y <= 160) && (y >= -160))
		{
			// S'il y a intersection :
			Material m; m.emission = Vector(16.0, 0, 0); m.reflect = 0; // seulement pour le matériau par défaut ; vous devrez utiliser le matériau de l'objet intersecté
			Intersection intersection;	// seulement par défaut ; vous devrez passer l'intersection trouvée par l'appel à Object::intersect()
			outColor = shade(ray, rayDepth, intersection, m, scene);
			depth = 1280;	// la profondeur devrait être mise à jour dans la méthode Object::intersect()
		}
	}
	else
	{

		Intersection intersect;
		intersect.depth = depth;
		double mini = depth;
		for (int i = 0; i < scene.objects.size(); i++) { //on itere sur chaque objets
			if ((scene.objects[i]->intersect(ray, intersect)) && (intersect.depth <= mini)) {
				mini = intersect.depth;
				outColor = shade(ray, rayDepth, intersect, scene.objects[i]->material, scene); //necessaire pour shade
			}
		}
		depth = mini;
	}

    // Décrémente la profondeur du rayon.
    rayDepth--;

    return false; 
}


Vector Raytracer::shade(Ray const &ray,
                 int &rayDepth,
                 Intersection const &intersection,
                 Material const &material,
                 Scene const &scene)
{

	Vector diffuse(0);
	Vector ambient(0);
	Vector specular(0);
	for (auto lightIter = scene.lights.begin(); lightIter != scene.lights.end(); lightIter++)
	{
		

		Vector n = intersection.normal.normalized(); // normale
		Vector l = (lightIter->position - intersection.position).normalized(); //light dir
		Vector r = (-2 * (l.dot(n)) * n + l).normalized(); //reflect dir
		Vector v = ray.direction.normalized(); //viewdir

		//les classiques intensités
		Vector I_a = (lightIter->ambient * material.ambient); //ambiant
		Vector I_d = (max(n.dot(l), 0.0) * material.diffuse * lightIter->diffuse); //diffuse
		Vector I_s = (pow(max(v.dot(r), 0.0), material.shininess) * material.specular * lightIter->specular); //specular

		//attenuation = 1.0 / (a0 + a1 * d + a2 * d * d)
		double d = (lightIter->position - intersection.position).length(); //distance lumiere
		double attenuation = 1/(lightIter->attenuation[0] + lightIter->attenuation[1] * d + lightIter->attenuation[2] * d * d);
		
		
		ambient += I_a;
		specular += I_s * attenuation;
		diffuse += I_d * attenuation;


		Intersection shadowIntersect;
		Ray shadowRay = Ray(intersection.position + (l * 0.00001), l);
		shadowIntersect.depth = (lightIter->position - intersection.position).length();

		//si c'est un ombre, on enleve la lumiere specular et diffuse
		for (int i = 0; i < scene.objects.size(); ++i) {
			if (scene.objects[i]->intersect(shadowRay, shadowIntersect)) { 
				specular += -I_s * attenuation;
				diffuse += -I_d * attenuation;
				break;
			}
		}		

	}

	Vector reflectedLight(0);
	if ((!(ABS_FLOAT(material.reflect) < 1e-6)) && (rayDepth < MAX_RAY_RECURSION))
	{
		// Calculez la couleur réfléchie en utilisant trace() de manière récursive.
		Vector n = intersection.normal.normalized(); // normale
		Vector l = ray.direction.normalized(); //light dir
		Vector r = (-2 * (l.dot(n)) * n + l).normalized(); //reflect dir
		Ray reflect = Ray((intersection.position + (r * 0.00001)), r); // r * 0.00001 important
		double maxdepth = 99999999999; // big nombre

		trace(reflect, rayDepth, scene, reflectedLight, maxdepth); //classique récurssivité 
	}
	return material.emission + ambient + diffuse + specular + material.reflect * reflectedLight;
}