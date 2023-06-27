#include "parser.hpp"
#include "raytracer.hpp"

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

int main(int argc, char **argv)
{
	string scenefile; // Si non spécifié, le programme rendra la scène par défaut qui ne contient qu'un cube.
	string outfilename("scenes/default_output.bmp"); // NOTE : on ne supporte que les sorties bmp
	string outfilename_depth("scenes/default_depth_output.bmp");

    // Non de fichier de scène spécifié
	if (argc > 1)
	{
		scenefile = string(argv[1]);
		int dot = scenefile.find_last_of('.');
		outfilename = scenefile.substr(0, dot) + "_output.bmp";
		outfilename_depth = scenefile.substr(0, dot) + "_depth_output.bmp";
	}

    std::cout << "Rendering " << (scenefile.empty() ? "default scene" : scenefile) << std::endl;
    std::cout << "Output to " << outfilename << std::endl;
   
	Raytracer raytracer;
	if (scenefile.empty())
	{
		// Rend la scène par défaut
		raytracer.render(
			outfilename.c_str(),
			outfilename_depth.c_str(),
			Scene()
			);
	}
	else
	{
		// Analyse le fichier de la scène
		Parser parser(new std::ifstream(scenefile.c_str()));
		if (!parser.parse()) {
			puts("Scene file can't be parsed. Use default scene.");
			raytracer.render(
				outfilename.c_str(),
				outfilename_depth.c_str(),
				Scene()
				);
		}
		else
		{
			// Rend la scène donnée avec le lancer de rayon
			raytracer.render(
				outfilename.c_str(),
				outfilename_depth.c_str(),
				parser.scene
				);
		}
	}
	// Décommentez si vous utilisez Visual Studio
	//system("pause");

	return 0;
}

