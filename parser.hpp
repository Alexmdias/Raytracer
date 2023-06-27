#ifndef PARSER_H
#define PARSER_H

#include "scene.hpp"

#include <climits>
#include <deque>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>


// Les différents types de token pouvant être lexés
enum TokenType {
	STRING,
	NUMBER,
	NAME,
	ARRAY_BEGIN,
	ARRAY_END,
	END_OF_FILE,
	ERROR
};


// Une classe pour représenter un token qui a été lu
class Token 
{
public:
	TokenType type;

	// Variables pour les différents types de données.
	double number;
	std::string string;

	// Quelques constructeurs pour assigner directement des valeurs.
	Token(TokenType type) : type(type) {}
	Token(double value) : type(NUMBER), number(value) {}
	Token(TokenType type, std::string value) : type(type), string(value) {}

	// Opérateur d'égalité.
	bool operator==(Token const &other) const;
};


// Permet l'écriture d'un token directement dans le flux de sortie.
std::ostream& operator<<(std::ostream &out, Token const &token);


// La classe du lexer.
class Lexer 
{
private:
	// Le flux d'entrée.
	std::istream *_input;

	// Fonction pour lire le flux d'entrée et retourner le token suivant.
	Token _processStream();

	// Un buffer temporaire pour les tokens qui n'ont pas encore été analysés.
	std::deque<Token> _buffer;

public:
	// Constructeur. Un flux d'entrée doit être fourni.
	Lexer(std::istream *input) : _input(input) {}

	// Regarde le prochain token mais ne le consomme pas.
	Token peek(unsigned int index = 0);

	// Accède au token suivant.
	Token next();

	// Passe un certain nombre de tokens.
	void skip(unsigned int count = 1);

	// Les fonctions suivantes produiront une exception std::string si
	// elles ne peuvent pas fonctionner comme demandé.

	// Récupère un nom de commande.
	std::string getName();

	// Récupère une liste de nombres. Min/max font référence
	// à la taille requise de la liste.
	std::vector<double> getNumbers(unsigned int min = 0,
		unsigned int max = UINT_MAX);

	// Récupère un unique nombre.
	double getNumber();

	// Récupère une unique string.
	std::string getString();

	// Récupère une liste de paramètres (i.e. strings mappés à des listes de nombres).
	// Min/max s'appliquent à chaque liste, tout comme getNumbers().
	ParamList getParamList(unsigned int min = 0,
		unsigned int max = UINT_MAX);
};


class Parser 
{
private:
    Lexer lexer; // Le lexer utilisé pour séparer le fichier en tokens.
    
    std::vector<Matrix> transformStack;  // Pile de transformations.

    // Les fonctions suivantes analysent toutes les commandes qui peuvent être
	// trouvées dans un fichier .ray.
    void parseDimensions();
    void parsePerspective();
    void parseLookAt();
    void parseMaterial();

    void parsePushMatrix();
    void parsePopMatrix();
    void parseTranslate();
    void parseRotate();
    void parseScale();

    void parseSphere();
    void parsePlane();
    void parseMesh();
    void parseConic();

    void parsePointLight();

    // Analyse les parties communes de chaque objet et met en place les objets dans la scène.
    void finishObject(Object *obj);

public:
    Scene scene; // La scène qui sera créée pendant l'analyse.

    Parser(char const * filename) : lexer(new std::ifstream(filename)) {}
    Parser(std::istream *input) : lexer(input) {}

    // Analyse le fichier ou le flux passé au constructeur.
	// Sauvegarde le résultat dans une scène.
	// Retourne false sur un échec ; erreur écrite dans std::cerr.
    bool parse();
};


#endif
