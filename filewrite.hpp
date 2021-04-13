#ifndef FILEWRITE_H
#define FILEWRITE_H
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "mesh.hpp"
#include "fields.hpp"
#include <cstring>
#include <sstream>
#include <string>
using namespace std;
using std::vector;
using std::string;


class fileWriter
{

public:
	fileWriter();

	virtual ~fileWriter();

	void writeUVWP(string&, int , Mesh&, Fields::vec3dField&, Fields::vec3dField&, 
			Fields::vec3dField&, Fields::vec3dField&);

	// int is the time step , vec filds - U, V, W, P, string - filename
};

#endif
