#pragma once

#include <fstream>

int gmshFracGrid (double * xy, const char * fname);
void addFracMaterial (int ncells, const char * fname);
