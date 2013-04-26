#ifndef __COMMAND_H__
#define __COMMAND_H__
 
#include <fstream>
#include "vl/VLd.h"
#include <vector>
#include "RealTimeIKui.h"
 
typedef void (*Command)(void*);
Vec3d CalculateC();
void computeJ();
void LoadModel(void*);
void LoadC3d(void*);
void Exit(void*);
void Solution(void*);
float CalculateFQ(Vec3d c);
#endif