#ifndef __COMMAND_H__
#define __COMMAND_H__
 
#include <fstream>
#include "vl/VLd.h"
#include <vector>
#include "RealTimeIKui.h"
 
typedef void (*Command)(void*);
void CalculateC();
void computeJ();
void LoadModel(void*);
void LoadC3d(void*);
void Exit(void*);
void Solution(void*);
double CalculateFQ(int frame);
double CalculateFQ(int handle, int frame);

Vec3d CalculateC(int handle, int frame);
Vecd CalculateCVec(int handle, int frame);
Vecd CalculateWeights(int frame);
bool KeepGoing();
void ComputeJ(int handle);
#endif