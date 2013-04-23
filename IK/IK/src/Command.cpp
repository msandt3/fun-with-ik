#ifndef __COMMAND_H__
#include "Command.h"
#endif //__COMMAND_H__

#ifndef __C3DFILEINFO_H__
#include "C3dFileInfo.h"
#endif	//__C3DFILEINFO_H__

#ifndef __ARTICULATEDBODY_H__
#include "ArticulatedBody.h"
#endif	//__ARTICULATEDBODY_H__

#ifndef RealTimeIKui_h
#include "RealTimeIKui.h"
#endif //RealTimeIKui_h

#ifndef __PHYLTERFLGLWINDOW_H__
#include "PhylterGLWindow.h"
#endif	//__PHYLTERFLGLWINDOW_H__

#ifndef	__PHOWARDDATA_H__
#include "PhowardData.h"
#endif

#ifndef __TRANSFORM_H__
#include "Transform.h"
#endif	//__TRANSFORM_H__



int readSkelFile( FILE* file, ArticulatedBody* skel );

extern RealTimeIKUI *UI;

void LoadModel(void *v)
{
  char *params = (char*)v;
  if(!params){
    params = (char*)fl_file_chooser("Open File?", "{*.skel}", "../src/skels" );
  }

  if(!params)
    return;

  FILE *file = fopen(params, "r");
    
  if(file == NULL){
    cout << "Skel file does not exist" << endl;
    return;
  }

  ArticulatedBody *mod = new ArticulatedBody();
  UI->mData->mModels.push_back(mod);
  UI->mData->mSelectedModel = mod;

  readSkelFile(file, mod);
  UI->CreateDofSliderWindow();

  mod->InitModel();
  UI->mGLWindow->mShowModel = true;
  UI->mShowModel_but->value(1);
  UI->mGLWindow->refresh();
  
  cout << "number of dofs in model: " << UI->mData->mModels[0]->GetDofCount() << endl;
}

void Solution(void *v)
{
    cout << "TODO: Solve inverse kinematics problem" << endl;
    bool test = UI->mData->mSelectedModel->mLimbs[0]->mTransforms[0]->IsDof();
	cout << test;

	TMat jacobian = TMat();
	jacobian.SetSize(3,9);
	//Calculate C
	Vec3d c = CalculateC();
	float error = 0.1f;

	while(CalculateFQ(c) < error){
		//Calculate Jacobian
		computeJ(jacobian);
		//Calculate Jacobian Transpose
		TMat transposeJ=trans(jacobian);
		//Calculate partial F/ partial q
		Vecd pFpq = 2*(transposeJ*c);
		cout << pFpq;
		//update q term
		c = CalculateC();
	}
		
}

void Exit(void *v)
{
  exit(0);
}

void LoadC3d(void *v)
{
  if(!UI->mData->mSelectedModel){
    cout << "Load skeleton first";
    return;
  }
  char *params = (char*)v;
  if(!params){
    params = fl_file_chooser("Open File?", "{*.c3d}", "mocap/" );
  }

  if(!params)
    return;
  
  char *c3dFilename = new char[80];
  
  // load single c3d file
 
  C3dFileInfo *openFile = new C3dFileInfo(params);
  openFile->LoadFile();
  UI->mData->mSelectedModel->mOpenedC3dFile = openFile;
  cout << "number of frames in c3d: " << openFile->GetFrameCount() << endl;

  UI->InitControlPanel();
  UI->mGLWindow->mShowConstraints = true;
  UI->mShowConstr_but->value(1);
}
Vec3d CalculateC(){
	Marker* mark=UI->mData->mSelectedModel->mHandleList[0];
	Vec3d pBar=UI->mData->mSelectedModel-> mOpenedC3dFile->GetMarkerPos(0,0); 
	Vec3d temp=mark->mGlobalPos-pBar;
	return temp;
}

TMat computeJ(TMat Jacobian){
	Marker* mark=UI->mData->mSelectedModel->mHandleList[0];
	TransformNode* node=UI->mData->mSelectedModel->mLimbs[mark->mNodeIndex];
	Mat4d parent=node->mParentTransform;
	Mat4d T=node->mTransforms[0]->GetTransform();
	Mat4d partRpartQ=node->mTransforms[1]->GetDeriv(0);
	Mat4d Rq= node->mTransforms[2]->GetTransform();
	Vec4d offset=Vec4d(mark->mOffset,1);
	Vec4d Ji=parent*T*partRpartQ*Rq*offset;
	int column=node->mTransforms[1]->GetDof(0)->mId;
	Jacobian[column]=Ji;
	return Jacobian;
}

float CalculateFQ(Vec3d c){
	float length = len(c);
	float error = length * length;
	return error;
}
