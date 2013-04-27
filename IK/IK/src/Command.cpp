#ifndef __COMMAND_H__
#include "Command.h"
#endif //__COMMAND_H__
 
#ifndef __C3DFILEINFO_H__
#include "C3dFileInfo.h"
#endif  //__C3DFILEINFO_H__
 
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
bool solve = false;
Vec3d c=Vec3d();
TMat Jacobian = TMat(3,9);
 
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
    //cout << "TODO: Solve inverse kinematics problem" << endl;
    bool test = UI->mData->mSelectedModel->mLimbs[0]->mTransforms[0]->IsDof();
 
	Jacobian.MakeZero();
	//Calculate C
	CalculateC();
	double error = 0.01;
	double alpha = 0.01;
	double FQ = CalculateFQ();
 
	solve=true;
	//while(FQ > error){
		//Calculate Jacobian
		//Jacobian.MakeZero();
		computeJ();
		//Calculate Jacobian Transpose
		TMat transposeJ=trans(Jacobian);
		//Calculate partial F/ partial q
		Vecd pFpq = 2*(transposeJ*c);
		Vecd qNew = Vecd(9);
		//update q term
		for(int i=0; i<pFpq.Elts(); i++){
			double qOld = UI->mData->mSelectedModel->mDofList.GetDof(i);
			//subtract partial from q
			qNew[i] = qOld - alpha * pFpq[i];
			//set the new q value
		}
		UI->mData->mSelectedModel->SetDofs(qNew);
		CalculateC();
		FQ = CalculateFQ();
	//}
	if(FQ > 0.0025)
		Fl::add_timeout(0.001,Solution);
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
 
 void CalculateC(){
	Marker* mark=UI->mData->mSelectedModel->mHandleList[0];
	Vec3d pBar=UI->mData->mSelectedModel-> mOpenedC3dFile->GetMarkerPos(0,0);
	Vec3d handlePos=mark->mGlobalPos;
	c=mark->mGlobalPos-pBar;
	
}
void computeJ(){
 
 
	Marker* mark=UI->mData->mSelectedModel->mHandleList[0];
	//node we're computing partials for
	TransformNode* node = UI->mData->mSelectedModel->mLimbs[mark->mNodeIndex];
	//remaining transformations
	int NeedOffset = 1;
	Vec4d Ji;
	Vec4d u = Vec4d(1,1,1,1);
	Mat4d um = Mat4d(vl_1);
 
	/** While there are still nodes to process **/
	while(node != NULL){
		Mat4d parent = node->mParentTransform;
		//loop over the transforms for this node
 
 
		for(int trans=0; trans<node->mTransforms.size(); trans++){
			Transform* current = node->mTransforms[trans];
			
 
			//determine if the current transform is a dof
			if(current->IsDof()){
				//loop over the DOF's in the transform
 
 
				for(int dof=0; dof<current->GetDofCount(); dof++){
					//compute partial derivative
					Mat4d partial = current->GetDeriv(dof);
					//cout << "Deriv is " << partial << "\n";
					Mat4d Jim = parent;

 
 
					//compute jacobian entry & u (as matrices)
					for(int i=0; i<node->mTransforms.size(); i++){
						if(i == trans){
							Jim *= partial;
						}else{
							Jim *= node->mTransforms[i]->GetTransform();
						}
						um *= node->mTransforms[i]->GetTransform();
					}
 
					
					//multiply by offset -- only if at the foot joint
					if(NeedOffset == 1){
						Ji = Jim * Vec4d(mark->mOffset,1);
						u = um * Vec4d(mark->mOffset,1);
					}
					else{
						u = u * um;
						Ji = Jim * u;
					}
 
 
					//compute row & column of entry
					int column = current->GetDof(dof)->mId;
 
					//set jacobian column
					for(int j=0; j<3; j++){\
						//Jacobian[j][column] = Ji[j];
						Jacobian[j][column] = Ji[j];
					}
 

					
				}//end DOF loop
 
 
 
			}//end dof check
 
 
		}//end transform loop
 
 
		NeedOffset = 0;
		node = node->mParentNode;
	}//end while
}
double CalculateFQ(){
	double length = len(c);
	double error = length * length;
	return error;
}