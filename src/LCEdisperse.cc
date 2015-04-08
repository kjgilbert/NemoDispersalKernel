/** $Id: LCEdisperse.cc,v 1.19 2015-02-05 16:12:38 fred Exp $
 *
 *  @file LCEdisperse.cc
 *  Nemo2
 *
 *   Copyright (C) 2006-2011 Frederic Guillaume
 *   frederic.guillaume@env.ethz.ch
 *
 *   This file is part of Nemo
 *
 *   Nemo is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   Nemo is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *  
 *  Created on @date 08.07.2004
 *  @author fred
 */

#include <stdlib.h>
#include <string>
#include <map>
#include <iterator>
#include <cmath>
#include "LCEdisperse.h"
#include "metapop.h"
#include "individual.h"
#include "Uniform.h"

using namespace std;
// ------------------------------------------------------------------------------

//                             LCE_Disperse_base

// ----------------------------------------------------------------------------------------
LCE_Disperse_base::LCE_Disperse_base () 
: LifeCycleEvent("", ""), _disp_model(-1), _disp_propagule_prob(-1.0), _PropaguleTargets()
, _fem_rate (-1), _mal_rate(-1), _isForward(1)
{
  _DispMatrix[0] = NULL;
  _DispMatrix[1] = NULL;
//  cout << "calling LCE_Disperse_base()\n";
}
// ----------------------------------------------------------------------------------------
LCE_Disperse_base::~LCE_Disperse_base ()
{
  if(NULL != _DispMatrix[0])
    delete _DispMatrix[0];
  
  if(NULL != _DispMatrix[1])
    delete _DispMatrix[1];
}
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::addParameters (string prefix, ParamUpdaterBase* updater)
{
  add_parameter(prefix + "_model",INT,false,true,1,4,updater);
  add_parameter(prefix + "_border_model",INT,false,true,1,3,updater);
  add_parameter(prefix + "_lattice_range",INT,false,true,1,2,updater);
  add_parameter(prefix + "_propagule_prob",DBL,false,true,0,1,updater);
  add_parameter(prefix + "_matrix",MAT,false,false,0,0,updater);
  add_parameter(prefix + "_matrix_fem",STR,false,false,0,0,updater);
  add_parameter(prefix + "_matrix_mal",STR,false,false,0,0,updater);
  add_parameter(prefix + "_rate",DBL,false,true,0,1,updater);
  add_parameter(prefix + "_rate_fem",DBL,false,true,0,1,updater);
  add_parameter(prefix + "_rate_mal",DBL,false,true,0,1,updater);
  // adding my own after here KJG:
  add_parameter(prefix + "_matrix_xy",MAT,false,false,0,0,updater); // this will be the x and y coordinates to take the dispersal function to the aimed patch once it finds the index of where the migrant will go
  add_parameter(prefix + "_kernel_sorted",MAT,false,false,0,0,updater); // this will be the 1-d array holding the sorted probabilities of dispersing to patches 1 through n, and corresponding to the x,y coordinates in the above matrix. not sure yet if I need to make one for x and one for y or if this will suffice

}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setBaseParameters(string prefix)
{
  _prefix = prefix;
  
  _npatch = _popPtr->getPatchNbr();
  
  _disp_model = (int)_paramSet->getValue(prefix + "_model"); // set dispersal model to 1-4 if a model has been fed in, not a disp. matrix
  
  _disp_propagule_prob = _paramSet->getValue(prefix + "_propagule_prob");
  
  for (unsigned int sex = 0; sex < 2; sex++) {
    if(_DispMatrix[sex]) { // false is zero, true is 1
      delete _DispMatrix[sex]; // so if _DispMatrix[0] equals 1, then delete _DispMatrix[0]?
      _DispMatrix[sex] = NULL;
    }
  }
  
  if(_paramSet->isSet(prefix + "_matrix")) { // set normal dispersal matrix here if meet criteria
    
    _DispMatrix[0] = new TMatrix();
    
    _paramSet->getMatrix(prefix + "_matrix",_DispMatrix[0]);
    
    //same dispersal matrix for males and females
    _DispMatrix[1] = new TMatrix(*_DispMatrix[0]);
    
//  } else { // else do sex-specific dispersal matrices

  } else if(_paramSet->isSet(prefix + "_matrix_xy")) { // this is true if the input includes the xy coordinate matrix
	
	// want to be sure then that the sorted kernel has also been given in the input file
	if(!_paramSet->isSet(prefix + "_kernel_sorted")) { // if not set, return error, i.e. if it is set, the ! should make it return false and not throw the error
	  error("Dispersal rate parameters not set!\n");
      return false;
    }
    
    // HERE MAKE/RUN A FUNCTION TO DO MY EDITED DISPERSAL METHOD
cout << "tesing";

	} else{ //continue onto next line of original code to do sex-specific dispersal matrices
    
    if(_paramSet->isSet(prefix + "_matrix_fem")) {
      
      _DispMatrix[FEM] = new TMatrix();
      
      _paramSet->getMatrix(prefix + "_matrix_fem",_DispMatrix[FEM]);
      
    }
    
    if(_paramSet->isSet(prefix + "_matrix_mal")) {
      
      _DispMatrix[MAL] = new TMatrix();
      
      _paramSet->getMatrix(prefix + "_matrix_mal",_DispMatrix[MAL]);
    
    }
  }  

// add this stuff to see if that causes problem of looking for other parameters that I had before I got the segmentation fault using the or statement below
  if( _paramSet->isSet(prefix + "_matrix_xy") ) 
  {
    if(  ( _paramSet->isSet(prefix + "_rate") ||
          (_paramSet->isSet(prefix + "_rate_fem") &&  _paramSet->isSet(prefix + "_rate_mal")) )
       || _paramSet->isSet(prefix + "_model") )
      warning("parameter \"dispersal_matrix\" takes precedence over parameters \"dispersal_rate\" and \"dispersal_model\"\n");
    
    _disp_model = 0;

  }


// below code sends to setReducedDispMatrix, don't want that for my edits, return to original state:
  if( _paramSet->isSet(prefix + "_matrix") || 
     ( _paramSet->isSet(prefix + "_matrix_fem") && _paramSet->isSet(prefix + "_matrix_mal") )  ) 

  {
    if(  ( _paramSet->isSet(prefix + "_rate") ||
          (_paramSet->isSet(prefix + "_rate_fem") &&  _paramSet->isSet(prefix + "_rate_mal")) )
       || _paramSet->isSet(prefix + "_model") )
      warning("parameter \"dispersal_matrix\" takes precedence over parameters \"dispersal_rate\" and \"dispersal_model\"\n");
    
    _disp_model = 0; // if not dispersal models 1-4, set to zero, because a matrix has been fed in from the input

    if(_DispMatrix[FEM]) {

      if(_DispMatrix[FEM]->length() != _npatch*_npatch) {
        error("the size of the female dispersal matrix is not equal to patch_number X patch_number (%i[%i,%i] != %i)!\n",
              _DispMatrix[FEM]->length(),_DispMatrix[FEM]->getNbRows(),_DispMatrix[FEM]->getNbCols(),_npatch*_npatch);
        return false;
      }
      if(_isForward) {
        if(!checkForwardDispersalMatrix(_DispMatrix[FEM])) return false;
      } else {
        if(!checkBackwardDispersalMatrix(_DispMatrix[FEM])) return false; // check that these both sum to 1, either across rows for forward or across columns for backwards
      }
      
    }
    
    if(_DispMatrix[MAL]) {
      
      if(_DispMatrix[MAL]->length() != _npatch*_npatch) {
        error("the size of the male dispersal matrix is not equal to patch_number X patch_number (%i[%i,%i] != %i)!\n",
                _DispMatrix[MAL]->length(),_DispMatrix[MAL]->getNbRows(),_DispMatrix[MAL]->getNbCols(),_npatch*_npatch);
        return false;
      }
      if(_isForward) {
        if(!checkForwardDispersalMatrix(_DispMatrix[MAL])) return false;
      } else {
        if(!checkBackwardDispersalMatrix(_DispMatrix[MAL])) return false;
      }
    }
    
    setReducedDispMatrix(); // calls on setReducedDispMatrix once has read in all matrices so it can order patches for optimal searching rather than searching all despite order of probabilities
    
  } else {
   if(!_paramSet->isSet(prefix + "_matrix_xy")){ // start if statement I added to only do this if not using my added dispersal kernel method
       
    if(!_paramSet->isSet(prefix + "_model")) {
      error("Dispersal model not set!\n");
      return false;
    }
    
    if(_paramSet->isSet(prefix + "_rate")) 
      
    {
      _fem_rate = _mal_rate = _paramSet->getValue(prefix + "_rate");
      
      if(!setDispMatrix()) return false;
    }
    
    else if(  _paramSet->isSet(prefix + "_rate_fem") &&  _paramSet->isSet(prefix + "_rate_mal")  ) 
      
    {
      _fem_rate = _paramSet->getValue(prefix + "_rate_fem");
      
      _mal_rate = _paramSet->getValue(prefix + "_rate_mal");
      
      if(!setDispMatrix()) return false;
    }
    
    else {
      error("Dispersal rate parameters not set!\n");
      return false;
    }

    if(_isForward) {
      if(!checkForwardDispersalMatrix(_DispMatrix[FEM])) return false;
      if(!checkForwardDispersalMatrix(_DispMatrix[MAL])) return false;
    } else {
      if(!checkBackwardDispersalMatrix(_DispMatrix[FEM])) return false;
      if(!checkBackwardDispersalMatrix(_DispMatrix[MAL])) return false;
    }
  } // end if statement I added to only do this if not using my added dispersal kernel method
  }  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::updateDispMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::updateDispMatrix()
{ //called by the 'execute' function when a change in patch number is detected 
  if ( getDispersalModel() == 0) {
    error("cannot update the dispersal matrix provided in input when number of populations changes.\n");
    return false;
  }
  
  return setDispMatrix();
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::swapPostDisp
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::swapPostDisp ( )
{
  Patch *patch;
  
  for(unsigned int i = 0; i < _npatch; i++) {
    patch = _popPtr->getPatch(i);
    patch->swap(FEM, PDISPx, OFFSx);
    patch->swap(MAL, PDISPx, OFFSx);
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::reset_counters
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::reset_counters()
{
  Patch *patch;
  for(unsigned int i = 0; i < _npatch; i++) {
    
    patch = _popPtr->getPatch(i);
    
    patch->reset_counters();
    
    patch->flush(PDISPx, _popPtr);
  }  
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::allocateDispMatrix
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::allocateDispMatrix (sex_t sex, unsigned int dim)
{  
  if(_DispMatrix[sex] != NULL)
    _DispMatrix[sex]->reset(dim,dim);
  else
    _DispMatrix[sex] = new TMatrix(dim,dim);
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::checkDispMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::checkForwardDispersalMatrix (TMatrix* mat)
{
//  unsigned int dim = _popPtr->getPatchNbr();
  double cntr;

  //add one patch for absorbing border model
//  if( get_parameter_value( _prefix + "_border_model") == 3 ) dim++;
//  
//  if(mat->length() != dim*dim) {
//    error("the size of the dispersal matrix is not equal to patch_number X patch_number (%i[%i,%i] != %i)!\n",
//          mat->length(),mat->getNbRows(),mat->getNbCols(),dim*dim);
//    return false;
//  }
  
  for(unsigned int i = 0; i < mat->getNbRows(); ++i) {
    cntr = 0;
    for(unsigned int j = 0; j < mat->getNbCols(); ++j)
      cntr += mat->get(i,j);
    if(cntr < 0.999999 || cntr > 1.000001) {
      error("the elements of row %i of the dispersal matrix do not sum to 1!\n",i+1);
      error("sum of row %i is: %f\n",i+1, cntr);
      return false;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::checkDispMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::checkBackwardDispersalMatrix (TMatrix* mat)
{
//  unsigned int dim = _popPtr->getPatchNbr();
  double cntr;
  
  //add one patch for absorbing border model
//  if( get_parameter_value( _prefix + "_border_model") == 3 ) dim++;
//  
//  if(mat->length() != dim*dim) {
//    error("The size of the dispersal matrix is not equal to patch_number X patch_number (%i[%i,%i] != %i)!\n",
//          mat->length(),mat->getNbRows(),mat->getNbCols(),dim*dim);
//    return false;
//  }
  
  for(unsigned int i = 0; i < mat->getNbCols(); ++i) {
    cntr = 0;
    for(unsigned int j = 0; j < mat->getNbRows(); ++j)
      cntr += mat->get(j,i);
    if(cntr < 0.999999 || cntr > 1.000001) {
      error("The elements of column %i of the dispersal matrix do not sum to 1!\n",i+1);
      return false;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setPropaguleTargets
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::setPropaguleTargets ( )
{
  unsigned int nb_patch = _popPtr->getPatchNbr();
  unsigned int tmp_array[nb_patch];
  unsigned int table_index, target_patch;
  unsigned int size, last;
  
  //shuffling algorithm:
  do {
    for(unsigned int i = 0; i < nb_patch; ++i)
      tmp_array[i] = i;
    
    size = nb_patch;
    
    for(unsigned int orig_patch = 0; orig_patch < nb_patch-1; ++orig_patch) {
      
      do{
        
        table_index = RAND::Uniform( size );
        
        target_patch = tmp_array[ table_index ];
        
      }while(target_patch == orig_patch);
      
      size--;
      
      last = tmp_array[size];
      
      tmp_array[table_index] = last;
      
      tmp_array[size] = target_patch;
    }
    //do this until the last element left is not the last patch:
  }while (tmp_array[0] == nb_patch-1);
  
  _PropaguleTargets.assign(nb_patch,0);
  
  unsigned int reverse_i = nb_patch;
  
  //we read the shuffled array in reverse order:
  for(unsigned int i=0; i < _PropaguleTargets.size(); i++) {
    _PropaguleTargets[i] = tmp_array[--reverse_i];
    
#ifdef _DEBUG_
    cout<<" -- Patch "<<i<<" : assigned Patch "<<_PropaguleTargets[i]<<endl;
#endif
    
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setDispMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setDispMatrix ()
{  
  
  switch ( getDispersalModel() ) {
    case 1:
      if( !setIsland_MigrantPool_Matrix() ) return false;
      break;
    case 2:
      if( !setIsland_PropagulePool_Matrix() ) return false;
      break;
    case 3:
      if( !setSteppingStone1DMatrix() ) return false;
      break;
    case 4:
      if( !setLatticeMatrix() ) return false;
      break;
    default: {
      error("Dispersal model '%i' not yet implemented\n",getDispersalModel());
      return false;
    }
  }
  
  return setReducedDispMatrix(); // also call on reduced disp matrix here if none of the other dispersal models has been set
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setIsland_MigrantPool_Matrix()  (set the Island dispersal matrix)
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setIsland_MigrantPool_Matrix()
{
#ifdef _DEBUG_
  cout<<"setIsland_MigrantPool_Matrix(_npatch="<<_npatch<<", _mal_rate="
      <<_mal_rate<<", _fem_rate="<<_fem_rate<<")"<<endl;
#endif
  allocateDispMatrix(MAL, _npatch);
  allocateDispMatrix(FEM, _npatch);
  
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  double pmal = 1 - _mal_rate;
  double pfem = 1 - _fem_rate;
  double mmal = _mal_rate/(_npatch-1);
  double mfem = _fem_rate/(_npatch-1);
  
  for (unsigned int i=0; i<_npatch; ++i){
    for (unsigned int j=0; j<_npatch; ++j){
      mmat->set(i,j, mmal);
      fmat->set(i,j, mfem);
    }
  }
  
  for (unsigned int i=0; i<_npatch; ++i){
    mmat->set(i,i, pmal);
    fmat->set(i,i, pfem);
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setIsland_PropagulePool_Matrix()  (set the Island dispersal matrix)
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setIsland_PropagulePool_Matrix()
{
  allocateDispMatrix(MAL, _npatch);
  allocateDispMatrix(FEM, _npatch);
  
  if( !_paramSet->isSet(_prefix + "_propagule_prob") ) {
    error("Missing parameter \"dispersal_propagule_prob\" with dispersal model 2!\n");
    return false;
  }
  
  setPropaguleTargets();
  
  double propagulePHI = getPropaguleProb();
  double c1 = (1 - _fem_rate), c2 = (_fem_rate*propagulePHI),
	c3 = (_fem_rate*(1.0 - propagulePHI)/(_npatch-2));
  
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  for (unsigned int i=0; i < _npatch; ++i){
		
    fmat->set(i, i, c1);
    
    for (unsigned int j=i+1; j < _npatch; ++j){
      fmat->set(i, j, c3);
      fmat->set(j, i, c3);
    }
    fmat->set(i, getPropaguleTarget(i), c2);
  }
  
  c1 = (1 - _mal_rate);
  c2 = (_mal_rate*propagulePHI);
  c3 = (_mal_rate*(1.0 - propagulePHI)/(_npatch-2));
  
  for (unsigned int i=0; i < _npatch; ++i){
    
    mmat->set(i, i, c1);
    
    for (unsigned int j=i+1; j< _npatch; ++j) {
      mmat->set(i, j, c3);
      mmat->set(j, i, c3);
    }
    mmat->set(i, getPropaguleTarget(i), c2);
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setSteppingStone1DMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setSteppingStone1DMatrix()
{
#ifdef _DEBUG_
  message("setSteppingStone1DMatrix()\n");
#endif
  int border_model = (unsigned int)get_parameter_value(_prefix + "_border_model");
  
  //check for the border model, the extra patch is the absorbing patch
  if(border_model == 3) {
    allocateDispMatrix(MAL, _npatch+1);
    allocateDispMatrix(FEM, _npatch+1);
  } else {    
    allocateDispMatrix(MAL, _npatch);
    allocateDispMatrix(FEM, _npatch);
  }
  
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  //philopatry:
  double pmal = 1 - _mal_rate, pfem = 1 - _fem_rate;
  //migration:
  double mmal = _mal_rate/2, mfem = _fem_rate/2;
  
  fmat->assign(0);
  mmat->assign(0);
  
  //diagonal:
  for (unsigned int i = 0; i < _npatch; ++i){
    fmat->set(i, i, pfem);
    mmat->set(i, i, pmal);
  }
  //around the diagonal
  for (unsigned int i = 0; i < _npatch-1; ++i){
    fmat->set(i, i+1, mfem);
    fmat->set(i+1, i, mfem);
    mmat->set(i, i+1, mmal);
    mmat->set(i+1, i, mmal);
  }
  
  if(border_model == 3) {
    //absorbing boders
    //emigrants who leave the population are put in the sink patch
    fmat->set(0, _npatch, mfem);
    mmat->set(0, _npatch, mmal);
    fmat->set(_npatch -1, _npatch, mfem);
    mmat->set(_npatch -1, _npatch, mmal);
    fmat->set(_npatch, _npatch, 1);
    mmat->set(_npatch, _npatch, 1);
    
    if(!_isForward) { 
      fmat->transpose();  mmat->transpose(); //this creates artificial immigration from sink
      //need to reset border patches as immigration from sink not allowed:
      fmat->set(0, 0, pfem + mfem); //the proportion not comming from the sink must be added 
      fmat->set(_npatch, 0, 0);  //no immigration from sink
      fmat->set(_npatch-1, _npatch-1, pfem + mfem);
      fmat->set(_npatch, _npatch-1, 0);
      
      mmat->set(0, 0, pmal + mmal);
      mmat->set(_npatch, 0, 0);
      mmat->set(_npatch-1, _npatch-1, pmal + mmal);
      mmat->set(_npatch, _npatch-1, 0);
    }
    
  } else if (border_model == 2) {
    //reflective borders, 
    //emigrants that cannot leave stay in place
    fmat->set(0, 0, pfem+mfem);
    mmat->set(0, 0, pmal+mmal);
    fmat->set(_npatch -1, _npatch -1, pfem+mfem);
    mmat->set(_npatch -1, _npatch -1, pmal+mmal);
    
    //no need to transpose for backward migration as the matrix is symmetrical
    
  } else { //is a torus by default
    //the 2 last elements, this is a ring population!
    fmat->set(0, _npatch -1, mfem);
    mmat->set(0, _npatch -1, mmal);
    fmat->set(_npatch -1, 0, mfem);
    mmat->set(_npatch -1, 0, mmal);
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeMatrix
// ----------------------------------------------------------------------------------------
/**Sets the dispersal matrices for the Lattice dispersal model.*/
bool LCE_Disperse_base::setLatticeMatrix()
{
#ifdef _DEBUG_
  message("setLatticeMatrix()\n");
#endif
  if(!_paramSet->isSet(_prefix + "_border_model")) {
    error("Missing parameter \"dispersal_border_model\" with dispersal model 4!\n");
    return false;
  }
  
  if(!_paramSet->isSet(_prefix + "_lattice_range") ) {
    error("Missing parameter \"dispersal_lattice_range\" with dispersal model 4!\n");
    return false;
  }
  
  unsigned int side = (unsigned int) sqrt((double)_npatch);
  
  if( side*side != _npatch ) {
    error("The number of patches is not a square number in the lattice dispersal model\n");
    return false;
  }
  
  /**Each matrix has 'patch number' x 'patch number' cells unless the lattice model is the absorbing
   boundaries model where we add the sink patch.*/
  if((unsigned int)get_parameter_value(_prefix + "_border_model") == 3) {
    allocateDispMatrix(MAL, _npatch+1);
    allocateDispMatrix(FEM, _npatch+1);
  } else {
    allocateDispMatrix(MAL, _npatch);
    allocateDispMatrix(FEM, _npatch);
  }
  
  /**The "dispersal_lattice_range" parameter defines the number of neighbouring patches to disperse into. 
   Option 1 sets this number to 4 (left and right, up and down patches) whereas option 2 allows to disperse to the 8 neighbouring
   patches, including the patches on the diagonals.*/
  unsigned int range = (unsigned int)get_parameter_value(_prefix + "_lattice_range");
  //philopatry:
  double pmal = 1 - _mal_rate, pfem = 1 - _fem_rate;
  //migration:
  double mmal = _mal_rate/(range == 1 ? 4 : 8), mfem = _fem_rate/(range == 1 ? 4 : 8);
  
  setBasicLatticeMatrix(side, pmal, pfem, mmal, mfem);
  
  
  return true;
  
  //  TMatrix* fmat = _DispMatrix[FEM];
  //  for (unsigned int i=0; i<_npatch; ++i){
  //    for (unsigned int j=0; j<_npatch; ++j)
  //      cout<<fmat->get(i,j)<<"\t";
  //    cout<<endl;
  //  }
} 
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setBasicLatticeMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setBasicLatticeMatrix(int side, double phi_mal, double phi_fem,
                                              double disp_mal, double disp_fem) 
{  
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  //init:
  fmat->assign(0.0);
  mmat->assign(0.0);
  
  TMatrix grid(side, side);
  
  int p = 0;
  for (int i = 0; i < side; ++i) {
    for (int j = 0; j < side ; ++j) {
      grid.set(i, j, p++);
    }
  }
  assert( p == side*side );
  
  //diagonal:
  for (unsigned int i = 0; i < _npatch; ++i){
    fmat->set(i, i, phi_fem);
    mmat->set(i, i, phi_mal);
  }
  
  int connect;
  
  for (int x = 0; x < side; ++x) {
    for (int y = 0; y < side; ++y) {
      
      p = grid.get(x, y);
      
      connect = p + 1; //patch to the right
    
      if (connect < (x + 1)*side) { //stay on the same row
        fmat->set(p, connect, disp_fem);
        mmat->set(p, connect, disp_mal);
      } 
      
      connect = p - 1; //patch to the left
      
      if (connect >= x*side ) {
        fmat->set(p, connect, disp_fem);
        mmat->set(p, connect, disp_mal);
      }

      connect = p + side; //patch one row up
      
      if (connect < (int)_npatch ) {
        fmat->set(p, connect, disp_fem);
        mmat->set(p, connect, disp_mal);
      }

      connect = p - side; //patch one row down
      
      if (connect >= 0 ) {
        fmat->set(p, connect, disp_fem);
        mmat->set(p, connect, disp_mal);
      }
      
      //diagonal steps:
      if((unsigned int)get_parameter_value(_prefix + "_lattice_range") == 2) {
        
        connect = p + side; //patches one row up
        
        if (connect < (int)_npatch ) { //we are not on the last row
          
          if (connect + 1 < (x + 2)*side) {
            fmat->set(p, connect + 1, disp_fem);
            mmat->set(p, connect + 1, disp_mal);
          }
          
          if (connect - 1 >= (x + 1)*side) {
            fmat->set(p, connect - 1, disp_fem);
            mmat->set(p, connect - 1, disp_mal);
          }
        }
        
        connect = p - side; //patches one row down
        
        if (connect >= 0) { //we are not on the first row
          
          if (connect + 1 < (x + 1)*side) {
            fmat->set(p, connect + 1, disp_fem);
            mmat->set(p, connect + 1, disp_mal);
          }
          
          if (connect - 1 >= (x - 1)*side) {
            fmat->set(p, connect - 1, disp_fem);
            mmat->set(p, connect - 1, disp_mal);
          }            
        }
      } //lattice range == 2
      
    } //y
  } //x
  
  switch((unsigned int)get_parameter_value(_prefix + "_border_model")) {
    case 1:
      return setLatticeTorrusMatrix((int)side, disp_mal, disp_fem, &grid);
    case 2:
      return setLatticeReflectingMatrix((int)side, &grid);
    case 3:
      return setLatticeAbsorbingMatrix();
    default:
      error("parameter \"%s_border_model\" accepts only three values: [1,2,3]\n", _prefix.c_str());
      break;
  }  
  
  return false;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeTorrusMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setLatticeTorrusMatrix(int side, double disp_mal, double disp_fem, TMatrix* grid) 
{
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];

  int x, y;
  
  x = 0; //first row, connect with upper row
  for (y = 0; y < side; ++y) {
    fmat->set( grid->get(x, y), grid->get( side-1, y ), disp_fem);
    mmat->set( grid->get(x, y), grid->get( side-1, y ), disp_mal);
  }
  
  x = side - 1; //last row, connect with bottom row
  for (y = 0; y < side; ++y) {
    fmat->set( grid->get(x, y), grid->get( 0, y ), disp_fem);
    mmat->set( grid->get(x, y), grid->get( 0, y ), disp_mal);
  }
  
  y = 0; //left column, connect with right-most column
  for (x = 0; x < side; ++x) {
    fmat->set( grid->get(x, y), grid->get( x, side - 1 ), disp_fem);
    mmat->set( grid->get(x, y), grid->get( x, side - 1 ), disp_mal);
  }

  y = side - 1; //right-most column, connect with left-most column
  for (x = 0; x < side; ++x) {
    fmat->set( grid->get(x, y), grid->get( x, 0 ), disp_fem);
    mmat->set( grid->get(x, y), grid->get( x, 0 ), disp_mal);
  }
  
  //connect along diagonals
  
  if((int)get_parameter_value(_prefix + "_lattice_range") == 2) {
    
    
    //CORNERS:  
    
    //first patch on first row, wrap around to top right border (last patch)
    fmat->set( grid->get(0, 0), grid->get( side - 1, side - 1 ), disp_fem);
    mmat->set( grid->get(0, 0), grid->get( side - 1, side - 1 ), disp_mal);
    //last patch to first patch
    fmat->set( grid->get( side - 1, side - 1 ), grid->get(0, 0), disp_fem);
    mmat->set( grid->get( side - 1, side - 1 ), grid->get(0, 0), disp_mal);
    
    //last patch on first row, wrap around to top left border
    fmat->set( grid->get(0, side - 1), grid->get( side - 1, 0 ), disp_fem);
    mmat->set( grid->get(0, side - 1), grid->get( side - 1, 0 ), disp_mal);
    //third corner to second corner
    fmat->set( grid->get( side - 1, 0 ), grid->get(0, side - 1), disp_fem);
    mmat->set( grid->get( side - 1, 0 ), grid->get(0, side - 1), disp_mal);
    
    
    //BORDERS:
    
    x = 0;//we are on the first row
    
    for (y = 0; y < side - 1; ++y) {
      //diagonal step to the right, move to the other side of the grid
      fmat->set( grid->get(x, y), grid->get( side-1, y + 1 ), disp_fem);
      mmat->set( grid->get(x, y), grid->get( side-1, y + 1 ), disp_mal);
    }
    for (y = 1; y < side; ++y) {
      //diagonal step to the left, move to the other side of the grid
      fmat->set( grid->get(x, y), grid->get( side-1, y - 1 ), disp_fem);
      mmat->set( grid->get(x, y), grid->get( side-1, y - 1 ), disp_mal);
    }
    
    x = side - 1; //last row
    
    for (y = 0; y < side - 1; ++y) {
      //diagonal step to the right, move to the first row
      fmat->set( grid->get(x, y), grid->get( 0, y + 1 ), disp_fem);
      mmat->set( grid->get(x, y), grid->get( 0, y + 1 ), disp_mal);
    }
    for (y = 1; y < side; ++y) {
      //diagonal step to the left, move to the first row
      fmat->set( grid->get(x, y), grid->get( 0, y - 1 ), disp_fem);
      mmat->set( grid->get(x, y), grid->get( 0, y - 1 ), disp_mal);
    }
        
    //in-between rows, connect first and last columns
    for (x = 0; x < side; ++x) {
      
      if(x > 0) {
        //diagonal step down, to previous row
        y = 0; //first colum
        fmat->set( grid->get(x, y), grid->get( x - 1, side - 1 ), disp_fem);
        mmat->set( grid->get(x, y), grid->get( x - 1, side - 1 ), disp_mal);
        y = side - 1; //last colum
        fmat->set( grid->get(x, y), grid->get( x - 1, 0 ), disp_fem);
        mmat->set( grid->get(x, y), grid->get( x - 1, 0 ), disp_mal);
      }
      if(x < side - 1){
        //diagonal step up, to next row
        y = 0;
        fmat->set( grid->get(x, y), grid->get( x + 1, side - 1 ), disp_fem);
        mmat->set( grid->get(x, y), grid->get( x + 1, side - 1 ), disp_mal);
        y = side - 1;
        fmat->set( grid->get(x, y), grid->get( x + 1, 0 ), disp_fem);
        mmat->set( grid->get(x, y), grid->get( x + 1, 0 ), disp_mal);
      }
    }
  }
//the matrix is symmetrical no need to transpose for backward migration
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeReflectingMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setLatticeReflectingMatrix(int side, TMatrix* grid)
{
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
//  double dm, df;
//  unsigned int range = (unsigned int)get_parameter_value(_prefix + "_lattice_range");
  
  vector<unsigned int> border_cells;
  
  for (unsigned int i = 0; i < _npatch; ++i) {
    
    for (int j = 0; j < side; ++j) 
      border_cells.push_back( j ); //first row
    
    for (int j = 1; j < side - 1; ++j) { //intermediate rows, only the two border cells
      border_cells.push_back( j*side );
      border_cells.push_back( (j + 1)*side - 1 );
    }
    
    for (unsigned int j = side*(side-1); j < _npatch; ++j) 
      border_cells.push_back( j ); //last row
    
  }
  
  double sum;
  // individuals who can't disperse past the border stay in place:
  for (unsigned int i = 0; i < border_cells.size(); ++i) {
    //sum of dispersal rates:
    sum = fmat->rowSum( border_cells[i] ) - fmat->get(border_cells[i], border_cells[i]);
    //difference gives 1 - m:
    fmat->set(border_cells[i], border_cells[i], 1 - sum);
    
    sum = mmat->rowSum( border_cells[i] ) - mmat->get(border_cells[i], border_cells[i]);
    mmat->set(border_cells[i], border_cells[i], 1 - sum);
  }
  
  if(!_isForward) {  fmat->transpose();  mmat->transpose(); }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeAbsorbingMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setLatticeAbsorbingMatrix()
{
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  fmat->set(_npatch, _npatch, 1.0);
  mmat->set(_npatch, _npatch, 1.0);
  
  double sum;
  
  if(_isForward) {
  //set the absorbing patch probs to 1 - sum(row)
    for(unsigned int i = 0; i < _npatch; ++i) {
      sum = 0;
      for (unsigned int j = 0; j < _npatch; ++j) {
        sum += fmat->get(i, j);
      }
      fmat->set(i, _npatch, 1.0 - sum);
    }
    
    for(unsigned int i = 0; i < _npatch; ++i) {
      sum = 0;
      for (unsigned int j = 0; j < _npatch; ++j) {
        sum += mmat->get(i, j);
      }
      mmat->set(i, _npatch, 1.0 - sum);
    }
  //backward migration:  
  }else {
    //the missing immigrant rate from non-existing patches must be added to the "philopatric" rate
    for(unsigned int i = 0; i < _npatch; ++i) {
      sum = 0;
      for (unsigned int j = 0; j < _npatch; ++j) {
        sum += fmat->get(j, i);
      }
      fmat->plus(i, i, 1.0 - sum);
    }
    
    for(unsigned int i = 0; i < _npatch; ++i) {
      sum = 0;
      for (unsigned int j = 0; j < _npatch; ++j) {
        sum += mmat->get(j, i);
      }
      mmat->plus(i, i, 1.0 - sum);
    }
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setReducedDispMatrix
// ----------------------------------------------------------------------------------------
/** The reduced dispersal matrix contains the indices of the patches to which each patch is
 connected. The connected patches are further ordered in descending order of the migration rates.
 This offers a double speed-up compared to the classical method.
*/
bool LCE_Disperse_base::setReducedDispMatrix() /// CAN USE THIS FOR DISP KERNEL MODIFICATIONS I'LL BE MAKING
{
  unsigned int border_model = (unsigned int)get_parameter_value(_prefix + "_border_model"); // check if there are reflecting or absorbing boundaries set from the input
  unsigned int num_patch = (border_model == 3 ? _npatch + 1 : _npatch); // SOMEHOW? check on this. this is getting the number of patches from that code?

  //multimap automatically orders the key values in ascending order
  multimap<double, unsigned int> ordered_rates_fem, ordered_rates_mal;
  typedef multimap<double, unsigned int>::const_iterator CI;

#ifdef _DEBUG_
  message("== Dispersal matrices ==\n");
    _DispMatrix[FEM]->show_up();
    _DispMatrix[MAL]->show_up();
#endif
  for (unsigned int sex = 0; sex < 2; sex++)
    if(_reducedDispMat[sex].size() != 0) _reducedDispMat[sex].clear();
  
  
  for (unsigned int i = 0; i < num_patch; ++i) { // go through all the patches
    
    _reducedDispMat[0].push_back(vector<unsigned int>());
    _reducedDispMat[1].push_back(vector<unsigned int>());
    
    ordered_rates_fem.clear();
    ordered_rates_mal.clear();
    
    if(_isForward) { // if doing forward migration
      
      for (unsigned int j = 0; j < num_patch; ++j)
        if(_DispMatrix[0]->get(i, j) != 0) ordered_rates_mal.insert(make_pair(_DispMatrix[0]->get(i, j), j));
      
      
      for (unsigned int j = 0; j < num_patch; ++j)      // go through all the patches
        if(_DispMatrix[1]->get(i, j) != 0) ordered_rates_fem.insert(make_pair(_DispMatrix[1]->get(i, j),j));

      
    } else { // otherwise we are doing backwards migration
      //backward migration matrices are read column-wise
      for (unsigned int j = 0; j < num_patch; ++j)
        if(_DispMatrix[0]->get(j, i) != 0) ordered_rates_mal.insert(make_pair(_DispMatrix[0]->get(j, i),j));
      
      for (unsigned int j = 0; j < num_patch; ++j)
        if(_DispMatrix[1]->get(j, i) != 0) ordered_rates_fem.insert(make_pair(_DispMatrix[1]->get(j, i),j));
    }
    
//    cout << "=== female ordered rates ===\n";
//    cout << "ordered_rates_fem.size = " << ordered_rates_fem.size()<<"\n";
//    for (CI p = ordered_rates_fem.begin(); p != ordered_rates_fem.end(); p++) {
//      cout << p->first <<", "<< p->second<<"\n";
//    }
    
    if(ordered_rates_fem.size() == 1) {
        _reducedDispMat[1][i].push_back(ordered_rates_fem.begin()->second);
      //store the patch indices in reverse order of the migration rates:
    } else {
        CI p;
      for (p = --ordered_rates_fem.end(); p != --ordered_rates_fem.begin(); --p) {
          _reducedDispMat[1][i].push_back(p->second);
      }
  
      //bug fix for the case of 2 patches only (for some reason the second patch recieves only one value...???)
      if(p == ordered_rates_fem.begin() && (_reducedDispMat[1][i].size() < ordered_rates_fem.size()) )
        _reducedDispMat[1][i].push_back(p->second);
    }

    
    if(ordered_rates_mal.size() == 1) {
      _reducedDispMat[0][i].push_back(ordered_rates_mal.begin()->second);
    } else {
      CI p;
      for (p = --ordered_rates_mal.end(); p != --ordered_rates_mal.begin(); --p) {
        _reducedDispMat[0][i].push_back(p->second);
      }
      
      if(p == ordered_rates_mal.begin() && (_reducedDispMat[0][i].size() < ordered_rates_mal.size()) )
        _reducedDispMat[0][i].push_back(p->second);
    }
  }
#ifdef _DEBUG_
    cout << "=== female reduced dispersal matrix ===\n";
  for (unsigned int i = 0; i < num_patch; ++i) {
    cout << "  [";
    for (unsigned int k = 0; k < _reducedDispMat[FEM][i].size(); k++) {
      cout << _reducedDispMat[FEM][i][k] <<" ";
    }
    cout<<"]\n";
  }
  
    cout << "=== male reduced dispersal matrix ===\n";
  for (unsigned int i = 0; i < num_patch; ++i) {
    cout << "  [";
    for (unsigned int k = 0; k < _reducedDispMat[MAL][i].size(); k++) {
      cout << _reducedDispMat[MAL][i][k] <<" ";
    }
    cout<<"]\n";
    
  }
#endif
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::Migrate Forward
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_base::getMigrationPatchForward (sex_t SEX, unsigned int LocalPatch)
{
  double sum = 0, random = RAND::Uniform(); // draw a random number
  unsigned int AimedPatch = 0;
  
  if(random > 0.999999) random = 0.999999;//this to avoid overflows when random == 1
  
  sum = _DispMatrix[SEX]->get(LocalPatch, _reducedDispMat[SEX][LocalPatch][AimedPatch]); // FIGURE OUT THIS LINE

  while (random > sum) { // find the aimed patch whose probability matches that of the random number drawn, i.e. the patch that will be migrated into
    AimedPatch++; // keep going through patches
    sum += _DispMatrix[SEX]->get(LocalPatch, _reducedDispMat[SEX][LocalPatch][AimedPatch]); // increase sum until hit the patch matching the drawn random number
  }

  return _reducedDispMat[SEX][LocalPatch][AimedPatch]; // FIGURE OUT THE DETAILS OF WHAT THIS RETURNS - this must be after finding the aimed patch, returning that patch's ID? what is the sex part?
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::Migrate Backward
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_base::getMigrationPatchBackward (sex_t SEX, unsigned int LocalPatch)
{
  double sum = 0, random = RAND::Uniform();
  unsigned int SourcePatch = 0;
  
  if(random > 0.999999) random = 0.999999;//this to avoid overflows when random == 1
    
  sum = _DispMatrix[SEX]->get(_reducedDispMat[SEX][LocalPatch][SourcePatch], LocalPatch);
  
  while (random > sum) {
    SourcePatch++;
    sum += _DispMatrix[SEX]->get(_reducedDispMat[SEX][LocalPatch][SourcePatch], LocalPatch);
  }
    
  return _reducedDispMat[SEX][LocalPatch][SourcePatch];
}
// ------------------------------------------------------------------------------

//                             LCE_Disperse_ConstDisp

// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp
// ----------------------------------------------------------------------------------------
LCE_Disperse_ConstDisp::LCE_Disperse_ConstDisp () 
: LifeCycleEvent ("disperse",""), doMigration(0), doPatchMigration(0)
{
  
  ParamUpdater< LCE_Disperse_ConstDisp > * updater = 
  new ParamUpdater< LCE_Disperse_ConstDisp > (&LCE_Disperse_ConstDisp::setParameters);
  
  LCE_Disperse_base::addParameters("dispersal", updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_ConstDisp::setParameters(string prefix)
{
  
  if(!LCE_Disperse_base::setBaseParameters(prefix)) return false;
  
  switch ( getDispersalModel() ) {
    case 0: //dispersal matrix given in input
    case 1:
    case 3:
    case 4:
      doMigration = &LCE_Disperse_ConstDisp::Migrate;
      break;
    case 2:
      doMigration = &LCE_Disperse_ConstDisp::Migrate_propagule;
      break;
    default: {
      error("\nDispersal model '%i' not yet implemented\n",getDispersalModel());
      return false;
    }
  }
  
  switch ((int)get_parameter_value(prefix + "_border_model")) {
    case -1:
    case 1:
    case 2:
      doPatchMigration = &LCE_Disperse_ConstDisp::MigratePatch;
      break;
    case 3:
      doPatchMigration = &LCE_Disperse_ConstDisp::MigratePatch_AbsorbingBorder;
      break;
    default:{
      error("\nDispersal border model '%i' not yet implemented\n",
            (int)get_parameter_value(prefix + "_border_model"));
      return false;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp::execute
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::execute ()
{
#ifdef _DEBUG_
  message("LCE_Disperse_ConstDisp::execute (Patch nb: %i offsprg nb: %i adlt nb: %i "
          ,_popPtr->getPatchNbr(), _popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif
  //check wether the number of patches changed during simulation:
  if(_npatch != _popPtr->getPatchNbr()) {
    _npatch = _popPtr->getPatchNbr();
    if(!updateDispMatrix()) fatal("bailing out\n");
  }
  
  reset_counters();
  
  (this->*doMigration)();
  
#ifdef _DEBUG_
  unsigned int c = 0;
  for(unsigned int i = 0; i < _npatch; i++)
    c += _popPtr->getPatch(i)->nbEmigrant;
  message("emigrants nb: %i)\n",c);
#endif
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp::Migrate
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::Migrate ()
{
  Patch *patch;
  
  for(unsigned int i = 0; i < _npatch; i++) {
    
    patch = _popPtr->getPatch(i);
    
    (this->*doPatchMigration)(FEM, i);
    
    (this->*doPatchMigration)(MAL, i);
    
    //set coloniser counter
    if(patch->get_isExtinct()) 
      patch->nbKolonisers = patch->size(PDISPx);
    else 
      patch->nbKolonisers = -1; //value used in stat_demo to check for colonization
    
  }//end for nb_patch
  
  //put back the indviduals into the offspring container
  swapPostDisp();
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::Migrate_propagule ()
{
  setIsland_PropagulePool_Matrix();
  Migrate();
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp::MigratePatch
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::MigratePatch (sex_t SEX, unsigned int LocalPatch)
{
  Patch* patch = _popPtr->getPatch(LocalPatch);
  unsigned int AimedPatch;
  unsigned int Limit = _npatch -1;
  
  while( patch->size(SEX, OFFSx) != 0 ) {
    
    do{
      AimedPatch = getMigrationPatchForward(SEX, LocalPatch);
    }while ( AimedPatch > Limit );

    _popPtr->move(SEX, OFFSx, LocalPatch, PDISPx, AimedPatch, 0);
    
    if(LocalPatch != AimedPatch) {
      patch->nbEmigrant++;
      _popPtr->getPatch(AimedPatch)->nbImigrant++;
    } else
      patch->nbPhilopat++;
    
  }//end while
}
// ----------------------------------------------------------------------------------------
// MigratePatch_AbsorbingBorder
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::MigratePatch_AbsorbingBorder (sex_t SEX, unsigned int LocalPatch)
{   
  Patch* patch = _popPtr->getPatch(LocalPatch);
  unsigned int AimedPatch; 
  
  while( patch->size(SEX, OFFSx) != 0 ) {
   
    do{
      AimedPatch = getMigrationPatchForward(SEX, LocalPatch);
    }while ( AimedPatch > _npatch );
    
    if(AimedPatch == _npatch) {
      _popPtr->recycle( patch->get(SEX, OFFSx, 0) );
      patch->remove(SEX, OFFSx, 0);
    } else
      _popPtr->move(SEX, OFFSx, LocalPatch, PDISPx, AimedPatch, 0);
    
    if(LocalPatch != AimedPatch) {
      patch->nbEmigrant++;
      if(AimedPatch < _npatch)
        _popPtr->getPatch(AimedPatch)->nbImigrant++;
    } else
      patch->nbPhilopat++;
    
  }//end while
}
// ----------------------------------------------------------------------------------------

//                             LCE_SeedDisp/

// ----------------------------------------------------------------------------------------
LCE_SeedDisp::LCE_SeedDisp () : LifeCycleEvent ("seed_disp","")
{
  ParamUpdater<LCE_SeedDisp> * updater =
  new ParamUpdater<LCE_SeedDisp> (&LCE_SeedDisp::setParameters);
  
//  cout << "calling LCE_SeedDisp()\n";
  
  set_event_name("seed_disp"); //this resets the paramset, erases the params added by LCE_Disperse_ConstDisp

  addParameters("seed_disp", updater);

//  get_paramset()->show_up();
}
// ------------------------------------------------------------------------------

//                             LCE_Disperse_EvolDisp/

// ----------------------------------------------------------------------------------------
LCE_Disperse_EvolDisp::LCE_Disperse_EvolDisp () 
: LifeCycleEvent("disperse_evoldisp",""), _fem_cost(-1.0), _mal_cost(-1.0),
_fdisp_trait_link(0), _mdisp_trait_link(0)
{
  ParamUpdater<LCE_Disperse_EvolDisp> * updater =
  new ParamUpdater<LCE_Disperse_EvolDisp> (&LCE_Disperse_EvolDisp::setParameters);

  LCE_Disperse_base::addParameters("dispersal", updater);
  
  add_parameter("dispersal_cost",DBL,false,true,0,1, updater);
  add_parameter("dispersal_cost_fem",DBL,false,true,0,1, updater);
  add_parameter("dispersal_cost_mal",DBL,false,true,0,1, updater);
  add_parameter("dispersal_fixed_trait",STR,false,false,0,0, updater);
  add_parameter("dispersal_fixed_rate",DBL,false,true,0,1, updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_EvolDisp::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_EvolDisp::setParameters ()
{
  //we do not call the LCE_Disperse_base::setParameters here because we don't use
  //dispersal matrices in input
  _npatch = _popPtr->getPatchNbr();
  
  _prefix = "dispersal";
  
  _disp_model = (int)_paramSet->getValue("dispersal_model");
  
  _disp_propagule_prob = _paramSet->getValue("dispersal_propagule_prob");
  
  if(_disp_model == -1) {
    error("dispersal model not specified!\n");
    return false;
  }
  
  if(_disp_model == 2 && _disp_propagule_prob == -1) {
    error("dispersal propagule probability is missing!\n");
    return false;
  }
  
  _fdisp_trait_link = _popPtr->getTraitIndex("fdisp");
  _mdisp_trait_link = _popPtr->getTraitIndex("mdisp");
  
  if(_paramSet->isSet("dispersal_cost")) {
    
    _fem_cost = _mal_cost = _paramSet->getValue("dispersal_cost");
    
  } else if(_paramSet->isSet("dispersal_cost_fem") && _paramSet->isSet("dispersal_cost_mal")) {
    
    _fem_cost = _paramSet->getValue("dispersal_cost_fem");
    
    _mal_cost = _paramSet->getValue("dispersal_cost_mal");
    
  } else {
    error("dispersal cost params are not set !\n");
    return false;
  }
  
  if(_paramSet->isSet("dispersal_fixed_trait")) {
    
    if(_paramSet->getArg("dispersal_fixed_trait").compare("female") == 0) {
      
      exec = &LCE_Disperse_EvolDisp::exec_evolmale;
      
    } else if(_paramSet->getArg("dispersal_fixed_trait").compare("male") == 0) {
      
      exec = &LCE_Disperse_EvolDisp::exec_evolfemale;
      
    } else {
      error("wrong argument value for \"dispersal_fixed_trait\"\n");
      return false;
    }
    
    _fixed_disp_rate = _paramSet->getValue("dispersal_fixed_rate");
    
  } else
    exec = &LCE_Disperse_EvolDisp::exec_evol2sex;
  
  switch ( getDispersalModel() ) {
    case 0:
    case 1:
      getAimedPatch = &LCE_Disperse_EvolDisp::Migrate_Island;
      break;
    case 2:
      getAimedPatch = &LCE_Disperse_EvolDisp::Migrate_Island_Propagule;
      break;
    case 3:
      getAimedPatch = &LCE_Disperse_EvolDisp::Migrate_SteppingStone1D;
      break;
    default: {
      /**@todo implement lattice dispersal models with evolving dispersal.*/
      error("\nDispersal model %i not yet implemented with evolving dispersal!\n",
            getDispersalModel());
      return false;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_EvolDisp::execute
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::execute ()
{  
#ifdef _DEBUG_
  message("LCE_Disperse_EvolDisp::execute (Patch nb: %i offsprg nbr: %i)\n"
          ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ));
#endif
  if( getDispersalModel() == 2 ) setPropaguleTargets();
  
  reset_counters();
  
  (this->*exec)();
  
  Patch *current_patch;  
  
  _npatch = _popPtr->getPatchNbr();
  
  for(unsigned int i = 0; i < _npatch; i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    //set coloniser counter
    if(current_patch->get_isExtinct()) 
      current_patch->nbKolonisers = current_patch->size(PDISPx);
    else 
      current_patch->nbKolonisers = -1;
    
  }//end_for
  
  swapPostDisp();
  
}
// ----------------------------------------------------------------------------------------
// exec_evolfemale
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::exec_evolfemale ()
{
  evoldisp(FEM, _fdisp_trait_link, _fem_cost);
  fixdisp(MAL, _fixed_disp_rate, _mal_cost);
}
// ----------------------------------------------------------------------------------------
// exec_evolmale
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::exec_evolmale ()
{
  evoldisp(MAL, _mdisp_trait_link, _mal_cost);
  fixdisp(FEM, _fixed_disp_rate, _fem_cost);
}
// ----------------------------------------------------------------------------------------
// exec_evol2sex
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::exec_evol2sex ()
{
  evoldisp(FEM, _fdisp_trait_link, _fem_cost);
  evoldisp(MAL, _mdisp_trait_link, _mal_cost);
}
// ----------------------------------------------------------------------------------------
// evoldisp
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::evoldisp (sex_t SEX, int trait_link, double cost)
{
  
  Patch *current_patch;
  unsigned int AimedPatch;
  
  for(unsigned int i = 0; i < _npatch; i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    while( current_patch->size(SEX, OFFSx) != 0 ) {
      
      if(RAND::Uniform() < *(double*)current_patch->get(SEX, OFFSx, 0)->getTraitValue(trait_link)) {
        //this one disperses
        AimedPatch = (this->*getAimedPatch)(i);
        
        current_patch->nbEmigrant++;
        
        if(RAND::Uniform() > cost) {
          //survives
          _popPtr->move(SEX, OFFSx, i, PDISPx, AimedPatch, 0);
          
          _popPtr->getPatch(AimedPatch)->nbImigrant++;
        } else {
          
          _popPtr->recycle( current_patch->get(SEX, OFFSx, 0) );
          
          current_patch->remove(SEX, OFFSx, 0);
        }
      } else {
        //no dispersal
        current_patch->move(SEX,OFFSx,PDISPx,0);
        current_patch->nbPhilopat++;
      }
      
    }//end while
  }//end for
}
// ----------------------------------------------------------------------------------------
// fixdisp
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::fixdisp (sex_t SEX, double rate, double cost)
{
  
  Patch *current_patch;
  unsigned int AimedPatch;
  
  for(unsigned int i = 0; i < _npatch; i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    while( current_patch->size(SEX, OFFSx) != 0 ) {
      
      if(RAND::Uniform() < rate) {
        //this one disperses
        AimedPatch = (this->*getAimedPatch)(i);
        
        current_patch->nbEmigrant++;
        
        if(RAND::Uniform() > cost) {
          //survives
          _popPtr->move(SEX, OFFSx, i, PDISPx, AimedPatch, 0);
          
          _popPtr->getPatch(AimedPatch)->nbImigrant++;
        } else {
          
          _popPtr->recycle( current_patch->get(SEX, OFFSx, 0) );
          
          current_patch->remove(SEX, OFFSx, 0);
        }
      } else {
        //no dispersal
        current_patch->move(SEX,OFFSx,PDISPx,0);
        current_patch->nbPhilopat++;
      }
      
    }//end while
  }//end for
}
// ----------------------------------------------------------------------------------------
// Migrate_Island_EvolDisp
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_EvolDisp::Migrate_Island(unsigned int home)
{
  unsigned int AimedPatch;
  //assign a Patch of arrival at random
  do{
    AimedPatch = RAND::Uniform(_npatch);
  }while(AimedPatch == home);
  
  return AimedPatch;
}
// ----------------------------------------------------------------------------------------
// Migrate_Island_Propagule_EvolDisp
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_EvolDisp::Migrate_Island_Propagule(unsigned int home)
{
  unsigned int AimedPatch, PropaguleTarget = getPropaguleTarget(home);
  
  if(!(RAND::Uniform() > getPropaguleProb()) )
    AimedPatch = PropaguleTarget;
  else
    do{
      AimedPatch = RAND::Uniform(_npatch);
    }while(AimedPatch == home || AimedPatch == PropaguleTarget);
  
  return AimedPatch;
}
// ----------------------------------------------------------------------------------------
// Migrate_SteppingStone1D_EvolDisp
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_EvolDisp::Migrate_SteppingStone1D(unsigned int home)
{
  int neighbours[2] = {(int)(home - 1),(int)(home + 1)};
  //if we are at one of the bound of the Patch array, we may migrate to the other bound:
  if(neighbours[0] < 0) neighbours[0] = _npatch -1;
  else if(neighbours[1] == (int)_npatch) neighbours[1] = 0;
  
  return(RAND::RandBool() ? neighbours[0] : neighbours[1]);
}

