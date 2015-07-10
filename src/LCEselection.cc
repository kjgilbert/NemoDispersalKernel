/**  $Id: LCEselection.cc,v 1.21 2015-04-01 14:25:16 fred Exp $
*
*  @file LCEselection.cc
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
*  created on @date 09.10.2007
*
*  @author fred
*/

#include <sstream>
#include "LCEselection.h"
#include "Uniform.h"
#include "utils.h"
#include "tstring.h"
#include "simenv.h"

// ------------------------------------------------------------------------------

//                             LCE_Selection_base/

// ----------------------------------------------------------------------------------------
// LCE_Selection_base::LCE_Selection_base
// ----------------------------------------------------------------------------------------
LCE_Selection_base::LCE_Selection_base ( ) : LifeCycleEvent("viability_selection", ""),
_selection_matrix(0), _gsl_selection_matrix(0), _diffs(0), _res1(0), _local_optima(0), _phe(0),
_selectTraitDimension(1), _base_fitness(1), _mean_fitness(0), _max_fitness(0), _scaling_factor(1), 
_is_local(0), _is_absolute(1), _eVariance(0), _getRawFitness(0), _getFitness(0), _stater(0)
{ 
  add_parameter("selection_trait",STR,false,false,0,0); //no updaters here, to keep it safe...
  
  ParamUpdater< LCE_Selection_base > * updater = 
    new ParamUpdater<LCE_Selection_base>(&LCE_Selection_base::set_fit_model);
  
  add_parameter("selection_fitness_model",STR,false,false,0,0, updater);
  
  updater = new ParamUpdater<LCE_Selection_base>(&LCE_Selection_base::set_sel_model);
  
  add_parameter("selection_model",STR,false,false,0,0, updater);
  add_parameter("selection_matrix",MAT,false,false,0,0, updater);
  add_parameter("selection_variance",DBL,false,false,0,0, updater);
  add_parameter("selection_correlation",DBL,false,false,0,0, updater);
  add_parameter("selection_trait_dimension",INT,false,false,0,0, updater);
  add_parameter("selection_base_fitness",DBL,false,true,0,1, updater);
  add_parameter("selection_lethal_equivalents",DBL,false,false,0,0, updater);
  add_parameter("selection_pedigree_F",MAT,false,false,0,0, updater);
  add_parameter("selection_randomize",BOOL,false,false,0,0, updater);
  add_parameter("selection_environmental_variance", DBL, false, false, 0, 0, updater);
  
  updater = new ParamUpdater<LCE_Selection_base>(&LCE_Selection_base::set_local_optima);
  add_parameter("selection_local_optima",DBL,false,false,0,0, updater);
  
  updater = new ParamUpdater<LCE_Selection_base>(&LCE_Selection_base::set_param_rate_of_change);
  add_parameter("selection_rate_environmental_change", DBL, false, false, 0, 0, updater);
  add_parameter("selection_std_rate_environmental_change", DBL, false, false, 0, 0, updater);
  add_parameter("selection_std_rate_set_at_generation", INT, false, false, 0, 0, updater);
  add_parameter("selection_std_rate_reference_patch", INT, false, false, 0, 0, updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::~LCE_Selection_base
// ----------------------------------------------------------------------------------------
LCE_Selection_base::~LCE_Selection_base ( )
{
  vector< TMatrix* >::iterator selIT = _selection_matrix.begin();
  for(; selIT != _selection_matrix.end(); ++selIT)
    if((*selIT)) delete (*selIT);
  _selection_matrix.clear();
  
  for(unsigned int i  = 0; i < _gsl_selection_matrix.size(); ++i)
    if(_gsl_selection_matrix[i]) gsl_matrix_free(_gsl_selection_matrix[i]);
  _gsl_selection_matrix.clear();
  
  if(_local_optima) delete _local_optima;
  if(_diffs) gsl_vector_free(_diffs);
  if(_res1) gsl_vector_free(_res1);
  if(_stater) delete _stater;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::loadStatServices
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::loadStatServices ( StatServices* loader )
{
  if(_stater == NULL) _stater = new LCE_SelectionSH(this);
  loader->attach(_stater);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::setParameters ( )
{   
  //selection may on more than on trait at a time
  vector< string >  traits = get_parameter("selection_trait")->getMultiArgs();
  
  _TraitIndices.clear();

  for(unsigned int i = 0; i < traits.size(); i++)  {

    if(_popPtr->getTraitIndex(traits[i].c_str()) == -1) {
      return error("cannot attach trait \"%s\" to life cycle event \"%s\", trait has not been initiated.\n",
            traits[i].c_str(), _event_name.c_str());
    } 
    else {    
      _TraitIndices.push_back(_popPtr->getTraitIndex(traits[i].c_str()));
    }
    
  }
  
  if(!set_fit_model()) return false;
  
  if(!set_sel_model()) return false;
  
  if(!set_local_optima()) return false;
  
  if(!set_param_rate_of_change()) return false;
  
  _mean_fitness = _max_fitness = 0;
  _scaling_factor = 1;
  resetCounters();

  return true; 
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::set_fit_model
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::set_fit_model()
{  
    
  if(_paramSet->isSet("selection_fitness_model")) {
    
    string fit_model = _paramSet->getArg("selection_fitness_model");
    
    if(fit_model.compare("absolute") == 0) {
      
      _is_local = false;
      _is_absolute = true;
      _getFitness = &LCE_Selection_base::getFitnessAbsolute;
      _setScalingFactor = &LCE_Selection_base::setScalingFactorAbsolute;
      
    } else if(fit_model.compare("relative_local") == 0) { 
      
      _is_local = true;      
      _is_absolute = false;
      _getFitness = &LCE_Selection_base::getFitnessRelative;
      _setScalingFactor = &LCE_Selection_base::setScalingFactorLocal;
      
    } else if(fit_model.compare("relative_global") == 0) {
            
      _is_local = false;
      _is_absolute = false;
      _getFitness = &LCE_Selection_base::getFitnessRelative;
      _setScalingFactor = &LCE_Selection_base::setScalingFactorGlobal;
      
    } else {
      error("Unknown fitness model \"%s\"", fit_model.c_str());
      return false;
    }  
  } //default case:
  else {
    _is_local = false;
    _is_absolute = true;
    _getFitness = &LCE_Selection_base::getFitnessAbsolute;
    _setScalingFactor = &LCE_Selection_base::setScalingFactorAbsolute;
  }

  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::set_sel_model
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::set_sel_model()
{  
  
  _selectTraitDimension = (int)get_parameter_value("selection_trait_dimension");
  
  _base_fitness = get_parameter_value("selection_base_fitness");

  if(get_parameter("selection_environmental_variance")->isSet())
    //the variable actually holds the standard dev...
    _eVariance =  sqrt(get_parameter_value("selection_environmental_variance"));
  else
    _eVariance = 0;
  
  _SelectionModels.clear();
  
  if(!_paramSet->isSet("selection_model")) {
    
    _getRawFitness.push_back(&LCE_Selection_base::getFitnessDirect);
    
    _SelectionModels.push_back("direct");
      
  } else {
    
    _SelectionModels = get_parameter("selection_model")->getMultiArgs();
  } //end_if isSet()

  if (_SelectionModels.size() != _TraitIndices.size()) {
    error("\"selection_trait\" and \"selection_model\" must have the same number of arguments.");
    return false;
  }

  
  string sel_model;
  
  _getRawFitness.clear();
  
  for(unsigned int t = 0; t < _SelectionModels.size(); t++) {

    sel_model = _SelectionModels[t];
  
    if(sel_model == "fix") {
      
      if(!get_parameter("selection_lethal_equivalents")->isSet()) {
        error("\"selection_lethal_equivalents\" parameter is missing with \"fix\" selection!\n");
        return false;
      } else
        _letheq = get_parameter_value("selection_lethal_equivalents");
      
      if(!get_parameter("selection_base_fitness")->isSet()) {
        warning("\"selection_base_fitness\" parameter is missing under fix selection model, setting it to 1.\n");
        _base_fitness = 1.0;
      }
      
      if(!get_parameter("selection_pedigree_F")->isSet()) {
      
        return error("\"selection_pedigree_F\" parameter is missing with \"fix\" selection!\n");
     
      } else {
        
        TMatrix tmp_mat;
        
        get_parameter("selection_pedigree_F")->getMatrix(&tmp_mat);
        
        if(tmp_mat.getNbCols() != 5) {
          return error("\"selection_pedigree_F\" must be an array of size 5.\n");
        }

        for(unsigned int i = 0; i < 5; i++)
          _Fpedigree[i] = tmp_mat.get(0,i);
      }
      
      for(unsigned int i = 0; i < 5; i++)
        _FitnessFixModel[i] = _base_fitness * exp( -_letheq * _Fpedigree[i] );
            
      _getRawFitness.push_back(&LCE_Selection_base::getFitnessFixedEffect);
      
    } else if(sel_model == "direct") {
      
       _getRawFitness.push_back(&LCE_Selection_base::getFitnessDirect);
      
    } else if(sel_model == "quadratic") {
      
      if(_selectTraitDimension == 1){
      
        if(!setSelectionMatrix()) return false; //this to set the selection variance params
        //now reset the fitness function:
         _getRawFitness.push_back(&LCE_Selection_base::getFitnessUnivariateQuadratic);
      
      } else return error("\"quadratic\" fitness model implemented for a single trait only.\n");
    

    } else if(sel_model == "gaussian") {
      
      if(!setSelectionMatrix()) return false;
      
    } else return error("wrong selection model, must be either \"fix\", \"direct\", or \"gaussian\".\n");
    
  }
  
  if(sel_model != "fix" && !_paramSet->isSet("selection_trait"))
    return error("trait under selection is not set, please add parameter \"selection_trait\"\n");

  return true;
}
// ----------------------------------------------------------------------------------------
// setSelectionMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::setSelectionMatrix()
{
  TMatrix tmp_mat; 
  unsigned int patchNbr = _popPtr->getPatchNbr();
  
  if(!get_parameter("selection_matrix")->isSet() && !get_parameter("selection_variance")->isSet()) 
    return error("\"selection_matrix\" or \"selection_variance\" must be set with selection model = \"gaussian\".\n");

  
  if(get_parameter("selection_variance")->isSet() && !get_parameter("selection_trait_dimension")->isSet())
    return error("parameter \"selection_trait_dimension\" is missing!\n");
   
  
  //clear the selection matrix container
  vector< TMatrix* >::iterator selIT = _selection_matrix.begin();
  
  for(; selIT != _selection_matrix.end(); ++selIT) if((*selIT)) delete (*selIT);  
  
  _selection_matrix.clear();
  
  
  if(get_parameter("selection_matrix")->isSet()) {
    
    //selection matrix provided, same selection surface in each patch
    _paramSet->getMatrix("selection_matrix", &tmp_mat);
    
    if(tmp_mat.getNbCols() != tmp_mat.getNbRows())
      return error("\"selection_matrix\" must be a square matrix!\n");

    _selectTraitDimension = tmp_mat.getNbCols();
    
    //we have one selection matrix per patch, copy it in the container for each patch
    for(unsigned int i = 0; i < patchNbr; ++i)
      _selection_matrix.push_back( new TMatrix(tmp_mat) );
    
  } else {
    //we have to check for spatial variation in variance and covariances
    _selectTraitDimension = (unsigned int)get_parameter_value("selection_trait_dimension");
    
    TMatrix var_spatmat, corr_spatmat;
    
    //setting variance spatial matrix:
    var_spatmat.reset(patchNbr, (unsigned)_selectTraitDimension);
    
    if(get_parameter("selection_variance")->isMatrix()) {
    
      _paramSet->getMatrix("selection_variance", &tmp_mat);

      if( !setSpatialMatrix("selection_variance","\"selection_trait_dimension\"", &tmp_mat, &var_spatmat,
                            (unsigned)_selectTraitDimension, patchNbr, _paramSet->isSet("selection_randomize") ) ) 
        return false;
    
    } else {
      
      var_spatmat.assign(get_parameter_value("selection_variance"));
    }
    
    //setting correlation spatial matrix:
    corr_spatmat.reset(patchNbr, (unsigned)_selectTraitDimension*(_selectTraitDimension-1)/2);
    
    if(get_parameter("selection_correlation")->isMatrix()) {
      
      _paramSet->getMatrix("selection_correlation", &tmp_mat);
      
      if( !setSpatialMatrix("selection_correlation","the num of correlation coefficients", &tmp_mat, &corr_spatmat, 
                            (unsigned)_selectTraitDimension*(_selectTraitDimension-1)/2, patchNbr, _paramSet->isSet("selection_randomize") ) )
        return false;
      
    } else {
      
      corr_spatmat.assign((get_parameter("selection_correlation")->isSet() ?
                           get_parameter_value("selection_correlation") : 0.0 ));
    }
    
    //set the selection matrix:
    tmp_mat.reset(_selectTraitDimension, _selectTraitDimension);
    double covar;
    unsigned int col;
    for( unsigned int p = 0; p < patchNbr; p++) {
      col = 0;
      for( int i = 0; i < _selectTraitDimension; i++) {
        tmp_mat.set(i, i, var_spatmat.get(p, i));
        for( int j = i+1; j < _selectTraitDimension; j++) {
          covar = corr_spatmat.get(p, col) * sqrt( var_spatmat.get(p, i) * var_spatmat.get(p, j) );
          tmp_mat.set(i, j, covar);
          tmp_mat.set(j, i, covar);
          col++;
        }
      }
      _selection_matrix.push_back(new TMatrix(tmp_mat));
    }
  }
  
  if(_selectTraitDimension > 1) {
    
    //selection on more than one trait
    
    //inversing the selection matrices:    
    if(_gsl_selection_matrix.size() != 0)
      for(unsigned int i = 0; i < _gsl_selection_matrix.size(); ++i)
        if(_gsl_selection_matrix[i] != NULL) gsl_matrix_free( _gsl_selection_matrix[i] );
    
    _gsl_selection_matrix.clear();
    
    for( unsigned int p = 0; p < patchNbr; p++) {
      _selection_matrix[p]->inverse();

      _gsl_selection_matrix.push_back( gsl_matrix_alloc(_selectTraitDimension, _selectTraitDimension) );
      
      _selection_matrix[p]->get_gsl_matrix(_gsl_selection_matrix[p]);
    }
    //allocate the vectors used by the fitness function:
    if(_diffs != NULL) gsl_vector_free(_diffs);
    _diffs = gsl_vector_alloc( _selectTraitDimension );
    
    if(_res1 != NULL) gsl_vector_free(_res1);
    _res1 = gsl_vector_alloc( _selectTraitDimension );
    
    if(_eVariance > 0)
       _getRawFitness.push_back(&LCE_Selection_base::getFitnessMultivariateGaussian_VE);
    else
       _getRawFitness.push_back(&LCE_Selection_base::getFitnessMultivariateGaussian);
    
  } else {
    //selection on one trait only
    _selection_variance.clear();
    for(unsigned int p = 0; p < patchNbr; p++)
      _selection_variance.push_back( _selection_matrix[p]->get(0, 0) );
    
    if(_eVariance > 0)
       _getRawFitness.push_back(&LCE_Selection_base::getFitnessUnivariateGaussian_VE);
    else
       _getRawFitness.push_back(&LCE_Selection_base::getFitnessUnivariateGaussian);
  }
  
  
  return true;
}

bool LCE_Selection_base::set_local_optima ()
{//set the traits' local optima, _selectTraitDimension must be set before that (= nbr of traits to select on)    
  string model = _paramSet->getArg("selection_model");
  
  if( model == "fix"  || model == "direct") return true;
  
  if( (model == "gaussian" || model == "quadratic")
     && !get_parameter("selection_local_optima")->isSet()) {
    
    return error("parameter \"selection_local_optima\" must be set to have Gaussian or quadratic selection.\n");
    
  } else {
    
    TMatrix tmp_mat;
    
    if(_local_optima == 0) _local_optima = new TMatrix();
    
    _local_optima->reset(_popPtr->getPatchNbr(), _selectTraitDimension);
    
    _paramSet->getMatrix("selection_local_optima", &tmp_mat);
    
    return setSpatialMatrix("selection_local_optima", "\"selection_trait_dimension\"", &tmp_mat, _local_optima, 
                            _selectTraitDimension, _popPtr->getPatchNbr(), _paramSet->isSet("selection_randomize"));
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::set_param_rate_of_change
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::set_param_rate_of_change ()
{
  _rate_of_change_is_std = false;
  _do_change_local_opt = false;
  _setNewLocalOptima = 0;
  _rate_of_change_local_optima.reset();
  
  
  if (!get_parameter("selection_rate_environmental_change")->isSet() 
      && !get_parameter("selection_std_rate_environmental_change")->isSet() )
    return true;
  
  if (get_parameter("selection_rate_environmental_change")->isSet() 
      && get_parameter("selection_std_rate_environmental_change")->isSet() ) {
    return error("both \"selection_rate_environmental_change\" and \"selection_std_rate_environmental_change\" are set, need only one.\n");
  }
  
  TMatrix tmpMat;
  
  if (get_parameter("selection_rate_environmental_change")->isSet()) {
    
    if(!get_parameter("selection_rate_environmental_change")->isMatrix()) {
      
      double val = get_parameter_value("selection_rate_environmental_change");
      
      tmpMat.reset(1, _selectTraitDimension);
      tmpMat.assign(val);
      
    } else {
      get_parameter("selection_rate_environmental_change")->getMatrix(&tmpMat);
    }
    
    _setNewLocalOptima = &LCE_Selection_base::changeLocalOptima;
  
  } else if (get_parameter("selection_std_rate_environmental_change")->isSet()){
    
    if(!get_parameter("selection_std_rate_environmental_change")->isMatrix()) {
      
      double val = get_parameter_value("selection_std_rate_environmental_change");
      
      tmpMat.reset(1, _selectTraitDimension);
      tmpMat.assign(val);
      
    } else {
      get_parameter("selection_std_rate_environmental_change")->getMatrix(&tmpMat);
    }
    
    _rate_of_change_is_std = true;
    
    if(get_parameter("selection_std_rate_set_at_generation")->isSet())
      _set_std_rate_at_generation = (unsigned int)get_parameter_value("selection_std_rate_set_at_generation");
    else
      _set_std_rate_at_generation = 1;
    
    //check if phenotypic SD is to be computed in a single reference patch
    //is -1 if parameter not set, which corresponds to the whole population then
    _std_rate_reference_patch = get_parameter_value("selection_std_rate_reference_patch");
    
    _setNewLocalOptima = &LCE_Selection_base::set_std_rate_of_change;
  }

  if(tmpMat.nrows() != 1)
    return error("The matrix of rate of change in local optima must be an array (nrows=1)\n");
  
  if((int)tmpMat.ncols() > _selectTraitDimension)
    return error("The matrix of rate of change in local optima must have at most as many columns as dimensions of the trait(s) under selection\n");
  
  _do_change_local_opt = true;

  _rate_of_change_local_optima.reset(1, _selectTraitDimension);
  
  unsigned int nVal = tmpMat.ncols();

  //copy values, with pattern propagation
  for (int i = 0; i < _selectTraitDimension; ++i) {
    _rate_of_change_local_optima.set(0, i, tmpMat.get(0, i % nVal));
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessUnivariateQuadratic
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessUnivariateQuadratic ( Individual* ind, unsigned int patch, unsigned int trait )
{
  register double res2, diff;
  
  _phe = (double*)ind->getTraitValue(trait);
  
  diff = _phe[0] - _local_optima->get(patch, 0);
  
  res2 = diff*diff / _selection_variance[patch];
  
  return 1 - res2;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessUnivariateGaussian
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessUnivariateGaussian ( Individual* ind, unsigned int patch, unsigned int trait )
{
  register double res2, diff;
  
  _phe = (double*)ind->getTraitValue(trait);
  
  diff = _phe[0] - _local_optima->get(patch, 0);
  
  res2 = diff*diff / _selection_variance[patch];

  return exp( -0.5 * res2 );
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessMultivariateGaussian
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessMultivariateGaussian ( Individual* ind, unsigned int patch, unsigned int trait )
{
  register double res2;
  
//  if(!ind) fatal("passing NULL ind ptr to LCE_Selection_base::getFitnessMultivariateGaussian!!!\n");
  
  _phe = (double*)ind->getTraitValue(trait);
  
  for( int i = 0; i < _selectTraitDimension; i++)
    gsl_vector_set(_diffs, i, _phe[i] - _local_optima->get(patch, i));
  
  //(diff)T * W * diff:  
  //right partial product:
  gsl_blas_dsymv(CblasUpper, 1.0, _gsl_selection_matrix[patch], _diffs, 0.0, _res1);
  //left product:
  gsl_blas_ddot(_diffs, _res1, &res2);
  
  return exp( -0.5 * res2 );
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessUnivariateGaussian_VE
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessUnivariateGaussian_VE ( Individual* ind, unsigned int patch, unsigned int trait )
{
  register double res2, diff;
  
  _phe = (double*)ind->getTraitValue(trait);
  
  //add the environmental variance here:
  diff = _phe[0] + RAND::Gaussian(_eVariance) - _local_optima->get(patch, 0);
  
  res2 = diff*diff / _selection_variance[patch];
  
  return exp( -0.5 * res2 );
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessMultivariateGaussian_VE
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessMultivariateGaussian_VE ( Individual* ind, unsigned int patch, unsigned int trait )
{
  register double res2;
  
  _phe = (double*)ind->getTraitValue(trait);
  
  //add the environmental variance here:
  for( int i = 0; i < _selectTraitDimension; i++)
    gsl_vector_set(_diffs, i, _phe[i]  + RAND::Gaussian(_eVariance) - _local_optima->get(patch, i));
  
  //(diff)T * W * diff:  
  //right partial product:
  gsl_blas_dsymv(CblasUpper, 1.0, _gsl_selection_matrix[patch], _diffs, 0.0, _res1);
  //left product:
  gsl_blas_ddot(_diffs, _res1, &res2);
  
  return exp( -0.5 * res2 );
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getMeanFitness
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessAbsolute (Individual* ind, unsigned int patch)
{
  double fitness = 1;
  //at this point, _getRawFitness.size() == _TraitIndices.size()
  //and we assume that the functions are "aligned" with the traits
  if((SIMenv::getCurrentGeneration() >= 5000 && SIMenv::getCurrentGeneration() <= 7500) && SIMenv::getCurrentGeneration() % 25 == 0) cout << SIMenv::getCurrentGeneration() << " " << patch << " ";
  
  for(unsigned int i = 0; i < _getRawFitness.size(); i++)
  {
    fitness *= (this->*_getRawFitness[i])(ind, patch, _TraitIndices[i]);
    if((SIMenv::getCurrentGeneration() >= 5000 && SIMenv::getCurrentGeneration() <= 7500) && SIMenv::getCurrentGeneration() % 25 == 0) cout << (this->*_getRawFitness[i])(ind, patch, _TraitIndices[i]) << " ";
  }
  if((SIMenv::getCurrentGeneration() >= 5000 && SIMenv::getCurrentGeneration() <= 7500) && SIMenv::getCurrentGeneration() % 25 == 0) cout << fitness << endl;
  return fitness;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getMeanFitness
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getMeanFitness (age_idx age)
{
  double mean = 0;
  Patch *patch;
//  age_idx age = (AGE == ADULTS ? ADLTx : OFFSx);
  
  for(unsigned int i = 0, npatch = _popPtr->getPatchNbr(); i < npatch; i++) {
    patch = _popPtr->getPatch(i);
    for(unsigned int j = 0, size = patch->size(FEM, age); j < size; j++)
      mean += getFitness( patch->get(FEM, age, j), i);
    for(unsigned int j = 0, size = patch->size(MAL, age); j < size; j++)
      mean += getFitness( patch->get(MAL, age, j), i);
  }
  return mean/_popPtr->size(age);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getMeanPatchFitness
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getMeanPatchFitness (age_idx age, unsigned int p)
{
  double mean = 0;
  Patch *patch = _popPtr->getPatch(p);
  
  for(unsigned int j = 0, size = patch->size(FEM, age); j < size; j++)
    mean += getFitness( patch->get(FEM, age, j), p);
  
  for(unsigned int j = 0, size = patch->size(MAL, age); j < size; j++)
    mean += getFitness( patch->get(MAL, age, j), p);
  
  return mean/patch->size(age);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::setScalingFactorLocal
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::setScalingFactorLocal (age_idx age, unsigned int p)
{
  _scaling_factor = 1; //this to have the raw mean fitness below
  _scaling_factor = 1.0/getMeanPatchFitness(age, p);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::setScalingFactorGlobal
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::setScalingFactorGlobal (age_idx age, unsigned int p)
{
  if(p != 0) return; //we compupte the mean fitness only once
  _scaling_factor = 1;
  _scaling_factor = 1.0/getMeanFitness(age);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::execute
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::execute ()
{
  Patch * patch;
  unsigned int popSize = _popPtr->size(OFFSPRG);
  
  resetCounters();

  if(_do_change_local_opt) {
  
    if(_popPtr->getCurrentGeneration() == 1) {
      set_local_optima(); //reset local optima to initial values
      if(_rate_of_change_is_std) //reset rate of change relative to SD of that replicate
        _setNewLocalOptima = &LCE_Selection_base::set_std_rate_of_change;
    }
    
    (this->*_setNewLocalOptima)();
  }
  
  for(unsigned int p = 0; p < _popPtr->getPatchNbr(); p++) {
    
    (this->*_setScalingFactor)(OFFSx, p);
    
    patch = _popPtr->getPatch(p);
    
    doViabilitySelection(FEM, OFFSx, patch, p);
    
    doViabilitySelection(MAL, OFFSx, patch, p);
  }

  setMeans(popSize);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::doViabilitySelection
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::doViabilitySelection (sex_t SEX, age_idx AGE, Patch* patch, unsigned int p)
{
  register Individual* ind;
  register double fitness;
  register unsigned int cat;
  
  for(unsigned int i = 0; i < patch->size(SEX, AGE); i++) {
  
    ind = patch->get(SEX, AGE, i);
    cat = ind->getPedigreeClass();
    fitness = getFitness( ind, p);
       
    _fitness[cat] += fitness;
    _ind_cntr[cat]++;
    
    if(RAND::Uniform() > fitness ) {
      //this individual dies
      patch->remove(SEX, AGE, i);
      
      _popPtr->recycle(ind);
    
      i--;
      
    } //else; this individual stays in the patch
    else {
      _survival[cat]++;
    }

  }
}
// ----------------------------------------------------------------------------------------
// LCE_Selection::set_std_rate_of_change
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::set_std_rate_of_change()
{
  if(_popPtr->getCurrentGeneration() == _set_std_rate_at_generation) {
    
    //reset rate_of_change matrix in each new replicate:
    TMatrix tmpMat;
    
    if(!get_parameter("selection_std_rate_environmental_change")->isMatrix()) {
      
      double val = get_parameter_value("selection_std_rate_environmental_change");
      
      tmpMat.reset(1, _selectTraitDimension);
      tmpMat.assign(val);
      
    } else {
      get_parameter("selection_std_rate_environmental_change")->getMatrix(&tmpMat);
    }
    
    unsigned int nVal = tmpMat.ncols();

    for (int i = 0; i < _selectTraitDimension; ++i) {
      _rate_of_change_local_optima.set(0, i, tmpMat.get(0, i % nVal));
    }
    
    double* SD = new double [_selectTraitDimension];
    
    for (int i = 0; i < _selectTraitDimension; ++i) {
      SD[i] = 0;
    }
    
    // check if SD is taken in a reference patch instead of the whole pop
    if (_std_rate_reference_patch > -1) {
      
      if( _popPtr->getPatch(_std_rate_reference_patch)->size(OFFSx) != 0 )
        addPhenotypicSD(_std_rate_reference_patch, SD);
    
    } else {
      // compute SD as the mean within-patch SD
      
      unsigned int cnt = 0;
      
      //get SD only in extant demes, to avoid nans
      for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {
        if (_popPtr->getPatch(i)->size(OFFSx) != 0 ) {
          cnt++;
          addPhenotypicSD(i, SD);
        }
      }
      
      //compute mean within-patch phenotypic standard deviation:
      for (int i = 0; i < _selectTraitDimension; ++i) {
        SD[i]/= cnt;
        
      }
    } //end if
    
    //multiply the per-trait rates of change by SD:
    for (int i = 0; i < _selectTraitDimension; ++i){
      _rate_of_change_local_optima.multi(0, i, SD[i]);
    }
    
    //log the rates in the simulation log file:
    SIMenv::MainSim->_FileServices.log("#selection_rate_environmental_change " +
                                       _rate_of_change_local_optima.to_string());
    
    //compute the change of local optima for current generation:
    changeLocalOptima();
    
    //now reset the function pointer to changeLocalOptima() for next generation:
    _setNewLocalOptima = &LCE_Selection_base::changeLocalOptima;
    
    delete [] SD;
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Selection::set_std_rate_of_change
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::addPhenotypicSD (unsigned int deme, double *stDev)
{
  Patch *patch = _popPtr->getPatch(deme);
  Individual* ind;
  unsigned int table_size = patch->size(OFFSx);
  
  double **phenot = new double* [_selectTraitDimension];

  for (int i = 0; i < _selectTraitDimension; ++i) {
    
    phenot[i] = new double [table_size];
    
  }
  
  unsigned int pos = 0;
  
  for (unsigned int j = 0; j < patch->size(FEM, OFFSx) && pos < table_size; ++j) {
    
    ind = patch->get(FEM, OFFSx, j);
    _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);
    
    for (int i = 0; i < _selectTraitDimension; ++i)
      phenot[i][pos] = _phe[i];
    
    pos++;
    
  }
  
  for (unsigned int j = 0; j < patch->size(MAL, OFFSx) && pos < table_size; ++j) {
    
    ind = patch->get(MAL, OFFSx, j);
    _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);
    
    for (int i = 0; i < _selectTraitDimension; ++i)
      phenot[i][pos] = _phe[i];
    
    pos++;
    
  }
    
  assert(pos == table_size);
  
  for (int i = 0; i < _selectTraitDimension; ++i) {
    stDev[i] += sqrt( my_variance_with_fixed_mean( phenot[i], table_size, my_mean(phenot[i], table_size) ) );
    delete [] phenot[i];
  }
  
  delete [] phenot;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection::changeLocalOptima
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::changeLocalOptima ()
{
  for (int i = 0; i < _selectTraitDimension; ++i)
    for (unsigned int j = 0; j < _local_optima->nrows(); ++j) {
      _local_optima->plus(j, i, _rate_of_change_local_optima.get(0, i));
    }
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::setStatRecorders
// ----------------------------------------------------------------------------------------
bool LCE_SelectionSH::setStatRecorders (string& token)
{
  string age_tag = token.substr(0,token.find_first_of("."));
  string sub_token;
  age_t AGE = ALL;
  
  if (age_tag.size() != 0 && age_tag.size() != string::npos) {
    
    if (age_tag == "adlt") AGE = ADULTS;
    
    else if (age_tag == "off") AGE = OFFSPRG;
    
    else age_tag = "";
    
  } else {
    age_tag = "";
  }
  
  if (age_tag.size() != 0) 
    sub_token = token.substr(token.find_first_of(".") + 1, string::npos);
  else
    sub_token = token;
  
  if(sub_token == "fitness") {
    add("Mean population fitness","fitness.mean",ALL,0,0,&LCE_SelectionSH::getMeanFitness,0,0,0);
    add("Mean population fitness","fitness.outb",ALL,0,0,0,&LCE_SelectionSH::getFitness,0,0);
    add("Mean population fitness","fitness.outw",ALL,1,0,0,&LCE_SelectionSH::getFitness,0,0);
    add("Mean population fitness","fitness.hsib",ALL,2,0,0,&LCE_SelectionSH::getFitness,0,0);
    add("Mean population fitness","fitness.fsib",ALL,3,0,0,&LCE_SelectionSH::getFitness,0,0);
    add("Mean population fitness","fitness.self",ALL,4,0,0,&LCE_SelectionSH::getFitness,0,0);
  } else if(sub_token == "survival") {
    add("Mean offspring survival","survival.outb",ALL,0,0,0,&LCE_SelectionSH::getSurvival,0,0);
    add("Mean offspring survival","survival.outw",ALL,1,0,0,&LCE_SelectionSH::getSurvival,0,0);
    add("Mean offspring survival","survival.hsib",ALL,2,0,0,&LCE_SelectionSH::getSurvival,0,0);
    add("Mean offspring survival","survival.fsib",ALL,3,0,0,&LCE_SelectionSH::getSurvival,0,0);
    add("Mean offspring survival","survival.self",ALL,4,0,0,&LCE_SelectionSH::getSurvival,0,0);
  } else if(sub_token == "fitness.prop") {
    add("Proportion of b/n demes outbreds","prop.outb",ALL,0,0,0,&LCE_SelectionSH::getPedProp,0,0);
    add("Proportion of w/n demes outbreds","prop.outw",ALL,1,0,0,&LCE_SelectionSH::getPedProp,0,0);
    add("Proportion of half-sib crossings","prop.hsib",ALL,2,0,0,&LCE_SelectionSH::getPedProp,0,0);
    add("Proportion of full-sib crossings","prop.fsib",ALL,3,0,0,&LCE_SelectionSH::getPedProp,0,0);
    add("Proportion of selfed progeny","prop.self",ALL,4,0,0,&LCE_SelectionSH::getPedProp,0,0);
    
  } else if(sub_token == "fitness.patch") {
    
    addMeanPerPatch(AGE);  
    
  } else if(sub_token == "fitness.var.patch") {
    
    addVarPerPatch(AGE);
    
  } else return false;
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::addMeanPerPatch
// ----------------------------------------------------------------------------------------
void LCE_SelectionSH::addMeanPerPatch (age_t AGE)
{
  unsigned int patchNbr = _pop->getPatchNbr();
  
  if (AGE == ALL) {
    addMeanPerPatch(ADULTS);
    addMeanPerPatch(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off."); //at this stage, AGE != ALL
  string name = suffix + "W.avg.p";
  string long_name = "Mean fitness in patch ";
  string patch;

  void (LCE_SelectionSH::* setter) (void) = (AGE == ADULTS ?
                                             &LCE_SelectionSH::setAdultTable : 
                                             &LCE_SelectionSH::setOffsprgTable);

  unsigned int int_agex = static_cast<age_idx> ((AGE == ADULTS ? ADLTx : OFFSx));

  //first patch, gets the data table setter:
  add(long_name + "1", name + "1", AGE, 0, int_agex, 
      0,0,&LCE_SelectionSH::getMeanPatchFitness, setter);

  for(unsigned int p = 1; p < patchNbr; p++) {
    patch = tstring::int2str(p+1);
    add(long_name + patch, name + patch, AGE, p, int_agex,
        0,0,&LCE_SelectionSH::getMeanPatchFitness,0);
  }
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::addVarPerPatch
// ----------------------------------------------------------------------------------------
void LCE_SelectionSH::addVarPerPatch (age_t AGE)
{
  unsigned int patchNbr = _pop->getPatchNbr();
  
  if (AGE == ALL) {
    addVarPerPatch(ADULTS);
    addVarPerPatch(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off."); //at this stage, AGE != ALL
  string name = suffix + "W.var.p";
  string long_name = "Var fitness in patch ";
  string patch;

  void (LCE_SelectionSH::* setter) (void) = (AGE == ADULTS ?
                                             &LCE_SelectionSH::setAdultTable : 
                                             &LCE_SelectionSH::setOffsprgTable);
  
  unsigned int int_agex = static_cast<age_idx> ((AGE == ADULTS ? ADLTx : OFFSx));
  
  //first patch, gets the data table setter:
  add(long_name + "1", name + "1", AGE, 0, int_agex, 
      0,0,&LCE_SelectionSH::getVarPatchFitness, setter);
  
  for(unsigned int p = 1; p < patchNbr; p++) {
    patch = tstring::int2str(p+1);
    add(long_name + patch, name + patch, AGE, p, int_agex,
        0,0,&LCE_SelectionSH::getVarPatchFitness,0);
  }
}
// ----------------------------------------------------------------------------------------
// setDataTable
// ----------------------------------------------------------------------------------------
void LCE_SelectionSH::setDataTable(age_t AGE) 
{
  if(_table_set_age == AGE 
     && _table_set_gen == _pop->getCurrentGeneration()
     && _table_set_repl == _pop->getCurrentReplicate()
     ) return;
  
  unsigned int patchNbr = _pop->getPatchNbr();

  if(_phenoTable.size() != patchNbr) {
    if(_phenoTable.size() < patchNbr) {
      while (_phenoTable.size() < patchNbr)
        _phenoTable.push_back(vector<double>());
    } else {
      while (_phenoTable.size() > patchNbr) {
        _phenoTable.pop_back();
      }
    }
  }

  for(unsigned int i = 0; i < patchNbr; ++i) {
    _phenoTable[i].assign(_pop->size(AGE, i),0);
  }
    
  Patch* patch;
  
  age_idx age = (AGE == ADULTS ? ADLTx : OFFSx);
  
  for(unsigned int i = 0, n; i < patchNbr; i++) {
    
    
    patch = _pop->getPatch(i);
    
    
    if( !_SHLinkedEvent->_is_absolute ) {
      if(_SHLinkedEvent->_is_local)
        _SHLinkedEvent->setScalingFactorLocal(age, i);
      else
        _SHLinkedEvent->setScalingFactorGlobal(age, i);
    }
    
    n = 0;
    
    for(unsigned int j = 0, size = patch->size(FEM, age); 
        j < size && n < _phenoTable[i].size(); 
        j++)
    {
      _phenoTable[i][n++] = _SHLinkedEvent->getFitness( patch->get(FEM, age, j), i);
    }
    
    for(unsigned int j = 0, size = patch->size(MAL, age); 
        j < size && n < _phenoTable[i].size(); 
        j++)
    {
      _phenoTable[i][n++] = _SHLinkedEvent->getFitness( patch->get(MAL, age, j), i);
    }
     
    
    if (n != _phenoTable[i].size()) {
      fatal("problem while recording fitness trait values; size counter doesn't match table size.\n");
    }
  }
  
  _table_set_age  = AGE;
  _table_set_gen  = _pop->getCurrentGeneration();
  _table_set_repl = _pop->getCurrentReplicate();
  
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::getMeanPatchFitness
// ----------------------------------------------------------------------------------------
double LCE_SelectionSH::getMeanPatchFitness (unsigned int i, unsigned int int_agex)
{
  age_idx age = static_cast<age_idx> (int_agex);
  unsigned int patch_size = _pop->getPatchPtr(i)->size(age);

  assert(patch_size == _phenoTable[i].size());
  
  if(patch_size == 0) return (nanf("NULL"));
  
  return getMeanPatchFitness(i);
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::getMeanPatchFitness
// ----------------------------------------------------------------------------------------
double LCE_SelectionSH::getMeanPatchFitness (unsigned int i)
{
  double mean = 0;
  unsigned int size = _phenoTable[i].size();
  
  for(unsigned int j = 0; j < size; j++)
    mean += _phenoTable[i][j];
  
  return mean/size;
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::getVarPatchFitness
// ----------------------------------------------------------------------------------------
double LCE_SelectionSH::getVarPatchFitness (unsigned int i, unsigned int int_agex)
{
  age_idx age = static_cast<age_idx> (int_agex);
  unsigned int patch_size = _pop->getPatchPtr(i)->size(age);
  
  assert(patch_size == _phenoTable[i].size());
  
  if(patch_size == 0) return nanf("NULL");
  
  double mean = getMeanPatchFitness(i);
  double var = 0;
  
  for(unsigned int j = 0; j < patch_size; j++) {
    
    var += pow(_phenoTable[i][j] - mean, 2.0);
  }
    
  return var/patch_size;
}
