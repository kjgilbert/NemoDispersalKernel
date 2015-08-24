/**  $Id: LCEselection.h,v 1.13 2015-04-01 14:15:04 fred Exp $
*
*  @file LCEselection.h
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

#ifndef LCE_SELECTION_H
#define LCE_SELECTION_H

#include <cmath>
#include "lifecycleevent.h"
#ifdef HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#endif

class LCE_SelectionSH;
//CLASS LCE_Selection_base
//
/**Base class performing (viability) selection on an arbitrary trait.
   Implements several types of selection models, from fixed inbreeding depression-type model to 
   Gaussian stabilizing selection on an adaptive landscape.
   This class is the base-class for the composite LCEs breed_selection and breed_selection_disperse.
*/
class LCE_Selection_base : public virtual LifeCycleEvent
{
  ///@name Gaussian selection variables
  ///@{
  
  vector< TMatrix* > _selection_matrix;
  vector< gsl_matrix* > _gsl_selection_matrix;
  gsl_vector *_diffs, *_res1;
  TMatrix *_local_optima;
  TMatrix _rate_of_change_local_optima;
  bool _do_change_local_opt, _rate_of_change_is_std;
  unsigned int _set_std_rate_at_generation;
  int _std_rate_reference_patch;
  /**Array to temporarily hold the phenotypic values of an individual, is never allocated.*/
  double *_phe;
  
  ///@}
  
  ///@name Fixed selection (lethal-equivalents-based) parameters
  ///@{
  
  /**Array of pedigree values used in the fixed-selection model.*/
  double _Fpedigree[5];
  
  /**Absolute fitness values of the five pedigree class for the fixed selection model (lethal equivalents model).
    Fitness is computed as:  _FitnessFixModel[i] = _base_fitness * exp( -_letheq * _Fpedigree[i] ); */
  double _FitnessFixModel[5];
  
  ///@}
  
  /**The StatHandler associated class is a friend.*/
  friend class LCE_SelectionSH;
  friend class LCE_Breed_Selection;
  
protected:
  
  double _letheq, _base_fitness;
  double _mean_fitness, _max_fitness, _scaling_factor;
  bool _is_local, _is_absolute;
  /**Fitness counters, one for each pedigree class.*/
  double _fitness[5], _survival[5], _ind_cntr[5]; 

  ///@name Quantitative traits parameters
  ///@{
  
  /**Number of quantitative traits under selection.*/
  int _selectTraitDimension;
  
  /**Patch-specific selection variance.*/
  vector< double > _selection_variance;
  
  /**Evironmental variance.*/
  double _eVariance;
  ///@}
  
  /**The indices of the traits under selection.*/
  vector< unsigned int > _TraitIndices;
  
  /**The selection models associated with each trait under selection.*/
  vector< string > _SelectionModels;  
  
  ///@name Pointers to functions, set in LCE_Selection_base::setParameters
  ///@{
  vector< double (LCE_Selection_base::* ) (Individual*, unsigned int, unsigned int) > _getRawFitness;
  double (LCE_Selection_base::* _getFitness) (Individual*, unsigned int);
  void (LCE_Selection_base::* _setScalingFactor) (age_idx, unsigned int);
  void (LCE_Selection_base::* _setNewLocalOptima) (void);
  ///@}
  
  LCE_SelectionSH* _stater;
    
public:
    
  LCE_Selection_base ( );
  
  virtual ~LCE_Selection_base ( );
  
  bool setSelectionMatrix();
  bool setLocalOptima ();  
  bool set_sel_model  ();
  bool set_fit_model  ();
  bool set_local_optima ();
  bool set_param_rate_of_change();
  
  void set_std_rate_of_change();
  void addPhenotypicSD (unsigned int deme, double *stDev);
  void changeLocalOptima ();

  /**Resets the fitness counters.*/
  void resetCounters  ()
  {
    for(unsigned int i = 0; i < 5; i++) {
      _fitness[i] = 0;
      _survival[i] = 0;
      _ind_cntr[i] = 0;
    }
  }
  /**Computes the average fitness of each pedigree class.
     Called after selection has been done in the whole population.
     @param tot_ind count of the total number of individuals produced*/
  void setMeans       (unsigned int tot_ind)
  {
    _mean_fitness = 0;
    for(unsigned int i = 0; i < 5; i++) {
      _mean_fitness += _fitness[i];
      _fitness[i] /= _ind_cntr[i];
      _survival[i] /= _ind_cntr[i];
      _ind_cntr[i] /= tot_ind;
    }
    _mean_fitness /= tot_ind;
  }
  
  
  /**Computes the mean fitness of the whole population for a given age class.
    @param age the age-class index*/
  double getMeanFitness         (age_idx age);
  
  /**Computes the mean fitness in a given patch for a given age class.
    @param age the age-class index
    @param p the patch index*/
  double getMeanPatchFitness    (age_idx age, unsigned int p);
  
  /**Sets the _mean_fitness variable to the value of the mean population fitness.
    @param age the age-class index*/
  double setMeanFitness         (age_idx age) { return (_mean_fitness = getMeanFitness(age));}

  /**Calls the fitness function according to the fitness model.
    The fitness model can be "absolute", "relative_local" or "relative_global".
    @param ind the focal indvidual, we want to know its fitness
    @param patch the index of the patch of the focal individual*/
  double getFitness (Individual* ind, unsigned int patch)
  {
    return (this->*_getFitness)(ind, patch);
  }
  
  /**Returns the fitness of an individual in the fixed selection model.*/
  double getFitnessFixedEffect (Individual* offsprg, unsigned int patch, unsigned int trait)
  {
    return _FitnessFixModel[ offsprg->getPedigreeClass() ];
  }
  
  /**Returns the fitness of an individual following the direct selection model.*/
  double getFitnessDirect (Individual* offsprg, unsigned int patch, unsigned int trait)
  { 
    return *(double*)offsprg->getTraitValue(trait); 
  }
  /**Quadratic fitness surface, approximates the Gaussian model for weak selection and/or small deviation from the optimum.*/
  double getFitnessUnivariateQuadratic ( Individual* ind, unsigned int patch, unsigned int trait);

  /**Returns the fitness of an individual following the Gaussian selection model with one trait under selection.*/
  double getFitnessMultivariateGaussian (Individual* offsprg, unsigned int patch, unsigned int trait);
  
  /**Returns the fitness of an individual following the Gaussian selection model with several traits under selection.*/
  double getFitnessUnivariateGaussian (Individual* offsprg, unsigned int patch, unsigned int trait);
  
  /**Returns the fitness of an individual following the Gaussian selection model with one trait under selection and environmental variance.*/
  double getFitnessMultivariateGaussian_VE (Individual* offsprg, unsigned int patch, unsigned int trait);
  
  /**Returns the fitness of an individual following the Gaussian selection model with several traits under selection and environmental variance.*/
  double getFitnessUnivariateGaussian_VE (Individual* offsprg, unsigned int patch, unsigned int trait);
  
  
  /**Returns the raw fitness of the individual, without adjustment (absolute fitness).
    Calls the fitness function according to the right selection model.
    @param ind the focal indvidual, we want to know its fitness
    @param patch the index of the patch of the focal individual*/
  double getFitnessAbsolute (Individual* ind, unsigned int patch);
  
  /**Returns the relative fitness of the individual, adjusted by a scaling factor.
    The scaling factor is set according to the relative fitness model (local or global)
    outside of this function.
    If the scaling factor is one, this is the absolute fitness of the individual.
    Calls the fitness function according to the right selection model.
    @param ind the focal indvidual, we want to know its fitness
    @param patch the index of the patch of the focal individual*/
  double getFitnessRelative (Individual* ind, unsigned int patch)
  {
    return getFitnessAbsolute(ind, patch) * _scaling_factor;
  }
  
  /**Sets the fitness scaling factor equal to the inverse of the mean local patch fitness.
    @param age the age-class index
    @param p the focal patch index*/
  void   setScalingFactorLocal  (age_idx age, unsigned int p);
  
  /**Sets the fitness scaling factor equal to the inverse of the mean population fitness.
    Function exits on p != 0; the mean population fitness is computed only once in the execute()
    procedure, at the beginning of the patch loop.
    @param age the age-class index
    @param p the focal patch index*/
  void   setScalingFactorGlobal (age_idx age, unsigned int p);
  
  /**Resets the fitness scaling factor equal to one.
    @param age the age-class index
    @param p the focal patch index*/
  void   setScalingFactorAbsolute (age_idx age, unsigned int p)
  { _scaling_factor = 1; }
  
  /**Selectively removes individuals in the population depending on their fitness.
    Calls the fitness function for each individual. Updates the fitness counters.
    The fitness scaling factor is set outside this function, in the execute() procedure.
    @param SEX the sex class
    @param AGE the age-class index
    @param patch the focal patch
    @param p the focal patch index*/
  void   doViabilitySelection   (sex_t SEX, age_idx AGE, Patch* patch, unsigned int p);
  
  ///@name Implementations
  ///@{
  virtual bool setParameters ();
  virtual void execute ();
  virtual void loadStatServices (StatServices* loader);
  virtual void loadFileServices (FileServices* loader) {}
  virtual LifeCycleEvent* clone () {return new LCE_Selection_base();}
  virtual age_t removeAgeClass   () {return NONE;}
  virtual age_t addAgeClass      () {return NONE;}
  virtual age_t requiredAgeClass () {return OFFSPRG;}
  ///@}
};


/**StatHandler class for the LCE_Selection class. Records the fitness stats.*/
class LCE_SelectionSH : public EventStatHandler< LCE_Selection_base, LCE_SelectionSH >
{
  
  vector< vector< double > > _phenoTable;
  unsigned int _table_set_gen, _table_set_age, _table_set_repl;
  
public:
  
  LCE_SelectionSH (LCE_Selection_base* event) :
  EventStatHandler< LCE_Selection_base, LCE_SelectionSH > (event),
  _table_set_gen(999999), _table_set_age(999999), _table_set_repl(999999)
  {}
  
  virtual ~LCE_SelectionSH() {}
  
  virtual bool setStatRecorders (string& token);
  
  void   addMeanPerPatch (age_t AGE);
  void   addVarPerPatch  (age_t AGE);
  void   setDataTable    (age_t AGE);
  void   setAdultTable   ( ) {setDataTable(ADULTS);}
  void   setOffsprgTable ( ) {setDataTable(OFFSPRG);}

  
  double getMeanFitness  () {return _SHLinkedEvent->_mean_fitness;}
  double getFitness  (unsigned int i) {return _SHLinkedEvent->_fitness[i];}
  double getSurvival (unsigned int i) {return _SHLinkedEvent->_survival[i];}
  double getPedProp  (unsigned int i) {return _SHLinkedEvent->_ind_cntr[i];}
  double getMeanPatchFitness (unsigned int i, unsigned int int_agex);
  double getMeanPatchFitness (unsigned int i);
  double getVarPatchFitness (unsigned int i, unsigned int int_agex);
};

#endif

