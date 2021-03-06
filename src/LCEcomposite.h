/**  $Id: LCEcomposite.h,v 1.10 2015-04-01 14:15:04 fred Exp $
*
*  @file LCEcomposite.h
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

#ifndef LCE_COMPOSITE_H
#define LCE_COMPOSITE_H

#include <cmath>
#include "LCEbreed.h"
#include "LCEdisperse.h"
#include "LCEselection.h"
#include "lifecycleevent.h"

//CLASS LCE_Breed_Disperse
//
/**Performs breeding and migration in one, migration rates are backward rates.
   Inherits parameters from LCE_Breed_base and LCE_Disperse_base. Population size is constant,
   There is no demographic stochasticity with this LCE. The number of colonizers of extinct patches
   can be set differently from the patch carrying capacities.*/
class LCE_Breed_Disperse : public virtual LCE_Breed_base, public virtual LCE_Disperse_base
{
  void (LCE_Breed_Disperse::* _breed_disperse)();
  
  /**Maximum size of a patch after colonisation.*/
  int _num_colonizers;
  
  /**Which part of the population is dispersing its gametes.*/
  sex_t _dispersing_sex;

  /**Patch-specific growth rates.*/
  double *_growthRates;
    
protected: //so that the derived class can use these function pointers
  Individual*  (LCE_Breed_Disperse::* _make_offspring)(sex_t SEX, Patch* patch, unsigned int LocalPatch);
  unsigned int (LCE_Breed_Disperse::* _get_numFemOffspring)(Patch* patch);
  unsigned int (LCE_Breed_Disperse::* _get_numMalOffspring)(Patch* patch);
  unsigned int (LCE_Breed_Disperse::* _get_patchFecundity )(Patch* patch, sex_t SEX);
  
public:
  
  LCE_Breed_Disperse();
  
  virtual ~LCE_Breed_Disperse (){if(_growthRates) delete [] _growthRates;}
  
  void do_breed_disperse             ();
  void do_breed_disperse_propagule   ();
  
  
  unsigned int numFemOffspring              (Patch *patch)
  {
    if (patch->size(FEM, ADLTx) == 0 && _dispersing_sex != FEM)//to avoid choosing a local female when there is none
      return 0;
    else
      return (this->*_get_patchFecundity)(patch, FEM);
  }
  
  /**The number of females produced in case a max number of colonizers was specified.*/
  unsigned int numFemOffspring_colonizers   (Patch *patch)
  {
    if(patch->get_isExtinct()) {
      //check if only one sex sends gametes around, capturing the case of pollen dispersal
      if(_dispersing_sex != FEM)
        return 0;
      else
        return _num_colonizers; 
    } else
      return (this->*_get_patchFecundity)(patch, FEM);
  }

  /**The number of males produced is always zero when the mating system is not random mating.*/
  unsigned int numMalOffspring_notrandom     (Patch *patch)
  {
    return 0;
  }
  
  unsigned int numMalOffspring_random        (Patch *patch)
  { //return 0 if only males gametes migrate and there is no females in the patch
    //to avoid choosing a local female when there is none
    //(we assume it doesn't make sense to have only female gametes migrating)
    if (patch->size(FEM, ADLTx) == 0 && _dispersing_sex != FEM)
      return 0;
    else
      return (this->*_get_patchFecundity)(patch, MAL);
  }
  
  unsigned int numMalOffspring_random_colonizers (Patch *patch)
  {
    if(patch->get_isExtinct()) {
      //check if only one sex sends gametes around, capturing the case of pollen dispersal
      if(_dispersing_sex != FEM)
        return 0;
      else
        return _num_colonizers; 
    } else
      return (this->*_get_patchFecundity)(patch, MAL);
  }
  
  ///@name Growth Functions
  ///@{
  /**The number of offspring produced corresponds to the carrying capacity of the patch.*/
  unsigned int instantGrowth  (Patch* patch, sex_t SEX)
  {
    return patch->get_K(SEX);
  }
  /**The number of offspring produced is given by the logistic growth function.*/
  unsigned int logisticGrowth (Patch* patch, sex_t SEX) 
  {
    double K = (double)patch->get_K(SEX);
    double r = _growthRates[ patch->getID() ];
    double N = (double)patch->size(ADLTx);
    return (unsigned int)ceil(N + r*N*((K-N)/K));
  }
  /**The number of offspring produced is drawn from a Poisson with mean equal to the logistic growth predicate.*/
  unsigned int stochasticLogisticGrowth (Patch* patch, sex_t SEX) 
  {
    return (unsigned int)RAND::Poisson((double)logisticGrowth(patch, SEX));
  }
  /**The number of offspring produced depends on the adult density. A minimum number of offspring is
     produced when less than half of the carrying capacity of adults is present. It is logistic otherwise.*/
  unsigned int conditionalLogisticGrowth (Patch* patch, sex_t SEX)
  {
    if (patch->size(SEX, ADLTx) < patch->get_K(SEX)/2) {
      return fixedFecundityGrowth(patch, SEX);
    } else {
      return logisticGrowth(patch, SEX);
    }
  }
  /**The number of offspring produced depends on the adult density, similar to 'conditionalLogisticGrowth'
      except that this time, the fecundities are drawn from Poisson distributions.*/
  unsigned int conditionalStochasticLogisticGrowth (Patch* patch, sex_t SEX)
  {
    if (patch->size(SEX, ADLTx) < patch->get_K(SEX)/2) {
      return stochasticFecundityGrowth(patch, SEX);
    } else {
      return stochasticLogisticGrowth(patch, SEX);
    }
  }
  /**The number of offspring produced is equal to the carrying capacity multiplied by the mean fecundity
     of the focal patch.*/
  unsigned int fixedFecundityGrowth (Patch* patch, sex_t SEX)
  {
    return patch->get_K(SEX)*LCE_Breed_base::getMeanFecundity(patch->getID());
  }
  /**The number of offspring produced is a random number drawn from a Poisson distribution with mean
   equal to the carrying capacity multiplied by the mean fecundity of the focal patch.*/
  unsigned int stochasticFecundityGrowth (Patch* patch, sex_t SEX)
  {
    return RAND::Uniform(patch->get_K(SEX)*LCE_Breed_base::getMeanFecundity(patch->getID()));
  }
  ///@}
  
  ///@name Mating Functions
  ///@{
  Individual* mate_random         (sex_t SEX, Patch *patch, unsigned int LocalPatch);
  Individual* mate_random_hermaphrodite (sex_t SEX, Patch *patch, unsigned int LocalPatch);
  Individual* mate_selfing        (sex_t SEX, Patch *patch, unsigned int LocalPatch);
  Individual* mate_full_selfing   (sex_t SEX, Patch *patch, unsigned int LocalPatch);
  Individual* mate_cloning        (sex_t SEX, Patch *patch, unsigned int LocalPatch);
  Individual* makeOffspring       (Individual* ind);
  void breed_disperse(sex_t SEX, Patch *patch, unsigned int LocalPatch, unsigned int size);
  
  Individual* get_parent          (sex_t SEX, sex_t DispSex, Patch* LocalPatch, unsigned int patchNbr);
  ///@}
  
  ///@name Implementations
  ///@{
//  virtual void init(Metapop* popPtr);  
  virtual bool setParameters ();
  virtual void execute ();
  virtual LifeCycleEvent* clone () {return new LCE_Breed_Disperse();}
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return NONE;}
  virtual age_t addAgeClass ( ) {return OFFSPRG;}
  virtual age_t requiredAgeClass () {return ADULTS;}
  ///@}
  
};

// Class LCE_Breed_Selection
// 
/**Composite LCE implementing breeding and viability selection on a given trait type.
Inherits from LCE_Breed_base and LCE_Selection_base.*/
class LCE_Breed_Selection : public virtual LCE_Breed_base, public virtual LCE_Selection_base
{
  
  void (LCE_Breed_Selection::* _breed_selection) (Patch* patch, unsigned int patchID, unsigned int *cntr);
  
  unsigned int _nb_trait;
  
  vector< unsigned int > _nonSelectedTraitIndices;
  
public:
    
  LCE_Breed_Selection ( ); // : LifeCycleEvent("breed_selection","") {}
  virtual ~LCE_Breed_Selection ( ) {}
  
  /**Performs viability selection and breeding at the same time. The selected trait is
     first inherited (with recombination and mutations) and survival checked before the 
     remaining traits are inherited. This helps save some computing time. The function
     returns a null pointer in case of selective death.
    @param ind the newborn, its new traits have not been computed yet
    @param natalpatch the index of the natal patch of the offfspring where selectin takes palce*/
  Individual* makeOffspringWithSelection (Individual* ind, unsigned int natalpatch);

  void do_breed_selection_FecFitness (Patch* patch, unsigned int patchID, unsigned int *cntr);
  void do_breed_selection_OffSurvival (Patch* patch, unsigned int patchID, unsigned int *cntr);
  
  ///@name Implementations
  ///@{
//  virtual void init (Metapop* popPtr);
  virtual bool setParameters ();
  virtual void  execute ();
  virtual LifeCycleEvent* clone ( ){ return new LCE_Breed_Selection(); }
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader );
  virtual age_t removeAgeClass ( ) {return NONE;}
  virtual age_t addAgeClass ( ) {return OFFSPRG;}
  virtual age_t requiredAgeClass () {return ADULTS;}
  ///@}  
};

//CLASS LCE_Breed_Selection_Disperse
//
/**Composite LCE performing breeding, migration and viability selection all in one. Migration rates are 
   backward rates and population sizes are constant.
   Inherits from LCE_Breed_Disperse and LCE_Breed_Selection (and thus from LCE_Breed_base, LCE_Selection_base,
   and LCE_Disperse_base as well).
   Fitness is always absolute here.
*/
class LCE_Breed_Selection_Disperse : public virtual LCE_Breed_Disperse, public virtual LCE_Breed_Selection
{
//  int _max_try;
  
  void (LCE_Breed_Selection_Disperse::* _breed_selection_disperse)(int* counter);
  
public:

  LCE_Breed_Selection_Disperse();
  
  virtual ~LCE_Breed_Selection_Disperse() {}
  
  void breed_selection_disperse(int* counter);
  void breed_selection_disperse_propagule(int* counter);
  
  void do_breed (sex_t SEX, unsigned int size, int* counter, Patch* patch, unsigned int patchNbr);
    
  ///@name Implementations
  ///@{
//  virtual void init(Metapop* popPtr);  
  virtual bool setParameters ();
  virtual void execute ();
  virtual LifeCycleEvent* clone () {return new LCE_Breed_Selection_Disperse();}
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader );
  virtual age_t removeAgeClass ( ) {return NONE;}
  virtual age_t addAgeClass ( ) {return OFFSPRG;}
  virtual age_t requiredAgeClass () {return ADULTS;}
  ///@}
};
#endif
