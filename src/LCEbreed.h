/**  $Id: LCEbreed.h,v 1.13 2015-02-05 16:11:33 fred Exp $
*
*  @file LCEbreed.h
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
*  Created on @date 07.07.2004
*
*  @author: fred
*/

#ifndef LCEBREED_H
#define LCEBREED_H
#include "lifecycleevent.h"
#include "Uniform.h"


// Class LCE_Breed_base
// 
/**Base class for the breeding (and mating) life cycle events. 
   This class registers the whole set of basic mating parameters. 
   Sets the function pointers for the mating system, fecundity and sex-ratio models,
   and breeding conditions checker (mostly checks for state of individual containers).*/

class LCE_Breed_base : public virtual LifeCycleEvent
{
  ///@name Parameters
  ///@{
  int _mating_system;
  unsigned int _mating_males;
  unsigned int _alpha_male;
  double _mating_proportion, _sd_fecundity;
  bool _do_inherit;
  ///@}
  Individual* (LCE_Breed_base::* MatingFuncPtr)   (Patch*, Individual*, unsigned int);
  Individual* (LCE_Breed_base::* DoBreedFuncPtr)  (Individual* mother, Individual* father, unsigned int LocalPatch);
  double (LCE_Breed_base::* FecundityFuncPtr)     (double mean);
  bool   (LCE_Breed_base::* CheckMatingConditionFuncPtr) (Patch* thePatch);
  sex_t  (LCE_Breed_base::* GetOffsprgSex)        ();
  
protected:
  
  TMatrix *_mean_fecundity;
  
public:

  LCE_Breed_base ();
  
  virtual ~LCE_Breed_base ( ) {if(_mean_fecundity) delete _mean_fecundity;}
  
  /**Calls the mating function according to the model chosen using the function pointer, used to get the father from the mother in a patch
    @param thePatch pointer to the focal patch where mating is taking place
    @param mother pointer to the mother, returned when mating is done by self-fertilization or cloning
    @param motherIndex index of the mother in the current patch female container, used in the \a polyginy an \a monoginy mating systems
    @return the pointer to the father following the mating scheme chosen
  */
  virtual Individual* getFatherPtr (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  { 
    return (this->*MatingFuncPtr)(thePatch, mother, motherIndex);
  }  
  ///@}
  ///@name Implementations
  ///@{
  //  virtual void init(Metapop* popPtr);
  virtual bool setParameters();
  ///@}
  
  ///@name Parameter setters/updaters
  ///@{
  bool setMatingSystem ();
  bool setFecundity ();
  bool setSexRatio ();
  ///@}
  ///@name Accessors
  ///@{
  double getMatingProportion ()                       {return _mating_proportion;}
  double getMeanFecundity    (unsigned int patch)     {return _mean_fecundity->get(0, patch);}
  int    getMatingSystem     ()                       {return _mating_system;}
  bool   doInheritance       ()                       {return _do_inherit;}
  double getPoissonFecundity (double mean)            {return RAND::Poisson(mean);}
  double getFixedFecundity   (double mean)            {return mean;}
  double getGaussianFecundity(double mean)            {
    double fec; do{fec = mean + RAND::Gaussian(_sd_fecundity);}while(fec < 0);
    return fec;}
  double getFecundity        (unsigned int patch)     {return (this->* FecundityFuncPtr)(_mean_fecundity->get(0, patch));}
  double getFecundity        (double mean)            {return (this->* FecundityFuncPtr)(mean);}
  sex_t  getOffsprgSex       ()                       {return (this->* GetOffsprgSex) ();}
  sex_t  getOffsprgSexRandom ()                       {return (sex_t)RAND::RandBool();}
  sex_t  getOffsprgSexFixed  ();
  sex_t  getOffsprgSexSelfing()                       {return FEM;}
  sex_t  getOffsprgSexCloning()                       {return FEM;}

// Kim adding:
  bool findMale ();

  ///@}
  /**Makes a new individual with the right parents.
     Calls IndFactory::makeNewIndividual. The sex of the offspring is determined by a call to 
     getOffsprgSex(). Recombination and mutation are done later, in the makeOffspring() procedure.
    @param mother pointer to the mother
    @param father pointer to the father
    @param LocalPatch index of the natal patch
    */
  Individual* breed           (Individual* mother, Individual* father, unsigned int LocalPatch);
  
  /**Makes a new individual by doing a deep copy of the mother (copies the mother's genes into the offspring).
    Calls IndFactory::getNewIndividual() and then copy the mother into the new offspring.
    @param mother pointer to the mother
    @param father pointer to the father (of no use here)
    @param LocalPatch index of the natal patch
    */
  Individual* breed_cloning   (Individual* mother, Individual* father, unsigned int LocalPatch);
  
  /**Last step of the breeding process, does inheritance and mutation of the parents' genes.
     Calls Individual::create(do_inherit, do_mutate) with do_inherit set following the local
     _do_inherit value. Updates the parent's fecundity counters.
     A breeding session looks like that:
     \code
     Individual* newind;
     Patch* natalPatch;
    
     newind = makeOffspring( do_breed( mother, father = getFatherPtr(natalPatch, mother, motherIndex), LocalPatch = natalPatch->getID() ) )
     
     natalPatch->add( newind->getSex(), OFFSx, newind );
     \endcode
    
    @param ind the offspring, as returned by the do_breed function.
    */
  Individual* makeOffspring   (Individual* ind);
  
  /**Calls the breeding function unsing its pointer.
     Used to distinguish cloning from other mating systems. 
    */
  Individual* do_breed        (Individual* mother, Individual* father, unsigned int LocalPatch)
  {
    return (this->* DoBreedFuncPtr)(mother, father, LocalPatch);
  }
  
  /**Checks if any mating will take place in the patch passed as argument.
    Is called prior to breeding in each patch. Calls the check function using its pointer.
    @param thePatch the focal patch
    */
  bool checkMatingCondition (Patch* thePatch)
  {
    return (this->* CheckMatingConditionFuncPtr) (thePatch);
  }
  
  
//****************************************************************************
/**
	I AM ADDING CODE HERE - KJG
	I think this spot should work okay because where I am adding the code in breed.cc is right next to the above function.
	...though these are all boolean, so not sure if there is anything important to do with classes there? :/
**/
  
/*  
    bool checkMatingCondition (Patch* thePatch)
  {
    return (this->* CheckMatingConditionFuncPtr) (thePatch);
  }
 */ 
  
  
//****************************************************************************

  
  /**Checks whether mating will take place in the current patch when mating is not selfing or cloning.
     Males and females must be present in the patch for mating to occur.
    @param thePatch the focal patch
    */
  bool checkNoSelfing (Patch* thePatch)
  {
    return (thePatch->size(FEM, ADLTx) != 0 && thePatch->size(MAL, ADLTx) != 0);
  }
  
  /**Checks whether mating will take place in the current patch when mating is polygynous.
    Males and females must be present in the patch for mating to occur.
    @param thePatch the focal patch
    */
  bool checkPolygyny (Patch* thePatch)
  {
    if(thePatch->size(FEM, ADLTx) == 0 || thePatch->size(MAL, ADLTx) == 0) return false;
    
//    if(thePatch->size(MAL, ADLTx) < _mating_males) _mating_males = thePatch->size(MAL, ADLTx);
    
    _alpha_male = (unsigned int)RAND::Uniform(thePatch->size(MAL, ADLTx));
    
    return true;
  }
  
  /**Checks whether mating will take place in the current patch when mating is selfing.
    Only females must be present.
    @param thePatch the focal patch
    */
  bool checkSelfing (Patch* thePatch)
  {
    if(thePatch->size(MAL, ADLTx) != 0) thePatch->flush(MAL, ADLTx, this->_popPtr);
    return (thePatch->size(FEM, ADLTx) != 0);
  }
  
  /**Checks whether mating will take place in the current patch when mating is cloning.
    Only females must be present. Males containers are flushed if not empty.
    @param thePatch the focal patch
    */
  bool checkCloning (Patch* thePatch)
  {    
    if(thePatch->size(MAL, ADLTx) != 0) thePatch->flush(MAL, ADLTx, this->_popPtr);
   
    return (thePatch->size(FEM, ADLTx) != 0);
  }
  
  
  ///@name Mating functions
  ///@{
  /**Returns a pointer to a male drawn randomly from a patch.
    @param thePatch the focal patch.
    @param mother the mother to mate with (unused here)
    @param motherIndex index of the mother in the patch adult female container (unused here)
  */
  Individual* RandomMating   (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  { return thePatch->get(MAL, ADLTx, RAND::Uniform(thePatch->size(MAL, ADLTx)) ); }
  
  /**Returns a pointer to the alpha male of the patch. The alpha male of a patch is set in the 
    LCE_Breed_base::checkPolygyny function called before mating.
    @param thePatch the focal patch.
    @param mother the mother to mate with (unused here).
    @param motherIndex index of the mother in the patch adult female container (unused here)
  */
  Individual* fullPolyginy (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  { return thePatch->get(MAL, ADLTx, _alpha_male); }
  
  /**Returns a pointer to one of the first _mating_males males of the patch.
    @param thePatch the focal patch.
    @param mother the mother to mate with (unused here).
    @param motherIndex index of the mother in the patch adult female container (unused here)
  */
  Individual* fullPolyginy_manyMales (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  { 
    if(thePatch->size(MAL,ADLTx) < _mating_males)
      return thePatch->get(MAL, ADLTx, RAND::Uniform( thePatch->size(MAL, ADLTx) ) );
    else
      return thePatch->get(MAL, ADLTx, RAND::Uniform( _mating_males ) ); 
  }
  
  /**Returns a pointer to a male from a patch chosen at random if _mating_proportion != 1, or the first male otherwise. 
    @param thePatch the focal patch.
    @param mother the mother to mate with (unused here).
    @param motherIndex index of the mother in the patch adult female container (unused here)
  */  
  Individual* partialPolyginy (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  {
    if(RAND::Uniform() > _mating_proportion)
      return RandomMating(thePatch, mother, 0);
    else
      return fullPolyginy(thePatch, 0, 0);
  }
  
  /**Returns a pointer to a male from a patch chosen at random if _mating_proportion != 1, or
    one of the _mating_males first males otherwise. 
    @param thePatch the focal patch.
    @param mother the mother to mate with (unused here).
    @param motherIndex index of the mother in the patch adult female container (unused here)
  */   
  Individual* partialPolyginy_manyMales  (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  {
    if(RAND::Uniform() > _mating_proportion)
      return RandomMating(thePatch, mother, 0);
    else
      return fullPolyginy_manyMales(thePatch, mother, 0);
  }
  
  /**Returns a pointer to a male with same index as mother (if available) from the focal patch.
     If the male is not available, one is drawn randomly from the patch.
    @param thePatch the focal patch.
    @param mother the mother to mate with (unused here)
    @param motherIndex index of the mother in the patch adult female container
    **/
  Individual* fullMonoginy (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  {
    if(thePatch->size(MAL, ADLTx) < motherIndex+1)
      return RandomMating(thePatch, mother, motherIndex);
    else
      return thePatch->get(MAL, ADLTx, motherIndex);
  }
  
  /**Returns a pointer to a male with same index as mother (if available) from the focal patch.
    If the male is not available or _mating_proportion != 1, one is drawn randomly from the patch.
    @param thePatch the focal patch.
    @param mother the mother to mate with (unused here)
    @param motherIndex index of the mother in the patch adult female container
    **/  
  Individual* partialMonoginy (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  {
    if(RAND::Uniform() > _mating_proportion || thePatch->size(MAL, ADLTx) < motherIndex+1)
      return RandomMating(thePatch, mother, motherIndex);
    else
      return thePatch->get(MAL, ADLTx, motherIndex);
  }
  
  /**Returns the mother pointer.
    @param thePatch the focal patch.
    @param mother the mother to mate with (returned here)
    @param motherIndex index of the mother in the patch adult female container (unused here)
    **/  
  Individual* fullSelfing  (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  {
      return mother;
  }
  
  /**Returns the mother pointer or a random female if _mating_proportion != 1.
    @param thePatch the focal patch.
    @param mother the mother to mate with (returned here)
    @param motherIndex index of the mother in the patch adult female container (unused here)
    **/  
  Individual* partialSelfing  (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  {
    unsigned int fem;
    if(RAND::Uniform() > _mating_proportion) {
      do {
        fem = RAND::Uniform(thePatch->size(FEM, ADLTx));
      } while(fem == motherIndex && thePatch->size(FEM, ADLTx) != 1);
      return thePatch->get(FEM, ADLTx, fem);
    }else
      return mother;
  }
  
  /**Returns a random female from the patch, will be the same mother with probability 1/N (Wright-Fisher model).
    @param thePatch the focal patch.
    @param mother the mother to mate with (unused here)
    @param motherIndex index of the mother in the patch adult female container (unused here)
    **/  
  Individual* random_hermaphrodite  (Patch* thePatch, Individual* mother, unsigned int motherIndex)
  {
    return thePatch->get(FEM, ADLTx, RAND::Uniform(thePatch->size(FEM, ADLTx)) );
  }
  ///@}
    
};

// Class LCE_Breed
// 
/**Implementation of the basic breeding and mating procedures, does not link to any trait.
   Individuals mate according to the mating system chosen. The mated adults
   are not removed from the population. The offspring containers are filled with the new generation. Note that
   they are first emptied if still containing offspring individuals when starting the breeding process.

   The population's age is set to \c ALL. The mating and realized fecundity counters of the reproducing
   males and females are updated.

   @see LCE_Breed_base
**/
class LCE_Breed : public virtual LCE_Breed_base
{
  
public:
    
  LCE_Breed ( ) : LifeCycleEvent("breed","") { }
  
  virtual ~LCE_Breed   ( ) { }
  
  ///@name Implementations
  ///@{
//  virtual void init (Metapop* popPtr);
  virtual bool setParameters ();
  virtual void  execute (); 
  
  virtual LifeCycleEvent* clone ( ) {return new LCE_Breed();}
  
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return OFFSPRG;}
  virtual age_t requiredAgeClass () {return ADULTS;}
  ///@}
};

#endif //LCEBREED_H
