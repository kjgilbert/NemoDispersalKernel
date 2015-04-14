/**  $Id: LCEbreed.cc,v 1.12 2015-02-05 16:11:33 fred Exp $
*
*  @file LCEbreed.cc
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
*  created on @date 07.07.2004
*
*  @author fred
*/
#include <deque>
#include "output.h"
#include "LCEbreed.h"
#include "individual.h"
#include "metapop.h"
#include "Uniform.h"
#include "simenv.h"

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

//                             LCE_Breed_base/

// ----------------------------------------------------------------------------------------
// LCE_Breed_base
// ----------------------------------------------------------------------------------------
LCE_Breed_base::LCE_Breed_base ()
: LifeCycleEvent("", ""),
_mating_system(0), _mating_proportion(1), _mean_fecundity(0), 
//_growthRates(0),
MatingFuncPtr(0),
DoBreedFuncPtr(0), FecundityFuncPtr(0), CheckMatingConditionFuncPtr(0), GetOffsprgSex(0)
//,GetPatchFecundityFuncPtr(0)
{
  ParamUpdater<LCE_Breed_base> * updater = new ParamUpdater<LCE_Breed_base> (&LCE_Breed_base::setMatingSystem);
  
  add_parameter("mating_system",INT,true,true,1,6, updater);
  add_parameter("mating_proportion",DBL,false,true,0,1,updater);
  add_parameter("mating_males",INT,false,false,0,0, updater);
//  add_parameter("growth_model", INT, false, true, 1, 7, updater);
//  add_parameter("growth_rate", DBL, false, false, 0, 0, updater);
  
  updater = new ParamUpdater< LCE_Breed_base > (&LCE_Breed_base::setFecundity);
  add_parameter("mean_fecundity",DBL,false,false,0,0, updater);
  add_parameter("fecundity_dist_stdev",DBL,false,false,0,0, updater);
  add_parameter("fecundity_distribution",STR,false,false,0,0, updater);
  
  updater = new ParamUpdater< LCE_Breed_base > (&LCE_Breed_base::setSexRatio);
  add_parameter("sex_ratio_mode",STR,false,false,0,0, updater);
  
  // Kim adding parameters here now
  add_parameter("breeding_connectivity_matrix",MAT,false,false,0,0,updater); // this will be the patch IDs of patches to potentially look for mates in for a given patch
  add_parameter("breeding_kernel",MAT,false,false,0,0,updater); // this will be the 1-d array holding the sorted probabilities of sending gametes to patches 1 through n, and corresponding to the IDs in the above matrix.
  add_parameter("self_if_alone", BOOL, false, false, 0, 0, updater); // ad a parameter to init file that I can set and if exists will say to self when no mates are found
  add_parameter("always_breed_window", BOOL, false, false, 0, 0, updater); // ad a parameter to init file that I can set and if exists will say to self when no mates are found
  add_parameter("never_breed_window", BOOL, false, false, 0, 0, updater); // ad a parameter to init file that I can set and if exists will say to self when no mates are found
 
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Breed_base::setParameters ()
{
  return ( setMatingSystem() && setFecundity() && setSexRatio() );
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::setMatingSystem
// ----------------------------------------------------------------------------------------
bool LCE_Breed_base::setMatingSystem ()
{
  _mating_system = (int)this->get_parameter_value("mating_system");
  
  if(get_parameter("mating_proportion")->isSet())
    _mating_proportion = this->get_parameter_value("mating_proportion");
  else
    _mating_proportion = 1;
  
  if(_paramSet->isSet("mating_males"))
    _mating_males = (int)_paramSet->getValue("mating_males");
  else
    _mating_males = 1;
  
  //set the mating functions ptr:
  CheckMatingConditionFuncPtr = &LCE_Breed_base::checkNoSelfing;
  DoBreedFuncPtr = &LCE_Breed_base::breed;
  _do_inherit = true; //is true unless the mating system is cloning
  
  switch(_mating_system) {
    //random mating:
    case 1:
    {
      MatingFuncPtr = &LCE_Breed_base::RandomMating;
      break;
    }
    //polygyny:
    case 2:
    {
      if(_mating_proportion == 1)
        if(_mating_males == 1)
          MatingFuncPtr = &LCE_Breed_base::fullPolyginy;
        else
          MatingFuncPtr = &LCE_Breed_base::fullPolyginy_manyMales;
      else
        if(_mating_males == 1)
          MatingFuncPtr = &LCE_Breed_base::partialPolyginy;
        else
          MatingFuncPtr = &LCE_Breed_base::partialPolyginy_manyMales;
        
      CheckMatingConditionFuncPtr = &LCE_Breed_base::checkPolygyny;
      break;
    }
    //monogamy:
    case 3:
    {
      if(_mating_proportion == 1)
        MatingFuncPtr = &LCE_Breed_base::fullMonoginy;
      else 
        MatingFuncPtr = &LCE_Breed_base::partialMonoginy;
      
      break;
    }
    //selfing:
    case 4:
    {
      if(_mating_proportion == 1)
        MatingFuncPtr = &LCE_Breed_base::fullSelfing;
      else
        MatingFuncPtr = &LCE_Breed_base::partialSelfing;
      
      CheckMatingConditionFuncPtr = &LCE_Breed_base::checkSelfing;
      break;
    }
    //cloning
    case 5:
    {
      if(_mating_proportion == 1)
        MatingFuncPtr = &LCE_Breed_base::fullSelfing;
      else
        MatingFuncPtr = &LCE_Breed_base::partialSelfing;
 
      CheckMatingConditionFuncPtr = &LCE_Breed_base::checkCloning;
      DoBreedFuncPtr = &LCE_Breed_base::breed_cloning;
      //      _do_inherit = false;
      break;
    }
    //random mating in Wright-Fisher model with hermaphrodites:
    case 6:
    {
      MatingFuncPtr = &LCE_Breed_base::random_hermaphrodite;
      CheckMatingConditionFuncPtr = &LCE_Breed_base::checkSelfing;
      break;
    }
      
  }
  
  //Growth model
//  unsigned int model;
//  if(_paramSet->isSet("growth_model"))
//    model = get_parameter_value("growth_model");
//  else 
//    model = 1;
//  
//  switch (model) {
//    case 1:
//      GetPatchFecundityFuncPtr = &LCE_Breed_base::instantGrowth;
//      break;
//    case 2:
//      GetPatchFecundityFuncPtr = &LCE_Breed_base::logisticGrowth;
//      break;
//    case 3:
//      GetPatchFecundityFuncPtr = &LCE_Breed_base::stochasticLogisticGrowth;
//      break;
//    case 4:
//      GetPatchFecundityFuncPtr = &LCE_Breed_base::conditionalLogisticGrowth;
//      break;
//    case 5:
//      GetPatchFecundityFuncPtr = &LCE_Breed_base::conditionalStochasticLogisticGrowth;
//      break;
//    case 6:
//      GetPatchFecundityFuncPtr = &LCE_Breed_base::fixedFecundityGrowth;
//      break;
//    case 7:
//      GetPatchFecundityFuncPtr = &LCE_Breed_base::stochasticFecundityGrowth;
//      break;
//    default:
//      GetPatchFecundityFuncPtr = &LCE_Breed_base::instantGrowth;
//      break;
//  }
//  
//  //growth rate
//  if (model > 1 && model < 6) {
//    
//    if(!_paramSet->isSet("growth_rate")) {
//      error("parameter \"growth_rate\" needs to be set\n");
//      return false;
//    }
//    
//    if(_growthRates) delete [] _growthRates;
//    _growthRates = new double [ _popPtr->getPatchNbr() ];
//    
//    if(_paramSet->isMatrix("growth_rate")) {
//      
//      TMatrix tmp;
//      
//      _paramSet->getMatrix("growth_rate", &tmp);
//      
//      if(tmp.getNbCols() != _popPtr->getPatchNbr()){
//        error("matrix argument to \"growth_rate\" has wrong number of elements,\
//              must equal the number of patches.\n");
//        return false;
//      }
//      
//      for (unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
//        _growthRates[i] = tmp.get(0, i);
//      }
//      
//      
//    } else { //not a matrix
//      
//      _growthRates[0] = get_parameter_value("growth_rate");
//      
//      for (unsigned int i = 1; i < _popPtr->getPatchNbr(); i++) {
//        _growthRates[i] = _growthRates[0];
//      }
//    }
//    
//    
//  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::setSexRatio
// ----------------------------------------------------------------------------------------
bool LCE_Breed_base::setSexRatio ()
{
  // SEX-RATIO
  
  if(get_parameter("sex_ratio_mode")->isSet()) {
    
    if(get_parameter("sex_ratio_mode")->getArg().compare("fixed") == 0) 
  
      GetOffsprgSex = &LCE_Breed_base::getOffsprgSexFixed;
    
    else if(get_parameter("sex_ratio_mode")->getArg().compare("random") != 0) {
    
      error("\"sex_ratio_mode\" parameter argument must be either \"fixed\" or \"random\".");
      return false;
    } else
      GetOffsprgSex = &LCE_Breed_base::getOffsprgSexRandom;
  } else {
    
    switch(_mating_system) {
      case 4: GetOffsprgSex = &LCE_Breed_base::getOffsprgSexSelfing;
        break;
      case 5: GetOffsprgSex = &LCE_Breed_base::getOffsprgSexCloning;
        break;
      case 6: GetOffsprgSex = &LCE_Breed_base::getOffsprgSexSelfing;
        break;
      default: GetOffsprgSex = &LCE_Breed_base::getOffsprgSexRandom;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::setFecundity
// ----------------------------------------------------------------------------------------
bool LCE_Breed_base::setFecundity ()
{
  // FECUNDITY
  if (_mean_fecundity) delete _mean_fecundity;
  
  _mean_fecundity = new TMatrix();
  
  if (get_parameter("mean_fecundity")->isMatrix()) {
    
    get_parameter("mean_fecundity")->getMatrix(_mean_fecundity);
    
    if(_mean_fecundity->getNbRows() != 1) {
      error("\"mean_fecundity\" accepts a single number or a single array with patch-specific values.\n");
      return false;
    }
    
    if(_mean_fecundity->getNbCols() > _popPtr->getPatchNbr()) {
      error("\"mean_fecundity\" accepts an array of max num patches in length.\n");
      return false;
    }
    else if (_mean_fecundity->getNbCols() < _popPtr->getPatchNbr()) 
    {
      unsigned int npat = _mean_fecundity->getNbCols();
      
      TMatrix tmp(*_mean_fecundity);
      
      _mean_fecundity->reset(1, _popPtr->getPatchNbr());
      
      for (unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
        _mean_fecundity->set(0, i, tmp.get(0 , i % npat));
      }
    }
  } else {
    _mean_fecundity->reset(1, _popPtr->getPatchNbr()); //single array
    _mean_fecundity->assign(get_parameter_value("mean_fecundity"));//will be set to -1 if param isn't set    
  }

  

  if(get_parameter("fecundity_distribution")->isSet()) {
    
    string dist = get_parameter("fecundity_distribution")->getArg();
    
    if( dist.compare("fixed") == 0 )
    
      FecundityFuncPtr = &LCE_Breed_base::getFixedFecundity;
    
    else if( dist.compare("poisson") == 0 )
      
      FecundityFuncPtr = &LCE_Breed_base::getPoissonFecundity;
    
    else if( dist.compare("normal") == 0 ) {
      
      FecundityFuncPtr = &LCE_Breed_base::getGaussianFecundity;
      
      if(get_parameter("fecundity_dist_stdev")->isSet())
        _sd_fecundity = get_parameter_value("fecundity_dist_stdev");
      else {
        error("parameter \"fecundity_dist_stdev\" is missing!\n");
        return false;
      }
      
    } else {
      error("unknown fecundity distribution parameter's argument!\n");
      return false;
    }
  
  } else { //default distribution is Poisson:
    FecundityFuncPtr = &LCE_Breed_base::getPoissonFecundity;
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::getOffsprgSexFixed
// ----------------------------------------------------------------------------------------
sex_t LCE_Breed_base::getOffsprgSexFixed  ()
{
  static bool sex = RAND::RandBool();
  sex ^= 1;
  return (sex_t)sex;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed_base::makeOffspring
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_base::makeOffspring(Individual* ind)
{
  unsigned int cat = ind->getPedigreeClass();
  
  ind->getMother()->DidHaveABaby(cat);
  if(cat!=4) ind->getFather()->DidHaveABaby(cat);

  return ind->create(doInheritance(), true);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::findMales      NOT WORKING, WHY? OR IS THIS THE RIGHT WAY TO TRY?
// ----------------------------------------------------------------------------------------
/*  bool LCE_Breed_base::findMale ()
{
  bool malePresent = 0;   // if no one in the focal patch, find a nearby male to mate with
  Patch* focalPatch;
   
  focalPatch = _popPtr->getPatch(i);
    
  aimedPatchMatrix[0] = new TMatrix();       // make an empty matrix
  _paramSet->getMatrix("breeding_connectivity_matrix",aimedPatchMatrix[0]);  // fill it with parameters from input file
  aimedPatchMatrix[1] = new TMatrix(*aimedPatchMatrix[0]); // don't know what this is doing? copied from LCEdisperse.

  while(malePresent = 0)     // do I need ";" here?
  do{          // find a patch in the breeding kernel that contains a male, as long as find at least one, exit loop and go to next step
    
    aimedList = aimedPatchMatrix[i]; // pull out the row of the matrix for the focal patch?
    
    for(unsigned int j = 1; j < size(aimedList), j++ ) { // j = 1 because have already checked focal patch in previous line and only continuing if no male there(spot one is j=0)

      fatherPatch = _popPtr->getPatch(j); // check the next patch in the list
        
      if( checkMatingCondition(fatherPatch) ) { // if there IS a male in the patch being checked (checkMatCond is boolean)
             //first patch you find with >0 males changes male_present to TRUE
        malePresent = 1;   // then change to true, and loop should stop searching, and we proceed onward to breeding 
        break;  // BREAK OUT OF THE FOR LOOP IF FIND ANY 1 MALE   -- maybe don't need it because of the while-do? 
      }
      
    } // end for loop, if enter the if statement, should leave for loop early with malePresent = 1;
      // if finds no males, should still have "male_present = 0" here
  }   // end do loop
  
  if(!malePresent) {
    return false;      // if found no males, return false
  }
  
  return true;
}  */
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::breed
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_base::breed(Individual* mother, Individual* father, unsigned int LocalPatch)
{
  return _popPtr->makeNewIndividual(mother, father, getOffsprgSex(), LocalPatch);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::breed_cloning
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_base::breed_cloning(Individual* mother, Individual* father, unsigned int LocalPatch)
{
  Individual *newind;
  
  if(mother == father) {
    newind = _popPtr->getNewIndividual();
    //cloning:
    (*newind) = (*mother); 
    
    newind->reset_counters();
    newind->setFather(NULL);
    newind->setFatherID(0);
    newind->setMother(mother);
    newind->setMotherID(mother->getID());
    newind->setIsSelfed(true);
    newind->setHome(LocalPatch);
    newind->setAge(0);
    _do_inherit = false;  
  } else {
    newind =  _popPtr->makeNewIndividual(mother, father, getOffsprgSex(), LocalPatch);
    _do_inherit = true;
  }
  
  return newind;
}
// ----------------------------------------------------------------------------------------

//                             LCE_Breed

// ----------------------------------------------------------------------------------------
// LCE_Breed::init
// ----------------------------------------------------------------------------------------
bool LCE_Breed::setParameters ( )
{ 
  return LCE_Breed_base::setParameters( ); 
}
// ----------------------------------------------------------------------------------------
// LCE_Breed::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed::execute()
{
  Patch* patch;
  Individual* mother;
  Individual* father;  // make a father variable of class "individual"
  Individual* NewOffsprg;
  unsigned int nbBaby;
  // Kim adding following params
  bool malePresent; // don't necessarily need this here if I create it with the declaration at the same time below
  bool breedWindow;
  Patch* fatherPatch;
  Patch* checkPatch;
  Patch* focalPatch;
  bool doSelfing;
  
#ifdef _DEBUG_
  message("LCE_Breed::execute (Patch nb: %i offsprg nb: %i adlt nb: %i)\n"
          ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS ));
#endif
  
  if(_popPtr->size(OFFSPRG) != 0) {
    warning("offspring containers not empty at time of breeding, flushing.\n");
    _popPtr->flush(OFFSx);
  }
  //because mean fecundity can be patch-specific, we have to check whether the patch number changed
  if (_mean_fecundity->getNbCols() != _popPtr->getPatchNbr()) LCE_Breed_base::setFecundity(); // find the mother's patch's fecundity


  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) { // which patch are we in
    
    patch = _popPtr->getPatch(i);  // current patch we're doing breeding for

// ARE THESE IN THE RIGHT SPOT? might be faster if I move them out of the for loop, but if I do that, after one loop will they always stay at that value?
    doSelfing = 0;   // set defaults
    breedWindow = 0;  
    if(_paramSet->isSet("self_if_alone")) doSelfing = 1;  // only change defaults if specified from init file
    if(_paramSet->isSet("always_breed_window")) breedWindow = 1;   // 

    if( (_paramSet->isSet("never_breed_window")) && !checkMatingCondition(patch) && !doSelfing) continue;
    	// if no males and we never want to use the breeding window and we're not allowing selfing, nemo behaves like the old version and exits loop and continue through code
		// this calls "checkNoSelfing" function in line 93 of this file which is defined in LCEbreed.h
		//	checkNoSelfing returns true if it counts >0 females and >0 males in a patch 	
/******  MY EDITED VERSION
// below should all be uncommented once working

    if( !checkMatingCondition(patch) ) {    // if breedWindow is true (=1) then will always use the breeding window regardless of any males in focal patch
		breedWindow = 1;           // to use below when I tell it to actually use the breeding window function, 1=true

		bool malePresent = 0;   // if no one in the focal patch, find a nearby male to mate with; false = 0, true = 1
		Patch* focalPatch;

		focalPatch = _popPtr->getPatch(i);

		// maybe not right? no female dispersal here, just male gametes
*		aimedPatchMatrix[0] = new TMatrix();       // make an empty matrix
*		_paramSet->getMatrix("breeding_connectivity_matrix",aimedPatchMatrix[0]);  // fill it with parameters from input file
*		aimedPatchMatrix[1] = new TMatrix(*aimedPatchMatrix[0]); // don't know what this is doing? copied from LCEdisperse.

		// find a patch in the breeding kernel that contains a male, as long as find at least one, exit loop and go to next step

*		aimedList = aimedPatchMatrix[i]; // pull out the row of the matrix for the focal patch?

		//CHECK - UNSAFE because might be approximate and not have given the number of elements -- hopefully elsewhere should have this infor when the matrix is created    
*		unsigned int lengthAimedList = sizeof(aimedList) /  sizeof(aimedList[0]);  // b/c all elements have the same size, this is the way to find the length, i.e. number of elements in the array

		for(unsigned int j = 0; j < lengthAimedList; j++ ) { // j = 1 because have already checked focal patch in previous line and only continuing if no male there(spot one is j=0)

			fatherPatch = _popPtr->getPatch(j); // check the next patch in the list

			if( checkMatingCondition(fatherPatch) ) { // if there IS a male in the patch being checked (checkMatCond is boolean)
					  //first patch you find with >0 males changes male_present to TRUE
				 malePresent = 1;   // then change to true, and loop should stop searching, and we proceed onward to breeding 
				 break;  // BREAK OUT OF THE FOR LOOP IF FIND ANY 1 MALE 
		     }

		} // end for loop, if enter the if statement, should leave for loop early with malePresent = 1;
		// if finds no males, should still have "male_present = 0" here
        if( !malePresent && !doSelfing ) continue;
     }
 
*/
    
    unsigned int cnt =0;
    for(unsigned int size = patch->size(FEM, ADLTx), indexOfMother = 0;
        indexOfMother < size;
        indexOfMother++)
    {
      mother = patch->get(FEM, ADLTx, indexOfMother);
      
      nbBaby = (unsigned int)mother->setFecundity( getFecundity(i) ) ; //allows for patch-specific fec
      cnt += nbBaby;
      //-----------------------------------------------------------------------
      while(nbBaby != 0) {
        
      // if we have reached this point, then we know there will be at least one male in at least one patch within the breeding kernel or that selfing will occur
      // want to implement "backwards migration" here to then correctly and randomly find the male that becomes the father
         // will have to on the fly convert the forward migration matrix input from the init file into the appropriate backwards migration matrix
         
        // use the corrected backwards migration matrix to find the correct father patch, then just continue to next function as it exists? 
/*


        if(breedWindow){ // then find the other patch that the father comes from
             //  males per patch in breeding window: 
           
  ** sizeof ....         unsigned int lengthBreedKernel = sizeof(breeding_kernel) /  sizeof(breeding_kernel[0]);  // b/c all elements have the same size, this is the way to find the length, i.e. number of elements in the array
       // NEW things here    
 *          unsigned int arrayNumMales[lengthBreedKernel]; // empty array to fill in number of males per patch
 *          double numerator[lengthBreedKernel];
 *          double normalBreedKernel[lengthBreedKernel];
*		aimedList = aimedPatchMatrix[i]; // pull out the row of IDS (NOT PROB) the matrix for the focal patch
  // normalize mating probabilities into backwards migration rates
           double denominator = 0;


           for(unsigned int k = 0; k < lengthBreedKernel; k++){ 

               checkPatch = _popPtr->getPatch(aimedList[k]); // check the patch being iterated - this should be the patch's ID number relative to the whole landscape

               arrayNumMales[k] = checkPatch->size(MAL, ADLTx);  // put that number in the respective spot in the new array
                              
               numerator[k] = (checkPatch->size(MAL, ADLTx))*(breeding_kernel[aimedList[k]]);
               
               denominator += numerator[k];
               
               
           }   // end for loop finding number of males per aimed patch
           
                                 
          // have to iterate through to divide an array by a single number
           for(unsigned int k = 0; k < lengthBreedKernel; k++){ 
            
             normalBreedKernel[k] = numerator[k] / denominator;  // numerator is an array, denominator is a number
           
           }

           
             // draw a random number between 0 and 1, see what patch that picks probability-wise
             // some spots in the array will be zero, so should never be picked because there are no males there
           
           double sumProbs = 0;
           double cumSums[lengthBreedKernel];
           double total = 0;
           
           for(unsigned int k = 0; k < lengthBreedKernel; k++) {
              
              total += breeding_kernel[k];
              cumSums[k] = total;
           
           }
           
           unsigned int c = 0;
           unsigned int fatherPatchID;
           double randNum = RAND::Uniform();  // maybe check that this isn't the default rand num generator, but is instead one made for nemo 
           
           while(c < lengthBreedKernel) {
           
              if(randNum < cumSums[c]) {
                 fatherPatchID = aimedList[c]; // to get the universal patch ID
                 break; // this breaks out of the whole while loop, c won't iterate up
              } 
              c++;
           }
           
           assert(c < lengthBreedKernel); // if somehow the function above doesn't work, this means it didn't find a box with the probability matching the rand number and c became greater than length of disp kernel      
               // if this happens, probably when one patch has super high prob vs others, and in that case cn add code to fix because that patch is probably the one to choose from
 
          fatherPatch = _popPtr->getPatch(fatherPatchID);
        
          father = this->getFatherPtr(fatherPatch, mother, indexOfMother);

// end of if, DELETE THINGS HERE          
          
        } else { // get the father from the focal patch - will not reach here if no selfing and no male present in focal patch
        
          // breeding window is not happening, = 0, one of 2 things happens, normal nemo with a male in the patch, or if no male then we self

            // focal patch has male, proceed with normal nemo:
          if( checkMatingCondition(patch) ) father = this->getFatherPtr(patch, mother, indexOfMother); // if there was a focal patch male, normal nemo mating within patch
          else if(doSelfing) father = mother;

        }  
        
*/
// comment next line out once I get my code working above
        father = this->getFatherPtr(patch, mother, indexOfMother); // probably will want to change this to a new function using my inputs
        
        NewOffsprg = makeOffspring( do_breed(mother, father, i) );
        
        patch->add(NewOffsprg->getSex(), OFFSx, NewOffsprg); // this is okay, b/c offspring is always added to FOCAL patch, not source patch
        
        nbBaby--;
      }//_END__WHILE
      
    } // end for loop going through all mothers in given patch, size = total number of mothers, i.e. number of females in the patch
  } // end for loop going through all patches, i
}
