/**  $Id: individual.cc,v 1.9 2011/11/29 12:49:23 freg Exp $
*
*  @file individual.cc
*  NEMO
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
*  Created on @date 05.08.2004
*
*  @author fred
*/

#include <iostream>
#include <string.h>
#include "individual.h"
#include "ttrait.h"
#include "output.h"

unsigned long Individual::currentID = 0;

// ----------------------------------------------------------------------------------------
// individual
// ----------------------------------------------------------------------------------------
Individual::Individual ( ) 
: _age(0), _sex(MAL),_motherID(0),_fatherID(0),_mother(NULL),_father(NULL),
_home(0),_pedigreeClass(0),_fecundity(0),_trait_nb(0)
{
  _id = currentID++;
  for(unsigned int i = 0; i < 5; i++) {
    _matings[i] = 0; _realizedFecundity[i] = 0;
  }
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
Individual * Individual::init ()
{
  _age = 0;
  _pedigreeClass = 0;
  _motherID = 0;
  _fatherID = 0;
  _mother = NULL;
  _father = NULL;
  _home = 0;
  for(unsigned int i = 0; i < 5; i++) {
    _matings[i] = 0; _realizedFecundity[i] = 0;
  }
  
  if(_trait_nb != Traits.size()){
    error("Individual::init: trait counter and table size differ, resetting\n");
    _trait_nb = Traits.size();
  }
  
  for(unsigned int i = 0; i < _trait_nb; i++)
    Traits[i]->init();

  return this;
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void Individual::reset ()
{
  _id = 0;
  _age = 0;
  _sex = MAL;
  _pedigreeClass = 0;
  _motherID = 0;
  _fatherID = 0;
  _mother = NULL;
  _father = NULL;
  _home = 0;
  for(unsigned int i = 0; i < 5; i++) {
    _matings[i] = 0; _realizedFecundity[i] = 0;
  }
   
  if(_trait_nb != Traits.size()){
    warning("Individual::reset: trait counter and table size differ, resetting\n");
    _trait_nb = Traits.size();
  }
}
// ----------------------------------------------------------------------------------------
// store_data
// ----------------------------------------------------------------------------------------
void Individual::store_data ( BinaryStorageBuffer* saver )
{
  saver->store(&_id, sizeof(unsigned long));
  saver->store(&_motherID, sizeof(unsigned long));
  saver->store(&_fatherID, sizeof(unsigned long));
  saver->store(&_sex, sizeof(sex_t));
  saver->store(&_home, sizeof(unsigned int));
  saver->store(&_matings, 2*sizeof(unsigned short));
  saver->store(&_realizedFecundity, 2*sizeof(unsigned short));
  //saver->store(&_pedigreeClass, 1);
  //saver->store(&_age, sizeof(unsigned short));
}
// ----------------------------------------------------------------------------------------
// retrieve_data
// ----------------------------------------------------------------------------------------
void Individual::retrieve_data ( BinaryStorageBuffer* reader )
{
  reader->read(&_id, sizeof(unsigned long));
  reader->read(&_motherID, sizeof(unsigned long));
  reader->read(&_fatherID, sizeof(unsigned long));
  reader->read(&_sex, sizeof(sex_t));
  reader->read(&_home, sizeof(unsigned int));
  reader->read(&_matings, 2*sizeof(unsigned short));
  reader->read(&_realizedFecundity, 2*sizeof(unsigned short));
//reader->read(&_pedigreeClass, 1);
//reader->read(&_age, sizeof(unsigned short));
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void Individual::show_up ()
{
  message("\n Individual ID: %i\n\
           age: %i\n\
           sex: %i\n\
        mother: %i\n\
        father: %i\n\
pedigree class: %i\n\
          home: %i\n\
 traits values: \n",_id,_age,_sex,_motherID,_fatherID, _pedigreeClass,_home);

for(unsigned int i = 0; i < _trait_nb; i++)
  Traits[i]->show_up();
}
// ----------------------------------------------------------------------------------------
// clone
// ----------------------------------------------------------------------------------------
Individual* Individual::clone ()
{
  Individual* myClone = new Individual();
  
  for(unsigned int i = 0; i < _trait_nb; i++)
    myClone->addTrait(Traits[i]->clone(), i);
  
  return myClone;
}
// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
Individual& Individual::operator=(const Individual& i)
{
  if(this != &i) {  
    
    if(Traits.size() != i.Traits.size()) fatal("Individual::operator=:not same number of traits in left and right sides of assignment\n");
    if(_trait_nb != i._trait_nb) {
      error("Individual::operator=:trait counters differ, restting\n");
      _trait_nb = i._trait_nb;
    }
    _sex = i._sex;
    _age = i._age;
    _motherID = i._motherID;
    _fatherID = i._fatherID;
    _mother = i._mother;
    _father = i._father;
    _home = i._home;
    _pedigreeClass = i._pedigreeClass;
    _fecundity = i._fecundity;
    for(unsigned int j = 0; j < 5; j++) {
      _matings[j] = i._matings[j];
      _realizedFecundity[j] = i._realizedFecundity[j];
    }

    for(unsigned int t = 0; t < _trait_nb; t++){ 
      if(Traits[t]->get_type().compare(i.Traits[t]->get_type()) != 0) 
        fatal("Individual::operator=: not same kinds of traits on left and right sides of assignment\n");
      (*Traits[t]) = (*i.Traits[t]);
    }
  }
  return *this;
}
// ----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool Individual::operator==(const Individual& i)
{
  if(this != &i) {
    if(Traits.size() != i.Traits.size()) return false;

    for(unsigned int t = 0; t < Traits.size(); t++)
      if((*Traits[t]) != (*i.Traits[t])) return false;

  //if(_sex != i._sex) return false;
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool Individual::operator!=(const Individual& i)
{
  if(!((*this) == i))
    return true;
  else
    return false;
}
// ----------------------------------------------------------------------------------------
// getPedigreeClass
// ----------------------------------------------------------------------------------------
unsigned int Individual::getPedigreeClass (Individual* mother, Individual* father)
{
  if(mother == father) return 4; //selfed

  if(mother->getHome() != father->getHome()) return 0; //outbred between patch
  
  unsigned int mm,mf,fm,ff;
  //mother's parents:
  mm = mother->getMotherID();
  mf = mother->getFatherID();
  //father's parents:
  fm = father->getMotherID();
  ff = father->getFatherID();

  if(mm != fm && mf != ff) return 1; //outbred within patch
    
  else if((mm == fm && mf != ff) || (mm != fm && mf == ff)) return 2; //half sibs
    
  else return 3; //full sibs
}

