/**  $Id: ttdeletmutations_bitstring.cc,v 1.23 2015-03-20 08:00:58 fred Exp $
 *
 *  @file ttdeletmutations_bitstring.cc
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
 *  Created on @date 05.08.2004
 *  @author fred
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <list>
#include "ttdeletmutations_bitstring.h"
#include "Uniform.h"
#include "output.h"

float** TTDeletMutations_bitstring::_effects = NULL;
void    TTDeletMutations_bitstring::set_effects (float** fx) {_effects = fx;}

// ------------------------------------------------------------------------------

//                          TProtoDeletMutations/

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TProtoDeletMutations_bitstring::TProtoDeletMutations_bitstring () 
: _nb_locus(0), _fitness_model(0), _mutation_model(0), _dominance_model(1), 
_fitness_scaling_factor(1), _init_freq(0),_mut_rate(0), _strength(0), 
_dominance(0), _continuous_effects(0), _viability_func_ptr(0), _inherit_func_ptr(0), _stats(0),
_writer(0), _reader(0), _effects(0)
{
  set_paramset("delet", false, this);
  
  add_parameter("delet_loci",INT,true,false,0,0);
  
  add_parameter("delet_mutation_rate",DBL,true,true,0,1, 0);
  add_parameter("delet_mutation_model",INT,false,true,1,2, 0);
  add_parameter("delet_init_freq",DBL,false,true,0,1, 0);
  
  //genetic map parameters:
  TTProtoWithMap::addGeneticMapParameters("delet");
  
  add_parameter("delet_fitness_model",INT,true,true,1,2, 0);
  add_parameter("delet_effects_distribution",STR,false,false,0,0, 0);
  add_parameter("delet_effects_mean",DBL,true,false,0,1, 0);
  add_parameter("delet_effects_dist_param1",DBL,false,false,0,0, 0);
  add_parameter("delet_effects_dist_param2",DBL,false,false,0,0, 0);
  add_parameter("delet_dominance_mean",DBL,true,true,0,1, 0);
  add_parameter("delet_fitness_scaling_factor",DBL,false,false,0,0, 0);
  add_parameter("delet_dominance_model", INT, false, true, 1, 2, 0);
  //for back compatibility:
  add_parameter("delet_sel_coef",DBL,false,true,0,1, 0);
  add_parameter("delet_dom_coef",DBL,false,true,0,1, 0);
  //output parameters:
  add_parameter("delet_save_genotype",BOOL,false,false,0,0);
  add_parameter("delet_genot_dir",STR,false,false,0,0);
  add_parameter("delet_genot_logtime",INT,false,false,0,0);
}
// ----------------------------------------------------------------------------------------
// copy cstor
// ----------------------------------------------------------------------------------------
TProtoDeletMutations_bitstring::TProtoDeletMutations_bitstring (const TProtoDeletMutations_bitstring& T) 
: _nb_locus(T._nb_locus), _fitness_model(T._fitness_model), _mutation_model(T._mutation_model),
_dominance_model(T._dominance_model),_fitness_scaling_factor(T._fitness_scaling_factor), 
_init_freq(T._init_freq), _mut_rate(T._mut_rate), _strength(T._strength), _dominance(T._dominance),
_continuous_effects(T._continuous_effects), _viability_func_ptr(T._viability_func_ptr), 
_inherit_func_ptr(T._inherit_func_ptr), _stats(0), _writer(0), _reader(0), _effects(0)
{
  _paramSet = new ParamSet( *(T._paramSet) ) ;
}
// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TProtoDeletMutations_bitstring::~TProtoDeletMutations_bitstring ()
{
  if(_stats != NULL)      {delete _stats;     _stats = NULL;}
  if(_reader)             {delete _reader;    _reader= NULL;}
  if(_writer != NULL)     {delete _writer;   _writer = NULL;} 
  if(_effects != NULL)    {
    delete [] _effects[0];
    delete [] _effects[1];
    delete [] _effects;
    _effects = NULL;
  }
  
}
// ----------------------------------------------------------------------------------------
// setParameters
// ----------------------------------------------------------------------------------------
bool TProtoDeletMutations_bitstring::setParameters ()
{  
  _nb_locus = (int)get_parameter_value("delet_loci");
  _mut_rate = get_parameter_value("delet_mutation_rate");
  _fitness_model = (int)get_parameter_value("delet_fitness_model");
  
  _inherit_func_ptr = &TProtoDeletMutations_bitstring::inherit_low;
  
  //mutation model
  if( get_parameter("delet_mutation_model")->isSet() )
    _mutation_model = (unsigned int)get_parameter_value("delet_mutation_model");
  else
    _mutation_model = 1;
  
  if( (_init_freq = get_parameter_value("delet_init_freq")) == -1.0)
    _init_freq = 0;
  
  //selection parameters:
  if(!setSelectionParameters()) return false;
  
  return TTProtoWithMap::setGeneticMapParameters ("delet");
}
// ----------------------------------------------------------------------------------------
// setSelectionParameters
// ----------------------------------------------------------------------------------------
bool TProtoDeletMutations_bitstring::setSelectionParameters ()
{
  //effects and dominance parameters:
  _strength = get_parameter_value("delet_effects_mean");
  _dominance = get_parameter_value("delet_dominance_mean");
  _dist_p1 = get_parameter_value("delet_effects_dist_param1");
  _dist_p2 = get_parameter_value("delet_effects_dist_param2");
  
  _dominance_model = (int)get_parameter_value("delet_dominance_model");
  if(_dominance_model == -1) _dominance_model = 1; //1 means continuous; 2 means constant
  
  if(get_parameter("delet_sel_coef")->isSet())
    _strength = get_parameter_value("delet_sel_coef");
  if(get_parameter("delet_dom_coef")->isSet())
    _dominance = get_parameter_value("delet_dom_coef");
  
  string distro = _paramSet->getArg("delet_effects_distribution");
  
  _continuous_effects = true;
  
  if(distro.empty() || distro.compare("constant") == 0) {
    
    _continuous_effects = false;
    
    _set_effects_func = 0;
    
  } else if(distro.compare("exponential") == 0) {
    
    _set_effects_func = &TProtoDeletMutations_bitstring::set_effects_exp;
    
  } else if(distro.compare("gamma") == 0) {
    
    _set_effects_func = &TProtoDeletMutations_bitstring::set_effects_gamma;
    //check params:
    if(_dist_p1 == -1) {
      error("missing parameter 1 (shape) of the gamma distribution for the deleterious mutation effects.");
      return false;
    }
    
    if(_dist_p2 == -1)
      _dist_p2 = _strength / _dist_p1; //mean of gamma == _strength = a*b with a == shape == _dist_p1
    
  } else if(distro.compare("lognormal") == 0) {
    
    _set_effects_func = &TProtoDeletMutations_bitstring::set_effects_lognorm;
    //check params:
    if(_dist_p1 == -1) {
      error("missing parameter 1 (mu) of the log-normal distribution for the deleterious mutation effects.");
      return false;
    }
    if(_dist_p2 == -1) {
      error("missing parameter 2 (sigma) of the log-normal distribution for the deleterious mutation effects.");
      return false;
    }
    
  } else {
    error("Deleterious effects distribution \"%s\" not implemented.",distro.c_str());
    return false;
  }
  
  if( (_fitness_scaling_factor = get_parameter_value("delet_fitness_scaling_factor")) == -1.0)
    _fitness_scaling_factor = 1.0;
  
  if(_continuous_effects) set_effects();
  
  return true;
}
// ----------------------------------------------------------------------------------------
// set_effects
// ----------------------------------------------------------------------------------------
void TProtoDeletMutations_bitstring::set_effects()
{ 
  if(_effects != NULL) {
    delete [] _effects[0];
    delete [] _effects[1];
    delete [] _effects;
  }
  
  _effects = new float* [2];
  _effects[0] = new float [_nb_locus];
  _effects[1] = new float [_nb_locus];
  
  double dom;
  double k = -1.0*log(2*_dominance) / _strength;
  
  for(unsigned int i = 0; i < _nb_locus; ++i){
    
    if(RAND::Uniform()>0.03) {   // remove "rand()" by Remi, better to use built-in RAND::Uniform()
      _effects[1][i] = 1;
      _effects[0][i] = 0;
    } else {
      do{
        _effects[1][i] = (float)(this->*_set_effects_func)(); //homozygote effect: s  ORIG NEMO
      }while(_effects[1][i] > 1); //truncate distribution		ORIG NEMO
      
      if(_dominance_model == 1)	// ORIG NEMO
        dom = exp(-1.0*_effects[1][i]*k)/2.0; //scaling of h on s    ORIG NEMO
      else dom = _dominance;	// ORIG NEMO
      
      _effects[0][i] = (float)(dom * _effects[1][i]); //heterozygote effect: hs
    }  // end new code added - this was changing the distribution of adding lethal mutations to if greater than 0.3, adds a lethal,  REMEMBER that at each del. locus there is a set size mutation for the whole course of the simulation
  }
  //set the TTDeletMutations global var:
  TTDeletMutations_bitstring::set_effects(_effects);
}
// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TProtoDeletMutations_bitstring::loadStatServices ( StatServices* loader ) 
{
  if(_stats != NULL) 
    delete _stats;
  
  _stats = new TTDeletMutBitstrSH(this);
  
  loader->attach(_stats);
}// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void TProtoDeletMutations_bitstring::loadFileServices  (FileServices* loader)
{ 
  // --- THE READER ---
  //always add the reader:
  if(_reader) delete _reader;
  _reader = new TTDeletMutBitstrFH(this);
  //set to read mode:
  _reader->set_isInputHandler(true);
  //attach to file manager:
  loader->attach_reader(_reader);
  
  //writer
  if(get_parameter("delet_save_genotype")->isSet()) {
    
    if(_writer == NULL) _writer = new TTDeletMutBitstrFH(this);
    
    Param* param = get_parameter("delet_genot_logtime");
    
    if(param->isMatrix()) {
      
      TMatrix temp;
      param->getMatrix(&temp);
      _writer->set_multi(true, true, 1, &temp, get_parameter("delet_genot_dir")->getArg());
      
    } else {      
      //           rpl_per, gen_per, rpl_occ, gen_occ, rank, path, self-ref
      _writer->set(true, true, 1, (param->isSet() ? (int)param->getValue() : 0),
                   0, get_parameter("delet_genot_dir")->getArg(),this);
    }
    
    loader->attach(_writer);
    
  } else if(_writer != NULL) {
    delete _writer;
    _writer = NULL;
  }
}
// ----------------------------------------------------------------------------------------
// hatch
// ----------------------------------------------------------------------------------------
TTDeletMutations_bitstring* TProtoDeletMutations_bitstring::hatch ( )
{
  TTDeletMutations_bitstring* new_trait = new TTDeletMutations_bitstring();
  
  new_trait->set_proto(this);
  new_trait->set_nb_locus(_nb_locus);
  new_trait->set_init_freq(_init_freq);
  new_trait->set_fitness_model(_fitness_model);
  new_trait->set_fitness_scaling_factor(_fitness_scaling_factor);
  new_trait->set_dominance(_dominance);
  new_trait->set_strength(_strength);
  new_trait->set_continuous_effects(_continuous_effects);
  new_trait->set_viability_func_ptr(_fitness_model, _continuous_effects);
  new_trait->set_mut_rate(_mut_rate,_nb_locus);
  new_trait->set_inherit_func_ptr(_inherit_func_ptr);
  new_trait->set_mutation_func_ptr(_mutation_model);
  
  return new_trait;
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TProtoDeletMutations_bitstring::inherit_free (sex_t SEX, bitstring* seq, bitstring** parent)
{ 
  for(unsigned int i = 0; i < _nb_locus; ) 
    seq[i] = (bool)(*parent[RAND::RandBool()])[i];  
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TProtoDeletMutations_bitstring::inherit_low (sex_t SEX, bitstring* seq, bitstring** parent)
{
  
  register unsigned int nbRec, prevLoc = 0, chrm_bloc;
  register bool flipper;
  
  vector< unsigned int >& recTable = _map.getRecLoci(SEX, _mapIndex);
  vector< bool > & firstRecPos = _map.getFirstRecPosition(SEX);
  
  nbRec = recTable.size();
  
  //  cout << "\nTProtoDeletMutations_bitstring::inherit_low -- "<<nbRec<<" recombination events\n";
  
  for(unsigned int c = 0, stride = 0, rec = 0; c < _numChromosome; ++c) {
    
    flipper = firstRecPos[c];
    
    chrm_bloc = stride + _numLociPerChrmsm[c]; //number of loci considered so far
    
    prevLoc = stride; //stride is the first locus of a chromosome
    
    //do recombination chromosome-wise
    for(; recTable[rec] < chrm_bloc && rec < nbRec; rec++) {
      
      //      cout << " -- copy from "<<prevLoc<<" to "<<recTable[rec]<<" (side="<<flipper<<")"<<endl;
      
      seq->copy( *parent[flipper], prevLoc, recTable[rec]);
      
      //      memcpy(&seq[prevLoc], &parent[flipper][prevLoc], (recTable[rec] - prevLoc));
      
      prevLoc = recTable[rec];
      
      flipper = !flipper;
    }
    
    //    cout << " -- copy from "<<prevLoc<<" to "<<chrm_bloc<<" (side="<<flipper<<")"<<endl;
    //copy what's left between the last x-over point and the end of the chrmsme
    seq->copy( *parent[flipper], prevLoc, chrm_bloc);
    
    //    memcpy(&seq[prevLoc], &parent[flipper][prevLoc], (chrm_bloc - prevLoc));
    
    stride += _numLociPerChrmsm[c];
    
  }
  //  cout << "parent chromosomes:\n --0:";
  //  for(unsigned int i = 0; i < _nb_locus; ++i)
  //    cout << (*parent[0])[i];
  //  cout << "\n --1:";
  //  for(unsigned int i = 0; i < _nb_locus; ++i)
  //    cout << (*parent[1])[i];
  //  cout << "\ngamete:\n --0:";
  //  for(unsigned int i = 0; i < _nb_locus; ++i)
  //    cout << (*seq)[i];
  //  cout << "\n";
}
// ----------------------------------------------------------------------------------------
// store_data
// ----------------------------------------------------------------------------------------
void TProtoDeletMutations_bitstring::store_data (BinaryStorageBuffer* saver)
{
  saver->store(&_nb_locus,sizeof(int));
  //store all bool on 1 byte
  char dummy = _continuous_effects;
  saver->store(&dummy, 1);
  if(_continuous_effects) {
    saver->store(_effects[0],_nb_locus * sizeof(float));
    saver->store(_effects[1],_nb_locus * sizeof(float));
  }
}
// ----------------------------------------------------------------------------------------
// retrieve_data
// ----------------------------------------------------------------------------------------
bool TProtoDeletMutations_bitstring::retrieve_data (BinaryStorageBuffer* reader)
{
  unsigned int dummy_int;
  reader->read(&dummy_int,sizeof(int));
  if(dummy_int != _nb_locus ){
    error("TProtoDeletMutations::retrieve_data:nb locus in file differ from parameter value!\n");
    _nb_locus = dummy_int;
  }
  
  char dummy_bool = 0;
  reader->read(&dummy_bool, 1);
  if(dummy_bool != _continuous_effects) {
    error("TProtoDeletMutations::retrieve_data:effects in file differ from parameter value!\n");
    _continuous_effects = dummy_bool;
  }
  
  if(_continuous_effects) {
    reader->read(_effects[0],_nb_locus * sizeof(float));
    reader->read(_effects[1],_nb_locus * sizeof(float));
    TTDeletMutations_bitstring::set_effects(_effects);
  }
  return true;
}

// ------------------------------------------------------------------------------

//                             TTDeletMutations/

// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTDeletMutations_bitstring& TTDeletMutations_bitstring::operator=(const TTrait& T)
{
  const TTDeletMutations_bitstring& TD = dynamic_cast<const TTDeletMutations_bitstring&> (T);
  if(this != &TD) {
    
//    _myProto = TD._myProto;
    _nb_locus = TD._nb_locus;
//    _fitness_model = TD._fitness_model;
//    _fitness_scaling_factor = TD._fitness_scaling_factor;
//    _init_freq = TD._init_freq;
//    _mut_rate = TD._mut_rate;
    //    _recomb_rate = TD._recomb_rate;
//    _strength = TD._strength;
//    _dominance = TD._dominance;
//    _continuous_effects = TD._continuous_effects;
//    _viability_func_ptr = TD._viability_func_ptr;
//    _inherit_func_ptr = TD._inherit_func_ptr;
    //deallocate any previous sequence memory:
    reset();
    //allocate sequence memory:
    init();
    //copy sequence:
    sequence[0]->copy(*TD.sequence[0]);
    sequence[1]->copy(*TD.sequence[1]);
    
    set_value();
  }
  return *this;
}
// ----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TTDeletMutations_bitstring::operator==(const TTrait& T)
{
  if(_type.compare(T.get_type()) != 0) return false;
  
  const TTDeletMutations_bitstring& TD = dynamic_cast<const TTDeletMutations_bitstring&> (T);
  if(this != &TD) {
    if(_nb_locus != TD._nb_locus) return false;
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TTDeletMutations_bitstring::operator!=(const TTrait& T)
{
  if(!((*this) == T))
    return true;
  else
    return false;
}
// ----------------------------------------------------------------------------------------
// set_viability_func_ptr
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::set_viability_func_ptr (unsigned int f_model, bool is_cont ) {  
  switch(f_model) {
    case 1:
    {
      if(is_cont)
        _viability_func_ptr = &TTDeletMutations_bitstring::viability_multi_continuous;
      else
        _viability_func_ptr = &TTDeletMutations_bitstring::viability_multi;
      break;
    }
    case 2:
    {
      if(is_cont)
        _viability_func_ptr = &TTDeletMutations_bitstring::viability_epist_continuous;
      else
        _viability_func_ptr = &TTDeletMutations_bitstring::viability_epist;
      break;
    }
  }
}
// ----------------------------------------------------------------------------------------
// set_mutation_func_ptr
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::set_mutation_func_ptr (unsigned int m_model)
{  
  switch(m_model) {
    case 1:
      _mutation_func_ptr = &TTDeletMutations_bitstring::mutate_noredraw;
      break;
    case 2:
      _mutation_func_ptr = &TTDeletMutations_bitstring::mutate_redraw;
      break;
  }
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::init ()
{  
  if(sequence[0] != NULL) fatal("TTDeletMutations_bitstring::init::sequence[0] is not NULL !\n"); 
  if(sequence[1] != NULL) fatal("TTDeletMutations_bitstring::init::sequence[1] is not NULL !\n");
  
  sequence[0] = new bitstring(_nb_locus);
  sequence[1] = new bitstring(_nb_locus);
  
  _htz = new bitstring(_nb_locus);
  _hmz = new bitstring(_nb_locus);
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::reset ()
{
  if(sequence[0] != NULL) delete sequence[0];
  if(sequence[1] != NULL) delete sequence[1];
  
  sequence[0] = NULL;
  sequence[1] = NULL; 
  
  if(_htz != NULL) delete _htz;
  _htz = 0;
  
  if(_hmz != NULL) delete _hmz;
  _hmz = 0;
}
// ----------------------------------------------------------------------------------------
// set_sequence
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::set_sequence(bitstring** seq)
{
  reset(); init();
  sequence[0]->copy(*seq[0]);
  sequence[1]->copy(*seq[1]);
}
// ----------------------------------------------------------------------------------------
// init_sequence
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::init_sequence ()
{
  if(sequence[0] == NULL) sequence[0] = new bitstring(_nb_locus);
  
  if(sequence[1] == NULL) sequence[1] = new bitstring(_nb_locus);
  
  unsigned int nb_mut, locus;
  
  sequence[0]->reset();
  sequence[1]->reset();
  
  if(_init_freq != 0) {
    
    nb_mut = (unsigned int)(_init_freq*_nb_locus*2);
    
    for(unsigned int i = 0; i < nb_mut; i++) {
      
      do {
        locus = RAND::Uniform( _nb_locus );
        //check if the locus is not already homozygote:
      } while( (*sequence[0])[ locus ] && (*sequence[1])[ locus ] );
      
      if((*sequence[0])[ locus ])
        sequence[1]->set(locus);
      else
        sequence[0]->set(locus);
    }
    
  }
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::inherit (TTrait* mother, TTrait* father)
{
  bitstring** mother_seq = (bitstring**)mother->get_sequence();
  bitstring** father_seq = (bitstring**)father->get_sequence();
  
  (_myProto->* _inherit_func_ptr) (FEM, sequence[FEM], mother_seq);
  
  (_myProto->* _inherit_func_ptr) (MAL, sequence[MAL], father_seq);
}

// ----------------------------------------------------------------------------------------
// mutate_redraw
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::mutate_redraw ()
{
  unsigned int NbMut, mutLocus;
  
  NbMut = (unsigned int)RAND::Poisson(_genomic_mut_rate);
  
  if( (int)_nb_hmz_mutations - (int)NbMut < 0 ) NbMut -= _nb_hmz_mutations;
  
  while(NbMut != 0) {  
    
    do { 
      mutLocus = RAND::Uniform( _nb_locus );
      //check if the locus is not already homozygote:
    } while( (*sequence[0])[mutLocus] && (*sequence[1])[mutLocus] );
    
    if ( !( (*sequence[0])[mutLocus] || (*sequence[1])[mutLocus] ) )
      sequence[RAND::RandBool()]->set(mutLocus);
    else if((*sequence[0])[mutLocus])
      sequence[1]->set(mutLocus);
    else
      sequence[0]->set(mutLocus);
    
    NbMut--;
  }
}
// ----------------------------------------------------------------------------------------
// mutate
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::mutate_noredraw ()
{
  unsigned int NbMut;
  
  NbMut = (unsigned int)RAND::Poisson(_genomic_mut_rate);
  
  while(NbMut != 0) {  
    
    sequence[RAND::RandBool()]->set(  RAND::Uniform( _nb_locus ) );
    
    NbMut--;
  }
}
// ----------------------------------------------------------------------------------------
// set_value
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::set_value ( )
{
  set_nb_mutations();
  set_nb_hmz_mutations();
  set_nb_htz_mutations();
  _phenotype = (this->*_viability_func_ptr)();
}
// ----------------------------------------------------------------------------------------
// set_nb_mutations
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::set_nb_mutations ()
{
  _nb_mutations = sequence[0]->count() + sequence[1]->count();
}
// ----------------------------------------------------------------------------------------
// set_nb_htz_mutations
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::set_nb_htz_mutations ()
{
  (*_htz) = *sequence[0] ^ *sequence[1];
  _nb_htz_mutations = _htz->count();
}
// ----------------------------------------------------------------------------------------
// set_nb_hmz_mutations
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::set_nb_hmz_mutations ()
{
  (*_hmz) = *sequence[0] & *sequence[1];
  _nb_hmz_mutations = _hmz->count();
}
// ----------------------------------------------------------------------------------------
// viability_multi
// ----------------------------------------------------------------------------------------
double TTDeletMutations_bitstring::viability_multi()
{
  //w = (1-s)^hmz * (1-hs)^htz
  return ( pow( 1.0 - _strength, (double)_nb_hmz_mutations ) *
          pow( 1.0 - (_dominance * _strength), (double)_nb_htz_mutations) *
          _fitness_scaling_factor );
}
// ----------------------------------------------------------------------------------------
// viability_epist
// ----------------------------------------------------------------------------------------
double TTDeletMutations_bitstring::viability_epist()
{
  //w = 1 - s(hmz + h*htz)
  return 1.0 - _strength * (_nb_hmz_mutations + _dominance * _nb_htz_mutations) * _fitness_scaling_factor;
}
// ----------------------------------------------------------------------------------------
// viability_multi
// ----------------------------------------------------------------------------------------
double TTDeletMutations_bitstring::viability_multi_continuous()
{
  //w = (1-s)^hmz * (1-hs)^htz
  double fitness = 1.0;
  
  for(unsigned int i = 0; i < _nb_locus; ++i)
    //_effects[0] stores hs, _effects[1] stores s
    fitness *= (1.0f - (*_hmz)[i]*_effects[1][i]) * (1.0f - (*_htz)[i]*_effects[0][i]);
  
  return fitness * _fitness_scaling_factor;
}
// ----------------------------------------------------------------------------------------
// viability_epist
// ----------------------------------------------------------------------------------------
double TTDeletMutations_bitstring::viability_epist_continuous()
{
  //w = 1 - s(hmz + h*htz)
  double fitness = 1.0;
  
  for(unsigned int i = 0; i < _nb_locus; ++i)
    fitness -= (*_hmz)[i]*_effects[1][i] + (*_htz)[i]*_effects[0][i];
  
  return fitness * _fitness_scaling_factor;
}
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::store_data (BinaryStorageBuffer* saver)
{
  size_t wnb = sequence[0]->nb_words();
  size_t bytes = sizeof(bitstring::_ul);
  
  for(size_t i = 0; i < wnb; i++) 
    saver->store(sequence[0]->getword_atIdx(i), bytes);
  
  for(size_t i = 0; i < wnb; i++) 
    saver->store(sequence[1]->getword_atIdx(i), bytes);
}
// ----------------------------------------------------------------------------------------
bool TTDeletMutations_bitstring::retrieve_data (BinaryStorageBuffer* reader)
{
  size_t wnb = sequence[0]->nb_words();
  size_t bytes = wnb * sizeof(bitstring::_ul);
  bitstring::_ul *srce = new bitstring::_ul [wnb];
  
  reader->read(srce, bytes); 
  sequence[0]->set_data(srce,wnb);
  
  reader->read(srce, bytes);
  sequence[1]->set_data(srce,wnb);
  
  return true;
}

// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TTDeletMutations_bitstring::show_up ()
{
  set_nb_mutations();
  set_nb_hmz_mutations();
  set_nb_htz_mutations();
  set_value();
  std::cout<<"\n  Trait type: delet"
  <<"\n       value: "<<*(double*)getValue()
  <<"\nnb mutations: "<<get_nb_mutations()
  <<"\n  nb hmz mut: "<<get_nb_hmz_mutations()
  <<"\n  nb htz mut: "<<get_nb_htz_mutations()
  <<"\n    sequence: "<<std::endl;
  cout<<"0: ";
  for(unsigned int i = 0; i < _nb_locus && i < 64; i++)
    cout<<(*sequence[0])[i];
  cout<<"\n1: ";
  for(unsigned int i = 0; i < _nb_locus && i < 64; i++)
    cout<<(*sequence[1])[i];
  cout<<endl;
}



// ----------------------------------------------------------------------------------------

//                             StatHandler

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTDeletMutBitstrSH::setStatRecorders(std::string& token)
{
#ifdef _DEBUG_
  message("-TTDeletMutBitstrSH::setStatRecorders ");
#endif
  if(token.compare("delet") == 0) {
    
    setStatsForDeletMutations(ALL);
    
  }else if(token.compare("adlt.delet") == 0) {
    
    setStatsForDeletMutations(ADULTS);
    
    
  }else if(token.compare("off.delet") == 0) {
    
    setStatsForDeletMutations(OFFSPRG);
    
    
  } else if(token.compare("viability") == 0) {
    
    setViabStats(ALL);
    
  } else if(token.compare("off.viability") == 0) {
    
    setViabStats(OFFSPRG);
    
  } else if(token.compare("adlt.viability") == 0) {
    
    setViabStats(ADULTS);
    
  } else if(token.compare("meanviab") == 0) {
    
    add("Offspring Viability","off.viab",OFFSPRG, static_cast<unsigned int> (OFFSx),
        0,0,&TTDeletMutBitstrSH::getMeanViability,0,0);
    add("Adults Viability","adlt.viab",ADULTS, static_cast<unsigned int> (ADLTx),
        0,0,&TTDeletMutBitstrSH::getMeanViability,0,0);
    
  } else
    return false;
  
  return true;
}
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrSH::setStatsForDeletMutations(age_t AGE)
{
  if(AGE == ALL) {
    setStatsForDeletMutations(OFFSPRG); setStatsForDeletMutations(ADULTS);
    return;
  }
  
  string prefix = (AGE == ADULTS ? "adlt." : "off.");
  void (TTDeletMutBitstrSH::* setter) () = (AGE == ADULTS ?
                                            &TTDeletMutBitstrSH::setAdultDeletStats :
                                            &TTDeletMutBitstrSH::setOffsprgDeletStats);
  
  add("Frequency of deleterious alleles", prefix + "delfreq",   AGE,0,0,
      &TTDeletMutBitstrSH::getDeletAllFreq,0,0,
      setter);
  
  add("Homozygosity of delet. allele", prefix + "delhmz",   AGE,0,0,
      &TTDeletMutBitstrSH::getDeletAllHmz,0,0,0);
  
  add("Heterozygosity of delet. allele", prefix + "delhtz",   AGE,0,0,
      &TTDeletMutBitstrSH::getDeletAllHtz,0,0,0);
  
  add("Fixed Delet Mutations - Global", prefix + "delfix",   AGE,0,0,
      &TTDeletMutBitstrSH::getFixedDeletLoci,0,0,0);
  
  add("Fixed Delet Mutations - Local", prefix + "delfixp",   AGE,0,0,
      &TTDeletMutBitstrSH::getFixedDeletLociPerPatch,0,0,0);
  
  add("Segregating Delet Mutations - Global", prefix + "delsegr",   AGE,0,0,
      &TTDeletMutBitstrSH::getSegregatingDeletLoci,0,0,0);
  
  add("Segregating Delet Mutations - Local",  prefix + "delsegrp",   AGE,0,0,
      &TTDeletMutBitstrSH::getSegregatingDeletLociPerPatch,0,0,0);
  
  add("Delet Fst", prefix + "delfst",   AGE,0,0,
      &TTDeletMutBitstrSH::getFst,0,0,0);
  
  add("Lethal Equivalents", prefix + "lethequ",   AGE,0,0,
      &TTDeletMutBitstrSH::getLethalEquivalents,0,0,0);
  
  add("Heterosis", "heterosis",   ADULTS,0,0,
      &TTDeletMutBitstrSH::getHeterosis,0,0,0);
  
  add("Genetic Load", "load",   ADULTS,0,0,
      &TTDeletMutBitstrSH::getLoad,0,0,0);
  
}
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrSH::setViabStats(age_t AGE)
{
  if(AGE == ALL) {
    setViabStats(OFFSPRG); setViabStats(ADULTS);
    return;
  }
  
  string prefix = (AGE == ADULTS ? "adlt." : "off.");
  void (TTDeletMutBitstrSH::* setter) () = (AGE == ADULTS ?
                                            &TTDeletMutBitstrSH::setAdultViab :
                                            &TTDeletMutBitstrSH::setOffsprgViab);
  
  add("Average Viability",       prefix + "viab",  AGE, 0, 0,
      &TTDeletMutBitstrSH::getMeanViability,0,0,
      setter);
  
  add("Outbreds Viability",      prefix + "viab.outb",  AGE, 0,0,0,
      &TTDeletMutBitstrSH::getViability,0,0);
  
  add("Inbreds Viability",       prefix + "viab.outw",  AGE, 1,0,0,
      &TTDeletMutBitstrSH::getViability,0,0);
  
  add("Half Sibs Viability",     prefix + "viab.hsib",  AGE, 2,0,0,
      &TTDeletMutBitstrSH::getViability,0,0);
  
  add("Full Sibs Viability",     prefix + "viab.fsib",  AGE, 3,0,0,
      &TTDeletMutBitstrSH::getViability,0,0);
  
  add("Selfed Viability",        prefix + "viab.self",  AGE, 4,0,0,
      &TTDeletMutBitstrSH::getViability,0,0);
  
  add("Outbred btw Proportion",  prefix + "prop.outb",  AGE, 0,0,0,
      &TTDeletMutBitstrSH::getSibProportions,0,0);
  
  add("Outbred wtn Proportion",  prefix + "prop.outw",  AGE, 1,0,0,
      &TTDeletMutBitstrSH::getSibProportions,0,0);
  
  add("Half Sibs Proportion",    prefix + "prop.hsibs",  AGE, 2,0,0,
      &TTDeletMutBitstrSH::getSibProportions,0,0);
  
  add("Full Sibs Proportion",    prefix + "prop.fsibs",  AGE, 3,0,0,
      &TTDeletMutBitstrSH::getSibProportions,0,0);
  
  add("Selfed Proportion",       prefix + "prop.self",  AGE, 4,0,0,
      &TTDeletMutBitstrSH::getSibProportions,0,0);
}
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrFH::FHwrite()
{
  if(!get_pop_ptr()->isAlive()) return;
  
  int nb_locus = this->_FHLinkedTrait->get_nb_locus();
  int patchNbr = get_pop_ptr()->getPatchNbr();
  bitstring** seq;
  Patch* current_patch;
  Individual *ind;
  
  std::string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + get_extension();
  
#ifdef _DEBUG_
  message("TTNeutralGenesFH::FHwrite (%s)\n",filename.c_str());
#endif
  
  ofstream FILE (filename.c_str(), ios::out);
  
  if(!FILE) fatal("could not open DELET output file!!\n");
  
  //FILE<<patchNbr<<" "<<nb_locus<<" "<<2<<" "<<1<<"\n";
  
  FILE<<"pop ";
  
  for(int i = 0; i < nb_locus; i++)
    FILE<<"loc"<<i+1<<" ";
  
  FILE<<"age sex ped origin"<<endl;
  
  if(_FHLinkedTrait->get_iscontinuous()){
    float *s = _FHLinkedTrait->get_s_continous();
    float *hs = _FHLinkedTrait->get_hs_continous();
    FILE<<"-1 ";
    for(int i = 0; i < nb_locus; i++)
      FILE<< s[i] << " " ;
    FILE<<"-1 -1 -1 -1"<< endl;
    
    FILE<<"-1 ";
    for(int i = 0; i < nb_locus; i++)
      FILE<< hs[i] << " " ;
    FILE<<"-1 -1 -1 -1"<< endl;
  }    
  
  for (int i = 0; i < patchNbr; ++i) {
    
    current_patch = get_pop_ptr()->getPatch(i);
    
    for (unsigned int j = 0, size = current_patch->size(FEM, OFFSx); j < size; ++j) {
      
      FILE<<i+1<<" ";
      ind = current_patch->get(FEM, OFFSx, j);
      seq = (bitstring**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
      
      for(int k = 0; k < nb_locus; ++k)
        FILE<<(int)((*seq[0])[k])<<(int)((*seq[1])[k])<<" ";
      
      FILE << OFFSPRG << " " << FEM << " " << ind->getPedigreeClass() << " " << ind->getHome()<<endl;
    }
    
    for (unsigned int j = 0, size = current_patch->size(MAL, OFFSx); j < size; ++j) {
      
      FILE<<i+1<<" ";
      ind = current_patch->get(MAL, OFFSx, j);
      seq = (bitstring**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
      
      for(int k = 0; k < nb_locus; ++k)
        FILE<<(int)((*seq[0])[k])<<(int)((*seq[1])[k])<<" ";
      
      FILE << OFFSPRG << " " << MAL << " "  << ind->getPedigreeClass() << " " << ind->getHome()<<endl;
    }
    
    for (unsigned int j = 0; j < current_patch->size(FEM, ADLTx); ++j) {
      
      FILE<<i+1<<" ";
      ind = current_patch->get(FEM, ADLTx, j);
      seq = (bitstring**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
      
      for(int k = 0; k < nb_locus; ++k)
        FILE<<(int)((*seq[0])[k])<<(int)((*seq[1])[k])<<" ";
      
      FILE << ADULTS << " " << FEM << " " << ind->getPedigreeClass() << " " << ind->getHome()<<endl;
    }
    
    for (unsigned int j = 0; j < current_patch->size(MAL, ADLTx); ++j) {
      
      FILE<<i+1<<" ";
      ind = current_patch->get(MAL, ADLTx, j);
      seq = (bitstring**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
      
      for(int k = 0; k < nb_locus; ++k)
        FILE<<(int)((*seq[0])[k])<<(int)((*seq[1])[k])<<" ";
      
      FILE << ADULTS << " " << MAL << " " << ind->getPedigreeClass() << " " << ind->getHome()<<endl;
    }
  }
  
  FILE.close();
  
}
// ----------------------------------------------------------------------------------------
// FHread
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrFH::FHread (string& filename)
{
  unsigned int nb_locus = _FHLinkedTrait->get_nb_locus();
  unsigned int patchNbr = _pop->getPatchNbr();
  TTDeletMutations_bitstring* trait;
  
  ifstream FILE(filename.c_str(),ios::in);
  
  if(!FILE) fatal("could not open FSTAT input file \"%s\"\n",filename.c_str());
  
  unsigned int genot, age, sex, ped, origin, xline = 0;
  int pop;
  age_idx agex;
  Individual *ind;
  unsigned char all0, all1;
  bitstring *seq[2];  
  seq[0] = new bitstring(nb_locus);
  seq[1] = new bitstring(nb_locus);
  
  double* effects[2];
  int lnbr = 2;
  
  effects[0] = new double [nb_locus];
  effects[1] = new double [nb_locus];
  
  cout<<"file state: "<<FILE<<endl;
  
  string str;
  for (unsigned int i = 0; i < nb_locus+5; ++i) {
    FILE>>str;
  }
  cout<<"str: "<<str<<endl;
  while(FILE>>pop) {
    //FILE;
    cout<<"pop "<<pop<<endl;
    if(pop > (int)patchNbr)
      fatal("Patch number found in file exceeds number of patches in the population.\n");
    
    if(pop == -1) {
      
      for(unsigned int i = 0; i < nb_locus; ++i) {
        FILE>>effects[xline][i];
      }
      
      for(unsigned int i = 0; i < 4; ++i)
        FILE>>str;
      
      cout<<"read effects line "<<xline<<" last was "<<effects[xline][nb_locus-1]<<endl;
      xline++;
    }
    else {
      
      for(unsigned int i = 0; i < nb_locus; ++i) {
        FILE>>genot;
        
        all0 = (unsigned char) genot/10;
        all1 = (unsigned char) genot%10;
        
        if(all0 <= 2) {
          if(all0) seq[0]->set(i);
        } else {
          error("in DELET input file at line %i, locus %i : \
                first allele value %d is greater than 2!\n", lnbr, i+1, all0);
          fatal("Please check the input file.\n");
        }
        
        if(all1 <= 2){
          if(all1) seq[1]->set(i);
        } else {
          error("in DELET input file at line %i, locus %i : \
                second allele value %i is greater than 2!\n", lnbr, i+1, all1);
          fatal("Please check the input file.\n");
        }
      }
      
      FILE >> age >> sex >> ped >> origin;
      
      message("age %i sex %i ped %i origin %i\n", age , sex , ped , origin);
      
      agex = (age == ADULTS ? ADLTx : OFFSx);
      
      ind = _pop->makeNewIndividual(0, 0, static_cast<sex_t> (sex), origin - 1);
      ind->setPedigreeClass((unsigned char)ped);
      ind->setAge(age);
      trait = dynamic_cast<TTDeletMutations_bitstring*> (ind->getTrait(_FHLinkedTraitIndex));
      trait->set_sequence(seq);
      
      ind->show_up();
      
      _pop->getPatch(pop-1)->add(static_cast<sex_t> (sex), agex, ind);
    }
    lnbr++;
    
  }
  cout<<"pop size after loading: "<<_pop->size()<<endl;
  FILE.close();  
  
  delete seq[0];
  delete seq[1];
  delete [] effects[0];
  delete [] effects[1];
}

