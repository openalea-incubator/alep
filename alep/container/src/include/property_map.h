/* -*-c++-*- 
 *------------------------------------------------------------------------------
 *                                                                              
 *        openalea.container: property map package                                     
 *                                                                              
 *        Copyright 2006 INRIA - CIRAD - INRA                      
 *                                                                              
 *        File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>         
 *                                                                              
 *        Distributed under the Cecill-C License.                               
 *        See accompanying file LICENSE.txt or copy at                          
 *            http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html       
 *                                                                              
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr                    
 *       
 *                                                                       
 *-----------------------------------------------------------------------------*/

#ifndef __CONTAINER_PROPERTY_MAP_H__
#define __CONTAINER_PROPERTY_MAP_H__


#include <utility>
#include "container/custom_map.h"

namespace container {
	typedef CustomMap<int,float> fmap;

	class PropertyMap : public fmap {
	public:
		PropertyMap () {}
		void iadd (float val) {
			for(fmap::iterator it=begin();it!=end();++it) {
				it->second+=val;
			}
		}

		void isub (float val) {
			for(fmap::iterator it=begin();it!=end();++it) {
				it->second-=val;
			}
		}

		void imul (float val) {
			for(fmap::iterator it=begin();it!=end();++it) {
				it->second*=val;
			}
		}

		void idiv (float val) {
			for(fmap::iterator it=begin();it!=end();++it) {
				it->second/=val;
			}
		}

	};

	class FactorMapIterator {
	private:
		typedef FactorMapIterator self;
	public:
		typedef fmap::iterator::value_type value_type;
	public:
		fmap::iterator current_it;
		float factor;
		int first;
		float second;
	public:
		FactorMapIterator (fmap::iterator it, float ref_factor) : current_it(it),
									  factor(ref_factor),
	       								  first(0),
									  second(1.) {
			first=current_it->first;
			second=current_it->second;
		}
		~FactorMapIterator () {}
		value_type operator* () {
			//return value_type(current_it->first,current_it->second * factor);
			return value_type(first,second);
		}
		const self* operator-> () const {
			return this;
		}
		self& operator++ () {
			current_it++;
			first=current_it->first;
			second=current_it->second*factor;
			return *this;
		}
		self operator++ (int) {
			self ret=*this;
			++current_it;
			first=current_it->first;
			second=current_it->second*factor;
			return ret;
		}
		bool operator== (const self& other) const {
			return current_it == other.current_it;
		}
		bool operator!= (const self& other) const {
			return current_it != other.current_it;
		}
	};

	class FactorMap {
	public:
		typedef fmap associated_map;
		typedef FactorMapIterator iterator;
		typedef fmap::key_iterator key_iterator;
	private:
		associated_map& ref_prop;
		float factor;
	public:
		FactorMap (associated_map& ref_map, float mul_factor) : ref_prop(ref_map),
									factor(mul_factor) {}
		~FactorMap () {}
		float getitem (int key) {
			return ref_prop.getitem(key)*factor;
		}
		void setitem (int key, float val) {
			ref_prop.setitem(key,val/factor);
		}
		iterator begin () {
			return FactorMapIterator(ref_prop.begin(),factor);
		}
		iterator end () {
			return FactorMapIterator(ref_prop.end(),factor);
		}
		key_iterator key_begin () {
			return ref_prop.key_begin();
		}
		key_iterator key_end () {
			return ref_prop.key_end();
		}

	};

	class DensityMap {
	};

	class QuantityAdaptor {
	};
};
#endif //__CONTAINER_PROPERTY_MAP_H__

