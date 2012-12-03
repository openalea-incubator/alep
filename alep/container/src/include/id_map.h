/* -*-c++-*- 
 *------------------------------------------------------------------------------
 *                                                                              
 *        openalea.container: utils package                                     
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

#ifndef __CONTAINER_IDMAP_H__
#define __CONTAINER_IDMAP_H__


#include <utility>
#include "container/id_generator.h"
#include "container/custom_map.h"

namespace container {
#define MAP CustomMap<int,D>

	template <typename D>
	class IdMap : public MAP {
	public:
		typedef typename MAP::value_type value_type;//type of data in this map
		typedef typename MAP::iterator iterator;//used to traverse all elements : iter of (int,D)
		typedef typename MAP::const_iterator const_iterator;//const version of previous iterator
		typedef typename MAP::key_iterator key_iterator;//iterator on keys in the map : iter of int
		typedef typename MAP::const_key_iterator const_key_iterator;//iterator on keys in the map : iter of int

	
	public:
		typedef IdGenerator::InvalidId InvalidId;//error thrown when playing with invalid ids

	private:
		IdGenerator id_generator;//generator used to create ids when not provided
		//set of methods that can not be used
		void swap (const IdMap<D>& map) {}
		D& operator[] (const int& key) {}

	public:
		//inheritance from map
		IdMap () {}
		//template <typename InputIterator>
		//IdMap(InputIterator f, InputIterator l); TODO

		IdMap (const MAP& ref) {
		//fill this map with values taken from ref
			for(iterator it=ref.begin();it!=ref.end();++it) {
				add(it->second,it->first);
			}
		}

		~IdMap () {}

		std::pair<iterator,bool> insert (const value_type& x) {
		//try to insert a new (int,D) pair in the map
		//if the key is already in the map throw IdError
			int key = id_generator.get_id(x->first);
			return MAP::insert(x);
		}

		void erase (iterator pos) {
		//remove element at position pos
		//release the corresponding id
			id_generator.release_id(pos->first);
			MAP::erase(pos);
		}

		void erase (const int& key) {
		//remove element which id is key
		//if key not in the map
		//throw IdError
			id_generator.release_id(key);
			MAP::erase(key);
		}

		void clear () {
		//clear all elements from the map
		//and release all used ids
			id_generator.clear();
			MAP::clear();
		}

		int add (const D& data) {
		//add a new element in the map
		//create a new id to be used as key
		//return the id used to store the data
			int key=id_generator.get_id();
			MAP::insert(typename IdMap<D>::value_type(key,data));
			return key;
		}

		int add (const D& data, const int& key) {
		//try to add a new element in th emap
		//using the provided id as a key
		//if key is already used throw IdError
			int used_key=id_generator.get_id(key);
			MAP::insert(typename IdMap<D>::value_type(used_key,data));
			return used_key;
		}

		void setitem (const int& key, const D& data) {
		//return the value stored with key
			try {
				id_generator.get_id(key);
			}
			catch (IdGenerator::InvalidId) {
			}
			MAP::setitem(key,data);
		}

	};

#undef MAP
};
#endif //__CONTAINER_IDMAP_H__

