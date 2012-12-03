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

#ifndef __CONTAINER_CUSTOM_MAP_H__
#define __CONTAINER_CUSTOM_MAP_H__


#include <iostream>
#include <utility>
#include <map>

namespace container {
#define MAP std::map<K,D>

	template <typename K, typename D>
	class MapKeyIterator {
	private:
		typedef typename MAP::iterator mapit;//for internal use only
		typedef MapKeyIterator<K,D> self;//for internal use only
	public:
		typedef K value_type;//the type of data returned when dereferenced
	public:
		mapit map_it;//map iterator
	public:
		MapKeyIterator() {}
		//blanck constructor used
		//because I don't know how to initialise
		//iterators in for loops otherwise
		//do not use
		//
		MapKeyIterator(const mapit& ref_it) : map_it(ref_it) {
		//copy ref_it inside self
		}

		~MapKeyIterator() {}

		D& data () const {
		//return the second argument of this iterator
			return map_it->second;
		}

		K operator* () {
		//dereference operator
		//return first argument of the pair owned
		//by map_it
			return map_it->first;
		}

		self& operator++ () {
		//post increment operator
			map_it++;
			return *this;
		}

		self operator++ (int) const {
		//preincrement operator
			self ret=*this;
			++map_it;
			return ret;
		}

		bool operator== (const self& other) const {
		//comparaison operator
		//use comparaison operator of owned iterator
			return map_it == other.map_it;
		}
		bool operator!= (const self& other) {
		//comparaison operator
		//use comparaison operator of owned iterator
			return map_it != other.map_it;
		}
	};

	template <typename K, typename D>
	class ConstMapKeyIterator {
	private:
		typedef typename MAP::const_iterator mapit;//for internal use only
		typedef ConstMapKeyIterator<K,D> self;//for internal use only
	public:
		typedef K value_type;//the type of data returned when dereferenced
	public:
		mapit map_it;//map iterator
	public:
		ConstMapKeyIterator() {}
		//blanck constructor used
		//because I don't know how to initialise
		//iterators in for llops otherwise
		//do not use
		//
		ConstMapKeyIterator(const mapit& ref_it) : map_it(ref_it) {
		//copy ref_it inside self
		}

		~ConstMapKeyIterator() {}

		const D& data () const {
		//return the second argument of this iterator
			return map_it->second;
		}

		K operator* () const {
		//dereference operator
		//return first argument of the pair owned
		//by map_it
			return map_it->first;
		}

		self& operator++ () {
		//post increment operator
			map_it++;
			return *this;
		}

		self operator++ (int) const {
		//preincrement operator
			self ret=*this;
			++map_it;
			return ret;
		}

		bool operator== (const self& other) const {
		//comparaison operator
		//use comparaison operator of owned iterator
			return map_it == other.map_it;
		}
		bool operator!= (const self& other) {
		//comparaison operator
		//use comparaison operator of owned iterator
			return map_it != other.map_it;
		}
	};

	template <typename K, typename D>
	class CustomMap : public MAP {
	public:
		typedef typename MAP::key_type key_type;//type of key in this map
		typedef D data_type;//type of data stored in this map
		typedef typename MAP::value_type value_type;//type of elements stored in this map

		typedef typename MAP::iterator iterator;//used to traverse all elements : iter of (int,D)
		typedef typename MAP::const_iterator const_iterator;//const version of previous iterator
		typedef MapKeyIterator<K,D> key_iterator;//iterator on keys in the map : iter of int
		typedef ConstMapKeyIterator<K,D> const_key_iterator;//iterator on keys in the map : iter of int

	
	public:
		//inheritance from map
		CustomMap () {}
		CustomMap (iterator ref_begin, iterator ref_end) : MAP(ref_begin,ref_end) {}
		~CustomMap () {}

		bool has_key (const int& key) {
		//test wether a key is inside the map
			return MAP::find(key)!=MAP::end();
		}

		//specific functions
		key_iterator key_begin () {
		//iterator on all key in the map
			return key_iterator(MAP::begin());
		}

		const_key_iterator key_begin () const {
		//iterator on all key in the map
			return const_key_iterator(MAP::begin());
		}

		key_iterator key_end () {
		//iterator on all keys in the map
			return key_iterator(MAP::end());
		}

		const_key_iterator key_end () const {
		//iterator on all keys in the map
			return const_key_iterator(MAP::end());
		}

		const D& getitem (const K& key) {
		//return the value stored with key
			typename MAP::iterator it=MAP::find(key);
			if(it==MAP::end()) {
				throw key;
			}
			else {
				return it->second;
			}
		}

		void setitem (const K& key, const D& data) {
		//store data in the map using key
			MAP::operator[](key)=data;
		}

		void delitem (const K& key) {
		//erase element at position key
			MAP::erase(key);
		}
		D pop (const K& key) {
		//remove and return value associated with key
			D ret=getitem(key);
			MAP::erase(key);
			return ret;
		}

		iterator popitem () {
			iterator it=MAP::begin();
			MAP::erase(it);
			return it;
		}

		CustomMap<K,D> copy () {
			return CustomMap(MAP::begin(),MAP::end());
		}

		const D& setdefault (const int& key, const D& default_value) {
			if(has_key(key)) {
				return getitem(key);
			}
			else {
				setitem(key,default_value);
				return default_value;
			}
		}

	};

#undef MAP
};
#endif //__CONTAINER_CUSTOM_MAP_H__

