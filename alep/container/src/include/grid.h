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

#ifndef __CONTAINER_GRID_H__
#define __CONTAINER_GRID_H__

#include <vector>
#include <exception>

namespace container {
	class IntRangeIterator {
	//an iterator that go from one integer value to another
	private:
		typedef IntRangeIterator self;//for internal use
	public:
		typedef int value_type;
	private:
		int current_value;//current integer value
	public:
		IntRangeIterator(int value=0) {
		//initialise the iterator at value
			current_value=value;
		}
		~IntRangeIterator() {}
		int operator* () const {
		//dereferencing operator
			return current_value;
		}
		
		self& operator++ () {
		//post increment operator
			current_value++;
			return *this;
		}

		self operator++ (int) {
		//pre increment operator
			self ret=*this;
			++current_value;
			return ret;
		}

		bool operator== (const IntRangeIterator& other) const {
		//comparison operator
		//used stored value to compare
			return current_value == *other;
		}

		bool operator!= (const IntRangeIterator& other) const {
		//comparison operator
		//used stored value to compare
			return current_value != *other;
		}

		//debug
		void state () const {}
	};

	class Grid {
	//regular multi dimensional grid
	public:
		typedef std::vector<int> coord_list;//type of coordinates
		typedef IntRangeIterator iterator;//iterator on indexes : iter of int
	public:
		struct OutOfBoundError : public std::exception {
		//small exception raised when attempting to access an invalid cell
			int cid;
			int cid_max;
			OutOfBoundError (int l, int lmax) : std::exception(), cid(l), cid_max(lmax) {};
			virtual ~OutOfBoundError() throw () {};
			OutOfBoundError* copy () { return new OutOfBoundError(cid,cid_max); }
			virtual void rethrow () { throw OutOfBoundError(cid,cid_max); }
		};
	private:
		coord_list max_coords;//number of elements in each dimension of the grid
		coord_list offset_values;//internal offset for coordinates of cases
	public:
		Grid() {}
		Grid(const coord_list& shape);
		//create a grid with shape[i] cases in the ieme dimension

		~Grid() {}

		int dim () const {
		//dimension of the space
			return max_coords.size();
		}

		void set_shape (const coord_list& shape);
		//set a new shape for the grid

		const coord_list& shape () const {
		//number of elements in each dimension
			return max_coords;
		}

		int size () const;
		//total number of cells in the grid

		iterator begin () const {
		//iterator on cell indexes
			return iterator(0);
		}

		iterator end () const {
		//iterator on cell indexes
			return iterator(size());
		}

		int index (const coord_list& coord) const;
		//compute index of a cell from its coordinates
		//if coord[i]>shape[i] throw OutOfBoundError

		coord_list coordinates (int ind) const;
		//compute coordinates of a cell from its index
		//if index<0 or index>size()
		//throw OutOfBoundError

		coord_list coordinates (const iterator& it) const {
		//compute coordinates of a cell from an iterator
		//on its index
		//if index<0 or index>size()
		//throw OutOfBoundError
			return coordinates(*it);
		}

		//debug
		void state () const;
	};
};
#endif //__CONTAINER_GRID_H__
