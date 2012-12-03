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

#ifndef __CONTAINER_INT_ITERATOR_H__
#define __CONTAINER_INT_ITERATOR_H__

namespace container {
	class IntIterator {
	//base class for all iterator that reference an integer
	public:
		typedef int value_type;
	public:
		IntIterator() {}
		~IntIterator() {}
		virtual int operator* () const {
		//dereferencing operator
			return 0;
		}
		
		virtual IntIterator& operator++ () {
		//post increment operator
			return *this;
		}

		virtual IntIterator operator++ (int) {
		//pre increment operator
		}

		bool operator== (const IntIterator& other) const {
		//comparison operator
		//used stored value to compare
			return this->operator*() == other.operator*() ;
		}

		bool operator!= (const IntIterator& other) const {
		//comparison operator
		//used stored value to compare
			return this->operator*() != other.operator*();
		}

		//debug
		void state () const {}
	};

	template<typename T_iterator>
	class IntCloneIterator : public IntIterator {
	private:
		T_iterator& current_it;
	public:
		IntCloneIterator (T_iterator& ref_it) : current_it(ref_it) {}
		~IntCloneIterator () {}
		int operator* () const {
			return *current_it;
		}
		IntIterator& operator++ () {
			current_it++;
			return *this;
		}
		IntIterator operator++ (int) {
			IntCloneIterator<T_iterator> ret=*this;
			++current_it;
			return ret;
		}
	};
};
#endif //__CONTAINER_INT_ITERATOR_H__
