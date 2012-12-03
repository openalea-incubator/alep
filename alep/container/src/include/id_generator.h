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

#ifndef __CONTAINER_UTILS_IDGENERATOR_H__
#define __CONTAINER_UTILS_IDGENERATOR_H__

#include <set>

namespace container {
	class IdGenerator {
	public:
		struct InvalidId : public std::exception {
		//small exception raised when attempting to play with invalid ids
			int id;
			InvalidId (int l) : std::exception(), id(l) {};
			virtual ~InvalidId() throw () {};
			InvalidId* copy () { return new InvalidId(id); }
			virtual void rethrow () { throw InvalidId(id); }
		};
	private:
		std::set<int> free_ids;
		//list of ids not used
		//all ids bigger than the last element of the list are free
	private:
		int id_max () {
		//maximum not yet used id
			return *free_ids.rbegin();
		}
	public:
		IdGenerator();
		~IdGenerator() {}

		int get_id ();
		//return a free id

		int get_id (int id);
		//return id if id is free
		//otherwise throw InvalidId

		void release_id (int id);
		//release id to be used later
		//if id is not used throw InvalidId

		void clear ();
		//clear all used ids

		//debug
		void state () const;
	};

};
#endif //__CONTAINER_UTILS_IDGENERATOR_H__
