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

#include <iostream>
#include "container/id_generator.h"

namespace container {
	IdGenerator::IdGenerator() {
		free_ids.insert(0);
	}

	int IdGenerator::get_id () {
		std::set<int>::iterator free_id=free_ids.begin();
		int ret=*free_id;
		free_ids.erase(free_id);
		if ( free_ids.empty() ) {
			free_ids.insert(ret+1);
		}
		return ret;
	}

	int IdGenerator::get_id (int id) {
		int idmax=id_max();
		if(id==idmax) {
			free_ids.erase(id);
			free_ids.insert(id+1);
			return id;
		}
		else {
			if (id>idmax) {
				for(int i=idmax+1;i<id;++i) {
					free_ids.insert(i);
				}
				free_ids.insert(id+1);
				return id;
			}
			else {
				std::set<int>::iterator id_it=free_ids.find(id);
				if (id_it==free_ids.end()) {
					throw IdGenerator::InvalidId(id);
				}
				else {
					free_ids.erase(id_it);
					return id;
				}
			}
		}
	}

	void IdGenerator::release_id (int id) {
		std::set<int>::iterator id_it=free_ids.find(id);
		if (id_it==free_ids.end()) {
			free_ids.insert(id);
		}
		else {
			throw IdGenerator::InvalidId(id);
		}
	}

	void IdGenerator::clear () {
		free_ids.clear();
		free_ids.insert(0);
	}

	void IdGenerator::state () const {
		std::set<int>::const_iterator it;
		for(it=free_ids.begin();it!=free_ids.end();++it) {
			std::cout << "val: " << *it << std::endl;
		}
	}
};
