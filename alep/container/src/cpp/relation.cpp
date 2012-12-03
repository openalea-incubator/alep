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
#include "container/relation.h"

namespace container {
	void Relation::state () const {
		std::cout << "size :" << size() << std::endl;
		std::cout << "links :";
		for(Relation::const_iterator it=begin();it!=end();++it) {
			int link_id=*it;
			std::cout << "link: " << link_id << " source: " << source(link_id) << "target: " << target(link_id) << std::endl;
		}
	}
};
