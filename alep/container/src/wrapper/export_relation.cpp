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

#include "container/relation.h"
#include "container/relation_iterators.h"
using namespace container;

#include <boost/python.hpp>
#include "export_iterator.h"
using namespace boost::python;

typedef PyCustomRange<Relation::const_iterator> relation_link_range;

relation_link_range relation_iter_links (const Relation& rel) {
	return relation_link_range(rel.begin(),rel.end());
}

void export_relation () {
	export_custom_range<Relation::const_iterator>("_PyConstRelationLinkRange");
	export_custom_range<Relation::iterator>("_PyRelationLinkRange");
	export_custom_range<LinkSourceIterator>("_PyRelationLinkSourceRange");
	export_custom_range<LinkTargetIterator>("_PyRelationLinkTargetRange");

	class_<Relation>("Relation", "dual relation between entities")
		.def("size",&Relation::size,"number of relations")
		.def("links",&relation_iter_links,"iterator on all link_ids")
		.def("source",(int (Relation::*) (int) const)& Relation::source,"source element of a link")
		.def("target",(int (Relation::*) (int) const)& Relation::target,"target element of a link")
		.def("add_link",(int (Relation::*) (int,int))& Relation::add_link,"add a new link between two elements")
		.def("add_link",(int (Relation::*) (int,int,int))& Relation::add_link,"add a new link between two elements")
		.def("remove_link",&Relation::remove_link,"remove the specified link")
		.def("clear",&Relation::clear,"remove all links")
		.def("state",&Relation::state,"debug function");

}
