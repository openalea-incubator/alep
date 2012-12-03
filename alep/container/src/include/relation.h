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

#ifndef __CONTAINER_RELATION_H__
#define __CONTAINER_RELATION_H__

#include <utility>
#include <map>
#include "container/id_map.h"

namespace container {
	class Relation {
	//maintain a set of links between elements
	private:
		typedef IdMap< std::pair<int,int> > idmap;
	public:
		typedef idmap::key_iterator iterator;//iterator on link ids : iter of int
		typedef idmap::const_key_iterator const_iterator;//iterator on link ids : iter of int

	public:
		struct InvalidLink : public idmap::InvalidId {
		//small exception raised when attempting to access an invalid link
			InvalidLink (int l) : idmap::InvalidId(l) {};
			virtual ~InvalidLink() throw () {};
			InvalidLink* copy () { return new InvalidLink(id); }
			virtual void rethrow () { throw InvalidLink(id); }
		};
	private:
		idmap link_elms;//map of link_id:(source,target)
	public:
		Relation() {}
		~Relation() {}
		bool has_link (int link_id) const {
		//test wether a link is inside this relation
			return link_elms.find(link_id)!=link_elms.end();
		}

		int size () const {
		//number of links in this object
			return link_elms.size();
		}

		iterator begin () {
		//iterator on all links
			return link_elms.key_begin();
		}

		const_iterator begin () const {
		//iterator on all links
			return link_elms.key_begin();
		}

		iterator end () {
		//iterator on all links
			return link_elms.key_end();
		}

		const_iterator end () const {
		//iterator on all links
			return link_elms.key_end();
		}

		int source (const_iterator link_it) const {
		//source element of a link
		//if the link is not inside this
		//raise InvalidLink
			if(link_it==end()) {
				throw InvalidLink(*link_it);
			}
			else {
				return link_it.data().first;
			}
		}

		int source (int link_id) const {
		//source element of a link
		//if the link is not inside this
		//raise InvalidLink
			return source(const_iterator(link_elms.find(link_id)));
		}

		int target (const_iterator link_it) const {
		//source element of a link
		//if the link is not inside this
		//raise InvalidLink
			if(link_it==end()) {
				throw InvalidLink(*link_it);
			}
			else {
				return link_it.data().second;
			}
		}

		int target (int link_id) const {
		//target element of a link
		//if the link is not inside this
		//raise InvalidLink
			return target(const_iterator(link_elms.find(link_id)));
		}

		int add_link (int source_elm, int target_elm) {
		//add a new link between source_elm and target_elm
		//return id used for this link
			return link_elms.add(std::make_pair(source_elm,target_elm));
		}

		int add_link (int source_elm, int target_elm, int link_id) {
		//try to add a new link between source_elm and target_elm
		//using link_id
		//if link_id is already used, throw InvalidLink
			try {
				return link_elms.add(std::make_pair(source_elm,target_elm),link_id);
			}
			catch (idmap::InvalidId e) {
				throw InvalidLink(link_id);
			}
		}

		void remove_link (int link_id) {
		//remove a link from this object
		//if the link does not exist throw InvalidLink
			try {
				return link_elms.erase(link_id);
			}
			catch (idmap::InvalidId e) {
				throw InvalidLink(link_id);
			}
		}

		void clear () {
		//remove all links from this relation
			link_elms.clear();
		}

		//debug
		void state () const;
	};
};
#endif //__CONTAINER_RELATION_H__
