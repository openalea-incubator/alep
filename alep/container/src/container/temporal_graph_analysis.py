# -*- python -*-
#
#       OpenAlea.Container
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s):  Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                        Frederic Boudon <frederic.boudon@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################
"""This module helps to analyse TemporalPropertyGraph from Spatial Images."""
 
import types
import numpy as np
from interface.property_graph import IPropertyGraph, PropertyError


def __normalized_parameters(func):
    def wrapped_function(graph, vertex_property, vids = None, rank = 1 , verbose = False):
        """
        :Parameters:
        - 'graph' : a TPG.
        - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
        - 'vids' : by default a vertex id or a list of vertex ids. If 'vids=None' the mean absolute deviation will be computed for all ids present in the graph provided.
        - 'rank' : neighborhood at distance 'rank' will be used.

        :Return:
        - a single value if vids is an interger, or a dictionnary of *keys=vids and *values= "result of applyed fucntion `func`"
        """
        # if a name is given, we use vertex_property stored in the graph with this name.
        if isinstance(vertex_property,str):
            vertex_property = graph.vertex_property(vertex_property)

        # -- If no vids provided we compute the function for all keys present in the vertex_property
        if vids==None:
            vids = vertex_property.keys()

        # if an instancemethod is given, we use create a dictionary for the vids base ont the method.
        if isinstance(vertex_property,types.MethodType):
            tmp_vertex_property = {}
            for vid in graph.vertices():
                tmp_vertex_property[vid] = vertex_property(vid)
            vertex_property = tmp_vertex_property

        if type(vids)==int:
            # for single id, compute single result
            return func(graph, vertex_property, vids, rank)
        else:
            # for set of ids, we compute a dictionary of resulting values.
            l={}
            for k in vids:
                if verbose and k%10==0: print k,'/',len(vids)
                l[k] = func(graph, vertex_property, k, rank, edge_type='s')
            return l

    return  wrapped_function


@__normalized_parameters
def laplacian(graph, vertex_property, vid, rank, edge_type):
    """
    Sub-function computing the laplacian between ONE vertex ('vid') and its neighbors at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = laplacian between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)
        vid_neighborhood.remove(vid)
    else: # if rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)-graph.neighborhood(vid,rank-1, edge_type)

    nb_neighborhood = len(vid_neighborhood)

    result = 0
    ivalue = vertex_property[vid]
    if nb_neighborhood != 0 : # if ==0 it's mean that there is no neighbors for the vertex vid.
        for i in vid_neighborhood:
            result = result + vertex_property[i]
        return ivalue - (result / float(nb_neighborhood))

@__normalized_parameters
def mean_abs_dev(graph, vertex_property, vid, rank, edge_type):
    """
    Sub-function computing the mean sum of absolute difference between ONE vertex ('vid') and its neighbors at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = the mean absolute deviation between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)
        vid_neighborhood.remove(vid)
    else: # if rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)-graph.neighborhood(vid,rank-1, edge_type)

    nb_neighborhood = len(vid_neighborhood)

    result = 0
    ivalue = vertex_property[vid]
    if nb_neighborhood != 0 : # if ==0 it's mean that there is no neighbors for the vertex vid.
        for i in vid_neighborhood:
            result = result + abs(ivalue - vertex_property[i])
        return result / float(nb_neighborhood)


@__normalized_parameters
def change(graph, vertex_property, vid, rank, edge_type):
    """
    Sub-function computing the laplacian between ONE vertex ('vid') and its neighbors at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = laplacian between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)
        vid_neighborhood.remove(vid)
    else: # if rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)-graph.neighborhood(vid,rank-1, edge_type)

    nb_neighborhood = len(vid_neighborhood)

    result = 0
    ivalue = vertex_property[vid]
    if nb_neighborhood != 0 : # if ==0 it's mean that there is no neighbors for the vertex vid.
        for i in vid_neighborhood:
            result = result + vertex_property[i]
        return result - ivalue


def __normalized_temporal_parameters(func):
    def wrapped_function(graph, vertex_property, vids = None, rank = 1, labels_at_t_n = True, check_full_lineage = True, verbose = False):
        """
        :Parameters:
        - 'graph' : a TPG.
        - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
        - 'vids' : by default a vertex id or a list of vertex ids. If 'vids=None' the mean absolute deviation will be computed for all ids present in the graph provided.
        - 'rank' : neighborhood at distance 'rank' will be used.

        :Return:
        - a single value if vids is an interger, or a dictionnary of *keys=vids and *values= Mean absolute deviation
        """
        # if a name is given, we use vertex_property stored in the graph with this name.
        if isinstance(vertex_property,str):
            vertex_property = graph.vertex_property(vertex_property)

        # -- If no vids provided we compute the function for all keys present in the vertex_property
        if vids==None:
            vids = vertex_property.keys()

        if isinstance(vids,int):
            vids=[vids] # for single id, compute single result

        # for a list of ids, we create a dictionary of resulting values from function `func`.
        l={}
        for k,vid in enumerate(vids):
            if verbose and k%10==0: print k,'/',len(vids)
            if check_full_lineage:
                if full_lineage(graph, vid, rank): # Check if ALL descendants up to `rank` exists !
                    l[vid] = func(graph, vertex_property, vid, rank)
            else:
                if full_lineage(graph, vid, 1): # Check if there is at least one descendants for `vid`!
                    l[vid] = func(graph, vertex_property, vid, rank)

        if labels_at_t_n:
            return l
        else:
            m={}
            print "You have asked for labels @ t_n+"+str(rank)
            for vid in l:
                if rank > 1:
                    vid_descendants = graph.descendants(vid ,rank)-graph.descendants(vid,rank-1)
                else:
                    vid_descendants = graph.children(vid)
                
                for id_descendant in vid_descendants:
                    m[id_descendant]=l[vid]
            return m

    return  wrapped_function

@__normalized_temporal_parameters
def temporal_change(graph, vertex_property, vid, rank):
    """
    Sub-function computing the temporal change between ONE vertex ('vid') and its descendants at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.
    
    :Return:
    - a single value = temporal change between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_descendants = graph.children(vid)
    else: # if rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
        vid_descendants = graph.descendants(vid,rank)-graph.descendants(vid,rank-1)

    nb_descendants = len(vid_descendants)
    descendants_value = 0
    vid_value = vertex_property[vid]
    for id_descendant in vid_descendants:
        descendants_value = descendants_value + vertex_property[id_descendant]

    return descendants_value - vid_value

@__normalized_temporal_parameters
def relative_temporal_change(graph, vertex_property, vid, rank):
    """
    Sub-function computing the relative temporal change between ONE vertex ('vid') and its descendants at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.
    
    :Return:
    - a single value = relative temporal change between vertex 'vid' and its neighbors at rank 'rank'.
    """
    return temporal_change(graph, vertex_property, vid, rank).values()[0] / float(vertex_property[vid])

def full_lineage(graph, vid, rank):
    """
    Check if lineage is complete over several ranks. 
    i.e. every decendants cells from `vid` have a lineage up to rank `rank`.
    Suppose that the lineage has been correctly done: a lineage is given to the graph only if we are sure to have all the daugthers from a mother.

    :Parameters:
    - `graph` TPG to browse;
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.
    """
    # -- We create a first list of descendants at rank 1.
    vids_descendants={}
    vids_descendants[1] = graph.descendants(vid,1)
    vids_descendants[1].remove(vid)

    # -- First thing you want to know is if there is a lineage at least at the first order !
    if len(vids_descendants[1]) < 1:
        return False

    # -- If rank == 1 (and there is a lineage but we already did the test!) we suppose that it has been correctly done !
    if rank == 1:
        return True

    # -- To check a full lineage @rank_n :
    for r in xrange(2,rank+1):
        len_lineage_at_r = 0
        for v in vids_descendants[r-1]: # for every mother cell @ rank_n-1 (n>2),
            if len(graph.descendants(v,1)) > 1: # we verify she has a daughter @ rank_n
                len_lineage_at_r += 1
        if len_lineage_at_r != len(vids_descendants[r-1]):
            #~ print len_lineage_at_r,'/', len(vids_descendants[r-1])
            return False
        else:
            # Retreive cells only at rank_n-1: (not all cells between rank 1 and `rank`)
            vids_descendants[r] = graph.descendants(vid,r)-graph.descendants(vid,r-1)

    return True


def time_point_property(graph, time_point, vertex_property):
    """
    Allow to extract a property 'vertex_property' from the temporal graph for one time-point.
    
    :Parameters:
    - `graph` TPG to browse;
    - `time_point` integer defining the time-point to consider;
    - `vertex_property` vertex property to extract;
    :Return:
    - dictionnary of vertex property extracted from the time-point 'time_point';
    """
    # if a name is given, we use vertex_property stored in the graph with this name.
    if isinstance(vertex_property,str):
        vertex_property = graph.vertex_property(vertex_property)
        
    if time_point not in graph.vertex_property('index').values():
        import warnings
        warnings.warn(str(time_point)+"not in"+str(graph))

    k=[i for i in graph.vertex_property('index') if graph.vertex_property('index')[i]==time_point]
    tmp={}
    for i in k:
        tmp[i]=vertex_property[i]
    
    return tmp


def cell_vtx_time_association(graph, return_cells_vertices_relations = False ):
    """
    Creates vrtx2vrtx dictionnary (v2v): associate the corresponding cell vertex over time.
    
    :INPUTS:
        .l21: t_n+1-> t_n (LienTissuTXT) cells lineage
        .vrtx2cell_1: dict at t_n *keys=vertex id ; *values=ids of the 4 associated cells
        .vrtx2cell_2: dict at t_n+1 *keys=vertex id ; *values=ids of the 4 associated cells
    
    :OUPTUT:
        .v2v: dict *keys=t_n+1 vertex number; *values=associated t_n vertex.
    """
    from openalea.image.algo.analysis import cells_vertices_relations
    
    cell2vertices, vtx2coords, vtx2cells = [],[],[]
    for t in xrange( len(graph.graph_property('cell_vertices_coord')) ):
        tmp_1,tmp_2,tmp_3 = cells_vertices_relations(graph.graph_property('cell_vertices_coord')[t])
        cell2vertices.append( tmp_3 )
        vtx2coords.append( tmp_2 )
        vtx2cells.append( tmp_1 )
    
    vtx_time_association = []
    for t in xrange( len(graph.graph_property('cell_vertices_coord'))-1 ):
        v2v = {} ##vertex t_n vers t_n+1
        ## Loop on the vertices label of t_n+1:
        for vtx in vtx2cells[t+1].keys():
            associated_cells = list( vtx2cells[t+1][vtx] )
            ## For the 4 (daugthers) cells associated to this vertex, we temporary replace it by it's mother label:
            for n,cell in enumerate(associated_cells):
                if cell == None: # Mean background...
                    associated_cells.remove(None)
                else:
                    ancestors = graph.ancestors(cell,1)
                    if len(ancestors) != 1: ## if the cell has ancestors (ancestor return a set containing the vertex you ask for)...
                        associated_cells[n] = list(ancestors-set([cell]))[0] ## ...we replace the daughters' label by the one from its mother.
                    else:
                        associated_cells[n] = 0 ## ...else we code by a 0 the absence of a mother in the lineage file (for one -or more- of the 4 daugthers)
            if 0 not in associated_cells: ## If the full topology around the vertex is known:
                associated_cells.sort()
                for k in vtx2cells[t].keys(): 
                    if len(set(vtx2cells[t][k])&set(associated_cells)) == len(associated_cells):
                        v2v[k] = vtx
        vtx_time_association.append(v2v)
    
    if return_cells_vertices_relations:
        return vtx_time_association, cell2vertices, vtx2coords, vtx2cells
    else:
        return vtx_time_association


def __strain_parameters(func):
    def wrapped_function(graph, vids = None, verbose = False):
        """
        :Parameters:
        - 'graph' : a TPG.
        - 'vids' : by default a vertex id or a list of vertex ids. If 'vids=None' the mean absolute deviation will be computed for all ids present in the graph provided.
        """
        # Check if cell vertices have been recovered and associated.
        try :
            graph.graph_property('cell_vertices_coord')
        except:
            import warnings
            warnings.warn("It seems that you don't have detected cell vertices...")
            return 0
        else:
            from openalea.container.temporal_graph_analysis import cell_vtx_time_association
            vtx_time_association, cell2vertices, vtx2coords, vtx2cells = cell_vtx_time_association(graph, return_cells_vertices_relations = True)

        strain_matrix = {}
        if vids == None:
            #~ for t in xrange(max(g.vertex_property('index').values())+1):
            for t in xrange( len(cell2vertices)-1 ):
                for n, cell in enumerate( cell2vertices[t] ):
                    N = len( cell2vertices[t][cell] )
                    # - We make sure that all vertex of the cell 'cell' have been associtates over time.
                    if sum( [(cell2vertices[t][cell][k] in vtx_time_association[t].keys()) for k in xrange(N)] ) == N :
                        # - Now we recover the coordinates @t_n and @t_n+1
                        xyz_t1 = np.array( [vtx2coords[t][cell2vertices[t][cell][k]] for k in xrange(N)] )
                        xyz_t2 = np.array( [vtx2coords[t][vtx_time_association[t][cell2vertices[t][cell][k]]] for k in xrange(N)] )
                        strain_matrix[cell] = func(graph, xyz_t1, xyz_t2, N)
        return strain_matrix

    return  wrapped_function


@__strain_parameters
def strain_matrix(graph, xyz_t1, xyz_t2, N, dimension = 2):
    """
    Suppose that you made sure vid has descendants.
    """
    from numpy.linalg import svd, lstsq

    # -- We start by making sure we can compute the strain matrix in the dimensionnality provided.
    if (dimension > 3) or (dimension < 1) :
        import warnings
        warnings.warn("You can only compute the strain in 1D, 2D or 3D!")
        return None

    ## Compute the centroids:
    c_t1=np.array((np.mean(xyz_t1[:,0]),np.mean(xyz_t1[:,1]),np.mean(xyz_t1[:,2])))
    c_t2=np.array((np.mean(xyz_t2[:,0]),np.mean(xyz_t2[:,1]),np.mean(xyz_t2[:,2])))
    ## Compute the centered matrix:
    c_xyz_t1=np.array(xyz_t1-c_t1)
    c_xyz_t2=np.array(xyz_t2-c_t2)
    ## Compute the Singular Value Decomposition (SVD) of centered coordinates:
    U_t1,D_t1,V_t1=svd(c_xyz_t1, full_matrices=False)
    U_t2,D_t2,V_t2=svd(c_xyz_t2, full_matrices=False)
    V_t1=V_t1.T ; V_t2=V_t2.T
    ## Projection of the vertices' xyz 3D co-ordinate into the 2D subspace defined by the 2 first eigenvector
    #(the third eigenvalue is really close from zero confirming the fact that all the vertices are close from the plane -true for external part of L1, not for inner parts of the tissue).
    c_xy_t1=np.array([np.dot(U_t1[k,0:dimension],np.diag(D_t1)[0:dimension,0:dimension]) for k in range(N)])
    c_xy_t2=np.array([np.dot(U_t2[k,0:dimension],np.diag(D_t2)[0:dimension,0:dimension]) for k in range(N)])
    ## Least-square estimation of A.
    #A is the transformation matrix in the regression equation between the centered vertices position of two time points:
    lsq=lstsq(c_xy_t1,c_xy_t2)
    A=lsq[0]

    return A


def strain_rate(strain_matrix, deltaT = 24):
    """
    Compute the strain rate: sr[c] = np.log(D_A[0])/deltaT , np.log(D_A[1])/deltaT
    Dimensionnality is imposed by the one of the strain matrix.
    """
    from numpy.linalg import svd
    sr = {}
    for c in strain_matrix:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q=svd(strain_matrix[c])
        # Compute Strain Rates :
        sr[c] = np.log(D_A)/deltaT

    return sr


def areal_strain_rate(strain_matrix, deltaT = 24):
    """
    Compute the areal strain rate: asr[c] = sum(np.log(D_A[0])/deltaT , np.log(D_A[1])/deltaT
    Dimensionnality is imposed by the one of the strain matrix.
    """
    from numpy.linalg import svd
    asr = {}
    for c in strain_matrix:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q=svd(strain_matrix[c])
        # Compute Strain Rates and Areal Strain Rate:
        asr[c] = sum(np.log(D_A)/deltaT)

    return asr

def anisotropy(strain_matrix, deltaT = 24):
    """
    Compute the "Growth" Anisotropy = (sr1-sr2)/(sr1+sr2)
    Dimensionnality is imposed by the one of the strain matrix.
    """
    from numpy.linalg import svd
    anisotropy = {}
    for c in strain_matrix:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q=svd(strain_matrix[c])
        # Compute Strain Rates :
        strain_rate = np.log(D_A)/deltaT
        anisotropy[c] = (strain_rate[0]-strain_rate[1])/(strain_rate[0]+strain_rate[1])

    return anisotropy

@__strain_parameters
def strain_cross(graph, xyz_t1, xyz_t2, N, dimension = 2):
    """
    Suppose that you made sure vid has descendants.
    """
    from numpy.linalg import svd, lstsq

    # -- We start by making sure we can compute the strain matrix in the dimensionnality provided.
    if (dimension > 3) or (dimension < 1) :
        import warnings
        warnings.warn("You can only compute the strain in 1D, 2D or 3D!")
        return None

    ## Compute the centroids:
    c_t1=np.array((np.mean(xyz_t1[:,0]),np.mean(xyz_t1[:,1]),np.mean(xyz_t1[:,2])))
    c_t2=np.array((np.mean(xyz_t2[:,0]),np.mean(xyz_t2[:,1]),np.mean(xyz_t2[:,2])))
    ## Compute the centered matrix:
    c_xyz_t1=np.array(xyz_t1-c_t1)
    c_xyz_t2=np.array(xyz_t2-c_t2)
    ## Compute the Singular Value Decomposition (SVD) of centered coordinates:
    U_t1,D_t1,V_t1=svd(c_xyz_t1, full_matrices=False)
    U_t2,D_t2,V_t2=svd(c_xyz_t2, full_matrices=False)
    V_t1=V_t1.T ; V_t2=V_t2.T
    ## Projection of the vertices' xyz 3D co-ordinate into the 2D subspace defined by the 2 first eigenvector
    #(the third eigenvalue is really close from zero confirming the fact that all the vertices are close from the plane -true for external part of L1, not for inner parts of the tissue).
    c_xy_t1=np.array([np.dot(U_t1[k,0:dimension],np.diag(D_t1)[0:dimension,0:dimension]) for k in range(N)])
    c_xy_t2=np.array([np.dot(U_t2[k,0:dimension],np.diag(D_t2)[0:dimension,0:dimension]) for k in range(N)])
    ## Compute the Singular Value Decomposition (SVD) of the least-square estimation of A.
    #A is the (linear) transformation matrix in the regression equation between the centered vertices position of two time points:
    lsq=lstsq(c_xy_t1,c_xy_t2)
    ##  Singular Value Decomposition (SVD) of A.
    R,D_A,Q=svd(lsq[0])
    Q=Q.T
    ##  Getting back in 3D: manually adding an extra dimension.
    R=np.hstack([np.vstack([R,[0,0]]),[[0],[0],[1]]])
    D_A=np.hstack([np.vstack([np.diag(D_A),[0,0]]),[[0],[0],[0]]])
    Q=np.hstack([np.vstack([Q,[0,0]]),[[0],[0],[1]]])
    ##  Getting back in 3D: strain of cell c represented at each time point.
    s_t1 = np.dot(np.dot(np.dot(np.dot(V_t1, R), D_A), R.T), V_t1.T)
    s_t2 = np.dot(np.dot(np.dot(np.dot(V_t2, Q), D_A), Q.T), V_t2.T)

    return s_t1,s_t2


#~ def display_lignage(graph, ini):
    
#~ def strain2D(graph, tp_1, tp_2):
    #~ """
    #~ Strain computation based on the 3D->2D->3D GOODALL method.
    #~ 
    #~ :INPUTS:
        #~ .t1: t_n Spatial Image containing cells (segmented image)
        #~ .t2: t_n+1 Spatial Image containing cells (segmented image)
        #~ .l12: lineage between t_n & t_n+1;
        #~ .l21: INVERTED lineage between t_n & t_n+1;
        #~ .deltaT: time interval between two time points;
        #~ 
    #~ :Variables:
        #~ .v2v_21: vertex (keys=t_n+1) to vertex (values=t_n) association.
        #~ .c2v_1: cells 2 vertex @ t_n
        #~ .v2b_1: vextex 2 barycenters @ t_n
        #~ .v2b_2: vextex 2 barycenters @ t_n+1
    #~ 
    #~ :OUTPUTS: (c= keys= mother cell number)
        #~ .sr[c]: Strain Rate = np.log(D_A[0])/deltaT , np.log(D_A[1])/deltaT
        #~ .asr[c]: Areal Strain Rate = (sr1+sr2)
        #~ .anisotropy[c]: Growth Anisotropy = (sr1-sr2)/(sr1+sr2)
        #~ .s_t1[c]: t_n strain cross in 3D (tensor)
        #~ .s_t2[c]: t_n+1 strain cross in 3D (tensor)
    #~ 
    #~ ########## Relationship between least-squares method and principal components: ##########
    #~ ## The first principal component about the mean of a set of points can be represented by that line which most closely approaches the data points 
    #~ #(as measured by squared distance of closest approach, i.e. perpendicular to the line).
    #~ ## In contrast, linear least squares tries to minimize the distance in the y direction only.
    #~ ## Thus, although the two use a similar error metric, linear least squares is a method that treats one dimension of the data preferentially, while PCA treats all dimensions equally.
    #~ #########################################################################################
    #~ """
    #~ from numpy.linalg import svd, lstsq
#~ 
    #~ ## Variable creation used to comput the strain.
    #~ v2v_12 = dict((v,k) for k, v in v2v_21.items())
    #~ lsq={}
    #~ s_t1,s_t2={},{}
    #~ sr={}
    #~ asr={}
    #~ anisotropy={}
#~ 
    #~ for c in l12.keys():
        #~ if c in c2v_1.keys():
            #~ if sum([(c2v_1[c][k] in v2v_12.keys()) for k in range(len(c2v_1[c]))])==len(c2v_1[c]):
                #~ N = len(c2v_1[c])
                #~ if N>2:
                    #~ ## Retreive positions of the vertices belonging to cell 'c':
                    #~ xyz_t1=np.array([v2b_1[c2v_1[c][k]] for k in range(N)])
                    #~ xyz_t2=np.array([v2b_2[v2v_12[c2v_1[c][k]]] for k in range(N)])
                    #~ ## Compute the centroids:
                    #~ c_t1=np.array((np.mean(xyz_t1[:,0]),np.mean(xyz_t1[:,1]),np.mean(xyz_t1[:,2])))
                    #~ c_t2=np.array((np.mean(xyz_t2[:,0]),np.mean(xyz_t2[:,1]),np.mean(xyz_t2[:,2])))
                    #~ ## Compute the centered matrix:
                    #~ c_xyz_t1=np.array(xyz_t1-c_t1)
                    #~ c_xyz_t2=np.array(xyz_t2-c_t2)
                    #~ ## Compute the Singular Value Decomposition (SVD) of centered coordinates:
                    #~ U_t1,D_t1,V_t1=svd(c_xyz_t1, full_matrices=False)
                    #~ U_t2,D_t2,V_t2=svd(c_xyz_t2, full_matrices=False)
                    #~ V_t1=V_t1.T ; V_t2=V_t2.T
                    #~ ## Projection of the vertices' xyz 3D co-ordinate into the 2D subspace defined by the 2 first eigenvector
                    #~ #(the third eigenvalue is really close from zero confirming the fact that all the vertices are close from the plane -true for external part of L1, not for inner parts of the tissue).
                    #~ c_xy_t1=np.array([np.dot(U_t1[k,0:2],np.diag(D_t1)[0:2,0:2]) for k in range(N)])
                    #~ c_xy_t2=np.array([np.dot(U_t2[k,0:2],np.diag(D_t2)[0:2,0:2]) for k in range(N)])
                    #~ ## Compute the Singular Value Decomposition (SVD) of the least-square estimation of A.
                    #~ #A is the (linear) transformation matrix in the regression equation between the centered vertices position of two time points:
                    #~ lsq[c]=lstsq(c_xy_t1,c_xy_t2)
                    #~ ##  Singular Value Decomposition (SVD) of A.
                    #~ R,D_A,Q=svd(lsq[c][0])
                    #~ Q=Q.T
                    #~ # Compute Strain Rates and Areal Strain Rate:
                    #~ sr[c] = np.log(D_A)/deltaT
                    #~ asr[c] = sum(sr[c])
                    #~ anisotropy[c]=((sr[c][0]-sr[c][1])/asr[c])
                    #~ ##  Getting back in 3D: manually adding an extra dimension.
                    #~ R=np.hstack([np.vstack([R,[0,0]]),[[0],[0],[1]]])
                    #~ D_A=np.hstack([np.vstack([np.diag(D_A),[0,0]]),[[0],[0],[0]]])
                    #~ Q=np.hstack([np.vstack([Q,[0,0]]),[[0],[0],[1]]])
                    #~ ##  Getting back in 3D: strain of cell c represented at each time point.
                    #~ s_t1[c] = np.dot(np.dot(np.dot(np.dot(V_t1, R), D_A), R.T), V_t1.T)
                    #~ s_t2[c] = np.dot(np.dot(np.dot(np.dot(V_t2, Q), D_A), Q.T), V_t2.T)
#~ 
    #~ return sr,asr,anisotropy,s_t1,s_t2



def triplot(graphs_list, values2plot, labels_list=None, values_name="",normed=False):
    """
    TO DO
    """
    import numpy as np
    if labels_list==None:
        labels_list=[]
        for g in graphs_list:
            labels_list.append(g.vertex_property('label'))
    
    values=[]
    abs_dev_values=[]
    laplacian_values=[]
    #-- if 'values2plot' is a string, it must be a property found in all graphs in the 'graph_list'.
    if type(values2plot)==type(str('str')):
        for g in graphs_list:
            if values2plot not in g.vertex_property_names():
                import sys
                sys.exit(1)
            else:
                if (values_name==""):
                    values_name=values2plot
                values.append(g.vertex_property(values2plot).values())
                abs_dev_values.append(dev_abs(g,values2plot,True))
                laplacian_values.append(laplacian(g,values2plot,True))

    import matplotlib.pyplot as plt
    fig = plt.figure()
    fig.subplots_adjust( wspace=0.13, left=0.05, right=0.95, top=0.95)
    main=fig.add_subplot(1,2,1)
    main.hist(values, bins=20,normed=normed,
        label=( ('t1, n='+str(len(values[0]))+', mean='+str(np.round(np.mean(values[0]), 2))) ,
        ('t2, n='+str(len(values[1]))+', mean='+str(np.round(np.mean(values[1]), 2))) ,
        ('t3, n='+str(len(values[2]))+', mean='+str(np.round(np.mean(values[2]), 2))) ), histtype='bar' )
    plt.title("L1 cells' "+values_name)
    if values_name=='volume':
        plt.xlabel('Volumes'+ r' ($\mu m^3$)')
    else:
        plt.xlabel(values_name)
    if normed:
            plt.ylabel('Frequency')
    else:
        plt.ylabel('Number of observations')
    plt.legend()
    
    dev=fig.add_subplot(2,2,2)
    dev.hist(abs_dev_values, bins=20,normed=normed,
        label=( ('t1, n='+str(len(abs_dev_values[0]))+', mean='+str(np.round(np.mean(abs_dev_values[0]), 2))) ,
        ('t2, n='+str(len(abs_dev_values[1]))+', mean='+str(np.round(np.mean(abs_dev_values[1]), 2))) ,
        ('t3, n='+str(len(abs_dev_values[2]))+', mean='+str(np.round(np.mean(abs_dev_values[2]), 2))) ), histtype='bar' )
    plt.title("L1 cells' absolute deviance from neighbors in "+values_name)
    plt.xlabel('Deviance from neighbors in volumes'+ r' ($\mu m^3$)')
    if normed:
            plt.ylabel('Frequency')
    else:
        plt.ylabel('Number of observations')
    plt.legend()
    
    lap=fig.add_subplot(2,2,4)
    lap.hist(laplacian_values, bins=20,normed=normed,
        label=( ('t1, n='+str(len(laplacian_values[0]))+', mean='+str(np.round(np.mean(laplacian_values[0]), 2))) ,
        ('t2, n='+str(len(laplacian_values[1]))+', mean='+str(np.round(np.mean(laplacian_values[1]), 2))) ,
        ('t3, n='+str(len(laplacian_values[2]))+', mean='+str(np.round(np.mean(laplacian_values[2]), 2))) ), histtype='bar' )
    plt.title("L1 cells' laplacian from neighbors in "+values_name)
    plt.xlabel('Laplacian from neighbors in volumes'+ r' ($\mu m^3$)')
    if normed:
            plt.ylabel('Frequency')
    else:
        plt.ylabel('Number of observations')
    plt.legend()
    plt.show()
