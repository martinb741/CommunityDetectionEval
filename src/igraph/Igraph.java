package igraph;

import java.util.concurrent.TimeUnit;

import com.google.common.base.Stopwatch;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.sun.jna.ptr.DoubleByReference;

import igraph.IgraphLibrary.igraph_arpack_options_t;
import igraph.IgraphLibrary.igraph_matrix_t;
import igraph.IgraphLibrary.igraph_t;
import igraph.IgraphLibrary.igraph_vector_ptr_t;
import igraph.IgraphLibrary.igraph_vector_t;
import librec.data.MatrixEntry;
import librec.data.SparseMatrix;

public class Igraph {
	
	
	public static SparseMatrix getCommunities(SparseMatrix userMatrix, String algorithm) {
		SparseMatrix membershipsMatrix = null;
		
		Stopwatch sw = Stopwatch.createStarted();

		// get igraph instance
		IgraphLibrary igraph = IgraphLibrary.INSTANCE;
		
		// get the number of users in the user adjacency matrix
		int numUsers = userMatrix.numRows();
		
		// create an undirected graph with numUser nodes
		igraph_t graph = new igraph_t();
		int directed = 0;
		igraph.igraph_empty(graph, numUsers, directed);
		
		// write edges into a map to avoid duplicates and self edges
		Multimap<Integer, Integer> edgeMap = HashMultimap.create();
		for (MatrixEntry e : userMatrix){
			int user1 = e.row();
			int user2 = e.column();
			if(!edgeMap.containsEntry(user2, user1)){
				edgeMap.put(user1, user2);
			}
		}
		
		// create an edge vector, initialize to 0 elements and reserve memory for the number of edges contained in the edge map
		igraph_vector_t edges = new igraph_vector_t();
		igraph.igraph_vector_init(edges, 0);
		int numEdges = edgeMap.size();
		igraph.igraph_vector_reserve(edges, numEdges*2);
		
		// create a weight vector
		igraph_vector_t weights = new igraph_vector_t();
		igraph.igraph_vector_init(weights, 0);
		igraph.igraph_vector_reserve(edges, numEdges);
		
		// append each edge to the vector
		for (int user1 : edgeMap.keySet()){
			for (int user2 : edgeMap.get(user1)){
				double weight = userMatrix.get(user1, user2);
				igraph.igraph_vector_push_back(edges, user1);
				igraph.igraph_vector_push_back(edges, user2);
				igraph.igraph_vector_push_back(weights, weight);
			}
		}
		
		// add edges to the graph
		igraph.igraph_add_edges(graph, edges, null);
		
		long graphConstructTime = sw.elapsed(TimeUnit.SECONDS);
		
        System.out.println("Time to construct graph: " + graphConstructTime);
		int numNodesInGraph = IgraphLibrary.INSTANCE.igraph_vcount(graph);
        System.out.println("Number of nodes in graph: " + numNodesInGraph);

		// run community detection
		igraph_vector_t membershipsVector = getCommunitiesFromGraph(graph, weights, algorithm);
		
		// fill membership information into the memberships matrix
		Table<Integer, Integer, Double> membershipsTable = HashBasedTable.create();
		Multimap<Integer, Integer> membershipsColMap = HashMultimap.create();

//		int membershipsVectorLength = igraph.igraph_vector_size(membershipsVector);
		int maxCommunity = 0;
		
		int membershipsVectorLength = igraph.igraph_vector_size(membershipsVector);
		
		if (membershipsVectorLength == numUsers){
		for (int user = 0; user < numUsers; user++){
			int community = (int) igraph.igraph_vector_e(membershipsVector, user);
			membershipsTable.put(user,community,1.0);
			membershipsColMap.put(community, user);
			if(community > maxCommunity){
				maxCommunity = community;
			}
		}
		}
		else{
			System.out.println("Length of memberships vector is " + membershipsVectorLength + " but should be " + numUsers + ".");
		}
		membershipsMatrix = new SparseMatrix(numUsers, maxCommunity+1, membershipsTable, membershipsColMap);
		
		sw.stop();
		long communityDetectTime = sw.elapsed(TimeUnit.SECONDS) - graphConstructTime;
		
//		int weightsVectorLength = igraph.igraph_vector_size(weights);
//        int numEdgesInGraph = IgraphLibrary.INSTANCE.igraph_ecount(graph);
//        int modularityVectorLength = igraph.igraph_vector_size(modularityVector);
//        int mergesMatrixRows = igraph.igraph_matrix_nrow(mergesMatrix);
//        int mergesMatrixCols = igraph.igraph_matrix_ncol(mergesMatrix);
        int numCommunities = membershipsMatrix.columns().size();
        
//        System.out.println("Number of adjacencies in matrix: " + userMatrix.size());
//        System.out.println("Number of edges in map: " + numEdges);
//        System.out.println("Length of weights vector: " + weightsVectorLength);
//        System.out.println("Number of edges in graph: " + numEdgesInGraph);
        System.out.println("Time to detect communities: " + communityDetectTime);
//        System.out.println("Dimensions of merges matrix: " + mergesMatrixRows + "x" + mergesMatrixCols);
//        System.out.println("Length of modularity vector: " + modularityVectorLength);
//        System.out.println("Length of memberships vector: " + membershipsVectorLength);
//        System.out.println("Highest community id: " + maxCommunity);
        System.out.println("Number of communties: " + numCommunities);
		
		return membershipsMatrix;
	}
	
	private static igraph_vector_t getCommunitiesFromGraph(igraph_t graph, igraph_vector_t weights, String algorithm){
		
		// get igraph instance
		IgraphLibrary igraph = IgraphLibrary.INSTANCE;
		
		// create a vector that will contain the membership information
		igraph_vector_t membershipsVector = new igraph_vector_t();
		igraph.igraph_vector_init(membershipsVector, 0);
		
		// declare some vectors and parameters that are used by some of the algorithms
		igraph_matrix_t mergesMatrix = new igraph_matrix_t();
		igraph.igraph_matrix_init(mergesMatrix, 0, 0);
		
		igraph_vector_t modularityVector = new igraph_vector_t();
		igraph.igraph_vector_init(modularityVector, 0);
		
		igraph_vector_t betweennessVector = new igraph_vector_t();
		igraph.igraph_vector_init(betweennessVector, 0);
		
		igraph_vector_t bridgesVector = new igraph_vector_t();
		igraph.igraph_vector_init(bridgesVector, 0);
		
		igraph_vector_t resultVector = new igraph_vector_t();
		igraph.igraph_vector_init(resultVector, 0);
		
		igraph_vector_t cSizeVector = new igraph_vector_t();
		igraph.igraph_vector_init(cSizeVector, 0);
		
		igraph_matrix_t multilevelMembershipsMatrix = new igraph_matrix_t();
		igraph.igraph_matrix_init(multilevelMembershipsMatrix, 0, 0);
		
		// select algorithm, set parameters and run community detection
		System.out.println("Community detection algorithm: " + algorithm);
		switch (algorithm){
			case "Spinglass":
				/*
				int igraph_community_spinglass(
				   const igraph_t *graph,
			       const igraph_vector_t *weights,
			       igraph_real_t *modularity,	// Pointer to a real number, if not NULL then the modularity score of the solution
			        							// will be stored here.
			       igraph_real_t *temperature,	// Pointer to a real number, if not NULL then the temperature at the end of the
			        							// algorithm will be stored here. 
			       igraph_vector_t *membership, // Pointer to an initialized vector or NULL. If not NULL then the result of the
			        							// clustering will be stored here.
			       igraph_vector_t *csize, 		// Pointer to an initialized vector or NULL. If not NULL then the sizes of the clusters
			        							// will stored here in cluster number order.
			       igraph_integer_t spins,		// Integer giving the number of spins, ie. the maximum number of clusters.
			        							// Default in Paper 25.
	       		   igraph_bool_t parupdate,		// A logical constant, whether to update all spins in parallel.
	       			 							// The default for this argument was FALSE (ie. 0) in the original code.
	       			 							// It is not implemented in the IGRAPH_SPINCOMM_INP_NEG implementation.					
			       igraph_real_t starttemp,		// Real number, the temperature at the start.
			        							// The value of this argument was 1.0 in the original code.
			       igraph_real_t stoptemp,		// Real number, the algorithm stops at this temperature.
			        							// The default was 0.01 in the original code. 
			       igraph_real_t coolfact,		// Real number, the coolinf factor for the simulated annealing.
												// The default was 0.99 in the original code. 
			       igraph_spincomm_update_t update_rule,	// enum: IGRAPH_SPINCOMM_UPDATE_SIMPLE=0, IGRAPH_SPINCOMM_UPDATE_CONFIG=1
			       											// If this is IGRAPH_SPINCOMM_UPDATE_SIMPLE then the random graph (ie. G(n,p)),
			       											// if it is IGRAPH_SPINCOMM_UPDATE then the configuration model is used.
			       											// The configuration means that the baseline for the clustering is a random
			       											// graph with the same degree distribution as the input graph. 
			       igraph_real_t gamma, 		// This defined the weight of the missing and existing links in the quality function
			        							// for the clustering. The default value in the original code was 1.0, which is equal
			        							// weight to missing and existing edges.
			       igraph_spinglass_implementation_t implementation,	// enum: IGRAPH_SPINCOMM_IMP_ORIG=0, IGRAPH_SPINCOMM_IMP_NEG=1
			       								// IGRAPH_SPINCOMM_IMP_ORIG selects the original implementation, this is faster,
			       								// IGRAPH_SPINCOMM_INP_NEG selects a new implementation by Vincent Traag that allows
			       								// negative edge weights. 
			       igraph_real_t gamma_minus);	// Parameter for the IGRAPH_SPINCOMM_IMP_NEG implementation. This specifies the balance
												// between the importance of present and non-present negative weighted edges in a community.
				*
				* O(?)
				*/
//				double temperature = 0;
				int spins = 25;
				int parupdate = 0;
				double starttemp = 1.0;
				double stoptemp = 0.01;
//				double coolfact = 0.99;
				double coolfact = 0.99;
				int updateRule = 1; // IGRAPH_SPINCOMM_UPDATE_CONFIG
				double gamma = 1.0;
				int spinglassImplementation = 0;
				double gamma_minus = 1; // not relevant for IGRAPH_SPINCOMM_IMP_ORIG
				
				// get the components of the graph and run spinglass on each component (for now just take the first component)
				igraph_vector_ptr_t componentsVector = new igraph_vector_ptr_t();
				igraph.igraph_vector_init(componentsVector, 0);
				
				int connectednessMode = 0;
				igraph_vector_t componentMembershipVector = new igraph_vector_t();
				igraph.igraph_vector_init(componentMembershipVector, 0);
				igraph_vector_t componentSizeVector = new igraph_vector_t();
				igraph.igraph_vector_init(componentSizeVector, 0);
				igraph.igraph_clusters(graph, componentMembershipVector, componentSizeVector, null, connectednessMode);
				int numComponents = igraph.igraph_vector_size(componentSizeVector);
				
				if(numComponents == 1){
				igraph.igraph_community_spinglass(graph, weights, null, null, membershipsVector, null, spins,
						parupdate, starttemp, stoptemp, coolfact, updateRule, gamma, spinglassImplementation, gamma_minus);
				}
				else{
					System.out.println("Graph contains " + numComponents + " components. Cannot run Spinglass.");
				}
				break;
			case "LE":
				/*
				int igraph_community_leading_eigenvector(
						const igraph_t *graph,
				        const igraph_vector_t *weights,
						igraph_matrix_t *merges,
						igraph_vector_t *membership,
						igraph_integer_t steps,				// The maximum number of steps to perform.
						igraph_arpack_options_t *options, 
						igraph_real_t *modularity,
						igraph_bool_t start,				// Boolean, whether to use the community structure given in the membership argument as a starting point.
						igraph_vector_t *eigenvalues,		// If not a null pointer, then the eigenvalues calculated along the community structure detection are stored here.
						igraph_vector_ptr_t *eigenvectors,	// If not a null pointer, then the eigenvectors that are calculated in each step of the algorithm, are stored here, in a pointer vector.
						igraph_vector_t *history,			// If not a null pointer, then a trace of the algorithm is stored here, encoded numerically.
				        igraph_community_leading_eigenvector_callback_t *callback,
				        void *callback_extra);
		        *
		        * O(n*n)
				*/
				int leSteps = 100000;
				igraph_arpack_options_t leOptions = new igraph_arpack_options_t();
				igraph.igraph_arpack_options_init(leOptions);
				int leStart = 0;
				igraph.igraph_community_leading_eigenvector(graph, weights, mergesMatrix, membershipsVector, leSteps, leOptions, null, leStart,
						null, null, null, null, null);
				break;
			case "Walktrap":
				/*
				int igraph_community_walktrap(
						  const igraph_t *graph, 
					      const igraph_vector_t *weights,
					      int steps,						// Integer constant, the length of the random walks. Paper uses t=2 and t=5.
					      igraph_matrix_t *merges,
					      igraph_vector_t *modularity, 
					      igraph_vector_t *membership);
				*
				* O(m*n)
			 	*/
				int walktrapSteps = 2;
				igraph.igraph_community_walktrap(graph, weights, walktrapSteps, mergesMatrix, modularityVector, membershipsVector);
				break;
			case "EB":
				/*
				int igraph_community_edge_betweenness(
						  const igraph_t *graph, 
					      igraph_vector_t *result,
					      igraph_vector_t *edge_betweenness,	// Pointer to an initialized vector or NULL. In the former case the edge#
					       										// betweenness of the removed edge is stored here. The vector will be resized
					       										// as needed.
					      igraph_matrix_t *merges,
					      igraph_vector_t *bridges,				// Pointer to an initialized vector of NULL. If not NULL then all edge removals
					       										// which separated the network into more components are marked here. 
					      igraph_vector_t *modularity,
					      igraph_vector_t *membership,
					      igraph_bool_t directed,
					      const igraph_vector_t *weights);
				*
				* O(n*n*n)
				*/
				int directed = 0;
				igraph.igraph_community_edge_betweenness(graph, resultVector, betweennessVector, mergesMatrix, modularityVector,
															bridgesVector, membershipsVector, directed, weights);
				break;
			case "FGMod":
				/* 
				 * int igraph_community_fastgreedy(
				 * 		const igraph_t *graph,				// The input graph. It must be a graph without multiple edges.
				 * 											// This is checked and an error message is given for graphs with multiple edges. 
				 *		const igraph_vector_t *weights,		// Potentially a numeric vector containing edge weights. Supply a null pointer here for unweighted graphs.
				 *											// The weights are expected to be non-negative.
				 *		igraph_matrix_t *merges,			// Pointer to an initialized matrix or NULL, the result of the computation is stored here.
				 *											// The matrix has two columns and each merge corresponds to one merge, the ids of the two
				 *											// merged components are stored. The component ids are numbered from zero and the first n
				 *											// components are the individual vertices, n is the number of vertices in the graph.
				 *											// Component n is created in the first merge, component n+1 in the second merge, etc.
				 *											// The matrix will be resized as needed. If this argument is NULL then it is ignored completely.
				 *		igraph_vector_t *modularity, 		// Pointer to an initialized vector or NULL pointer, in the former case the modularity scores along
				 *											// the stages of the computation are recorded here. The vector will be resized as needed. 
				 *		igraph_vector_t *membership);		// Pointer to a vector. If not a null pointer, then the membership vector corresponding to the
				 *											// best split (in terms of modularity) is stored here.
				 *
				 * O(m)
				 */
				igraph.igraph_community_fastgreedy(graph, weights, mergesMatrix, modularityVector, membershipsVector);
				break;
			case "MLMod":
				/*
				int igraph_community_multilevel(
						  const igraph_t *graph,
						  const igraph_vector_t *weights,
						  igraph_vector_t *membership,		// The membership vector, the result is returned here.
						  igraph_matrix_t *memberships,		// Numeric matrix that will contain the membership vector after each level.
						  igraph_vector_t *modularity);		// Numeric vector that will contain the modularity score after each level
				*
				* O(n) - O(m)
				*/
				igraph.igraph_community_multilevel(graph, weights, membershipsVector, multilevelMembershipsMatrix, modularityVector);
				break;
			case "LP":
				/*
				int igraph_community_label_propagation(const igraph_t *graph,
                        igraph_vector_t *membership,
                        const igraph_vector_t *weights,
                        const igraph_vector_t *initial,	// The initial state. If NULL, every vertex will have a different label at the
                         								// beginning. Otherwise it must be a vector with an entry for each vertex.
                        igraph_vector_bool_t *fixed, 	// Boolean vector denoting which labels are fixed. Of course this makes sense
                         								// only if you provided an initial state, otherwise this element will be ignored.
                        igraph_real_t *modularity);		// The modularity score of the detected community structure is stored here.
				*
				* O(m+n)
				*/
				igraph.igraph_community_label_propagation(graph, membershipsVector, weights, null, null, null);
				break;
			case "Infomap":
				/*
				int igraph_community_infomap(const igraph_t * graph,
                        const igraph_vector_t *e_weights,	// Numeric vector giving the weights of the edges.
                        const igraph_vector_t *v_weights,	// Numeric vector giving the weights of the vertices.
                        int nb_trials,						// The number of attempts to partition the network
                         									// (can be any integer value equal or larger than 1)
                        igraph_vector_t *membership,
                        igraph_real_t *codelength);			// If not NULL the code length of the partition is stored here.
				*
				* O(?)
				*/
				int nbTrials = 5;
				DoubleByReference codelengthRef = new DoubleByReference();
//				igraph_vector_t nodeWeights = new igraph_vector_t();
//				igraph.igraph_vector_init(nodeWeights, igraph.igraph_vector_size(weights));
				igraph.igraph_community_infomap(graph, weights, null, nbTrials, membershipsVector, codelengthRef);
				break;
		}
		
		
		
		return membershipsVector;
	}
}
