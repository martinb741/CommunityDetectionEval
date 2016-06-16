package igraph;

import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import igraph.IgraphLibrary.igraph_t;
import igraph.IgraphLibrary.igraph_vector_t;

public class IgraphTest {

	public static void main(String[] args) {
		
		IgraphLibrary inst = IgraphLibrary.INSTANCE;
		
        PointerByReference verPtr = new PointerByReference();
        IntByReference majPtr = new IntByReference();
        IntByReference minPtr = new IntByReference();
        IntByReference sminPtr = new IntByReference();
        IgraphLibrary.INSTANCE.igraph_version(verPtr, majPtr, minPtr, sminPtr);
        String ver = verPtr.getValue().getString(0);
        int maj = majPtr.getValue();
        int min = minPtr.getValue();
        int smin = sminPtr.getValue();
        System.out.println("Version String: " + ver);
        System.out.println("Major Version: " + maj);
        System.out.println("Minor Version: " + min);
        System.out.println("Subminor Version: " + smin);
        
        System.out.println("---");
        
        // create an undirected graph with 5 nodes
        
        igraph_t graph = new igraph_t();
        int retvalue = IgraphLibrary.INSTANCE.igraph_empty(graph, 5, 0);
        System.out.println("Return value: " + retvalue);
        int numNodes = IgraphLibrary.INSTANCE.igraph_vcount(graph);
        System.out.println("Number of nodes in graph: " + numNodes);
        int numEdges = IgraphLibrary.INSTANCE.igraph_ecount(graph);
        System.out.println("Number of edges in graph: " + numEdges);
        
        System.out.println("---");
        
        // add edges to graph one by one
        
        IgraphLibrary.INSTANCE.igraph_add_edge(graph, 0, 1);
        IgraphLibrary.INSTANCE.igraph_add_edge(graph, 1, 2);
        IgraphLibrary.INSTANCE.igraph_add_edge(graph, 2, 3);

        numNodes = IgraphLibrary.INSTANCE.igraph_vcount(graph);
        System.out.println("Number of nodes in graph: " + numNodes);
        numEdges = IgraphLibrary.INSTANCE.igraph_ecount(graph);
        System.out.println("Number of edges in graph: " + numEdges);

        System.out.println("---");
        
        // create vector
        
        igraph_vector_t vector = new igraph_vector_t();
//        retvalue = IgraphLibrary.INSTANCE.igraph_vector_init(vector, 4);
        retvalue = IgraphLibrary.INSTANCE.igraph_vector_init(vector, 0);
        System.out.println("Init return value: " + retvalue);
        retvalue = IgraphLibrary.INSTANCE.igraph_vector_reserve(vector, 4);
        System.out.println("Reserve return value: " + retvalue);
        
        int vectorSize = IgraphLibrary.INSTANCE.igraph_vector_size(vector);
        System.out.println("Vector size: " + vectorSize);
        
        for (int i = 0; i < vectorSize; i++){
        	double element = IgraphLibrary.INSTANCE.igraph_vector_e(vector, i);
        	System.out.println("Element " + i + ": " + element);
        }
        
//        IgraphLibrary.INSTANCE.igraph_vector_set (vector, 0, 0);
//        IgraphLibrary.INSTANCE.igraph_vector_set (vector, 1, 2);
//        IgraphLibrary.INSTANCE.igraph_vector_set (vector, 2, 3);
//        IgraphLibrary.INSTANCE.igraph_vector_set (vector, 3, 2);
        
        IgraphLibrary.INSTANCE.igraph_vector_push_back(vector, 0);
        IgraphLibrary.INSTANCE.igraph_vector_push_back (vector, 2);
        IgraphLibrary.INSTANCE.igraph_vector_push_back (vector, 3);
        IgraphLibrary.INSTANCE.igraph_vector_push_back (vector, 2);
        
        vectorSize = IgraphLibrary.INSTANCE.igraph_vector_size(vector);
        System.out.println("Vector size: " + vectorSize);
        
        for (int i = 0; i < vectorSize; i++){
        	double element = IgraphLibrary.INSTANCE.igraph_vector_e(vector, i);
        	System.out.println("Element " + i + ": " + element);
        }
        
        System.out.println("---");
        
        // Add edges to graph using vector
        
    	retvalue = IgraphLibrary.INSTANCE.igraph_add_edges(graph, vector, null);
        System.out.println("Return value: " + retvalue);
        
        numEdges = inst.igraph_ecount(graph);
        System.out.println("Number of edges in graph: " + numEdges);
        
        System.out.println("---");
        
        
        
        
	}

}
