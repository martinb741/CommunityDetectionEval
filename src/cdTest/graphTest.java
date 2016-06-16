package cdTest;

import i5.las2peer.services.ocd.graphs.CustomGraph;
import y.base.Edge;
import y.base.Node;

public class graphTest {

	public static void main(String[] args) {
		CustomGraph graph = new CustomGraph();

		Node node1 = graph.createNode();
		Node node2 = graph.createNode();
		Node node3 = graph.createNode();
		Node node4 = graph.createNode();
		graph.setNodeName(node1, "1");
		graph.setNodeName(node2, "2");
		graph.setNodeName(node3, "3");
		graph.setNodeName(node4, "4");
		
		Edge edge1 = graph.createEdge(node1, node2);
		Edge edge2 = graph.createEdge(node2, node3);
		Edge edge3 = graph.createEdge(node3, node1);
		Edge edge4 = graph.createEdge(node2, node4);
		graph.setEdgeWeight(edge1, 2);
		graph.setEdgeWeight(edge2, 1);
		graph.setEdgeWeight(edge3, 1);
		graph.setEdgeWeight(edge4, 4);

		System.out.println(graph.toString());
	}

}
