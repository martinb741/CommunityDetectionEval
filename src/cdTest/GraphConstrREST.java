package cdTest;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import javax.ws.rs.client.Client;
import javax.ws.rs.client.Entity;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;

import org.apache.commons.dbutils.DbUtils;
import org.glassfish.jersey.apache.connector.ApacheConnectorProvider;
import org.glassfish.jersey.client.ClientConfig;
import org.glassfish.jersey.client.ClientProperties;
import org.glassfish.jersey.client.JerseyClientBuilder;

import com.google.common.base.Stopwatch;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.LinkedListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.Table;

import greedyFiltering.GFSorting;
import greedyFiltering.GFSparseMatrix;
import greedyFiltering.GreedyFiltering;
import librec.data.MatrixEntry;
import librec.data.SparseMatrix;
import librec.data.SparseVector;
import librec.data.VectorEntry;

public class GraphConstrREST {

	public static void main(String[] args) {
		SparseMatrix ratings = null;
		
		boolean proxy = false;
		int k = 10;
		int mu = 300;
		
//		for (int maxIds : new int[] {10, 100, 1000,2000,4000,6000,8000,10000,12000, 20000, 40000, 60000, 80000, 100000, 120000}){
//		for (int maxIds : new int[] {20000, 40000, 60000, 80000, 100000, 120000}){
//		for (int maxIds : new int[] {1000,4000,8000,12000,20000,40000,60000,80000,100000,120000}){
		for (int maxIds : new int[] {60000,80000,100000,120000}){
			System.out.println("======");
			try {

				// Get ratings matrix from database
				ratings = loadData(maxIds);
				
//				String nodeEdgeList = getUserNodeEdgeList(ratings);
//				SparseMatrix userMatrix = getUserMatrix(ratings);
				
			} catch (Exception e1) {
				e1.printStackTrace();
//				System.out.println("Exception: " + e1.getMessage());
			}
			
			// Build user graph from ratings matrix, store in a sparse adjacency matrix
			SparseMatrix userMatrix = getUserMatrix3(ratings, k, mu);
			
			// Create node edge list that can be processed by WebOCD from that adjacency matrix
			String nodeEdgeList = getUserNodeEdgeList(userMatrix);
			
			// Send graph to WebOCD
			int numUsers = ratings.numRows();
			uploadGraph(nodeEdgeList, numUsers, k, mu, proxy);
			
		}
		
	}
	
	private static void uploadGraph(String nodeEdgeList, int numUsers, int k, int mu, boolean proxy) {
		Stopwatch sw = Stopwatch.createStarted();
		
		Client client;
		
		if (proxy){
	        ClientConfig config = new ClientConfig();
	        config.connectorProvider(new ApacheConnectorProvider());
	        config.property(ClientProperties.PROXY_URI, "http://localhost:8082");
			client = JerseyClientBuilder.newClient(config);
		}
		else{
			client = JerseyClientBuilder.newClient();
		}
		
		String name = "Graph-"+numUsers+"u-"+k+"k-"+mu+"mu";
		
		WebTarget target = client.target("http://localhost:8080").path("ocd").path("graphs")
				.queryParam("name", name)
				.queryParam("creationType", "REAL_WORLD")
				.queryParam("inputFormat", "NODE_WEIGHTED_EDGE_LIST")
				.queryParam("doMakeUndirected", "TRUE");
		Response response = target.request(MediaType.TEXT_XML).post(Entity.entity(nodeEdgeList, MediaType.TEXT_PLAIN));
		
		sw.stop();
		
		long uploadTime = sw.elapsed(TimeUnit.SECONDS);

		System.out.println("Response Status: " + response.getStatus());
		System.out.println("Response Text: " + response.readEntity(String.class));
		System.out.println("Upload time: " + uploadTime);
	}

	private static String getUserNodeEdgeList(SparseMatrix adjacencies){
		
		System.out.println("Constructing list of nodes and edges ...");

		Stopwatch sw = Stopwatch.createStarted();

		StringBuilder list = new StringBuilder();
		
		int numUsers = adjacencies.numRows();
		
		// Add each user to the list
		for (int u = 0; u < numUsers; u++){
			list.append(u + "\n");
		}		
		
		for (MatrixEntry e : adjacencies){
			list.append(e.row() + " " + e.column() + " " + e.get() + "\n");
		}
			
		sw.stop();
		
		long computeTime = sw.elapsed(TimeUnit.SECONDS);
		
		System.out.println("Node edge list computation time: " + computeTime);
		
		return list.toString();
	}
	
	private static SparseMatrix getUserMatrix(SparseMatrix ratings){
		
		System.out.println("Computing user matrix ...");

		Stopwatch sw = Stopwatch.createStarted();
		
		SparseMatrix userMatrix = null;

		BiMap<Integer, Integer> userIds = HashBiMap.create();
		Table<Integer, Integer, Double> userAdjacencyTable = HashBasedTable.create();
		Multimap<Integer, Integer> colMap = HashMultimap.create();
		
		int numUsers = 0;
		int numAdjacencies = 0;
		
		if (ratings == null){
			return null;
		}
		
		// Add a table entry for each user
		for (int u = 0; u < ratings.numRows(); u++){
			numUsers++;
			userIds.put(u, u);
		}		
		
		for (int u1 = 0; u1 < ratings.numRows()-1; u1++){
			for (int u2 = u1+1; u2 < ratings.numRows(); u2++){
				double sim = getSimilarity(ratings.row(u1), ratings.row(u2));
				if (sim != 0){
					userAdjacencyTable.put(u1,u2,sim);
					colMap.put(u2, u1);
					numAdjacencies++;
				}
			}
		}
		
		userMatrix = new SparseMatrix(ratings.numRows(), ratings.numRows(), userAdjacencyTable, colMap);
		
		sw.stop();
		
		long computeTime = sw.elapsed(TimeUnit.SECONDS);

		System.out.println("Users: " + numUsers + " / Adjacencies: " + numAdjacencies + " / Time: " + computeTime);
		
		return userMatrix;
	}
	
	private static SparseMatrix getUserMatrix2(SparseMatrix ratings, int mu){
		
		System.out.println("Computing user matrix ...");

		if (ratings == null){
			return null;
		}
		
		Stopwatch sw = Stopwatch.createStarted();
		
		int numUsers = ratings.numRows();
		int numItems = ratings.numColumns();
		
		// Matrix that will be returned from this method, contains the user adjacencies
		SparseMatrix adjacencyMatrix = null;

		// Table from which the adjacency matrix will be created
		Table<Integer, Integer, Double> adjacencyTable = HashBasedTable.create();
		// Column map to build a matrix with both compressed column and compressed row storage
		Multimap<Integer, Integer> adjacencyColMap = HashMultimap.create();
		
		// Maps each user to a dimension-feature map that maps each item that was rated by the user to a feature derived from the rating.
		// The inner map must be sorted by the feature value for the KNN algorithm to work
//		Map<Integer, Map<Integer,Double>> featureMap = new HashMap<Integer, Map<Integer,Double>>();

		Map<Integer, double[]> featureMap = new HashMap<Integer, double[]>();
		Map<Integer, int[]> dimensionMap = new HashMap<Integer, int[]>();

		int numAdjacencies = 0;
		
		
		// Get user features and store in featureMap
		
		for (int u = 0; u < ratings.numRows(); u++){
			SparseVector userRatings = ratings.row(u);
			Map<Integer, Double> userFeatureMap = new HashMap<Integer, Double>();
			
			double normFactor = 0;
			
			for (VectorEntry ratingEntry : userRatings){
				double rating = ratingEntry.get();
				normFactor += rating * rating;
			}
			normFactor = Math.sqrt(normFactor);
			
			for (VectorEntry ratingEntry : userRatings){
				int dimension = ratingEntry.index();
				double rating = ratingEntry.get();
				double feature = rating / normFactor;
				userFeatureMap.put(dimension, feature);
			}
			userFeatureMap = MapUtils.sortByValueDesc(userFeatureMap);
			double[] userFeatureArray = new double[userFeatureMap.size()];
			int[] userDimensionArray = new int[userFeatureMap.size()];
			
			int index = 0;
			for (Map.Entry<Integer, Double> e : userFeatureMap.entrySet()){
				userFeatureArray[index] = e.getValue();
				userDimensionArray[index] = e.getKey();
			}
//			featureMap.put(u, userFeatureMap);
			featureMap.put(u, userFeatureArray);
			dimensionMap.put(u, userDimensionArray);
		}
		
		
		// Perform Greedy Filtering
		
		// Maps items to a list of users whose feature vector prefix contains this item (dimension) 
		LinkedListMultimap<Integer, Integer> L = LinkedListMultimap.create();
		
		int c = 0; // counts iterations
		
		// Contains the users whose vectors are to be processed
		HashSet<Integer> R = new HashSet<Integer>();
		for (int u = 0; u < numUsers; u++){
			R.add(u);
		}
		
		do{
			// Find candidates
			for (int u : R){
				int dimension = dimensionMap.get(u)[c];
				L.put(dimension, u);
//				System.out.println("L: added (dim " + dimension + ", user " + u + ")");
			}
			// Check stop condition
			HashSet<Integer> removeFromR = new HashSet<Integer>();
			for (int u : R){
				int m = 0;
				for (int j = 0; j < c; j++){
					int dimension = dimensionMap.get(u)[j];
					m += L.get(dimension).size();
				}
				if (m > mu || c >= dimensionMap.get(u).length - 1){
					removeFromR.add(u);
				}
			}
			R.removeAll(removeFromR);
			c++;
		}
		while(!R.isEmpty());
		
		
		// Brute-force search in L
		
		HashMultimap<Integer, Integer> userPairs = HashMultimap.create();
		for (int i = 0; i < numItems; i++){
			// Store each pair of users that need to be compared into the userPairs map with the lower user id as the key
			// ...
			List<Integer> users = L.get(i);
			
			// can be made more efficient by not iterating through the whole list in the inner loop
			for (int user1 : users){
				for (int user2 : users){
					if (user1 < user2){
						userPairs.put(user1, user2);
					}
					else if (user2 < user2){
						userPairs.put(user2, user1);
					}
				}
			}
		}
		
		
//		for (int u1 = 0; u1 < ratings.numRows()-1; u1++){
//			for (int u2 = u1+1; u2 < ratings.numRows(); u2++){
//				double sim = getSimilarity(ratings.row(u1), ratings.row(u2));
//				if (sim != 0){
//					adjacencyTable.put(u1,u2,sim);
//					adjacencyColMap.put(u2, u1);
//					numAdjacencies++;
//				}
//			}
//		}
		
		adjacencyMatrix = new SparseMatrix(ratings.numRows(), ratings.numRows(), adjacencyTable, adjacencyColMap);
		
		sw.stop();
		
		long computeTime = sw.elapsed(TimeUnit.SECONDS);
		
		// Collect debug information
		
		// ratings: Average number of ratings per user
		double avgRatings = 0;
		for (int u = 0; u < numUsers; u++){
			avgRatings += (double) ratings.row(u).sum() / (double) numUsers;
		}
		
		// featureMap: number of users, average number of features per user
		int numFeatureUsers = featureMap.size();
		double avgFeatures = 0;
		for (Map.Entry<Integer, double[]> e : featureMap.entrySet()){
			avgFeatures += (double) e.getValue().length / (double) numFeatureUsers;
		}
		// dimensionMap: number of users, average number of features per user
		int numDimensionUsers = dimensionMap.size();
		double avgDimensions = 0;
		for (Map.Entry<Integer, int[]> e : dimensionMap.entrySet()){
			avgDimensions += (double) e.getValue().length / (double) numDimensionUsers;
		}
		// L: number of user-dimension pairs
		int numUserDimensionPairs = L.entries().size();
		// userPairs: number of user pairs
		int numUserPairs = userPairs.entries().size();

		System.out.println("Users: " + numUsers + " / "
				+ "Avg ratings/user: " + avgRatings + " / "
				+ "Adjacencies: " + numAdjacencies + " / "
				+ "Time: " + computeTime);
		System.out.println("Users in featureMap: " + numFeatureUsers + " / "
				+ "Avg features/user: " + avgFeatures + " / "
				+ "Users in dimensionMap: " + numDimensionUsers + " / "
				+ "Avg dimensions/user: " + avgDimensions + " / "
				+ "Pairs in L: " + numUserDimensionPairs + " / "
				+ "User Pairs: " + numUserPairs
				);
		
		return adjacencyMatrix;
	}
	
	private static SparseMatrix getUserMatrix3(SparseMatrix ratings, int k, int mu){
		
		System.out.println("Computing user matrix ...");

		if (ratings == null){
			return null;
		}
		
		Stopwatch sw = Stopwatch.createStarted();
		
		// Matrix that will be returned from this method, contains the user adjacencies
		SparseMatrix adjacencyMatrix = null;
		// Table to create the matrix from
		Table<Integer, Integer, Double> adjacencyTable = HashBasedTable.create();
		// Column map to allow compressed row and column storage
		Multimap<Integer, Integer> adjacencyColMap = HashMultimap.create();
		
		int numUsers = ratings.numRows();
		int numItems = ratings.numColumns();
		int numAdjacencies = 0;
		
		
		// Construct GFSparseMatrix from ratings matrix
		
		int numElements = ratings.size();
		GFSparseMatrix inputMatrix = new GFSparseMatrix(numUsers, numElements, numItems);

		int vectorCounter = 0;
		int elementCounter = 0;
		
		// Get value frequencies (number of vectors containing each value
		int[] frequencies = new int[numItems];
		for (int u = 0; u < numUsers; u++){
			for (VectorEntry element : ratings.row(u)){
				int dimension = element.index();
				frequencies[dimension]++;
			}
		}

		for (int u = 0; u < numUsers; u++){
			inputMatrix.setVectorIndex(vectorCounter, elementCounter); 
			
			// Get max element value
			double maxValue = 0;
			for (VectorEntry element : ratings.row(u)){
				double value = element.get();
				if (maxValue < value){
					maxValue = value;
				}
			}
			
			for (VectorEntry element : ratings.row(u)){
				int dimension = element.index();
				double value = element.get();
				
				double tfidf = (0.5 + (0.5 * value / maxValue)) * Math.log((double) numItems / (double) frequencies[dimension]);
				
				inputMatrix.setDimension(elementCounter, dimension);
				inputMatrix.setValue(elementCounter, tfidf);
				
				elementCounter++;
			}

			vectorCounter++;
		}
		inputMatrix.setVectorIndex(vectorCounter, elementCounter);
		
		// Sort each vector by value
		GFSorting.sortByValue(inputMatrix);
		
		// Run KNNGraphConstruction - kNNGraphConstruction(SparseMatrix inputMatrix, double[][] kNNGraph, int k, int mu)
		double[][] kNNGraph = new double[numUsers][];
		for (int i=0; i<numUsers; i++) {
			kNNGraph[i] = new double[k*2];
			for (int j=0; j<k*2; j++) {
				kNNGraph[i][j] = -1;
			}
		}
		GreedyFiltering.kNNGraphConstruction(inputMatrix, kNNGraph, k, mu);
		
		
		// Convert double[][] containing KNN graph into a sparse adjacency matrix
		// kNNGraph has dimensions [numUsers][k*2]
		// e.g. kNNGraph[5][8]=17, kNNGraph[5][9]=0.35 means that user 5 has neighbor 17 with similarity 0.35
		
		for (int u1 = 0; u1 < numUsers; u1++){
			for (int idx = 0; idx < k*2; idx += 2){
				int u2 = (int) kNNGraph[u1][idx];
				if (u2 != -1){
					double sim = kNNGraph[u1][idx+1];
					adjacencyTable.put(u1, u2, sim);
					adjacencyColMap.put(u2, u1);
				}
			}
		}
		
		adjacencyMatrix = new SparseMatrix(numUsers, numUsers, adjacencyTable, adjacencyColMap);
		sw.stop();
		
		long computeTime = sw.elapsed(TimeUnit.SECONDS);
		
		// Collect debug information
		
		// ratings: Average number of ratings per user
		double avgRatings = 0;
		for (int u = 0; u < numUsers; u++){
			avgRatings += (double) ratings.row(u).sum() / (double) numUsers;
		}
		
		// adjacencyMatrix: Average number of neighbors
		double avgNeighbors = (double) adjacencyMatrix.size() / (double) numUsers;
		
		System.out.println("Users: " + numUsers + " / "
				+ "Avg ratings/user: " + avgRatings + " / "
				+ "Avg neighbors/user: " + avgNeighbors + " / "
				+ "Adjacencies: " + numAdjacencies + " / "
				+ "Time: " + computeTime);
		
		return adjacencyMatrix;
	}
	
	private static double getSimilarity(SparseVector u1, SparseVector u2){
		
		for (VectorEntry ve : u1){
			if (ve.get() != 0 && u2.get(ve.index()) != 0){
				return 1;
			}
		}
		
		return 0;
	}
	
	private static SparseMatrix loadData(int maxIds) throws Exception {
		
		System.out.println("Loading ratings ...");
		
		Stopwatch sw = Stopwatch.createStarted();

		// store rate data as {user, item, rate} matrix
		SparseMatrix rateMatrix;
		// store time data as {user, item, timestamp} matrix
//		SparseMatrix timeMatrix;

		List<Double> ratingScale;
		Multiset<Double> scaleDist = HashMultiset.create();

//		private BiMap<String, Integer> userIds, itemIds;
		BiMap<Integer, Integer> userIds = HashBiMap.create(), itemIds = HashBiMap.create();

		TimeUnit timeUnit = TimeUnit.SECONDS;

		long minTimestamp, maxTimestamp;
		
		
		// Table {row-id, col-id, rate}
		Table<Integer, Integer, Double> dataTable = HashBasedTable.create();
		// Table {row-id, col-id, timestamp}
		Table<Integer, Integer, Long> timeTable = null;
		// Map {col-id, multiple row-id}: used to fast build a rating matrix
		Multimap<Integer, Integer> colMap = HashMultimap.create();

		DatabaseManager dbm = new DatabaseManager();
		Connection conn = null;
		PreparedStatement stmnt = null;
		ResultSet rs = null;

		minTimestamp = Long.MAX_VALUE;
		maxTimestamp = Long.MIN_VALUE;

		try {
			conn = dbm.getConnection();

			if (maxIds > 0){
				stmnt = conn.prepareStatement("SELECT UserId, ItemId, Time, Rating FROM Rating WHERE ItemId < ? AND UserId < ?;");
				stmnt.setInt(1, maxIds);
				stmnt.setInt(2, maxIds);
			}
			else{
				stmnt = conn.prepareStatement("SELECT UserId, ItemId, Time, Rating FROM Rating;");
			}
			rs = stmnt.executeQuery();
			while(rs.next()){
				
				int user = rs.getInt("UserId");
				int item = rs.getInt("ItemId");
				
				// convert time to milliseconds
				long mms = rs.getTimestamp("Time").getTime();
				long timestamp = timeUnit.toMillis(mms);

				Double rate = rs.getDouble("Rating");
				
				scaleDist.add(rate);

				int row = userIds.containsKey(user) ? userIds.get(user) : userIds.size();
				userIds.put(user, row);

				int col = itemIds.containsKey(item) ? itemIds.get(item) : itemIds.size();
				itemIds.put(item, col);

				dataTable.put(row, col, rate);
				colMap.put(col, row);
				
				if (timeTable == null)
					timeTable = HashBasedTable.create();

				if (minTimestamp > timestamp)
					minTimestamp = timestamp;

				if (maxTimestamp < timestamp)
					maxTimestamp = timestamp;

				timeTable.put(row, col, timestamp);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		finally{
			DbUtils.closeQuietly(stmnt);
			DbUtils.closeQuietly(conn);
			DbUtils.closeQuietly(rs);
		}

		ratingScale = new ArrayList<>(scaleDist.elementSet());
		Collections.sort(ratingScale);

		int numRows = userIds.size(), numCols = itemIds.size();

		double minRate = ratingScale.get(0).doubleValue();
		double epsilon = minRate == 0.0 ? ratingScale.get(1).doubleValue() - minRate : 0;
		if (epsilon > 0) {
			for (int i = 0, im = ratingScale.size(); i < im; i++) {
				double val = ratingScale.get(i);
				ratingScale.set(i, val + epsilon);
			}
			for (int row = 0; row < numRows; row++) {
				for (int col = 0; col < numCols; col++) {
					if (dataTable.contains(row, col))
						dataTable.put(row, col, dataTable.get(row, col) + epsilon);
				}
			}
		}

		rateMatrix = new SparseMatrix(numRows, numCols, dataTable, colMap);
		
//		timeMatrix = new SparseMatrix(numRows, numCols, timeTable, colMap);

		sw.stop();
		
		long loadTime = sw.elapsed(TimeUnit.SECONDS);

		System.out.println("maxIds: " + maxIds + " / "
				+ "users: " + numRows + " / "
				+ "items: " + numCols + " / "
				+ "ratings: " + dataTable.size() + " / "
				+ "time: " + loadTime
//				+ "ratings/sec: " + dataTable.size() / loadTime 
				);
		
		dataTable = null;
		timeTable = null;

		return rateMatrix;

	}


}
