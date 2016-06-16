package cdTest;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.apache.commons.dbutils.DbUtils;
import org.la4j.matrix.Matrix;

import com.google.common.base.Stopwatch;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.Table;

import greedyFiltering.GFSorting;
import greedyFiltering.GFSparseMatrix;
import greedyFiltering.GreedyFiltering;
import i5.las2peer.services.ocd.algorithms.ClizzAlgorithm;
import i5.las2peer.services.ocd.algorithms.MergingOfOverlappingCommunitiesAlgorithm;
import i5.las2peer.services.ocd.algorithms.RandomWalkLabelPropagationAlgorithm;
import i5.las2peer.services.ocd.algorithms.SpeakerListenerLabelPropagationAlgorithm;
import i5.las2peer.services.ocd.algorithms.SskAlgorithm;
import i5.las2peer.services.ocd.algorithms.WeightedLinkCommunitiesAlgorithm;
import i5.las2peer.services.ocd.algorithms.utils.OcdAlgorithmException;
import i5.las2peer.services.ocd.graphs.Cover;
import i5.las2peer.services.ocd.graphs.CustomGraph;
import librec.data.MatrixEntry;
import librec.data.SparseMatrix;
import librec.data.VectorEntry;
import y.base.Edge;
import y.base.Node;

public class GraphOCDLib {

	public static void main(String[] args) {
		SparseMatrix ratings = null;
		
		boolean checkDeterministic = false;
		
//		int k = 10;
		int mu = 300;
		
//		int[] ks = new int[] {2,5,10,20,100};
		int[] ks = new int[] {2};

//		String[] ocdAlgorithms = new String[] {"SLPA","CLIZZ","DMID","LC","MONC","SSK","DOCA",""};
//		String[] ocdAlgorithms = new String[] {"SLPA"};
//		String[] ocdAlgorithms = new String[] {"CLIZZ"};
		String[] ocdAlgorithms = new String[] {"DMID"};
//		String[] ocdAlgorithms = new String[] {"LC"};
//		String[] ocdAlgorithms = new String[] {"MONC"};
//		String[] ocdAlgorithms = new String[] {"SSK"};
//		String[] ocdAlgorithms = new String[] {"DOCA"};
		
		int[] maxIdsArray = new int[] {80000, 100000, 120000};
//		int[] maxIdsArray = new int[] {100,1000,4000,8000,12000,20000,40000,60000,80000,100000,120000};

		for (int maxIds : maxIdsArray){
			System.out.println("======");
			try {

				// Get ratings matrix from database
				ratings = loadData(maxIds);
				
			} catch (Exception e1) {
				e1.printStackTrace();
//				System.out.println("Exception: " + e1.getMessage());
			}
			
			for (int k : ks){
				// Build user graph from ratings matrix, store in a sparse adjacency matrix
				SparseMatrix userMatrix = getUserMatrix(ratings, k, mu);
				
				// Construct CustomGraph from user adjacency matrix
				CustomGraph graph = buildGraph(userMatrix);
				
				for (String algorithm : ocdAlgorithms){
					if (checkDeterministic){
						Cover cover1 = getCover(graph, algorithm);
						Cover cover2 = getCover(graph, algorithm);
						Cover cover3 = getCover(graph, algorithm);
						if (compareCovers(cover1, cover2) == 0 && compareCovers(cover2, cover3) == 0){
							System.out.println("Community structures are identical.");
						}
						else{
							System.out.println("Community structures are different.");
						}
					}
					else{
						try{
							getCover(graph, algorithm);
						}
						catch (Exception e){
							System.out.println(e.getMessage());
						}
					}
				}
			}
		}
	}
	
	private static int compareCovers(Cover cover1, Cover cover2) {
		int numCommunities1 = cover1.communityCount();
		int numCommunities2 = cover2.communityCount();
		if (numCommunities1 != numCommunities2){
			System.out.println("Different number of communities: " + numCommunities1 + "/" + numCommunities2);
			return 1;
		}
		Matrix memberships1 = cover1.getMemberships();
		Matrix memberships2 = cover2.getMemberships();
		int numUsers = memberships1.rows();
		for (int user = 0; user < numUsers; user++){
			for (int community = 0; community < numCommunities1; community++){
				double membership1 = memberships1.get(user, community);
				double membership2 = memberships2.get(user, community);
				if (membership1 != membership2){
					System.out.println("Community membership for user " + user + " and community " + community + " different.");
					return 1;
				}
			}
		}
		return 0;
	}

	private static Cover getCover(CustomGraph graph, String algorithm) {
		Stopwatch sw = Stopwatch.createStarted();
		
		Cover cover = null;
		
		int numCommunities = 0;
		
		Map<String, String> parameters = new HashMap<String, String>();

		try {
			switch (algorithm){
			case "SLPA":
				/* Parameters and default values
				 *  memorySize				100
				 *  probabilityThreshold	0.15
				 */
				SpeakerListenerLabelPropagationAlgorithm slpaAlgo = new SpeakerListenerLabelPropagationAlgorithm();
				cover = slpaAlgo.detectOverlappingCommunities(graph);
				break;
			case "CLIZZ":
				/* Parameters and default values
				 *  membershipsIterationBound	1000
				 *  membershipsPrecisionFactor	0.001
				 *  influenceFactor				0.9
				 */
				ClizzAlgorithm clizzAlgo = new ClizzAlgorithm();
				
				parameters.put("membershipsIterationBound", "1000");
				parameters.put("membershipsPrecisionFactor", "0.001");
				parameters.put("influenceFactor", "1000");
				clizzAlgo.setParameters(parameters);
				
				cover = clizzAlgo.detectOverlappingCommunities(graph);
				break;
			case "DMID":
				/* Parameters and default values
				 *  leadershipIterationBound	1000
				 *	leadershipPrecisionFactor	0.001
				 *	profitabilityDelta			0.1
				 */
				RandomWalkLabelPropagationAlgorithm dmidAlgo = new RandomWalkLabelPropagationAlgorithm();
				cover = dmidAlgo.detectOverlappingCommunities(graph);
				break;
			case "LC":
				/* No parameters
				 */
				WeightedLinkCommunitiesAlgorithm lcAlgo = new WeightedLinkCommunitiesAlgorithm();
				cover = lcAlgo.detectOverlappingCommunities(graph);
				break;
			case "MONC":
				/* No parameters
				 */
				MergingOfOverlappingCommunitiesAlgorithm moncAlgo = new MergingOfOverlappingCommunitiesAlgorithm();
				cover = moncAlgo.detectOverlappingCommunities(graph);
				break;
			case "SSK":
				/* Parameters and default values
				 *  leadershipIterationBound	1000
				 *  membershipsIterationBound	1000
				 *  leadershipPrecisionFactor	0.001
				 *  membershipsPrecisionFactor	0.001
				 */
				SskAlgorithm sskAlgo = new SskAlgorithm();
				
				parameters.put("leadershipIterationBound", "50");
				parameters.put("membershipsIterationBound", "50");
				parameters.put("leadershipPrecisionFactor", "0.1");
				parameters.put("membershipsPrecisionFactor", "0.1");
				sskAlgo.setParameters(parameters);

				cover = sskAlgo.detectOverlappingCommunities(graph);
				break;
			}
		} catch (InterruptedException | OcdAlgorithmException e) {
			e.printStackTrace();
		}
		
		sw.stop();
		
		long computeTime = sw.elapsed(TimeUnit.SECONDS);
		numCommunities = cover.communityCount();
		double avgCommunitySize = 0;
		for (int community = 0; community < numCommunities; community++){
			avgCommunitySize += (double) cover.getCommunitySize(community) / (double) numCommunities;
		}
		
		System.out.println("Algorithm: " + algorithm
				+ " / NumCommunities: "	+ numCommunities
				+ " / Avg community size: "	+ avgCommunitySize
				+ " / Time: " + computeTime);

		return cover;
	}

	private static CustomGraph buildGraph(SparseMatrix matrix){
		CustomGraph graph = new CustomGraph();
		
		BiMap<Integer, Node> userNodeMap = HashBiMap.create();
		
		int numUsers = matrix.numRows();
		
		for (int user = 0; user < numUsers; user++){
			Node node = graph.createNode();
			userNodeMap.put(user, node);
		}
		
		for (MatrixEntry e : matrix){
			int user1 = e.row();
			int user2 = e.column();
			double weight = e.get();
			Node node1 = userNodeMap.get(user1);
			Node node2 = userNodeMap.get(user2);
			
			Edge edge = graph.createEdge(node1, node2);
			graph.setEdgeWeight(edge, weight);
		}
		
		return graph;
	}
	
	private static SparseMatrix getUserMatrix(SparseMatrix ratings, int k, int mu){
		
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
