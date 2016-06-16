package greedyFiltering;

import java.io.*;
import java.util.*;

/*
 * This program is "GF Single Machine" 
 * which is used to generate an approximate k-NN graph on a single machine.
 * If you want to construct a k-NN graph using multiple machines,
 * we recommend you to use "GF MapReduce"
 *    
 * If you want to use GF in your research,
 * please cite the following paper:
 * 
 * Greedy Filtering: "A Scalable Algorithm for K-Nearest Neighbor Graph Construction" 
 * (by Youngki Park, Sungchan Park, Woosung Jung, and Sang-goo Lee). 
 * In Proceedings of the 19th International Conference on 
 * Database Systems for Advanced Applications (DASFAA'14), 
 * pp. 327-341, 2014. 
 */

/*
 * The process of this program:
 * loop
 *   1) Read one line from config.txt
 *   2) Read the sparse (input) matrix
 *   3) Construct a k-NN graph
 *   4) Write results to an output file
 *   5) Print the (estimated) accuracy of the constructed k-NN graph
 * end
 */

/*
 * The format of config.txt:
 * inputFileName outputFileName numVector k mu
 * - The first "numVector" lines of "inputFileName" is read
 * - The constructed k-NN graph is written into "outputFileName"
 * - k is used to construct a k-NN graph
 * - mu is be used to control the accuracy level of greedy filtering
 *   (please refer to the abovementioned paper)
 */

public class Main {	
	public static void main(String args[]) {		
		try {
			BufferedReader reader = new BufferedReader(new FileReader("config.txt"));
			String line;
			
			while ((line = reader.readLine()) != null) {
				// 1) Read one line from config.txt
				StringTokenizer tokenizer = new StringTokenizer(line, " ");
				String inputFileName = tokenizer.nextToken();
				String outputFileName = tokenizer.nextToken();
				int numVector = Integer.parseInt(tokenizer.nextToken());
				int k = Integer.parseInt(tokenizer.nextToken());
				int mu = Integer.parseInt(tokenizer.nextToken());
				int numElement = IOProcessing.getNumElement(inputFileName);
				int numVerification;				
				
				// 2) Read the sparse (input) matrix
				GFSparseMatrix inputMatrix = new GFSparseMatrix();
				inputMatrix.numVector = numVector;
				inputMatrix.vectorIndex = new int[numVector+1];
				inputMatrix.dimension = new int[numElement];
				inputMatrix.value = new double[numElement];
				inputMatrix.numDimension = IOProcessing.readSparseMatrix(inputFileName, inputMatrix);
				
				System.out.println("Input: " + inputFileName);
				System.out.println("Output: " + outputFileName);
				System.out.println("# Vectors: " + numVector);
				System.out.println("# Dimensions: " + inputMatrix.numDimension);
				System.out.println("# Elements: " + numElement);				
				System.out.println("k = " + k + ", mu: " + mu);
				
				// 3) Construct a k-NN graph
				double[][] kNNGraph = new double[numVector][];
				for (int i=0; i<numVector; i++) {
					kNNGraph[i] = new double[k*2];
					for (int j=0; j<k*2; j++) {
						kNNGraph[i][j] = -1;
					}
				}				
				GreedyFiltering.kNNGraphConstruction(inputMatrix, kNNGraph, k, mu);
				
				// 4) Write results to "outputFileName"
				IOProcessing.writekNNGraph(outputFileName, kNNGraph, k, numVector);
								
				// 5) Print the (estimated) accuracy of the constructed k-NN graph				
				if (numVector < 1000) numVerification = numVector;				
				else numVerification = 1000;
				AccuracyEstimation.printAccuracy(inputMatrix, kNNGraph, k, numVerification);
			}			
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}