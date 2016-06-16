package greedyFiltering;

import java.io.*;
import java.util.*;

/* 
 * Note: we assume that 
 * 1) The values of each vector are weighted by TF-IDF
 * 2) The vaules of each vector are L2-normalized
 * 3) The elements of each vector are sorted in descending order of their values
 * 
 * The format of inputFileName:
 * 
 * [a dummy value]	[a dimension number] [its corresponding value] [a dimension number] [its corresponding value] ...
 * [a dummy value]	[a dimension number] [its corresponding value] [a dimension number] [its corresponding value] ...
 * ...
 * 
 * The range of the dimension numbers should be [0, #dimensions-1]
 */

public class IOProcessing {
	public static int getNumElement(String inputFileName) {
		int numElement = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(inputFileName));		
			String line;
			
			// Update the max dimension			
			while ((line = reader.readLine()) != null) {
				StringTokenizer tokenizer = new StringTokenizer(line, " \t");
				tokenizer.nextToken(); // ignore the first token (a vector number)
				while (tokenizer.hasMoreTokens()) {
					for (int i=0; i<2; i++)
						tokenizer.nextToken();
					numElement++;					
				}				
			}
		} catch(Exception e) {
			e.printStackTrace();
		}
		return numElement;
	}
	
	public static int readSparseMatrix(String inputFileName, GFSparseMatrix inputMatrix) {		
		int max_dimension = -1;
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(inputFileName));
			String line;
									
			// Read the sparse matrix
			int vectorCounter = 0;
			int elementCounter = 0;
						
			while ((line = reader.readLine()) != null) {				
				StringTokenizer tokenizer = new StringTokenizer(line, " \t");
				tokenizer.nextToken(); // ignore the first token (a vector number)

				inputMatrix.vectorIndex[vectorCounter] = elementCounter; 
				
				while (tokenizer.hasMoreTokens()) {
					int dimension = Integer.parseInt(tokenizer.nextToken());
					double value = Double.parseDouble(tokenizer.nextToken());
					
					inputMatrix.dimension[elementCounter] = dimension;
					inputMatrix.value[elementCounter] = value;
					
					if (max_dimension < dimension)
						max_dimension = dimension;
					
					elementCounter++;
				}

				vectorCounter++;
			}
			inputMatrix.vectorIndex[vectorCounter] = elementCounter;
			
			reader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}		
		
		return max_dimension+1;
	}
	
	public static void writekNNGraph(String outputFileName, double[][] kNNGraph, int k, int numVector) {
		System.out.println("Writing a k-NN graph...");
		
		try {
			FileWriter writer = new FileWriter(outputFileName, false);
			String writeLine;
			
			for (int i=0; i<numVector; i++) {
				for (int j=0; j<k*2; j+=2) {
					writeLine = i + ", " + (int)Math.round(kNNGraph[i][j]) + ", " + kNNGraph[i][j+1] + "\r\n"; 
					writer.write(writeLine);
				}
			}		
			
			writer.close();				
		} catch (Exception e) {
			e.printStackTrace();
		}	
	}
}
