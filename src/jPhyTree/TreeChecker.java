package jPhyTree;

import java.io.*;
import java.lang.*;
import java.util.*;

import edu.uci.ics.jung.graph.util.*;

/**
 * Class: TreeChecker
 * Constructor: None
 * Last Edited: September 13, 2012
 * ----
 * This class embodies the majority of the algorithmic work in regards to building
 * the tree.
 *
 */
public class TreeChecker {

	private static final double THRESHOLD_PVALUE = VCFConstants.THRESHOLD_PVALUE;
	
	/* Private Instance Variables */
	private static int numRows;
	private static int numCols;
	private static Map<Integer, Double> rowToMutRateMap;
	private static Map<String, Double> codeToMutRateMap;	
	private static Set<Integer> minCover;
	private static Double min;
	
	/*===============================================================
	 * Public Methods
	 */
	
	/**
	 * Function: main(String[] args)
	 * Usage: (Main Method)
	 * ----
	 * Wrapper for running the main program. What is called
	 * when TreeChecker is called from command line.
	 * 
	 * @param args	Input arguments from command line
	 */
	public static void main(String[] args) {
		if (args.length == 0) {
			/**
			 * Taken out for time testing
			 */
//			System.out.println("Must pass in an input matrix file!");
		}
		else {
			ArrayList<ArrayList<Integer>> matrixPrime = checkIfTree(args[0]);
			/**
			 * Taken out for time testing
			 */
//			if (matrixPrime != null) System.out.println("This can be a PhyTree!");
//			else System.out.println("This cannot be a PhyTree!");
		}
	}

	/**
	 * Function: checkIfTree(String matrixFile)
	 * Usage: ArrayList<ArrayList<Integer>> matrixPrime = TreeChecker.checkIfTree(matrixFile)
	 * ----
	 * This method checks if the given matrix file can be turned into
	 * a phylogenetic tree.
	 * 
	 * This method utilizes several static methods of the class to
	 * generate a new tree. It creates an instance of the class
	 * and refers to it for all operations. Matrices are stored
	 * as 2D ArrayLists of Integer objects.
	 * 
	 * @param	matrixFile	A file consisting of 0s and 1s which define the input matrix
	 * @return 				true if file can be made into a matrix, false otherwise 
	 */
	public static ArrayList<ArrayList<Integer>> checkIfTree(String matrixFile) {
		TreeChecker treeCheckInstance = new TreeChecker();
		ArrayList<ArrayList<Integer>> matrix = treeCheckInstance.processMatrixFile(matrixFile);
		//Can probably just put the static method inside one another
		ArrayList<ArrayList<Integer>> transMatrix = transposeMatrix(matrix);
		ArrayList<Integer> binaryCodeList = generateBinaryCodeList(transMatrix);
		HashMap<Integer, ArrayList<Integer>> codeToColumnHash = treeCheckInstance.getHashFromTransMatrix(transMatrix);
		ArrayList<ArrayList<Integer>> matrixPrime = treeCheckInstance.generateMatrixPrime(binaryCodeList, codeToColumnHash, transMatrix);
		/**
		 * Taken out for time testing
		 */
//		printMatrix(matrixPrime);
		HashMap<Pair<Integer>, Integer> LFunctionMap = treeCheckInstance.generateLFunctionMap(matrixPrime);
		HashMap<Integer, Integer> LColFuncMap = treeCheckInstance.generateLColFuncMap(LFunctionMap, matrixPrime);
		boolean isTree = treeCheckInstance.isPhyTree(LFunctionMap, LColFuncMap);
		//If is false, print conflicts
		//if (!isTree) findConflicts(matrixPrime);
		if (isTree) return matrixPrime;
		else return null;
	}
	
	/**
	 * Function: getConflicts(ArrayList<ArrayList<Integer>> matrixPrime,
	 * 			 ArrayList<ArrayList<Integer>> noConflictMatrixPrime)
	 * Usage: TreeChecker.getConflicts(matrixPrime, noConflictMatrixPrime)
	 * ----
	 * Takes in matrixPrime and noConflictMatrixPrime and finds the differences.
	 * These binary codes are conflicts, which are returned in a set.
	 * 
	 * @param matrixPrime
	 * @param noConflictMatrixPrime matrixPrime with conflicts removed
	 * @return A set of conflicting binary codes
	 */
	public static Set<ArrayList<Integer>> getConflicts(ArrayList<ArrayList<Integer>> matrixPrime, 
			ArrayList<ArrayList<Integer>> noConflictMatrixPrime){
		Set<ArrayList<Integer>> allCodes = new HashSet<ArrayList<Integer>>(transposeMatrix(matrixPrime));
		//System.out.println(allCodes.toString());
		Set<ArrayList<Integer>> noConflictCodes = new HashSet<ArrayList<Integer>>(transposeMatrix(noConflictMatrixPrime));
		allCodes.removeAll(noConflictCodes);
		Set<ArrayList<Integer>> conflicts = allCodes;
		/**
		 * Taken out for time testing
		 */
//		System.out.println("Final Matrices: ");
//		TreeChecker.printMatrix(matrixPrime);
//		System.out.println("-----");
//		TreeChecker.printMatrix(noConflictMatrixPrime);
		//System.out.println("-----");
		//System.out.println(conflicts.toString());
		return conflicts;
	}
	
	/**
	 * Function: findConflicts(ArrayList<ArrayList<Integer>> matrixPrime)
	 * Usage: ArrayList<Integer> conflicts = findConflicts(matrixPrime)
	 * ----
	 * Returns a list of columns which are in conflict
	 * @param matrixPrime	
	 * @return	An ArrayList of conflicting columns in matrixPrime
	 */
	public static ArrayList<Integer> findConflicts(ArrayList<ArrayList<Integer>> matrixPrime) {
		ArrayList<ArrayList<Integer>> conflictMatrix = new ArrayList<ArrayList<Integer>>();
		//Makes the conflict matrix as according to Gusfield's algorithm
		for (int i = 0; i < matrixPrime.size(); i++){
			conflictMatrix.add(new ArrayList<Integer>());
			int counter = 0;
			for (int j = 0; j < matrixPrime.get(i).size(); j++){
				if (matrixPrime.get(i).get(j) == 1){
					conflictMatrix.get(i).add(counter);
					counter = j + 1;
				} else conflictMatrix.get(i).add(-1);
			}
		}
		/**
		 * Taken out for time testing
		 */
//		System.out.println("ConflictMatrix:");
//		printMatrix(conflictMatrix);
//		System.out.println("--");
		ArrayList<Integer> conflictCols = new ArrayList<Integer>();
		/*
		 * Vertex Cover Initialization - Creating Conflict Graph
		 * 1. Create Nodes
		 * 2. Create Edges by finding conflicts
		 */
		Set<Integer> nodes = new HashSet<Integer>();
		Set<Pair<Integer>> edges = new HashSet<Pair<Integer>>(); 
		ArrayList<Integer> bCodeList = generateBinaryCodeList(transposeMatrix(matrixPrime));
		for (int i = 0; i < numRows; i++) nodes.add(i);
		for (int j = 0; j < conflictMatrix.get(0).size(); j++){
//			if (colHasConflicts(conflictMatrix, j)) {
				edges.addAll(findEdges(matrixPrime, j, bCodeList));
//			}
		}
		/**
		 * Taken out for time testing
		 */
//		System.out.println("Nodes: " + nodes.toString());
//		System.out.println("Edges: " + edges.toString());
		Set<Integer> vertexCover = new HashSet<Integer>();
		//Set<Integer> minCover = null;
//		for (int k = 1; k < nodes.size(); k++){
//			vertexCover = recVertexCover(edges, k, new HashSet<Integer>());
//			if (vertexCover != null) break;
//		}
		minCover = null;
		min = -1.00;
		//long startTime = System.nanoTime();	
		//System.out.println("Start: " + Long.toString(startTime));
		//recVertexCover(edges, nodes, vertexCover);
		//Do approximate vertex cover to find which conflicts to remove
		minCover = approxVertexCover(edges, nodes, vertexCover);
		//long endTime = System.nanoTime();
		//System.out.println("End: " + Long.toString(endTime));
		//System.out.println("Took "+(endTime - startTime) + " ns for vertex cover");
		//Add the result of approximate vertex cover
		conflictCols.addAll(minCover);
		//printMatrix(conflictMatrix);
		Set<Integer> nonConflictCols = getNonConflictCols(bCodeList);
		conflictCols.removeAll(nonConflictCols);
		/**
		 * Taken out for time testing
		 */
//		System.out.print("The following columns had conflicts: ");
//		for (int i = 0; i < conflictCols.size(); i++){
//			System.out.print(conflictCols.get(i) + " ");
//		}
//		System.out.println();
		return conflictCols;
	}
	
	private static Set<Integer> getNonConflictCols(ArrayList<Integer> bCodeList) {
		Set<Integer> nonConflictCols = new HashSet<Integer>();
		for (int i = 0; i < bCodeList.size(); i++){
			String codeStr = Integer.toBinaryString(bCodeList.get(i));
			int counter = 0;
			for (int j = 0; j < codeStr.length(); j++){
				if (codeStr.charAt(j) == '1') counter++;
			}
			if (counter == 1) nonConflictCols.add(i);
		}
		return nonConflictCols;
	}

	/**
	 * Function: getCFMatrixPrime(String matrixFile)
	 * Usage: ArrayList<ArrayList<Integer>> CFMatrixPrime = TreeChecker.getCFMatrixPrime(matrixFile)
	 * ----
	 * Takes the file with the matrix and returns matrixPrime without the conflicts.
	 * @param matrixFile	The pathname to the text file containing the matrix
	 * @return	matrixPrime with conflicts removed
	 */
	public static ArrayList<ArrayList<Integer>> getCFMatrixPrime(String matrixFile){
		ArrayList<ArrayList<Integer>> matrixPrime = getMatrixPrime(matrixFile);
		ArrayList<ArrayList<Integer>> matrixPrimeTrans = transposeMatrix(matrixPrime);
		ArrayList<Integer> conflicts = findConflicts(matrixPrime);
		Collections.sort(conflicts);
		int counter = 0;
		for (Integer conflict: conflicts){
			//System.out.println("Removing code for column " + conflict.toString() + ": ");
			ArrayList<Integer> curr = matrixPrimeTrans.get(conflict.intValue() - counter);
			//for (int i = 0; i < curr.size(); i++) System.out.print(curr.get(i));
			//System.out.println();
			matrixPrimeTrans.remove(conflict.intValue() - counter);
			counter++;
		}
		//printMatrix(matrixPrimeTrans);
		return transposeMatrix(matrixPrimeTrans);
	}
	
	/**
	 * Function: getMatrixPrime(String matrixFile)
	 * Usage: ArrayList<ArrayList<Integer>> matrixPrime = TreeCheckInstance.getMatrixPrime(matrixFile)
	 * ----
	 * Public static method which will return M' given an input matrix
	 * 
	 * Makes a new instance of TreeChecker and runs through necessary steps
	 * to get M'.
	 * 
	 * @param matrixFile	Input matrix file
	 * @return				matrixPrime (M')
	 */
	public static ArrayList<ArrayList<Integer>> getMatrixPrime(String matrixFile){
		TreeChecker treeCheckInstance = new TreeChecker();
		ArrayList<ArrayList<Integer>> matrix = treeCheckInstance.processMatrixFile(matrixFile);
		ArrayList<ArrayList<Integer>> transMatrix = treeCheckInstance.transposeMatrix(matrix);
		ArrayList<Integer> binaryCodeList = treeCheckInstance.generateBinaryCodeList(transMatrix);
		HashMap<Integer, ArrayList<Integer>> codeToColumnHash = treeCheckInstance.getHashFromTransMatrix(transMatrix);
		return treeCheckInstance.generateMatrixPrime(binaryCodeList, codeToColumnHash, transMatrix);
	}
	
	/**
	 * Function: getMutMap()
	 * Usage: Map<String, Double> mutMap = TreeChecker.getMutMap()
	 * ----
	 * Returns a map of codes to number of mutations for the code.
	 * @return	A map of GATK strings to number of mutations for that code.
	 */
	public static Map<String, Double> getMutMap(){
		return new HashMap<String, Double>(codeToMutRateMap);
	}
	
	/**
	 * Function: getNumTotalEdges()
	 * Usage: int numEdges = getNumTotalEdges()
	 * ----
	 * Returns the number of edges in the matrix. This is equivalent to 
	 * the number of binary codes in the map.
	 * @return Number of edges in the matrix
	 */
	private static int getNumTotalEdges(){
		//return getMutMap().keySet().size();
		return numCols+1;
	}
	
	private static double getTotalMutations(){
		Map<String, Double> mutMap = getMutMap();
		double counter = 0.0;
		for (Double value: mutMap.values()){
			counter += value;
		}
		return counter;
	}
	
	public static int getConflictThreshold(int n, int k){
		double pValue = 0.0;
		int x = n / k;
		for (; x > 1; x--){
			pValue = Math.pow((n - x * k)/(double) n, k - 1);
			if (pValue > THRESHOLD_PVALUE) break;
		}
		return x;
	}
	
	/**
	 * Function: approxVertexCover(Set<Pair<Integer>> edges, Set<Integer> nodes, Set<Integer> currCover)
	 * Usage: Set<Integer> vertexCover = approxVertexCover(edges, nodes, currCover)
	 * ----
	 * @param edges	A set of pairs of integers. Each pair represents an edge.
	 * @param nodes	A set of integers
	 * @param currCover	The set of nodes which make the vertex cover
	 * @return	The vertex cover
	 */
	private static Set<Integer> approxVertexCover(Set<Pair<Integer>> edges, Set<Integer> nodes, Set<Integer> currCover){
		/*
		 * 1. Sort mutation groups by size.
		 * 2. Find largest mutation (that hasn't been removed).
		 * -> Find all edges with this group.
		 * -> Store connecting nodes in currCover
		 * -> Remove the node's connecting edges and corresponding nodes
		 * 3. Repeat 2 until no edges left 
		 */
		//Initialize edges and nodes
		Set<Pair<Integer>> allEdges = new HashSet<Pair<Integer>>(edges);
		Set<Integer> allNodes = new HashSet<Integer>(nodes);
		
		int numEdgesInTree = getNumTotalEdges();
		int totalMut = (int) getTotalMutations();
		/**
		 * Taken out for time testing
		 */
//		System.out.println("Total Mutations: " + totalMut);
		int conflictThreshold = getConflictThreshold(totalMut, numEdgesInTree);
		//System.out.println("Conflict Threshold: " + conflictThreshold);
		//while there are still conflicting edges
		while (!allEdges.isEmpty()){
			//find largest node
			int maxNode = -1;
			double maxMutRate = 0.0;
			for (Integer node: allNodes){
				double currMutRate = rowToMutRateMap.get(node+1);
				if (currMutRate > maxMutRate){
					maxMutRate = currMutRate;
					maxNode = node;
				}
			}
			/* recalculate threshold */
			if (maxNode == 113){
			}
			conflictThreshold = getConflictThreshold(totalMut, numEdgesInTree);
			/**
			 * Taken out for time testing
			 */
//			System.out.println("Conflict Threshold: " + conflictThreshold);
			/*
			 * if the highest mut rate is less than
			 * conflictThreshold, no more nodes
			 * can be added to maximal independent set,
			 * meaning remaining nodes should be put
			 * in vertex cover
			 */
			if (maxMutRate < conflictThreshold) {
				currCover.addAll(allNodes);
				break;
			}
			//remove largest node
			allNodes.remove(maxNode);
			numEdgesInTree++;
			//find all nodes which conflict with largest
			for (Pair<Integer> edge: allEdges){
				if (edge.getFirst().equals(maxNode)){
					int conflictNode = edge.getSecond().intValue();
					currCover.add(conflictNode);
					//allEdges.remove(edge);
				} else if (edge.getSecond().equals(maxNode)){
					int conflictNode = edge.getFirst().intValue();
					currCover.add(conflictNode);
					//allEdges.remove(edge);
				}
			}
			//System.out.println("Edges left after initial trimming: " + allEdges.toString());
			//Remove these nodes
			allNodes.removeAll(currCover);
			//System.out.println("Remaining Nodes: " + allNodes.toString());
			//Remove any edges 
			Set<Pair<Integer>> edgesToRemove = new HashSet<Pair<Integer>>();
			for (Pair<Integer> edge: allEdges){
				int first = edge.getFirst().intValue();
				int second = edge.getSecond().intValue();
				if (currCover.contains(first) || currCover.contains(second)){
					//allEdges.remove(edge);
					edgesToRemove.add(edge);
				}
			}
			allEdges.removeAll(edgesToRemove);
			//System.out.println("Edges left after final cuts: " + allEdges.toString());
		}
		return currCover;
	}

//	private static void recVertexCover(Set<Pair<Integer>> edges, Set<Integer> nodes, Set<Integer> currCover){
//		boolean isVertexCover = vertexCoverCheck(new HashSet<Pair<Integer>>(edges), currCover);
//		double weight = coverWeight(currCover);
//		if (isVertexCover && (weight < min || min < 0)){
//			minCover = new HashSet<Integer>(currCover);
//			min = weight;
//		} else {
//			for (Integer node: nodes){
//				currCover.add(node);
//				Set<Integer> newNodes = new HashSet<Integer>(nodes);
//				newNodes.remove(node);
//				recVertexCover(edges, newNodes, currCover);
//				currCover.remove(node);
//			}
//		}
//		return;
//	}
	
//	private static double coverWeight(Set<Integer> currCover) {
//		double sum = 0;
//		for (Integer node: currCover) sum += rowToMutRateMap.get(node+1);
//		return sum;
//	}

//	private static boolean vertexCoverCheck(Set<Pair<Integer>> edges,
//			Set<Integer> currCover) {
//		for (Pair<Integer> edge: edges){
//			int first = edge.getFirst();
//			int second = edge.getSecond();
//			if (!currCover.contains(first) && !currCover.contains(second)) return false;
//		}
//		return true;
//	}

//	private static Set<Integer> recVertexCover(Set<Pair<Integer>> edges,
//			int k, Set<Integer> currVertexCover) {
//		// TODO Auto-generated method stub
//	/*
//	 * 1. For i = 1 to numNodes-1, do VertexCover
//	 * 2. if edges is empty, return currentSet of nodes
//	 * 3. if k = 0; didn't find one in this branch
//	 * 4. Else, add endpoint with lowest mutation rate
//	 * 5. Let newEdges be all of the edges in edges without newest endpoint
//	 * 6. if VertexCover(edges', k-1, currNodeSet).notEmpty; found one
//	 * 7. Otherwise, didn't find one in this branch.
//	 */
//		if (edges.isEmpty()) return currVertexCover;
//		else if (k == 0) return null;
//		else{
////			Set<Integer> nodesInEdges = new HashSet<Integer>();
////			for (Pair<Integer> edge: edges){
////				nodesInEdsges.add(edge.getFirst());
////				nodesInEdges.add(edge.getSecond());
////			}
//			ArrayList<Pair<Integer>> edgesList = new ArrayList<Pair<Integer>>();
//			edgesList.addAll(edges);
//			Pair<Integer> currEdge = edgesList.get(0);
//			//Need to +1 because map is indexed starting at 1s
//			Integer endpointX = (rowToMutRateMap.get(currEdge.getFirst() + 1) <= rowToMutRateMap.get(currEdge.getSecond() + 1) ? 
//					currEdge.getFirst() : currEdge.getSecond());
//			currVertexCover.add(endpointX);
//			Set<Pair<Integer>> edgesPrime = new HashSet<Pair<Integer>>();
//			for (Pair<Integer> edge: edges){
//				if (edge.getFirst() != endpointX && edge.getSecond() != endpointX)
//					edgesPrime.add(edge);	
//			}
//			Set<Integer> result = recVertexCover(edgesPrime, k-1, currVertexCover);
//			if (result != null) return result;
//			else {
//				currVertexCover.remove(endpointX);
//				Integer endpointY = endpointX == currEdge.getFirst() ? currEdge.getSecond(): currEdge.getFirst();
//				currVertexCover.add(endpointY);
//				edgesPrime = new HashSet<Pair<Integer>>();
//				for (Pair<Integer> edge: edges){
//					if (edge.getFirst() != endpointY && edge.getSecond() != endpointY)
//						edgesPrime.add(edge);	
//				}
//				return recVertexCover(edgesPrime, k-1, currVertexCover);
//			}
//		}
//	}

	
	/**
	 * Function: findEdges(ArrayList<ArrayList<Integer>> matrixPrime, int j, ArrayList<Integer> bCodeList)
	 * Usage: ArrayList<Pair<Integer>> edges = findEdges(matrixPrime, j, bCodeList)
	 * ----
	 * Given a specific column, checks all previous columns in matrixPrime to see if the pair conflict.
	 * Returns a list of all edges which conflict with the given column.
	 *
	 * @param matrixPrime
	 * @param j	The current column to check for conflicts
	 * @param bCodeList	The list of binary codes from matrixPrime
	 * @return	A list of conflict edges
	 */
	private static ArrayList<Pair<Integer>> findEdges(
			ArrayList<ArrayList<Integer>> matrixPrime, int j, ArrayList<Integer> bCodeList) {
		ArrayList<Pair<Integer>> newEdges = new ArrayList<Pair<Integer>>();
		int bCode1 = bCodeList.get(j);
		for (int i = 0; i < j; i++){
			int bCode2 = bCodeList.get(i);
			//if not subsets or disjoint, conflict
			if ((bCode1 | bCode2) != bCode2 && (bCode1 & bCode2) != 0){
				newEdges.add(new Pair<Integer>(i, j));
			}
		}
		return newEdges;
	}

	/**
	 * Function: colHasConflicts(ArrayList<ArrayList<Integer>> conflictMatrix, int j)
	 * Usage: boolean hasConflict = colHasConflicts(conflictMatrix, int j)
	 * ----
	 * Given the conflictMatrix and a column, checks whether the column has any conflicts.
	 * 
	 * @param conflictMatrix	A matrix which is marked with any present conflicts
	 * @param j					The column being checked for conflicts
	 * @return					True if the column does have conflicts; else, false
	 */
	private static boolean colHasConflicts(ArrayList<ArrayList<Integer>> conflictMatrix, int j) {
		int colValue = -1;
		for (int i = 0; i < conflictMatrix.size(); i++){
			int currValue = conflictMatrix.get(i).get(j);
			if (currValue != -1 && colValue == -1){
				colValue = currValue;
			} else if (currValue != -1 && colValue != -1){
				if (colValue != currValue) return true;
			}
		}
		return false;
	}

	/**
	 * Function: printMatrix(ArrayList<ArrayList<Integer>> matrix)
	 * Usage: printMatrix(matrix)
	 * ----
	 * Static method for the class that prints out a matrix
	 * 
	 * Prints out to std.out
	 * 
	 * @param matrix	The 2D matrix which is to be printed
	 */
	public static void printMatrix(ArrayList<ArrayList<Integer>> matrix){
		for (int i = 0; i < matrix.size(); i++){
			char index = (char) ('A' + i);
			System.out.print(index);
			System.out.print(":\t "); 
			for (int j = 0; j < matrix.get(i).size(); j++){
				System.out.print(matrix.get(i).get(j));
				if (j != matrix.get(i).size() - 1) System.out.print("\t|");
				else System.out.println();
			}
		}
	}
	
	/**
	 * Class: decreasingOrderComparator
	 * Purpose: Normal comparisons are done in ascending order. To achieve
	 * descending order, made this short comparator.
	 */
	public class decreasingOrderComparator implements Comparator<Integer>{
		public int compare(Integer int1, Integer int2){
			if (int1 > int2) return -1;
			else if (int1 < int2) return 1;
			else return 0;
		}
	}
		
	/*================================================================
	 * Private Methods
	 */
	
	/**
	 * Uses the LFunctionMap and LColFuncMap to determine if a matrix can be built
	 * 
	 * Returns a boolean
	 * 
	 * @param lFunctionMap	HashMap created by generateLFunctionMap
	 * @param lColFuncMap	Hashmap created by generateLColFuncMap
	 * @return				true if input matrix can be a PhyTree; else, false
	 */
	private boolean isPhyTree(
			HashMap<Pair<Integer>, Integer> lFunctionMap,
			HashMap<Integer, Integer> lColFuncMap) {
		for (Pair<Integer> key : lFunctionMap.keySet()){
			int j = key.getSecond();
			if (lFunctionMap.get(key) != lColFuncMap.get(j)) return false;
		}
		return true;
	}

	/**
	 * Takes in the LFunction HashMap and M' and returns the L function for columns.
	 * 
	 * The L Function for column j returns the largest L(i, j) such that M'(i, j) == 1.
	 * 
	 * @param lFunctionMap	The HashMap created by generateLFunctionMap
	 * @param matrixPrime	M'
	 * @return				A HashMap of M' columns to respective L(j) value
	 */
	private HashMap<Integer, Integer> generateLColFuncMap(
			HashMap<Pair<Integer>, Integer> lFunctionMap,
			ArrayList<ArrayList<Integer>> matrixPrime) {
		int numRows = matrixPrime.size();
		int numCols = matrixPrime.get(0).size();
		HashMap<Integer, Integer> colFuncMap = new HashMap<Integer, Integer>();
		for (int j = 0; j < numCols; j++){
			int currLValue = 0;
			for (int i = 0; i < numRows; i++){
				Pair<Integer> currPair = new Pair<Integer>(i, j);
				if (lFunctionMap.containsKey(currPair)){
					Integer lValue = lFunctionMap.get(currPair);
					if (lValue.intValue() > currLValue) currLValue = lValue.intValue();
				}
			}
			colFuncMap.put(j, currLValue);
		}
		return colFuncMap;
	}

	/**
	 * Returns a HashMap which simulates the L function in Gusfield's paper
	 * 
	 * The L function takes an (i,j) coordinate-pair from M' and maps it to
	 * the largest index k < j such that M'(i, k) == 1. If there is no such
	 * index, L(i,j) is set to 0.
	 * 
	 * @param matrixPrime	M'
	 * @return				A HashMap that takes coordinates i and j and returns L(i, j)
	 */
	private HashMap<Pair<Integer>, Integer> generateLFunctionMap(
			ArrayList<ArrayList<Integer>> matrixPrime) {
		HashMap<Pair<Integer>, Integer> funcMap = new HashMap<Pair<Integer>, Integer>();
		int numRows = matrixPrime.size();
		int numCols = matrixPrime.get(0).size();
		for (int j = 0; j < numCols; j++){
			for (int i = 0; i < numRows; i++){
				if (matrixPrime.get(i).get(j) == 1){
					Pair<Integer> newPair = new Pair<Integer>(i, j);
					funcMap.put(newPair, getFuncValue(matrixPrime, i, j));
				} 
			}
		}
		return funcMap;
	}

	/**
	 * Function: getFuncValue(ArrayList<ArrayList<Integer>> matrixPrime, int i, int j)
	 * Usage: int value = getFuncValue(matrixPrime, i, j)
	 * ----
	 * Finds the value of the L-Function for this (i,j) coordinate in matrixPrime
	 * 
	 * @param matrixPrime
	 * @param i	The row
	 * @param j	The col
	 * @return	L(i,j)
	 */
	private int getFuncValue(ArrayList<ArrayList<Integer>> matrixPrime, int i, int j) {
		if (j == 0) return 0;
		for (int currCol = j-1; currCol >= 0; currCol--){
			if (matrixPrime.get(i).get(currCol) == 1) return currCol + 1;
		}
		return 0;
	}
	
	/**
	 * Returns M' from Gusfield 1991 Paper - 
	 * "Efficient Algorithms for Inferring Evolutionary Trees"
	 * 
	 * M' (or MatrixPrime) is a transformed matrix from the original input version
	 * M' has all of its columns sorted in decreasing order by the column's binary 
	 * code. Repeated columns (which would appear right beside each other) have been
	 * deleted. If one wishes to see which codes correspond with which column, use
	 * the code-To-Col HashMap.
	 * 
	 * @param binaryCodeList	The list of binary codes corresponding to the column values
	 * @param codeToColumnHash	The hash of each code to its respective column values
	 * @param tMat				The transposed matrix
	 * @return					M' - a 2D ArrayList of Integers
	 */
	private ArrayList<ArrayList<Integer>> generateMatrixPrime(
			ArrayList<Integer> binaryCodeList,
			HashMap<Integer, ArrayList<Integer>> codeToColumnHash, ArrayList<ArrayList<Integer>> tMat) {
		Comparator<Integer> decreaseOrderComp = new decreasingOrderComparator();
		Set<Integer> tempSet = new TreeSet<Integer>(decreaseOrderComp);
		tempSet.addAll(binaryCodeList);
		ArrayList<Integer> sortedCodes = new ArrayList<Integer>(tempSet);
		return generateMatrixFromSortedCodes(sortedCodes, codeToColumnHash, tMat);
	}

	/**
	 * Function: generateMatrixFromSortedCodes(ArrayList<Integer> sortedCodes, HashMap<Integer, ArrayList<Integer>> codeToColumnHash, ArrayList<ArrayLIst<Integer>> tMat)
	 * Usage: ArrayList<ArrayList<Integer>> matrixPrime = generateMatrixFromSortedCodes(sortedCodes, codeToColumnHash, tMat)
	 * ----
	 * Given the sorted binary codes, a map of which column they line up to, and the transposed matrix,
	 * returns matrixPrime.
	 * @param sortedCodes	The binary codes sorted in a list.
	 * @param codeToColumnHash	A map of binary codes to their respective columns
	 * @param tMat	The input matrix transposed
	 * @return	matrixPrime
	 */
	private ArrayList<ArrayList<Integer>> generateMatrixFromSortedCodes(
			ArrayList<Integer> sortedCodes,
			HashMap<Integer, ArrayList<Integer>> codeToColumnHash,
			ArrayList<ArrayList<Integer>> tMat) {
		ArrayList<ArrayList<Integer>> matP = new ArrayList<ArrayList<Integer>>();
		for (int i = 0; i < sortedCodes.size(); i++){
			ArrayList<Integer> cols = codeToColumnHash.get(sortedCodes.get(i));
			matP.add(tMat.get(cols.get(0)));
		}
		return transposeMatrix(matP);
	}

	/**
	 * Returns a hash of binary codes to a list of their respective column numbers.
	 * 
	 * @param tMat	The transposed matrix
	 * @return		A HashMap of integers to an ArrayList of integers
	 */
//Might want to make this a public static method to see column repeats in other classes
	private HashMap<Integer, ArrayList<Integer>> getHashFromTransMatrix(
			ArrayList<ArrayList<Integer>> tMat){
		HashMap<Integer, ArrayList<Integer>> codeToColHash = 
			new HashMap<Integer, ArrayList<Integer>>();
		ArrayList<Integer> binaryCodeList = generateBinaryCodeList(tMat);
		for (int i = 0; i < binaryCodeList.size(); i++){
			if (codeToColHash.containsKey(binaryCodeList.get(i))){
				ArrayList<Integer> currCols = codeToColHash.get(binaryCodeList.get(i));
				currCols.add(i);
				codeToColHash.put(binaryCodeList.get(i), currCols);
			} else {
				ArrayList<Integer> newColsList = new ArrayList<Integer>();
				newColsList.add(i);
				codeToColHash.put(binaryCodeList.get(i), newColsList);
			}
		}
		return codeToColHash;
	}
	
	/**
	 * Returns an ArrayList of all the transposed matrix rows as binary codes.
	 * 
	 * Each "code" is an integer who's binary representation matches
	 * the corresponding row's values. The codes are aggregated and
	 * returned.
	 * 
	 * @param tMat	The transposed matrix (columns of original are rows of this matrix)
	 * @return		An ArrayList of integers representing binary codes
	 */
	private static ArrayList<Integer> generateBinaryCodeList(ArrayList<ArrayList<Integer>> tMat){
		ArrayList<Integer> bcList = new ArrayList<Integer>();
		for (int i = 0; i < tMat.size(); i++){
			int code = getBinaryCode(tMat, i);
			bcList.add(code);
		}
		return bcList;
	}
	
	/**
	 * Takes in a transposed matrix and row number, returning the
	 * row's values as an integer who has a matching binary representation.
	 * 
	 * The function takes the String representation of the ArrayList,
	 * takes out non-digit characters, and reparses the now fully-formed integer.
	 * 
	 * Example: If a row's representation is 4 - then it has a binary code of
	 * ...100 (with however many leading zeros as needed)
	 * 
	 * Note: "row" refers to row of the transposed matrix. 
	 * This corresponds to a column of the input matrix.
	 * 
	 * @param tMat	The transposed matrix
	 * @param row	The specific row in tMat
	 * @return		The values of the row encoded in an Integer
	 */
	private static int getBinaryCode(ArrayList<ArrayList<Integer>> tMat, int row){
		ArrayList<Integer> currRow = tMat.get(row);
		String rowStr = currRow.toString().replaceAll("[^0-9]", "");
		return Integer.parseInt(rowStr, 2);
	}
	
	/**
	 * Returns the transpose of a matrix.
	 * 
	 * Returns a matrix such that every (i, j)
	 * value in original matrix is replaced by the
	 * value at (j, i) essentially "flipping" the matrix.
	 * 
	 * @param matrix	Input 2D matrix	
	 * @return			Newly transposed matrix
	 */
	private static ArrayList<ArrayList<Integer>> transposeMatrix(ArrayList<ArrayList<Integer>> matrix){
		int numRows = matrix.size();
		int numCols = matrix.get(0).size();
		ArrayList<ArrayList<Integer>> newMat = new ArrayList<ArrayList<Integer>>();
		//Put this in second for-loop to not have to loop twice?
		for(int i = 0; i < numCols; i++) newMat.add(new ArrayList<Integer>());
		for(int i = 0; i < numRows; i++){
			for (int j = 0; j < numCols; j++){
				newMat.get(j).add(matrix.get(i).get(j));
			}
		}
		return newMat;
	}
	
	/**
	 * Takes in an inputFile, parses it, and builds the matrix
	 * 
	 * Uses a BufferedReader to read each line in one-by-one.
	 * Each line is processed into its own ArrayList of Integers,
	 * and then added as a row to the overall 2D ArrayList.
	 * 
	 * @param inputFile	Matrix input file
	 * @return 			A 2D ArrayList which contains the input matrix
	 */
	private ArrayList<ArrayList<Integer>> processMatrixFile(String inputFile){
		ArrayList<ArrayList<Integer>> matrix = new ArrayList<ArrayList<Integer>>();
		rowToMutRateMap = new HashMap<Integer, Double>();
		codeToMutRateMap = new HashMap<String, Double>();
		BufferedReader input = null;
		try{
			input = new BufferedReader(new FileReader(inputFile));
			processFirstLine(input.readLine());
			String currLine;
			int counter = 0;
			while ((currLine = input.readLine()) != null) {
				matrix.add(processLine(currLine, counter));
				counter++;
			}
		} catch (IOException e) {
			System.out.println("Unforeseen IOException!");
		} finally {
			try {
				if (input != null) input.close();
			} catch (IOException e){
				System.out.println("Unforeseen IOException when closing file!");
			}
		}
		//Transposed for the new type of input. CHECK FOR DOUBLE TRANSPOSITIONS.
		matrix = transposeMatrix(matrix);
		//printMatrix(matrix);
		//System.out.println(rowToMutRateMap.toString());
		return matrix;
	}
	
	/**
	 * Function: processFirstLine(String line)
	 * Usage: processFirstLine(line)
	 * ----
	 * Parses the first line from the matrix text file to set the
	 * numRows and numCols
	 * @param line	First line of the matrix text file
	 */
	private void processFirstLine(String line) {
		String[] numsAsStr = line.split(" ");
		numRows = Integer.parseInt(numsAsStr[0]);
		numCols = Integer.parseInt(numsAsStr[1]);
		/**
		 * Taken out for time testing
		 */
//		System.out.println("Rows: " + numRows + " Cols: " + numCols);
	}

//	/**
//	 * Process a single line from input matrix text file
//	 * 
//	 * Takes in a string line of boolean values separated by spaces
//	 * and converts it into an ArrayList of the values as Integers.
//	 * 
//	 * @param currLine	Input of a single string line
//	 * @return			An ArrayList of the input values
//	 */
//	private ArrayList<Integer> processLine(String currLine){
//		ArrayList<Integer> newRow = new ArrayList<Integer>();
//		String[] numsAsStr = currLine.split(" ");
//		for (int i = 0; i < numsAsStr.length; i++) 
//			newRow.add(Integer.valueOf(numsAsStr[i]));
//		return newRow;
//	}
	
	/**
	 * Function: processLine(String currLine, int counter)
	 * Usage: processLine(currLine, counter)
	 * ----
	 * Parses each line of the matrix text file and returns
	 * a row of the matrix
	 */
	private ArrayList<Integer> processLine(String currLine, int counter){
		ArrayList<Integer> newRow = new ArrayList<Integer>();
		for (int i = 0; i < numCols; i++){
			newRow.add(Character.getNumericValue(currLine.charAt(i))); 
		}
		//System.out.println(currLine);
		String[] numsAsStr = currLine.split(" ");
		double mutRate = Double.parseDouble(numsAsStr[1]);
		//System.out.println(mutRate);
		rowToMutRateMap.put(counter + 1, mutRate);
		codeToMutRateMap.put(numsAsStr[0], mutRate);
		return newRow;
	}
}
