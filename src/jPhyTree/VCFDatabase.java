/**
 * Class: VCFDatabase
 * Constructor: VCFDatabase(String TESTFILE)
 * Last Edited: September 13, 2012
 * ----
 * This class is a container for VCFEntries.
 * The constructor takes the path name of the VCF file,
 * constructs a VCFEntry for each record, and stores it
 * internally. One can get back VCFEntries using their
 * GATK codes.
 */
package jPhyTree;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

public class VCFDatabase {
	
	/* Private Instance Variables */
	private static final double THRESHOLD = VCFConstants.THRESHOLD;
	private static final double COVERAGE = VCFConstants.COVERAGE;
	private ArrayList<VCFEntry> database;
	private int germlineCounter; 
	
	/**
	 * Function: VCFDatabase(String TESTFILE)
	 * Usage: (Constructor)
	 * ----
	 * This is the constructor for the VCFDatabase.
	 * @param TESTFILE	The VCF file path with which to build the database
	 */
	public VCFDatabase(String TESTFILE){
		BufferedReader rd;
		database = new ArrayList<VCFEntry>();
		try{
			rd = new BufferedReader(new FileReader(TESTFILE));
			String currLine = rd.readLine();
			while (currLine.substring(0, 1).equals("#")) currLine = rd.readLine();
			int allCounter = 0;
			int validCounter = 0;
			//germlineCounter = 0;
			while (currLine != null){
				VCFEntry entry = new VCFEntry(currLine);
				currLine = rd.readLine();
				allCounter++;
				//Makes sure entries are legitimate.
				boolean isLegitimate = true;
				int numSamples = entry.getNumSamples();
				for (int i = 0; i < numSamples; i++){
					if (entry.getAlleleFreq(i).equals("./.")){
						isLegitimate = false;
						break;
					}
				}
				if (!isLegitimate) continue;
				/* TO FILTER OUT ALL SHARED SAMPLES */
				boolean isGermline = true;
				for (int i = 0; i < numSamples; i++){
					if(!entry.getGenotype(i).equals("0/1") && !entry.getGenotype(i).equals("1/1")){
						isGermline = false;
						break;
					}
				}
				if (isGermline) {
					//germlineCounter++;
					continue;
				}
				/* TO FILTER OUT SNVs with low coverage in all SAMPLES*/
				int totalCoverage = 0;
				for (int i = 0; i < numSamples; i++){
					totalCoverage += entry.getReadDepth(i);
				}
				if ((2.0*totalCoverage)/numSamples <= COVERAGE){
					isLegitimate = false;
				}
				if (!isLegitimate) continue;
				database.add(entry);
				validCounter++;
			}
			/**
			 * Taken out for time testing
			 */
//			System.out.println("There were " + allCounter + " entries in the file.");
//			System.out.println("Of those, " + validCounter + " were valid, and " + (allCounter - validCounter) + " were not.");
			rd.close();
		} catch (IOException e){
			System.out.println("File Reading Error!");
		}
	}
	
	/**
	 * Function: getEntriesByGATK(String inputCode)
	 * Usage: ArrayList<VCFEntry> entries = db.getEntriesByGATK(inputCode)
	 * ----
	 * Searches the database for all VCFEntries with matching GATK codes
	 * and returns an ArrayList of such entries
	 * 
	 * @param inputCode	The GATK code with which to filter the entries
	 * @return	An ArrayList of VCFEntries with a GATK equivalent to inputCode
	 */
	public ArrayList<VCFEntry> getEntriesByGATK(String inputCode){
		ArrayList<VCFEntry> list = new ArrayList<VCFEntry>();
		for (VCFEntry entry: database){
			if (entry.getGATK().equals(inputCode)) list.add(entry);
		}
		return list;
	}
	
	/**
	 * Function: getSortedEntriesByGATK(String inputCode, String destCode)
	 * Usage: ArrayList<VCFEntry> entries = db.getSortedEntriesByGATK(inputCode, destCode)
	 * ----
	 * Returns all VCFEntries with matching GATK codes to inputCode in an ArrayList. 
	 * The entries are sorted by how likely they are to be converted to destCode in 
	 * ascending order.
	 * 
	 * @param inputCode	The GATK code that should match valid entries
	 * @param destCode	The GATK code to convert inputCode into
	 * @param pw 
	 * @return	A sorted ArrayList of VCFEntries
	 */
	public ArrayList<VCFEntry> getSortedEntriesByGATK(String inputCode, String destCode, PrintWriter pw){
		ArrayList<VCFEntry> entryList = getEntriesByGATK(inputCode);
		//index -> if 0 to 1
		Map<Integer, Boolean> codeDifferenceMap = initCodeDiffMap(inputCode, destCode);
		Set<Integer> mismatchIndices = new HashSet<Integer>(codeDifferenceMap.keySet());
		ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap = new ArrayList<MyEntry<Double, ArrayList<VCFEntry>>>();
		for (VCFEntry entry: entryList){
			double total = -1.0;
			for (Integer index: mismatchIndices){
				//boolean is0to1 = codeDifferenceMap.get(index);
				double currProb = entry.getSumProb(index);
				//if (codeDifferenceMap.get(index)) currProb = 1.0 - currProb;
				total = (total < 0.0 ? currProb : total * currProb);
			}
			if (total < 0.0) total = 0.0;
			
			if (entryListContainsKey(probsToEntryMap, total)) {
				ArrayList<VCFEntry> currEntries = entryListGet(probsToEntryMap, total);
				currEntries.add(entry);
				entryListPut(probsToEntryMap, total, currEntries);
				//probsToEntryMap.put(total, entry);
			} else {
				ArrayList<VCFEntry> newEntries = new ArrayList<VCFEntry>();
				newEntries.add(entry);
				entryListPut(probsToEntryMap, total, newEntries);
			}
		}
		ArrayList<Double> probs = new ArrayList<Double>(new HashSet<Double>(entryListKeySet(probsToEntryMap)));
		Collections.sort(probs);
		//String divider = "-----\n";
		//String outputStr = "Input: " + inputCode + " Dest: " + destCode + "\n";
		//System.out.println(outputStr);
		//pw.write(divider);
		//pw.write(outputStr);
		//pw.write(divider);
		ArrayList<VCFEntry> sortedList = new ArrayList<VCFEntry>();
		for (int i = 0; i < probs.size(); i++){
			//sortedList.add(probsToEntryMap.get(probs.get(i)));
			ArrayList<VCFEntry> newEntryList = entryListGet(probsToEntryMap, probs.get(i));
			if (newEntryList.size() > 1){
				newEntryList.size();
			}
			for(int j = 0; j < newEntryList.size(); j++){
				sortedList.add(newEntryList.get(j));
				//System.out.println(newEntryList.get(j).toString());
			}
//			if (i >= 25){
//				System.out.println("Entries: " + i);
//			}
			//sortedList.addAll();
			//String probLine = i+1 + ". " + "Prob: " + probs.get(i) + "\n";
			//pw.write(probLine);
			//System.out.println(probLine);
			//VCFEntry entry = probsToEntryMap.get(probs.get(i));
			//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + "\n";
			//pw.write(entryLine);
			//System.out.println(entryLine);
		}
		//pw.close();
		//System.out.println(sortedList.toString());
		return sortedList;
	}
	
	private ArrayList<Double> entryListKeySet(
			ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap) {
		// TODO Auto-generated method stub
		ArrayList<Double> newKeyList = new ArrayList<Double>();
		for (MyEntry<Double, ArrayList<VCFEntry>> entry: probsToEntryMap){
			newKeyList.add(entry.getKey());
		}
		return newKeyList;
	}

	private void entryListPut(
			ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap,
			double total, ArrayList<VCFEntry> currEntries) {
		MyEntry<Double, ArrayList<VCFEntry>> newEntry = new MyEntry<Double, ArrayList<VCFEntry>>(total, currEntries);
		probsToEntryMap.add(newEntry);
		// TODO Auto-generated method stub
		
	}

	private ArrayList<VCFEntry> entryListGet(
			ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap,
			double total) {
		// TODO Auto-generated method stub
		ArrayList<VCFEntry> returnList = new ArrayList<VCFEntry>();
		for (MyEntry<Double, ArrayList<VCFEntry>> currEntry : probsToEntryMap){
			//ArrayList<VCFEntry> currVCFEntries = currEntry.getValue();
			if (currEntry.getKey().doubleValue() == (total)) returnList.addAll(currEntry.getValue());
		}
		return returnList;
	}

	private boolean entryListContainsKey(
			ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap,
			double total) {
		for (MyEntry<Double, ArrayList<VCFEntry>> currEntry: probsToEntryMap){
			if (currEntry.getKey().equals(total)) return true;
		}
		return false;
	}

	/**
	 * Function: getValidEntries(String inputCode, String destCode, ArrayList<VCFEntry> validEntries)
	 * Usage: int numValid = db.getValidEntries(String inputCode, String destCode, ArrayList<VCFEntry> validEntries)
	 * ----
	 * Finds all valid entries which have the GATK code equivalent to inputCode and can be converted to
	 * destCode. These valid entries are stored in validEntries, which should be an empty ArrayList that is
	 * passed into the function. The function returns the number of valid entries stored in validEntries.
	 * 
	 * @param inputCode	The GATK code to convert
	 * @param destCode	The GATK code to convert to
	 * @param editDistance 
	 * @param validEntries	The ArrayList of entries in which valid entries are stored
	 * @param failCodes 
	 * @param totalMut 
	 * @param pw 
	 * @return	The number of validEntries found as an int
	 */
	public int getValidEntries(String inputCode, String destCode, double editDistance, ArrayList<VCFEntry> validEntries, ArrayList<VCFEntry> failCodes, Integer totalMut){
		//ArrayList<VCFEntry> entries = getSortedEntriesByGATK(inputCode, destCode, pw);
		ArrayList<VCFEntry> entries = getEntriesByGATK(inputCode);
		//System.out.println("For code " + inputCode + " and dest " + destCode + " we found " + entries.size() + " entries.");
		//System.out.println("For code " + inputCode + " we found " + entries.size() + " entries.");
		if (totalMut != null) totalMut = entries.size();
		//For each mismatch in an entry 
		//boolean is0to1 = checkIf0to1(inputCode, destCode);
		Map<Integer, Boolean> mismatchMap = new HashMap<Integer, Boolean>();
		for (int j = 0; j < inputCode.length(); j++){
			if (inputCode.charAt(j) != destCode.charAt(j)) {
				if (inputCode.charAt(j) == '0') mismatchMap.put(j, true);
				else mismatchMap.put(j, false);
			}
		}
		boolean isAll0to1 = true;
		for (Boolean is0to1 : mismatchMap.values()) if (!is0to1) isAll0to1 = false;
		int counter = 0;
		for (int i = 0; i < entries.size(); i++){
			//Note this currently only works for edit distance == 1
			//WELL NOT ANYMORE!
			//int sampleIndex = mismatchIndex(inputCode, destCode);
			if (isAll0to1){
				//check each probability. If all pass, entries are valid
				boolean canConvertAll = true;
				for (Integer indexKey : mismatchMap.keySet()){
					double indexProb = entries.get(i).getSumProb(indexKey);
					if (indexProb >= THRESHOLD){
						canConvertAll = false;
						failCodes.add(entries.get(i));
						break;
					}
				}
				if (canConvertAll){
					counter++;
					validEntries.add(entries.get(i));
				}
			} else {
				//for (int j = entries.size() - 1; j >= 0; j--){
				VCFEntry entry = entries.get(i);
				failCodes.add(entry);
				//}
			}
		}
//			//If all 0 -> 1 keep, else fail
//			double entryProb = entries.get(i).getSumProb(sampleIndex);
//			if (entryProb < THRESHOLD){
//				counter++;
//				VCFEntry entry = entries.get(i);
//				validEntries.add(entry);
//				//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + destCode + "T\n";
//				//System.out.println(entryLine);
//				//Check these
//				//pw.write(entryLine);
//			} else {
//				VCFEntry entry = entries.get(i);
//				failCodes.add(entry);
//				//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + inputCode + "\n";
//				//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + destCode + "F\n";
//				//System.out.println(entryLine);
//				//pw.write(entryLine);
//				//break;
//			}
//		}
//		if (is0to1){
//			//pw.write("--Testing Codes--\n");
//			//ArrayList<VCFEntry> currFailCodes = new ArrayList<VCFEntry>();
//			for (int i = 0; i < entries.size(); i++){
//				//Note this currently only works for edit distance == 1
//				int sampleIndex = mismatchIndex(inputCode, destCode);
//				double entryProb = entries.get(i).getSumProb(sampleIndex);
//				if (entryProb < THRESHOLD){
//					counter++;
//					VCFEntry entry = entries.get(i);
//					validEntries.add(entry);
//					//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + destCode + "T\n";
//					//System.out.println(entryLine);
//					//Check these
//					//pw.write(entryLine);
//				} else {
//					VCFEntry entry = entries.get(i);
//					failCodes.add(entry);
//					//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + inputCode + "\n";
//					//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + destCode + "F\n";
//					//System.out.println(entryLine);
//					//pw.write(entryLine);
//					//break;
//				}
//			}
//			//if 
//			//pw.write("--Testing Ended--\n");
//		} else {
//			for (int i = entries.size() - 1; i >= 0; i--){
//				VCFEntry entry = entries.get(i);
//				failCodes.add(entry);
//				//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + inputCode + "\n";
//				//pw.write(entryLine);
////				int sampleIndex = mismatchIndex(inputCode, destCode);
////				double entryProb = entries.get(i).getSumProb(sampleIndex);
////				if (entryProb >= THRESHOLD){
////					counter++;
////					validEntries.add(entries.get(i));
////				} else break;
//			}
//		}
		return counter;
	}
	
	/**
	 * Function: mismatchIndex(String inputCode, String destCode)
	 * Usage: int index = mismatchIndex(inputCode, destCode)
	 * ----
	 * Returns the index of the first difference between the two
	 * strings. NOTE: IF EDIT DISTANCE NEEDS TO BE >1, THIS MUST
	 * BE CHANGED.
	 * 
	 * @param inputCode	The first GATK code to compare
	 * @param destCode The second GATK code to compare
	 * @return	the first index at which the characters in the strings do not match
	 */
	private int mismatchIndex(String inputCode, String destCode) {
		for (int i = 0; i < inputCode.length(); i++){
			if (inputCode.charAt(i) != destCode.charAt(i)){
				return i;
			}
		}
		return -1;
	}

	/**
	 * Function: checkIf0to1(String inputCode, String destCode)
	 * Usage: boolean isOto1 = checkIf0to1(inputCode, destCode)
	 * ----
	 * Finds if the difference between the inputCode and the destCode
	 * is from 0 to 1 or from 1 to 0. NOTE: THIS ONLY WORKS IF 
	 * EDIT DISTANCE == 1. OTHERWISE, MUST BE CHANGED.
	 * @param inputCode	The original code
	 * @param destCode	The code to change into
	 * @return	true if going to 0 to 1; else, false
	 */
	private boolean checkIf0to1(String inputCode, String destCode) {
		for (int i = 0; i < inputCode.length(); i++){
			if (inputCode.charAt(i) != destCode.charAt(i)){
				return (inputCode.charAt(i) == '0');
			}
		}
		return true;
	}
	
	/**
	 * Function: initCodeDiffMap(String inputCode, String destCode)
	 * Usage: Map<Integer, Boolean> diffMap = initCodeDiffMap(inputCode, destCode)
	 * ----
	 * Generates a map of mismatched indices between the two input parameters and
	 * whether the mismatch is from 0 to 1. In other words, every key is an index
	 * where the characters between the strings do not match. Every value to a key
	 * is a boolean that tells whether the original string's mismatched character 
	 * was a 0.
	 * @param inputCode	The original GATK code
	 * @param destCode	The GATK code to convert to
	 * @return	The map with the mismatched indices to booleans 
	 */
	private Map<Integer, Boolean> initCodeDiffMap(String inputCode,
			String destCode) {
		Map<Integer, Boolean> result = new HashMap<Integer, Boolean>();
		for (int i = 0; i < inputCode.length(); i++){
			char inputChar = inputCode.charAt(i);
			char destChar = destCode.charAt(i);
			if (inputChar != destChar) {
				if (inputChar == '0') result.put(i, true);
				else result.put(i, false);
			}
		}
		return result;
	}

	/**
	 * Function: generateMatrix(String outputFile)
	 * Usage: generateMatrix(outputFile)
	 * ----
	 * Generates the matrix of GATK codes and mutation numbers needed
	 * to build a phylogenetic tree. The matrix is written to 
	 * the file specified by the path outputFile.
	 * @param outputFile	The file which to write the matrix
	 */
	public void generateMatrix(String outputFile){
		Map<String, Integer> GATKCounter = new HashMap<String, Integer>();
		int codeLength = -1;
		for (VCFEntry entry: database){
			String code = entry.getGATK();
			if (codeLength == -1) codeLength = code.length();
			if (GATKCounter.containsKey(code)){
				GATKCounter.put(code, GATKCounter.get(code) + 1);
			} else {
				GATKCounter.put(code, 1);
			}
		}
		//Quick-fix to remove the entry for all 1's
		String imposCode = "";
		for (int i = 0; i < codeLength; i++){
			imposCode += "1";
		}
		if (GATKCounter.containsKey(imposCode)) GATKCounter.remove(imposCode);
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(outputFile));
			int rows = GATKCounter.size();
			int cols = database.get(0).getGATK().length();
			pw.write(rows + " " + cols + "\n");
			ArrayList<String> keySet = new ArrayList<String>(GATKCounter.keySet());
			for (int i = 0; i < rows; i++){
				String key = largestKey(keySet, GATKCounter);
				if (i < rows - 1){
					pw.write(key + " " + GATKCounter.get(key) + "\n");
				} else {
					pw.write(key + " " + GATKCounter.get(key));
				}
				keySet.remove(key);
			}
			pw.close();
		} catch (IOException e) {
			//e.printStackTrace();
			System.out.println("Caught Unknown IO Exception!");
		}
	}

	/**
	 * Function: largestKey(ArrayList<String> keySet, Map<String, Integer> GATKCounter)
	 * Usage: String key = largestKey(keySet, GATKCounter)
	 * ----
	 * Finds the largest key in the map. The keySet may only be a subset of the keys in
	 * GATKCounter. This function is generally used to sort the map.
	 * @param keySet	The remaining keys to sort
	 * @param GATKCounter	The map of GATK codes to mutation numbers
	 * @return	The largest key that is in keySet
	 */
	private String largestKey(ArrayList<String> keySet,
			Map<String, Integer> GATKCounter) {
		// TODO Auto-generated method stub
		int maxCount = 0;
		String topKey = "";
		for (String key: keySet){
			int currCount = Integer.parseInt(key, 2);
			if (currCount > maxCount){
				topKey = key;
				maxCount = currCount;
			}
		}
		return topKey;
	}

	public void generateGATKFile(String outputFile) {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new FileWriter(outputFile));
		} catch (IOException e){
			System.out.println("IOEXCEPTION OCCURRED!");
		}
		String header = "Chromosome\tLocation\tGATK\n";
		pw.print(header);
		for (int i = 0; i < database.size(); i++){
			VCFEntry entry = database.get(i);
			String outputStr = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + entry.getGATK() + "\n";
			pw.print(outputStr);
		}
		pw.close();
		// TODO Auto-generated method stub
	}
	
	//Taken from StackOverflow
	//http://stackoverflow.com/questions/3110547/java-how-to-create-new-entry-key-value
	final class MyEntry<K, V> implements Map.Entry<K, V> {
	    private final K key;
	    private V value;

	    public MyEntry(K key, V value) {
	        this.key = key;
	        this.value = value;
	    }

	    @Override
	    public K getKey() {
	        return key;
	    }

	    @Override
	    public V getValue() {
	        return value;
	    }

	    @Override
	    public V setValue(V value) {
	        V old = this.value;
	        this.value = value;
	        return old;
	    }
	}

}
