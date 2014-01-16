package jPhyTree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class HelloWorld {

	private static VCFDatabase db;
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//System.out.println("Hello, world!");
		db = new VCFDatabase("testCases/simulation_vcfs/tree_6_12.raw.vcf");
		OutputRatio(4, "000110", "000011");
//		BufferedReader br = null;
//		try{
//			br = new BufferedReader(new FileReader("Allele Ratios.txt"));
//			String currLine;
//			while ((currLine = br.readLine()) != null){
//			//int counter = 0;
//				//System.out.println(currLine);
//				String[] lineParts = currLine.split("\\s+");
//				
//				if (Character.isDigit(lineParts[0].charAt(0))) OutputRatio(Integer.parseInt(lineParts[0]), lineParts[1], lineParts[2]);
//				else System.out.println(currLine);
//			}
//		} catch (IOException e) {
//			System.out.println("Unforeseen IOException!");
//		}
//		db.generateMatrix("output");
	}
	private static void OutputRatio(int sample, String codeA, String codeB) {
		// TODO Auto-generated method stub
//		int sample = 1;
//		String codeA = "110";
//		String codeB = "011";
		ArrayList<VCFEntry> test = db.getEntriesByGATK(codeA);
		double totRatio = 0.0;
		int numEntries = 0;
		for (VCFEntry entry: test){
			int majAllele = entry.getAlleleCount(sample, 0);
			int readDepth = entry.getReadDepth(sample);
			//System.out.println(entry.toString());
			totRatio += (double) majAllele / readDepth;
			numEntries++;
		}
		double ratioA = totRatio / numEntries;
		totRatio = 0.0;
		numEntries = 0;
		test = db.getEntriesByGATK(codeB);
		for (VCFEntry entry: test){
			int majAllele = entry.getAlleleCount(sample, 0);
			int readDepth = entry.getReadDepth(sample);
			//System.out.println(entry.toString());
			totRatio += (double) majAllele / readDepth;
			numEntries++;
		}
		double ratioB = totRatio / numEntries;
		System.out.println(sample + "\t" + codeA + "\t" + codeB + "\t" + ratioA / ratioB);
	}

}
