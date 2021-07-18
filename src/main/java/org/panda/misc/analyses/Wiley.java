package org.panda.misc.analyses;

import org.panda.causalpath.run.JasonizeResultGraphsRecursively;
import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.Overlap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Wiley
{
	static final String DATA_DIR = "/Users/ozgun/Documents/Data/Wiley/";
	static final String BASE = "/Users/ozgun/Documents/Analyses/Wiley/";

	static Set<String> idMemory = new HashSet<>();
	static final String SITE_NAME = "Residue";
	static final String SYM_NAME = "Gene";
	static final String VAL_PREFIX = "Fold";

	public static void main(String[] args) throws IOException
	{
//		convertData();
//		generateAnalysisFolders();
//		analyzeAgreement();
//		writeEGFDownstream();
//		writeTimecourseDosesIntersectionGraph();
		writeConflicting();
	}

	static void convertData() throws IOException
	{
		convertData(DATA_DIR + "timecourse.csv", BASE + "data-timecourse.csv");
		convertData(DATA_DIR + "doses-run1.csv", BASE + "data-doses-run1.csv");
		convertData(DATA_DIR + "doses-run2.csv", BASE + "data-doses-run2.csv");
	}

	static void convertData(String inFile, String outFile) throws IOException
	{
		idMemory.clear();
		String[] header = FileUtil.readHeader(inFile);
		int siteIndex = ArrayUtil.indexOf(header, SITE_NAME);
		int symIndex = ArrayUtil.indexOf(header, SYM_NAME);
		List<Integer> valInds = new ArrayList<>();
		for (int i = 0; i < header.length; i++)
		{
			if (header[i].startsWith(VAL_PREFIX)) valInds.add(i);
		}

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tModification\tEffect");
		valInds.forEach(i -> FileUtil.tab_write(header[i], writer));

		FileUtil.linesTabbedSkip1(inFile).filter(t -> !t[symIndex].isEmpty()).forEach(t ->
		{
			String sym = t[symIndex];
			List<String> sites = parseSites(t[siteIndex]);
			String id = getID(sym, sites);

			FileUtil.lnwrite(id + "\t" + sym + "\t" + CollectionUtil.merge(sites, "|") + "\tP\t", writer);

			valInds.forEach(i ->
			{
				double v = Double.valueOf(t[i]);
				if (v < 1) v = -(1 / v);
				FileUtil.tab_write(v, writer);
			});
		});
		writer.close();
	}

	static List<String> parseSites(String s)
	{
		return Arrays.stream(s.split("s|t|y")).filter(site -> !site.isEmpty()).collect(Collectors.toList());
	}

	static String getID(String sym, List<String> sites)
	{
		String pref = sym + "-" + CollectionUtil.merge(sites, "-");

		String id = pref + "-P";
		int repCounter = 0;

		while (idMemory.contains(id))
		{
			id = pref + "-" + (++repCounter) + "-P";
		}

		idMemory.add(id);
		return id;
	}

	static void generateAnalysisFolders()
	{
		generateAnalysisFolders("data-timecourse.csv");
		generateAnalysisFolders("data-doses-run1.csv");
		generateAnalysisFolders("data-doses-run2.csv");
	}

	static void generateAnalysisFolders(String dataFile)
	{
		String base = dataFile.substring(dataFile.indexOf("-") + 1, dataFile.lastIndexOf("."));
		base = BASE + base + "/";

		String[] header = FileUtil.readHeader(BASE + dataFile);
		for (int i = 5; i < header.length; i++)
		{
			String valCol = header[i];
			String valDir = valCol.substring(valCol.indexOf(" "));
			valDir = valDir.replaceAll(" ", "");
			valDir = valDir.replaceAll("@", "dose");

			for (double thr = 1.25; thr <= 2; thr+=0.25)
			{
				String dir = base + valDir + "/fc" + thr + "/";
				FileUtil.mkdirs(dir);
				FileUtil.writeStringToFile(getParameterFileText(dataFile, valCol, thr), dir + "parameters.txt");
			}
		}
	}

	static String getParameterFileText(String dataFile, String valCol, double thr)
	{
		return
			"proteomics-values-file = ../../../" + dataFile + "\n" +
			"id-column = ID\n" +
			"symbols-column = Symbols\n" +
			"sites-column = Sites\n" +
			"modification-column = Modification\n" +
			"effect-column = Effect\n" +
			"\n" +
			"value-transformation = as-is-single-value\n" +
			"\n" +
			"threshold-for-data-significance = " + thr + " phosphoprotein\n" +
			"color-saturation-value = 5\n" +
			"show-all-genes-with-proteomic-data = true\n" +
			"\n" +
			"calculate-network-significance = true\n" +
			"use-network-significance-for-causal-reasoning = true\n" +
			"permutations-for-significance = 100000\n" +
			"fdr-threshold-for-network-significance = 0.1\n" +
			"\n" +
			"value-column = " + valCol;
	}

	static void analyzeAgreement()
	{
		Map<String, String> tcMap = FileUtil.readMap(BASE + "data-timecourse.csv", "\t", "ID", "Fold 4 min");
		Map<String, String> r1Map = FileUtil.readMap(BASE + "data-doses-run1.csv", "\t", "ID", "Fold @1");
		Map<String, String> r2Map = FileUtil.readMap(BASE + "data-doses-run2.csv", "\t", "ID", "Fold @1");

		System.out.println("Presence of measurements:");
		CollectionUtil.printNameMapping("Timecourse", "Doses Run 1", "Doses Run 2");
		CollectionUtil.printVennCounts(tcMap.keySet(), r1Map.keySet(), r2Map.keySet());
		printSetStats(tcMap.keySet(), r1Map.keySet(), r2Map.keySet());


//		Set<String> common = CollectionUtil.getIntersection(tcMap.keySet(), r1Map.keySet(), r2Map.keySet());
		Set<String> common = CollectionUtil.getIntersection(tcMap.keySet(), r1Map.keySet());

		BufferedWriter writer = FileUtil.newBufferedWriter(BASE + "common.txt");
		for (String id : common)
		{
			double x = Double.valueOf(tcMap.get(id));
			double y = Double.valueOf(r1Map.get(id));

			FileUtil.writeln(x + "\t" + y, writer);
		}
		FileUtil.closeWriter(writer);

		for (double i = 1.5; i <= 2; i+=0.25)
		{
			System.out.println("\nthr = " + i);
			double thr = i;

			Set<String> tcpSet = tcMap.keySet().stream().filter(common::contains).filter(id -> discretize(tcMap.get(id), thr).equals("+")).collect(Collectors.toSet());
			Set<String> tcnSet = tcMap.keySet().stream().filter(common::contains).filter(id -> discretize(tcMap.get(id), thr).equals("-")).collect(Collectors.toSet());
			Set<String> r1pSet = r1Map.keySet().stream().filter(common::contains).filter(id -> discretize(r1Map.get(id), thr).equals("+")).collect(Collectors.toSet());
			Set<String> r1nSet = r1Map.keySet().stream().filter(common::contains).filter(id -> discretize(r1Map.get(id), thr).equals("-")).collect(Collectors.toSet());

			CollectionUtil.printVennCounts(tcpSet, tcnSet, r1pSet, r1nSet);
			printSetStats(tcpSet, tcnSet, r1pSet, r1nSet);
			System.out.println("Positive overlap");
			printRandomOverlapProb(common, tcpSet, r1pSet);
			System.out.println("Negative overlap");
			printRandomOverlapProb(common, tcnSet, r1nSet);
			System.out.println("tc+ r1- overlap");
			printRandomOverlapProb(common, tcpSet, r1nSet);
			System.out.println("tc- r1+ overlap");
			printRandomOverlapProb(common, tcnSet, r1pSet);


//			CollectionUtil.printVennCounts(tcSet, r1Set, r2Set);
//			System.out.println("Overlap shrinkage = " + CollectionUtil.getOverlapShrinkage(tcSet, r1Set, r2Set));
		}
	}

	static String discretize(String valStr, double thr)
	{
		double fc = Double.valueOf(valStr);
		if (Math.abs(fc) < thr) return "0";
		return fc > 0 ? "+" : "-";
	}

	static void printSetStats(Set<String> s1, Set<String> s2, Set<String> s3)
	{
		System.out.println("s1.size() = " + s1.size());
		System.out.println("s2.size() = " + s2.size());
		System.out.println("s3.size() = " + s3.size());
		System.out.println("s1-s2 = " + CollectionUtil.getIntersection(s1, s2).size());
		System.out.println("s1-s3 = " + CollectionUtil.getIntersection(s1, s3).size());
		System.out.println("s2-s3 = " + CollectionUtil.getIntersection(s2, s3).size());
		System.out.println("s1-s2-s3 = " + CollectionUtil.getIntersection(s1, s2, s3).size());
	}

	static void printSetStats(Set<String> s1, Set<String> s2, Set<String> s3, Set<String> s4)
	{
		System.out.println("s1.size() = " + s1.size());
		System.out.println("s2.size() = " + s2.size());
		System.out.println("s3.size() = " + s3.size());
		System.out.println("s4.size() = " + s4.size());
		System.out.println("s1-s2 = " + CollectionUtil.getIntersection(s1, s2).size());
		System.out.println("s1-s3 = " + CollectionUtil.getIntersection(s1, s3).size());
		System.out.println("s1-s4 = " + CollectionUtil.getIntersection(s1, s4).size());
		System.out.println("s2-s3 = " + CollectionUtil.getIntersection(s2, s3).size());
		System.out.println("s2-s4 = " + CollectionUtil.getIntersection(s2, s4).size());
		System.out.println("s3-s4 = " + CollectionUtil.getIntersection(s3, s4).size());
		System.out.println("s2-s3-s4 = " + CollectionUtil.getIntersection(s2, s3, s4).size());
		System.out.println("s1-s3-s4 = " + CollectionUtil.getIntersection(s1, s3, s4).size());
		System.out.println("s1-s2-s4 = " + CollectionUtil.getIntersection(s1, s2, s4).size());
		System.out.println("s1-s2-s3 = " + CollectionUtil.getIntersection(s1, s2, s3).size());
		System.out.println("s1-s2-s3-s4 = " + CollectionUtil.getIntersection(s1, s2, s3, s4).size());
	}

	static void printRandomOverlapProb(Set<String> bg, Set<String> s1, Set<String> s2)
	{
		int n = bg.size();
		int o = CollectionUtil.countOverlap(s1, s2);
		double p = Overlap.calcCoocPval(n, o, s1.size(), s2.size());
		System.out.println("p = " + p);
	}

	static void writeEGFDownstream() throws IOException
	{
		writeEGFDownstreamRecursive(BASE + "timecourse");
		writeEGFDownstreamRecursive(BASE + "doses-run1");
		writeEGFDownstreamRecursive(BASE + "doses-run2");

		JasonizeResultGraphsRecursively.generate(BASE, BASE, Collections.singleton("EGF-downstream"), BASE + "EGF-downstream/", "causative.json");
	}

	static void writeEGFDownstreamRecursive(String base) throws IOException
	{
		FileUtil.processDirsRecursive(new File(base), dir ->
		{
			if (Files.exists(Paths.get(dir + "/results.txt")))
			{
				writeEGFDownstream(dir.getPath());
			}
		});
	}


	static void writeEGFDownstream(String dir) throws IOException
	{
		CausalPathSubnetwork.writeDownstreamSubgraph(dir, Collections.singleton("EGF"), "EGF-downstream");
	}

	static void writeTimecourseDosesIntersectionGraph() throws IOException
	{
		CausalPathSubnetwork.writeIntersectionSubgraph(BASE + "/timecourse/4min/fc2.0",
			"causative", BASE + "doses-run1/dose1/fc2.0/causative.sif",
			"tc-doses-common");
		JasonizeResultGraphsRecursively.generate(BASE, BASE, Collections.singleton("tc-doses-common"), BASE + "common/", "causative.json");
	}

	static void writeConflicting() throws IOException
	{
		JasonizeResultGraphsRecursively.generate(BASE, BASE, Collections.singleton("conflicting"), BASE + "conflicts/", "causative.json");
	}
}
