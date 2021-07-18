package org.panda.misc.analyses;

import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.utility.*;
import org.panda.utility.statistics.BoxPlot;
import org.panda.utility.statistics.TTest;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class CPTACLSCC
{
//	public static final String BASE = "C:\\Users\\Owner\\Documents\\Analyses\\CPTAC-LSCC\\"; // Laptop
	public static final String DIR = "/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/"; // MA

	public static void main(String[] args) throws IOException
	{
//		convertData();
//		printColNames();
//		prepareCosmicGraph();
//		prepareGOIGraph();
//		generateSigNeighGraphs();
//		generatePartialSigGraph();
//		separateEdgeTypes();
//		compareClusters();
//		printResultSimilarities();
//		normalTumorImmuneComparison();
//		separateTurmors();
//		printImmunePCADimComparisons();

		exploreResultOverlaps();
	}

	private static void convertData() throws IOException
	{
		Map<String[], Map<String, Double>> pp = processPhosphoproteomicData();
		Map<String, Map<String, Double>> tp = processProteomicData();

		List<String> samples = pp.values().iterator().next().keySet().stream().sorted().collect(Collectors.toList());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "data.txt"));

		writer.write("ID\tSymbols\tSites\tEffect");
		samples.forEach(s -> FileUtil.tab_write(s, writer));

		pp.forEach((s, map) ->
		{
			FileUtil.lnwrite(ArrayUtil.merge("\t", s), writer);
			samples.forEach(c -> FileUtil.tab_write(map.get(c), writer));
		});

		tp.forEach((s, map) ->
		{
			FileUtil.lnwrite(s + "\t" + s + "\t\t", writer);
			samples.forEach(c -> FileUtil.tab_write(map.get(c), writer));
		});

		writer.close();
	}


	private static Map<String, Map<String, Double>> processProteomicData() throws IOException
	{
		String inFile = DIR + "lscc-v3.2-proteome-ratio-norm-NArm.gct";

		Map<String, String[]> rowMap = new HashMap<>();
		Map<String, Integer> specMap = new HashMap<>();

		Files.lines(Paths.get(inFile)).skip(83).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = t[2];
			Integer spec = Integer.valueOf(t[4]);

			if (!specMap.containsKey(sym) || specMap.get(sym) < spec)
			{
				specMap.put(sym, spec);
				rowMap.put(sym, t);
			}
		});

		System.out.println("specMap.size() = " + specMap.size());

//		printIsoformFreq(rowMap);

		Map<String, Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(inFile)).skip(2).findFirst().get().split("\t");

		rowMap.forEach((sym, t) ->
		{
			HashMap<String, Double> valMap = new HashMap<>();
			for (int i = 18; i < header.length; i++)
			{
				valMap.put(header[i], t[i].equals("NA") ? Double.NaN : Double.valueOf(t[i]));
			}
			map.put(sym, valMap);
		});

		return map;
	}

	private static void printIsoformFreq(Map<String, String[]> rowMap)
	{
		TermCounter tc = new TermCounter();
		for (String[] t : rowMap.values())
		{
			int ind = t[1].indexOf(" isoform ");
			if (ind > 0)
			{
				String num = t[1].substring(ind + 9, t[1].length());
				tc.addTerm(num.split(" ")[0]);
			}
		}

		tc.print();
	}

	private static Map<String[], Map<String, Double>> processPhosphoproteomicData() throws IOException
	{
		String inFile = DIR + "lscc-v3.2-phosphoproteome-ratio-norm-NArm.gct";

		Map<String[], Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(inFile)).skip(2).findFirst().get().split("\t");

		Files.lines(Paths.get(inFile)).skip(83).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = t[2];
			List<String> sites = new ArrayList<>();
			for (String s : t[12].split(" "))
			{
				if (s.startsWith("S") || s.startsWith("T") || s.startsWith("Y"))
				{
					sites.add(s.substring(0, s.length() - 1));
				}
			}

			HashMap<String, Double> valMap = new HashMap<>();
			for (int i = 25; i < header.length; i++)
			{
				valMap.put(header[i], t[i].equals("NA") ? Double.NaN : Double.valueOf(t[i]));
			}

			String id = sym + "-" + CollectionUtil.merge(sites, "-");

			String site = CollectionUtil.merge(sites, "|");

			map.put(new String[]{id, sym, site, ""}, valMap);
		});

		return map;
	}

	private static void printColNames() throws IOException
	{
		String inFile = DIR + "lscc-v3.2-proteome-ratio-norm-NArm.gct";
		String dir = null;
		String[] header = Files.lines(Paths.get(inFile)).skip(2).findFirst().get().split("\t");
		String[] tp53 = Files.lines(Paths.get(inFile)).filter(l -> l.startsWith("TP53.mutation.status")).findFirst().get().split("\t");
		String[] type = Files.lines(Paths.get(inFile)).skip(8).findFirst().get().split("\t");
		String[] cluster = Files.lines(Paths.get(inFile)).skip(59).findFirst().get().split("\t");

		// For correlation
//		for (int i = 18; i < header.length; i++)
//		{
//			if (type[i].equals("Tumor"))
//			{
//				System.out.println("value-column = " + header[i]);
//			}
//		}


		// Tumor vs normal
//		Set<String> normals = Arrays.stream(header).filter(s -> s.endsWith(".N")).collect(Collectors.toSet());
//
//		for (int i = 18; i < header.length; i++)
//		{
//			if (type[i].equals("Tumor") && normals.contains(header[i] + ".N"))
//			{
//				System.out.println("test-value-column = " + header[i]);
//				System.out.println("control-value-column = " + header[i] + ".N");
//			}
//		}

		// Print distribution of TP53 status to clusters
//		TermCounter tc = new TermCounter();
//		for (int i = 0; i < cluster.length; i++)
//		{
//			if (cluster[i].startsWith("C"))
//			{
//				tc.addTerm(cluster[i] + (tp53[i].equals("1") ? "+" : "-"));
//			}
//		}
//		tc.print();

//		dir = BASE + "clusters/one-off/P53-intact-vs-normals/";
//		Files.createDirectories(Paths.get(dir));
//		BufferedWriter writer5 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
//		Files.lines(Paths.get(BASE + "clusters/parameters.template.txt")).map(l -> l.endsWith("-mean") ? l + "-paired" : l).forEach(l -> FileUtil.writeln(l, writer5));
//		for (int i = 16; i < header.length; i++)
//		{
//			if (type[i].equals("NAT"))
//			{
//				String tum = header[i].substring(0, header[i].length() - 2);
//				int tInd = ArrayUtil.indexOf(header, tum);
//				if (tp53[tInd].equals("0"))
//				{
//					FileUtil.lnwrite("control-value-column = " + header[i], writer5);
//					FileUtil.lnwrite("test-value-column = " + header[tInd], writer5);
//				}
//			}
//		}
//		writer5.close();
//
//		dir = BASE + "clusters/one-off/P53-mut-in-C3/";
//		Files.createDirectories(Paths.get(dir));
//		BufferedWriter writer6 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
//		Files.lines(Paths.get(BASE + "clusters/parameters.template.txt")).forEach(l -> FileUtil.writeln(l, writer6));
//		for (int i = 16; i < header.length; i++)
//		{
//			if (type[i].equals("Tumor") && cluster[i].equals("C3"))
//			{
//				if (tp53[i].equals("1")) FileUtil.lnwrite("test-value-column = " + header[i], writer6);
//				else FileUtil.lnwrite("control-value-column = " + header[i], writer6);
//			}
//		}
//		writer6.close();



		// Determine normals paired to tumor clusters
		String[] correspNormal = new String[cluster.length];
		for (int i = 0; i < cluster.length; i++)
		{
			String name = header[i];
			if (name.endsWith(".N"))
			{
				name = name.substring(0, name.length() - 2);
				int ind = ArrayUtil.indexOf(header, name);
				if (ind > 0)
				{
					String c = cluster[ind];
					if (!c.equalsIgnoreCase("NA"))
					{
						correspNormal[i] = c;
					}
				}
			}
			if (correspNormal[i] == null) correspNormal[i] = "NA";
		}

//		generateBoxPlot("MCM3", inFile, cluster, correspNormal);


		// Clusters

//		List<String> clusters = Arrays.stream(cluster).filter(c -> !c.equalsIgnoreCase("NA") && !c.equals("NMF.cluster")).distinct().sorted().collect(Collectors.toList());
//		for (int i = 0; i < clusters.size(); i++)
//		{
//			String testC = clusters.get(i);
//
//			dir = DIR + "clusters/against-others/" + testC + "/";
//			Files.createDirectories(Paths.get(dir));
//			BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
//			Files.lines(Paths.get(DIR + "clusters/parameters.template.txt")).forEach(l -> FileUtil.writeln(l, writer1));
//			for (int j = 18; j < header.length; j++)
//			{
//				if (cluster[j].equals(testC)) FileUtil.lnwrite("test-value-column = " + header[j], writer1);
//				else if (!cluster[j].equalsIgnoreCase("NA")) FileUtil.lnwrite("control-value-column = " + header[j], writer1);
//			}
//			writer1.close();
//
//			dir = DIR + "clusters/against-normals/" + testC + "/";
//			Files.createDirectories(Paths.get(dir));
//			BufferedWriter writer3 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
//			Files.lines(Paths.get(DIR + "clusters/parameters.template.txt")).map(l -> l.endsWith("-mean") ? l + "-paired" : l).forEach(l -> FileUtil.writeln(l, writer3));
//			List<String> tumors = new ArrayList<>();
//			List<String> normals = new ArrayList<>();
//			for (int j = 18; j < header.length; j++)
//			{
//				if (cluster[j].equals(testC)) tumors.add(header[j]);
//			}
//			for (int j = 18; j < header.length; j++)
//			{
//				if (header[j].endsWith(".N") && tumors.contains(header[j].substring(0, header[j].length() - 2))) normals.add(header[j]);
//			}
//			Set<String> unmatchedT = tumors.stream().filter(t -> !normals.contains(t + ".N")).collect(Collectors.toSet());
//			tumors.removeAll(unmatchedT);
//
//			tumors.stream().sorted().forEach(s -> FileUtil.lnwrite("test-value-column = " + s, writer3));
//			normals.stream().sorted().forEach(s -> FileUtil.lnwrite("control-value-column = " + s, writer3));
//			writer3.close();
//
//			for (int j = i+1; j < clusters.size(); j++)
//			{
//				String ctrlC = clusters.get(j);
//
//				dir = DIR + "clusters/pairs/" + testC + "-vs-" + ctrlC + "/";
//				Files.createDirectories(Paths.get(dir));
//				BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
//				Files.lines(Paths.get(DIR + "clusters/parameters.template.txt")).forEach(l -> FileUtil.writeln(l, writer2));
//				for (int k = 18; k < header.length; k++)
//				{
//					if (cluster[k].equals(testC)) FileUtil.lnwrite("test-value-column = " + header[k], writer2);
//					else if (cluster[k].equals(ctrlC)) FileUtil.lnwrite("control-value-column = " + header[k], writer2);
//				}
//				writer2.close();
//
//				dir = DIR + "clusters/normal-pairs/" + testC + "-vs-" + ctrlC + "/";
//				Files.createDirectories(Paths.get(dir));
//				BufferedWriter writer4 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
//				Files.lines(Paths.get(DIR + "clusters/parameters.template.txt")).forEach(l -> FileUtil.writeln(l, writer4));
//				for (int k = 18; k < header.length; k++)
//				{
//					if (correspNormal[k].equals(testC)) FileUtil.lnwrite("test-value-column = " + header[k], writer4);
//					else if (correspNormal[k].equals(ctrlC)) FileUtil.lnwrite("control-value-column = " + header[k], writer4);
//				}
//				writer4.close();
//			}
//		}

		// Mutations

		Set<String> oncogenes = new HashSet<>(Arrays.asList("PIK3CA", "KRAS", "NOTCH1", "HRAS", "EGFR", "IL21R", "EGFL6",
			"LMO2", "BIRC6", "BRAF", "ARAF", "ERBB2"));

		Files.lines(Paths.get(inFile)).map(l -> l.split("\t")).filter(t -> t[0].endsWith(".mutation.status")).forEach(t ->
		{
			String gene = t[0].substring(0, t[0].indexOf("."));
			String directory = DIR + "mutations/" + gene + "/";
			FileUtil.mkdirs(directory);
			BufferedWriter writer = FileUtil.newBufferedWriter(directory + "parameters.txt");
			FileUtil.lines(DIR + "mutations/parameters.template.txt").forEach(l -> FileUtil.writeln(l, writer));

			FileUtil.lnwrite("gene-activity = " + gene + " " + (oncogenes.contains(gene) ? "a" : "i") + "\n", writer);

			for (int i = 18; i < header.length; i++)
			{
				if (!type[i].equals("Tumor")) continue;

				if (t[i].equals("1")) FileUtil.lnwrite("test-value-column = " + header[i], writer);
				else if (t[i].equals("0")) FileUtil.lnwrite("control-value-column = " + header[i], writer);
			}
			FileUtil.closeWriter(writer);
		});
	}

	private static void generateBoxPlot(String id, String origDataSheet, String[] cluster, String[] correspNormal) throws IOException
	{
		List<String[]> rows = Files.lines(Paths.get(origDataSheet)).map(l -> l.split("\t")).filter(t -> t.length > 2 && t[2].equals(id)).collect(Collectors.toList());
		for (String[] vals : rows)
		{
			String[]header = new String[]{"N1", "C1", "N2", "C2", "N3", "C3"};
			List<Double>[] cols = new List[header.length];
			for (int i = 0; i < cols.length; i++)
			{
				cols[i] = new ArrayList<>();
			}
			for (int i = 16; i < vals.length; i++)
			{
				double v = Double.valueOf(vals[i]);
				if (!Double.isNaN(v))
				{
					if (cluster[i].startsWith("C"))
					{
						int ind = ArrayUtil.indexOf(header, cluster[i]);
						cols[ind].add(v);
					}
					else if (correspNormal[i].startsWith("C"))
					{
						int ind = ArrayUtil.indexOf(header, correspNormal[i].replace("C", "N"));
						cols[ind].add(v);
					}
				}
			}
			BoxPlot.write(DIR + "box-plots/" + id + "-" + vals[0] + ".txt", header, cols);
		}
	}

	private static String[] readClusters(String[] samples, String annotFile) throws IOException
	{
		String[] header = Files.lines(Paths.get(annotFile)).findFirst().get().split("\t");
		int idInd = ArrayUtil.lastIndexOf(header, "id");
		int cInd = ArrayUtil.lastIndexOf(header, "NMF.consensus");

		Map<String, String> map = Files.lines(Paths.get(annotFile)).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[idInd], t -> t[cInd]));

		String[] clust = new String[samples.length];
		for (int i = 0; i < clust.length; i++)
		{
			if (map.containsKey(samples[i]))
			{
				clust[i] = "C" + map.get(samples[i]);
			}
			else clust[i] = "NA";
		}
		return clust;
	}

	private static void prepareCosmicGraph() throws IOException
	{
		CausalPathSubnetwork.writeGOINeighForCorrBased(DIR + "correlation-tumor", COSMIC, StreamDirection.UPSTREAM);
		CausalPathSubnetwork.writeGOINeighForCompBased(DIR + "tumor-vs-normal", COSMIC, StreamDirection.UPSTREAM);

//		CausalPathSubnetwork.writeGOINeighForCorrBased(BASE + "correlation-tumor", Collections.singleton("RUNX1"), StreamDirection.BOTHSTREAM);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "cosmic.highlight"));
		COSMIC.stream().sorted().forEach(g -> FileUtil.writeln("node\t" + g, writer));
		writer.close();
	}


	private static void prepareGOIGraph() throws IOException
	{
		CausalPathSubnetwork.writeGOINeighForCorrBased(DIR + "clusters/against-others/C3", new HashSet<String>(Arrays.asList(new String[]{"KEAP1", "CUL3", "NFE2L2"})), StreamDirection.BOTHSTREAM, "KEAP1-CUL3-NFE2L2");
	}

	private static void generateSigNeighGraphs() throws IOException
	{
		CausalPathSubnetwork.generateNeighborhoodSubgraphsForSignificantsRecursively(DIR + "immune-PCA-axis", 0.1);
//		CausalPathSubnetwork.writeSignifNeighForCompBased(BASE + "clusters/one-off/P53-intact-vs-normals", StreamDirection.DOWNSTREAM, 0.1);

//		CausalPathSubnetwork.writeSignifNeighForCorrBased(BASE + "correlation-tumor", StreamDirection.DOWNSTREAM);
	}

	private static void generatePartialSigGraph() throws IOException
	{
		// Update the below line and DEFINERS
		String dir = DIR + "clusters/pairs/C2-vs-C3/";
		String[] roots = DEFINERS.split("\n\n");

		for (int i = 0; i < roots.length; i++)
		{
			Set<String> goi = new HashSet<>(Arrays.asList((roots[i]).split("\n")));
			String subname = "Part-" + (i + 1);
			SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/causative-" + subname + ".sif", StreamDirection.DOWNSTREAM);
			Set<String> ids = CausalPathSubnetwork.getIDsAtTheNeighborhood(dir + "/results.txt", goi, StreamDirection.DOWNSTREAM);
			System.out.println("ids.size() = " + ids.size());
			CausalPathSubnetwork.writeSubsetFormat(dir + "/causative.format", dir + "/causative-" + subname + ".format", goi, ids);
		}
	}

	private static void separateEdgeTypes() throws IOException
	{
		String file = DIR + "clusters/pairs/C2-vs-C3/causative-Part-2";
		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(file + "-phospho.sif"));
		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(file + "-express.sif"));

		Files.lines(Paths.get(file + ".sif")).forEach(l ->
		{
			BufferedWriter writer = l.contains("phospho") ? writer1 : writer2;
			FileUtil.writeln(l, writer);
		});

		writer1.close();
		writer2.close();
	}

	private static void printResultSimilarities() throws IOException
	{
		double score = getSimilarityScore(DIR + "clusters/against-normals/C3", DIR + "clusters/against-normals/C2");
		System.out.println("score = " + score);

		Set<String> sif = SIFFileUtil.getRelationsAsString(DIR + "clusters/pairs/C2-vs-C3/causative.sif");
		System.out.println("sif.size() = " + sif.size());
	}

	private static double getSimilarityScore(String dir1, String dir2) throws IOException
	{
		Set<String> signedEdges1 = readSignedEdges(dir1 + "/results.txt");
		Set<String> signedEdges2 = readSignedEdges(dir2 + "/results.txt");

		Set<String> diff1 = new HashSet<>();
		Set<String> diff2 = new HashSet<>();
		Set<String> comm = new HashSet<>();
		Set<String> opp = new HashSet<>();

		for (String edge : signedEdges1)
		{
			if (signedEdges2.contains(edge))
			{
				comm.add(edge);
			}
			else
			{
				String oEdge = opposite(edge);
				if (!signedEdges1.contains(oEdge) && signedEdges2.contains(oEdge))
				{
					opp.add(edge);
				}
				else
				{
					diff1.add(edge);
				}
			}
		}

		for (String edge : signedEdges2)
		{
			if (!comm.contains(edge))
			{
				String oEdge = opposite(edge);
				if (!opp.contains(oEdge))
				{
					diff2.add(edge);
				}
			}
		}

		return (comm.size() - opp.size()) / (double)(diff1.size() + diff2.size() + comm.size() + opp.size());
	}

	private static String opposite(String edge)
	{
		if (edge.endsWith("-")) return edge.substring(0, edge.length() - 1) + "+";
		else return edge.substring(0, edge.length() - 1) + "-";
	}

	private static Set<String> readSignedEdges(String file) throws IOException
	{
		return Files.lines(Paths.get(file)).skip(1).map(l -> l.split("\t"))
			.map(t -> ArrayUtil.getString(" ", t[0], t[1], t[2], (t[8].startsWith("-") ? "-" : "+")))
			.collect(Collectors.toSet());
	}

	//-------

	private static void compareClusters() throws IOException
	{
		vennCompareTwoCPComparisons(DIR + "clusters/against-normals/C3/", DIR + "clusters/pairs/C1-vs-C3/", 0.01);
	}

	private static void vennCompareTwoCPComparisons(String dir1, String dir2, double pvalThr) throws IOException
	{
		Set<String> ids1 = readMeasuredIDs(dir1);
		Set<String> ids2 = readMeasuredIDs(dir2);
		Set<String> ids = CollectionUtil.getIntersection(ids1, ids2);
		Map<String, Double> vc1 = readCPValChanges(dir1, pvalThr);
		Map<String, Double> vc2 = readCPValChanges(dir2, pvalThr);

		int[] cnt = new int[4];
		int diff1Ind = 0;
		int diff2Ind = 1;
		int intersectInd = 2;
		int oppositeInd = 3;

		for (String id : ids)
		{
			if (vc1.containsKey(id))
			{
				double v1 = vc1.get(id);

				if (vc2.containsKey(id))
				{
					double v2 = vc2.get(id);

					if (v1 * v2 < 0)
					{
						cnt[oppositeInd]++;
//						System.out.println("id = " + id);
					}
					else cnt[intersectInd]++;
				}
				else cnt[diff1Ind]++;
			}
			else if (vc2.containsKey(id)) cnt[diff2Ind]++;
		}

		System.out.println("cnt[diff1Ind] = " + cnt[diff1Ind]);
		System.out.println("cnt[diff2Ind] = " + cnt[diff2Ind]);
		System.out.println("cnt[intersectInd] = " + cnt[intersectInd]);
		System.out.println("cnt[oppositeInd] = " + cnt[oppositeInd]);
	}

	private static Map<String, Double> readCPValChanges(String dir, double pvalThr) throws IOException
	{
		return Files.lines(Paths.get(dir + "value-changes.txt")).map(l -> l.split("\t"))
			.filter(t -> t.length > 2).filter(t -> !t[0].equals("Row ID")).filter(t -> Double.valueOf(t[2]) <= pvalThr)
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));
	}
	private static Set<String> readMeasuredIDs(String dir) throws IOException
	{
		return Files.lines(Paths.get(dir + "value-changes.txt")).map(l -> l.split("\t"))
			.filter(t -> t.length > 2).filter(t -> !t[0].equals("Row ID")).filter(t -> !t[2].equals("NaN"))
			.map(t -> t[0]).collect(Collectors.toSet());
	}

	private static void exploreResultOverlaps() throws IOException
	{
		String base = DIR + "clusters/against-normals/";

		File[] d = new File(base).listFiles();

		List<String> dirs = Arrays.stream(d).map(f -> f.toString().substring(f.toString().lastIndexOf("/") + 1)).filter(s -> !s.startsWith(".")).collect(Collectors.toList());

		List<Set<String>> relsList = dirs.stream().map(dir -> readSIFRelations(base + dir)).collect(Collectors.toList());

		CollectionUtil.printNameMapping(dirs.toArray(new String[0]));
		CollectionUtil.printVennCounts(relsList.toArray(new Collection[0]));

		Set<String> intersect = new HashSet<>(relsList.get(0));
		relsList.forEach(l -> intersect.retainAll(l));

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "temp/common.sif"));
		intersect.forEach(rel -> FileUtil.writeln(rel, writer));
		writer.close();

	}

	private static Set<String> readSIFRelations(String dir)
	{
		try
		{
			return Files.lines(Paths.get(dir + "/causative.sif")).map(l -> l.split("\t"))
				.filter(t -> t.length > 2).map(t -> t[0] + "\t" + t[1] + "\t" + t[2] + (t.length > 4 ? "\t\t" + t[4] : "")).collect(Collectors.toSet());
		}
		catch (IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}

	final static String DEFINERS = "NELFB\n" +
		"SUPT16H\n" +
		"CDK12\n" +
		"SUPT5H\n" +
		"GTF2F1\n" +
		"FANCI\n" +
		"POLR2K\n" +
		"ELOC\n" +
		"ELOB\n" +
		"CCNT2\n" +
		"NELFA\n" +
		"CTDP1\n" +
		"GTF2H5\n" +
		"GTF2H4\n" +
		"GTF2H3\n" +
		"MNAT1\n" +
		"GTF2H1\n" +
		"CCNH\n" +
		"SSRP1\n" +
		"POLR2F\n" +
		"GTF2F2\n" +
		"MDC1\n" +
		"POLR2I\n" +
		"POLR2H\n" +
		"POLR2E\n" +
		"POLR2J\n" +
		"BRCA1\n" +
		"ELOA\n" +
		"ELL\n" +
		"POLR2C\n" +
		"TCEA1\n" +
		"ERCC3\n" +
		"SUPT4H1\n" +
		"POLR2L\n" +
		"CCNT1\n" +
		"NELFE\n" +
		"POLR2A\n" +
		"CCNK\n" +
		"\n" +
		"E2F3\n" +
		"TFDP2\n" +
		"CCNB1\n" +
		"MYC\n" +
		"CDK1\n" +
		"TFDP1\n" +
		"RBL2\n" +
		"RBL1\n" +
		"CDK2\n" +
		"E2F4\n" +
		"MAX\n" +
		"E2F1\n" +
		"\n" +
		"RUVBL2\n" +
		"TP53\n" +
		"MYBL2\n" +
		"FOXM1\n" +
		"SP2\n" +
		"JUN\n" +
		"TRRAP\n";

	final static Set<String> COSMIC = new HashSet<>(Arrays.asList(("ARAF\n" +
		"BIRC6\n" +
		"CPEB3\n" +
		"CSMD3\n" +
		"CUL3\n" +
		"EED\n" +
		"EGFR\n" +
		"EPHA3\n" +
		"GPC5\n" +
		"GRIN2A\n" +
		"HIF1A\n" +
		"KRAS\n" +
		"LEPROTL1\n" +
		"MALAT1\n" +
		"MB21D2\n" +
		"MYCL\n" +
		"N4BP2\n" +
		"NOTCH1\n" +
		"PTPN13\n" +
		"PTPRD\n" +
		"PTPRT\n" +
		"RAD21\n" +
		"RB1\n" +
		"RBM10\n" +
		"RFWD3\n" +
		"SIRPA\n" +
		"STRN\n" +
		"TP53\n" +
		"USP44\n" +
		"ZNF479\n" +
		"AKT1\n" +
		"ALK\n" +
		"BAP1\n" +
		"BRAF\n" +
		"CCDC6\n" +
		"CD74\n" +
		"DDR2\n" +
		"DROSHA\n" +
		"EGFR\n" +
		"EML4\n" +
		"ERBB2\n" +
		"ERBB4\n" +
		"EZR\n" +
		"FGFR2\n" +
		"HIP1\n" +
		"KDR\n" +
		"KEAP1\n" +
		"KIF5B\n" +
		"LRIG3\n" +
		"MAP2K1\n" +
		"MAP2K2\n" +
		"NFE2L2\n" +
		"NKX2-1\n" +
		"NRG1\n" +
		"PIK3CB\n" +
		"PTPN13\n" +
		"RET\n" +
		"ROS1\n" +
		"SDC4\n" +
		"SLC34A2\n" +
		"SMARCA4\n" +
		"SOX2\n" +
		"STK11\n" +
		"TFG\n" +
		"TP63\n" +
		"TPM3\n" +
		"TPR\n" +
		"RASA1\n" +
		"KLF5\n" +
		"EP300\n" +
		"CREBBP\n" +
		"MAPK1").split("\n")));
}
