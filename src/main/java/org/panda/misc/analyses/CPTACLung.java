package org.panda.misc.analyses;

import org.panda.misc.CausalPathSubnetwork;
import org.panda.utility.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class CPTACLung
{
	public static final String DIR = "/Users/ozgun/Documents/Analyses/CPTAC-LUAD/";

	public static void main(String[] args) throws IOException
	{
//		convertData();
//		printColNames();
//		prepareCosmicGraph();
		generateSigNeighGraphs();
//		generatePartialSigGraph();
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
		String inFile = DIR + "luad-v3.1-proteome-ratio-norm-NArm.gct";

		Map<String, String[]> rowMap = new HashMap<>();
		Map<String, Integer> specMap = new HashMap<>();

		Files.lines(Paths.get(inFile)).skip(81).map(l -> l.split("\t")).forEach(t ->
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
			for (int i = 16; i < header.length; i++)
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
		String inFile = DIR + "luad-v3.1-phosphoproteome-ratio-norm-NArm.txt";

		Map<String[], Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(inFile)).skip(2).findFirst().get().split("\t");

		Files.lines(Paths.get(inFile)).skip(81).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = t[2];
			List<String> sites = new ArrayList<>();
			for (String s : t[9].split(" "))
			{
				if (s.startsWith("S") || s.startsWith("T") || s.startsWith("Y"))
				{
					sites.add(s.substring(0, s.length() - 1));
				}
			}

			HashMap<String, Double> valMap = new HashMap<>();
			for (int i = 22; i < header.length; i++)
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
		String inFile = DIR + "luad-v3.1-proteome-ratio-norm-NArm.gct";
		String[] header = Files.lines(Paths.get(inFile)).skip(2).findFirst().get().split("\t");
		String[] type = Files.lines(Paths.get(inFile)).skip(8).findFirst().get().split("\t");
		String[] cluster = Files.lines(Paths.get(inFile)).skip(29).findFirst().get().split("\t");

		// For correlation
//		for (int i = 16; i < header.length; i++)
//		{
//			if (type[i].equals("Tumor"))
//			{
//				System.out.println("value-column = " + header[i]);
//			}
//		}

		// Tumor vs normal
//		for (int i = 16; i < header.length; i++)
//		{
//			if (type[i].equals("Tumor"))
//			{
//				System.out.println("test-value-column = " + header[i]);
//			}
//			else if (type[i].equals("NAT"))
//			{
//				System.out.println("control-value-column = " + header[i]);
//			}
//		}

		// Clusters

		List<String> clusters = Arrays.stream(cluster).filter(c -> c.startsWith("C")).distinct()
			.sorted().collect(Collectors.toList());
		for (int i = 0; i < clusters.size(); i++)
		{
			String testC = clusters.get(i);
			String dir;

//			dir = DIR + "clusters/against-others/" + testC + "/";
//			Files.createDirectories(Paths.get(dir));
//			BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
//			Files.lines(Paths.get(DIR + "clusters/parameters.template.txt")).forEach(l -> FileUtil.writeln(l, writer1));
//			for (int j = 16; j < header.length; j++)
//			{
//				if (cluster[j].equals(testC)) FileUtil.lnwrite("test-value-column = " + header[j], writer1);
//				else if (cluster[j].startsWith("C")) FileUtil.lnwrite("control-value-column = " + header[j], writer1);
//			}
//			writer1.close();

			dir = DIR + "clusters/against-normals/" + testC + "/";
			Files.createDirectories(Paths.get(dir));
			BufferedWriter writer3 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
			Files.lines(Paths.get(DIR + "clusters/parameters.template.txt")).forEach(l -> FileUtil.writeln(l, writer3));
			Set<String> tumors = new HashSet<>();
			Set<String> normals = new HashSet<>();
			for (int j = 16; j < header.length; j++)
			{
				if (cluster[j].equals(testC)) tumors.add(header[j]);
			}
			for (int j = 16; j < header.length; j++)
			{
				if (header[j].endsWith(".N") && tumors.contains(header[j].substring(0, header[j].length() - 2))) normals.add(header[j]);
			}
			tumors.stream().sorted().forEach(s -> FileUtil.lnwrite("test-value-column = " + s, writer3));
			normals.stream().sorted().forEach(s -> FileUtil.lnwrite("control-value-column = " + s, writer3));
			writer3.close();

//			for (int j = i+1; j < clusters.size(); j++)
//			{
//				String ctrlC = clusters.get(j);
//				dir = DIR + "clusters/pairs/" + testC + "-vs-" + ctrlC + "/";
//				Files.createDirectories(Paths.get(dir));
//				BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(dir + "parameters.txt"));
//				Files.lines(Paths.get(DIR + "clusters/parameters.template.txt")).forEach(l -> FileUtil.writeln(l, writer2));
//				for (int k = 16; k < header.length; k++)
//				{
//					if (cluster[k].equals(testC)) FileUtil.lnwrite("test-value-column = " + header[k], writer2);
//					else if (cluster[k].equals(ctrlC)) FileUtil.lnwrite("control-value-column = " + header[k], writer2);
//				}
//				writer2.close();
//			}
		}

		// Mutations

//		Set<String> oncogenes = new HashSet<>(Arrays.asList("KRAS", "EGFR", "IL21R", "EGFL6", "LMO2", "BIRC6", "BRAF", "ARAF", "ERBB2"));
//
//		Files.lines(Paths.get(inFile)).map(l -> l.split("\t")).filter(t -> t[0].endsWith(".mutation.status")).forEach(t ->
//		{
//			String gene = t[0].substring(0, t[0].indexOf("."));
//			String dir = DIR + "mutations/" + gene + "/";
//			FileUtil.mkdirs(dir);
//			BufferedWriter writer = FileUtil.newBufferedWriter(dir + "parameters.txt");
//			FileUtil.lines(DIR + "mutations/parameters.template.txt").forEach(l -> FileUtil.writeln(l, writer));
//
//			FileUtil.lnwrite("gene-activity = " + gene + " " + (oncogenes.contains(gene) ? "a" : "i") + "\n", writer);
//
//			for (int i = 16; i < header.length; i++)
//			{
//				if (!type[i].equals("Tumor")) continue;
//
//				if (t[i].equals("1")) FileUtil.lnwrite("test-value-column = " + header[i], writer);
//				else if (t[i].equals("0")) FileUtil.lnwrite("control-value-column = " + header[i], writer);
//			}
//			FileUtil.closeWriter(writer);
//		});
	}

	private static void prepareCosmicGraph() throws IOException
	{
//		CausalPathSubnetwork.writeGOINeighForCorrBased(DIR + "correlation-tumor", COSMIC, StreamDirection.UPSTREAM);
//		CausalPathSubnetwork.writeGOINeighForCompBased(DIR + "tumor-vs-normal", COSMIC, StreamDirection.UPSTREAM);

//		CausalPathSubnetwork.writeGOINeighForCorrBased(DIR + "correlation-tumor", Collections.singleton("RUNX1"), StreamDirection.BOTHSTREAM);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "cosmic.highlight"));
		COSMIC.stream().sorted().forEach(g -> FileUtil.writeln("node\t" + g, writer));
		writer.close();
	}

	private static void generateSigNeighGraphs() throws IOException
	{
//		CausalPathSubnetwork.generateNeighborhoodSubgraphsForSignificantsRecursively(DIR + "correlation-tumor-rnaseq", 0.1);
		CausalPathSubnetwork.writeSignifNeighForCorrBased(DIR + "correlation-tumor-rnaseq", StreamDirection.DOWNSTREAM);
	}

	private static void generatePartialSigGraph() throws IOException
	{
		String dir = DIR + "clusters/against-normals/C3";
		Set<String> goi = CausalPathSubnetwork.combine(CausalPathSubnetwork.getSignificantGenes(dir + "/significance-pvals.txt", 0.1));
		goi.removeAll(Arrays.asList("CDK1", "CDK2", "RUNX1", "CBFB", "MYC", "MAX"));
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/causative-sig-partial-neigh.sif", StreamDirection.DOWNSTREAM);
		Set<String> ids = CausalPathSubnetwork.getIDsAtTheNeighborhood(dir + "/results.txt", goi, StreamDirection.DOWNSTREAM);
		System.out.println("ids.size() = " + ids.size());
		CausalPathSubnetwork.writeSubsetFormat(dir + "/causative.format", dir + "/causative-sig-partial-neigh.format", goi, ids);
	}

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
		"TPR").split("\n")));
}
