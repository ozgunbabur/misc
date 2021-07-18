package org.panda.misc.analyses;

import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.StreamDirection;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Bouhaddou
{
	private static final String BASE = "/Users/ozgun/Documents/Analyses/Bouhaddou/";
	private static double pThr = 0.1;
	public static void main(String[] args) throws IOException
	{
		convertData();
//		subgraph();
	}

	private static void convertData() throws IOException
	{
		Map<String[], Map<String, Double>> pp = processPhosphoproteomicData();
		Map<String, Map<String, Double>> tp = processProteomicData();

		List<String> samples = pp.values().iterator().next().keySet().stream().sorted().collect(Collectors.toList());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(BASE + "data.txt"));

		writer.write("ID\tSymbols\tSites\tEffect");
		samples.forEach(s -> FileUtil.tab_write(s, writer));

		pp.forEach((s, map) ->
		{
			FileUtil.lnwrite(ArrayUtil.merge("\t", s), writer);
			samples.forEach(c -> FileUtil.tab_write(map.get(c), writer));
		});

		tp.forEach((s, map) ->
		{
			FileUtil.lnwrite(s.replaceAll(" ", "_") + "\t" + s + "\t\t", writer);
			samples.forEach(c -> FileUtil.tab_write(map.get(c), writer));
		});

		writer.close();
	}


	private static Map<String, Map<String, Double>> processProteomicData() throws IOException
	{
		String inFile = BASE + "totprot.csv";

		Map<String, Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");
		int ctrlPInd = ArrayUtil.indexOf(header, "Ctrl_24Hr.adj.pvalue");
		int ctrlFCInd = ArrayUtil.indexOf(header, "Ctrl_24Hr.log2FC");
		int symInd = ArrayUtil.indexOf(header, "Gene_Name");
		int valStartInd = ArrayUtil.indexOf(header, "Inf_00Hr.adj.pvalue");

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).filter(t -> !t[symInd].isEmpty()).forEach(t ->
		{
			// If the control case is missing then skip it
			if (t[ctrlPInd].isEmpty()) return;

			if (t[symInd].isEmpty()) return;

			// Determine the change direction of the control case
			int ctrlChange = 0;
			if (Double.valueOf(t[ctrlPInd]) <= pThr)
			{
				ctrlChange = t[ctrlFCInd].startsWith("-") ? -1 : 1;
			}

			List<String> geneList = Arrays.stream(t[symInd].split(",|;")).filter(l -> !l.equals("NA")).distinct().sorted().collect(Collectors.toList());

			String sym = CollectionUtil.merge(geneList, " ");

			HashMap<String, Double> valMap = new HashMap<>();

			for (int i = valStartInd; i < valStartInd + 6; i++)
			{
				double sp = t[i].isEmpty() ? Double.NaN : Double.valueOf(t[i]);
				if (sp == 0) sp = 1E-10;

				if (((int) Math.signum(sp)) == ctrlChange) sp = 1;

				if (t[i-7].startsWith("-")) sp = -sp;
				valMap.put(header[i].substring(0, header[i].indexOf(".")), sp);
			}

			map.put(sym, valMap);

		});

		return map;
	}

	private static Map<String[], Map<String, Double>> processPhosphoproteomicData() throws IOException
	{
		String inFile = BASE + "phospho.csv";

		Map<String[], Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");
		int ctrlPInd = ArrayUtil.indexOf(header, "Ctrl_24Hr.adj.pvalue");
		int ctrlFCInd = ArrayUtil.indexOf(header, "Ctrl_24Hr.log2FC");
		int symInd = ArrayUtil.indexOf(header, "Gene_Name");
		int siteInd = ArrayUtil.indexOf(header, "H.s.Ph");
		int valStartInd = ArrayUtil.indexOf(header, "Inf_00Hr.adj.pvalue");


		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).filter(t -> !t[symInd].isEmpty()).forEach(t ->
		{
			// If the control case is missing then skip it
			if (t[ctrlPInd].isEmpty()) return;

			if (t[symInd].isEmpty()) return;

			// Determine the change direction of the control case
			int ctrlChange = 0;
			if (Double.valueOf(t[ctrlPInd]) <= pThr)
			{
				ctrlChange = t[ctrlFCInd].startsWith("-") ? -1 : 1;
			}

			String sym = t[symInd];

			List<String> sites = Arrays.stream(t[siteInd].split(";")).map(s -> s.substring(s.indexOf("_") + 1)).collect(Collectors.toList());
			String siteStr = CollectionUtil.merge(sites, "|");

			String[] fourCols = {null, sym, siteStr, ""};
			setID(fourCols);


			HashMap<String, Double> valMap = new HashMap<>();

			for (int i = valStartInd; i < valStartInd + 6; i++)
			{
				double sp = t[i].isEmpty() ? Double.NaN : Double.valueOf(t[i]);
				if (sp == 0) sp = 1E-10;

				if (((int) Math.signum(sp)) == ctrlChange) sp = 1;

				if (t[i-7].startsWith("-")) sp = -sp;
				valMap.put(header[i].substring(0, header[i].indexOf(".")), sp);
			}


			map.put(fourCols, valMap);
		});

		return map;
	}

	private static void setID(String[] fourCols)
	{
		String[] genes = fourCols[1].split(" ");
		String[] sites = fourCols[2].split(" ");

		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < genes.length; i++)
		{
			String s = genes[i] + "-" + sites[i].replaceAll("\\|", "-");
			if (sb.length() > 0) sb.append("_");
			sb.append(s);
		}
		fourCols[0] = sb.toString();
	}


	private static void subgraph() throws IOException
	{
		String dir = BASE + "fdr0.1/";
//		Set<String> seed = SIFFileUtil.getNeighborNodes(dir + "causative.sif", Collections.singleton("PDPK1"), StreamDirection.BOTHSTREAM);
		Set<String> seed = new HashSet<>();
		seed.add("PDPK1");
		CausalPathSubnetwork.writeGOINeighForCompBased(dir, seed, StreamDirection.BOTHSTREAM, "focus-PDPK1-2-neigh");
	}
}
