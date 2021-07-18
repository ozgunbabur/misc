package org.panda.misc.analyses;

import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.utility.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Klann
{
	public static final String BASE = "/Users/ozgun/Documents/Analyses/Klann/";
	public static void main(String[] args) throws IOException
	{
//		convertData();
		subgraph();
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
			FileUtil.lnwrite(s + "\t" + s + "\t\t", writer);
			samples.forEach(c -> FileUtil.tab_write(map.get(c), writer));
		});

		writer.close();
	}


	private static Map<String, Map<String, Double>> processProteomicData() throws IOException
	{
		String inFile = BASE + "totprot.csv";

		Map<String, Map<String, Double>> map = new HashMap<>();

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = t[1];
			if (sym.contains("GN="))
			{
				int beginIndex = sym.indexOf("GN=") + 3;
				sym = sym.substring(beginIndex, sym.indexOf(" ", beginIndex));
			}
			else return;

			HashMap<String, Double> valMap = new HashMap<>();

			double sp = Double.valueOf(t[t.length-1]);
			if (t[t.length-2].startsWith("-")) sp = -sp;
			valMap.put("SignedP", sp);

			map.put(sym, valMap);

		});

		return map;
	}

	private static Map<String[], Map<String, Double>> processPhosphoproteomicData() throws IOException
	{
		String inFile = BASE + "phospho.csv";

		Map<String[], Map<String, Double>> map = new HashMap<>();

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			List<List<String>> sitesList = parseSites(t[0], t[1]);

			String syms = t[2].replaceAll(";", " ");

			if (sitesList.size() != syms.split(" ").length)
			{
				System.out.println("Genes and sites don't match: " + Arrays.asList(t));
				return;
			}

			HashMap<String, Double> valMap = new HashMap<>();

			if (t[4].startsWith("#")) return;

			double sp = Double.valueOf(t[4]);
			if (t[3].startsWith("-")) sp = -sp;
			valMap.put("SignedP", sp);

			List<String> collect = sitesList.stream().map(l -> CollectionUtil.merge(l, "|")).collect(Collectors.toList());
			String sites = CollectionUtil.merge(collect, " ");

			String[] fourCols = {null, syms, sites, ""};
			setID(fourCols);

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

	private static List<List<String>> parseSites(String ups, String cell)
	{
		String[] ids = ups.split(";");

		List<String> sitesList = new ArrayList<>();

		for (int i = 0; i < ids.length; i++)
		{
			int beginIndex = cell.indexOf(ids[i]);
			int endIndex = i < ids.length - 1 ? cell.indexOf(ids[i+1]) : cell.length();
			String sub = cell.substring(beginIndex, endIndex);
			sitesList.add(sub);
		}

		return sitesList.stream().map(Klann::parseOneProteinSites).collect(Collectors.toList());
	}

	private static List<String> parseOneProteinSites(String string)
	{
		List<String> sites = new ArrayList<>();

		int beginIndex = string.indexOf("Phospho [") + 9;
		if (beginIndex > 9)
		{
			string = string.substring(beginIndex, string.indexOf("]", beginIndex));

			String[] t = string.split(";");

			for (int i = 0; i < t.length; i++)
			{
				if (t[i].contains("(")) sites.add(t[i].substring(0, t[i].indexOf("(")));
			}
		}

		return sites;
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
