package org.panda.misc.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.network.UniProt;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Fabian
{
	public static final String BASE = "/Users/ozgun/Documents/Analyses/Fabian/";
	public static void main(String[] args) throws IOException
	{
//		convertData();
		addModeratedValues();
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
			FileUtil.lnwrite(s.replaceAll(" ", "-") + "\t" + s + "\t\t", writer);
			samples.forEach(c -> FileUtil.tab_write(map.get(c), writer));
		});

		writer.close();
	}


	private static Map<String, Map<String, Double>> processProteomicData() throws IOException
	{
		String inFile = BASE + "20190313_230046_20190313_SanderCollab_Proteome_Alex.txt";

		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");
		int geneInd = ArrayUtil.indexOf(header, "PG.Genes");


		Map<String, Map<String, Double>> map = new HashMap<>();

		Files.lines(Paths.get(inFile)).skip(1).filter(l -> !l.isEmpty()).map(l -> l.split("\t")).forEach(t ->
		{
			if (t.length <= geneInd) return;
			List<String> genes = Arrays.stream(t[geneInd].split(";")).distinct().filter(g -> !g.isEmpty()).collect(Collectors.toList());

			HashMap<String, Double> valMap = new HashMap<>();

			for (int i = 0; i < 12; i++)
			{
				valMap.put(header[i], Double.valueOf(t[i]));
			}

			map.put(CollectionUtil.merge(genes, " "), valMap);

		});

		return map;
	}

	private static Map<String[], Map<String, Double>> processPhosphoproteomicData() throws IOException
	{
		String inFile = BASE + "phospho.csv";

		String[] sampleName = new String[]{"Ctrl_1", "Drug1_1", "Drug2_1", "Drug3_1", "Ctrl_2", "Drug1_2", "Drug2_2", "Drug3_2", "Ctrl_3", "Drug1_3", "Drug2_3", "Drug3_3"};

		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");
		int upInd = ArrayUtil.indexOf(header, "T: Proteins");
		int siteInd = ArrayUtil.indexOf(header, "T: Positions within proteins");
		int aaInd = ArrayUtil.indexOf(header, "C: Amino acid");

		Map<String[], Map<String, Double>> map = new HashMap<>();

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String[] up = t[upInd].split(";");
			String[] sn = t[siteInd].split(";");

			if (up.length != sn.length)
			{
				System.err.println("UP IDs and site numbers don't match: " + Arrays.asList(t));
				return;
			}

			List<String> genes = new ArrayList<>();
			List<String> sites = new ArrayList<>();

			for (int i = 0; i < up.length; i++)
			{
				String gene = HGNC.get().getSymbol(up[i]);
				if (gene != null && !genes.contains(gene))
				{
					genes.add(gene);
					sites.add(t[aaInd] + sn[i]);
				}
			}

			if (genes.isEmpty())
			{
				System.err.println("Could ot recognize any gene: " + Arrays.asList(t));
				return;
			}

			HashMap<String, Double> valMap = new HashMap<>();


			for (int i = 0; i < 12; i++)
			{
				valMap.put(sampleName[i], Double.valueOf(t[i]));
			}


			String[] fourCols = {null, CollectionUtil.merge(genes, " "), CollectionUtil.merge(sites, " "), ""};
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


	private static void addModeratedValues() throws IOException
	{
		Map<String, String> mapD1 = readModeratedValues(BASE + "D1.csv");
		Map<String, String> mapD2 = readModeratedValues(BASE + "D2.csv");
		Map<String, String> mapD3 = readModeratedValues(BASE + "D3.csv");

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(BASE + "data-moderated.txt"));

		Files.lines(Paths.get(BASE + "data.txt")).forEach(l ->
		{
			String[] t = l.split("\t");


			if (t[0].equals("ID")) FileUtil.write(l + "\tD1\tD2\tD3", writer);
			else FileUtil.lnwrite(l + "\t" + mapD1.get(t[0]) + "\t" + mapD2.get(t[0]) + "\t" + mapD3.get(t[0]), writer);
		});

		writer.close();
	}

	private static Map<String, String> readModeratedValues(String file) throws IOException
	{
		return FileUtil.readMap(file, ",", "ID", "SignedP");
	}
}
