package org.panda.misc.analyses;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

public class PlateletIL6
{
	public static final String BASE = "/Users/ozgun/Documents/Analyses/Platelet-IL6/";

	public static void main(String[] args) throws IOException
	{
		convertData();
	}

	public static void convertData() throws IOException
	{
		Map<String, String> linesMap = new HashMap<>();
		Map<String, Double> pMap = new HashMap<>();

		String inFile = BASE + "Quant_Analysis_3v3_paired.csv";

		String[] header = Files.lines(Paths.get(inFile)).skip(5).findFirst().get().split("\t");
		int geneInd = ArrayUtil.indexOf(header, "UniProt Gene Name");
		int siteInd = ArrayUtil.indexOf(header, "Site List");
		int pInd = ArrayUtil.indexOf(header, "PValue");
		int fdrInd = ArrayUtil.indexOf(header, "FDR");
		int fcInd = ArrayUtil.indexOf(header, "FC");

		Files.lines(Paths.get(inFile)).skip(6).map(l -> l.split("\t")).filter(t -> t.length > geneInd).forEach(t ->
		{
			String sym = t[geneInd].split(" ")[0];

			String sites = t[siteInd].replaceAll("; ", "|");

			double p = Double.valueOf(t[pInd]);
			double f = Double.valueOf(t[fdrInd]);
			if (t[fcInd].startsWith("-"))
			{
				p = -p;
				f = -f;
			}

			String id = sym.replaceAll(" ", "-") + "-" + sites.replaceAll("\\|| ", "-");

			String line = id + "\t" + sym + "\t" + sites + "\t\t" + p + "\t" + f;

			if (!linesMap.containsKey(id) || Math.abs(pMap.get(id)) > Math.abs(p))
			{
				linesMap.put(id, line);
				pMap.put(id, p);
			}
		});
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(BASE + "data.txt"));
		writer.write("ID\tSymbols\tSites\tEffect\tSignedP\tSignedFDR");
		linesMap.keySet().stream().sorted(Comparator.comparingDouble(o -> Math.abs(pMap.get(o)))).forEach(id -> FileUtil.lnwrite(linesMap.get(id), writer));
		writer.close();
	}

	public static void convertDataFirstVersion_deleteIfNoLongerNeeded() throws IOException
	{
		Map<String, String> linesMap = new HashMap<>();
		Map<String, Double> pMap = new HashMap<>();

		Files.lines(Paths.get(BASE + "IL6-preresultsN3-forOzgun-2-12282020.csv")).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String sym = t[0];
			if (sym.contains(";")) return;

			String sites = t[1].replaceAll("; ", "|");

			String tempS = "";
			int sNum = sym.split(" ").length;
			for (int i = 0; i < sNum; i++)
			{
				tempS += " " + sites;
			}
			sites = tempS.trim();

			double p = Double.valueOf(t[9]);
			if (t[8].startsWith("-")) p = -p;

			String id = sym.replaceAll(" ", "-") + "-" + sites.replaceAll("\\|| ", "-");

			String line = id + "\t" + sym + "\t" + sites + "\t\t" + p;

			if (!linesMap.containsKey(id) || Math.abs(pMap.get(id)) > Math.abs(p))
			{
				linesMap.put(id, line);
				pMap.put(id, p);
			}
		});
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(BASE + "data.txt"));
		writer.write("ID\tSymbols\tSites\tEffect\tSignedP");
		linesMap.keySet().stream().sorted(Comparator.comparingDouble(o -> Math.abs(pMap.get(o)))).forEach(id -> FileUtil.lnwrite(linesMap.get(id), writer));
		writer.close();
	}
}
