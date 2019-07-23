package org.panda.misc.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.MGIVertebrateHomology;
import org.panda.resource.SiteMappingRatToHuman;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class ShowLing
{
	private static final double LOG2 = Math.log(2);

	public static void main(String[] args) throws IOException
	{
		String dataDir = "/home/ozgun/Data/ShowLing/";
		String outDir = "/home/ozgun/Analyses/ShowLing/";

		convertFile(dataDir + "Exp1.csv", outDir + "Exp1/data.txt");
		convertFile(dataDir + "Exp2.csv", outDir + "Exp2/data.txt");
		convertFile(dataDir + "Exp3.csv", outDir + "Exp3/data.txt");
		convertFile(dataDir + "Exp4.csv", outDir + "Exp4/data.txt");
	}

	private static void convertFile(String inFile, String outFile) throws IOException
	{
		String[] header = Files.lines(Paths.get(inFile)).skip(8).findFirst().get().split("\t");

		int symbolsInd = 2;
		int upInd = 5;
		int sitesInd = 35;
		int sampleStart = 41;
		int sampleEnd = 50;

		Set<String> IDMemory = new HashSet<>();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect");

		for (int i = sampleStart; i <= sampleEnd; i++)
		{
			writer.write("\t" + header[i]);
		}

		Files.lines(Paths.get(inFile)).skip(9).map(l -> l.split("\t")).forEach(t ->
		{
			// we don't expect missing values
			if (t.length < sampleEnd || t[sampleEnd].isEmpty()) return;

			String[] syms = t[symbolsInd].split("; ");
			String[] ups = parseUP(t[upInd]);
			String[] sites = t[sitesInd].split("; ");

			if (syms[0].isEmpty()) return;
			if (sites[0].isEmpty()) return;

			// look for human mapping

			Set<String> hSyms = MGIVertebrateHomology.get().getCorrespondingHumanSymbols(syms[0], MGIVertebrateHomology.Organism.RAT);

			if (hSyms.isEmpty()) return;

			Map<String, List<String>> mapping = SiteMappingRatToHuman.get().mapToHumanSite(ups[0], sites);

			if (mapping.isEmpty()) return;

			List<String> symbolList = new ArrayList<>();
			List<String> siteList = new ArrayList<>();

			for (String hUP : mapping.keySet())
			{
				String hSym = HGNC.get().getSymbol(hUP);

				if (hSym == null || !hSyms.contains(hSym)) continue;

				symbolList.add(hSym);
				String siteStr = CollectionUtil.merge(mapping.get(hUP).stream()
					.sorted((s1, s2) -> new Integer(s1.substring(1)).compareTo(Integer.valueOf(s2.substring(1))))
					.collect(Collectors.toList()), "|");
				siteList.add(siteStr);
			}

			if (!symbolList.isEmpty())
			{
				assert symbolList.size() == siteList.size();
				String id = constructID(symbolList, siteList, IDMemory);
				String symStr = CollectionUtil.merge(symbolList, " ");
				String siteStr = CollectionUtil.merge(siteList, " ");
				FileUtil.lnwrite(id + "\t" + symStr + "\t" + siteStr + "\t", writer);

				for (int i = sampleStart; i <= sampleEnd; i++)
				{
					FileUtil.tab_write((Math.log(Double.valueOf(t[i].replaceAll(",", ""))) / LOG2), writer);
				}
			}
		});

		writer.close();
	}

	private static String[] parseUP(String value)
	{
		String[] t = value.split("; ");
		String[] up = new String[t.length];
		for (int i = 0; i < t.length; i++)
		{
			up[i] = t[i].split("\\|")[1];
		}
		return up;
	}

	private static String constructID(List<String> symbols, List<String> sites, Set<String> idMem)
	{
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < symbols.size(); i++)
		{
			sb.append(symbols.get(i)).append("_").append(sites.get(i).replaceAll("\\|", "-"));
			if (i < symbols.size() - 1) sb.append("_");
		}
		String id = sb.toString();

		if (idMem.contains(id))
		{
			System.err.println("Encountered the same ID = " + id);

			String alt = id;
			int i = 0;
			while (idMem.contains(alt)) alt = id + "-repeat-" + (++i);

			id = alt;
		}

		idMem.add(id);
		return id;
	}
}
