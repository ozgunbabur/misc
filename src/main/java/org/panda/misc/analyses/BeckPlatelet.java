package org.panda.misc.analyses;

import org.panda.resource.HGNC;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.channels.WritableByteChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;

/**
 * For the analysis of the data int he paper:
 * https://www.ncbi.nlm.nih.gov/pubmed/28060719
 *
 * @author Ozgun Babur
 */
public class BeckPlatelet
{
	public static final String DIR = "/home/babur/Documents/Analyses/BeckProteomics/";

	public static void main(String[] args) throws IOException
	{
//		prepareDataFile();
//		listUnmatchingRelations();
		printDifferentialGeneNames();
	}

	public static void prepareDataFile() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "data.txt"));
		writer.write("ID\tSymbols\tSites\tEffect\t10s\t20s\t30s");
		Files.lines(Paths.get(DIR + "TableS3.txt")).skip(2).map(l -> l.split("\t")).filter(t -> !t[2].startsWith("n"))
			.forEach(t ->
		{
			String sym = HGNC.get().getSymbol(t[2]);

			if (sym == null) return;

			String[] sites = t[4].split(" or | and | and/or ");
			Double[] v = new Double[]{Double.valueOf(t[5]), Double.valueOf(t[8]), Double.valueOf(t[11])};
			for (int i = 0; i < v.length; i++)
			{
				if (v[i] < 1) v[i] = -(1 / v[i]);
			}

			List<String> siteList = new ArrayList<>();
			for (String site : sites)
			{
				siteList.add("S" + site);
				siteList.add("T" + site);
				siteList.add("Y" + site);
			}

			String id = sym;
			for (String site : sites)
			{
				id += "_" + site;
			}

			FileUtil.lnwrite(id + "\t" + sym + "\t" + CollectionUtil.merge(siteList, "|") + "\t\t" +
				ArrayUtil.getString("\t", v), writer);

		});

		writer.close();
	}

	private static void listUnmatchingRelations() throws IOException
	{
		String dir = DIR + "run-no-site-match/";
		Map<String, Set<Integer>> observedSites = new HashMap<>();
		Files.lines(Paths.get(dir + "value-changes.txt")).map(l -> l.split("\t")[0]).forEach(s ->
		{
			String[] t = s.split("_");
			if (!observedSites.containsKey(t[0])) observedSites.put(t[0], new HashSet<>());
			for (int i = 1; i < t.length; i++)
			{
				observedSites.get(t[0]).add(Integer.valueOf(t[i]));
			}
		});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "unmatching.txt"));
		writer.write("Regulator\tRelation\tTarget\tAffected sites\tObserved sites");

		Files.lines(Paths.get(dir + "causative.sif")).sorted().map(l -> l.split("\t")).forEach(t ->
		{
			List<Integer> sites = new ArrayList<>(observedSites.get(t[2]));
			Collections.sort(sites);
			FileUtil.lnwrite(ArrayUtil.getString("\t", t[0], t[1], t[2]) + "\t" + (t.length > 4 ? t[4] : "") + "\t" + ArrayUtil.getString(";", sites.toArray()), writer);
		});

		writer.close();
	}

	private static void printDifferentialGeneNames() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "for-ppi.format"));
		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});
		Files.lines(Paths.get(DIR + "data.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			Double v = Double.valueOf(t[4]);
			if (Math.abs(Double.valueOf(t[5])) > Math.abs(v)) v = Double.valueOf(t[5]);
			if (Math.abs(Double.valueOf(t[6])) > Math.abs(v)) v = Double.valueOf(t[6]);

			FileUtil.writeln("node\t" + t[1] + "\tcolor\t" + vtc.getColorInString(v), writer);
		});
		writer.close();
	}
}
