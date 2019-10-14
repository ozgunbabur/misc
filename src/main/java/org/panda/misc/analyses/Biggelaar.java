package org.panda.misc.analyses;

import ch.qos.logback.core.rolling.helper.FileStoreUtil;
import org.panda.resource.HGNC;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class Biggelaar
{
	public static final String DIR = "/home/ozgun/Analyses/Biggelaar/";

	public static void main(String[] args) throws IOException
	{
//		test();
		convertDataset();
	}

	public static void convertDataset() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "data.txt"));
		writer.write("ID\tSymbols\tSites\tEffect");

		String[] header = Files.lines(Paths.get(DIR + "table1.txt")).findFirst().get().split("\t");
		Map<String, int[]> colGroups = getColGroups(header);
		List<String> cols = new ArrayList<>(colGroups.keySet());
		cols.sort(String::compareTo);

		cols.forEach(c -> FileUtil.tab_write(c, writer));

		Files.lines(Paths.get(DIR + "table1.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String UPs = t[1];
			String sym = t[3].split(";")[0];
			int uInd = getRelatedUniprotIndex(sym, UPs);
			if (uInd < 0) return;

			String pep = t[7];
			String aa = t[4];

			int firstPos = Integer.valueOf(t[5].split(";")[uInd]);

			List<String> positions = getPositions(pep, aa, firstPos);

			String id = sym + "_" + CollectionUtil.merge(positions, "_");

			FileUtil.lnwrite(id + "\t" + sym + "\t" + CollectionUtil.merge(positions, "|") + "\t", writer);

			cols.forEach(col -> FileUtil.tab_write(getMostDistant(t, colGroups.get(col)), writer));
		});

		writer.close();

		cols.forEach(col -> FileUtil.mkdirs(DIR + col));
		for (String col : cols)
		{
			BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(DIR + col + "/parameters.txt"));

			writer2.write("proteomics-values-file = ../data.txt\n" +
				"id-column = ID\n" +
				"symbols-column = Symbols\n" +
				"sites-column = Sites\n" +
				"effect-column = Effect\n" +
				"value-transformation = arithmetic-mean\n" +
				"value-column = " + col + "\n" +
				"threshold-for-data-significance = 0.75 phosphoprotein\n" +
				"calculate-network-significance = true\n" +
				"permutations-for-significance = 10000\n" +
				"fdr-threshold-for-network-significance = 0.1\n" +
				"use-network-significance-for-causal-reasoning = true\n" +
				"minimum-potential-targets-to-consider-for-downstream-significance = 5\n" +
				"color-saturation-value = 3\n" +
				"do-site-matching = true\n" +
				"show-all-genes-with-proteomic-data = true\n");

			writer2.close();
		}
	}

	private static double getMostDistant(String[] row, int[] inds)
	{
		double v = Double.NaN;

		for (int i = 0; i < inds.length; i++)
		{
			double d = Double.valueOf(row[inds[i]]);
			if (!Double.isNaN(d))
			{
				if (Double.isNaN(v)) v = d;
				else if (Math.abs(d) > Math.abs(v)) v = d;
			}
		}

		return v;
	}

	private static int getRelatedUniprotIndex(String sym, String UPs)
	{
		String up = HGNC.get().getUniProt(sym);
		if (up != null)
		{
			String[] upArray = UPs.split(";");
			return ArrayUtil.indexOf(upArray, up);
		}
		return -1;
	}

	private static List<String> getPositions(String markedPep, String aa, int first)
	{
		markedPep = markedPep.replaceAll("\\(ox\\)", "");
		markedPep = markedPep.replaceAll("\\(ac\\)", "");
		markedPep = markedPep.replaceAll("_", "");

		String[] t = markedPep.split("\\(ph\\)");

		List<String> pos = new ArrayList<>();

		int ind = 0;

		for (int i = 0; i < t.length - 1; i++)
		{
			if (t[i].endsWith(aa))
			{
				ind = i;
				break;
			}
		}

		pos.add(aa + first);
		int loc = first;

		for (int i = ind + 1; i < t.length - 1; i++)
		{
			loc += t[i].length();
			if (t[i].endsWith(aa)) pos.add(aa + loc);
		}

		return pos;
	}

	private static Map<String, int[]> getColGroups(String[] header)
	{
		Map<String, int[]> map = new HashMap<>();

		for (int i = 0; i < header.length; i++)
		{
			if (header[i].startsWith("Ratio "))
			{
				String time = header[i].substring(header[i].indexOf("(") + 1, header[i].lastIndexOf(" "));
				String type = header[i].substring(header[i].indexOf(" ") + 1, header[i].indexOf(" norm")).replaceAll("/", "-");
				String col = time + "-" + type;

				if (!map.containsKey(col)) map.put(col, new int[]{i, 0, 0});
				else
				{
					int[] p = map.get(col);
					if (p[1] == 0) p[1] = i;
					else p[2] = i;
				}
			}
		}
		return map;
	}

	public static void test() throws IOException
	{
//		Set<String> items = new HashSet<>();
//		Files.lines(Paths.get(DIR + "table1.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
//		{
//			String pep = t[7];
//
//			while (pep.contains("("))
//			{
//				items.add(pep.substring(pep.indexOf("(") + 1, pep.indexOf(")")));
//				pep = pep.substring(pep.indexOf(")") + 1);
//			}
//		});
//
//		System.out.println("items = " + items);

		String up = HGNC.get().getUniProt("NCOA7");
		System.out.println("up = " + up);
	}
}
