package org.panda.misc;

import org.panda.utility.ValToColor;
import org.panda.utility.statistics.Summary;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * Created by babur on 2/23/16.
 */
public class PanCanMutexSIFGenerator
{
	static Set<String> avoid = new HashSet<>(Arrays.asList(
		"PanCan",
		"UVM-keep",
		"mutations-only",
		"whole",
//		"outliers-excluded",
//		"no-network",
		"PC2v7",
		"fries1K",
		"fries290K",
		"fries290K-leidos-PC2v7",
		"fries290K-PC2v7",
		"fries1K-PC2v7"
	));

	public static void main(String[] args) throws IOException
	{
		prepare("/home/babur/Documents/mutex/TCGA", "/home/babur/Documents/DARPA/BigMech/no-network.sif");
	}

	private static void prepare(String dir, String outSIFFile) throws IOException
	{
		double thr = 0.01;
		int freqThr = 1;
		Map<String, Integer> cnt = new HashMap<>();
		Map<String, Map<String, Integer>> linkMap = new HashMap<>();

		for (File f : new File(dir).listFiles())
		{
			if (avoid.contains(f.getName())) continue;

			if (f.isDirectory())
			{
				List<List<String>> results = new ArrayList<>();
				readRecursive(f, thr, results);
				updateGeneCounts(results, cnt);
				updateLinkCounts(results, linkMap);
			}
		}

		List<String> frequent = new ArrayList<>();
		for (String gene : cnt.keySet())
		{
			if (cnt.get(gene) >= freqThr) frequent.add(gene);
		}

		int maxNodeFreq = Summary.max(cnt.values());



		Set<String> edges = new HashSet<>();
		BufferedWriter writer = new BufferedWriter(new FileWriter(outSIFFile));

		int maxEdgeFreq = 0;
		for (String g1 : frequent)
		{
			for (String g2 : frequent)
			{
				if (g1.compareTo(g2) < 0)
				{
					Integer c = linkMap.get(g1).get(g2);
					if (c != null)
					{
						String edge = g1 + "\tinteracts-with\t" + g2;
						writer.write(edge + "\n");
						edges.add(edge.replaceAll("\t", " "));

						if (c > maxEdgeFreq) maxEdgeFreq = c;
					}
				}
			}
		}

		writer.close();

		ValToColor nodeColor = new ValToColor(new double[]{0, Math.log(maxNodeFreq)}, new Color[]{Color.LIGHT_GRAY, Color.BLACK});
		ValToColor edgeColor = new ValToColor(new double[]{0, Math.log(maxEdgeFreq)}, new Color[]{Color.LIGHT_GRAY, Color.BLACK});

		writer = new BufferedWriter(new FileWriter(outSIFFile.substring(0, outSIFFile.lastIndexOf(".")) + ".format"));

		writer.write("node\tall-nodes\tborderwidth\t3\n");
		writer.write("node\tall-nodes\tcolor\t255 255 255\n");

		for (String gene : frequent)
		{
			writer.write("node\t" + gene + "\tbordercolor\t" +
				nodeColor.getColorInString(Math.log(cnt.get(gene))) +"\n");
			writer.write("node\t" + gene + "\ttooltip\t" + cnt.get(gene) + "\n");
		}

		for (String edge : edges)
		{
			String g1 = edge.split(" ")[0];
			String g2 = edge.split(" ")[2];
			writer.write("edge\t" + edge + "\twidth\t" + (linkMap.get(g1).get(g2) >= freqThr ? "2" : 1) + "\n");
			writer.write("edge\t" + edge + "\tcolor\t" +
				edgeColor.getColorInString(Math.log(linkMap.get(g1).get(g2))) + "\n");
		}

		writer.close();
	}

	private static void readRecursive(File dir, double thr, List<List<String>> results) throws FileNotFoundException
	{
		List<List<String>> list = read(dir.getPath(), thr);
		results.addAll(list);

		for (File f : dir.listFiles())
		{
			if (avoid.contains(f.getName())) continue;

			if (f.isDirectory()) readRecursive(f, thr, results);
		}
	}

	private static void updateGeneCounts(List<List<String>> results, Map<String, Integer> cnt)
	{
		Set<String> genes = MutexReader_old.convertToSet(results);
		for (String gene : genes)
		{
			if (!cnt.containsKey(gene)) cnt.put(gene, 1);
			else cnt.put(gene, cnt.get(gene) + 1);
		}
	}

	private static void updateLinkCounts(List<List<String>> results, Map<String, Map<String, Integer>> linkMap)
	{
		Map<String, Set<String>> map = new HashMap<>();

		for (List<String> result : results)
		{
			for (String g1 : result)
			{
				for (String g2 : result)
				{
					if (!g1.equals(g2))
					{
						if (!map.containsKey(g1)) map.put(g1, new HashSet<>());
						map.get(g1).add(g2);
					}
				}
			}
		}

		for (String g1 : map.keySet())
		{
			if (!linkMap.containsKey(g1)) linkMap.put(g1, new HashMap<>());

			for (String g2 : map.get(g1))
			{
				if (!linkMap.get(g1).containsKey(g2)) linkMap.get(g1).put(g2, 1);
				else linkMap.get(g1).put(g2, linkMap.get(g1).get(g2) + 1);
			}
		}
	}

	private static List<List<String>> read(String dir, double thr) throws FileNotFoundException
	{
		File f = new File(dir + "/ranked-groups.txt");
		if (!f.exists()) return Collections.emptyList();

		return MutexReader_old.readMutexResults(dir, thr, false);
	}
}
