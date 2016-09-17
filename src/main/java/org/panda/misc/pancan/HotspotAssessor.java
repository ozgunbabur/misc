package org.panda.misc.pancan;

import org.panda.resource.tcga.MutTuple;
import org.panda.resource.tcga.MutationReader;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.StringUtil;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class HotspotAssessor
{
	public static Map<String, Double> getHotspotPvalues(String maffile) throws IOException
	{
		MutationReader mr = new MutationReader(maffile);
		String[] samples = mr.getSamples().toArray(new String[0]);

		Map<String, Double> hotP = new HashMap<>();

		Set<String> genes = mr.getGenes();
		Progress p = new Progress(genes.size(), "Calculating hotspot p-values");
		for (String gene : genes)
		{
			p.tick();
			List<MutTuple>[] muts = mr.getMutations(gene, samples);

			if (muts != null)
			{
				Map<Integer, Integer> locCounts = getLocCounts(muts);
				if (locCounts.size() > 1)
				{
					double pval = calcHotspotPvalue(locCounts);
					hotP.put(gene, pval);
				}
			}
		}

		return hotP;
	}

	public static Map<String, List<Integer>> getHotspotLocs(String maffile, double pvalThr) throws IOException
	{
		MutationReader mr = new MutationReader(maffile);
		String[] samples = mr.getSamples().toArray(new String[0]);

		Map<String, List<Integer>> hotLocMap = new HashMap<>();

		Set<String> genes = mr.getGenes();
		Progress p = new Progress(genes.size(), "Calculating hotspot locations");
		for (String gene : genes)
		{
			p.tick();
			List<MutTuple>[] muts = mr.getMutations(gene, samples);

			if (muts != null)
			{
				Map<Integer, Integer> locCounts = getLocCounts(muts);
				if (locCounts.size() > 1)
				{
					if (gene.equals("APOB"))
					{
						System.out.print("");
					}
					int thr = calcHotspotThreshold(locCounts, pvalThr);
					List<Integer> locs = locCounts.keySet().stream().filter(l -> locCounts.get(l) >= thr).sorted()
						.collect(Collectors.toList());
					if (!locs.isEmpty())
					{
						hotLocMap.put(gene, locs);
					}
				}
			}
		}

		return hotLocMap;
	}



	private static Map<Integer, Integer> getLocCounts(List<MutTuple>[] muts)
	{
		Map<Integer, Integer> cnt = new HashMap<>();

		for (List<MutTuple> mut : muts)
		{
			if (mut == null) continue;

			for (MutTuple m : mut)
			{
				int loc = StringUtil.findFirstInt(m.value);
				if (loc < 0) continue;

				if (cnt.containsKey(loc)) cnt.put(loc, cnt.get(loc) + 1);
				else cnt.put(loc, 1);
			}
		}
		return cnt;
	}

	private static double calcHotspotPvalue(Map<Integer, Integer> cnt)
	{
		if (cnt.isEmpty()) return 1;

		int total = cnt.values().stream().mapToInt(Integer::valueOf).sum();
		int max = cnt.values().stream().mapToInt(Integer::valueOf).max().getAsInt();
		int distinct = cnt.size();

		return calcHotspotPvalue(distinct, total, max);
	}

	private static int calcHotspotThreshold(Map<Integer, Integer> cnt, double pvalThr)
	{
		if (cnt.isEmpty()) return Integer.MAX_VALUE;

		int total = cnt.values().stream().mapToInt(Integer::valueOf).sum();
		int distinct = cnt.size();

		return calcHotspotThreshold(distinct, total, pvalThr);
	}

	static Random rand = new Random();
	static final int maxTrials = 1000000;
	static final int maxHit = 100;

	private static int calcHotspotThreshold(int categories, int n, double pvalThr)
	{
		int trials = maxTrials / 100;
		List<Integer> list = new ArrayList<>(trials);
		for (int i = 0; i < trials; i++)
		{
			list.add(randomRun(categories, n));
		}
		Collections.sort(list);
		int index = (int) (list.size() * (1 - pvalThr));
		return list.get(index);
	}

	private static double calcHotspotPvalue(int categories, int n, int observedHigh)
	{
		int hit = 0;
		int trials = 0;

		while (hit < maxHit && trials < maxTrials)
		{
			if (randomRun(categories, n) >= observedHigh) hit++;
			trials++;
		}
		return hit / (double) trials;
	}

	private static int randomRun(int categories, int n)
	{
		int[] arr = new int[categories];
		for (int i = 0; i < n; i++)
		{
			arr[rand.nextInt(categories)]++;
		}
		return Summary.max(arr);
	}



	public static void main(String[] args) throws IOException
	{
//		Map<String, Double> pvals = getHotspotPvalues("/home/babur/Documents/TCGA/PanCan/mutation.maf");
//		BufferedWriter w1 = Files.newBufferedWriter(Paths.get("/home/babur/Documents/PanCan/hotspot-pvals.txt"));
//		w1.write("Gene\tP-val");
//		pvals.keySet().stream().sorted((g1, g2) -> pvals.get(g1).compareTo(pvals.get(g2))).forEach(g ->
//			FileUtil.lnwrite(g + "\t" + pvals.get(g), w1));
//		w1.close();

//		Map<String, List<Integer>> locMap = getHotspotLocs("/home/babur/Documents/TCGA/PanCan/mutation.maf", 0.05);
//		BufferedWriter w2 = Files.newBufferedWriter(Paths.get("/home/babur/Documents/PanCan/hotspot-locs.txt"));
//		locMap.keySet().stream().forEach(g ->
//		{
//			FileUtil.lnwrite(g, w2);
//			locMap.get(g).stream().forEach(l -> FileUtil.write("\t" + l, w2));
//		});
//		w2.close();

		System.out.println(calcHotspotPvalue(90, 99, 5));
	}
}
