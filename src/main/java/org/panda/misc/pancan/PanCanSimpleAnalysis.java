package org.panda.misc.pancan;

import org.panda.utility.Progress;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * For finding genes that are more likely to be mutated in samples with fewer alterations.
 * Created by babur on 3/8/16.
 */
public class PanCanSimpleAnalysis
{
	public static final int ITERATION = 100000;

	public static void main(String[] args) throws FileNotFoundException
	{
		findGenesOverlapWithSilentSamples();
	}

	private static void findGenesOverlapWithSilentSamples() throws FileNotFoundException
	{
		Map<String, boolean[]> data = readData();
		System.out.println("Number of genes: " + data.keySet().size());
		int[] freq = calcAltFreqs(data);

		Map<String, Double> pvals = new HashMap<>();

		Map<Integer, List<Integer>> randomTrialMap = new HashMap<>();

		Progress prg = new Progress(data.keySet().size(), "Calculating p-values");
		for (String gene : data.keySet())
		{
			pvals.put(gene, calcSilenceTendency(data.get(gene), freq, randomTrialMap));
			prg.tick();
		}

		List<String> genes = new ArrayList<>(pvals.keySet());
		Collections.sort(genes, (o1, o2) -> pvals.get(o1).compareTo(pvals.get(o2)));

		for (String gene : genes)
		{
			Double p = pvals.get(gene);
			if (p > 0.05) break;
			System.out.println(gene + "\t" + p);
		}
	}

	private static double calcSilenceTendency(boolean[] alts, int[] freq, Map<Integer, List<Integer>> randomTrialMap)
	{
		int sum = 0;
		int cnt = 0;
		for (int i = 0; i < alts.length; i++)
		{
			if (alts[i])
			{
				sum += freq[i];
				cnt++;
			}
		}

		if (!randomTrialMap.containsKey(cnt)) randomTrialMap.put(cnt, generateRandomDistribution(freq, cnt));

		int pass = 0;
		for (Integer val : randomTrialMap.get(cnt))
		{
			if (val <= sum) pass++;
			else break;
		}

		return pass / (double) ITERATION;
	}

	private static List<Integer> generateRandomDistribution(int[] freq, int altCount)
	{
		List<Integer> vals = new ArrayList<>();
		List<Integer> list = new ArrayList<>();
		for (int i : freq)
		{
			list.add(i);
		}

		for (int i = 0; i < ITERATION; i++)
		{
			Collections.shuffle(list);
			int sum = 0;
			for (int j = 0; j < altCount; j++)
			{
				sum += list.get(j);
			}
			vals.add(sum);
		}
		Collections.sort(vals);
		return vals;
	}

	private static int[] calcAltFreqs(Map<String, boolean[]> dataMap)
	{
		int[] f = new int[dataMap.values().iterator().next().length];

		for (boolean[] b : dataMap.values())
		{
			for (int i = 0; i < b.length; i++)
			{
				if (b[i]) f[i]++;
			}
		}
		return f;
	}

	private static Map<String, boolean[]> readData() throws FileNotFoundException
	{
		Map<String, boolean[]> dataMap = new HashMap<>();
		Scanner sc = new Scanner(new File("/home/babur/Documents/mutex/TCGA/PanCan/mutations-only/whole/DataMatrix" +
			".txt"));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			String gene = token[0];
			boolean[] b = new boolean[token.length - 1];
			for (int i = 0; i < b.length; i++)
			{
				b[i] = !token[i + 1].equals("0");
			}
			dataMap.put(gene, b);
		}
		sc.close();
		return dataMap;
	}
}
