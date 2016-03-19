package org.babur.misc;

import org.cbio.causality.idmapping.CancerGeneCensus;
import org.cbio.causality.util.DiscretePvalHisto;
import org.cbio.causality.util.FishersCombinedProbability;
import org.cbio.causality.util.Summary;

import java.io.*;
import java.util.*;

/**
 * Created by babur on 2/23/16.
 */
public class MutexResultAggregator
{
	public static void main(String[] args) throws IOException
	{
		prepare();
	}
	private static void prepare() throws IOException
	{
		Set<String> dirs = getDirs();
		Set<Map<String, Double>> results = readResults(dirs);
		Set<String> genes = getGenes(results);
		Map<String, Double> scores = calcAggregateScores(results, genes);
		List<String> geneList = new ArrayList<>(genes);
		Collections.sort(geneList, ((o1, o2) -> scores.get(o1).compareTo(scores.get(o2))));

		BufferedWriter writer = new BufferedWriter(new FileWriter("pancan-mut-scores.txt"));
		writer.write("Score\tGene");
		for (String gene : geneList)
		{
			writer.write("\n" + gene + "\t" + scores.get(gene) + "\t" + (CancerGeneCensus.isCancerGene(gene) ? "" :
				"X"));
		}
		writer.close();
		printHistogram(scores);
	}

	private static Set<String> getGenes(Set<Map<String, Double>> results)
	{
		Set<String> genes = new HashSet<>();
		for (Map<String, Double> map : results)
		{
			genes.addAll(map.keySet());
		}
		return genes;
	}

	private static Set<String> getDirs()
	{
		Set<String> dirs = new HashSet<>();
		for (int i = 1; i <= 10; i++)
		{
			dirs.add("/home/babur/Documents/mutex/TCGA/PanCan/mutations-only/" + i + "/no-network");
		}
		return dirs;
	}

	private static Set<Map<String, Double>> readResults(Set<String> dirs) throws FileNotFoundException
	{
		Set<Map<String, Double>> set = new HashSet<>();
		for (String dir : dirs)
		{
			Map<String, Double> scoreMap = readScores(dir);
			if (scoreMap != null) set.add(scoreMap);
		}
		return set;
	}

	private static Map<String, Double> calcAggregateScores(Set<Map<String, Double>> set, Set<String> genes)
	{
		Map<String, Double> scores = new HashMap<>();
		for (String gene : genes)
		{
			scores.put(gene, calcAggregateScore(set, gene));
		}
		return scores;
	}

	private static double calcAggregateScore(Set<Map<String, Double>> set, String gene)
	{
		List<Double> list = new ArrayList<>();

		for (Map<String, Double> map : set)
		{
			if (map.containsKey(gene)) list.add(map.get(gene));
		}

		double[] scores = new double[list.size()];
		for (int i = 0; i < scores.length; i++)
		{
			scores[i] = list.get(i);
		}

//		return FishersCombinedProbability.pValue(scores);
		return Summary.min(scores);
	}

	private static Map<String, Double> readScores(String dir) throws FileNotFoundException
	{
		Map<String, Double> map = new HashMap<>();
		File file = new File(dir + "/ranked-groups.txt");
		if (!file.exists())
		{
			System.err.println("File not found: " + file.getPath());
			return null;
		}

		Scanner sc = new Scanner(file);
		String header = sc.nextLine();
		boolean hasQval = header.split("\t").length > 2;
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");

			for (String gene : token[hasQval ? 2 : 1].split(" "))
			{
				if (!map.containsKey(gene)) map.put(gene, Double.parseDouble(token[0]));
			}
		}
		sc.close();
		return map;
	}

	private static void printHistogram(Map<String, Double> scores)
	{
		DiscretePvalHisto h = new DiscretePvalHisto(new ArrayList<>(scores.values()), 0.05);
		h.plot();
	}
}
