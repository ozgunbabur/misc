package org.panda.misc.analyses;

import org.panda.misc.MutexReader;
import org.panda.resource.CancerGeneBushman;
import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;

import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class BMEGSampleMutexRun
{
	public static void main(String[] args)
	{
		printResults();
	}

	public static void printResults()
	{
		String dir = "/home/babur/Documents/mutex/BMEG/TCGA";
		double scoreThrTight = 0.01;
		double scoreThrRelax = 0.01;
		int minTypeCnt = 1;

		Set<MutexReader.Group> groups = MutexReader.readMutexResultsRecursive(dir, new HashSet<>());

		Map<String, Double> geneScores = MutexReader.convertGroupsToGeneBestScores(groups);

		Map<String, List<String>> geneToType = new HashMap<>();
		Set<String> types = new HashSet<>();
		Map<String, List<String>> pairCounts = new HashMap<>();

		groups.stream().filter(g -> g.score < scoreThrRelax).sorted((g1, g2) -> new Double(g1.score).compareTo(g2.score)).forEach(g ->
		{
			String type = g.fromDir.split("/")[7];
			g.genes.stream().forEach(gene ->
			{
				if (!geneToType.containsKey(gene)) geneToType.put(gene, new ArrayList<>());
				if (!geneToType.get(gene).contains(type)) geneToType.get(gene).add(type);
				types.add(type);

				g.genes.stream().filter(g2 -> gene.compareTo(g2) < 0).forEach(gene2 ->
				{
					String key = gene + " " + gene2;
					if (!pairCounts.containsKey(key)) pairCounts.put(key, new ArrayList<>());
					if (!pairCounts.get(key).contains(type)) pairCounts.get(key).add(type);
				});
			});

			System.out.println(type + ": " + g.genes);
		});

//		Set<String> cancerGenes = new HashSet<>();
//		cancerGenes.addAll(CancerGeneCensus.get().getAllGenes());
//		cancerGenes.addAll(OncoKB.get().getAllGenes());
//		cancerGenes.addAll(CancerGeneBushman.get().getAllGenes());

		geneScores.keySet().stream().filter(k -> geneScores.get(k) < scoreThrRelax).filter(k -> geneToType.get(k).size() >= minTypeCnt || geneScores.get(k) < scoreThrTight)
			.sorted((k1, k2) -> geneScores.get(k1).compareTo(geneScores.get(k2)))
			.forEach(k -> System.out.println((OncoKB.get().isCancerGene(k) ? "O" :
				CancerGeneCensus.get().isCancerGene(k) ? "C" : CancerGeneBushman.get().isCancerGene(k) ? "B" : "") +
				"\t" + k + "\t" + geneToType.get(k)));

		System.out.println("\ntypes = " + types.size() + "\n");

		pairCounts.keySet().stream().filter(k -> pairCounts.get(k).size() > 1)
			.sorted((k1, k2) -> new Integer(pairCounts.get(k2).size()).compareTo(pairCounts.get(k1).size()))
			.forEach(k -> System.out.println(k + "\t" + pairCounts.get(k)));
	}


}
