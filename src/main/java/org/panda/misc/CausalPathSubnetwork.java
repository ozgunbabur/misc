package org.panda.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;
import org.panda.utility.graph.UndirectedGraph;
import org.panda.utility.statistics.FDR;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Imagine you have a comparison-based CausalPath result with a lot of relations. It is not easily readable. So we
 * visualize a subset of those relations. Since some of those subset graphs are so small, we can fit in further PPI
 * relations to the visualizations with the same proteomic data overlay on the interacting proteins. The result is a
 * mixture of CausalPath graph and PPI networks, with all nodes showing significant proteomic data on them.
 */
public class CausalPathSubnetwork
{
	public static void main(String[] args) throws IOException
	{
		generateNeighborhoodSubgraphsForSignificantsRecursively("/home/ozgun/Analyses/Jaffe", 0.1);
//		if (true) System.exit(0);

//		String dir = "/home/ozgun/Analyses/Aslan-platelet/with-feedback-more-hyp-added/";
//		String dir = "/home/ozgun/Analyses/Hisham/Proteome/primed-vs-naive-phospho/";
//		String dir = "/home/ozgun/Analyses/Hisham/Proteome/mouse-primed-vs-naive-expression-rnaseq/";

//		generateSubsetWithPPI(dir + "causative.sif", Arrays.asList("IKBKB", "BIN2"), loadPPIGraph());

//		SIFFileUtil.writeNeighborhood(dir + "causative.sif", getSignificantGenes(
//			dir + "significance-pvals.txt", 0.1), dir + "causative-neighborhood-of-enriched-fdr0.1.sif");

//		SIFFileUtil.writeIntersection("/home/ozgun/Analyses/Hisham/Proteome/primed-vs-naive-expression-rnaseq/causative.sif",
//			"/home/ozgun/Analyses/Hisham/Proteome/primed-vs-naive-expression/causative.sif",
//			"/home/ozgun/Analyses/Hisham/Proteome/primed-vs-naive-expression-rnaseq/intersection.sif");
	}

	public static UndirectedGraph loadPPIGraph()
	{
		UndirectedGraph graph = (UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.INTERACTS_WITH);
		graph.merge((UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.IN_COMPLEX_WITH));
		return graph;
	}

	public static void generateSubsetWithPPI(String sifFile, Collection<String> seed, UndirectedGraph ppiGraph)
		throws IOException
	{
		StringBuilder sb = new StringBuilder(sifFile.substring(0, sifFile.lastIndexOf("."))).append("-subgraph");
		seed.stream().sorted().forEach(gene -> sb.append("-").append(gene));
		sb.append(".sif");

		SIFFileUtil.generateSubsetAndAddPPI(sifFile, seed, sb.toString(), ppiGraph);
	}



	public static Set<String>[] getSignificantGenes(String sigFile, double fdrThr) throws IOException
	{
		Map<String, Double> actMap = new HashMap<>();
		Map<String, Double> inhMap = new HashMap<>();

		Files.lines(Paths.get(sigFile)).skip(2).map(l -> l.split("\t")).forEach(t ->
		{
			actMap.put(t[0], Double.valueOf(t[2]));
			inhMap.put(t[0], Double.valueOf(t[3]));
		});

		return new Set[]{new HashSet<>(FDR.select(actMap, null, fdrThr)), new HashSet<>(FDR.select(inhMap, null, fdrThr))};
	}

	private static Set<String> combine(Set<String>[] sets)
	{
		Set<String> comb = new HashSet<>(sets[0]);
		comb.addAll(sets[1]);
		return comb;
	}

	public static void generateNeighborhoodSubgraphsForSignificantsRecursively(String dir, double netSig) throws IOException
	{
		String sifFile = dir + "/causative.sif";
		String sigFile = dir + "/significance-pvals.txt";
		if (Files.exists(Paths.get(sifFile)) && Files.exists(Paths.get(sigFile)))
		{
			String subFile = sifFile.substring(0, sifFile.lastIndexOf(".")) + "-subset.sif";
			Set<String> genes = combine(getSignificantGenes(sigFile, netSig));

			if (!genes.isEmpty())
			{
				SIFFileUtil.writeNeighborhood(sifFile, genes, subFile);
				FileUtil.copyFile(dir + "/causative.format", subFile.substring(0, subFile.lastIndexOf(".")) + ".format");
				System.out.println(dir + ": " + genes);
			}
		}
		else
		{
			for (File subdir : new File(dir).listFiles())
			{
				if (subdir.isDirectory())
				{
					generateNeighborhoodSubgraphsForSignificantsRecursively(subdir.getPath(), netSig);
				}
			}
		}
	}

	public static void printSignificantGenesRecursively(String dir, double netSig) throws IOException
	{
		String sigFile = dir + "/significance-pvals.txt";
		if (Files.exists(Paths.get(sigFile)))
		{
			Set<String>[] genes = getSignificantGenes(sigFile, netSig);

			Set<String> comb = combine(genes);

			List<String> sorted = comb.stream().sorted().map(g -> genes[0].contains(g) && genes[1].contains(g) ? g + "(a/i)" : genes[0].contains(g) ? g + "(a)" : g + "(i)").collect(Collectors.toList());

			if (!sorted.isEmpty())
			{
				System.out.println(dir.substring(dir.lastIndexOf("/") + 1) + "\t" + CollectionUtil.merge(sorted, ", "));
			}
		}
		else
		{
			for (File subdir : new File(dir).listFiles())
			{
				if (subdir.isDirectory())
				{
					printSignificantGenesRecursively(subdir.getPath(), netSig);
				}
			}
		}
	}
}
