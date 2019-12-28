package org.panda.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.*;
import org.panda.utility.graph.UndirectedGraph;
import org.panda.utility.statistics.FDR;

import java.io.BufferedWriter;
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
//		generateNeighborhoodSubgraphsForSignificantsRecursively("/home/ozgun/Analyses/Jaffe", 0.1);
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

//		writeSignifNeighForCorrBased("/Users/ozgun/Documents/Analyses/CPTAC-LUAD/correlation-tumor", StreamDirection.DOWNSTREAM);
		writeSignifNeighForCompBased("/Users/ozgun/Documents/Analyses/SUNIL-EXP12/netsig_Molm_Par_vs_AllEarly", StreamDirection.DOWNSTREAM, 0.1);
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

	public static Set<String> combine(Set<String>[] sets)
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
			writeSignifNeighForCompBased(dir, StreamDirection.DOWNSTREAM, netSig);
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
		printSignificantGenesRecursively(dir, dir, netSig);
	}

	public static void printSignificantGenesRecursively(String dir, String root, double netSig) throws IOException
	{
		String sigFile = dir + "/significance-pvals.txt";
		if (Files.exists(Paths.get(sigFile)))
		{
			Set<String>[] genes = getSignificantGenes(sigFile, netSig);

			Set<String> comb = combine(genes);

			List<String> sorted = comb.stream().sorted().map(g -> genes[0].contains(g) && genes[1].contains(g) ? g + "(a/i)" : genes[0].contains(g) ? g + "(a)" : g + "(i)").collect(Collectors.toList());

			if (!sorted.isEmpty())
			{
				System.out.println(dir.replace(root, "") + "\t" + CollectionUtil.merge(sorted, ", "));
			}
		}
		else
		{
			for (File subdir : new File(dir).listFiles())
			{
				if (subdir.isDirectory())
				{
					printSignificantGenesRecursively(subdir.getPath(), root, netSig);
				}
			}
		}
	}

	// Comparison based results

	public static void writeSignifNeighForCompBased(String dir, StreamDirection d, double fdrThr) throws IOException
	{
		if (!Files.exists(Paths.get(dir + "/results.txt"))) return;

		System.out.println("dir = " + dir);
		Set<String> goi = combine(getSignificantGenes(dir + "/significance-pvals.txt", fdrThr));
		System.out.println("signif = " + goi);
		System.out.println("signif size = " + goi.size());
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/causative-sig-neigh.sif", d);
		Set<String> ids = getIDsAtTheNeighborhood(dir + "/results.txt", goi, d);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(dir + "/causative.format", dir + "/causative-sig-neigh.format", goi, ids);
	}

	public static void writeGOINeighForCompBasedRecursively(String dir, Set<String> goi, StreamDirection d) throws IOException
	{
		String sifFile = dir + "/causative.sif";
		if (Files.exists(Paths.get(sifFile)))
		{
			writeGOINeighForCompBased(dir, goi, d);
		}
		else
		{
			for (File subdir : new File(dir).listFiles())
			{
				if (subdir.isDirectory())
				{
					writeGOINeighForCompBasedRecursively(subdir.getPath(), goi, d);
				}
			}
		}

	}

	public static void writeGOINeighForCompBased(String dir, Set<String> goi, StreamDirection d) throws IOException
	{
		writeGOINeighForCompBased(dir, goi, d, "causative-goi-neigh");
	}

	public static void writeGOINeighForCompBased(String dir, Set<String> goi, StreamDirection d, String outSIFNoExt) throws IOException
	{
		if (!Files.exists(Paths.get(dir + "/results.txt"))) return;

		System.out.println("dir = " + dir);
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/" + outSIFNoExt + ".sif", d);
		Set<String> ids = getIDsAtTheNeighborhood(dir + "/results.txt", goi, d);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(dir + "/causative.format", dir + "/" + outSIFNoExt + ".format", goi, ids);
	}

	// Correlation-based results

	/**
	 * @deprecated Provide the FDR threshold and use the next method.
	 */
	public static void writeSignifNeighForCorrBased(String dir, StreamDirection d) throws IOException
	{
		writeSignifNeighForCorrBased(dir, d, 0.1);
	}

	public static void writeSignifNeighForCorrBased(String dir, StreamDirection d, double fdrThr) throws IOException
	{
		System.out.println("dir = " + dir);
		Set<String> goi = getDownstreamEnrichedForCorrelation(dir + "/significance-pvals.txt", fdrThr);
		System.out.println("signif size = " + goi.size());
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/causative-sig-neigh.sif", d);
		Set<String> ids = getIDsAtTheNeighborhood(dir + "/results.txt", goi, d);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(dir + "/causative.format", dir + "/causative-sig-neigh.format", ids);
	}

	public static void writeGOINeighForCorrBased(String dir, Set<String> goi, StreamDirection d) throws IOException
	{
		System.out.println("dir = " + dir);
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/causative-goi-neigh.sif", d);
		Set<String> ids = getIDsAtTheNeighborhood(dir + "/results.txt", goi, d);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(dir + "/causative.format", dir + "/causative-goi-neigh.format", ids);
	}

	public static void writeSubsetFormat(String inFile, String outFile, Set<String> goi, Set<String> ids) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		Files.lines(Paths.get(inFile)).forEach(l ->
		{
			String[] t = l.split("\t");

			if (!t[2].equals("rppasite") || ids.contains(t[3].split("\\|")[0]))
			{
				FileUtil.writeln(l, writer);
			}
			else if (goi != null && goi.contains(t[1])) FileUtil.writeln(l, writer);
		});

//		if (goi != null) goi.forEach(g -> FileUtil.writeln("node\t" + g + "\tborderwidth\t2", writer));

		writer.close();
	}

	private static void writeSubsetFormat(String inFile, String outFile, Set<String> ids) throws IOException
	{
		writeSubsetFormat(inFile, outFile, null, ids);
	}

	private static Set<String> getDownstreamEnrichedForCorrelation(String sigFile, double fdrThr) throws IOException
	{
		Map<String, Double> pMap = Files.lines(Paths.get(sigFile)).skip(2).map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));
		return new HashSet<>(FDR.select(pMap, null, fdrThr));
	}

	public static Set<String> getIDsAtTheNeighborhood(String resultFile, Set<String> goi, StreamDirection d) throws IOException
	{
		String[] header = Files.lines(Paths.get(resultFile)).findFirst().get().split("\t");
		int sInd = ArrayUtil.indexOf(header, "Source");
		int tInd = ArrayUtil.indexOf(header, "Target");
		int sDI = ArrayUtil.indexOf(header, "Source data ID");
		int tDI = ArrayUtil.indexOf(header, "Target data ID");

		Set<String> ids = new HashSet<>();

		Files.lines(Paths.get(resultFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if ((d == StreamDirection.DOWNSTREAM && goi.contains(t[sInd])) ||
				(d == StreamDirection.UPSTREAM && goi.contains(t[tInd])) ||
				(d == StreamDirection.BOTHSTREAM && (goi.contains(t[tInd]) || goi.contains(t[sInd]))))
			{
				ids.add(t[sDI]);
				ids.add(t[tDI]);
			}
		});
		return ids;
	}
}
