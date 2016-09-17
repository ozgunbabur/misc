package org.panda.misc;

import com.sun.org.apache.xpath.internal.SourceTree;
import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.CancerGeneBushman;
import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.FileUtil;
import org.panda.utility.graph.Graph;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Written to help Sikander's mom and potentially many other cancer patients who got their tumors sequenced.
 *
 * @author Ozgun Babur
 */
public class CancerGeneNetworkAnalysis
{
	public static void generateNetwork(Set<String> mutated, String sifNoExt) throws IOException
	{
		Map<String, Set<String>> cancerGenes = loadCancerGenes();
		System.out.println("cancerGenes.keySet() = " + cancerGenes.keySet().size());
		Map<String, Graph> graphs = loadGraphs();

		mutated.stream().filter(cancerGenes::containsKey).peek(System.out::print).forEach(g -> System.out.println("\t" + cancerGenes.get(g)));

		Set<String> genesInGraph = new HashSet<>();
		Set<String> edges = new HashSet<>();

		findEdges(mutated, cancerGenes.keySet(), graphs, genesInGraph, edges);
//		findEdges(genesInGraph, genesInGraph, graphs, new HashSet<>(), edges);

		BufferedWriter sifWriter = Files.newBufferedWriter(Paths.get(sifNoExt + ".sif"));
		BufferedWriter fmtWriter = Files.newBufferedWriter(Paths.get(sifNoExt + ".format"));

		fmtWriter.write("node\tall-nodes\tcolor\t255 255 255");
		fmtWriter.write("\tnode\tall-nodes\tbordercolor\t0 0 0");

		for (String edge : edges)
		{
			sifWriter.write(edge + "\n");
		}

		mutated.stream().filter(cancerGenes::containsKey).forEach(g ->
		{
			FileUtil.writeln(g, sifWriter);
			genesInGraph.add(g);
		});

		for (String gene : genesInGraph)
		{
			if (mutated.contains(gene))
			{
				FileUtil.lnwrite("node\t" + gene + "\tbordercolor\t200 0 0", fmtWriter);
				FileUtil.lnwrite("node\t" + gene + "\tborderwidth\t2", fmtWriter);
			}
			if (cancerGenes.containsKey(gene))
			{
				Set<String> res = cancerGenes.get(gene);
				FileUtil.lnwrite("node\t" + gene + "\ttooltip\t" + res, fmtWriter);

				Color bgcolor = res.contains("OncoKB-actionable") ? new Color(255, 120, 120) :
					res.contains("OncoKB") ? new Color(150, 255, 150) :
					res.contains("COSMIC") ? new Color(150, 150, 255) :
					res.contains("Mutex") ? new Color(255, 150, 255) :
					res.contains("Cooc") ? new Color(150, 255, 255) :
						new Color(255, 255, 150);

				FileUtil.lnwrite("node\t" + gene + "\tcolor\t" + bgcolor.getRed() + " " + bgcolor.getGreen() + " " +
					bgcolor.getBlue(), fmtWriter);
			}
		}

		sifWriter.close();
		fmtWriter.close();
	}

	private static void findEdges(Set<String> set1, Set<String> set2, Map<String, Graph> graphs,
		Set<String> genesInGraph, Set<String> edges)
	{
		for (String gene1 : set1)
		{
			for (String type : graphs.keySet())
			{
				Graph graph = graphs.get(type);

				if (graph.isDirected())
				{
					for (String gene2 : graph.getDownstream(gene1))
					{
						if (set2.contains(gene2))
						{
							edges.add(gene1 + "\t" + type + "\t" + gene2);
							genesInGraph.add(gene1);
							genesInGraph.add(gene2);
						}
					}
					for (String gene2 : graph.getUpstream(gene1))
					{
						if (set2.contains(gene2))
						{
							edges.add(gene2 + "\t" + type + "\t" + gene1);
							genesInGraph.add(gene1);
							genesInGraph.add(gene2);
						}
					}
				}
				else
				{
					for (String gene2 : graph.getNeighbors(gene1))
					{
						if (set2.contains(gene2))
						{
							String g1 = gene1;
							String g2= gene2;
							if (g2.compareTo(g1) < 0)
							{
								String temp = g1;
								g1 = g2;
								g2 = temp;
							}

							edges.add(g1 + "\t" + type + "\t" + g2);
							genesInGraph.add(gene1);
							genesInGraph.add(gene2);
						}
					}
				}
			}
		}
	}

	private static Map<String, Set<String>> loadCancerGenes() throws IOException
	{
		Map<String, Set<String>> genes = new HashMap<>();

		for (String res : CancerGeneBushman.get().getAllResources())
		{
			for (String gene : CancerGeneBushman.get().getSubset(res))
			{
				if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
				genes.get(gene).add(res);
			}
		}

		for (String gene : CancerGeneCensus.get().getAllSymbols())
		{
			if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
			genes.get(gene).add("COSMIC");
		}

		for (String gene : OncoKB.get().getAllSymbols())
		{
			if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
			String level = OncoKB.get().getLevel(gene);
			genes.get(gene).add(level != null && !level.isEmpty() ? "OncoKB-actionable" : "OncoKB");
		}

		Files.lines(Paths.get("/home/babur/Documents/PanCan/pancan.txt")).skip(1).map(l -> l.split("\t"))
			.filter(t -> Double.parseDouble(t[2]) < 0.1).map(t -> t[0]).forEach(gene ->
		{
			if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
			genes.get(gene).add("Mutex");
		});

		Files.lines(Paths.get("/home/babur/Documents/PanCan/cooc-genes-right.txt")).skip(1).map(l -> l.split("\t"))
			.filter(t -> Double.parseDouble(t[2]) < 0.1).map(t -> t[0]).forEach(gene ->
		{
			if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
			genes.get(gene).add("Cooc");
		});

		return genes;
	}

	private static Map<String, Graph> loadGraphs()
	{
		Map<String, Graph> map = new HashMap<>();
		map.put(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), PathwayCommons.get().getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF));
//		map.put(SIFEnum.CONTROLS_EXPRESSION_OF.getTag(), PathwayCommons.get().getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));
		map.put(SIFEnum.IN_COMPLEX_WITH.getTag(), PathwayCommons.get().getGraph(SIFEnum.IN_COMPLEX_WITH));
		return map;
	}

	public static void main(String[] args) throws IOException
	{
		Set<String> sikGenes = Files.lines(Paths.get("/home/babur/Downloads/EUP_protein_variants.txt")).skip(1)
			.map(l -> l.split(",")).filter(t -> t[3].equals("germline")).map(t -> t[0]).filter(g -> !g.isEmpty())
			.collect(Collectors.toSet());
		System.out.println("sikGenes.size() = " + sikGenes.size());

		generateNetwork(sikGenes, "the-network");
	}
}
