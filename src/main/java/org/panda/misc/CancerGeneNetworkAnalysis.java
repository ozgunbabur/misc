package org.panda.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.GraphList;
import org.panda.utility.statistics.FDR;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Written to help Sikander's mom and potentially many other cancer patients who got their tumors sequenced.
 *
 * @deprecated see the new project cancer-network under PathwayAndDataAnalysis in GitHub
 * @author Ozgun Babur
 */
public class CancerGeneNetworkAnalysis
{
	Map<String, Set<Info>> altered;

	SIFEnum[] edgeTypes;

	public CancerGeneNetworkAnalysis()
	{
		altered = new HashMap<>();
	}

	public void addGeneAlteration(String gene, String letter)
	{
		addGeneAlteration(gene, letter, null, null, null);
	}

	public void addGeneAlteration(String gene, String letter, String tooltip, String bgColor, String borderColor)
	{
		if (!altered.containsKey(gene)) altered.put(gene, new HashSet<>());
		altered.get(gene).add(new Info(letter, tooltip, bgColor, borderColor));
	}

	public void setEdgeTypes(SIFEnum[] edgeTypes)
	{
		this.edgeTypes = edgeTypes;
	}

	public void generateNetwork(String sifNoExt) throws IOException
	{
		Map<String, Set<String>> cancerGenes = loadCancerGenes();
		System.out.println("cancerGenes.keySet() = " + cancerGenes.keySet().size());
		Map<String, Graph> graphs = loadGraphs();

		Set<String> genes = altered.keySet();
		genes.stream().filter(cancerGenes::containsKey).peek(System.out::print).forEach(g -> System.out.println("\t" + cancerGenes.get(g)));

		Set<String> enriched = findEnrichedNeighborhoods(graphs.get(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag()) != null ?
			graphs.get(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag()) : graphs.values().iterator().next(), cancerGenes.keySet());
		System.out.println("enriched = " + enriched);

		Set<String> genesInGraph = new HashSet<>();
		Set<String> edges = new HashSet<>();

		findEdges(genes, cancerGenes.keySet(), graphs, genesInGraph, edges);
//		findEdges(genesInGraph, genesInGraph, graphs, new HashSet<>(), edges);

		writeRelatedGeneJSON(sifNoExt, cancerGenes, graphs, genes);

		BufferedWriter sifWriter = Files.newBufferedWriter(Paths.get(sifNoExt + ".sif"));
		BufferedWriter fmtWriter = Files.newBufferedWriter(Paths.get(sifNoExt + ".format"));

		fmtWriter.write("node\tall-nodes\tcolor\t255 255 255");
		fmtWriter.write("\tnode\tall-nodes\tbordercolor\t0 0 0");

		for (String edge : edges)
		{
			sifWriter.write(edge + "\n");
		}

		genes.stream()./*filter(cancerGenes::containsKey).*/forEach(g ->
		{
			FileUtil.writeln(g, sifWriter);
			genesInGraph.add(g);
		});

		for (String gene : genesInGraph)
		{
			if (genes.contains(gene))
			{
				if (cancerGenes.containsKey(gene))
				{
					FileUtil.lnwrite("node\t" + gene + "\tbordercolor\t200 0 0", fmtWriter);
					FileUtil.lnwrite("node\t" + gene + "\tborderwidth\t2", fmtWriter);
				}

				for (Info info : altered.get(gene))
				{
					FileUtil.lnwrite("node\t" + gene + "\trppasite\t" + info, fmtWriter);
				}
			}
			if (cancerGenes.containsKey(gene))
			{
				Set<String> res = cancerGenes.get(gene);
				FileUtil.lnwrite("node\t" + gene + "\ttooltip\t" + res, fmtWriter);

				Color bgcolor = new Color(200, 255, 180);
//				Color bgcolor = res.contains("OncoKB-actionable") ? new Color(255, 120, 120) :
//					res.contains("OncoKB") ? new Color(150, 255, 150) :
//					res.contains("COSMIC") ? new Color(150, 150, 255) :
//					res.contains("Mutex") ? new Color(255, 150, 255) :
//					res.contains("Cooc") ? new Color(150, 255, 255) :
//						new Color(255, 255, 150);

				FileUtil.lnwrite("node\t" + gene + "\tcolor\t" + bgcolor.getRed() + " " + bgcolor.getGreen() + " " +
					bgcolor.getBlue(), fmtWriter);

//				if (genes.contains(gene))
//				{
//					FileUtil.lnwrite("node\t" + gene + "\thighlight\ton", fmtWriter);
//					FileUtil.lnwrite("node\t" + gene + "\thighlightcolor\t220 220 255", fmtWriter);
//				}
				if (enriched.contains(gene))
				{
					FileUtil.lnwrite("node\t" + gene + "\thighlight\ton", fmtWriter);
				}
			}
		}

		sifWriter.close();
		fmtWriter.close();
	}

	private void writeRelatedGeneList(String sifNoExt, Map<String, Set<String>> cancerGenes, Map<String, Graph> graphs, Set<String> genes) throws IOException
	{
		BufferedWriter geneWriter = Files.newBufferedWriter(Paths.get(sifNoExt + ".genes"));
		GraphList glist = new GraphList("merged");
		graphs.values().stream().forEach(glist::addGraph);
		for (String gene : genes)
		{
			Set<String> up = glist.getUpstream(gene);
			Set<String> dw = glist.getDownstream(gene);
			Set<String> ppi = glist.getNeighbors(gene);
			up.retainAll(cancerGenes.keySet());
			dw.retainAll(cancerGenes.keySet());
			ppi.retainAll(cancerGenes.keySet());
			ppi.removeAll(up);
			ppi.removeAll(dw);

			Set<String> neigh = new HashSet<>();
			neigh.addAll(up);
			neigh.addAll(dw);
			neigh.addAll(ppi);

			if (!(up.isEmpty() && dw.isEmpty() && ppi.isEmpty()))
			{
				geneWriter.write("\nGene: " + gene);

				for (Info info : altered.get(gene))
				{
					geneWriter.write("\nAlteration: " + info.getDocText());
				}

				if (!up.isEmpty()) geneWriter.write("\n\nUpstream: " + ArrayUtil.getString(", ", up.toArray()));
				if (!dw.isEmpty()) geneWriter.write("\nDownstream: " + ArrayUtil.getString(", ", dw.toArray()));
				if (!ppi.isEmpty()) geneWriter.write("\nComplex partner: " + ArrayUtil.getString(", ", ppi.toArray()));

				geneWriter.write("\n\nPathways: " + glist.getMediatorsInString(gene, neigh));

				geneWriter.write("\nComment: ");
				geneWriter.write("\n\nReference: ");
				geneWriter.write("\nQuote: ");
				geneWriter.write("\nComment: ");
				geneWriter.write("\n\nTCGA-evidence: ");
				geneWriter.write("\n\n#################################\n");
			}
		}

		geneWriter.close();
	}

	private void writeRelatedGeneJSON(String sifNoExt, Map<String, Set<String>> cancerGenes, Map<String, Graph> graphs, Set<String> genes) throws IOException
	{
		List<Map> json = new ArrayList<>();
		GraphList glist = new GraphList("merged");
		graphs.values().stream().forEach(glist::addGraph);
		for (String gene : genes)
		{
			Set<String> up = glist.getUpstream(gene);
			Set<String> dw = glist.getDownstream(gene);
			Set<String> ppi = glist.getNeighbors(gene);
			up.retainAll(cancerGenes.keySet());
			dw.retainAll(cancerGenes.keySet());
			ppi.retainAll(cancerGenes.keySet());
			ppi.removeAll(up);
			ppi.removeAll(dw);

			if (!(up.isEmpty() && dw.isEmpty() && ppi.isEmpty()))
			{
				Map genMap = new LinkedHashMap<>();
				json.add(genMap);

				genMap.put("gene-symbol", gene);

				List<String> altList = new ArrayList<>();
				for (Info info : altered.get(gene))
				{
					altList.add(info.getDocText());
				}
				genMap.put("alterations", altList);

				if (!up.isEmpty()) genMap.put("upstream", new ArrayList<>(up));
				if (!dw.isEmpty()) genMap.put("downstream", new ArrayList<>(dw));
				if (!ppi.isEmpty()) genMap.put("complex-partner", new ArrayList<>(ppi));
			}
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(sifNoExt + ".json"));
//		JsonUtils.writePrettyPrint(writer, json);
		writer.close();
	}

	private void findEdges(Set<String> set1, Set<String> set2, Map<String, Graph> graphs,
		Set<String> genesInGraph, Set<String> edges)
	{
		for (String gene1 : set1)
		{
			for (String type : graphs.keySet())
			{
				Graph graph = graphs.get(type);

				if (graph.isDirected())
				{
					for (String gene2 : ((DirectedGraph) graph).getDownstream(gene1))
					{
						if (set2.contains(gene2))
						{
							edges.add(gene1 + "\t" + type + "\t" + gene2 + "\t" + graph.getMediatorsInString(gene1, gene2));
							genesInGraph.add(gene1);
							genesInGraph.add(gene2);
						}
					}
					for (String gene2 : ((DirectedGraph) graph).getUpstream(gene1))
					{
						if (set2.contains(gene2))
						{
							edges.add(gene2 + "\t" + type + "\t" + gene1 + "\t" + graph.getMediatorsInString(gene2, gene1));
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

							edges.add(g1 + "\t" + type + "\t" + g2 + "\t" + graph.getMediatorsInString(g1, g2));
							genesInGraph.add(gene1);
							genesInGraph.add(gene2);
						}
					}
				}
			}
		}
	}

	private Set<String> findEnrichedNeighborhoods(Graph graph, Set<String> cancerGenes)
	{
		Map<String, Double>[] scores = graph.isDirected() ?
			((DirectedGraph) graph).getEnrichmentScores(altered.keySet(), null, DirectedGraph.NeighborType.DOWNSTREAM, 1, 5) :
			graph.getEnrichmentScores(altered.keySet(), null, 1, 5);
//		return scores[0].keySet().stream().filter(g -> scores[0].get(g) < 0.05).filter(cancerGenes::contains).collect(Collectors.toSet());
		for (int i = 0; i < scores.length; i++)
		{
			scores[i] = scores[i].keySet().stream().filter(cancerGenes::contains).collect(Collectors.toMap(g -> g, scores[i]::get));
		}
		return new HashSet<>(FDR.select(scores[0], scores[1], 0.1));
	}

	private Map<String, Set<String>> loadCancerGenes() throws IOException
	{
		Map<String, Set<String>> genes = new HashMap<>();

//		for (String res : CancerGeneBushman.get().getAllResources())
//		{
//			for (String gene : CancerGeneBushman.get().getSubset(res))
//			{
//				if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
//				genes.get(gene).add(res);
//			}
//		}

		for (String gene : CancerGeneCensus.get().getAllSymbols())
		{
			if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
			genes.get(gene).add("COSMIC");
		}

		for (String gene : OncoKB.get().getAllSymbols())
		{
			if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
		}

//		Files.lines(Paths.get("/home/babur/Documents/PanCan/pancan.txt")).skip(1).map(l -> l.split("\t"))
//			.filter(t -> Double.parseDouble(t[2]) < 0.1).map(t -> t[0]).forEach(gene ->
//		{
//			if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
//			genes.get(gene).add("Mutex");
//		});
//
//		Files.lines(Paths.get("/home/babur/Documents/PanCan/cooc-genes-right.txt")).skip(1).map(l -> l.split("\t"))
//			.filter(t -> Double.parseDouble(t[2]) < 0.1).map(t -> t[0]).forEach(gene ->
//		{
//			if (!genes.containsKey(gene)) genes.put(gene, new HashSet<>());
//			genes.get(gene).add("Cooc");
//		});

		return genes;
	}

	private Map<String, Graph> loadGraphs()
	{
		Map<String, Graph> map = new HashMap<>();

		for (SIFEnum type : edgeTypes)
		{
			map.put(type.getTag(), PathwayCommons.get().getGraph(type));
		}

		return map;
	}

	static class Info
	{
		public static final String DEFAULT_BACKGROUND_COLOR = "255 255 255";
		public static final String DEFAULT_BORDER_COLOR = "0 0 0";
		public static final String DEFAULT_TOOLTIP = "";

		public Info(String letter, String tooltip, String bgColor, String borderColor)
		{
			this.letter = letter;
			this.tooltip = tooltip;
			this.bgColor = bgColor;
			this.borderColor = borderColor;

			if (tooltip == null) this.tooltip = DEFAULT_TOOLTIP;
			if (bgColor == null) this.bgColor = DEFAULT_BACKGROUND_COLOR;
			if (borderColor == null) this.borderColor = DEFAULT_BORDER_COLOR;
		}

		String letter;
		String bgColor;
		String borderColor;
		String tooltip;

		@Override
		public String toString()
		{
			return ArrayUtil.getString("|", tooltip, letter, bgColor, borderColor);
		}

		public String getDocText()
		{
			if (letter.equals("m"))
			{
				if (tooltip.toLowerCase().endsWith("mutation"))
				{
					return tooltip;
				}
				else
				{
					return tooltip + " Mutation";
				}
			}
			else if (tooltip.contains("-")) return "Copy number loss";
			else return "Copy number gain";
		}
	}


	public void readCNVkit(String filename, double thrLog2, double saturationVal) throws IOException
	{
		ValToColor vtc = new ValToColor(new double[]{-saturationVal, 0, saturationVal},
			new Color[]{new Color(100, 100, 255), Color.WHITE, new Color(255, 100, 100)});

		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t"))
			.filter(t -> Math.abs(Double.valueOf(t[4])) >= thrLog2).filter(t -> !t[3].equals("-")).forEach(t ->
		{
			Double log2 = Double.valueOf(t[4]);
			String color = vtc.getColorInString(log2);

			for (String gene : t[3].split(","))
			{
				addGeneAlteration(gene, "c", "log2 = " + log2, color, null);
			}
		});
	}

	public void readReadGeneTrailsCNV(String filename, double saturationVal) throws IOException
	{
		ValToColor vtc = new ValToColor(new double[]{-saturationVal, 0, saturationVal},
			new Color[]{new Color(100, 100, 255), Color.WHITE, new Color(255, 100, 100)});

		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t"))
			.forEach(t ->
		{
			Double log2 = Math.log(Double.valueOf(t[7])) / Math.log(2);
			String color = vtc.getColorInString(log2);

			addGeneAlteration(t[4], "c", "log2 = " + log2, color, null);
		});
	}

	public void readMutect(String filename) throws IOException
	{
		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t")).
			forEach(t -> addGeneAlteration(t[5], "m", t[6], Info.DEFAULT_BACKGROUND_COLOR,
				t[6].startsWith("Splice") ? "255 0 0" : Info.DEFAULT_BORDER_COLOR));
	}

	public void readGeneTrailsMutations(String filename) throws IOException
	{
		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t")).
			forEach(t -> addGeneAlteration(t[8], "m", t[7], Info.DEFAULT_BACKGROUND_COLOR,
				t[16].toLowerCase().contains("splice") ? "255 0 0" : Info.DEFAULT_BORDER_COLOR));
	}

	public Set<String> readMutexGenes(String study)
	{
		Map<String, Double> mutexScores = MutexReader.readBestScoresRecursive("/home/babur/Documents/mutex/TCGA/" + study);
		return mutexScores.keySet().stream().filter(g -> mutexScores.get(g) < 0.05).collect(Collectors.toSet());
	}

	public Set<String> readPanCanMutexGenes() throws IOException
	{
		Map<String, Double> scoreMap = Files.lines(Paths.get("/home/babur/Documents/PanCan/pancan.txt")).skip(1)
			.map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[3])));

		return scoreMap.keySet().stream().filter(g -> scoreMap.get(g) < 0.1).collect(Collectors.toSet());
	}

	public void reportOverlapWithMutex(String study) throws IOException
	{
		Set<String> genes = altered.keySet();

		Set<String> specific = readMutexGenes(study);

		CollectionUtil.printNameMapping("altered", "mutex " + study);
		CollectionUtil.printVennSets(genes, specific);

		System.out.println("\n");

		Set<String> pancan = readPanCanMutexGenes();

		CollectionUtil.printNameMapping("altered", "pancan mutex");
		CollectionUtil.printVennSets(genes, pancan);
	}

	public static void main(String[] args) throws IOException
	{
		smmartP1Run();
	}

	private static void sikanderRun() throws IOException
	{
		CancerGeneNetworkAnalysis cgna = new CancerGeneNetworkAnalysis();
		Files.lines(Paths.get("/home/babur/Downloads/EUP_protein_variants.txt")).skip(1)
			.map(l -> l.split(",")).filter(t -> t[3].equals("germline")).map(t -> t[0]).filter(g -> !g.isEmpty())
			.forEach(g -> cgna.addGeneAlteration(g, "m"));

		cgna.generateNetwork("the-network");
	}

	private static void smmartP1Run() throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/SMMART/Patient1/";

		CancerGeneNetworkAnalysis cgna = new CancerGeneNetworkAnalysis();
		cgna.setEdgeTypes(new SIFEnum[]{SIFEnum.CONTROLS_STATE_CHANGE_OF, SIFEnum.IN_COMPLEX_WITH});

//		cgna.readMutect(dir + "cfDNA-mutect.csv");
//		cgna.readCNVkit(dir + "cfDNA-CNA.csv", 1, 2);

		cgna.readGeneTrailsMutations(dir + "irb16113_reported_mutations_by_mrn_Subject 101_deID-CompBio.csv");
		cgna.readReadGeneTrailsCNV(dir + "irb16113_reported_copy_number_by_mrn_Subject 101_deID-CompBio.csv", 2);

		cgna.generateNetwork(dir + "network-GeneTrails");
//		cgna.reportOverlapWithMutex("BRCA");
	}
}
