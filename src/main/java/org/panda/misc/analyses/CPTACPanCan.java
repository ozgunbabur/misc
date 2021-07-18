package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.causalpath.run.JasonizeResultGraphsRecursively;
import org.panda.misc.causalpath.CausalPathFormatHelper;
import org.panda.misc.causalpath.CausalPathShowPathsToCancer;
import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.resource.*;
import org.panda.utility.*;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.UndirectedGraph;
import org.panda.utility.statistics.FDR;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class CPTACPanCan
{
	public static final String DATA_BASE = "/Users/ozgun/Documents/Data/CPTAC-PanCan/";
	public static final String OUT_BASE = "/Users/ozgun/Documents/Analyses/CPTAC-PanCan/";
	public static final String MAP_FILE = DATA_BASE + "var_map_full.tsv";
	public static final String DIFFEXP_FILE = DATA_BASE + "full_diffexp_results_signatures.tsv";
	public static final String DENDRO_FILE = DATA_BASE + "full_de_cohort_cov_dendrogram.tsv";
	public static final String MISSING_SITES_FILE = DATA_BASE + "full_de_cohort_cov_with_missing_sites.tsv";
	public static final String MUTLIST_FILE = DATA_BASE + "full_de_cohort_cov_for_mutation_list.tsv";

//	public static final String PHOSPHO_FILE = DATA_BASE + "phospho_ardnmf_X.tsv";
//	public static final String ACETYL_FILE = DATA_BASE + "acetyl_ardnmf_X.tsv";
//	public static final String PROT_FILE = DATA_BASE + "prot_ardnmf_X.tsv";
//	public static final String RNA_FILE = DATA_BASE + "rna_ardnmf_X.tsv";
//	public static final String PROT_OUT = OUT_BASE + "prot-data-reg.txt";
//	public static final String RNA_OUT = OUT_BASE + "rna-data-reg.txt";
//	public static final String DELIM = ",";

	public static final String PHOSPHO_FILE = DATA_BASE + "phosphoproteome_X.tsv";
	public static final String ACETYL_FILE = DATA_BASE + "acetylome_X.tsv";
	public static final String PROT_FILE = DATA_BASE + "proteome_X.tsv";
	public static final String RNA_FILE = DATA_BASE + "tumor_rna_tpm_norm_reg_X.tsv";
	public static final String PROT_OUT = OUT_BASE + "prot-data.txt";
	public static final String RNA_OUT = OUT_BASE + "rna-data.txt";
	public static final String DELIM = "\t";

	static List<String> KEGG_NAMES = Arrays.asList("KEGG_FATTY_ACID_METABOLISM", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_FOCAL_ADHESION");


	public static void main(String[] args) throws IOException
	{
//		convertData();
		convertDiffExpData();

//		convertSingledOutDiffExpData(DATA_BASE + "lowtreg_157", "lowtreg");
//		convertSingledOutDiffExpData(DATA_BASE + "high_vs_low_treg", "high_vs_low_treg");
//		summarizeGraphs();

//		prepareAnalysisFoldersAgainstOthers();
		prepareAnalysisFoldersAgainstOthersDiffExp();
//		printRelationNumbers();

//		generateOncogenicProcessSubgraphs();
//		generateHallmarkSubgraphs();
//		generateHATHDACSubgraphs();
//		generateAcetylFocusedSubgraphs();
//		generateNetworkSignificanceSubgraphs();
//		generateMetabolicViews();
//		generateGOIFocusCorrelationBased();
//		generateRBUpstream();


//		printOncoKBCoverage();
//		printAllSampleNames();
//		printOncogenicFrequencies();
//		printOverlapOfPCWithKEGGGeneSets();

//		temo();
	}


	private static void temo() throws IOException
	{
		String base = OUT_BASE + "clusters/metabolic-sifs-diffexp";

		Map<String, Map<String, Integer>> cntMap = new HashMap<>();

		for (File dir : new File(base).listFiles())
		{
			if (dir.isDirectory())
			{
				for (File file : dir.listFiles())
				{
					if (file.getName().endsWith(".format"))
					{
						int cnt = (int) FileUtil.countLines(file.getPath());

						if (!cntMap.containsKey(file.getName())) cntMap.put(file.getName(), new HashMap<>());
						cntMap.get(file.getName()).put(dir.getName(), cnt);
					}
				}
			}
		}

		for (String pName : cntMap.keySet())
		{
			System.out.println(pName);
			Map<String, Integer> map = cntMap.get(pName);
			map.keySet().stream().sorted(Comparator.comparing(map::get).reversed()).forEach(c -> System.out.println(c + "\t" + map.get(c)));
		}
	}

	private static void printOncogenicFrequencies()
	{
		String base = OUT_BASE + "clusters/against-others-mapping-based";
		Map<String, Set<String>> geneToClusters = new HashMap<>();
		Map<String, Set<String>> clusterToGenes = new HashMap<>();
		Map<String, Set<String>> regulatorToTargets = new HashMap<>();

		for (File dir : new File(base).listFiles())
		{
			if (dir.isDirectory())
			{
				String cluster = dir.getName();

				String filename = dir.getPath() + "/sitespec/oncogenic-changes.sif";

				if (Files.exists(Paths.get(filename)))
				{
					Set<String> genes = FileUtil.lines(filename).map(l -> l.split("\t")).filter(t -> t.length > 2)
//						.map(t -> Arrays.asList(t[0], t[2])).flatMap(Collection::stream).collect(Collectors.toSet());
						.peek(t ->
						{
							if (!regulatorToTargets.containsKey(t[0])) regulatorToTargets.put(t[0], new HashSet<>());
							regulatorToTargets.get(t[0]).add(t[2]);
						})
						.map(t -> t[0]).collect(Collectors.toSet());


					clusterToGenes.put(cluster, genes);

					genes.forEach(gene ->
					{
						if (!geneToClusters.containsKey(gene)) geneToClusters.put(gene, new HashSet<>());
						geneToClusters.get(gene).add(cluster);
					});
				}
			}
		}

		geneToClusters.keySet().stream().sorted((o1, o2) ->
		{
			int c = Integer.compare(geneToClusters.get(o2).size(), geneToClusters.get(o1).size());
			if (c == 0)
			{
				c = o1.compareTo(o2);
			}
			return c;
		}).forEach(gene ->
		{
			System.out.println(gene + "\t" + geneToClusters.get(gene).stream().sorted().collect(Collectors.toList()) +
				"\t" + regulatorToTargets.getOrDefault(gene, Collections.emptySet()).stream().sorted().collect(Collectors.toList()));
		});

		// Generate cluster similarity graph

		Map<String, Integer> edgeCnts = new HashMap<>();
		for (Set<String> clusters : geneToClusters.values())
		{
			for (String c1 : clusters)
			{
				for (String c2 : clusters)
				{
					if (c1.compareTo(c2) < 0)
					{
						String key = c1 + "\t" + c2;
						edgeCnts.put(key, edgeCnts.getOrDefault(key, 0) + 1);
					}
				}
			}
		}

		System.out.println();
		Integer max = edgeCnts.values().stream().max(Integer::compareTo).get();
		ValToColor vtc = new ValToColor(new double[]{0, max}, new Color[]{Color.WHITE, Color.BLACK});

		BufferedWriter sifWriter = FileUtil.newBufferedWriter(OUT_BASE + "cluster-similarity.sif");
		edgeCnts.keySet().stream().filter(k -> edgeCnts.get(k) > 3).forEach(k -> FileUtil.writeln(k.replaceFirst("\t", "\tinteracts-with\t"), sifWriter));
		FileUtil.closeWriter(sifWriter);
		BufferedWriter fmtWriter = FileUtil.newBufferedWriter(OUT_BASE + "cluster-similarity.format");
		edgeCnts.keySet().forEach(k -> FileUtil.writeln("edge\t" + k.replaceFirst("\t", " interacts-with ") + "\tcolor\t" + vtc.getColorInString(edgeCnts.get(k)), fmtWriter));
		FileUtil.closeWriter(fmtWriter);

		Set<Set<String>> components = new HashSet<>();
		edgeCnts.keySet().stream().sorted((o1, o2) -> Integer.compare(edgeCnts.get(o2), edgeCnts.get(o1))).forEach(k ->
		{
			String c1 = k.substring(0, k.indexOf("\t"));
			String c2 = k.substring(k.indexOf("\t") + 1);
			Set<Set<String>> overlapping = components.stream().filter(s -> s.contains(c1) || s.contains(c2)).collect(Collectors.toSet());

			if (overlapping.isEmpty())
			{
				Set<String> set = new HashSet<>(Arrays.asList(c1, c2));
				components.add(set);
			}
			else if (overlapping.size() == 1)
			{
				Set<String> set = overlapping.iterator().next();
				set.add(c1);
				set.add(c2);
			}
			else
			{
				// print components before merge
				for (Set<String> set : overlapping)
				{
					System.out.println(set);
				}

				components.removeAll(overlapping);
				Set<String> set = new HashSet<>();
				overlapping.forEach(set::addAll);
				set.add(c1);
				set.add(c2);
				components.add(set);
			}
		});
		System.out.println();
	}

	private static void printAllSampleNames() throws IOException
	{
		String[] header = Files.lines(Paths.get(PROT_OUT)).findFirst().get().split("\t");
		for (int i = 5; i < header.length; i++)
		{
			System.out.println("value-column = " + header[i]);
		}
	}
	private static void printOncoKBCoverage() throws IOException
	{
		String dir = OUT_BASE + "/clusters/against-others-diffexp";

		Set<String> found = new HashSet<>();

		for (File cDir : new File(dir).listFiles())
		{
			String cluster = cDir.getName();
			if (cluster.startsWith(".")) continue;

			for (File tDir : cDir.listFiles())
			{
				String tName = tDir.getName();

				if (tName.equals("phospho")) found.addAll(readAllGenes(tDir.getPath() + "/oncogenic-changes-phospho.sif"));
				else if (tName.equals("totprot") || tName.equals("rnaseq")) found.addAll(readAllGenes(tDir.getPath() + "/oncogenic-changes-express.sif"));
			}
		}

		System.out.println("found.size() = " + found.size());

		Set<String> oncogenes = new HashSet<>(OncoKB.get().getOncogenes());
		Set<String> tumSupps = new HashSet<>(OncoKB.get().getTumorSuppressors());
		Set<String> both = CollectionUtil.getIntersection(oncogenes, tumSupps);
		System.out.println("both.size() = " + both.size());
		CollectionUtil.removeIntersection(oncogenes, tumSupps);

		System.out.println("in OncoKB = " + oncogenes.size());

		CollectionUtil.printVennCounts(oncogenes, found);
	}

	private static void generateGOIFocusCorrelationBased() throws IOException
	{
		CausalPathSubnetwork.writeGOINeighForCorrBased(OUT_BASE + "correlation-based-reg/sitespec", Collections.singleton("RB1"), StreamDirection.BOTHSTREAM, "RB1-neighborhood", false);
		CausalPathSubnetwork.writeGOINeighForCorrBased(OUT_BASE + "correlation-based-reg/sitespec", new HashSet<>(Arrays.asList("GRK1", "GRK2", "GRK3", "GRK4", "GRK6", "GRK7")), StreamDirection.BOTHSTREAM, "GRK-neighborhood", false);
	}

	private static void generateOncogenicProcessSubgraphs() throws IOException
	{
		CausalPathShowPathsToCancer.generateRecursivelyComparison(OUT_BASE + "clusters/against-others-diffexp", null);
		CausalPathShowPathsToCancer.generateRecursivelyComparison(OUT_BASE + "clusters/against-others-mapping-based", null);
		CausalPathShowPathsToCancer.generateRecursivelyCorrelation(OUT_BASE + "correlation-based");
		CausalPathShowPathsToCancer.generateRecursivelyCorrelation(OUT_BASE + "correlation-based-reg");
	}



	private static void generateHallmarkSubgraphs() throws IOException
	{
		generateHallmarkSubgraphsClusters(OUT_BASE + "clusters/against-others-diffexp/");
		generateHallmarkSubgraphsClusters(OUT_BASE + "clusters/against-others-mapping-based/");
		generateHallmarkSubgraphsCorrelation(OUT_BASE + "correlation-based/");
		generateHallmarkSubgraphsCorrelation(OUT_BASE + "correlation-based-reg/");
	}

	private static void generateHallmarkSubgraphsClusters(String base) throws IOException
	{
		String enrichmentFilename = DATA_BASE + "W_fgsea.tsv";

		String[] header = ("no\t" + Files.lines(Paths.get(enrichmentFilename)).findFirst().get()).split("\t");
		int pInd = ArrayUtil.indexOf(header, "\"padj\"");
		int nameInd = ArrayUtil.indexOf(header, "\"pathway\"");
		int cInd = ArrayUtil.indexOf(header, "\"id\"");

		Files.lines(Paths.get(enrichmentFilename)).skip(1).map(l -> l.split("\t"))
			.filter(t -> t.length > cInd)
			.filter(t -> Double.valueOf(t[pInd]) < 0.1 && t[nameInd].startsWith("\"HALLMARK_"))
			.forEach(t ->
			{
				String hallmark = t[nameInd].replaceAll("\"", "");
				String cluster = t[cInd];
				String[] dirs = new String[]{base + cluster + "/all", base + cluster + "/sitespec", base + cluster + "/rnaseq", base + cluster + "/totprot"};
				try
				{
					for (String dir : dirs)
					{
						CausalPathSubnetwork.writeGOINeighForCompBased(dir, MSigDB.get().getGeneSet(hallmark), StreamDirection.BOTHSTREAM, hallmark);
					}
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			});

		for (File cDir : new File(base).listFiles())
		{
			for (String name : KEGG_NAMES)
			{
				String[] dirs = new String[]{cDir.getPath() + "/all", cDir.getPath() + "/sitespec", cDir.getPath() + "/rnaseq", cDir.getPath() + "/totprot"};
				for (String dir : dirs)
				{
					CausalPathSubnetwork.writeGOINeighForCompBased(dir, MSigDB.get().getGeneSet(name), StreamDirection.BOTHSTREAM, name);
				}
			}
		}
	}

	private static void generateRBUpstream() throws IOException
	{
		generateRBUpstream(OUT_BASE + "clusters/against-others-diffexp");
	}

	private static void generateRBUpstream(String base) throws IOException
	{
		for (File clus : new File(base).listFiles())
		{
			String dir = clus.getPath() + "/sitespec";
			CausalPathSubnetwork.writeGOINeighForCompBased(dir, Collections.singleton("RB1"), StreamDirection.UPSTREAM, "RB1-upstream");
		}
	}

	private static void generateHallmarkSubgraphsCorrelation(String base)
	{
		Map<String, Set<String>> geneSets = MSigDB.get().getSetsNameFiltered(name -> name.startsWith("HALLMARK_"));
		KEGG_NAMES.forEach(n -> geneSets.put(n, MSigDB.get().getGeneSet(n)));

		geneSets.forEach((name, geneSet) ->
		{
			try
			{
				String dir = base + "sitespec";
				CausalPathSubnetwork.writeGOINeighForCorrBased(dir, geneSet, StreamDirection.BOTHSTREAM, name);
				dir = base + "rnaseq";
				CausalPathSubnetwork.writeGOINeighForCorrBased(dir, geneSet, StreamDirection.BOTHSTREAM, name);
				dir = base + "totprot";
				CausalPathSubnetwork.writeGOINeighForCorrBased(dir, geneSet, StreamDirection.BOTHSTREAM, name);
			}
			catch (Exception e)
			{
				e.printStackTrace();
			}
		});
	}

	private static void generateHATHDACSubgraphs() throws IOException
	{
		generateHATHDACSubgraphsClusters(OUT_BASE + "/clusters/against-others-diffexp");
		generateHATHDACSubgraphsClusters(OUT_BASE + "/clusters/against-others-mapping-based");
		generateHATHDACSubgraphsCorrelation(OUT_BASE + "correlation-based/");
		generateHATHDACSubgraphsCorrelation(OUT_BASE + "correlation-based-reg/");
	}

	private static void generateHATHDACSubgraphsClusters(String dir) throws IOException
	{
		Set<String> goi = GO.get().getGenesContainingKeywordInTermNames(
			new HashSet<>(Arrays.asList("histone acetyltransferase activity", "histone deacetylase activity")),
			new HashSet<>(Arrays.asList("regulation of")));

		System.out.println("goi.size() = " + goi.size());
		StreamDirection direction = StreamDirection.BOTHSTREAM;
		String graphName = "GENESET_HAT_HDAC";

		for (File cDir : new File(dir).listFiles())
		{
			String cluster = cDir.getName();
			if (cluster.startsWith(".")) continue;

			String d = cDir.getPath() + "/phospho";
			CausalPathSubnetwork.writeGOINeighForCompBased(d, goi, direction, graphName);
			d = cDir.getPath() + "/rnaseq";
			CausalPathSubnetwork.writeGOINeighForCompBased(d, goi, direction, graphName);
			d = cDir.getPath() + "/totprot";
			CausalPathSubnetwork.writeGOINeighForCompBased(d, goi, direction, graphName);
		}
	}
	private static void generateHATHDACSubgraphsCorrelation(String base2) throws IOException
	{
		Set<String> goi = GO.get().getGenesContainingKeywordInTermNames(
			new HashSet<>(Arrays.asList("histone acetyltransferase activity", "histone deacetylase activity")),
			new HashSet<>(Arrays.asList("regulation of")));

		System.out.println("goi.size() = " + goi.size());
		StreamDirection direction = StreamDirection.BOTHSTREAM;
		String graphName = "GENESET_HAT_HDAC";

		String dir = base2 + "sitespec";
		CausalPathSubnetwork.writeGOINeighForCorrBased(dir, goi, direction, graphName);
		dir = base2 + "rnaseq";
		CausalPathSubnetwork.writeGOINeighForCorrBased(dir, goi, direction, graphName);
		dir = base2 + "totprot";
		CausalPathSubnetwork.writeGOINeighForCorrBased(dir, goi, direction, graphName);
	}

	private static void generateAcetylFocusedSubgraphs() throws IOException
	{
		generateAcetylFocusedSubgraphsClusters(OUT_BASE + "clusters/against-others-diffexp");
		generateAcetylFocusedSubgraphsClusters(OUT_BASE + "clusters/against-others-mapping-based");
		generateAcetylFocusedSubgraphsCorrelation(OUT_BASE + "correlation-based");
		generateAcetylFocusedSubgraphsCorrelation(OUT_BASE + "correlation-based-reg");
	}

	private static void generateAcetylFocusedSubgraphsClusters(String base) throws IOException
	{
		String sifName = "acetyl-focused";

		FileUtil.processDirsRecursive(new File(base), dir ->
		{
			String filename = dir.getPath() + "/results.txt";
			if (new File(filename).exists())
			{
				BufferedWriter writer = FileUtil.newBufferedWriter(dir.getPath() + "/" + sifName + ".sif");
				Set<String> ids = new HashSet<>();

				FileUtil.lines(filename).skip(1).map(l -> l.split("\t")).filter(t -> t.length > 7 &&
					(isIDAcetylPeptide(t[4]) || isIDAcetylPeptide(t[7]))).forEach(t ->
				{
					ids.add(t[4]);
					ids.add(t[7]);

					FileUtil.writeln(ArrayUtil.merge("\t", t[0], t[1], t[2], "", t[3]), writer);
				});
				writer.close();

				CausalPathSubnetwork.writeSubsetFormat(dir.getPath() + "/causative.format", dir.getPath() + "/" + sifName +".format", null, ids);
			}
		});
	}

	private static boolean isIDAcetylPeptide(String id)
	{
		String[] t = id.split("-");
		for (int i = 2; i < t.length; i++)
		{
			if (t[i].equals("A")) return true;
		}
		return false;
	}

	private static void generateAcetylFocusedSubgraphsCorrelation(String base) throws IOException
	{
		String sifName = "acetyl-focused";

		FileUtil.processDirsRecursive(new File(base), dir ->
		{
			String filename = dir.getPath() + "/results.txt";
			if (new File(filename).exists())
			{
				BufferedWriter writer = FileUtil.newBufferedWriter(dir.getPath() + "/" + sifName + ".sif");
				Set<String> ids = new HashSet<>();

				FileUtil.lines(filename).skip(1).map(l -> l.split("\t")).filter(t -> t.length > 5 &&
					(isIDAcetylPeptide(t[4]) || isIDAcetylPeptide(t[5]))).forEach(t ->
				{
					ids.add(t[4]);
					ids.add(t[5]);

					FileUtil.writeln(ArrayUtil.merge("\t", t[0], t[1], t[2], "", t[3]), writer);
				});
				writer.close();

				CausalPathSubnetwork.writeSubsetFormat(dir.getPath() + "/causative.format", dir.getPath() + "/" + sifName +".format", null, ids);
			}
		});
	}

	private static void generateNetworkSignificanceSubgraphs() throws IOException
	{
		double netSig = 0.1;
		CausalPathSubnetwork.generateNeighborhoodSubgraphsForSignificantsRecursively(OUT_BASE + "clusters/against-others-diffexp", netSig);
		CausalPathSubnetwork.generateNeighborhoodSubgraphsForSignificantsRecursively(OUT_BASE + "clusters/against-others-mapping-based", netSig);
		CausalPathSubnetwork.writeSignifNeighForCorrBasedRecursive(OUT_BASE + "correlation-based", StreamDirection.DOWNSTREAM, netSig);
		CausalPathSubnetwork.writeSignifNeighForCorrBasedRecursive(OUT_BASE + "correlation-based-reg", StreamDirection.DOWNSTREAM, netSig);
	}

	private static void generateMetabolicViews() throws IOException
	{
		generateMetabolicViews(OUT_BASE + "prot-data-diffexp.txt", OUT_BASE + "rna-data-diffexp.txt", OUT_BASE + "clusters/metabolic-sifs-diffexp");
//		generateMetabolicViews(OUT_BASE + "prot-data.txt", OUT_BASE + "rna-data.txt", OUT_BASE + "graphs/clusters/all-inclusive/metabolic");
	}

	private static void generateMetabolicViews(String protFile, String rnaFile, String outBase) throws IOException
	{
		Map<String, Set<String>> geneSets = MSigDB.get().getSetsNameFiltered(name -> KEGG_NAMES.contains(name));

//		geneSets.forEach((name, genes) ->
//		{
//			Set<String> inHGNC = genes.stream().filter(g -> HGNC.get().getSymbol(g) != null).collect(Collectors.toSet());
//			System.out.println(name + "\t" + inHGNC.size() + "/" + genes.size());
//		});

		DirectedGraph catGraph = loadCatalysisPrecedesGraph();
		DirectedGraph stcGraph = loadStateChangeGraph();

		String[] header = FileUtil.readHeader(protFile);

		for (int i = 5; i < header.length; i++)
		{
			String cluster = header[i];

			String dir = outBase + "/" + cluster;
			Files.createDirectories(Paths.get(dir));

			for (String name : KEGG_NAMES)
			{
				Set<String> genes = geneSets.get(name);

				String outFmt = dir + "/" + name + ".format";
				BufferedWriter fmtWriter = FileUtil.newBufferedWriter(outFmt);
				List<String> fmtLines = CausalPathFormatHelper.getFormatLinesAdjSignedP(genes, protFile, rnaFile, "ID", "Symbols", "Sites", "Modification", "Effect", cluster, 0.1);
				fmtLines.addAll(CausalPathFormatHelper.getFormatLinesForOncogenicEffect(genes));
				fmtLines.forEach(l -> FileUtil.writeln(l, fmtWriter));
				fmtWriter.close();
				Set<String> changedNodes = fmtLines.stream().map(l -> l.split("\t")).filter(t -> t[2].equals("color") || t[2].equals("rppasite")).map(t -> t[1]).collect(Collectors.toSet());

				DirectedGraph graph = name.endsWith("ADHESION") ? stcGraph : catGraph;
//				Set<String> selectNodes = name.endsWith("ADHESION") ? changedNodes : genes;
				Set<String> selectNodes = genes;
//				DirectedGraph subgraph = name.endsWith("ADHESION") ? graph.getNeighborhoodInTheInducedSubgraph(genes, changedNodes) : graph.getInducedSubgraphWithoutDisconnectedNodes(selectNodes);
				DirectedGraph subgraph = graph.getNeighborhoodInTheInducedSubgraph(genes, changedNodes);

				String outSIF = dir + "/" + name + ".sif";
				BufferedWriter sifWriter = FileUtil.newBufferedWriter(outSIF);
				subgraph.write(sifWriter);
				sifWriter.close();
			}
		}
	}

	private static void printOverlapOfPCWithKEGGGeneSets()
	{
		for (String keggName : KEGG_NAMES)
		{
			System.out.println("keggName = " + keggName);
			Set<String> geneSet = MSigDB.get().getGeneSet(keggName);
			PCPathwayHGNC pcp = new PCPathwayHGNC();
			Map<String, Double> pvals = pcp.calculateEnrichment(geneSet, 5, 200);
			Map<String, Double> qvals = FDR.getQVals(pvals, null);
			pvals.keySet().stream().sorted(Comparator.comparing(pvals::get)).filter(k -> pvals.get(k) < 0.05).forEach(k ->
				System.out.println(k + "\t" + pcp.getPathwayName(k) + "\t" + pvals.get(k) + "\t" + qvals.get(k) + "\t" +
					pcp.getOverlappingGenes(k, geneSet)));

		}
	}

	private static DirectedGraph loadCatalysisPrecedesGraph()
	{
		String allFile = "/Users/ozgun/Documents/Data/PathwayCommonsV12/PathwayCommons12.All.hgnc.sif";
		String keggFile = "/Users/ozgun/Documents/Data/PathwayCommonsV12/KEGG-catalysis-precedes.sif";
		String humancycFile = "/Users/ozgun/Documents/Data/PathwayCommonsV12/HumanCyc-catalysis-precedes.sif";
		String reactomeFile = "/Users/ozgun/Documents/Data/PathwayCommonsV12/Reactome-catalysis-precedes.sif";
		DirectedGraph graph = new DirectedGraph("PC", SIFEnum.CATALYSIS_PRECEDES.getTag());
		graph.load(keggFile, Collections.singleton(SIFEnum.CATALYSIS_PRECEDES.getTag()));
		graph.load(humancycFile, Collections.singleton(SIFEnum.CATALYSIS_PRECEDES.getTag()));
		graph.load(reactomeFile, Collections.singleton(SIFEnum.CATALYSIS_PRECEDES.getTag()));
//		graph.load(allFile, Collections.singleton(SIFEnum.CATALYSIS_PRECEDES.getTag()));
		return graph;
	}

	private static UndirectedGraph loadInteractionGraph()
	{
		String file = "/Users/ozgun/Documents/Data/PathwayCommonsV12/PathwayCommons12.All.hgnc.sif";
		UndirectedGraph graph = new UndirectedGraph("PC", SIFEnum.INTERACTS_WITH.getTag());
		graph.load(file, new HashSet<>(Arrays.asList(SIFEnum.INTERACTS_WITH.getTag(), SIFEnum.IN_COMPLEX_WITH.getTag())));
		return graph;
	}

	private static DirectedGraph loadStateChangeGraph()
	{
		String file = "/Users/ozgun/Documents/Data/PathwayCommonsV12/PathwayCommons12.All.hgnc.sif";
		DirectedGraph graph = new DirectedGraph("PC", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		graph.load(file, new HashSet<>(Arrays.asList(SIFEnum.INTERACTS_WITH.getTag(), SIFEnum.IN_COMPLEX_WITH.getTag())));
		return graph;
	}

	private static void printRelationNumbers() throws IOException
	{
		String dir = OUT_BASE + "/clusters/against-others";

		System.out.println("Cluster\tExplaining phospho\tExplaining global prot\tExplaining mRNA");
		for (File cDir : new File(dir).listFiles())
		{
			String cluster = cDir.getName();
			if (cluster.startsWith(".")) continue;

			int pCnt = 0;
			int tCnt = 0;
			int rCnt = 0;

			for (File tDir : cDir.listFiles())
			{
				String tName = tDir.getName();

				if (tName.equals("phospho")) pCnt = readResultCount(tDir);
				else if (tName.equals("totprot")) tCnt = readResultCount(tDir);
				else if (tName.equals("rnaseq")) rCnt = readResultCount(tDir);
			}

			System.out.println(cluster + "\t" + pCnt + "\t" + tCnt + "\t" + rCnt);
		}
	}
//	private static int readResultCount(File dir) throws IOException
//	{
//		String filename = dir.getPath() + "/causative.sif";
//		if (Files.exists(Paths.get(filename)))
//			return (int) Files.lines(Paths.get(filename)).map(l -> l.split("\t")).filter(t -> t.length > 2).map(t -> t[0] + "\t" + t[1] + "\t" + t[2]).distinct().count();
//		return 0;
//	}
	private static int readResultCount(File dir) throws IOException
	{
		String filename = dir.getName().equals("phospho") ? dir.getPath() + "/oncogenic-changes-phospho.sif" : dir.getPath() + "/oncogenic-changes-express.sif";
		if (Files.exists(Paths.get(filename)))
			return (int) Files.lines(Paths.get(filename)).map(l -> l.split("\t")).filter(t -> t.length > 2).map(t -> t[0] + "\t" + t[1] + "\t" + t[2]).distinct().count();
		return 0;
	}
	private static Set<String> readAllGenes(String file) throws IOException
	{
		if (Files.exists(Paths.get(file)))
			return Files.lines(Paths.get(file)).map(l -> l.split("\t"))
				.map(t -> t.length > 2 ? Arrays.asList(t[0], t[2]) : Arrays.asList(t[0]))
				.flatMap(Collection::stream).collect(Collectors.toSet());
		return Collections.emptySet();
	}

	private static void convertDiffExpData() throws IOException
	{
//		String base = OUT_BASE + "clusters/signatures/";
//		String inFile = DIFFEXP_FILE;
//		String outFile = OUT_BASE + "data-signatures.tsv";

//		String base = OUT_BASE + "clusters/dendrograms/";
//		String inFile = DENDRO_FILE;
//		String outFile = OUT_BASE + "data-dendrogram.tsv";

		String base = OUT_BASE + "clusters/missing-sites/";
		String inFile = MISSING_SITES_FILE;
		String outFile = OUT_BASE + "data-missing-sites.tsv";

		new File(base).mkdirs();

		Map<String, String> idToGene = readIDToGeneMapping();

		String[] header = FileUtil.lines(inFile).findFirst().get().split("\t");
		int idInd = ArrayUtil.indexOf(header, "index");
		int modInd = ArrayUtil.indexOf(header, "feature");
		int pInd = ArrayUtil.indexOf(header, "adj.P.Val");
		int fcInd = ArrayUtil.indexOf(header, "logFC");
		int clusterInd = ArrayUtil.indexOf(header, "id");

		Map<String, Map<String, Double>> idToVals = new HashMap<>();
		Map<String, String[]> idToProps = new HashMap<>();

		Set<String> mem = new HashSet<>();

		FileUtil.lines(inFile).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[idInd];

			String[] props = idToProps.get(id);

			if (props == null)
			{
				String sym = idToGene.get(id);
				if (HGNC.get().getSymbol(sym) == null) return;
				String mod;
				switch (t[modInd])
				{
					case "phosphoproteome":
						mod = "P";
						break;
					case "acetylome":
						mod = "A";
						break;
					case "transcriptome":
						mod = "R";
						break;
					case "proteome":
						mod = "G";
						break;
					default:
						mod = "";
				}

				if (mod.equals("P") || mod.equals("A"))
				{
					List<String> siteList = getSitesFromID(id);

					String site = CollectionUtil.merge(siteList, "|");

					String cpID = sym + "-" + mod + "-" + site.replaceAll("\\|", "-");

					while (mem.contains(cpID))
					{
						String ss = cpID.substring(cpID.lastIndexOf("-") + 1);
						cpID = Character.isDigit(ss.charAt(0)) ? cpID.substring(0, cpID.lastIndexOf("-") + 1) + (Integer.valueOf(ss) + 1) : cpID + "-2";
					}
					mem.add(cpID);

					props = new String[]{cpID, sym, site, mod, ""};
				}
				else
				{
					props = new String[]{sym + "-" + mod, sym, "", mod, ""};
				}
				idToProps.put(id, props);
			}

			double p = Double.valueOf(t[pInd]);
			if (t[fcInd].startsWith("-")) p = -p;

			if (!idToVals.containsKey(id)) idToVals.put(id, new HashMap<>());
			idToVals.get(id).put(t[clusterInd], p);
		});

		List<String> clusters = idToVals.keySet().stream().map(idToVals::get).map(Map::keySet).flatMap(Collection::stream).distinct().collect(Collectors.toList());
//		clusters.sort(String::compareTo);
		clusters.sort(Comparator.comparing(Double::valueOf));

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tFeature\tEffect");
		clusters.forEach(c -> FileUtil.tab_write(c, writer));
		idToProps.forEach((id, props) ->
		{
			FileUtil.lnwrite(ArrayUtil.merge("\t", props), writer);
			Map<String, Double> valMap = idToVals.get(id);
			clusters.forEach(c -> FileUtil.tab_write(valMap.getOrDefault(c, Double.NaN), writer));
		});
		writer.close();

		// Convert transcriptome

//		idToVals.clear();
//		FileUtil.lines(inFile).skip(1).map(l -> l.split("\t")).forEach(t ->
//		{
//			if (!t[modInd].equals("transcriptome")) return;
//
//			String id = t[0];
//			double p = Double.valueOf(t[pInd]);
//			if (t[fcInd].startsWith("-")) p = -p;
//
//			if (!idToVals.containsKey(id)) idToVals.put(id, new HashMap<>());
//			idToVals.get(id).put(t[clusterInd], p);
//		});
//
//		BufferedWriter writerT = FileUtil.newBufferedWriter(OUT_BASE + "rna-data-diffexp.txt");
//		writerT.write("Gene");
//		clusters.forEach(c -> FileUtil.tab_write(c, writerT));
//		idToVals.forEach((id, valMap) ->
//		{
//			String sym = idToGene.get(id);
//			if (HGNC.get().getSymbol(sym) == null) return;
//			FileUtil.lnwrite(sym, writerT);
//			clusters.forEach(c -> FileUtil.tab_write(valMap.get(c), writerT));
//		});
//		writerT.close();
	}

	private static void convertSingledOutDiffExpData(String srcDir, String outNamePrefix) throws IOException
	{
		Map<String, String> idToGene = readIDToGeneMapping();
		String phosphoFile = srcDir + "/phosphoproteome_de.tsv";
		String acetylFile = srcDir + "/acetylome_de.tsv";
		String globalProtFile = srcDir + "/proteome_de.tsv";
		String rnaFile = srcDir + "/transcriptome_de.tsv";
		int[] pInd = new int[2];

		BufferedWriter writer1 = FileUtil.newBufferedWriter(OUT_BASE + outNamePrefix + "-prot-data.tsv");
		writer1.write("ID\tSymbols\tSites\tModification\tEffect\tSignedAdjP");

		Set<String> idMem = new HashSet<>();

		// phosphoproteomics

		String[] header = FileUtil.lines(phosphoFile).map(l -> l.replaceAll("\"", "")).findFirst().get().split("\t");
		pInd[0] = ArrayUtil.indexOf(header, "adj.P.Val") + 1;
		pInd[1] = ArrayUtil.indexOf(header, "logFC") + 1;
		FileUtil.lines(phosphoFile).skip(1).map(l -> l.replaceAll("\"", "")).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = idToGene.get(t[0]);
			List<String> siteList = getSitesFromID(t[0]);
			String id = sym + "-" + CollectionUtil.merge(siteList, "-") + "-P";
			id = getUniqueID(id, idMem);
			double val = Double.valueOf(t[pInd[0]]);
			if (t[pInd[1]].startsWith("-")) val = -val;

			FileUtil.lnwrite(id + "\t" + sym + "\t" + CollectionUtil.merge(siteList, "|") + "\t" + "P" + "\t\t" + val, writer1);
		});

		// acetylomics

		header = FileUtil.lines(acetylFile).map(l -> l.replaceAll("\"", "")).findFirst().get().split("\t");
		pInd[0] = ArrayUtil.indexOf(header, "adj.P.Val") + 1;
		pInd[1] = ArrayUtil.indexOf(header, "logFC") + 1;
		FileUtil.lines(acetylFile).skip(1).map(l -> l.replaceAll("\"", "")).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = idToGene.get(t[0]);
			List<String> siteList = getSitesFromID(t[0]);
			String id = sym + "-" + CollectionUtil.merge(siteList, "-") + "-A";
			id = getUniqueID(id, idMem);
			double val = Double.valueOf(t[pInd[0]]);
			if (t[pInd[1]].startsWith("-")) val = -val;

			FileUtil.lnwrite(id + "\t" + sym + "\t" + CollectionUtil.merge(siteList, "|") + "\t" + "A" + "\t\t" + val, writer1);
		});

		// global protein

		header = FileUtil.lines(globalProtFile).map(l -> l.replaceAll("\"", "")).findFirst().get().split("\t");
		pInd[0] = ArrayUtil.indexOf(header, "adj.P.Val") + 1;
		pInd[1] = ArrayUtil.indexOf(header, "logFC") + 1;
		FileUtil.lines(globalProtFile).skip(1).map(l -> l.replaceAll("\"", "")).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = idToGene.get(t[0]);
			double val = Double.valueOf(t[pInd[0]]);
			if (t[pInd[1]].startsWith("-")) val = -val;

			FileUtil.lnwrite(sym + "\t" + sym + "\t\t\t\t" + val, writer1);
		});

		writer1.close();

		// RNA

		BufferedWriter writer2 = FileUtil.newBufferedWriter(OUT_BASE + outNamePrefix + "-rna-data.tsv");
		writer2.write("Gene\tSignedAdjP");

		header = FileUtil.lines(rnaFile).map(l -> l.replaceAll("\"", "")).findFirst().get().split("\t");
		pInd[0] = ArrayUtil.indexOf(header, "adj.P.Val") + 1;
		pInd[1] = ArrayUtil.indexOf(header, "logFC") + 1;
		FileUtil.lines(rnaFile).skip(1).map(l -> l.replaceAll("\"", "")).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = idToGene.get(t[0]);
			double val = Double.valueOf(t[pInd[0]]);
			if (t[pInd[1]].startsWith("-")) val = -val;

			FileUtil.lnwrite(sym + "\t" + val, writer2);
		});
	}

	private static String getUniqueID(String proposed, Set<String> memory)
	{
		if (memory.contains(proposed))
		{
			int i = 1;
			String newID = proposed + "-" + (i++);
			while (memory.contains(newID))
			{
				newID = proposed + "-" + (i++);
			}
			memory.add(newID);
			return newID;
		}
		else
		{
			memory.add(proposed);
			return proposed;
		}
	}

	private static void prepareAnalysisFoldersAgainstOthers() throws IOException
	{
		String base = OUT_BASE + "clusters/against-others-mapping-based/";
		new File(base).mkdirs();

		Map<String, String> clusMap = readClusterMappings();
		List<String> clusters = clusMap.values().stream().distinct().sorted().collect(Collectors.toList());

		for (int i = 0; i < clusters.size(); i++)
		{
			String c1 = clusters.get(i);
			String dir = base + c1 + "/";
			new File(dir).mkdirs();

			String valPointers = getParamsValPointersAgainstOthers(clusMap, c1);

			String innerDir = dir + "sitespec/";
			new File(innerDir).mkdirs();
			BufferedWriter writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(SITESPEC_PARAMS + valPointers);
			writer.close();
			innerDir = dir + "totprot/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(TOTPROT_PARAMS + valPointers);
			writer.close();
			innerDir = dir + "rnaseq/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(RNA_PARAMS + valPointers);
			writer.close();
			innerDir = dir + "all/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(ALL_PARAMS + valPointers);
			writer.close();
		}
	}

	private static void prepareAnalysisFoldersAgainstOthersDiffExp() throws IOException
	{
//		String base = OUT_BASE + "clusters/signatures/";
//		String base = OUT_BASE + "clusters/dendrogram/";
		String base = OUT_BASE + "clusters/missing-sites/";
		new File(base).mkdirs();

//		List<String> clusters = Arrays.asList(FileUtil.lines(OUT_BASE + "data-signatures.tsv").findFirst().get().split("\t"));
//		List<String> clusters = Arrays.asList(FileUtil.lines(OUT_BASE + "data-dendrogram.tsv").findFirst().get().split("\t"));
		List<String> clusters = Arrays.asList(FileUtil.lines(OUT_BASE + "data-missing-sites.tsv").findFirst().get().split("\t"));
		clusters = clusters.subList(5, clusters.size());

		for (int i = 0; i < clusters.size(); i++)
		{
			String c1 = clusters.get(i);
			String dir = base + c1 + "/";
			new File(dir).mkdirs();

			String valPointers = "value-column = " + c1;

			String innerDir = dir + "sitespec/";
			new File(innerDir).mkdirs();
			BufferedWriter writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(correctForDE(SITESPEC_PARAMS) + valPointers);
			writer.close();
			innerDir = dir + "globprot/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(correctForDE(TOTPROT_PARAMS) + valPointers);
			writer.close();
			innerDir = dir + "rnaseq/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(correctForDE(RNA_PARAMS) + valPointers);
			writer.close();
			innerDir = dir + "all/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(correctForDE(ALL_PARAMS) + valPointers);
			writer.close();
		}
	}

	private static void prepareAnalysisFoldersPairwise() throws IOException
	{
		String base = OUT_BASE + "clusters/pairwise/";
		new File(base).mkdirs();

		Map<String, String> samToClus = readClusterMappings();
		List<String> clusters = samToClus.values().stream().distinct().sorted().collect(Collectors.toList());

		for (int i = 0; i < clusters.size() - 1; i++)
		{
			String c1 = clusters.get(i);
			for (int j = i + 1; j < clusters.size(); j++)
			{
				String c2 = clusters.get(j);

				String dir = base + c1 + "-vs-" + c2 + "/";
				new File(dir).mkdirs();

				String valPointers = getParamsValPointersPairwise(samToClus, c1, c2);

				String innerDir = dir + "phospho/";
				new File(innerDir).mkdirs();
				BufferedWriter writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
				writer.write(SITESPEC_PARAMS + valPointers);
				writer.close();
				innerDir = dir + "totprot/";
				new File(innerDir).mkdirs();
				writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
				writer.write(TOTPROT_PARAMS + valPointers);
				writer.close();
				innerDir = dir + "rnaseq/";
				new File(innerDir).mkdirs();
				writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
				writer.write(RNA_PARAMS + valPointers);
				writer.close();
				innerDir = dir + "all/";
				new File(innerDir).mkdirs();
				writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
				writer.write(ALL_PARAMS + valPointers);
				writer.close();
			}
		}
	}

	private static String getParamsValPointersPairwise(Map<String, String> samToClus, String c1, String c2)
	{
		StringBuilder sb = new StringBuilder();

		samToClus.forEach((s, c) ->
		{
			if (c.equals(c1)) sb.append("\ntest-value-column = ").append(s);
			else if (c.equals(c2)) sb.append("\ncontrol-value-column = ").append(s);
		});
		return sb.toString();
	}

	private static String getParamsValPointersAgainstOthers(Map<String, String> clusMap, String c1)
	{
		StringBuilder sb = new StringBuilder();

		clusMap.forEach((s, c) ->
		{
			if (c.equals(c1)) sb.append("\ntest-value-column = ").append(s);
			else sb.append("\ncontrol-value-column = ").append(s);
		});
		return sb.toString();
	}

	private static Map<String, String> readIDToGeneMapping() throws IOException
	{
		return Files.lines(Paths.get(MAP_FILE)).skip(1)
			.map(l -> l.split("\t"))
//			.peek(t -> t[0] = t[0].startsWith("ENSG") ? t[0].substring(0, t[0].indexOf(".")) : t[0])
			.collect(Collectors.toMap(t -> t[0], t -> t[2]));//, (o, o2) -> o));
	}

	private static Map<String, String> readClusterMappings() throws IOException
	{
		return Files.lines(Paths.get(DATA_BASE + "mappings.tsv")).skip(1)
			.map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> t[2]));
	}

	private static Map<String, String> readSampleToStudyMap() throws IOException
	{
		return Files.lines(Paths.get(DATA_BASE + "mappings/grade_info.tsv")).skip(1)
			.map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> t[2]));
	}

	private static Set<String> readCompleteSampleIDs() throws IOException
	{
		String header = Files.lines(Paths.get(DATA_BASE + "imputed_v_res_f/pan_rna_X.tsv")).findFirst().get();
		header = header.substring(header.indexOf("\t") + 1);
		return new HashSet<>(Arrays.asList(header.split("\t")));
	}


	private static void convertData() throws IOException
	{
		Map<String, String> idToGene = readIDToGeneMapping();

		System.out.print("Reading phosphoproteome...");
		Map<String[], Map<String, Double>> pp = processSiteSpecificData(PHOSPHO_FILE, "P", idToGene);
		System.out.print("Done\nReading acetylome...");

 		List<String> samples = pp.values().iterator().next().keySet().stream().sorted().collect(Collectors.toList());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(PROT_OUT));

		writer.write("ID\tSymbols\tSites\tModification\tEffect");
		samples.forEach(s -> FileUtil.tab_write(s, writer));

		pp.forEach((s, map) ->
		{
			FileUtil.lnwrite(ArrayUtil.merge("\t", s), writer);
			samples.forEach(c -> FileUtil.tab_write(map.get(c), writer));
		});

		pp = null;
		System.gc();

		Map<String[], Map<String, Double>> ap = processSiteSpecificData(ACETYL_FILE, "A", idToGene);
		System.out.print("Done\nReading proteome...");

		ap.forEach((s, map) ->
		{
			FileUtil.lnwrite(ArrayUtil.merge("\t", s), writer);
			samples.forEach(c -> FileUtil.tab_write(map.getOrDefault(c, Double.NaN), writer));
		});

		ap = null;
		System.gc();

		Map<String, Map<String, Double>> tp = processProteomicData(idToGene);
		System.out.print("Done\nReading RNA...");
		tp.forEach((s, map) ->
		{
			FileUtil.lnwrite(s + "\t" + s + "\t\t\t", writer);
			samples.forEach(c -> FileUtil.tab_write(map.get(c), writer));
		});

		writer.close();

		tp = null;
		System.gc();

		convertRNAData(idToGene);
		System.out.println("Done");
	}

	private static void convertRNAData(Map<String, String> idToSym) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(RNA_OUT));
		String header = Files.lines(Paths.get(RNA_FILE)).findFirst().get().replaceAll(DELIM, "\t");

		writer.write("Gene" + header);

		Files.lines(Paths.get(RNA_FILE)).skip(1).forEach(l ->
		{
			String id = l.substring(0, l.indexOf(DELIM));
			l = l.substring(l.indexOf(DELIM) + 1);
			String sym = idToSym.get(id);
			if (sym == null) throw new RuntimeException("Missing ID --> Symbol mapping!: ID = " + id);
			if (HGNC.get().getSymbol(sym) == null) return;
			FileUtil.lnwrite(sym + "\t" + l.replaceAll(DELIM, "\t"), writer);
		});

		writer.close();
	}

	private static Map<String, Map<String, Double>> processProteomicData(Map<String, String> idToSym) throws IOException
	{
		Map<String, String[]> rowMap = new HashMap<>();

		Files.lines(Paths.get(PROT_FILE)).skip(1).map(l -> l.split(DELIM)).forEach(t ->
		{
			if (countNonEmpty(t) < 51) return;

			String id = t[0];
			String sym = idToSym.get(id);
			if (HGNC.get().getSymbol(sym) == null) return;

			if (sym == null) throw new RuntimeException("Missing ID --> Symbol mapping!: ID = " + id);

			rowMap.put(sym, t);
		});

		Map<String, Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(PROT_FILE)).findFirst().get().split(DELIM);

		rowMap.forEach((sym, t) ->
		{
			HashMap<String, Double> valMap = new HashMap<>();
			for (int i = 1; i < header.length; i++)
			{
				valMap.put(header[i], t.length <= i || t[i].isEmpty() || t[i].equals("NA") ? Double.NaN : Double.valueOf(t[i]));
			}
			map.put(sym, valMap);
		});

		return map;
	}

	private static Map<String[], Map<String, Double>> processSiteSpecificData(String inFile, String mod,
		Map<String, String> idToSym) throws IOException
	{
		Map<String[], Map<String, Double>> map = new HashMap<>();

		Set<String> mem = new HashSet<>();

		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split(DELIM);

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split(DELIM)).forEach(t ->
		{
			if (countNonEmpty(t) < 51) return;

			String id = t[0];
			String sym = idToSym.get(id);
			if (HGNC.get().getSymbol(sym) == null) return;

			if (sym == null) throw new RuntimeException("Missing ID --> Symbol mapping!: ID = " + id);

			List<String> sites = getSitesFromID(id);

			HashMap<String, Double> valMap = new HashMap<>();
			for (int i = 1; i < header.length; i++)
			{
				valMap.put(header[i], t.length <= i || t[i].isEmpty() || t[i].equals("NA") ? Double.NaN : Double.valueOf(t[i]));
			}

			String site = CollectionUtil.merge(sites, "|");

			id = sym + "-" + site.replaceAll("\\|", "-") + "-" + mod;

			while (mem.contains(id))
			{
				String s = id.substring(id.lastIndexOf("-")+1);
				id = Character.isDigit(s.charAt(0)) ? id.substring(0, id.lastIndexOf("-") + 1) + (Integer.valueOf(s) + 1) : id + "-2";
			}
			mem.add(id);

			map.put(new String[]{id, sym, site, mod, ""}, valMap);
		});

		return map;
	}

	private static int countNonEmpty(String[] t)
	{
		return (int) Arrays.stream(t).filter(s -> !s.isEmpty()).count();
	}

	private static List<String> getSitesFromID(String id)
	{
		id = id.split("_")[2];
		id = id.substring(0, id.length() - 1);
		return Arrays.asList(id.split("s|t|y|k"));
	}

	private static String correctForDE(String s)
	{
		return s.replaceFirst("value-transformation = significant-change-of-mean", "value-transformation = signed-p-values")
			.replaceAll("fdr-threshold-for-data", "threshold-for-data")
			.replaceFirst("minimum-sample-size = 3\n", "")
//			.replaceFirst("calculate-network-significance = true", "calculate-network-significance = false")
			;
	}

	private static void summarizeGraphs() throws IOException
	{
		String dir = OUT_BASE + "clusters/specific/high_vs_low_treg/all";
		CausalPathSubnetwork.writeSiteSpecAndStrongExpCompBased(dir, "summary");
		JasonizeResultGraphsRecursively.generate(dir, dir, Collections.singleton("summary"), dir, "causative.json");
//		CausalPathSubnetwork.writeSignifNeighForCompBased(dir, StreamDirection.BOTHSTREAM, 0.1);
//		JasonizeResultGraphsRecursively.generate(dir, dir, Collections.singleton("causative-sig-neigh"), dir, "causative.json");
	}

//	public static final String PARAM_START = "proteomics-values-file = ../../../../data-signatures.tsv\n" +
//	public static final String PARAM_START = "proteomics-values-file = ../../../../data-dendrogram.tsv\n" +
	public static final String PARAM_START = "proteomics-values-file = ../../../../data-missing-sites.tsv\n" +
		"id-column = ID\n" +
		"symbols-column = Symbols\n" +
		"sites-column = Sites\n" +
		"feature-column = Feature\n" +
		"effect-column = Effect\n" +
		"\n" +
		"value-transformation = significant-change-of-mean\n" +
		"fdr-threshold-for-data-significance = 0.1 protein\n" +
		"fdr-threshold-for-data-significance = 0.1 phosphoprotein\n" +
		"fdr-threshold-for-data-significance = 0.1 acetylprotein\n" +
		"fdr-threshold-for-data-significance = 0.1 rna\n" +
		"\n" +
		"minimum-sample-size = 3\n" +
		"color-saturation-value = 15\n" +
		"\n" +
		"calculate-network-significance = true\n" +
		"permutations-for-significance = 10000\n" +
		"fdr-threshold-for-network-significance = 0.1\n" +
		"use-network-significance-for-causal-reasoning = true\n" +
		"\n";

	public static final String SITESPEC_PARAMS = PARAM_START +
		"relation-filter-type = site-specific-only\n" +
		"\n";

	public static final String TOTPROT_PARAMS = PARAM_START +
		"relation-filter-type = expression-only\n" +
		"data-type-for-expressional-targets = protein\n" +
		"\n";

	public static final String ALL_PARAMS = PARAM_START +
		"data-type-for-expressional-targets = rna\n" +
		"\n";

	public static final String RNA_PARAMS = ALL_PARAMS +
		"relation-filter-type = expression-only\n" +
		"\n";
}
