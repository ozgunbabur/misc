package org.panda.misc.analyses;

import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.misc.MutexReader;
import org.panda.resource.GO;
import org.panda.resource.HGNC;
import org.panda.resource.proteomics.RPPAIDMapper;
import org.panda.utility.*;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.GeneSetEnrichment;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class AnilsPancan
{
	public static final String DIR = "/media/ozgun/6TB/TCGA-pancan/";
	public static final String MIX = "Mix";

	public static void main(String[] args) throws IOException
	{
//		generateSampelToTissueMapping();
//		generateComputeDirs();
//		generateSIFFormatFileForTP53PancanMutex();

//		annotateGeneLisFor3GoTerms();

//		printRasFamilyEnrichment();

//		convertDataCentricToGeneCentric();

		generateApoptosisSubgraphs();
	}

	private static void generateSampelToTissueMapping() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "sample-to-tissue-mapping.txt"));
		writer.write("Sample\tGroup");
		Files.lines(Paths.get(DIR + "pancan_samples.txt")).skip(1).map(l -> l.split("\t"))
			.forEach(t ->
			{
				String sample = t[1];
				String group = t[2];
				if (!t[3].equals("Not_Applicable")) group += "-" + t[3];

				FileUtil.lnwrite(sample + "\t" + group, writer);
			});
		writer.close();
	}

	private static void generateComputeDirs() throws IOException
	{
		String dir = DIR + "whole/types";

		for (File sub : new File(dir).listFiles())
		{
			Files.createDirectories(Paths.get(sub.getPath() + "/no-network"));
			BufferedWriter writer = Files.newBufferedWriter(Paths.get(sub.getPath() + "/no-network/parameters.txt"));

			writer.write(
				"search-on-signaling-network = false\n" +
				"first-level-random-iteration = 10000\n" +
				"max-group-size = 5\n" +
				"score-cutoff = 0.05\n" +
				"data-file = ../DataMatrix.txt\n" +
				"minimum-alteration-count-threshold = 2\n" +
				"gene-limit = 500");

			writer.close();
		}
	}

	private static void generateSIFFormatFileForTP53PancanMutex() throws IOException
	{
		String mutexRoot = DIR + "whole/types";
		Set<MutexReader.Group> groups = MutexReader.readMutexResultsRecursive(
			mutexRoot, new HashSet<>());

		Map<String, Color> colorMap = Files.lines(Paths.get(DIR + "TCGAcolors.csv"))
			.map(l -> l.split(",")).collect(Collectors.toMap(t -> t[1], t -> Color.decode(t[2])));

		Map<String, String> geneToStudy = new HashMap<>();

		groups.forEach(g -> g.shortenLoc(mutexRoot + "/"));
		groups.forEach(g -> g.fromDir = g.fromDir.split("/")[0].split("-")[0]);

		int i = 0;
		for (MutexReader.Group group : groups)
		{
			if (!group.genes.contains("TP53")) continue;
			if (group.score > 0.01) continue;
			i++;

			Memory.printIfNew(group.fromDir);

			for (String gene : group.genes)
			{
				if (!geneToStudy.containsKey(gene)) geneToStudy.put(gene, group.fromDir);
				else if (!geneToStudy.get(gene).equals(group.fromDir)) geneToStudy.put(gene, MIX);
			}
		}

		System.out.println("total groups = " + i);
		System.exit(0);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "network.format"));
		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tbordercolor\t0 0 0");

		geneToStudy.forEach((g, s) -> FileUtil.lnwrite("node\t" + g + "\tcolor\t" + (s.equals(MIX) ? "255 255 255" :
			colorMap.get(s).getRed() + " " + colorMap.get(s).getGreen() + " " + colorMap.get(s).getBlue()), writer));

		geneToStudy.forEach((g, s) -> {
			if (!s.equals(MIX)) FileUtil.lnwrite("node\t" + g + "\ttooltip\t" + s, writer);
		});

		writer.close();

		cropSIFFile(DIR + "TP53-neighborhood.sif", DIR + "TP53-neigh-cropped.sif", geneToStudy.keySet());
	}

	private static void cropSIFFile(String inFile, String outFile, Set<String> genes) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		Files.lines(Paths.get(inFile)).forEach(l ->
		{
			String[] t = l.split("\t");
			if (genes.contains(t[0]) && genes.contains(t[2]))
			{
				FileUtil.writeln(l, writer);
			}
		});

		genes.forEach(g -> FileUtil.lnwrite(g, writer));

		writer.close();
	}

	private static void annotateGeneLisFor3GoTerms() throws IOException
	{
		List<String> genes = Files.lines(Paths.get("/home/ozgun/Analyses/MCL1/GOI.txt")).collect(Collectors.toList());

		Set<String> motilityGenes = new HashSet<>(GO.get().getGenesOfTerm("GO:2000147"));
		motilityGenes.addAll(GO.get().getGenesOfTerm("GO:0016477"));
		Set<String> adhesionGenes = GO.get().getGenesOfTerm("GO:0030155");
		Set<String> lipidGenes = GO.get().getGenesOfTerm("GO:0006631");

		genes.forEach(gene ->
		{
			System.out.print("\n" + gene + "\t");
			if (lipidGenes.contains(gene)) System.out.print("1"); else System.out.print("0");
			System.out.print("\t");
			if (adhesionGenes.contains(gene)) System.out.print("1"); else System.out.print("0");
			System.out.print("\t");
			if (motilityGenes.contains(gene)) System.out.print("1"); else System.out.print("0");
		});


		Map<String, Double> pvals = GO.get().calculateEnrichment(genes, null, 200, Integer.MAX_VALUE);

		List<String> select = FDR.select(pvals, null, 0.1);

		select.forEach(go ->
		{
			String name = GO.get().getNameOfTerm(go);

			Set<String> hit = new HashSet<>(genes);
			hit.retainAll(GO.get().getGenes(go));

			System.out.println(pvals.get(go) + "\t" + go + " " + name + "\t" + hit);
		});
	}

	private static void printRasFamilyEnrichment() throws IOException
	{
		Map<String, Set<String>> geneSets = Files.lines(Paths.get("/home/ozgun/Analyses/MCL1/Ras-family.txt"))
			.map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> new HashSet<>(Arrays.asList(t[1].split("; ")))));

		geneSets.put("All", geneSets.values().stream().flatMap(Collection::stream).collect(Collectors.toSet()));

		Set<String> selected = Files.lines(Paths.get("/home/ozgun/Analyses/MCL1/GOI.txt")).collect(Collectors.toSet());

//		selected.retainAll(geneSets.get("Rho"));
//		System.out.println("selected.size() = " + selected.size());

		Map<String, Double> pvals = GeneSetEnrichment.calculateEnrichment(selected, HGNC.get().getAllSymbols(), 1, 1000, geneSets);

		pvals.keySet().stream().sorted(Comparator.comparing(pvals::get)).forEach(k -> System.out.println(k + "\t" + pvals.get(k)));
	}

	private static void convertDataCentricToGeneCentric() throws IOException
	{
		String dcSIFFile = "/home/ozgun/Downloads/Top20HCC1954TS48h3D.sif";
		String dataFile = "/home/ozgun/Downloads/Top20HCC1954TS48h3D.txt";
		String gcSIFFile = "/home/ozgun/Downloads/Top20HCC1954TS48h3D-gc.sif";
		String formatFile = "/home/ozgun/Downloads/Top20HCC1954TS48h3D-gc.format";

		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(gcSIFFile));
		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(formatFile));
		writer2.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer2.write("node\tall-nodes\tbordercolor\t50 50 50");

		Files.lines(Paths.get(dcSIFFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String[] ab1 = RPPAIDMapper.get().getPlatformLine(t[0]);
			String rel = t[1];
			String[] ab2 = RPPAIDMapper.get().getPlatformLine(t[2]);

			for (String source : ab1[1].split(" "))
			{
				for (String target : ab2[1].split(" "))
				{
					FileUtil.writeln(source + "\t" + rel + "\t" + target, writer1);
				}
			}
		});

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});

		Files.lines(Paths.get(dataFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String[] ab = RPPAIDMapper.get().getPlatformLine(t[0]);
			String[] genes = ab[1].split(" ");
			for (int i = 0; i < genes.length; i++)
			{
				String gene = genes[i];

				if (ab.length < 3 || ab[2].isEmpty())
				{
					FileUtil.lnwrite("node\t" + gene + "\tcolor\t" + vtc.getColorInString(Double.valueOf(t[1])),
						writer2);
				}
				else
				{
					String site = ab[2].split(" ")[i].replaceAll("\\|", "_");
					FileUtil.lnwrite("node\t" + gene + "\trppasite\t" + site + "|p|" +
						vtc.getColorInString(Double.valueOf(t[1])) + "|" +
						(ab[3].equals("i") ? "180 0 20" : ab[3].equals("a") ? "0 180 20" : "20 20 20") +
						"|" + t[1], writer2);
				}
			}
		});

		writer1.close();
		writer2.close();
	}

	private static void generateApoptosisSubgraphs() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/CPTAC-LUAD/mutations/TP53/";

		Set<String> apopt = new HashSet<>(Arrays.asList(("BCL2\n" +
			"MCL1\n" +
			"BCL2L1\n" +
			"BCL2L2\n" +
			"BCL2A1\n" +
			"BAD\n" +
			"PMAIP1\n" +
			"HRK\n" +
			"BCL2L11\n" +
			"BID\n" +
			"BBC3\n" +
			"BAX\n" +
			"BAK1\n" +
			"XIAP\n" +
			"BIRC2\n" +
			"BIRC3\n" +
			"BIRC5\n" +
			"DIABLO").split("\n")));

//		CausalPathSubnetwork.writeGOINeighForCorrBased(dir, apopt, StreamDirection.BOTHSTREAM, "apoptosis");

//		dir = "/Users/ozgun/Documents/Analyses/CPTAC-LUAD/tumor-vs-normal/";
		CausalPathSubnetwork.writeGOINeighForCompBased(dir, apopt, StreamDirection.BOTHSTREAM, "causative-apoptosis-neigh");
	}
}
