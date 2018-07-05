package org.panda.misc.siffile;

import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;
import org.panda.resource.network.SignedPC;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.PhosphoGraph;
import org.panda.utility.statistics.FDR;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Ozgun Babur
 */
public class SignedSIFQuery
{
	public static void pathsBetweenFromPC(Set<String> genes, String outFile) throws IOException
	{
		Map<SignedType, DirectedGraph> graphs = SignedPC.get().getAllGraphs();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (SignedType type : graphs.keySet())
		{
			DirectedGraph graph = graphs.get(type);

			for (String source : genes)
			{
				Set<String> down = new HashSet<>(graph.getDownstream(source));
				down.retainAll(genes);

				for (String target : down)
				{
					writer.write(source + "\t" + type.getTag() + "\t" + target + "\t" +
						graph.getMediatorsInString(source, target));

					if (graph instanceof PhosphoGraph && ((PhosphoGraph) graph).hasSites(source, target))
					{
						writer.write("\t" + CollectionUtil.merge(((PhosphoGraph) graph).getSites(source, target), ";"));
					}

					writer.write("\n");
				}
			}
		}

		writer.close();
	}

	public static void neighborhoodFromFile(String inFile, Set<String> genes, String outFile) throws IOException
	{
		Set<String> rows = Files.lines(Paths.get(inFile)).filter(l -> l.contains("\t")).collect(Collectors.toSet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (String row : rows)
		{
			String[] t = row.split("\t");
			if (t.length > 3 && (genes.contains(t[0]) || genes.contains(t[2])))
			{
				writer.write(row + "\n");
			}
		}

		writer.close();

		String outFmt = outFile.substring(0, outFile.lastIndexOf(".")) + ".format";
		File oFFile = new File(outFmt);
		if (oFFile.exists()) oFFile.delete();

		File fmtFile = new File(inFile.substring(0, inFile.lastIndexOf(".")) + ".format");
		if (fmtFile.exists())
		{
			FileUtil.copyFile(fmtFile.getPath(), outFmt);
		}
	}

	public static void pathsBetweenFromFile(String inFile, Set<String> genes, String outFile) throws IOException
	{
		Set<String> rows = Files.lines(Paths.get(inFile)).filter(l -> l.contains("\t")).collect(Collectors.toSet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (String row : rows)
		{
			String[] t = row.split("\t");
			if (t.length > 2 && (genes.contains(t[0]) && genes.contains(t[2])))
			{
				writer.write(row + "\n");
			}
		}

		writer.close();
	}

	public static void streamOverGenes(String inFile, Set<String> genes, String outFile, int dwDist, int upDist) throws IOException
	{
		Map<String[], String> rows = Files.lines(Paths.get(inFile)).filter(l -> l.contains("\t"))
			.collect(Collectors.toMap(l -> l.split("\t"), l -> l));

		Set<String> upBorder = new HashSet<>(genes);
		Set<String> dwBorder = new HashSet<>(genes);
		Set<String> visitedUp = new HashSet<>();
		Set<String> visitedDw = new HashSet<>();

		Set<String[]> result = new HashSet<>();

		while (dwDist > 0 && !dwBorder.isEmpty())
		{
			Set<String> nextRound = new HashSet<>();
			for (String[] t : rows.keySet())
			{
				if (dwBorder.contains(t[0]) && !visitedUp.contains(t[2]) && !visitedDw.contains(t[2]) && !dwBorder.contains(t[2]))
				{
					result.add(t);
					nextRound.add(t[2]);
				}
			}
			visitedDw.addAll(dwBorder);
			nextRound.removeAll(visitedDw);
			dwBorder = nextRound;
			dwDist--;
		}

		while (upDist > 0 && !upBorder.isEmpty())
		{
			Set<String> nextRound = new HashSet<>();
			for (String[] t : rows.keySet())
			{
				if (upBorder.contains(t[2]) && !visitedUp.contains(t[0]) && !visitedDw.contains(t[0]) && !upBorder.contains(t[0]))
				{
					result.add(t);
					nextRound.add(t[0]);
				}
			}
			visitedUp.addAll(upBorder);
			nextRound.removeAll(visitedUp);
			upBorder = nextRound;
			upDist--;
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (String[] t : result)
		{
			writer.write(rows.get(t) + "\n");
		}

		writer.close();
	}

	public static Set<String> getNodesNameStarting(String sifFile, String... prefix) throws IOException
	{
		return getNodeStreamInSIFFile(sifFile).filter(g -> startsWithOneOfThese(g, prefix)).collect(Collectors.toSet());
	}

	public static boolean startsWithOneOfThese(String s, String... q)
	{
		return Arrays.stream(q).filter(s::startsWith).findAny().isPresent();
	}

	public static Stream<String> getNodeStreamInSIFFile(String file) throws IOException
	{
		return Files.lines(Paths.get(file)).map(l -> l.split("\t")).filter(t -> t.length > 2)
			.map(t -> new String[]{t[0], t[2]}).flatMap(Arrays::stream);
	}

	public static List<String> chooseRelevantLinesFromFormatFile(String formatFile, Set<String> dataIDs) throws IOException
	{
		List<String> chosen = new ArrayList<>();
		Set<String> totals = dataIDs.stream().filter(id -> !id.contains("-")).collect(Collectors.toSet());

		Files.lines(Paths.get(formatFile)).forEach(l ->
		{
			String[] t = l.split("\t");
			if (t[1].equals("all-nodes") || t[1].equals("all-edges"))
			{
				chosen.add(l);
			}
			else if (t[2].equals("borderwidth") || t[2].equals("bordercolor") || t[2].equals("color"))
			{
				chosen.add(l);
			}
			else if (t[2].equals("rppasite"))
			{
				String[] tt = t[3].split("\\|");
				if (dataIDs.contains(tt[0]))
				{
					chosen.add(l);
				}
			}
		});
		return chosen;
	}

	public static Set<String> getDataIDsFromDataCentricSIF(String dataCentricFile, String subsetFile) throws IOException
	{
		Set<String> pairs = getGenePairs(subsetFile);
		Set<String> ids = new HashSet<>();
		Files.lines(Paths.get(dataCentricFile)).map(l -> l.split("\t")).forEach(t ->
		{
			String key = extractGeneFromDataID(t[0]) + " " + extractGeneFromDataID(t[2]);
			if (pairs.contains(key))
			{
				ids.add(t[0]);
				ids.add(t[2]);
			}
		});
		return ids;
	}

	public static String extractGeneFromDataID(String id)
	{
		String[] t = id.split("-");
		if (t.length == 1) return id;
		if ((t[1].startsWith("S") || t[1].startsWith("T") || t[1].startsWith("Y")) && t[1].length() > 1 && Character.isDigit(t[1].charAt(1)))
		{
			return t[0];
		}
		if (t[1].equals("rna") || t[1].equals("cna") || t[1].equals("mut"))
		{
			return t[0];
		}

		return t[0] + "-" + t[1];
	}

	public static Set<String> getGenePairs(String sifFile) throws IOException
	{
		return Files.lines(Paths.get(sifFile)).map(l -> l.split("\t")).map(t -> t[0] + " " + t[2])
			.collect(Collectors.toSet());
	}

	static Set<String> breastCancerGenes = new HashSet<>(Arrays.asList(
		("AKT1\n" +
		"APOBEC3B\n" +
		"ARID1A\n" +
		"ARID1B\n" +
		"BAP1\n" +
		"BARD1\n" +
		"BRCA1\n" +
		"BRCA2\n" +
		"BRIP1\n" +
		"CASP8\n" +
		"CCND1\n" +
		"CDH1\n" +
		"CDKN1B\n" +
		"CHEK2\n" +
		"CTCF\n" +
		"EP300\n" +
		"ERBB2\n" +
		"ESR1\n" +
		"ETV6\n" +
		"FBLN2\n" +
		"FEN1\n" +
		"FLNA\n" +
		"FOXA1\n" +
		"GATA3\n" +
		"IRS4\n" +
		"KEAP1\n" +
		"MAP2K4\n" +
		"MAP3K1\n" +
		"MAP3K13\n" +
		"NCOR1\n" +
		"NOTCH1\n" +
		"NTRK3\n" +
		"PALB2\n" +
		"PBRM1\n" +
		"PIK3CA\n" +
		"POLQ\n" +
		"PPM1D\n" +
		"RB1\n" +
		"SALL4\n" +
		"SMARCD1\n" +
		"TBX3\n" +
		"TP53\n" +
		"ZMYM3").split("\n")));

	public static void main(String[] args) throws IOException
	{
		// For Anil
//		String dir = "/home/babur/Documents/Analyses/JQ1/";
//		Set<String> genes = Files.lines(Paths.get(dir + "platform.txt")).skip(1)
//			.map(l -> l.split("\t")[2].split(" ")).flatMap(Arrays::stream).collect(Collectors.toSet());
//		pathsBetweenFromPC(genes, dir + "base-graph.sif");

//		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/correlation-based-phospho-0.1/";
//		Set<String> genes = Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient115/Galaxy156-[OncotatorMAF].maf"))
//			.filter(l -> !l.startsWith("#")).skip(1).map(l -> l.split("\t"))
//			.filter(t -> t.length > 8 && (t[8].equals("Missense_Mutation") || t[8].equals("Nonsense_Mutation")|| t[8].equals("Splice_Site") || t[8].startsWith("In_Frame") || t[8].startsWith("Frame_Shift")))
//			.map(t -> t[0]).distinct().peek(System.out::println).collect(Collectors.toSet());
//		neighborhoodFromFile(dir + "causative.sif", genes, dir + "causative-115-wes.sif");

//		String dir = "/home/babur/Documents/Analyses/Aslan/platelet-second-pass/fdr0.1-nositematch/";
//		neighborhoodFromFile(dir + "causative.sif", Collections.singleton("PRKCQ"), dir + "causative-PRKCQ.sif");

		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer77/subtypes/PAM50/LuminalAB-vs-Basal-like/";
//		Set<String> genes = loadSignificantGenes(dir + "significance-pvals.txt");
		Set<String> genes = new HashSet<>(Arrays.asList("ESR1"));

//		Set<String> genes = new HashSet<>(OncoKB.get().getAllSymbols());
//		genes.addAll(CancerGeneCensus.get().getAllSymbols());

		streamOverGenes(dir + "causative.sif", genes, dir + "causative-ESR1-neighborhood.sif", 5, 0);
		FileUtil.writeLinesToFile(chooseRelevantLinesFromFormatFile(dir + "causative.format",
			getDataIDsFromDataCentricSIF(dir + "causative-data-centric.sif", dir + "causative-ESR1-neighborhood.sif")),
			dir + "causative-ESR1-neighborhood.format");

//		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/correlation-based/";
//		String inFile = dir + "causative-data-centric.sif";
//		pathsBetweenFromFile(inFile, getNodesNameStarting(inFile, "BRAF", "AKT", "PPP"), dir + "causative-data-centric-subset.sif");

	}

	private static Set<String> loadSignificantGenes(String sigFile) throws IOException
	{
		Map<String, Double> pvals = Files.lines(Paths.get(sigFile)).skip(2).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));

		List<String> select = FDR.select(pvals, null, 0.1);
		return new HashSet<>(select);
	}
}
