package org.panda.misc.analyses;

import org.biopax.paxtools.examples.GOUnificationXREFtoRelationshipXREFConverter;
import org.panda.misc.MutexReader;
import org.panda.utility.FileUtil;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class AnilsPancan
{
	public static final String DIR = "/media/babur/6TB1/TCGA-pancan/";
	public static final String MIX = "Mix";

	public static void main(String[] args) throws IOException
	{
//		generateSampelToTissueMapping();
//		generateComputeDirs();
		generateSIFFormatFileForTP53PancanMutex();
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

		for (MutexReader.Group group : groups)
		{
			if (!group.genes.contains("TP53")) continue;
			if (group.score > 0.01) continue;

			for (String gene : group.genes)
			{
				if (!geneToStudy.containsKey(gene)) geneToStudy.put(gene, group.fromDir);
				else if (!geneToStudy.get(gene).equals(group.fromDir)) geneToStudy.put(gene, MIX);
			}
		}

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
}
