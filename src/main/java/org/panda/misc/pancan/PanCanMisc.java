package org.panda.misc.pancan;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Small things about pan-cancer that does not worth to generate a class for it.
 * @author Ozgun Babur
 */
public class PanCanMisc
{
	public static void main(String[] args) throws IOException
	{
		orderByMutationRate();
	}
	public static void orderByMutationRate() throws IOException
	{
		String base = "/home/babur/Documents/PanCan/tissue-unnormalized-results/";
		List<String> genes = Files.lines(Paths.get(base + "pancan.txt"))
			.skip(1).map(l -> l.split("\t")[0]).collect(Collectors.toList());

		genes = GOAssociationAnalysis.getSortedToMutationFreq(genes);

		Set<String> canGen = PanCanResultLoader.readCancerGenes();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(base + "sorted-to-freq.txt"));
		writer.write("Gene");
		for (String gene : genes)
		{
			writer.write("\n" + gene + "\t");
			if (canGen.contains(gene)) writer.write("X");
		}
		writer.close();
	}
}
