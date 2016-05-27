package org.panda.misc.altmatrix;

import org.panda.misc.PanMutSigReader;
import org.panda.resource.tcga.MutationReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Converts the PanCan MAF to mutation matrix.
 *
 * @author Ozgun Babur
 */
public class PanCanMAFToMatrix
{
	public static void main(String[] args) throws IOException
	{
//		convertToMatrix();
		separateToChunks();
	}


	public static void separateToChunks() throws IOException
	{
		String base = "/home/babur/Documents/mutex/TCGA/PanCan/";
		String wholeDir = base + "1/1";

		for (int i = 2; i <= 20; i++)
		{
			String dir = base + i;
			Files.createDirectories(Paths.get(dir));

			AlterationMatrixSeparator.separate(dir, wholeDir, i);
		}
	}

	public static void convertToMatrix() throws IOException
	{
		String pancanMAF = "/home/babur/Documents/TCGA/PanCan/mutation.maf";
		String outDir = "/home/babur/Documents/mutex/TCGA/PanCan/1/1";

		MutationReader mr = new MutationReader(pancanMAF);
		System.out.println("mr = " + mr);

		Map<String, Double> pvals = PanMutSigReader.getBestPValues();

		List<String> select = pvals.keySet().stream().sorted((s1, s2) -> pvals.get(s1).compareTo(pvals.get(s2)))
			.collect(Collectors.toList());

		String[] samples = mr.getSamples().stream().sorted().collect(Collectors.toList()).toArray(new String[0]);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outDir + "/DataMatrix.txt"));
		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}
		for (String gene : select)
		{
			boolean[] mut = mr.getGeneAlterationArray(gene, samples);
			if (mut != null)
			{
				writer.write("\n" + gene);
				for (boolean m : mut)
				{
					writer.write("\t" + (m ? "1" : "0"));
				}
			}
		}

		writer.close();
	}


}
